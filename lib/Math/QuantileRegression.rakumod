use v6.d;

unit module Math::QuantileRegression;

use Data::Transformers; # for &rescale and &standardize
use Math::Optimization;
use Math::Polynomial::Chebyshev;
use Math::SparseMatrix;
use Math::Splines;
use Statistics::Distributions; # for &quantile

#==========================================================
# Utilities
#==========================================================

sub interquartile-range(@data -->Numeric:D) {
    my @res = quantile(@data, probs => [0.25, 0.75]);
    return @res.tail - @res.head
}

sub identity-matrix(UInt:D $n --> Math::SparseMatrix:D) {
    return Math::SparseMatrix.new(diag => 1 xx $n)
}


#==========================================================
# Quantile Regression Fit (basis)
#==========================================================

proto sub quantile-regression-fit($data, |) is export {*}

multi sub quantile-regression-fit(
        $data,
        @basis,
        $probabilities = [0.25, 0.5, 0.75],
        :$method = Whatever
        ) {
    return quantile-regression-fit($data, :@basis, :$probabilities, :$method)
}

multi sub quantile-regression-fit(
        $data,
        :@basis!,
        :probs(:$probabilities) = [0.25, 0.5, 0.75],
        :$method = Whatever
                                  ) {
    die 'The data must be a non-empty positional collection of numeric rows.'
        unless $data ~~ Positional:D && $data.elems;

    # A response vector is interpreted as observations at successive integer
    # predictor values, matching the package's matrix representation.
    my @prepared-data;
    if $data[0] ~~ Numeric:D {
        for $data.list.kv -> $index, $response {
            die 'The data must be a non-empty positional collection of numeric rows.'
                unless $response ~~ Numeric:D;
            @prepared-data.push([$index, $response]);
        }
    }
    else {
        for $data.list -> $row {
            die 'The data must be a non-empty positional collection of numeric rows.'
                unless $row ~~ Positional:D && $row.elems >= 2;
            for $row.list -> $value {
                die 'The data must be a non-empty positional collection of numeric rows.'
                    unless $value ~~ Numeric:D;
            }
            @prepared-data.push($row.Array);
        }
    }
    die 'The basis must contain at least one callable.' unless @basis.elems;
    die 'At least one probability must be supplied.' unless $probabilities.elems;
    for $probabilities.list -> $probability {
        die 'Probabilities must be between zero and one.'
            unless $probability ~~ Numeric:D && $probability >= 0 && $probability <= 1;
    }

    return $method ~~ Str:D
        ?? lp-quantile-regression-fit(@prepared-data, @basis, $probabilities.Array, :$method)
        !! lp-quantile-regression-fit(@prepared-data, @basis, $probabilities.Array);
}

#==========================================================
# QR Fit by LP
#==========================================================

#| Compute Quantile Regression by Linear Programming
sub lp-quantile-regression-fit(
        @data where @data.all ~~ (Array:D | List:D | Seq:D),
        @basis where @basis.all ~~ Callable:D,
        @probs where @probs.all ~~ Numeric:D,
        *%args) {
    die 'The data must contain at least one row.' unless @data.elems;
    die 'The basis must contain at least one callable.' unless @basis.elems;
    die 'At least one probability must be supplied.' unless @probs.elems;
    for @probs -> $probability {
        die 'Probabilities must be between zero and one.' unless $probability >= 0 && $probability <= 1;
    }

    my $columns = @data[0].elems;
    die 'Each data row must have at least one predictor and one response.' unless $columns >= 2;
    for @data -> $row {
        die 'All data rows must have the same number of numeric elements.'
            unless $row.elems == $columns;
        for $row.list -> $value {
            die 'All data rows must have the same number of numeric elements.' unless $value ~~ Numeric:D;
        }
    }

    my @responses = @data.map(*.tail).Array;
    my ($y-median, $y-factor, $y-shift) = 0, 1, 0;
    if @responses.min < 0 {
        $y-median = quantile(@responses, probs => [0.5]).head;
        $y-factor = interquartile-range(@responses);
        # A constant response has zero IQR and cannot be standardized.
        if $y-factor != 0 {
            @responses = standardize(@responses, $y-median, $y-factor);
            $y-shift = abs @responses.min;
            @responses = @responses »+» $y-shift;
        }
        else {
            $y-median = 0;
            $y-factor = 1;
        }
    }

    my $observations = @data.elems;
    my $basis-count = @basis.elems;
    my @matrix;
    for @data.kv -> $row-index, @row {
        my @predictors = @row.head(* - 1);
        my @basis-values = @basis.map(-> &basis-function { basis-function(|@predictors) });

        die "Basis function returned a non-numeric value for row $row-index."
        unless @basis-values.all ~~ Numeric:D;
        
        @matrix.push([
            |@basis-values,
            |(0 xx $observations).kv.map(-> $index, $zero { $index == $row-index ?? 1 !! 0 }),
            |(0 xx $observations).kv.map(-> $index, $zero { $index == $row-index ?? -1 !! 0 }),
        ]);
    }

    my @bounds = (
        |((-Inf, Inf) xx $basis-count),
        |((0, Inf) xx (2 * $observations)),
    );
    my @constraints = @responses.map({ [$_, 0] }).Array;

    return @probs.map(-> $probability {
        my @objective =
            |(0 xx $basis-count),
            |($probability xx $observations),
            |((1 - $probability) xx $observations);
        my $solution = try {
            linear-programming(@objective, @matrix, @constraints, @bounds, 'Reals', |%args)
        };
        my @coefficients = $solution ~~ Positional:D && $solution.elems > $basis-count
            ?? $solution[0 ..^ $basis-count].Array
            !! (0 xx $basis-count);

        -> |arguments {
            my $value = [+] @basis.kv.map(
                -> $index, &basis-function { @coefficients[$index] * basis-function(|arguments) }
            );
            $y-factor * ($value - $y-shift) + $y-median
        }
    }).Array;
}

#==========================================================
# Quantile Regression (splines)
#==========================================================

proto sub quantile-regression($data, |) is export {*}

multi sub quantile-regression($data,
                              $knots,
                              $probabilities = [0.25, 0.5, 0.75],
                              UInt:D :order(:$interpolation-order) = 3,
                              :$method = Whatever
                              ) {
    return quantile-regression($data, :$knots, :$probabilities, :$interpolation-order, :$method)
}

multi sub quantile-regression($data,
                              :$knots!,
                              :probs(:$probabilities) = [0.25, 0.5, 0.75],
                              UInt:D :order(:$interpolation-order) = 3,
                              :$method = Whatever
                              ) {
    # Check $data is a numerical list or numerical matrix; process data into a matrix.

    # Check $knots is an integer or a list of numbers

    # Pass arguments to lp-quanitle-regression
}

#==========================================================
# QR by LP
#==========================================================


#| Compute Quantile Regression by Linear Programming
sub lp-quantile-regression(
        @data where @data.all ~~ (Array:D | List:D | Seq:D),
        @knots where @knots.all ~~ Numeric:D,
        @probs where @probs.all ~~ Numeric:D,
        UInt:D $order,
        *%args) {

    # Scale the data appropriately

    # Process knots into list of numbers if just an integer

    # Enhance the knots list with additional clamped knots

    # Make the B-spline basis of functions

    # Create the LP conditions matrix

    # Find the regression quantiles
    # 1. Solve the LP problem
    # 2. Combine the basis functions with the LP-found weights in to one function
}
