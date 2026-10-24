use v6.d;

unit module Math::QuantileRegression;

use Data::Transformers; # for &rescale and &standardize
use Math::Optimization;
use Math::Polynomial::Chebyshev;
use Math::SparseMatrix;
use Math::Splines;
use Statistics::Distributions; # for &quantile

# The implementations in this pacakge closely follow the Wolfram Language (aka Mathematica)
# implementations of the package "QuantileRegression" -- see:
# - https://github.com/antononcube/WL-QuantileRegression-paclet
# - https://resources.wolframcloud.com/PacletRepository/resources/AntonAntonov/QuantileRegression/

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

        -> +@args {
            my $value = [+] @basis.kv.map(
                -> $index, &basis-function { @coefficients[$index] * basis-function(|@args) }
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
    die 'The data must be a non-empty positional collection.'
    unless $data ~~ Positional:D && $data.elems;

    # A response vector is interpreted as observations at successive integer
    # predictor values. Matrix data uses its first two columns, as in the
    # spline formulation this routine implements.
    my @prepared-data;
    if $data[0] ~~ Numeric:D {
        for $data.list.kv -> $index, $response {
            die 'A response vector must contain only numeric values.' unless $response ~~ Numeric:D;
            @prepared-data.push([$index, $response]);
        }
    }
    else {
        for $data.list -> $row {
            die 'Each data row must contain numeric predictor and response values.'
                unless $row ~~ Positional:D && $row.elems >= 2
                    && $row[0] ~~ Numeric:D && $row[1] ~~ Numeric:D;
            @prepared-data.push([$row[0], $row[1]]);
        }
    }

    my @x = @prepared-data.map(*[0]);
    die 'Spline regression requires at least two distinct predictor values.'
    unless @x.min < @x.max;

    my @requested-knots = do given $knots {
        when Int:D {
            die 'The number of spline pieces must be positive.' unless $_ > 0;
            (0 .. $_).map({ @x.min + (@x.max - @x.min) * $_ / $knots }).Array
        }
        when Positional:D {
            die 'The knots must be numeric.' unless $_.all ~~ Numeric:D;
            $_.list.sort.Array
        }
        default {
            die 'The knots must be a positive integer or a positional collection of numbers.'
        }
    };
    die 'At least two knots must be supplied.' unless @requested-knots.elems >= 2;

    return $method ~~ Str:D
        ?? lp-quantile-regression(@prepared-data, @requested-knots, $probabilities.Array,
            $interpolation-order, :$method)
        !! lp-quantile-regression(@prepared-data, @requested-knots, $probabilities.Array,
            $interpolation-order);
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
    die 'The data must contain at least one row.' unless @data.elems;
    die 'At least two knots must be supplied.' unless @knots.elems >= 2;
    die 'Probabilities must be between zero and one.'
    unless @probs.elems && @probs.all ~~ { $_ >= 0 && $_ <= 1 };

    my @sorted-knots = @knots.sort;
    die 'Spline knots must span a non-zero interval.' unless @sorted-knots.head < @sorted-knots.tail;

    # Clamp the endpoints with additional copies, as in the reference
    # implementation. Math::Splines expects the spline degree (not order + 1).
    my @clamped-knots =
        |(@sorted-knots.head xx $order),
        |@sorted-knots,
        |(@sorted-knots.tail xx $order);
    
    my $basis-max-index = @clamped-knots.elems - $order - 2;

    die 'The supplied knots and interpolation order do not define a B-spline basis.'
    if $basis-max-index < 0;

    # Math::Splines handles the right endpoint for a domain ending at one.
    # A linear change of parameter leaves the B-spline basis unchanged while
    # making its endpoint convention work for arbitrary user knots.
    my $knot-min = @clamped-knots.head;
    my $knot-max = @clamped-knots.tail;
    my @unit-knots = @clamped-knots.map({ ($_ - $knot-min) / ($knot-max - $knot-min) }).Array;
    my @basis = (0 .. $basis-max-index).map(-> $index {
        -> $x {
            my $unit-x = ($x - $knot-min) / ($knot-max - $knot-min);
            b-spline-basis($order, $index, $unit-x, knots => @unit-knots)
        }
    }).Array;

    die 'The supplied knots and interpolation order produced no B-spline basis functions.'
    unless @basis.elems;

    # The basis-fit implementation creates the residual-variable LP, performs
    # the response scaling from the reference algorithm, and returns weighted
    # basis callables.
    return lp-quantile-regression-fit(@data, @basis, @probs, |%args)
}
