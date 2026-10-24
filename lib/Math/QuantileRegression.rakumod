use v6.d;

unit module Math::QuantileRegression;

use Data::Transformations; # for &rescale
use Math::Optimization;
use Math::Polynomial::Chebyshev;
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
    Math::SparseMatrix.new(diag => 1 xx $n);
}

#==========================================================
# Quantile Regression Fit (basis)
#==========================================================

proto sub quantile-regression-fit($data, |) is export {*}

multi sub quantile-regression-fit(
        $data,
        $basis,
        $probabilities = [0.25, 0.5, 0.75],
        :$method = Whatever
        ) {
    return quantile-regression-fit($data, :$basis, :$probabilities, :$method)
}

multi sub quantile-regression-fit(
        $data,
        :$basis!,
        :probs(:$probabilities) = [0.25, 0.5, 0.75],
        :$method = Whatever
                                  ) {
    # Scale the data appropriately

    # Create the LP conditions matrix

    # Find the regression quantiles
    # 1. Solve the LP problem
    # 2. Combine the basis functions with the LP-found weights in to one function
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
    # ...
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
