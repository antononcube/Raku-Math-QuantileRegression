#!/usr/bin/env raku
use v6.d;

use Math::QuantileRegression;
use Data::Generators;

my @x-values = -3, -2.99 ... 3.0;
my @dist-data = @x-values.map({
    [$_, exp( - $_ ** 2) + random-variate(NormalDistribution.new(0, 0.15 * sqrt(abs(1.5 - $_) / 1.5)))]
});

my @probs = 0.1, 0.25, 0.5, 0.75, 0.9;

my $start = now;
my @q-funcs = quantile-regression(@dist-data, :4knots, :@probs);
my $end = now;

say "Quantile regression time {$end - $start} s";

my $start2 = now;
my @variances = @dist-data.map({ @q-funcs[*-1]($_.head) - @q-funcs[0]($_.head)});
my $end2 = now;
say "Variance computation time {$end2 - $start2} s";
