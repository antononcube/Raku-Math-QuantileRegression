#!/usr/bin/env raku
use v6.d;

use Math::QuantileRegression;
use Data::Importers;
use Data::TypeSystem;

my $url = "https://raw.githubusercontent.com/antononcube/MathematicaVsR/master/Data/MathematicaVsR-Data-Atlanta-GA-USA-Temperature.csv";
#my $url = $*CWD.add(<.. MathematicaVsR Data MathematicaVsR-Data-Atlanta-GA-USA-Temperature.csv>);
say (:$url);
say "\$url as file exist: {$url.f}";
my @dsTemperature = data-import($url, headers => 'auto');
@dsTemperature .= deepmap({ $_ ~~ / ^ \d+ '-' / ?? Date.new($_) !! $_.Num });
@dsTemperature .= grep( 2020 ≤ *<Date>.year ≤ 2023);
say deduce-type(@dsTemperature);

my @temp-data = @dsTemperature.map(*<AbsoluteTime Temperature>);
say deduce-type(@temp-data);

my @probs = 0.01, |(0.1, 0.3 ... 0.9), 0.99;

my $start = now;
my @q-funcs = quantile-regression(@temp-data, :4knots, :@probs);
my $end = now;

say "Quantile regression time {$end - $start} s";

my $start2 = now;
my @variances = @temp-data.map({ @q-funcs[*-1]($_.head) - @q-funcs[0]($_.head)});
my $end2 = now;
say "Variance computation time {$end2 - $start2} s";
