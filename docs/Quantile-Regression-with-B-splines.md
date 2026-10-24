# Quantile Regression with B-splines

Anton Antonov    
[RakuForPrediction at WordPress](https://rakuforprediction.wordpress.com)   
October 2026

## Introduction

This notebook has explanations and examples of using the function `quantile-regression` provided by the Raku package ["QuantileRegression"](https://github.com/antononcube/Raku-Math-QuantileRegression), [AAp1].

The function `quantile-regression` computes [Quantile Regression (QR)](https://en.wikipedia.org/wiki/Quantile_regression) fits over a time series, a list of numbers or a list of numeric pairs using B-spline bases.

**Remark:** This notebook closely follows the page of the Wolfram Function Repository entry [`QuantileRegression`](https://resources.wolframcloud.com/FunctionRepository/resources/QuantileRegression/), [AAp2].

---

## `quantile-regression`

`quantile-regression(data, knots, probs)`    
does quantile regression over the times series or data array data using the knots specification knots for the probabilities probs.

`quantile-regression(data, knots, probs, opts)`    
does quantile regression with the options opts.

### Details & Options

- `quantile-regression` works on numpy arrays, lists of numbers, and lists of numeric pairs.


- The curves computed with quantile regression are called **regression quantiles**.


- The regression quantiles corresponding to the specified probabilities are linear combinations of B-splines generated over the specified knots.


- In other words, `quantile-regression` computes fits using a B-spline functions basis. The basis is specified with the `knots` argument and the option `order`.


- The named argument `order` is 3 by default.


---

## Setup

```raku
use Data::Generators;
use Statistics::Distributions;
use Data::TypeSystem;

use Math::QuantileRegression;

use Text::Plot;
use JavaScript::Google::Charts;
```

```raku
my $backgroundColor = '#1A1A1A';
my $format = 'html'; 
my $titleTextStyle = %(color => 'Gray', fontSize => 16);
my $chartArea = {
    :$backgroundColor,
    left => 50,
    right => 100,
    top => 100, 
    bottom => 50
}
```

```raku
sub qr-plot(
        @data where @data.all ~~ Positional:D,
        @q-funcs where @q-funcs.all ~~ Callable:D,
        @probs where @probs.all ~~ Numeric:D,
        Int:D :$width = 800,
        Int:D :$height = 600,
        :$pointColor = 'SteelBlue',
        :$format = 'jupyter',
        *%args
        --> Str:D
        ) {
    die 'There must be one quantile function for each probability.'
    unless @q-funcs.elems == @probs.elems;
    die 'Data must consist of numeric predictor-response pairs.'
    unless @data.grep(-> $row {
            $row.elems >= 2 && $row[0] ~~ Numeric:D && $row[1] ~~ Numeric:D
    }).elems == @data.elems;

    my @column-names = 'x', 'data', |@probs.map(*.Str);
    my @records = @data.map(-> $row {
        my $x = $row[0];
        my @values = $x, $row[1], |@q-funcs.map(-> &fit { &fit($x) });
        (@column-names.Array Z=> @values).Hash
    }).Array;

    my %series =
        0 => { type => 'scatter', pointSize => 4, lineWidth => 0, color => $pointColor},
        |(1 .. @q-funcs.elems).map(-> $index {
            $index => { type => 'line', pointSize => 0, lineWidth => 2 }
        });

    js-google-charts('Combo', @records,
        column-names => @column-names,
        :$format, :$width, :$height,
        title => 'Quantile regression with B-splines',
        hAxis => { title => 'X Axis Label' },
        vAxis => { title => 'Y Axis Label' },
        seriesType => 'line',
        series => %series,
        :!png-button,
        |%args
    )
}
```

---

## Basic Examples

Make a random signal:

```raku
srand(3232);
my $n = 200;
my @rand-data = (^$n) Z random-real([0, 100], $n);

deduce-type(@rand-data)
```

Compute `quantile-regression` with five knots for the probabilities $0.25$ and $0.75$:

```raku
my @q-funcs = quantile-regression(@rand-data, knots => 5, probs => [0.25, 0.75]);

@q-funcs.elems
```

Here is a plot of the original data and the obtained regression quantiles:

```raku
#% html
qr-plot(@rand-data, @q-funcs, [0.25, 0.75], 
    :$backgroundColor,
    :$titleTextStyle,
    :$chartArea,
    :$format, 
    div-id => 'first');
```

Find the fraction of the data points that are under the second regression quantile:

```raku
@rand-data.grep({ @q-funcs[1]($_.head) > $_.tail }).elems / @rand-data.elems
```

---

## Scope

Here is a quantile regression computation over a numerical vector:

```raku
#% html

my @vec-data = ((1...200) <</>> 6)».sin <<+>> random-real([-0.5, 0.5], 200);

say deduce-type(@vec-data);

my @q-funcs = quantile-regression(@vec-data, :12knots, probs => 0.5);
my @vec-data2 = (^@vec-data.elems) Z @vec-data;
say deduce-type(@vec-data2);

qr-plot(@vec-data2, @q-funcs, [0.5,], 
    :$backgroundColor,
    :$titleTextStyle,
    :$chartArea,
    legend => {
        position => 'right',
        textStyle => {
            color => '#dad5d5',
            fontSize => 12
        }
    },
    colors => [
        '#4285F4',
        '#EA4335'
    ],
    :$format, 
    div-id=>'scope-1'
);
```

The second argument—the knots specification—can be an integer specifying the number of knots or a list of numbers specifying the knots of the B-spline basis:

```raku
my @q-funcs2 = quantile-regression(@rand-data, knots => @rand-data.map(*.head).pick(5).sort, probs => 0.5);
@q-funcs2[0](100)
```

---

## Options

### Order (interpolation order)

```raku
#% html
my @q-funcs = [0, 1, 3].map({ quantile-regression(@rand-data, :5knots, probs => 0.5, order => $_).head });

qr-plot(@rand-data, @q-funcs, [0, 1, 3], 
    title => 'Interpolation orders',
    :$titleTextStyle, 
    :$chartArea, 
    legend => {
        position => 'right',
        textStyle => {
            color => '#dad5d5',
            fontSize => 12
        }
    }, 
    :$backgroundColor, 
    :$format, 
    div-id => 'order');
```

---

## Applications

### Fit for heteroscedastic data

Here is heteroscedastic data (the variance is not constant with respect to the predictor variable):

```raku
#% html
my @x-values = -3, -2.99 ... 3.01;
my @dist-data = @x-values.map({
    [$_, exp( - $_ ** 2) + random-variate(NormalDistribution.new(0, 0.15 * sqrt(abs(1.5 - $_) / 1.5)))]
});

js-google-charts('Scatter', 
    @dist-data,
    title => 'Distribution data',
    :$titleTextStyle,
    series => { 0 => {color => 'SteelBlue', pointSize => 2} },
    :$backgroundColor, 
    :$chartArea, 
    :$format, 
    div-id => 'dist-data')
```

Find quantile regression fits:

```raku
my @probs = 0.1, 0.25, 0.5, 0.75, 0.9;
my @q-funcs = quantile-regression(@dist-data, :4knots, :@probs);

@q-funcs.elems
```


Plot the data and the regression quantiles:


```raku
#% html
qr-plot(@dist-data, @q-funcs, @probs,
    title => 'Distribution data fits',
    :$titleTextStyle, 
    :$chartArea,
    pointColor => 'Gray',
    legend => {
        position => 'right',
        textStyle => {
            color => '#dad5d5',
            fontSize => 12
        }
    }, 
    :$backgroundColor, 
    :$format, 
    div-id => 'dist-data-fits');
```

Note that the regression quantiles clearly outline the heteroscedastic nature of the data.