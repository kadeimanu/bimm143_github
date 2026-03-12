# Class 4 R Intro
Kaliyah Adei-Manu

## Making basic plots in R

``` r
x <- 1:50
plot (x)
```

![](Untitled_files/figure-commonmark/unnamed-chunk-1-1.png)

Make a Sin Plot

``` r
plot (x, sin(x))
```

![](Untitled_files/figure-commonmark/unnamed-chunk-2-1.png)

Let’s make it nicer:

``` r
plot (x, sin(x), typ="l", col="blue", lwd=2)
```

![](Untitled_files/figure-commonmark/unnamed-chunk-3-1.png)

``` r
log(10, base =10)
```

    [1] 1
