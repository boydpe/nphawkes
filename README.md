
<!-- README.md is generated from README.Rmd. Please edit that file -->

# nphawkes

<!-- badges: start -->

<!-- badges: end -->

The nphawkes package can be used to study Hawkes processes using
nonparametric procedures. Data may be applied to the Model Independent
Stochastic Declustering (MISD) algorithm to probabilistically classify
events as background or triggered events. The conditional intensity of
the process can be calculated, and residual analysis may be exectued via
thinning, superpositioning, or super-thinning. This package may also be
used to create plots such as histogram estimators of the triggering
functions, line plots of the conditional intensity, and histograms of
the residual process. Temporal and spatio-temporal Hawkes processes can
be used, along with marked versions of each.

## Installation

<!-- You can install the released version of nphawkes from [CRAN](https://CRAN.R-project.org) with: -->

You can install the resleased version of nphawkes from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("boydpe/nphawkes")
```

## Example

We will demonstrate the use of the package with simulated earthquake
data available within the package.

``` r
library(nphawkes)
load(file = "data/eq.RData")
head(quake)
#>   X       lon         lat         t        m event_type     time       date
#> 1 1 0.6944789 0.698142489  8.966234 3.756671 background 23:11:23 2010-01-09
#> 2 2 0.6834387 0.686690810 10.223219 3.569237 aftershock  5:21:26 2010-01-11
#> 3 3 0.7150710 0.674706658 12.579882 3.038466 aftershock  13:55:2 2010-01-13
#> 4 4 0.9108068 0.629251374 19.361115 3.672588 background   8:40:0 2010-01-20
#> 5 5 0.4389249 0.002338193 23.427427 3.047991 background 10:15:30 2010-01-24
#> 6 6 0.3801922 0.120843447 29.829610 3.042964 background 19:54:38 2010-01-30
```

We will first employ the MISD algorithm to decluster the earthquake
catalog using time, space, and marks (magnitude). The temporal bins
(time\_breaks) will be defined, in terms of days, and we will use
roughly uniform bins for space and marks by setting the `space_quantile`
and `mark_quantile` inputs equal to `TRUE`.

``` r
out = nph(dates = quake$date, ref_date = "2010-01-01",
          lat = quake$lat, lon = quake$lon,
          marks = quake$m,
          time_breaks = c(0, 1, 3, 5, 7, 14, 31, 1500),
          #mark_breaks = c(3, 3.1, 3.25, 3.4, 6),
          space_breaks = c(0, 5, 10, 20, 50, 100),
          time_quantile = FALSE,
          space_quantile = FALSE,
          mark_quantile = TRUE,
          time_of_day = quake$time)

out1 = nph(dates = quake$date, ref_date = "2010-01-01",
          #lat = quake$lat, lon = quake$lon,
          marks = quake$m, 
          time_breaks = c(0, 1, 3, 5, 7, 14, 31, 1500),
          #mark_breaks = c(3, 3.1, 3.25, 3.4, 6),
          time_quantile = FALSE,
          #space_quantile = TRUE,
          mark_quantile = TRUE,
          time_of_day = quake$time)

out2 = nph(dates = quake$date, ref_date = "2010-01-01",
          lat = quake$lat, lon = quake$lon,
          #marks = quake$m,
          time_breaks = c(0, 1, 3, 5, 7, 14, 31, 1500),
          #mark_breaks = c(3, 3.1, 3.25, 3.4, 6),
          time_quantile = FALSE,
          space_quantile = TRUE,
          #mark_quantile = TRUE,
          time_of_day = quake$time)

out4 = nph(dates = quake$date, ref_date = "2010-01-01",
          lat = quake$lat, lon = quake$lon,
          marks = quake$m,
          time_breaks = c(0, 5, 7, 14, 31, 1500),
          #mark_breaks = c(3, 3.1, 3.25, 3.4, 6),
          #space_breaks = c(0, 5, 10, 20, 50, 100),
          time_quantile = FALSE,
          space_quantile = TRUE,
          mark_quantile = TRUE,
          h_length = 4,
          k_length = 4,
          #mark_quantile = FALSE,
          time_of_day = quake$time)
```

We can view the background rate and percentage of mass that lies on the
diagonal of the resulting probability matrix.

``` r
out$br
#> [1] 0.4709993
out$perc_diag
#> [1] 0.9999999
```

We can then inspect histogram estimators of the triggering functions.

``` r
trig_plots(model = out4, g_max = 32,
           k_max = 4, k_min = 3,  
           h_max = 22)
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

    #> TableGrob (1 x 3) "arrange": 3 grobs
    #>   z     cells    name           grob
    #> 1 1 (1-1,1-1) arrange gtable[layout]
    #> 2 2 (1-1,2-2) arrange gtable[layout]
    #> 3 3 (1-1,3-3) arrange gtable[layout]

We’ll now calculate the conditional intensity of the process.

``` r
ci = cond_int(model = out)
head(ci)
#>       times         lat       lon    marks  cond_int       Date
#> 1  8.966238 0.698142489 0.6944789 3.756671 0.4756867 2010-01-09
#> 2 10.223218 0.686690810 0.6834387 3.569237 0.4709993 2010-01-11
#> 3 12.579884 0.674706658 0.7150710 3.038466 0.4709993 2010-01-13
#> 4 19.361111 0.629251374 0.9108068 3.672588 0.4709993 2010-01-20
#> 5 23.427431 0.002338193 0.4389249 3.047991 0.4709993 2010-01-24
#> 6 29.829606 0.120843447 0.3801922 3.042964 0.4709993 2010-01-30
```

We can perfrom super-thinning as a residual method, thinning events in
areas of high conditional intenstiy and superpositioning events in areas
of low intensity. We’ll use the median value of the conditional
intensity as the threshhold value. We’ll use the `sim_grid = TRUE`
argument since the geographical coordinates are simulated within the
\([0\times1]\times[0\times1]\) grid.

``` r
st = super_thin(K = "median_ci", model = out,
                method = "superthin", sim_grid = TRUE)
```
