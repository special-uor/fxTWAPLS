Benchmarks
================

## `cv.w`

``` r
times_df <- reconbio::par_benchmark(CPUS = c(1, 2, 4, 8, 12),
                                    FUN = reconbio::cv.w, 
                                    plot = TRUE, 
                                    modern_tax = taxa, 
                                    modern_climate = modern_pollen$Tmin, 
                                    nPLS = 5, 
                                    trainfun = reconbio::WAPLS.w, 
                                    predictfun = reconbio::WAPLS.predict.w)
```

| CPUs | Time \[s\] |
| ---: | ---------: |
|    1 |   1318.007 |
|    2 |    693.736 |
|    4 |    378.848 |
|    8 |    259.177 |
|   12 |    240.585 |

![](inst/figures/benchmarks-cv.w-output-1.png)<!-- -->

## `cv.pr.w`

``` r
times_df <- reconbio::par_benchmark(CPUS = c(1, 2, 4, 8, 12),
                                    FUN = reconbio::cv.pr.w, 
                                    plot = TRUE, 
                                    modern_tax = taxa, 
                                    modern_climate = modern_pollen$Tmin, 
                                    nPLS = 5, 
                                    trainfun = reconbio::WAPLS.w, 
                                    predictfun = reconbio::WAPLS.predict.w,
                                    pseduo = pseduo_Tmin)
```

| CPUs | Time \[s\] |
| ---: | ---------: |
|    1 |   1465.608 |
|    2 |    756.902 |
|    4 |    413.547 |
|    8 |    326.225 |
|   12 |    231.002 |

![](inst/figures/benchmarks-cv.pr-w-output-1.png)<!-- -->

## `get_pseudo`

``` r
times_df <- reconbio::par_benchmark(CPUS = c(1, 2, 4, 8, 12),
                                    FUN = reconbio::get_pseduo, 
                                    plot = TRUE, 
                                    dist = dist, 
                                    x = modern_pollen$Tmin)
```

| CPUs | Time \[s\] |
| ---: | ---------: |
|    1 |    581.771 |
|    2 |    292.639 |
|    4 |    164.557 |
|    8 |    119.054 |
|   12 |    119.959 |

![](inst/figures/benchmarks-get_pseudo-pseudo-1.png)<!-- -->
