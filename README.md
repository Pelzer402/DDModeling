

# DDModeling

<!-- badges: start -->
<!-- badges: end -->
DDModeling is an R packet. Its goal is to give researchers access to an easy-to-use toolkit for the simulation and adaptation of drift-diffusion models in cognitive psychology.

## Installation

You can install the released version of DDModeling from [GitHub](https://github.com/)

``` r
# install.packages("devtools")
devtools::install_github("Pelzer402/DDModeling")
```
## Example

This is a basic example to show you how to generate a drift diffusion model and conduct a simulation with it


```r
library(DDModeling)
# First specify a model
DSTP_M1 <-  DDModel(model="DSTP",task="flanker",
                    CDF_perc = c(0.1,0.3,0.5,0.7,0.9),CAF_perc = c(0.0,0.2,0.4,0.6,0.8,1.0))
# Now simulate some data
R1 <- Sim_DDModel(model = DSTP_M1,trials = 10000)
```

## More Information
For more information visit the projects website [here](https://pelzer402.github.io/DDModeling/)
