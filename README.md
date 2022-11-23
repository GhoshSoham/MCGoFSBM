
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Monte Carlo Goodness of Fit tests for Stochastic Block Models

## Intended Use of the package

<!-- badges: start -->

There are various statistical model exists in the literature to fit a
network data. But a natural statistical question that arise is how good
is the fit. A chi-square goodness of fit test conditioning on sufficient
statistic was proposed (Karwa et al.) for ERSBMs. My package is based on
that. It will take Adjacency Matrix of an undirected, unweighted graph
with no self loop and also the block assignment of each nodes as input.
It will fit ERSBM model on that and after performing a chi-square
goodness of fit test will give the value of the test statistic of the
observed graph and p-value as an output.  
<!-- badges: end -->

## Installation

In order to install the package “MCGoFSBM”, the library “devtools” need
to be installed there. There after, You can install the development
version of the package from
[GitHub](https://github.com/GhoshSoham/MCGoFSBM) using the following
command in console:

``` r
# install.packages("devtools")
devtools::install_github("GhoshSoham/MCGoFSBM")
```

## Left for the project

## Other

<!-- ## Example -->
<!-- This is a basic example which shows you how to solve a common problem: -->
<!-- ```{r example} -->
<!-- library(MCGoFSBM) -->
<!-- ## basic example code -->
<!-- ``` -->
<!-- What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so: -->
<!-- ```{r cars} -->
<!-- summary(cars) -->
<!-- ``` -->
<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. You could also use GitHub Actions to re-render `README.Rmd` every time you push. An example workflow can be found here: <https://github.com/r-lib/actions/tree/v1/examples>. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
