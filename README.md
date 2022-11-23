
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Monte Carlo Goodness of Fit tests for Stochastic Block Models

## Intended use of the package

<!-- badges: start -->

There are various statistical models exist in the literature to fit
network data. But a natural statistical question arises whether the fit
is good or bad. A chi-square goodness of fit test conditioning on
sufficient statistics was proposed (Karwa et al.) for several popular
stochastic block models (SBM). My package now focuses only on a special
case of those i.e, ERSBM. Also, we have mainly focused on a simple
undirected, unweighted network with no self-loop. The main function
`goftest` take the adjacency matrix and the block assignment of each
node as input. It fits the ERSBM model on the corresponding graph of the
adjacency matrix and after performing a chi-square goodness of fit test,
gives the value of the test statistic and p-value as outputs.  
<!-- badges: end -->

## Installation

To install the package “MCGoFSBM”, the library “devtools” need to be
installed there. Thereafter, You can install the development version of
the package from [GitHub](https://github.com/GhoshSoham/MCGoFSBM) using
the following commands in the console:

``` r
# install.packages("devtools")
devtools::install_github("https://github.com/GhoshSoham/MCGoFSBM")
```

## Left for the project

Till now, I have assumed that the block assignment of each node is known
to users. I am planning to update the function when the block assignment
is not known to the user but has an idea about how many blocks are
there. As the code is not optimized fully, I am planning to implement
some loops in cpp. Apart from that, I am yet to add some more
compatibility checks, formally cite the references, and run some formal
tests to ensure correctness. I also have planned to add some examples
with the functions in the documentation.

## Other

When I submitted the proposal, I proposed that I will focus on different
stochastic block models discussed in (Karwa et al.) with known block
assignment rather than the case when block assignment is not known for
the nodes. But later I thought it will be better to focus on a
particular model i.e, ERSBM now. If time permits, I will upgrade the
package for other block models discussed in the paper.

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
