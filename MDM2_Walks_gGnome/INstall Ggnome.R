# install gGnome Package

install.packages('devtools')
install.packages('testthat')

## allows dependencies that throw warnings to install
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS = TRUE)


devtools::install_github('mskilab/gGnome')3
