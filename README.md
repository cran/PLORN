
<!-- README.md is generated from README.Rmd. Please edit that file -->

# PLORN

<!-- badges: start -->
<!-- badges: end -->

The goal of PLORN is to provide the functions to construct a prediction model of environments using noisy omics data linked with the environments based on PLORN algorithm. 

## Installation

You can install the development version of PLORN from
[GitHub](https://github.com/takakoizumi/PLORN) with:  
\#install.packages(“devtools”)  
devtools::install\_github(“takakoizumi/PLORN”)

## Example

install.packages(“PLORN”)  
library(PLORN)  
\# basic example code  
data(Pinus)  
x.train &lt;- Pinus\[\[1\]\]  
x.test &lt;- Pinus\[\[2\]\]  
y &lt;- Pinus\[\[3\]\]  
cor(y, plorn(x.train, y, newx = x.test, n.pred = 100))
