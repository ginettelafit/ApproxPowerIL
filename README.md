# ApproxPowerIL

A Shiny application and R package to perform power analysis to select the number of persons for multilevel models with autocorrelated errors using asymptotic approximations of the information matrix

The repository contains functions used in the following [preprint](https://psyarxiv.com/mnce4/):

Lafit... Complete

Users can download the app and run locally on their computer by executing the following commands in R or Rstudio. 

```
# Check if R packages are installed
list.of.packages = c("tidyverse","gridExtra","formattable","htmltools",
"shiny","DT","ggplot2","plyr","dplyr","tidyr","shinyjs","shinythemes","viridis","ploty","remotes")
new.packages = list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load packages
library(tidyverse)
library(gridExtra)
library(formattable)
library(htmltools)
library(shiny)
library(DT)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(shinyjs)
library(shinythemes)
library(viridis)
library(plotly)
library(remotes)

# Install package from GitHub repository
remotes::install_github("ginettelafit/ApproxPowerIL", force = T)

# Load package ApproxPowerIL
library(ApproxPowerIL)

# Using Gist: users can launch this app with:
shiny::runGist('302737dc046b89b7f09d15843389161c')

```

## How the app works in a nutshell

The shiny app focuses on a set of research questions regarding multilevel models with autocorrelated errors. The application allows conducting power analysis to select the number of persons using asymptotic approximations of the information matrix of the standard errors of fixed effects estimates for five two-level models commonly used in intensive longitudinal designs. For these models, when the number of persons and the number of repeated measurement occasions per person is sufficiently large we can obtain a closed-form expression of the standard error of the fixed effect estimates by deriving and approximating the information matrix of the fixed effect estimates of two-level models with autocorrelated errors. Using this closed-form expression, we use the analytic approach to compute statistical power for fixed effects in multilevel models with AR(1) errors, without having to fit the model to large amounts of simulated datasets.


We distinguish five different two-level models where the within-person errors follow an AR(1) process: 

* Model 1 includes a Level 2 or time-invariant continuous predictor that allows investigating changes in the mean level of an outcome due to changes in the person-level predictor. 
* Model 2 includes a Level 2 binary predictor to estimate differences in the mean level of the outcome between two groups of persons.
* Model 3 estimates the fixed within-person association between a Level 1 or time-varying predictor and the outcome. 
* Model 4 extends model 3 by including a Level 2 continuous predictor to investigate cross-level interaction effects between the time-varying predictor and the time-invariant predictor. 
* Model 5 extends the third model by including a binary Level 2 predictor to estimate group differences in the within-persons associations between the time-varying predictor and the outcome. 


[Lafit et al. (2022)](https://psyarxiv.com/mnce4/) includes a the analytic derivation of the information matrix of Models 1 to 5.

