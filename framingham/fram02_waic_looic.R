##
## Script name: fram02_waic_looic.R
##
## Purpose of script: calculates DIC and WAIC for the VCM model and the full BCTM for standardized y
##
## Author: Manuel Carlan
##
## Date Created: 2020-10-5
##
## Email: mcarlan@uni-goettingen.de
##
## ---------------------------

library(dplyr)
# rm(list = ls())

source("bctm_utils.R")
source("bctm_design_funs2.R")
source("bctm_design.R")
source("bctm_fun.R")

source("nuts/nuts_utils.R")
source("nuts/nuts.R")

packages <- c("Rcpp", "RcppArmadillo", "RcppEigen", "splines", "mgcv", "Matrix", "MCMCpack", 
              "tidyverse", "profvis",  "tictoc", "scales", "metR",
              "doParallel", "scam", "mvtnorm", "MCMCpack", "mcmcplots")
load_inst(packages)

sourceCpp("rcpp/posterior_grad_xx2.cpp")

data("Cholesterol", package="qrLMM")

data <- Cholesterol

# rescale for better sampling
data$age <- rescale(Cholesterol$age, to = c(0, 1))
data$year <- rescale(Cholesterol$year, to = c(0, 1))

#standardize
data$cholst <- scale(data$cholst)


seed <- 123

##########################################################################################################################################################
# model 1 - varying coefficient for age -----------------------------------------------------------------------------------------------------------------------------
##########################################################################################################################################################
object_vcm <- bctm(cholst ~  hy_sm(cholst, data=data) +
              hyx_vcm(cholst, by=age, center=T, data=data, add_to_diag=10e-4) +
              hx_lin(age) + hx_lin(sex) + hx_lin(year),
               family = "gaussian", data=data,
               iterations = 2000,
               hyperparams=list(a=2, b=0.5), nuts_settings=list(adapt_delta = 0.95, max_treedepth=12), seed = seed)



# mcmcplots::traplot(object_vcm$samples$beta)
# information criteria (we refer to WAIC2 in the paper)
object_vcm$IC

# $DIC
# DIC       pD     Dbar     Dhat 
# 2764.529   12.640 2751.888 2739.248 
# 
# $WAIC1
# WAIC1      llpd    pwaic1 
# 2764.502 -1369.637    12.614 
# 
# $WAIC2
# WAIC2      llpd    pwaic2 
# 2765.079 -1369.637    12.902 
# 
# $LOOIC
# $LOOIC[[1]]
# Estimate        SE
# elpd_loo -1382.61528 25.706790
# p_loo       12.97818  1.098911
# looic     2765.23055 51.413579


##########################################################################################################################################################
# model 2 - full BCTM for age ----------------------------------------------------------------------------------------------------------------------------
##########################################################################################################################################################

object_te <- bctm(cholst ~  hyx_sm(cholst, age, data=data, q=c(10,10), add_to_diag=10e-4) +
                    hx_lin(sex) + hx_lin(year),
               family = "gaussian", data=data,
               iterations = 2000,
               hyperparams=list(a=2, b=0.5), nuts_settings=list(adapt_delta = 0.90, max_treedepth=12), seed = seed)

# mcmcplots::traplot(object_te$samples$beta)

# get information criteria
object_te$IC

# $DIC
# DIC       pD     Dbar     Dhat 
# 2734.190   12.855 2721.334 2708.479 
# 
# $WAIC1
# WAIC1      llpd    pwaic1 
# 2734.054 -1354.307    12.720 
# 
# $WAIC2
# WAIC2      llpd    pwaic2 
# 2734.524 -1354.307    12.955 
# 
# $LOOIC
# $LOOIC[[1]]
#            Estimate        SE
# elpd_loo -1367.31304 25.057703
# p_loo       13.00599  0.824181
# looic     2734.62608 50.115405
