##
## Script name: vet03_NPO_trt.R
##
## Purpose of script: caclulates DIC/WAIC VA PO model with simple shift for treatment and NPO model with interaction for treatment 
##
## Author: Manuel Carlan
##
## Date Created: 2020-11-12
##
## Email: mcarlan@uni-goettingen.de
##
## ---------------------------

library(dplyr)
# rm(list = ls())
path <- "D:/GitHub/bctm_paper/code/"
setwd(path)

source("bctm_utils.R")
source("bctm_design_funs2.R")
source("bctm_design.R")
source("bctm_fun.R")

source("nuts/nuts_utils.R")
source("nuts/nuts.R")
source("nuts/adnuts_helper.R")

packages <- c("Rcpp", "RcppArmadillo", "RcppEigen", "splines", "mgcv", "Matrix", "MCMCpack", 
              "tidyverse", "profvis",  "tictoc", "scales", "metR", "caret",
              "doParallel", "scam", "mvtnorm", "MCMCpack", "mcmcplots")
load_inst(packages)

sourceCpp("rcpp/posterior_grad_xx2.cpp")

# get data
data(veteran, package="survival")
data <- veteran[!as.logical(veteran$prior),]
n <- nrow(data)

# create dummies
dummies <- predict(dummyVars(~ celltype, data = data), newdata = data)
colnames(dummies) <- c("squam", "small", "adeno", "large")
data <- cbind(data[,c("time", "status", "karno", "trt")], dummies)

# rescale karnofsky score to [0,1] for faster sampling
data$karno <- rescale(data$karno, c(0,1))
# dummy coding
data$trt <- ifelse(data$trt == 2, 1, 0)


its <- 2000

# logistic family for proportional odds

seed <- 123

# model 1 - po with treatment shift----------------------------------------------------------------------------------------------------------------
object_po <- bctm(time ~ hy_sm(time,data=data) + 
                    hx_lin(trt),  
                  family = "logistic", data=data, iterations = 2000, cens=as.logical(data$status),
                  hyperparams=list(a=2, b=0.5), nuts_settings=list(adapt_delta = 0.90, max_treedepth=12), seed = seed)


# traplot(object_po$samples$beta)

object_po$IC
# $DIC
# DIC      pD    Dbar    Dhat 
# 191.034  10.174 180.860 170.686 
# 
# $WAIC1
# WAIC1    llpd  pwaic1 
# 190.231 -86.218   8.898 
# 
# $WAIC2
# WAIC2    llpd  pwaic2 
# 190.995 -86.218   9.279 
# 
# $LOOIC
# $LOOIC[[1]]
# Estimate       SE
# elpd_loo -95.576588 10.75271
# p_loo      9.358656  1.22139
# looic    191.153176 21.50542


object_npo <- bctm(time ~ hy_sm(time,data=data) +
                     hyx_vcm(time, by= trt, data=data, center=T, add_to_diag=10e-4),  
                   family = "logistic", data=data, cens=as.logical(data$status), iterations = 2000,
                   hyperparams=list(a=2, b=0.5), nuts_settings=list(adapt_delta = 0.9, max_treedepth=12), seed = seed)

object_npo$IC
# $DIC
# DIC      pD    Dbar    Dhat 
# 192.185   9.466 182.719 173.253 
# 
# $WAIC1
# WAIC1    llpd  pwaic1 
# 191.851 -87.533   8.393 
# 
# $WAIC2
# WAIC2    llpd  pwaic2 
# 192.577 -87.533   8.756 
# 
# $LOOIC
# $LOOIC[[1]]
# Estimate        SE
# elpd_loo -96.378257 10.499887
# p_loo      8.845715  1.186272
# looic    192.756514 20.999774