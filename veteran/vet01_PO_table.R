##
## Script name: vet01_PO_table.R
##
## Purpose of script: calculates contents of Table 1 (VA lung cancer) in "Bayesian Conditional Transformation Models"
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
              "tidyverse", "profvis",  "tictoc", "scales", "metR", "caret",
              "doParallel", "scam", "mvtnorm", "MCMCpack", "mcmcplots")
load_inst(packages)

sourceCpp("rcpp/posterior_grad_xx2.cpp")


data(veteran, package="survival")

data2 <- veteran[!as.logical(veteran$prior),]

dummies <- predict(dummyVars(~ celltype, data = data2), newdata = data2)
colnames(dummies) <- c("squam", "small", "adeno", "large")
data2 <- cbind(data2[,c("time", "status", "karno", "trt")], dummies)

# data2$karno <- rescale(data2$karno, c(0,1))
data2$trt <- ifelse(data2$trt == 2, 1, 0)

# expand.grid for
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

n <- nrow(data)


family <- "logistic"

seed <- 123

its <- 2000


# model 1 - po with censoring-------------------------------------------------------------------------------------------------------------------------
object <- bctm(time ~ hy_sm(time,data=data2, q=22, add_to_diag=0) + hx_lin(karno) + hx_lin(adeno) + hx_lin(small) + hx_lin(squam),  
            family = family, data=data2,
            cens=as.logical(data2$status),
            iterations = 2000,
            hyperparams=list(a=2, b=0.5), nuts_settings=list(adapt_delta = 0.90, max_treedepth=12), seed = seed)


# traplot(object$samples$beta)

burnin <- object$mcmc$burnin

# indices of exponentiated coefficients
exp_ident <- object$model$exp_ident


# beta samples
beta_samples <- object$samples$beta[(burnin+1):its,]

# beta_tilde samples
bt_samples <- beta_samples
bt_samples[,exp_ident] <- exp(beta_samples[,exp_ident])


# same model with MLT
 library(mlt)
 data3 <- data2
 data3$time <- with(data2, survival::Surv(time, status))
 var_t <- numeric_var("time", support = c(0, 500), bounds = c(0, 600))
 b_t <- Bernstein_basis(var_t, order = 10, ui = "increasing")
 b_R <- as.basis(~ karno  + adeno + small + squam, data = data3, remove_intercept = TRUE, negative=F)
 ctm_t <- ctm(response = b_t, shifting = b_R, todistr = "Logistic")
 mlt_t <- mlt(ctm_t, data = data3, scale = TRUE)

 
coefs <- tibble(Parameter = c("score", "adeno vs. large", "small vs. large", "squamous vs. large"),
   BCTM= round(apply(bt_samples, 2, median), 3)[-(1:22)],
                 MLT =round(coef(mlt_t), 3)[-(1:11)],
                 MPT = c(-0.055, 1.303, 1.362, -0.173) # values from paper
 )
 
sds <- tibble(BCTM_sd=round(apply(bt_samples, 2, sd), 3)[-(1:22)],
              MLT_sd= round(sqrt(as.vector(diag(vcov(mlt_t)))), 3)[-(1:11)],
              MPT_sd=   c(0.010, 0.559, 0.527, 0.580)
)


table1 <- bind_cols(coefs, sds, )[,c(1,2,5, 3,6, 4,7)]
table1
# Parameter            BCTM BCTM_sd    MLT MLT_sd    MPT MPT_sd
# <chr>               <dbl>   <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
#    1 score              -0.055   0.011 -0.057  0.011 -0.055  0.01 
# 2 adeno vs. large     1.37    0.55   1.36   0.561  1.30   0.559
# 3 small vs. large     1.47    0.51   1.46   0.533  1.36   0.527
# 4 squamous vs. large -0.147   0.605 -0.188  0.598 -0.173  0.580

#xtable::xtable(table1, digits=3)
 