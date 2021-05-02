##
## Script name: fram01_plots.R
##
## Purpose of script: estimates and plots the VCM model and the full BCTM for unstandardized y
##
## Author: Manuel Carlan
##
## Date Created: 2020-10-7
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
              "tidyverse", "profvis",  "tictoc", "scales", "metR",
              "doParallel", "scam", "mvtnorm", "MCMCpack", "mcmcplots")
load_inst(packages)





sourceCpp("rcpp/posterior_grad_xx2.cpp")

data("Cholesterol", package="qrLMM")

data <- Cholesterol

# scale for better sampling
data$age <- rescale(Cholesterol$age, to = c(0, 1))
data$year <- rescale(Cholesterol$year, to = c(0, 1))

seed <- 123

# function that scales prediction variable according to the scaling in training data
scale_pred <- function(x, var){
  (x - attr(var, "scaled:center"))/
    attr(var, "scaled:scale")
}


##########################################################################################################################################################
# model 1 - varying coefficient for age -----------------------------------------------------------------------------------------------------------------------------
##########################################################################################################################################################
object_vcm <- bctm(cholst ~  hy_sm(cholst, data=data) +
                     hyx_vcm(cholst, by=age, center=T, data=data, q=20, add_to_diag=10e-4) +
                     hx_lin(age) + hx_lin(sex) + hx_lin(year),
                   family = "gaussian", data=data,
                   iterations = 2000,
                   hyperparams=list(a=2, b=0.5), nuts_settings=list(adapt_delta = 0.9, max_treedepth=12), seed = seed)


# indices of parameters that are exp-transformed
exp_ident <-object_vcm$model$exp_ident

# beta samples
betas <- object_vcm$samples$beta

#beta_tilde samples
bts <- betas
bts[,exp_ident] <- exp(bts[,exp_ident])

# prediction grid for y
y_pred <- seq(min(data$cholst),max(data$cholst), length=100)

# prediction for year and age(rescaled to [0,1])
year_pred <- 0.4
age_pred <- rescale(seq(30,60, 5), c(0, 1))

# prediction df
pred_grid <-expand.grid(cholst=y_pred, age=age_pred, year = 0.4, sex=0)

# 
predvars1 <-object_vcm$predvars[["hy_sm(cholst)"]]
predvars2 <-object_vcm$predvars[["hyx_vcm(cholst, by=age)"]]

B1_pred <-predvars1$B(pred_grid)
Bp1_pred <-predvars1$Bp(pred_grid)
B2_pred <-predvars2$B(pred_grid)
Bp2_pred <-predvars2$Bp(pred_grid)

# prediction design matrix
X_pred <- cbind(1, B1_pred, B2_pred, pred_grid[,c("age", "year", "sex")]) %>% as.matrix
# derivative of prediction design matrix
Xp_pred <-  cbind(0, Bp1_pred, Bp2_pred, 0, 0, 0)%>% as.matrix


# posterior mean beta_tilde
bt <- object_vcm$beta_tilde

# estimated conditional transformation function
h_est <- X_pred%*%bt
# estimated derivative of conditional transformation function
hp_est <- Xp_pred%*%bt


# estimated conditional density
d_est <- dnorm(h_est)*hp_est

# plot data
pred_dat<- data.frame(cholst=pred_grid$cholst, d_est, h_est, pred_grid)


library(viridis)
p1 <- ggplot(pred_dat) +
  geom_line(aes(x=cholst, y=d_est, group=age),show.legend = F, size=0.5)+
  geom_ribbon(aes(x=cholst, ymin=0, ymax=d_est, fill=as.factor(age),colour =as.factor(age)),
              alpha=0.5, inherit.aes = FALSE)+
  ylab("Density")  + xlab("Cholesterol")+# scale_colour_discrete(name="Age")   +
  scale_color_manual(values=rep("black", 8, name="Age"), guide=F)+
  scale_fill_viridis(option = "A", discrete=T, name="Age", labels = seq(30,60, 5))+ #+ scale_color_viridis(option = "E", discrete=T, guide=F)
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "white", fill=NA),
        legend.box.background = element_rect(color = NA, colour="white"),
        legend.background = element_rect(color=NA),
        rect = element_rect(fill = "transparent"), 
        legend.justification = c(1, 1), legend.position = c(0.9, 0.9))
p1
ggsave("D:/GitHub/bctm_paper/manuscript/figs/fram_densities_vcm.pdf", plot =p1, height=4, width=8, units="in", bg="transparent")

##########################################################################################################################################################
# model 2 - full BCTM for age ----------------------------------------------------------------------------------------------------------------------------
##########################################################################################################################################################

object_te <- bctm(cholst ~  hyx_sm(cholst, age, data=data, q=c(10,10), add_to_diag=10e-4) +
                    hx_lin(sex) + hx_lin(year),
                  family = "gaussian", data=data,
                  iterations = 2000,
                  hyperparams=list(a=2, b=0.5), nuts_settings=list(adapt_delta = 0.90, max_treedepth=12), seed = seed)

# prediction grid for y
y_pred <- seq(min(data$cholst),max(data$cholst), length=100)

# age predictions (rescaled to [0,1] in next step)
age_pred <- seq(30, 60, length=20)
year_pred <- 0.4

# grid for construction of tensor spline
pred_grid <- expand.grid(cholst=y_pred, age = rescale(age_pred, to=c(0,1)))

predvars <-object_te$predvars[["hyx_sm(cholst, age)"]]

#  tensor prediction matrix
B_pred <-predvars$B(pred_grid)
Bp_pred <-predvars$Bp(pred_grid)

# design matrices for predictions
X_pred <- cbind(1, B_pred, 0, year_pred)
Xp_pred <- cbind(0, Bp_pred, 0, 0)

# vector of posterior means of reparametrized basis coefficients beta_tilde 
bt <- object_te$beta_tilde

# estimated conditional transformation function
h_est <- X_pred%*%bt

# estimated conditional density
d_est <- dnorm(X_pred%*%bt)*Xp_pred%*%bt

pred_dat<- cbind(d_est, h_est, pred_grid)


p2 <- ggplot(pred_dat, aes(x = cholst, y = age, z = d_est))+
  geom_contour_filled(binwidth=0.001)+
  geom_text_contour(aes(z = d_est))+
  theme_bw()+
  theme(legend.position="none") +
  scale_fill_viridis_d(option="A") +
  xlab("Cholesterol")+ylab("Age")  +
  scale_x_continuous( expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), n.breaks=6, labels=function(x) rescale(x, to=range(Cholesterol$age)))

p2

# ggsave("D:/GitHub/bctm_paper/manuscript/figs/fram_contours.pdf", plot =p2, height=4, width=8, units="in", bg="transparent")
