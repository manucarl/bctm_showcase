## ---------------------------
##
## Script name: leuk01_plots.R
##
## Purpose of script: creates survivor function and spatial plot for leukemia data in "Bayesian Conditional Transformation Models"
##
## Author: Manuel Carlan
##
## Date Created: 2020-6-5
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

packages <- c("spBayesSurv", "survival", "BayesX", "scam", "Matrix", "Rcpp", "RcppArmadillo", "MCMCpack", "sf", "rgeos", "ggplot2")
load_inst(packages)


sourceCpp("rcpp/posterior_grad_xx2.cpp")

# get data
data("LeukSurv", package="spBayesSurv")
data <- LeukSurv[order(LeukSurv$district), ]
n <- nrow(data) 

#get boundary file
nwengland <- read.bnd(system.file("otherdata/nwengland.bnd", package = "spBayesSurv"))

# construct map in graph format for neighbour matrix
nmat <- bnd2gra(nwengland)

# function that scales prediction variable according to the scaling in training data
scale_pred <- function(x, var){
  (x - attr(var, "scaled:center"))/
    attr(var, "scaled:scale")
}

# scale for faster sampling
data$age <- scale(data$age)
data$wbc <- scale(data$wbc)

# minimum-extreme-value distribution for proportional hazards
family <- "mev"


seed <- 123
its <- 2000
# proportional hazards model with spatial effect -----------------------------------------------------------------------------------------------
object_ph <- bctm(time ~ hy_sm(time,data=data,  center=T)+
                    hx_spat(district, data=data, nmat=nmat)+
                    hx_lin(age) +hx_lin(sex) + hx_lin(wbc) + hx_lin(tpi) , #cens = as.logical(data$cens),
                  family = family, data=data, iterations = its, intercept=F, # remove intercept 
                  hyperparams=list(a=2, b=0.5), nuts_settings=list(adapt_delta = 0.80, max_treedepth=10), seed = seed)

# mcmcplots::traplot(object_ph$samples$beta)
object <- object_ph

# get model design matrix
# X <- object$X
# Xp <- object$Xp

burnin <- object$mcmc$burnin

# beta samples
betas <- object$samples$beta[burnin:its,]
colnames(betas) <- colnames(X)

# traplot(betas)

# indicators for exponentiation
exp_ident <- object$model$exp_ident


# beta_tilde samples
bt_samples <-  betas
bt_samples[,exp_ident] <- exp(bt_samples[,exp_ident])

# posterior mean
beta_tilde <- colMeans(bt_samples)


nwenglandsp <- bnd2sp(nwengland)
newenglandsf <-st_as_sf(nwenglandsp)

label_inds <-object$model$label_inds


# get beta_tilde vector of reparametrized basis coefficients
bt <- object$beta_tilde
xlabels <- c("hx_lin(tpi)", "hx_lin(age)", "hx_lin(sex)", "hx_lin(wbc)")

# create Table2 ---------------------------------------------------------------------------------------------------------------
bt[label_inds%in%xlabels]

beta_samples <- object$samples$beta[(burnin+1):its,]
bt_samples <- beta_samples
bt_samples[,exp_ident] <- exp(beta_samples[,exp_ident])

xsamples <- bt_samples[,label_inds%in%xlabels]
table2 <- cbind("Mean" = colMeans(xsamples), "s.d." = apply(xsamples, 2, sd), q=t(apply(xsamples , 2, quantile, c(0.5, 0.025, 0.975))))
rownames(table2) <- c("tpi", "age", "sex", "wbc")
table2

# Mean       s.d.        50%         2.5%      97.5%
# tpi 0.48866545 0.03420332 0.48796338  0.416508340 0.55748311
# age 0.08145468 0.06430202 0.08154937 -0.039386159 0.19997948
# sex 0.16455972 0.03459662 0.16480263  0.093313158 0.23242449
# wbc 0.02312324 0.00937159 0.02331872  0.004997466 0.04218026

# xtable::xtable(table2, digits=3)



# create plots -----------------------------------------------------------------------------------------------------------------
time_grid <- seq(0, 3000, length=300)

# adapt pred variables that are standardized in training data 
age_pred <-scale_pred(49, data$age)
wbc_pred <- scale_pred(38.6, data$wbc)

# construct design matrix for baseline transformation
pred_B <- lapply(object$predvars, "[[", "B")
B0pred <- pred_B[[1]](list(time=time_grid))


# get samples of transformation function for each quantile = 0.05, 0.5, 0.95 of tpi
effsamples095 <- cbind(
  h0=B0pred,
  age = age_pred,
  sex = 0,
  wbc = wbc_pred,
  tpi = quantile(LeukSurv$tpi, 0.99)
) %*% t(bt_samples[, which(label_inds != "hs(district)")]) 

effsamples050 <- 
  cbind(
    B0pred,
    age = age_pred,
    sex = 0,
    wbc = wbc_pred,
    tpi = quantile(LeukSurv$tpi, 0.5)
  ) %*% t(bt_samples[, which(label_inds != "hs(district)")])


effsamples005 <-
  cbind(
    h0=B0pred,
    age = age_pred,
    sex = 0,
    wbc = wbc_pred,
    tpi = quantile(LeukSurv$tpi, 0.01)
  ) %*% t(bt_samples[, which(label_inds != "hs(district)")]) 

effsamples005[1,] <- effsamples050[1,] <- effsamples095[1,] <- -Inf
# get data for plot A
gg <-
  lapply(list(
    "5%" = effsamples005,
    "50%" = effsamples050,
    "95%" = effsamples095
  ), function(x) {
    tibble(
      "pq975" = 1 - pmev(apply(x, 1, quantile, 0.975)),
      "pqmed" = 1 - pmev(apply(x, 1,  median)),
      "pq025" = 1 - pmev(apply(x, 1, quantile, 0.025)),
      time_grid
    )
  }) %>%
  bind_rows(.id = "q")

# panel A: survivor functions
p1 <- gg %>%
  ggplot(aes(x = time_grid, y = pqmed, col = q)) +
  geom_line() +
  theme_bw() +
  geom_ribbon(aes(ymin = pq025, ymax = pq975),
              alpha = 0.1,
              linetype = 0) +
  ylab("Survival") + xlab("Time") +
  theme(
    legend.justification = c(1, 1),
    legend.position = c(0.95, 0.95),
    panel.background = element_rect(fill = "transparent"),
    plot.background = element_rect(fill = "transparent", color = NA),
    legend.background = element_rect(fill = "transparent"),
    legend.box.background = element_rect(fill = "transparent")
  ) +
  scale_color_viridis_d(option = "C", begin=0.9, end=0.1) +
  labs(color = "Townsend score (quantiles)")

p1


# panel B: map ----------------------------------------------
library(viridis)
p2 <- ggplot(data = newenglandsf) +
  theme_void() +
  geom_sf(aes(fill = bt[label_inds == "hs(district)"])) +
  scale_fill_viridis(option = "A",
                     begin = 1,
                     end = 0.3) +
  theme(
    legend.justification = c(1, 1),
    legend.position = c(0.99, 0.99),
    panel.border = element_rect(colour = "white", fill = NA),
    legend.box.background = element_rect(color = NA, colour = "white"),
    legend.background = element_rect(color = NA),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent"),
    rect = element_rect(fill = "transparent") # all rectangles
  ) +
  labs(fill = "") #+ panel_border(remove = TRUE, color="white")
p2

p <-cowplot::plot_grid(p1, p2, labels=c("A","B"), ncol=2)#+ panel_border(remove = TRUE, color="white")
p
# ggsave("D:/GitHub/bctm_paper/manuscript/figs/leuk_ph.pdf", plot=p,height=3.5, width=9, units="in", bg = "transparent")


