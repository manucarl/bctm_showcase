##
## Script name: nuts.R
##
## Purpose of script: implements the No-U-Turn sampler by Hoffmann and Gelman (2014) with dual averaging and mass matrix adaption for 
##                    Bayesian Conditional transformation models; imitates the output style of
#                     rstan (Stan Development Team, 2020) and uses mass matrix adaptation from adnuts (Monnahan and Kristensen, 2018)
##
## Author: Manuel Carlan
##
## Date Created: 2020-10-7
##
## Email: mcarlan@uni-goettingen.de
##
## ---------------------------


NUTS <- function(n_iter, xx, f, gr, ll, start, warmup = floor(n_iter/2),thin=1,
                 seed = NULL,  chain = 1, nuts_settings, fixed = NULL,
                 prior_settings){
  

  its <- n_iter
  if(!is.null(seed)) set.seed(seed)
  if(chain != 1) stop("currently no parallel chains allowed")
  n_coef <- length(start)
  # supplement missing defaults
  if(!is.null(nuts_settings)){
    default_control <- list(adapt_delta = 0.8, metric = NULL, step_size = NULL, adapt_M = TRUE, max_treedepth=12, init_buffer = 75, term_buffer = 50, window = 25)
    nuts_settings <- modifyList(default_control, nuts_settings)
  }
  
  # number of coefficient  
  n_coef <- length(start)
  
  max_td <- nuts_settings$max_treedepth
  adapt_delta <- nuts_settings$adapt_delta
  adapt_M <- nuts_settings$adapt_M
  M <- nuts_settings$metric
  
  
  if(is.null(M)) M <- rep(1, n_coef)
  if(!is.vector(M)) stop("only diagonal mass matrices (with positive diagonal entries) allowed")
  
  
  ## BCTM uses stan default values for sampler warmup (https://mc-stan.org/docs/2_26/reference-manual/hmc-algorithm-parameters.html)
  interval1 <- nuts_settings$init_buffer 
  interval2 <- nuts_settings$term_buffer
  interval3 <- nuts_settings$window
  
  aws <- interval2 # adapt window size
  anw <- interval1 + interval2 # adapt next window
  
  if(warmup < (interval1+interval2+interval3) & adapt_M){
    warning("Specified warmup-up phase is too short for adaption.")
    adapt_M <- FALSE
  }
  
  rotation <- rotate_rescale(f = f,  gr = gr, ll = ll ,xx = xx, M=M, beta_y = start)
  

  n_coef <- xx[["n_coef"]]
  n_pen_grps <- xx[["npen"]]
  
  
  Ks_f <- xx[["Ks_f"]]
  
  hyper_a <- xx[["hyperparams"]][["a"]]
  hyper_b <- xx[["hyperparams"]][["b"]]
  
  ranks <- xx$ranks
  pen_ident <- xx$pen_ident
  
  # Smats <- vector("list", 2)
  # n_pen_grpss <- xx$npen
  # Ks_new <- vector("list", npen)
  
  step_size <- nuts_settings$step_size
  
  
  f2 <- rotation$f2 # posterior 
  gr2 <- rotation$gr2
  ll2 <- rotation$ll2
  
  # obtain rotated and rescaled parameters
  beta_cur <- rotation$beta_cur
  
  # square metric
  M_sq <- rotation$M_sq
  
  # sampler params as in rstan
  sampler_params <-  matrix(0, its, 6)            
  
  
  xx$Xpt <-  t(xx$Xp)

  
  ## how many steps were taken at each iteration, useful for tuning
  # eps_start <- 0.01
  
  # should dual averaging be applied (default TRUE)
  dual_averaging <- is.null(step_size)              
  
  
  if(dual_averaging){
    H_bar <- eps_out <-  eps_bar <- rep(NA, length = warmup+1)
    step_size <- eps_out[1] <- eps_bar[1] <- find_reasonable_epsilon(beta = beta_cur, f = f2, gr = gr2, xx = xx)
    mu <- log(10*step_size)
    
    H_bar[1] <- 0
    gamma <- 0.05
    t0 <- 10
    kappa <- 0.75
    
  } else {
    ## dummy values to return
    eps_out <- eps_bar <- H_bar <- NULL
  }
  
  # xx <- as.environment(xx)
  

  mbm <- microbenchmark::microbenchmark(
    "posterior" = {
      b <-  f(runif(n_coef), xx)
    },
    "gradient" = {
      b <- gr(runif(n_coef), xx)
    })
  print(mbm)
  j_out <-  rep(0, n_iter)
  
  
  # message('')
  # message(paste('Starting NUTS at', start))
  
  Ks_t <- Ks_f
  # names(Ks_t) <- names(pen_ident)
  
  K_inds <- xx$K_inds
  labels <- xx$labels
  npen <- xx$npen
  
  
  ## beta are the position parameters
  ## r are the momentum parameters
  beta_out <- matrix(0, nrow=n_iter, ncol=n_coef)
  colnames(beta_out) <- colnames(xx$X)
  
  log_liks <- log_posteriors <- rep(0, len=n_iter)
  tau2_out <- matrix(1, nrow=n_iter, ncol=n_pen_grps)
  colnames(tau2_out) <- unlist(  lapply(K_inds, function(x) paste0("tau2_", x)))
  # which(eff_pen[1] == name_groups)
  # lapply(  names(pen_ident))
  # inds[[1]] <- 
  pb <- progress::progress_bar$new(format = "[:bar] :current/:total (:percent)", total = n_iter)
  pb$tick(0)
  
  tau2 <- rep(0, npen)
  
  # inds <- 1:n_coef
  # Ks <- vector(mode="list", length=n_tau)
  
  start <- Sys.time()
  
  for(iter in 1:n_iter){
    
    # sourceCpp(paste0(sourcepath,"rcpp/gauss_hmc_update5.cpp"))
    beta_minus <- beta_plus <- beta_cur
    beta_out[iter,] <- M_sq*beta_cur 
    
    log_liks[iter] <- if(iter == 1) ll2(beta_cur, xx) else log_liks[iter-1]
    log_posteriors[iter] <- if(iter == 1) ll2(beta_cur, xx) else log_posteriors[iter-1]
    
    r_cur <- r_plus <- r_minus <-  rnorm(n_coef, 0, 1)
    H_prime <- calculate_hamiltonian(beta = beta_cur,  xx = xx, r = r_cur, f = f2)
    
    ## slicing step
    log_u <- log(runif(1)) + calculate_hamiltonian(beta=beta_cur, xx=xx,r=r_cur, f=f2)
    j <- 0
    n <- 1
    s <- 1
    div <- 0

    info <- as.environment(list(n.calls=0, div=0))
    
    while(s==1) {
      
      #print(j)
      
      # choose a direction v
      v <- sample(c(1,-1), 1)
      
      if(v==1){
        temp <- build_tree(beta = beta_plus,  xx = xx,r = r_plus, log_u = log_u, v = v,
                           j = j, step_size = step_size, H_prime = H_prime,
                           f = f2, gr =gr2, info = info)
        beta_plus <- temp$beta_plus
        r_plus <- temp$r_plus
      } else {
        
        temp <- build_tree(beta = beta_minus,  xx = xx,r = r_minus, log_u = log_u, v = v,
                           j = j, step_size = step_size, H_prime = H_prime,
                           f = f2, gr =gr2, info=info)
        beta_minus <- temp$beta_minus
        r_minus <- temp$r_minus
      }
      
      
      if(!is.finite(temp$s)) temp$s <- 0
      if(temp$s==1) {
        if(runif(1) <= temp$n/n){
          beta_cur <- temp$beta_prime
          log_liks[iter] <- ll2(beta_cur, xx)
          log_posteriors[iter] <- f2(beta_cur, xx)
          
          # beta_out is on the rotated scale
          beta_out[iter,] <- M_sq*beta_cur
        }
      }
      # end if
      n <- n + temp$n
      s <- temp$s*check_nuts(beta_plus, beta_minus, r_plus, r_minus)
      
      #print(paste0("s: ", s))
      j <- j+1
      
      if(!is.finite(s)) s <- 0
      
      if(j >= max_td) {
        warning(paste0("max_treedepth(", max_td,  ") reached"))
        break
      }
    }
    
    j_out[iter] <- j-1
    alpha2 <- temp$alpha/temp$n_alpha
    if(!is.finite(alpha2)) alpha2 <- 0
    
    
    if(dual_averaging){
      if(iter <= warmup){
        
        H_bar[iter+1] <- (1-1/(iter + t0))*H_bar[iter] + (adapt_delta - alpha2)/(iter + t0)
        # set logeps and logepsbar
        logeps <- mu - sqrt(iter)*H_bar[iter+1]/gamma
        eps_out[iter+1] <- exp(logeps)
        
        logepsbar <- iter^(-kappa)*logeps + (1-iter^(-kappa))*log(eps_bar[iter])
        eps_bar[iter+1] <- exp(logepsbar)
        
        step_size <- eps_out[iter+1]
      } else {
        step_size <- eps_bar[warmup]
        # step_size <- eps_out[warmup]
        if(iter== warmup +1) print(paste("step size after warmup: ",step_size))
        
      }
    }
    
    
    # this is from adnuts ----------------------------------------------
    if(adapt_M & .slow_phase(iter, warmup, interval1, interval3)){
      
      if(iter== interval1){
        
        m1 <- beta_out[iter,]
        s1 <- rep(0, len=n_coef)
        k <- 1
      } else if(iter==anw){
        
        M <- as.numeric(s1/(k-1)) # estimated variance
        
        rotation <- rotate_rescale(f=f,  xx=xx,gr=gr, ll=ll, M=M,  beta_y=beta_out[iter,])
        f2 <- rotation$f2
        gr2 <- rotation$gr2
        ll2 <- rotation$ll2
        
        M_sq <- rotation$M_sq
        
        beta_cur <- rotation$beta_cur 
        ## Reset the running variance calculation
        k <- 1
        s1 <- rep(0, n_coef)
        
        m1 <- beta_out[iter,]
        
        aws <- 2*aws
        anw <- .compute_next_window(iter, anw, warmup, interval1, aws, interval3)
        
        step_size <- find_reasonable_epsilon(beta = beta_cur,  xx = xx,f = f2,  gr = gr2)
        
      } else {
        k <- k+1; m0 <- m1; s0 <- s1
        m1 <- m0+(beta_out[iter,]-m0)/k
        s1 <- s0+(beta_out[iter,]-m0)*(beta_out[iter,]-m1)
      }
    }
    #---------------------------------------------------------------------------
    

    # update tau2s and multiply with corresponding precision matrix
    for(i in 1:npen){
      grp <- pen_ident[[i]] 
      par <- beta_out[iter, grp ]
      
      tau2_out[iter,i] <-  tau2 <- rinvgamma(1, hyper_a + 0.5*ranks[[i]], hyper_b + as.vector(0.5*t(par)%*%(Ks_f[[i]]%*%par)))
      
      Ks_t[[i]] <- Ks_f[[i]]/ tau2
      
    }
    
    
    
    S <- as.matrix(bdiag(lapply(K_inds, function(x) Reduce("+", Ks_t[x]))))
    xx$S <- S
    
    
    sampler_params[iter,] <-  c(alpha2, step_size, j, info$n.calls, info$div, f2(beta_cur, xx))
    
    pb$tick(1)
    if(iter==warmup) time.warmup <- difftime(Sys.time(), start, units='secs')
    .print.mcmc.progress(iter, n_iter, warmup, chain)
  } 
  
  beta_out <- beta_out[seq(1, nrow(beta_out), by=thin),]
  warmup <- warmup/thin
  colnames(sampler_params) <- c(NULL, "accept_stat", "step_size", "treedepth", "n_leapfrog", "divergent", "energy")
  sampler_params <- as_tibble(sampler_params)
  
  
  # this is from adnuts/Rstan--------------------------------------------------------------------------------------------
  sampler_params <- sampler_params[seq(1, nrow(sampler_params), by=thin),]
  ndiv <- sum(sampler_params[-(1:warmup),5])
  
  if(ndiv>0)    message(paste0("There were ", ndiv, " divergent transitions after warmup"))
  
  msg <- paste0("Final acceptance ratio=", sprintf("%.2f", colMeans(sampler_params[-(1:warmup),1])))
  
  if(dual_averaging) msg <- paste0(msg,", and target=", adapt_delta)
  message(msg)
  if(dual_averaging) message(paste0("Final step size=", round(step_size, 3),
                                    "; after ", warmup, " warmup iterations"))
  time.total <- difftime(Sys.time(), start, units='secs')
  .print.mcmc.timing(time.warmup=time.warmup, time.total=time.total)
  #-----------------------------------------------------------------------------------------------------------------
  
  list(beta = beta_out, tau2 = tau2_out, log_liks = log_liks, lp = log_posteriors, sampler_params = sampler_params,
       time.total = time.total, time.warmup = time.warmup,
       warmup = warmup, max_treedepth = max_td)
}

# log likelihoods (only used for calculation of ICs)
ll_gauss <- function(param, xx){
  bt <- param
  exp_ident <- xx$exp_ident+1
  bt[exp_ident] <- exp(bt[exp_ident])
  sum(dnorm(xx$X%*%bt, log=T)) + sum(log(xx$Xp%*%bt))
}


ll_logit <- function(param, xx){
  bt <- param
  exp_ident <- xx$exp_ident+1
  bt[exp_ident] <- exp(bt[exp_ident])
  sum(dlogis(xx$X%*%bt, log=T)) + sum(log(xx$Xp%*%bt))
}

ll_mev <- function(param, xx){
  bt <- param
  exp_ident <- xx$exp_ident+1
  bt[exp_ident] <- exp(bt[exp_ident])
  sum(bctm_dmev(xx$X%*%bt, log=T)) + sum(log(xx$Xp%*%bt))
}

# minimum-extreme-value density
bctm_dmev <- function(x, log = F) {
  
  ret <- x - exp(x)
  if (!log) return(exp(ret))
  ret
}

# minimum-extreme-value cdf

bctm_pmev <- function(x) {
  
  1 - exp(-exp(x))
  
}



calculate_hamiltonian <- function(beta, xx, r, f) f(beta, xx)- 0.5*sum(r^2)


check_nuts = function(beta_plus, beta_minus, r_plus, r_minus){
  (crossprod(beta_plus - beta_minus, r_minus)) >= 0 && (crossprod(beta_plus - beta_minus, r_plus) >= 0)
}





rotate_rescale <- function(f, gr, ll, xx, M,  beta_y){
  
  if(!is.vector(M)) stop("only diagonal mass matrices (with positive diagonal entries) allowed")
  M_sq <- sqrt(M)
  
  # apply opposite transformation to the parameters
  ll2 <- function(beta, xx) ll(M_sq * beta, xx)
  f2 <- function(beta, xx) f(M_sq * beta, xx)
  gr2 <- function(beta, xx) as.vector(gr(M_sq * beta, xx) ) * M_sq
  
  # apply transformation (rotation and scaling) to simplify kinetic energy
  beta_cur <- (1/M_sq) * beta_y
  
  list(f2 = f2, gr2 = gr2, ll2 = ll2, beta_cur = beta_cur, M_sq = M_sq)
}


# this is the FindReasonableEpsilon() in Hoffman and Gelman (2014)
find_reasonable_epsilon <- function(beta, f, gr, xx){
  
  #initialize
  step_size <- 1
  r <- rnorm(length(beta))
  
  ##leapfrog step
  lf <- leapfrog_step(beta, r, step_size, gr, xx)
  
  beta_prime <- lf$beta_prime
  r_prime <- lf$r_prime
  
  H <- calculate_hamiltonian(beta = beta,  xx=xx, r=r, f=f)
  H_prime <- calculate_hamiltonian(beta = beta_prime,  xx=xx, r=r_prime, f=f)
  a <- 2*(exp(H_prime)/exp(H) > 0.5)-1
  
  if(!is.finite(a)) a <- -1
  k <- 1
  
  while (!is.finite(H) | !is.finite(H_prime) | a*H_prime-a*H > -a*log(2)) {
    step_size <- (2^a)*step_size
    
    lf <- leapfrog_step(beta, r, step_size, gr, xx)
    beta_prime <- lf$beta_prime
    r_prime <- lf$r_prime
    H_prime <- calculate_hamiltonian(beta = beta_prime,  xx = xx, r = r_prime, f = f)
    k <- k + 1
    if(k > 500) stop("Could not find reasonable epsilon in 500 iterations")
  }
  # if(verbose) message(paste("Reasonable epsilon=", step_size, "found after", k, "steps"))
  return(invisible(step_size))
}


# generic leapfrog update
leapfrog_step = function(beta, r, step_size, gr, xx){
  r_prime <- r + 0.5 * step_size * gr(beta, xx)
  beta_prime <- beta + step_size * r_prime
  r_prime <- r_prime + 0.5 * step_size * gr(beta_prime, xx)
  
  list(beta_prime = beta_prime, r_prime = r_prime)
}





build_tree <- function(beta, xx, r, log_u, v, j, step_size, H_prime, f, gr,
                      delta_max=1000, info = environment() ){
  
  # print(diag(S))
  if(j==0){
    
    
    lf <- leapfrog_step(beta, r, v*step_size, gr, xx)
    beta <- lf$beta_prime
    r <- lf$r_prime
    
    
    H <- calculate_hamiltonian(beta=beta,  xx=xx, r=r, f=f)
    n <- log_u <= H
    s <- log_u < delta_max + H
    if(!is.finite(H) | s == 0){
      info$div <- 1
      s <- 0
    }
    
    log_alpha <- H-H_prime
    alpha <- min(exp(log_alpha), 1)

    
    return(list(beta_minus=beta, beta_plus=beta, beta_prime=beta, r_minus=r,
                r_plus=r,  s=s, n=n, alpha=alpha, n_alpha=1))
  } else {
    ## recursion - build left and right branches
    branch1 <- build_tree(beta=beta,  xx=xx, r=r, log_u=log_u, v=v, j=j-1, step_size=step_size,
                          H_prime=H_prime, f=f, gr=gr, info=info)
    
    beta_minus <- branch1$beta_minus
    beta_plus <- branch1$beta_plus
    beta_prime <- branch1$beta_prime
    r_minus <- branch1$r_minus
    r_plus <- branch1$r_plus
    alpha <- branch1$alpha
    n_alpha <- branch1$n_alpha
    s <- branch1$s
    if(!is.finite(s)) s <- 0
    nprime <- branch1$n
    
    if(s==1){
      if(v== -1){
        branch2 <- build_tree(beta = beta_minus,  xx = xx, r = r_minus, log_u = log_u, v=v, j=j-1, step_size=step_size, H_prime = H_prime, f=f, gr=gr, info=info)
        beta_minus <- branch2$beta_minus
        r_minus <- branch2$r_minus
      } else {
        branch2 <- build_tree(beta=beta_plus,  xx=xx, r =r_plus, log_u = log_u, v=v, j=j-1, step_size=step_size, H_prime = H_prime, f=f, gr=gr, info=info)
        beta_plus <- branch2$beta_plus
        r_plus <- branch2$r_plus
      }
      
      nprime <- branch2$n+ branch1$n
      if(!is.finite(nprime)) nprime <- 0
      
      ## acceptance step
      if(nprime>0){
        if(runif(1) <= branch2$n/nprime){
          beta_prime <- branch2$beta_prime
      alpha <- branch1$alpha+branch2$alpha
      n_alpha <- branch1$n_alpha+branch2$n_alpha
        }
      }
      # check if proposal is valid
      b <- check_nuts(beta_plus = beta_plus, beta_minus=beta_minus, 
                      r_plus = r_plus,  r_minus =  r_minus)
      s <- branch2$s*b
    }
    return(list(beta_minus =beta_minus, beta_plus = beta_plus, beta_prime=beta_prime,
                r_minus = r_minus, r_plus = r_plus,
                s = s, n = nprime,
                alpha = alpha, n_alpha=n_alpha))
  }
}








