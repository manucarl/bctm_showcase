## Author: Manuel Carlan


test_nuts <- function(beta_plus, beta_minus, r_plus, r_minus){
  beta_temp <- beta_plus - beta_minus
  temp <- (crossprod(beta_temp,r_minus) >= 0) *
    (crossprod(beta_temp, r_plus) >= 0)
  return(temp)
}

# implements Euclidian metric (qncour (2017)) by rotating and rescaling the target space, i.e. 
rotate_rescale <- function(f, gr, xx, M,  beta_y){
  
  if(!is.vector(M)) stop("only diagonal mass matrices (with positive diagonal entries) allowed")
  M_sq <- sqrt(M)
  
  # apply opposite transformation to the parameters
  f2 <- function(beta, xx) f(M_sq * beta, xx)
  gr2 <- function(beta, xx) as.vector(gr(M_sq * beta, xx) ) * M_sq
  
  # apply transformation (rotation and scaling) to simplify kinetic energy
  beta_x <- (1/M_sq) * beta_y
  
  list(f2 = f2, gr2 = gr2, beta_x = beta_x, M_sq = M_sq)
}


# this ist the FindReasonableEpsilon function in Hoffman et al(2014)
find_reasonable_epsilon <- function(beta, f, gr, xx){
  
  step_size <- 1
  r <- rnorm(n = length(beta))
  ## Do one leapfrog step
  r_star <- r + (step_size/2)*gr(beta, xx)
  beta_star <- beta + step_size * r_star
  r_star <- r_star + (step_size/2)*gr(beta_star,  xx)
  
  H1 <- calculate_hamiltonian(beta = beta,  xx=xx, r=r, f=f)
  H2 <- calculate_hamiltonian(beta = beta_star,  xx=xx, r=r_star, f=f)
  a <- 2*(exp(H2)/exp(H1)>.5)-1
  ## If jumped into bad region, a can be NaN so setup algorithm to keep
  ## halving step_size instead of throwing error
  if(!is.finite(a)) a <- -1
  k <- 1
  ## Similarly, keep going if there are infinite values
  while (!is.finite(H1) | !is.finite(H2) | a*H2-a*H1 > -a*log(2)) {
    step_size <- (2^a)*step_size
    ## Do one leapfrog step
    r_star <- r+(step_size/2)*gr(beta, xx)
    beta_star <- beta+step_size*r_star
    r_star <- r_star+(step_size/2)*gr(beta_star, xx)
    H2 <- .calculate.H(beta=beta_star,  xx=xx, r=r_star, f=f)
    k <- k+1
    if(k>500) {
      stop("Problem.")
    }
  }
  # if(verbose) message(paste("Reasonable epsilon=", step_size, "found after", k, "steps"))
  return(invisible(step_size))
}

calculate_hamiltonian <- function(beta, xx, r, f) f(beta, xx) - (1/2)*sum(r^2)



# helper functions for mass matrix adaption and mcmc output from package adnuts
.compute_next_window <- function(i, anw, warmup, w1, aws, w3){
  anw <- i+aws
  if(anw== (warmup-w3) ) return(anw)
  ## Check that the next anw is not too long. This will be the anw for the
  ## next time this is computed. If the next one is too long, extend this
  ## one to the very end.
  nwb <- anw+2*aws
  if(nwb >= warmup-w3){
    ## if(i != warmup-w3)
    ##   message(paste("Extending last slow window from", anw, "to", warmup-w3))
    anw <- warmup-w3
  }
  return(anw)
}




.slow_phase <- function(i, warmup, w1, w3){
  ## After w1, before start of w3
  x1 <- i>= w1 # after initial fast window
  x2 <- i<= (warmup-w3) # but before last fast window
  x3 <- i < warmup # definitely not during sampling
  return(x1 & x2 & x3)
}



.print.mcmc.progress <- function(iteration, iter, warmup, chain){
  i <- iteration
  refresh <- max(10, floor(iter/10))
  if(i==1 | i==iter | i %% refresh ==0){
    i.width <- formatC(i, width=nchar(iter))
    out <- paste0('Chain ',chain,', Iteration: ', i.width , "/", iter, " [",
                  formatC(floor(100*(i/iter)), width=3), "%]",
                  ifelse(i <= warmup, " (Warmup)", " (Sampling)"))
    message(out)
  }
}


.print.mcmc.timing <- function(time.warmup, time.total){
  x <- ' Elapsed Time: '
  message(paste0(x, sprintf("%.1f", time.warmup), ' seconds (Warmup)'))
  message(paste0(x, sprintf("%.1f", time.total-time.warmup), ' seconds (Sampling)'))
  message(paste0(x, sprintf("%.1f", time.total), ' seconds (Total)'))
}