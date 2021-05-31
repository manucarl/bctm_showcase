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