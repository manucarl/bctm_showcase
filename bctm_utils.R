# Helper functions. To be integrated somewhere else.
## Author: Manuel Carlan


# installs packages in string if not installed yet and loads them
load_inst <- function(packages){
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  invisible(lapply(packages, library, character.only = TRUE))
  
}





# expand.grid but for dfs and not recursive
expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))




#  helper  for building a custom prediction environment
custom_fun_env <- function (f, ...) 
{
  environment(f) <- new.env(parent = globalenv())
  dots <- list(...)
  nm <- names(dots)
  sapply(nm, function(x) {
    assign(x, dots[[x]], envir = environment(f))
  })
  return(f)
}



# constructs individual smooths -------------------------------------------------------------------------------------------------------
construct_individual_designs <- function(y, x, knots, degree, q, derivs = c(0,0), center=T){

  # degree  of spline
  B1 <- splineDesign(knots[[1]], y, ord = degree[1] + 1, derivs=derivs[1], outer.ok = F )
  B2 <- splineDesign(knots[[2]], x, ord = degree[1] + 1, derivs=derivs[2], outer.ok = F)

  list(B1 = B1, B2 = B2)
}

