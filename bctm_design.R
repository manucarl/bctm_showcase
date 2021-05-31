bctm_setup <- function (formula, data, intercept, specials = c("hx_lin", "hyx_sm,", "hyx_vcm", "hy_sm", "hx_sm", "hx_spat")){
  
  
  
  trms <- terms.formula(formula, specials = specials)
  trafos <- trms[attr(trms, "order") == 1]
  labels <- attr(trms, "term.labels")
  allvars <- all.vars(formula)
  
  trafo_designs <- sapply(1:sum(attr(trms, "order") == 
                                  1), function(i) {
                                    with(data, eval(trafos[i][[3]]))
                                  }, simplify = FALSE)
  
  
  label_inds <- lapply(trafo_designs, "[[", "label_inds")
  # names(trafo_designs) <- labels[attr(trms, "order") ==  1]
  labels <- unlist(lapply(trafo_designs, "[[", "label"))
  names(trafo_designs) <-labels
  
  Bs <- lapply(trafo_designs, "[[", "B")
  Bps <- lapply(trafo_designs, "[[", "Bp")
  # names(Bs) <- names(Bps) <-  first_labels
  
  
  interactions <- labels[attr(trms, "order") != 1]
  # construct_interactions <- function(term, Bs) {
  #   
  #   terms <- strsplit(interactions, ":")[[1]]
  #   marginals <- Bs[terms]
  #   marginals_p <- Bps[terms]
  #   
  #   label <- paste(sapply(Xs, attr, which = "label"), 
  #                  collapse = ":")
  #   
  #   B <- marginals[[1]]%.%marginals[[2]]
  #   
  #   Bp <- marginals_p[[1]]%.%marginals[[2]]
  #   
  #   colnames(B) <- colnames(Bp) <- namesB <- paste(rep(colnames(Xs[[1]]), 
  #                                                      e = NCOL(Xs[[2]])), rep(colnames(Xs[[2]]), t = NCOL(Xs[[1]])), 
  #                                                  sep = ":")
  #   
  #   list(B = B, Bp = Bp, predvars )
  #   
  # }
  
  # interactionDesigns <- sapply(interactions, construct_interactions, trafo_designs)
  
  y <- with(data, eval(attr(trms, "variables")))[[attr(trms, "response")]]
  n <- length(y)
  
  
  exp_inds <- lapply(trafo_designs, "[[", "exp_ident")
  
  pen_inds <-  unlist(lapply(trafo_designs, "[[", "pen_ind"))
  
  Ks <- lapply(trafo_designs, "[[", "K")
  
  
  
  pen_groups <- rep(unlist(labels[pen_inds]), times = sapply(Ks, "length")[pen_inds])
  
  
  predvars <- lapply(trafo_designs, "[[", "predvars")
  
  
  # if (attr(trms, "intercept") == 1) {
  if (intercept) {
    
    X <- matrix(1, nrow = n, ncol = 1)
    Xp <- matrix(0, nrow = n, ncol = 1)
    label_inds <- c("int", label_inds)
    labels <- c("int", labels)
    pen_inds <- c("int" = FALSE, pen_inds)
    colnames(X) <- "int"
    eff_inds <- c("int")
    Kint <- 1
  }  else {
    Kint <- NULL
    X <- Xp <- matrix(NA, nrow = n, ncol = 0)
    eff_inds <- c()
  }
  
  
  X <- cbind(X, do.call(cbind, Bs))
  Xp <- cbind(Xp, do.call(cbind, Bps))
  
  eff_inds <- unlist(c(eff_inds, rep(names(Bs), sapply(Bs, NCOL))))
  names(eff_inds) <- eff_inds
  
  eff_ident <- lapply(labels, function(x) unname(which(eff_inds == x)))
  names(eff_ident) <- labels
  
  
  pen_ident <- eff_ident[pen_groups]
  unpen_ident <- eff_ident[!pen_inds]
  
  response <- all.vars(formula)[1]
  y <- data[,response]
  attr(y, "response") <-  response
  
  
  ret <- list(y = y, X = X, Xp = Xp, Ks = c("int" = Kint, Ks), label_inds = label_inds, labels = labels,
              exp_inds = exp_inds, pen_ident = pen_ident, unpen_ident = unpen_ident, pen_inds = pen_inds,
              eff_ident = eff_ident, trafo_designs = trafo_designs, formula = formula, predvar = predvars, specials = specials
  )
  ret
  
}
