## Implementation of different (conditional) transformation function types
## Author: Manuel Carlan
## ---------------------------


# nonlinear interactions that are monotone in direction of the first variable ---------------------------------------------------------------------------
hyx_sm <- function(..., q = c(10,10) , data, knots = list(NULL, NULL), 
                   pen_order = c(2,2), degree = c(3, 3), center = T, add_to_diag = 10e-6 )
{
  
  if (!((pen_order[1] == 2) & pen_order[2] ==2))  stop("right now only default penalties possible")
  if (!(degree[1] == 3 & degree[2] ==3))  stop("right now only cubic splines possible")
  if (!center)  stop("right now only center == TRUE allowed")
  
  sspec <- s(..., k=q, bs="tesmi1")

  sm <- smoothCon(sspec,data=data,
                    scale.penalty=F, absorb.cons=center, knots=NULL)[[1]]
  
  q1 <- q[1]
  q2 <- q[2]
  
  B01 <- sm$X
  
  I <- diag(q2)
  
  K <- vector("list", 2)
  P <- diff(diag(q1 - 1), difference = 1)
  Pm1 <- matrix(0, q1 - 1, q1)
  Pm1[2:(q1 - 1), 2:q1] <- P
  K[[1]] <- Pm1 %x% I
  K[[1]] <- t(K[[1]])%*%K[[1]]

  
  I2 <- diff(diag(q2), difference = 1)
  I21 <- diff(diag(q2), difference = 2)
  I1 <- diag(q1)
  K[[2]] <- matrix(0, q2 - 2 + (q1 - 1) * (q2 - 1), q1 * q2)
  K[[2]][1:(q2 - 2), ] <- t(I1[1, ]) %x% I21
  K[[2]][(q2 - 1):nrow(K[[2]]), ] <- I1[2:q1, ] %x% I2
  K[[2]] <- t(K[[2]])%*%K[[2]]
  
  Zc <- diag(q1 * q2)
  Zc <- Zc[, -q2]
  D1 <- t(diff(diag(q2)))
  Zc[1:q2, 1:(q2 - 1)] <- D1
  K1 <- t(Zc)%*%K[[1]]%*%Zc
  K2 <- t(Zc)%*%K[[2]]%*%Zc
  
  
  diag(K1) <- diag(K1) + add_to_diag
  diag(K2) <- diag(K2) + add_to_diag
  
  vars <- as.list(substitute(list(...)))[-1]
  
  term1 <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  term2 <- deparse(vars[[2]], backtick = TRUE, width.cutoff = 500)
  
  
  knots <- sm$knots
  
  Bp01 <- te.BD_pya(y=data[,term1], x=data[,term2], knotsy=knots[[1]], knotsx=knots[[2]], q1=q1, q2=q2, center=T )

  exp_ident <- sm$p.ident
  label <- paste("hyx_sm(", vars[[1]], ", ", vars[[2]], ")", 
                 sep = "")
  colnames(B01) <-  paste(label, ".", 1:NCOL(B01), sep = "")
  
  # if(!is.null(add_to_diag)) diag(K1) <- diag(K1) + add_to_diag;   diag(K2) <- diag(K2) + add_to_diag
  
  rank <- lapply(list(K1, K2), rankMatrix)  
  attr(K1, "label") <-  attr(K2, "label") <- label
  attr(K1, "rank") <- rank[[1]]
  attr(K2, "rank") <- rank[[2]]
  
  
  pred_B <- function(newdata) {
    PredictMat(sm, newdata)%*%Zc
  }
  pred_B <- custom_fun_env(pred_B, sm=sm, Zc = Zc)
  
  
  pred_Bp <- function(newdata) {
    te.BD_pya(y=newdata[,term1], x=newdata[,term2], knotsy=knots[[1]], knotsx=knots[[2]], q1=q1, q2=q2, center=T )
  }
  pred_Bp <- custom_fun_env(pred_Bp, term1=term1, term2=term2, knots=knots, q1=q1, q2=q2)
  
  ret <- list(B=B01, Bp = Bp01, K = list(K1, K2), label = label,label_inds = rep(label, NCOL(B01)), knots = knots, exp_ident = exp_ident, pen_ind = TRUE,
              rank = rank, predvars = list(B = pred_B, Bp = pred_Bp))
  ret
}



# nonlinear monotone transformations with linear VCM interaction -----------------------------------------------------------------------------------
hyx_vcm <- function(..., by =NA, q =20,  knots = NULL, data,
                    pen_order = 2, degree = 3, center = F, add_to_diag = 10e-6 )
{
  
  if (pen_order != 2)  stop("right now only default penalties possible")
  if (degree != 3)  stop("right now only cubic splines possible")

  
  by.var <- deparse(substitute(by), backtick = TRUE, width.cutoff = 500)
  sspec <- s(..., k = q , bs = "ps")
  sspec$mono <- 1
  
  sm <- smoothCon(sspec,
                  data = data,
                  scale.penalty=F,
                  absorb.cons = center ,
                  knots = knots)[[1]]
  
  B <- sm$X
  K <- sm$S[[1]]
  
  (sm$deriv <- 1)
  Bp <- PredictMat(sm, data)
  
  x <- data[,by.var]
  B <- B*x
  Bp <- Bp*x
  
  vars <- as.list(substitute(list(...)))[-1]
  
  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  
  
  if (!is.null(add_to_diag))
    diag(K) <- diag(K) + add_to_diag
  
  
  exp_ident <- c(center, rep(T, ncol(B)-1))
  
  rank <- rankMatrix(K)
  
  label <- paste("hyx_vcm(", term, ", by=", by.var,")", 
                 sep = "")
  colnames(B) <-  paste(label, ".", 1:NCOL(B), sep = "")
  
  
  
  pred_B <- function(newdata) {
    sm$deriv <- 0
    PredictMat(sm, newdata)*newdata[,by.var]
  }
  pred_B <- custom_fun_env(pred_B,sm=sm, by.var=by.var)
  
  
  pred_Bp <- function(newdata) {
    sm$deriv <- 1
    PredictMat(sm, newdata)*newdata[,by.var]
  }
  pred_Bp <- custom_fun_env(pred_Bp, sm=sm, by.var=by.var)
  
  ret <- list(B=B, Bp = Bp, K = list(K), label = label, label_inds = rep(label, NCOL(B)), knots = knots, exp_ident = exp_ident, pen_ind = TRUE,
              rank = rank, predvars = list(B = pred_B, Bp = pred_Bp))
  # attr(ret, "pen_ind") <- F
  
  ret
}



# nonlinear monotone transformations ------------------------------------------------------------------------------------------
hy_sm <- function(..., q = 20,  knots = NULL, data,
                  pen_order = 2, degree = 3, center = T, add_to_diag = 10e-6 )
{
  
  if (pen_order != 2)
    stop("right now only default penalties possible")
  if (degree != 3)
    stop("right now only cubic splines possible")
  # if (!center)
    # stop("right now only center == TRUE allowed")
  
  
  
  sspec <- s(..., k = q , bs = "ps")
  sspec$mono <- 1
  
  sm <- smoothCon(sspec,
                  data = data,
                  scale.penalty=F,
                  absorb.cons = center ,
                  knots = knots)[[1]]
  
  B <- sm$X
  K <- sm$S[[1]]
  

  sm$deriv <- 1
  Bp <- PredictMat(sm, data)
  
  
  B <- B
  Bp <- Bp
  
  vars <- as.list(substitute(list(...)))[-1]
  
  term <- deparse(vars[[1]], backtick = TRUE, width.cutoff = 500)
  
  
  if (!is.null(add_to_diag))
    diag(K) <- diag(K) + add_to_diag
  
  
  exp_ident <- c(center, rep(T, ncol(B)-1))

  rank <- rankMatrix(K)
  
  label <- paste("hy_sm(", term, ")",
                 sep = "")
  colnames(B) <-  paste(label, ".", 1:NCOL(B), sep = "")
  
  
  
  pred_B <- function(newdata) {
    sm$deriv <- 0
    PredictMat(sm, data.frame(newdata))
  }
  pred_B <- custom_fun_env(pred_B,  sm=sm)
  
  
  pred_Bp <- function(newdata) {
    sm$deriv <- 1
    PredictMat(sm, data.frame(newdata))
  }
  pred_Bp <- custom_fun_env(pred_Bp,  sm=sm)
  
  ret <- list(B=B, Bp = Bp, K = list(K), label = label, label_inds = rep(label, NCOL(B)), knots = knots, exp_ident = exp_ident, pen_ind = TRUE,
              rank = rank, predvars = list(B = pred_B, Bp = pred_Bp))
  # attr(ret, "pen_ind") <- F
  
  ret
}



# simple linear effects ----------------------------------------------------------------------------------------------------------
hx_lin  <- function(x, scale = F, center = F){
  
  x <- scale(x, center, scale)
  B <- matrix(x, ncol=1)
  label <- paste("hx_lin(", match.call()$x, ")", sep = "")
  colnames(B) <- label
  
  Bp <- matrix(0, length(x), ncol=1)
  K <- structure(matrix(0), rank = 0, label = label)
  
  rng <- range(x)
  
  unq <- !duplicated(x)
  xu <- x[unq]
  construct_Bpred <- function(xnew) {
    # if (any(xnew < rng[1] | xnew > rng[2])) {
    #   warning("predictions outside fitted range for ", 
    #           label, ".")
    # }
    unlist(xnew)
  }
  construct_Bpred <- custom_fun_env(construct_Bpred, rng = rng, label = label, xu = xu)
  
  construct_Bpredp <- function(xnew) {
    # if (any(xnew < rng[1] | xnew > rng[2])) {
    #   warning("predictions outside fitted range for ", 
    #           label, ".")
    # }
    rep(0,length(unlist(xnew)))
  }
  construct_Bpredp <- custom_fun_env(construct_Bpredp, rng = rng, label = label, xu = xu)
  ret <- list(B = B, Bp = Bp, label = label, label_inds = rep(label, NCOL(B)), K = list(K), exp_ident = F, pen_ind = FALSE,
              predvars = list(B = construct_Bpred, Bp = construct_Bpredp))
  # attr(ret, "pen_ind") <- T
  ret
}



# spatial effect
hx_spat <- function(s, nmat = NULL, data){
  
  # if (center)  stop("right vcm-trafo cannot be centered")
  
  if(is.null(nmat)) stop("gra file needed for nmat")
  var <- deparse(match.call()$s)
  
  regions <- as.numeric(attr(nmat,"dimnames")[[1]])
  
  nreg <- length(regions)
  
  kspatmat <- matrix(unlist(nmat), nreg,nreg)
  # attr(nmat, "class") <- "mat"
  
  # regions <- as.numeric(attr(nmat,"dimnames")[[1]])[-1]

  
  B <- matrix(0,nrow(data),nreg)
  for(j in 1:nreg) {
    B[which(data[,var]==regions[j]),j] <- 1
  }
  
 
 Bp <- matrix(0, nrow(B), ncol(B))
  # K0 <- matrix(0, m+2, m+2)
  
  # class(nmat) <- "matr"
  # K <- as.matrix(nmat)[-1,-1]
 K <- matrix(nmat, nreg, nreg)
 
  class(K) <- "matrix"
  rank <- rankMatrix(K)

  exp_ident <- rep(F, ncol(B))
  
  
  label <- paste("hs(", deparse(match.call()$s),")", 
                 sep = "")
  colnames(B) <-  paste(label, ".", 1:NCOL(B), sep = "")
  
  
  # 
  pred_B <- function(newdata) {
    # B_pred<- matrix(0, nrow(newdata), nreg)
    newdata
  }
  pred_B <- custom_fun_env(pred_B, nreg=nreg)


  pred_Bp <- function(newdata) {
    Bp_pred<- matrix(0, nrow(newdata), nreg)
  }
  pred_Bp <- custom_fun_env(pred_Bp, nreg=nreg)
  
  ret <- list(B=B, Bp = Bp, K = list(K), label = label, label_inds = rep(label, NCOL(B)), knots = NULL, exp_ident = exp_ident, pen_ind = TRUE,
              rank = rank, predvars = list(B = pred_B, Bp = pred_Bp))
  # attr(ret, "pen_ind") <- F
  
  ret
}
