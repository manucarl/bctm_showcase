// posterior and gradient functions for Gaussian, logistic and minimum-extreme-value reference distribution

// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;




//[[Rcpp::export]]
double posterior_gauss(const arma::vec& params, const Rcpp::List& xx){
  
  uvec ind_exp = xx["exp_ident"];
  
  uvec m_pen_ident = xx["m_pen_ident"];
  
  
  mat X = xx["X"];
  mat Xp = xx["Xp"];
  mat XtX = xx["XtX"];
  
  mat S = xx["S"];
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  vec sparams = params(m_pen_ident);
  double pd = as_scalar(-0.5*bt.t()*(XtX*bt)+ sum(log(Xp*bt)) - 0.5*sparams.t()*(S*sparams));
  // double pd = as_scalar(-0.5*bt_vec.t()*XtX*bt_vec+ sum(log(Xp*bt_vec)));
  
  return pd;
}

// [[Rcpp::export]]
arma::vec gradf_gauss(const arma::vec &params, const Rcpp::List &xx) {
  
  mat X = xx["X"];
  mat Xp = xx["Xp"];
  mat XtX = xx["XtX"];
  
  mat S = xx["S"];
  uvec ind_exp = xx["exp_ident"];
  uvec m_pen_ident = xx["m_pen_ident"];
  
  arma::vec bt_vec= params;
  bt_vec(ind_exp) = exp(params(ind_exp));
  
  arma::vec w = 1/(Xp*bt_vec);
  
  int qb = bt_vec.n_elem;
  
  
  vec c = vec(qb);
  c.ones();
  c(ind_exp) = bt_vec(ind_exp);
  
  mat Xpt = Xp.t();
  
  
  
  arma::vec p1 = -(XtX*bt_vec)%c;
  
  
  arma::vec p2 = (Xpt*w)%c;
  
  vec p3 = vec(qb);
  p3.zeros();
  
  p3(m_pen_ident) = S*params(m_pen_ident);
  arma::vec res = p1 + p2 -p3;
  
  return res;
}


// [[Rcpp::export]]
double posterior_logit(const arma::vec &params, const Rcpp::List &xx) {
  
  uvec ind_exp = xx["exp_ident"];
  uvec m_pen_ident = xx["m_pen_ident"];
  
  arma::vec bt= params;
  bt.elem(ind_exp) = exp(params.elem(ind_exp));
  
  mat X = xx["X"];
  mat Xp = xx["Xp"];
  mat XtX = xx["XtX"];
  
  mat S = xx["S"];
  arma::vec bt_vec= params;
  bt_vec.elem(ind_exp) = exp(params.elem(ind_exp));
  vec sparams = params(m_pen_ident);
  
  arma::vec eXbt = exp(X*bt_vec);
  arma::vec likelihoods = log(eXbt/pow(eXbt+1, 2.0))+ log(Xp*bt_vec);
  double ll = sum(likelihoods)- as_scalar(0.5*sparams.t()*(S*sparams));
  
  return ll;
}





// [[Rcpp::export]]
vec gradf_logit(const arma::vec &params, const Rcpp::List &xx) {
  
  mat X = xx["X"];
  mat Xp = xx["Xp"];
  mat XtX = xx["XtX"];
  mat S = xx["S"];
  uvec ind_exp = xx["exp_ident"];
  uvec m_pen_ident = xx["m_pen_ident"];
  
  arma::vec bt_vec= params;
  bt_vec.elem(ind_exp) = exp(params.elem(ind_exp));
  
  arma::vec w = 1/(Xp*bt_vec);
  
  int qb = bt_vec.n_elem;
  
  
  vec diagvec = vec(qb);
  diagvec.ones();
  diagvec(ind_exp) = bt_vec(ind_exp);

  arma::vec eh = exp(X*bt_vec);
  
  vec p1 = trans(X.each_row()%diagvec.t())*(-(eh-1)/(eh+1));
  vec p2 = trans(Xp.each_row()%diagvec.t())*w;
  
  vec p3 = vec(qb);
  p3.zeros();
  p3(m_pen_ident) = S*params(m_pen_ident);
  
  
  arma::vec res = p1 + p2 -p3;
  
  return res;
}


// [[Rcpp::export]]
double posterior_logit_cens(const arma::vec &params, const Rcpp::List &xx) {
  
  mat X1 = xx["X1"];
  mat X2 = xx["X2"];
  mat Xp1 = xx["Xp1"];
  
  mat S = xx["S"];
  uvec ind_exp = xx["exp_ident"];
  uvec m_pen_ident = xx["m_pen_ident"];
  
  arma::vec bt_vec= params;
  bt_vec.elem(ind_exp) = exp(params.elem(ind_exp));
  
  vec eXbt = exp( X1*bt_vec);
  
  vec ll1 = log(eXbt/pow(eXbt+1, 2.0))+ log(Xp1*bt_vec);
  

  vec sparams = params(m_pen_ident);
  vec ll2 = log(1 - plogis(as<NumericVector>(wrap(X2*bt_vec)) ));
  
  double lp = sum(ll1) + sum(ll2) - as_scalar(0.5*sparams.t()*(S*sparams));

  return lp;
}


// [[Rcpp::export]]
double ll_logit_cens(const arma::vec &params, const Rcpp::List &xx) {
  
  mat X1 = xx["X1"];
  mat X2 = xx["X2"];
  mat Xp1 = xx["Xp1"];
  
  uvec ind_exp = xx["exp_ident"];
  uvec m_pen_ident = xx["m_pen_ident"];
  
  arma::vec bt_vec= params;
  bt_vec.elem(ind_exp) = exp(params.elem(ind_exp));
  
  vec eXbt = exp( X1*bt_vec);
  
  vec ll1 = log(eXbt/pow(eXbt+1, 2.0))+ log(Xp1*bt_vec);
  
  // Rcout << "The value is " << ll1;
  
  // vec ll2 = log(1 - plogis(X2*bt_vec, true, false ));
  //NumericVector p=X2*bt_vec;
  vec ll2 = log(1 - plogis(as<NumericVector>(wrap(X2*bt_vec)) ));
  
  double ll = sum(ll1) + sum(ll2);
  
  return ll;
}


// [[Rcpp::export]]
vec gradf_logit_cens(const arma::vec &params, const Rcpp::List &xx) {
  
  mat X1 = xx["X1"];
  mat X2 = xx["X2"];
  mat Xp = xx["Xp1"];
  mat S = xx["S"];
  int p = xx["p"];
  
  uvec ind_exp = xx["exp_ident"];
  uvec m_pen_ident = xx["m_pen_ident"];
  
  arma::vec bt_vec= params;
  bt_vec.elem(ind_exp) = exp(params.elem(ind_exp));
  
  arma::vec w = 1/(Xp*bt_vec);
  
  vec diagvec = vec(p);
  diagvec.ones();
  diagvec(ind_exp) = bt_vec(ind_exp);
  
  mat C = diagmat(diagvec);
  
  arma::vec eh = exp(X1*bt_vec);
  
  arma::vec p1 = trans(X1*C)*(-(eh-1)/(eh+1));
  arma::vec p2 = trans(Xp*C)*w;
  vec p3 = vec(p);
  
  p3.zeros();
  p3(m_pen_ident) = S*params(m_pen_ident);
  arma::vec res1 = p1 + p2 - p3;
  
  vec par =  -dlogis(as<NumericVector>(wrap(X2*bt_vec)))/(1-plogis(as<NumericVector>(wrap(X2*bt_vec))));
  vec res2 = C*trans(X2) * par;
  return res1 + res2;
}

// [[Rcpp::export]]
vec rowSums(const mat & X){
  int nRows = X.n_rows;
  vec out(nRows);
  for(int i = 0; i < nRows; i++){
    out(i) = sum(X.row(i));
  }
  return(out);
}

// [[Rcpp::export]]
double posterior_mev(const arma::vec &params, const Rcpp::List &xx) {
  
  mat X = xx["X"];
  mat Xp = xx["Xp"];
  
  mat S = xx["S"];
  uvec ind_exp = xx["exp_ident"];
  uvec m_pen_ident = xx["m_pen_ident"];
  
  arma::vec bt_vec= params;
  bt_vec.elem(ind_exp) = exp(params.elem(ind_exp));
  vec sparams = params(m_pen_ident);
  vec Xbt = X*bt_vec;
  
  arma::vec likelihoods = Xbt - exp(Xbt) + log(Xp*bt_vec);
  double ll = sum(likelihoods)- as_scalar(0.5*sparams.t()*(S*sparams));
  
  return ll;
}



// [[Rcpp::export]]
vec pmev(vec x) {  
 vec ret = 1 - exp(-exp(x));

  return ret;
}

// [[Rcpp::export]]
vec ldmev(vec x) {  
  vec ret = x - exp(x);
  return ret;
}

// [[Rcpp::export]]
vec dmev(vec x) {  
  return exp(x);
}

// [[Rcpp::export]]
double posterior_mev_cens(const arma::vec &params, const Rcpp::List &xx) {
  
  
  mat X1 = xx["X1"];
  mat X2 = xx["X2"];
  mat Xp1 = xx["Xp1"];
  
  mat S = xx["S"];
  uvec ind_exp = xx["exp_ident"];
  uvec m_pen_ident = xx["m_pen_ident"];
  
  arma::vec bt_vec= params;
  bt_vec.elem(ind_exp) = exp(params.elem(ind_exp));
  vec sparams = params(m_pen_ident);
  
  arma::vec ll1 = X1*bt_vec - exp(X1*bt_vec) + log(Xp1*bt_vec);
  arma::vec ll2 = log(1 - pmev(as<NumericVector>(wrap(X2*bt_vec)) ));
  double lp = sum(ll1) + sum(ll2)- as_scalar(0.5*sparams.t()*(S*sparams));
  
  return lp;
}


// [[Rcpp::export]]
vec gradf_mev_cens(const arma::vec &params, const Rcpp::List &xx){
  
  
  mat X1 = xx["X1"];
  mat X2 = xx["X2"];
  mat Xp = xx["Xp1"];
  int p = xx["p"];
  
  mat S = xx["S"];
  uvec ind_exp = xx["exp_ident"];
  uvec m_pen_ident = xx["m_pen_ident"];
  
  vec bt_vec= params;
  bt_vec.elem(ind_exp) = exp(params.elem(ind_exp));
  
  arma::vec w = 1/(Xp*bt_vec);
  
  
  vec diagvec = vec(p);
  diagvec.ones();
  diagvec(ind_exp) = bt_vec(ind_exp);
  mat C = diagmat(diagvec);
  
  // mat C = diagmat(diagvec);
  
  arma::vec eh = exp(X1*bt_vec);
  

  
  mat CXt = trans(X1.each_row()%diagvec.t());
  mat XpC = trans(Xp.each_row()%diagvec.t());
  // Rcout << "The value is " <<eye(n,n);
  // rowvec  = rowvec(1-eh);
  vec p1 = rowSums(CXt.each_row()%(1-eh).t() + XpC.each_row()%w.t());
  
  vec par =  -dmev(as<NumericVector>(wrap(X2*bt_vec)))/(1-pmev(as<NumericVector>(wrap(X2*bt_vec))));
  vec p2 = C*trans(X2) * par;
  
  vec p3 = vec(p);
  p3.zeros();
  p3(m_pen_ident) = S*params(m_pen_ident);
  
  return  p1 + p2 - p3;
}







// [[Rcpp::export]]
vec gradf_mev(const arma::vec &params, const Rcpp::List &xx){
  
  
  mat X = xx["X"];
  mat Xp = xx["Xp"];
  int p = xx["p"];
  
  mat S = xx["S"];
  uvec ind_exp = xx["exp_ident"];
  uvec m_pen_ident = xx["m_pen_ident"];
  
  vec bt_vec= params;
  bt_vec.elem(ind_exp) = exp(params.elem(ind_exp));
  
  arma::vec w = 1/(Xp*bt_vec);
  
  
  vec diagvec = vec(p);
  diagvec.ones();
  diagvec(ind_exp) = bt_vec(ind_exp);
  
  // mat C = diagmat(diagvec);
  
  arma::vec eh = exp(X*bt_vec);
  
  vec p3 = vec(p);
  p3.zeros();
  p3(m_pen_ident) = S*params(m_pen_ident);
  
  mat CXt = trans(X.each_row()%diagvec.t());
  mat XpC = trans(Xp.each_row()%diagvec.t());
  // Rcout << "The value is " <<eye(n,n);
  // rowvec  = rowvec(1-eh);
  vec res = rowSums(CXt.each_row()%(1-eh).t() + XpC.each_row()%w.t())- p3;
  
  return  res;
}







// [[Rcpp::depends(RcppArmadillo)]]