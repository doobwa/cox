# Wrappers for C++ code

require(inline)
require(Rcpp)
require(RcppArmadillo)
inc <- paste(readLines("cox.cpp"),collapse="\n")
src <- '
Rcpp::IntegerVector y(y_r);
Rcpp::List X(X_r);
Rcpp::NumericVector params(params_r);
double lp = llk(y,X,params);
return wrap(lp);
'
cox.llk.cpp <-  cxxfunction(signature(y_r="integer",X_r="list",params_r="numeric"),includes=inc, body = src, plugin="RcppArmadillo")

src <- '
Rcpp::IntegerVector y(y_r);
Rcpp::List X(X_r);
Rcpp::NumericVector params(params_r);
arma::colvec g = grad(y,X,params);
return wrap(g);
'
cox.grad.cpp <-  cxxfunction(signature(y_r="integer",X_r="list",params_r="numeric"),includes=inc, body = src, plugin="RcppArmadillo")
