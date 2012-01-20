// Compute the likelihood for the Cox model: 
//
// \prod_{i=1}^n beta x_i / \sum_j beta x_j
//
// y:      vector of responses, where y_i is in [1,...,J_i]
// X:      an R List object, where each element is a matrix of J_i rows and P columns
// params: vector of length P
//
// Returns: a double representing the loglikelihood of the data given the parameters

double llk(Rcpp::IntegerVector y,Rcpp::List & X, Rcpp::NumericVector params) {

  // Initialize parameter vector as an arma::colvec object
  arma::colvec beta(params.begin(), params.size(), false);

  double ll = 0;
  // Compute the loglikelihood of each observation
  for (int i = 0; i < y.size(); i++) {

    // Create the covariate matrix for this observation
    Rcpp::NumericMatrix xim = X[i];
    arma::mat xi = Rcpp::as<arma::mat>(xim);

    // Grab the row of covariates corresponding to the observed event
    arma::rowvec xiy = xim.row(y[i] - 1);  // y has 1-based indexing

    // Contribution of observed event
    ll += arma::as_scalar(xiy * beta);

    // Contribution of the normalizing constant
    ll -= arma::as_scalar(log(sum(exp(xi * beta))));
  }
  return ll;
}

// Compute the gradient for the Cox model.  
// Arguments are the same as for llk().  See notes for a derivation.
//
// Returns: an arma::colvec object of length P representing the gradient.

arma::colvec grad(Rcpp::IntegerVector y,Rcpp::List & X, Rcpp::NumericVector params) {

  arma::colvec beta(params.begin(), params.size(), false);

  // Initialize gradient vector
  arma::vec g = beta * 0;
  int P = params.size();

  // Loop through each observation
  for (int i = 0; i < y.size(); i++) {
    Rcpp::NumericMatrix xim = X[i];
    arma::mat xi = Rcpp::as<arma::mat>(xim);
    arma::colvec b = exp(xi * beta);

    // For each parameter, add the contribution of observation i
    for (int p = 0; p < P; p++) {

      // Observed event
      g[p] += arma::as_scalar(xi(y[i]-1,p));
      arma::colvec xip = xi.col(p);

      // Normalizing constant
      g[p] -= arma::as_scalar(sum(xip % b)/sum(b));
    }
  }
  return g;
}
