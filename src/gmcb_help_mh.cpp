#include <RcppArmadillo.h>
using namespace Rcpp;

// Taken from Dirk Eddelbuettel's answer here: 
// https://stackoverflow.com/questions/28442582/reproducing-r-rep-with-the-times-argument-in-c-and-rcpp
NumericVector rep_times_mh(NumericVector x, NumericVector y) {
  int n = y.size();
  NumericVector out(sum(y));
  int ind=0;
  for (int i=0; i < n; ++i) {
    int p = y[i];
    std::fill(out.begin()+ind, out.begin()+ind+p, x[i]);
    ind += p;
  }
  return out;
}

// log posterior conditional for alpha.b
double log_post_cond_alpha_b_mh(double alpha_b, NumericVector b, NumericVector rep_gamma,
                             NumericVector lambda, double k1, double k2) {
  if (alpha_b < k1 || alpha_b > k2) {
    return -std::numeric_limits<double>::infinity();
  } else {
    int pq = b.length();
    double exponent_term = -0.5*Rcpp::sum(lambda*Rcpp::pow(Rcpp::abs(b), alpha_b)/rep_gamma);
    double pq_term = pq * (std::log(alpha_b) - 1/alpha_b*std::log(2.0) - std::lgamma(1.0/alpha_b));
    return pq_term + 1.0/alpha_b*sum(log(lambda) - log(rep_gamma)) + exponent_term;
  }
}

// log posterior conditional for alpha.d
double log_post_cond_alpha_d_mh(double alpha_d, NumericVector delta, NumericVector tau,
                             NumericVector gamma, double k1, double k2) {
  if (alpha_d < k1 || alpha_d > k2) {
    return -std::numeric_limits<double>::infinity();
  } else {
    int q = gamma.length();
    NumericVector reps = wrap(arma::linspace(1, q-1, q-1));
    NumericVector gamma_minus1 = NumericVector::import(gamma.begin() + 1, gamma.end());
    NumericVector calc_gamma = rep_times_mh(gamma_minus1, reps);
    
    double alpha_only_term = q*(q-1)/2*(std::log(alpha_d) - 1.0/alpha_d*std::log(2.0) - std::lgamma(1.0/alpha_d));
    double gamma_term = Rcpp::sum(reps * Rcpp::log(gamma_minus1))/alpha_d;
    double tau_term = 1.0/alpha_d * Rcpp::sum(Rcpp::log(tau));
    double exp_term = -0.5*Rcpp::sum(tau/calc_gamma * Rcpp::pow(Rcpp::abs(delta), alpha_d));
    
    return alpha_only_term - gamma_term + tau_term + exp_term;
  }
}

// convert a q(q-1)/2 vector into a unit lower triangular matrix
arma::mat delta_to_matrix_inner_mh(NumericVector &delta, int q) {
  arma::vec delta_arma(delta.begin(), delta.size(), false);
  arma::mat out(q, q, arma::fill::eye);
  arma::uvec upper_indices = trimatu_ind( size(out) , 1); // don't include indices of diagonal
  out.elem(upper_indices) = -delta_arma;
  return out.t();
}

// calculate the precision or covariance matrix from delta and gamma
arma::mat sigmainv_calc_mh(NumericVector &delta, NumericVector &gamma,
                        bool cov) {
  int q = gamma.length();
  arma::mat unit_t = delta_to_matrix_inner_mh(delta, q); //T
  
  // construct D^(-1)
  arma::vec gamma_arma = as<arma::vec>(gamma);
  // arma::mat d = diagmat(gamma_arma);
  // arma::mat dinv = diagmat(1/gamma_arma);
  
  arma::mat out(q, q);
  if (cov) { // return covariance matrix
    arma::mat u = arma::inv(trimatl(unit_t));
    out = u * arma::diagmat(gamma_arma) * u.t(); 
  } else { // return precision matrix
    out = unit_t.t() * arma::diagmat(1/gamma_arma) * unit_t;
  }
  return out;
}

// log posterior conditional of B
double b_logpostcond(int index, const arma::mat &yty, const arma::mat &xtx, 
                     const arma::mat &ytx, NumericVector b,
                     double bij, NumericVector delta, NumericVector gamma,
                     double lambdaij, double alpha_b) {
  
  NumericVector b_test = b; // Rcpp does not pass by value, even though it may look like it in argument list
  b_test(index) = bij; // index will have a value between 0 and pq - 1, inclusive
  
  int gamma_index = index/xtx.n_cols;
  
  // reshape b for matrix calculations
  arma::vec b_arma(b.begin(), b.size(), false); 
  arma::mat b_mat = reshape(b_arma, xtx.n_cols, yty.n_cols);
  
  arma::mat omega = sigmainv_calc_mh(delta, gamma, false);
  
  return (-0.5 * (arma::trace((yty - 2.0 * ytx * b_mat + b_mat.t() * xtx * b_mat) * omega) + 
          lambdaij * std::pow(std::abs(bij), alpha_b)/gamma(gamma_index)));
}

// log posterior conditional of B
double b_logpostcond_direct(int index, const arma::mat &x, 
                            const arma::mat &y, NumericVector b,
                            double bij, NumericVector delta, NumericVector gamma,
                            double lambdaij, double alpha_b) {
  
  NumericVector b_test = b; // Rcpp does not pass by value, even though it may look like it in argument list
  b_test(index) = bij; // index will have a value between 0 and pq - 1, inclusive
  
  int gamma_index = index/x.n_cols;
  
  // reshape b for matrix calculations
  arma::vec b_arma(b.begin(), b.size(), false); 
  arma::mat b_mat = reshape(b_arma, x.n_cols, y.n_cols);
  
  arma::mat omega = sigmainv_calc_mh(delta, gamma, false);
  arma::mat ymxb = y - x * b_mat;
  
  return (-0.5 * (arma::trace(ymxb * omega * ymxb.t()) + 
          lambdaij * std::pow(std::abs(bij), alpha_b)/gamma(gamma_index)));
}

// // log posterior conditional of B
// double b_logpostcond_mixedcase_pln(int p, int q, int n, int index, 
//                                    const arma::mat &x, const arma::mat &y,
//                                    const arma::mat &ytx, const arma::mat &xtx,
//                                    NumericVector b,
//                                    double bij, NumericVector delta, NumericVector gamma,
//                                    double lambdaij, double alpha_b) {
//   
//   NumericVector b_test = b; // Rcpp does not pass by value, even though it may look like it in argument list
//   b_test(index) = bij; // index will have a value between 0 and pq - 1, inclusive
//   
//   int gamma_index = index/p;
//   
//   // reshape b for matrix calculations
//   arma::vec b_arma(b.begin(), b.size(), false); 
//   arma::mat b_mat = reshape(b_arma, p, q);
//   
//   arma::mat omega = sigmainv_calc_mh(delta, gamma, false);
//   
//   arma::mat y_omega_yt = y * omega * y.t();
//   arma::mat omega_bt_xty = omega * b_mat.t() * ytx.t();
//   arma::mat xtx_term = b_mat * omega * b_mat.t() * xtx;
//   
//   double traceterm = arma::trace(y_omega_yt) - 2 * arma::trace(omega_bt_xty) + 
//     arma::trace(xtx_term);
//   
//   return (-0.5 * (traceterm + 
//           lambdaij * std::pow(std::abs(bij), alpha_b)/gamma(gamma_index)));
// }
// 
// // log posterior conditional of B
// double b_logpostcond_mixedcase_qln(int p, int q, int n, int index, 
//                                    const arma::mat &x, const arma::mat &y, 
//                                    const arma::mat &yty, NumericVector b,
//                                    double bij, NumericVector delta, NumericVector gamma,
//                                    double lambdaij, double alpha_b) {
//   
//   NumericVector b_test = b; // Rcpp does not pass by value, even though it may look like it in argument list
//   b_test(index) = bij; // index will have a value between 0 and pq - 1, inclusive
//   
//   int gamma_index = index/p;
//   
//   // reshape b for matrix calculations
//   arma::vec b_arma(b.begin(), b.size(), false); 
//   arma::mat b_mat = reshape(b_arma, p, q);
//   
//   arma::mat omega = sigmainv_calc_mh(delta, gamma, false);
//   
//   arma::mat ytyomega = yty * omega;
//   
//   arma::mat xb = x * b_mat;
//   arma::mat yomega = y * omega;
//   arma::mat yomegaxb = yomega * xb.t();
//   
//   arma::mat xbomegaxbt = xb * omega * xb.t();
//   
//   double traceterm = arma::trace(ytyomega) - 2 * arma::trace(yomegaxb) +
//     arma::trace(xbomegaxbt);
//   
//   return (-0.5 * (traceterm + 
//           lambdaij * std::pow(std::abs(bij), alpha_b)/gamma(gamma_index)));
// }

// log posterior conditional for delta_c
double delta_logpostcond(int c, int j, NumericVector deltac, double deltacj,
                         const arma::mat &yctyc, const arma::mat &yctycm1, 
                         const arma::mat &ycm1tycm1, const arma::mat &yctx, 
                         const arma::mat &ycm1tx, const arma::mat&xtx,
                         arma::mat &b,
                         double gammac, double alphad, double taucj) {
  
  // note that the c, j are C++ indices
  
  // convert to arma objects for computation
  deltac(j) = deltacj;
  arma::colvec deltac_arma(deltac.begin(), deltac.size(), false); 
  
  double l2term = -0.5/gammac * arma::accu(yctyc + b.col(c).t() * xtx * b.col(c) +
                                           deltac_arma.t() * (ycm1tycm1 - ycm1tx * b.cols(0, c - 1) - 
                                           b.cols(0, c - 1).t() * ycm1tx.t() + 
                                           b.cols(0, c - 1).t() * xtx * b.cols(0, c-1)) * deltac_arma -
                                           2.0 * yctx * b.col(c) - 2.0 * (yctycm1 - yctx * b.cols(0, c - 1)) * deltac_arma +
                                           2.0 * deltac_arma.t() * (ycm1tx - b.cols(0, c - 1).t() * xtx) * b.col(c));
  double prior_term = -0.5/gammac * taucj * std::pow(std::abs(deltacj), alphad);
  
  return l2term + prior_term;
}  

// log posterior conditional for delta_c - direct calculation
double delta_logpostcond_direct(int c, int j, NumericVector deltac, double deltacj,
                                const arma::mat &y, const arma::mat &x, 
                                arma::mat &b,
                                double gammac, double alphad, double taucj) {
  
  // note that the c, j are C++ indices
  
  // convert to arma objects for computation
  deltac(j) = deltacj;
  arma::colvec deltac_arma(deltac.begin(), deltac.size(), false); 
  
  double l2term = -0.5/gammac * arma::accu(arma::square(y.col(c) - x * b.col(c) -
                                           (y.cols(0, c - 1) - x * b.cols(0, c - 1)) * deltac_arma));
  
  double prior_term = -0.5/gammac * taucj * std::pow(std::abs(deltacj), alphad);
  
  return l2term + prior_term;
}  

// log posterior conditional for delta_c
double delta_meanzero_logpostcond(int c, int j, NumericVector deltac, double deltacj,
                                  const arma::mat &yctyc, const arma::mat &yctycm1, 
                                  const arma::mat &ycm1tycm1,
                                  double gammac, double alphad, double taucj) {
  
  // note that the c, j are C++ indices
  
  // convert to arma objects for computation
  deltac(j) = deltacj;
  arma::colvec deltac_arma(deltac.begin(), deltac.size(), false); 
  
  double l2term = -0.5/gammac * arma::accu(yctyc + deltac_arma.t() * ycm1tycm1 * deltac_arma 
                                             - 2.0 * yctycm1 * deltac_arma);
  double prior_term = -0.5/gammac * taucj * std::pow(std::abs(deltacj), alphad);
  
  return l2term + prior_term;
}  

// log posterior conditional for delta_c - direct calculation
double delta_meanzero_logpostcond_direct(int c, int j, NumericVector deltac, double deltacj,
                                         const arma::mat &y, 
                                         double gammac, double alphad, double taucj) {
  
  // note that the c, j are C++ indices
  
  // convert to arma objects for computation
  deltac(j) = deltacj;
  arma::colvec deltac_arma(deltac.begin(), deltac.size(), false); 
  
  double l2term = -0.5/gammac * arma::accu(arma::square(y.col(c) - y.cols(0, c - 1) * deltac_arma));
  
  double prior_term = -0.5/gammac * taucj * std::pow(std::abs(deltacj), alphad);
  
  return l2term + prior_term;
}  

std::map<std::string,arma::mat> mh_precompute(int p, int q, int n, const arma::mat &x, 
                                              const arma::mat &y) {
  std::map<std::string, arma::mat>m;
  
  if (p < n) {
    arma::mat xtx = x.t() * x;
    arma::mat ytx = y.t() * x;
    
    m.insert(std::pair<std::string, arma::mat>("xtx", xtx));
    m.insert(std::pair<std::string, arma::mat>("ytx", ytx));
    
    if (q < n) {
      arma::mat yty = y.t() * y;
      
      m.insert(std::pair<std::string, arma::mat>("yty", yty));
      
      return(m);
    } else {
      arma::mat y_firstn = y.cols(0, n - 1); // select the first n columns
      arma::mat yty_firstn = y_firstn.t() * y_firstn;
      
      m.insert(std::pair<std::string, arma::mat>("yty", yty_firstn));
      
      return(m);
    }
    
  } else {
    
    if (q < n) {
      arma::mat yty = y.t() * y;
      
      m.insert(std::pair<std::string, arma::mat>("yty", yty));
      
      return(m);
    } else {
      arma::mat out;
      out.ones(1,1);
      
      m.insert(std::pair<std::string, arma::mat>("placeholder", out));
      
      return(m);
    }
  }
}