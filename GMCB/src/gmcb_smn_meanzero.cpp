#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include "retstable.h"
#include "gmcb_help_smn.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// with debugging output
// [[Rcpp::export]]
List gmcb_smn_meanzero_debug(const arma::mat &y, 
                             const arma::rowvec &d_init, 
                             const arma::rowvec &gamma_init, 
                             const double &alpha_d_init,
                             const arma::rowvec &tau_init,
                             const arma::vec &gamma_prior, 
                             const arma::vec &alpha_prior, 
                             const arma::mat &tau_prior, 
                             const int &iter, 
                             const double &alpha_d_scale, 
                             const List &pos) {
  
  Timer timer;
  timer.step("start"); // timer start
  
  // for ease of use, extracting the prior parameters
  double a = gamma_prior(0);
  double bg = gamma_prior(1);
  double k1 = alpha_prior(0);
  double k2 = alpha_prior(1);
  
  // Rcout << "prior parameters extracted" << std::endl;
  
  // dimensions
  int n = y.n_rows;
  int q = y.n_cols;
  int lowertcount = q*(q-1)/2;
  int iter_1 = iter + 1;
  
  // Rcout << "dimensions defined" << std::endl;
  
  // to hold output
  arma::mat tau(iter_1, lowertcount);
  arma::mat delta(iter_1, lowertcount);
  arma::mat gamma(iter_1, q);
  arma::vec alpha_d(iter_1);
  
  // Rcout << "output matrices created" << std::endl;
  
  // initialize values
  tau.row(0) = tau_init; 
  delta.row(0) = d_init;
  gamma.row(0) = gamma_init;
  alpha_d(0) = alpha_d_init;
  
  // Rcout << "initialization values saved" << std::endl;
  
  // acceptance indicator for MH steps
  arma::vec alpha_d_accept(iter);
  
  // debugging output, produced based on size of p and q relative to n
  std::map<std::string,arma::mat> debug = debugmatrices_meanzero(q, n, iter, pos);
  
  // computational aids
  arma::mat yty;
  if (q <= n) {
    yty = y.t() * y;
  } else {
    arma::mat y_firstn = y.cols(0, n - 1); // select the first n columns
    yty = y_firstn.t() * y_firstn;
  }
  
  timer.step("pre-sampling preparation complete");
  
  // sampler
  for (int i = 1; i < iter_1; i++) { // start at 1, since initializing values in row 0

    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    // sample epsilon, latent variables for delta
    double half_alpha_d = alpha_d(i-1)/2.0;
    
    arma::vec epsilon(lowertcount);
    for (int c = 1; c < q; c++) {
      IntegerVector pos_j = pos[c-1];
      for (int k = 0; k < c; k++) {
        int c_index = pos_j(k) - 1;
        double epsilon_tilt = std::pow(tau(i-1,c_index)/(2.0*gamma(i-1, c)), 2.0/alpha_d(i-1)) * std::pow(delta(i-1, c_index), 2.0)/2.0;
        
        epsilon(c_index) = retstable_LD(half_alpha_d, std::pow(2.0, half_alpha_d), epsilon_tilt);
      }
    }
    
    debug["epsilon_variates"].row(i-1) = arma::conv_to< arma::rowvec >::from(epsilon);
    
    // sample delta
    for (int c = 1; c < q; c++) {
      // Rcout << "c = " << c << std::endl;
      
      arma::uvec pos_j = pos[c-1];
      arma::uvec pos_j_c = pos_j - arma::ones<arma::uvec>(c); // convert indices from R indices to C++ indices
      
      arma::vec tau_c = arma::conv_to< arma::colvec >::from(tau.submat(i-1, pos_j_c(0), i-1, pos_j_c(pos_j_c.n_elem - 1)));
      arma::vec epsilon_c = epsilon(pos_j_c);
      
      if (c > n - 1) { // not the same c as in R version of the function, so require c + 1 > n 
        arma::vec phiinv = arma::pow(2.0*gamma(i-1,c)/tau_c, 2.0/alpha_d(i-1)) / epsilon_c;
        
        arma::mat wc_gammac = 1/sqrt(gamma(i-1,c)) * y.cols(0, c-1);
        arma::vec zc_gammac = 1/sqrt(gamma(i-1,c)) * y.col(c);
        
        int length_pos_j = pos_j.n_elem;
        arma::vec u(length_pos_j);
        for (int k = 0; k < length_pos_j; k++) {
          u(k) = sqrt(phiinv(k)) * norm_rand();
        }
        
        arma::vec z(n);
        for (int k = 0; k < n; k++) {
          z(k) = norm_rand();
        }
        
        arma::vec v = wc_gammac * u + z;
        arma::mat i_n(n, n, arma::fill::eye);
        // arma::mat phiinv_mat = arma::diagmat(phiinv);
        arma::mat w = solve(wc_gammac * arma::diagmat(phiinv) * wc_gammac.t() + i_n, zc_gammac - v, arma::solve_opts::likely_sympd);
        
        delta.submat(i, pos_j_c(0), i, pos_j_c(pos_j_c.n_elem - 1)) = 
          arma::conv_to< arma::rowvec >::from(u + arma::diagmat(phiinv) * wc_gammac.t() * w);
        
        // indexing for saving the debugging random variates
        arma::uvec npos1 = pos[n-1]; // to re-index, find the element of list pos that corresponds to n + 1 - this index is n in R, so n - 1 in C++
        arma::uvec npos1_c = npos1 - arma::ones<arma::uvec>(n); // indices are R indices, so subtract 1 to obtain the indices for C++
        debug["delta_u_gn"].submat(i-1, pos_j_c(0) - npos1_c(0), 
                         i-1, pos_j_c(pos_j_c.n_elem - 1) - npos1_c(0)) = arma::conv_to< arma::rowvec >::from(u);
        debug["delta_z_gn"].submat(i-1, (c - n) * n, i-1, (c - n) * n + (n-1)) = arma::conv_to< arma::rowvec >::from(z);
        
        
      } else {
        
        arma::vec gammacphi = gamma(i-1,c) * epsilon_c % arma::pow(tau_c / (2.0 * gamma(i-1,c)), 2.0/alpha_d(i-1));
        
        // arma::mat gammacphi_mat = arma::diagmat(gammacphi);
        
        arma::mat wctwc = yty.submat(0, 0, c-1, c-1);
        arma::vec wctzc = yty.submat(0, c, c-1, c);
        
        arma::mat postcond_prec_mat = wctwc + arma::diagmat(gammacphi);
        arma::mat postcond_prec_chol = arma::chol(postcond_prec_mat); // upper triangular form
        arma::vec postcond_mean_v = arma::solve(trimatl(postcond_prec_chol.t()), wctzc);
        arma::vec postcondmean = arma::solve(trimatu(postcond_prec_chol), postcond_mean_v);
        
        arma::vec randomvariates(c);
        for (int k = 0; k < c; k++) {
          randomvariates(k) = norm_rand();
        }
        arma::vec delta_w = solve(trimatu(1/std::sqrt(gamma(i-1,c)) * postcond_prec_chol),
                            randomvariates);
        
        delta.submat(i, pos_j_c(0), i, pos_j_c(pos_j_c.n_elem - 1)) = 
          arma::conv_to< arma::rowvec >::from(postcondmean + delta_w);
        
        if (q > n) {
          debug["delta_z_len"].submat(i-1, pos_j_c(0), i-1, pos_j_c(pos_j_c.n_elem - 1)) = 
            arma::conv_to< arma::rowvec >::from(randomvariates);
        } else {
          debug["delta_z"].submat(i-1, pos_j_c(0), i-1, pos_j_c(pos_j_c.n_elem - 1)) = 
            arma::conv_to< arma::rowvec >::from(randomvariates);
        }
      }
    }
    
    // sampling q(q-1)/2 taus
    for (int c = 1; c < q; c++) { // c indexes the gamma value in R indexing
      IntegerVector pos_j = pos[c-1]; // vector of tau indices
      
      NumericMatrix tau_postcond_par(4, c); // posterior conditional parameters for tau_c
      NumericVector tau_postcond_weight(c);
      
      for (int k = 0; k < c; k++) {
        int c_index = pos_j(k) - 1;
        
        tau_postcond_par(0,k) = tau_prior(0, c_index) + 1.0/alpha_d(i-1);
        tau_postcond_par(1,k) = tau_prior(1, c_index) +
          0.5/gamma(i-1,c)*std::pow(std::abs(delta(i,c_index)), alpha_d(i-1));
        tau_postcond_par(2,k) = tau_prior(2, c_index) + 1.0/alpha_d(i-1);
        tau_postcond_par(3,k) = tau_prior(3, c_index) +
          0.5/gamma(i-1,c)*std::pow(std::abs(delta(i,c_index)), alpha_d(i-1));
        
        double tau_intermediate_weight1 =
          std::pow(tau_prior(1, c_index), tau_prior(0, c_index))/std::tgamma(tau_prior(0, c_index)) *
          std::tgamma(tau_postcond_par(0,k))/std::pow(tau_postcond_par(1,k), tau_postcond_par(0,k));
        
        double tau_intermediate_weight2 =
          std::pow(tau_prior(3, c_index), tau_prior(2, c_index))/std::tgamma(tau_prior(2, c_index)) *
          std::tgamma(tau_postcond_par(2,k))/std::pow(tau_postcond_par(3,k), tau_postcond_par(2,k));
        
        tau_postcond_weight(k) =
          tau_intermediate_weight1/(tau_intermediate_weight1 + tau_intermediate_weight2);
        
        debug["tau_variates"](i-1,c_index) = unif_rand();
      }
      
      for (int k = 0; k < c; k++) {
        int c_index = pos_j(k) - 1;
        
        if (debug["tau_variates"](i-1,c_index) <= tau_postcond_weight(k)) {
          tau(i, c_index) = rgamma(1, tau_postcond_par(0,k), 1.0/tau_postcond_par(1,k))(0);
        } else {
          tau(i, c_index) = rgamma(1, tau_postcond_par(2,k), 1.0/tau_postcond_par(3,k))(0);
        }
      }
    }
    
    // sampling gamma
    for (int c = 0; c < q; c++) {
      if (c == 0) {
        // compute shape
        double shape = n/2.0 + a;
        
        //compute likelihood portion of rate
        double l2term = arma::accu(yty.submat(c, c, c, c));
        
        // compute rate
        double rate = bg + 0.5 * l2term;
        
        // sample
        gamma(i,c) = 1.0/rgamma(1, shape, 1.0/rate)(0);
      } else {
        // extract delta_c and tau_c and format
        arma::uvec temp_pos_j = pos[c-1]; // R indices
        arma::uvec pos_j = temp_pos_j - arma::ones<arma::uvec>(c); // convert to C++ indices
        arma::vec delta_c = arma::conv_to< arma::colvec >::from(delta.submat(i, pos_j(0), i, pos_j(pos_j.n_elem - 1)));
        
        arma::vec tau_c = arma::conv_to< arma::colvec >::from(tau.submat(i, pos_j(0), i, pos_j(pos_j.n_elem - 1)));
        
        // compute shape
        double shape = n/2.0 + a + c/alpha_d(i-1);
        
        // compute likelihood portion of rate
        double l2term;
        if (c > n - 1) {
          // direct calculation
          arma::vec ycm = y.col(c) - y.cols(0, c-1) * delta_c;
          l2term = dot(ycm, ycm);
        } else {
          l2term = arma::accu(yty.submat(c, c, c, c)  +
            delta_c.t() * yty.submat(0, 0, c-1, c-1) * delta_c -
            2.0 * yty.submat(c, 0, c, c - 1) * delta_c);
        }
        
        // compute delta prior portion of rate
        double delta_prior_term = arma::sum(tau_c % arma::pow(arma::abs(delta_c), alpha_d(i-1)));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * delta_prior_term;
        
        // sample
        gamma(i,c) = 1.0/rgamma(1, shape, 1.0/rate)(0);
      }
    }
    
    // sampling alpha.d
    double candidate_alpha_d = alpha_d(i-1) + alpha_d_scale * norm_rand();
    double num_alpha_d = log_post_cond_alpha_d_smn(candidate_alpha_d, arma::conv_to< arma::colvec >::from(delta.row(i)),
                                               arma::conv_to< arma::colvec >::from(tau.row(i)), 
                                               arma::conv_to< arma::colvec >::from(gamma.row(i)), k1, k2);
    double denom_alpha_d = log_post_cond_alpha_d_smn(alpha_d(i-1), arma::conv_to< arma::colvec >::from(delta.row(i)),
                                                 arma::conv_to< arma::colvec >::from(tau.row(i)), 
                                                 arma::conv_to< arma::colvec >::from(gamma.row(i)), k1, k2);
    
    double logmhratio_alpha_d = num_alpha_d - denom_alpha_d;
    double random_u_alpha_d = NA_REAL;
    if (logmhratio_alpha_d >= 0) {
      alpha_d(i) = candidate_alpha_d;
      alpha_d_accept(i-1) = true;
    } else {
      random_u_alpha_d = unif_rand();
      if (random_u_alpha_d < std::exp(logmhratio_alpha_d)) {
        alpha_d(i) = candidate_alpha_d;
        alpha_d_accept(i-1) = true;
      } else {
        alpha_d(i) = alpha_d(i-1);
        alpha_d_accept(i-1) = false;
      }
    }
    
    debug["alpha_d_previous_state"].row(i-1) = alpha_d(i-1);
    debug["alpha_d_proposal"].row(i-1) = candidate_alpha_d;
    debug["alpha_d_log_ratio"].row(i-1) = logmhratio_alpha_d;
    debug["alpha_d_u"].row(i-1) = random_u_alpha_d;
    
  }
  
  timer.step("sampling complete");
  
  List mcmc = List::create(Named("delta") = delta, 
                           Named("gamma") = gamma, 
                           Named("tau") = tau, 
                           Named("alpha_d") = alpha_d);
  
  // outputting as list preserves the matrix structure 
  List debug_out;
  debug_out["epsilon_variates"] = debug["epsilon_variates"];
  debug_out["tau_variates"] = debug["tau_variates"];
  debug_out["alpha_d_previous_state"] = debug["alpha_d_previous_state"];
  debug_out["alpha_d_proposal"] = debug["alpha_d_proposal"];
  debug_out["alpha_d_log_ratio"] = debug["alpha_d_log_ratio"];
  debug_out["alpha_d_u"] = debug["alpha_d_u"];
  
  if (q > n) {
    debug_out["delta_u_gn"] = debug["delta_u_gn"];
    debug_out["delta_z_gn"] = debug["delta_z_gn"];
    debug_out["delta_z_len"] = debug["delta_z_len"];
  } else {
    debug_out["delta_z"] = debug["delta_z"];
  }
  
  List acceptances = List::create(Named("alpha_accept_rate") = arma::mean(alpha_d_accept));
  
  NumericVector res(timer);
  for (int i = 0; i < res.size(); i++) {
    res[i] = res[i]/1e9; // convert from nanoseconds to seconds
  }
  
  return List::create(Named("mcmc") = mcmc, 
                      Named("debug_out") = debug_out,
                      Named("acceptances") = acceptances,
                      Named("timing") = res);
}

// no debug output
// [[Rcpp::export]]
List gmcb_smn_meanzero_nodebug(const arma::mat &y, 
                               const arma::rowvec &d_init, 
                               const arma::rowvec &gamma_init, 
                               const double &alpha_d_init,
                               const arma::rowvec &tau_init,
                               const arma::vec &gamma_prior, 
                               const arma::vec &alpha_prior, 
                               const arma::mat &tau_prior, 
                               const int &iter, 
                               const double &alpha_d_scale, 
                               const List &pos) {
  
  Timer timer;
  timer.step("start"); // timer start
  
  // for ease of use, extracting the prior parameters
  double a = gamma_prior(0);
  double bg = gamma_prior(1);
  double k1 = alpha_prior(0);
  double k2 = alpha_prior(1);
  
  // Rcout << "prior parameters extracted" << std::endl;
  
  // dimensions
  int n = y.n_rows;
  int q = y.n_cols;
  int lowertcount = q*(q-1)/2;
  int iter_1 = iter + 1;
  
  // Rcout << "dimensions defined" << std::endl;
  
  // to hold output
  arma::mat tau(iter_1, lowertcount);
  arma::mat delta(iter_1, lowertcount);
  arma::mat gamma(iter_1, q);
  arma::vec alpha_d(iter_1);
  
  // Rcout << "output matrices created" << std::endl;
  
  // initialize values
  tau.row(0) = tau_init; 
  delta.row(0) = d_init;
  gamma.row(0) = gamma_init;
  alpha_d(0) = alpha_d_init;
  
  // Rcout << "initialization values saved" << std::endl;
  
  // acceptance indicator for MH steps
  arma::vec alpha_d_accept(iter);
  
  // computational aids
  arma::mat yty;
  if (q <= n) {
    yty = y.t() * y;
  } else {
    arma::mat y_firstn = y.cols(0, n - 1); // select the first n columns
    yty = y_firstn.t() * y_firstn;
  }
  
  timer.step("pre-sampling preparation complete");
  
  // sampler
  for (int i = 1; i < iter_1; i++) { // start at 1, since initializing values in row 0
    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    // sample epsilon, latent variables for delta
    double half_alpha_d = alpha_d(i-1)/2.0;
    
    arma::vec epsilon(lowertcount);
    for (int c = 1; c < q; c++) {
      IntegerVector pos_j = pos[c-1];
      for (int k = 0; k < c; k++) {
        int c_index = pos_j(k) - 1;
        double epsilon_tilt = std::pow(tau(i-1,c_index)/(2.0*gamma(i-1, c)), 2.0/alpha_d(i-1)) * std::pow(delta(i-1, c_index), 2.0)/2.0;
        
        epsilon(c_index) = retstable_LD(half_alpha_d, std::pow(2.0, half_alpha_d), epsilon_tilt);
      }
    }
    
    // sample delta
    for (int c = 1; c < q; c++) {
      arma::uvec pos_j = pos[c-1];
      arma::uvec pos_j_c = pos_j - arma::ones<arma::uvec>(c);
      
      arma::vec tau_c = arma::conv_to< arma::colvec >::from(tau.submat(i-1, pos_j_c(0), i-1, pos_j_c(pos_j_c.n_elem - 1)));
      arma::vec epsilon_c = epsilon(pos_j_c);
      
      if (c > n - 1) {
        arma::vec phiinv = arma::pow(2.0*gamma(i-1,c)/tau_c, 2.0/alpha_d(i-1)) / epsilon_c;
        
        arma::mat wc_gammac = 1/sqrt(gamma(i-1,c)) * y.cols(0, c-1);
        arma::vec zc_gammac = 1/sqrt(gamma(i-1,c)) * y.col(c);
        
        int length_pos_j = pos_j.n_elem;
        arma::vec u(length_pos_j);
        for (int k = 0; k < length_pos_j; k++) {
          u(k) = sqrt(phiinv(k)) * norm_rand();
        }
        
        arma::vec z(n);
        for (int k = 0; k < n; k++) {
          z(k) = norm_rand();
        }
        
        arma::vec v = wc_gammac * u + z;
        arma::mat i_n(n, n, arma::fill::eye);
        // arma::mat phiinv_mat = diagmat(phiinv);
        arma::mat w = solve(wc_gammac * diagmat(phiinv) * wc_gammac.t() + i_n, zc_gammac - v, arma::solve_opts::likely_sympd);
        
        delta.submat(i, pos_j_c(0), i, pos_j_c(pos_j_c.n_elem - 1)) = 
          arma::conv_to< arma::rowvec >::from(u + diagmat(phiinv) * wc_gammac.t() * w);
        
      } else {
        
        arma::vec gammacphi = gamma(i-1,c) * epsilon_c % arma::pow(tau_c / (2.0 * gamma(i-1,c)), 2.0/alpha_d(i-1));
        
        // arma::mat gammacphi_mat = diagmat(gammacphi);
        
        arma::mat wctwc = yty.submat(0, 0, c-1, c-1);
        arma::vec wctzc = yty.submat(0, c, c-1, c);
        
        arma::mat postcond_prec_mat = wctwc + arma::diagmat(gammacphi);
        arma::mat postcond_prec_chol = chol(postcond_prec_mat); // upper triangular form
        arma::vec postcond_mean_v = solve(trimatl(postcond_prec_chol.t()), wctzc);
        arma::vec postcondmean = solve(trimatu(postcond_prec_chol), postcond_mean_v);
        
        arma::vec randomvariates(c);
        for (int k = 0; k < c; k++) {
          randomvariates(k) = norm_rand();
        }
        arma::vec delta_w = solve(trimatu(1/std::sqrt(gamma(i-1,c)) * postcond_prec_chol),
                            randomvariates);
        
        delta.submat(i, pos_j_c(0), i, pos_j_c(pos_j_c.n_elem - 1)) = 
          arma::conv_to< arma::rowvec >::from(postcondmean + delta_w);
        
      }
    }
    
    // sampling q(q-1)/2 taus
    for (int c = 1; c < q; c++) { // c indexes the gamma value in R indexing
      IntegerVector pos_j = pos[c-1]; // vector of tau indices
      
      NumericMatrix tau_postcond_par(4, c); // posterior conditional parameters for tau_c
      
      NumericVector tau_postcond_weight(c);
      NumericVector tau_variates(c);
      
      for (int k = 0; k < c; k++) {
        int c_index = pos_j(k) - 1;
        
        tau_postcond_par(0,k) = tau_prior(0, c_index) + 1.0/alpha_d(i-1);
        tau_postcond_par(1,k) = tau_prior(1, c_index) +
          0.5/gamma(i-1,c)*std::pow(std::abs(delta(i,c_index)), alpha_d(i-1));
        tau_postcond_par(2,k) = tau_prior(2, c_index) + 1.0/alpha_d(i-1);
        tau_postcond_par(3,k) = tau_prior(3, c_index) +
          0.5/gamma(i-1,c)*std::pow(std::abs(delta(i,c_index)), alpha_d(i-1));
        
        double tau_intermediate_weight1 =
          std::pow(tau_prior(1, c_index), tau_prior(0, c_index))/std::tgamma(tau_prior(0, c_index)) *
          std::tgamma(tau_postcond_par(0,k))/std::pow(tau_postcond_par(1,k), tau_postcond_par(0,k));
        
        double tau_intermediate_weight2 =
          std::pow(tau_prior(3, c_index), tau_prior(2, c_index))/std::tgamma(tau_prior(2, c_index)) *
          std::tgamma(tau_postcond_par(2,k))/std::pow(tau_postcond_par(3,k), tau_postcond_par(2,k));
        
       tau_postcond_weight(k) =
          tau_intermediate_weight1/(tau_intermediate_weight1 + tau_intermediate_weight2);
        
        tau_variates(k) = unif_rand();
      }
      
      for (int k = 0; k < c; k++) {
        int c_index = pos_j(k) - 1;
        
        if (tau_variates(k) <= tau_postcond_weight(k)) {
          tau(i, c_index) = rgamma(1, tau_postcond_par(0,k), 1.0/tau_postcond_par(1,k))(0);
        } else {
          tau(i, c_index) = rgamma(1, tau_postcond_par(2,k), 1.0/tau_postcond_par(3,k))(0);
        }
      }
    }
    
    // sampling gamma
    for (int c = 0; c < q; c++) {
      if (c == 0) {
        // compute shape
        double shape = n/2.0 + a;
        
        //compute likelihood portion of rate
        double l2term = arma::accu(yty.submat(c, c, c, c));
        
        // compute rate
        double rate = bg + 0.5 * l2term;
        
        // sample
        gamma(i,c) = 1.0/rgamma(1, shape, 1.0/rate)(0);
      } else {
        // extract delta_c and tau_c and format
        arma::uvec temp_pos_j = pos[c-1]; // R indices
        arma::uvec pos_j = temp_pos_j - arma::ones<arma::uvec>(c); // convert to C++ indices
        arma::vec delta_c = arma::conv_to< arma::colvec >::from(delta.submat(i, pos_j(0), i, pos_j(pos_j.n_elem - 1)));
        
        arma::vec tau_c = arma::conv_to< arma::colvec >::from(tau.submat(i, pos_j(0), i, pos_j(pos_j.n_elem - 1)));
        
        // compute shape
        double shape = n/2.0 + a + c/alpha_d(i-1);
        
        // compute likelihood portion of rate
        double l2term;
        if (c > n - 1) {
          // direct calculation
          arma::vec ycm = y.col(c) - y.cols(0, c-1) * delta_c;
          l2term = dot(ycm, ycm);
        } else {
          l2term = arma::accu(yty.submat(c, c, c, c)  +
            delta_c.t() * yty.submat(0, 0, c-1, c-1) * delta_c -
            2.0 * yty.submat(c, 0, c, c - 1) * delta_c);
        }
        
        // compute delta prior portion of rate
        double delta_prior_term = arma::sum(tau_c % arma::pow(arma::abs(delta_c), alpha_d(i-1)));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * delta_prior_term;
        
        // sample
        gamma(i,c) = 1.0/rgamma(1, shape, 1.0/rate)(0);
      }
    }
    
    // sampling alpha.d
    double candidate_alpha_d = alpha_d(i-1) + alpha_d_scale * norm_rand();
    double num_alpha_d = log_post_cond_alpha_d_smn(candidate_alpha_d, arma::conv_to< arma::colvec >::from(delta.row(i)),
                                               arma::conv_to< arma::colvec >::from(tau.row(i)), 
                                               arma::conv_to< arma::colvec >::from(gamma.row(i)), k1, k2);
    double denom_alpha_d = log_post_cond_alpha_d_smn(alpha_d(i-1), arma::conv_to< arma::colvec >::from(delta.row(i)),
                                                 arma::conv_to< arma::colvec >::from(tau.row(i)), 
                                                 arma::conv_to< arma::colvec >::from(gamma.row(i)), k1, k2);
    
    double logmhratio_alpha_d = num_alpha_d - denom_alpha_d;
    double random_u_alpha_d = NA_REAL;
    if (logmhratio_alpha_d >= 0) {
      alpha_d(i) = candidate_alpha_d;
      alpha_d_accept(i-1) = true;
    } else {
      random_u_alpha_d = unif_rand();
      if (random_u_alpha_d < std::exp(logmhratio_alpha_d)) {
        alpha_d(i) = candidate_alpha_d;
        alpha_d_accept(i-1) = true;
      } else {
        alpha_d(i) = alpha_d(i-1);
        alpha_d_accept(i-1) = false;
      }
    }
    
  }
  
  timer.step("sampling complete");
  
  List mcmc = List::create(Named("delta") = delta, 
                           Named("gamma") = gamma, 
                           Named("tau") = tau, 
                           Named("alpha_d") = alpha_d);
  
  List acceptances = List::create(Named("alpha_accept_rate") = arma::mean(alpha_d_accept));
  
  NumericVector res(timer);
  for (int i = 0; i < res.size(); i++) {
    res[i] = res[i]/1e9; // convert from nanoseconds to seconds
  }
  
  return List::create(Named("mcmc") = mcmc,
                      Named("acceptances") = acceptances,
                      Named("timing") = res);
}

