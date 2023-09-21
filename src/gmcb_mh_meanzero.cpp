#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include "gmcb_help_mh.h"

// [[Rcpp::depends(RcppArmadillo)]]

// 6-1-23: This is the working version that does not require initializing values for lambda and tau

using namespace Rcpp;

// debugging output
// [[Rcpp::export]]
List gmcb_mh_meanzero_debug(const arma::mat &y, 
                   const NumericVector &d_init, const NumericVector &gamma_init, 
                   const double &alpha_d_init, const NumericVector &gamma_prior, 
                   const NumericVector &alpha_prior, 
                   const NumericMatrix &tau_prior, const int &iter, 
                   const NumericVector &d_scale,
                   const double &alpha_d_scale, 
                   const List &pos) {
  
  Timer timer;
  timer.step("start"); // timer start
  
  // for ease of use, extracting the prior parameters
  double a = gamma_prior(0);
  double bg = gamma_prior(1);
  double k1 = alpha_prior(0);
  double k2 = alpha_prior(1);
  
  // dimensions
  int n = y.n_rows;
  int q = y.n_cols;
  int lowertcount = q*(q-1)/2;
  int iter_1 = iter + 1;
  
  // to hold output
  NumericMatrix tau(iter_1, lowertcount);
  NumericMatrix delta(iter_1, lowertcount);
  NumericMatrix gamma(iter_1, q);
  NumericVector alpha_d(iter_1);
  
  // initialize values
  tau(0, _) = rep(NA_REAL, lowertcount); // R version does not specify a starting value since it is not necessary
  delta(0, _) = d_init;
  gamma(0, _) = gamma_init;
  alpha_d(0) = alpha_d_init;
  
  // acceptance indicator for MH steps
  NumericMatrix delta_accept(iter, lowertcount);
  NumericVector alpha_d_accept(iter);
  
  // for storing random variates used in debugging
  NumericMatrix tau_variates(iter, lowertcount);
  
  NumericMatrix delta_previous_state(iter, lowertcount);
  NumericMatrix delta_proposal(iter, lowertcount);
  NumericMatrix delta_log_ratio(iter, lowertcount);
  NumericMatrix delta_u(iter, lowertcount);
  
  NumericVector alpha_d_previous_state(iter);
  NumericVector alpha_d_proposal(iter);
  NumericVector alpha_d_log_ratio(iter);
  NumericVector alpha_d_u(iter);
  
  // computational aids
  arma::mat yty;
  if (q < n) {
    yty = y.t() * y;
  } else {
    arma::mat y_firstn = y.cols(0, n - 1); // only perform the calculation for the first n columns
    yty = y_firstn.t() * y_firstn; 
  }
  
  timer.step("pre-sampling preparation complete");
  
  // sampler
  for (int i = 1; i < iter_1; i++) { // start at 1, since initializing values in row 0
    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
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
          0.5/gamma(i-1,c)*std::pow(std::abs(delta(i-1,c_index)), alpha_d(i-1));
        tau_postcond_par(2,k) = tau_prior(2, c_index) + 1.0/alpha_d(i-1);
        tau_postcond_par(3,k) = tau_prior(3, c_index) + 
          0.5/gamma(i-1,c)*std::pow(std::abs(delta(i-1,c_index)), alpha_d(i-1));
        
        double tau_intermediate_weight1 = 
          std::pow(tau_prior(1, c_index), tau_prior(0, c_index))/std::tgamma(tau_prior(0, c_index)) *
          std::tgamma(tau_postcond_par(0,k))/std::pow(tau_postcond_par(1,k), tau_postcond_par(0,k));
        
        double tau_intermediate_weight2 = 
          std::pow(tau_prior(3, c_index), tau_prior(2, c_index))/std::tgamma(tau_prior(2, c_index)) *
          std::tgamma(tau_postcond_par(2,k))/std::pow(tau_postcond_par(3,k), tau_postcond_par(2,k));
        
        tau_postcond_weight(k) = 
          tau_intermediate_weight1/(tau_intermediate_weight1 + tau_intermediate_weight2);
        
        tau_variates(i-1,c_index) = unif_rand();
      }
      
      for (int k = 0; k < c; k++) {
        int c_index = pos_j(k) - 1;
        
        if (tau_variates(i-1,c_index) <= tau_postcond_weight(k)) {
          tau(i, c_index) = rgamma(1, tau_postcond_par(0,k), 1.0/tau_postcond_par(1,k))(0);
        } else {
          tau(i, c_index) = rgamma(1, tau_postcond_par(2,k), 1.0/tau_postcond_par(3,k))(0);
        }
      }
    }
    
    // sampling delta
    // construct arma version of delta(i-1,_)
    NumericVector delta_im1 = delta(i-1, _);
    
    for (int c = 1; c < q; c++) { // corresponds to 2:q in R
      IntegerVector ones(c, 1);
      IntegerVector temp_pos_j = pos[c-1];
      IntegerVector pos_j = temp_pos_j - ones;
      NumericVector delta_c = delta_im1[pos_j];
      
      for (int j = 0; j < c; j++) {
        int c_index = pos_j(j);
        double candidate_deltac = delta(i-1, c_index) + d_scale(c_index) * norm_rand();
        
        double num;
        double denom;
        if (c > n - 1) { // the c here is the C++ index
          num = delta_meanzero_logpostcond_direct(c, j, delta_c, candidate_deltac, y, 
                                                  gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
          
          denom = delta_meanzero_logpostcond_direct(c, j, delta_c, delta(i-1, c_index), y,
                                                    gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
        } else {
          num = delta_meanzero_logpostcond(c, j, delta_c, candidate_deltac, yty.submat(c, c, c, c), 
                                           yty.submat(c, 0, c, c-1), yty.submat(0, 0, c-1, c-1),
                                           gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
          
          denom = delta_meanzero_logpostcond(c, j, delta_c, delta(i-1, c_index), yty.submat(c, c, c, c), 
                                             yty.submat(c, 0, c, c-1), yty.submat(0, 0, c-1, c-1),
                                             gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
        }
        
        double logmhratio = num - denom;
        
        double random_u = NA_REAL;
        if (logmhratio >= 0) {
          delta_c(j) = candidate_deltac;
          delta_accept(i-1, c_index) = true;
        } else {
          random_u = unif_rand();
          if (random_u < std::exp(logmhratio)) {
            delta_c(j) = candidate_deltac;
            delta_accept(i-1, c_index) = true;
          } else {
            delta_c(j) = delta(i-1, c_index);
            delta_accept(i-1, c_index) = false;
          }
        }
        
        delta_log_ratio(i-1, c_index) = logmhratio;
        delta_proposal(i-1, c_index) = candidate_deltac;
        delta_previous_state(i-1, c_index) = delta(i-1, c_index);
        delta_u(i-1, c_index) = random_u;
        delta(i, c_index) = delta_c(j);
      }
    }
    
    // sampling gamma
    NumericVector delta_i = delta(i, _);
    arma::colvec delta_arma(delta_i.begin(), delta_i.size(), false);
    
    NumericVector tau_i = tau(i, _);
    
    for (int c = 0; c < q; c++) { 
      if (c == 0) {
        // compute shape
        double shape = n/2.0 + a;
        
        //compute likelihood portion of rate
        double l2term = arma::accu(yty.submat(c, c, c, c)); // always have yty for c = 0
        
        // compute rate
        double rate = bg + 0.5 * l2term;
        
        // sample
        gamma(i,c) = 1.0/rgamma(1, shape, 1.0/rate)(0);
      } else {
        // extract delta_c and tau_c and format
        IntegerVector ones(c, 1);
        IntegerVector temp_pos_j = pos[c - 1]; // remember these are R indices
        IntegerVector pos_j = temp_pos_j - ones; // convert position to C++ indices
        NumericVector delta_c = delta_i[pos_j];
        arma::colvec deltac_arma(delta_c.begin(), delta_c.size(), false);
        NumericVector tau_c = tau_i[pos_j];
        
        // compute shape
        double shape = n/2.0 + a + c/alpha_d(i-1);
        
        // compute likelihood portion of rate
        double l2term;
        if (c > n - 1) {
          // direct calculation
          arma::vec ycm = y.col(c) - y.cols(0, c-1) * deltac_arma;
          l2term = dot(ycm, ycm);
        } else {
          l2term = arma::accu(yty.submat(c, c, c, c)  +
            deltac_arma.t() * yty.submat(0, 0, c-1, c-1) * deltac_arma -
            2.0 * yty.submat(c, 0, c, c - 1) * deltac_arma);
        }
        
        // compute delta prior portion of rate
        NumericVector delta_prior_term = tau_c * Rcpp::pow(Rcpp::abs(delta_c), alpha_d(i-1));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * Rcpp::sum(delta_prior_term);
        
        // sample
        gamma(i,c) = 1.0/rgamma(1, shape, 1.0/rate)(0);
      }
    }
    
    // sampling alpha.d
    double candidate_alpha_d = alpha_d(i-1) + alpha_d_scale * norm_rand();
    double num_alpha_d = log_post_cond_alpha_d_mh(candidate_alpha_d, delta_i, 
                                               tau_i, gamma(i,_), k1, k2);
    double denom_alpha_d = log_post_cond_alpha_d_mh(alpha_d(i-1), delta_i, 
                                                 tau_i, gamma(i,_), k1, k2);
    
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
    
    alpha_d_previous_state(i-1) = alpha_d(i-1);
    alpha_d_proposal(i-1) = candidate_alpha_d;
    alpha_d_log_ratio(i-1) = logmhratio_alpha_d;
    alpha_d_u(i-1) = random_u_alpha_d;
    
  }  
  
  timer.step("sampling complete");
  
  List mcmc = List::create(Named("delta") = delta, 
                           Named("gamma") = gamma, 
                           Named("tau") = tau, 
                           Named("alpha_d") = alpha_d);
  
  List debug_out = List::create(Named("tau_variates") = tau_variates,
                                Named("delta_previous_state") = delta_previous_state, 
                                Named("delta_proposal") = delta_proposal,
                                Named("delta_log_ratio") = delta_log_ratio, 
                                Named("delta_u") = delta_u,
                                Named("alpha_d_previous_state") = alpha_d_previous_state, 
                                Named("alpha_d_proposal") = alpha_d_proposal, 
                                Named("alpha_d_log_ratio") = alpha_d_log_ratio,
                                Named("alpha_d_u") = alpha_d_u);
  
  NumericVector delta_accept_rate = Rcpp::colMeans(delta_accept);
  double alpha_accept_rate = Rcpp::mean(alpha_d_accept);
  
  List acceptances = List::create(Named("delta_accept_rate") = delta_accept_rate,
                                  Named("alpha_accept_rate") = alpha_accept_rate);
  
  NumericVector res(timer);
  for (int i = 0; i < res.size(); i++) {
    res[i] = res[i]/1e9; // convert from nanoseconds to seconds
  }
  
  return List::create(Named("mcmc") = mcmc, 
                      Named("debug_out") = debug_out,
                      Named("acceptances") = acceptances,
                      Named("timing") = res);
}

// no debugging output
// [[Rcpp::export]]
List gmcb_mh_meanzero_nodebug(const arma::mat &y, 
                     const NumericVector &d_init, const NumericVector &gamma_init, 
                     const double &alpha_d_init, const NumericVector &gamma_prior, 
                     const NumericVector &alpha_prior, 
                     const NumericMatrix &tau_prior, const int &iter, 
                     const NumericVector &d_scale,
                     const double &alpha_d_scale, 
                     const List &pos) {
  Timer timer;
  timer.step("start"); // timer start
  
  // for ease of use, extracting the prior parameters
  double a = gamma_prior(0);
  double bg = gamma_prior(1);
  double k1 = alpha_prior(0);
  double k2 = alpha_prior(1);
  
  // dimensions
  int n = y.n_rows;
  int q = y.n_cols;
  int lowertcount = q*(q-1)/2;
  int iter_1 = iter + 1;
  
  // to hold output
  NumericMatrix tau(iter_1, lowertcount);
  NumericMatrix delta(iter_1, lowertcount);
  NumericMatrix gamma(iter_1, q);
  NumericVector alpha_d(iter_1);
  
  // initialize values
  tau(0, _) = rep(NA_REAL, lowertcount); // R version does not specify a starting value since it is not necessary
  delta(0, _) = d_init;
  gamma(0, _) = gamma_init;
  alpha_d(0) = alpha_d_init;
  
  // acceptance indicator for MH steps
  NumericMatrix delta_accept(iter, lowertcount);
  NumericVector alpha_d_accept(iter);
  
  // computational aids
  arma::mat yty;
  if (q < n) {
    yty = y.t() * y;
  } else {
    arma::mat y_firstn = y.cols(0, n - 1); // only perform the calculation for the first n columns
    yty = y_firstn.t() * y_firstn; 
  }
  
  timer.step("pre-sampling preparation complete");

  // sampler
  for (int i = 1; i < iter_1; i++) { // start at 1, since initializing values in row 0
    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
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
          0.5/gamma(i-1,c)*std::pow(std::abs(delta(i-1,c_index)), alpha_d(i-1));
        tau_postcond_par(2,k) = tau_prior(2, c_index) + 1.0/alpha_d(i-1);
        tau_postcond_par(3,k) = tau_prior(3, c_index) + 
          0.5/gamma(i-1,c)*std::pow(std::abs(delta(i-1,c_index)), alpha_d(i-1));
        
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
    
    // sampling delta
    // construct arma version of delta(i-1,_)
    NumericVector delta_im1 = delta(i-1, _);
    
    for (int c = 1; c < q; c++) { // corresponds to 2:q in R
      IntegerVector ones(c, 1);
      IntegerVector temp_pos_j = pos[c-1];
      IntegerVector pos_j = temp_pos_j - ones;
      NumericVector delta_c = delta_im1[pos_j];
      
      for (int j = 0; j < c; j++) {
        int c_index = pos_j(j);
        double candidate_deltac = delta(i-1, c_index) + d_scale(c_index) * norm_rand();
        
        double num;
        double denom;
        if (c > n - 1) { // the c here is the C++ index
          num = delta_meanzero_logpostcond_direct(c, j, delta_c, candidate_deltac, y, 
                                                  gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
          
          denom = delta_meanzero_logpostcond_direct(c, j, delta_c, delta(i-1, c_index), y,
                                                    gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
        } else {
          num = delta_meanzero_logpostcond(c, j, delta_c, candidate_deltac, yty.submat(c, c, c, c), 
                                           yty.submat(c, 0, c, c-1), yty.submat(0, 0, c-1, c-1),
                                           gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
          
          denom = delta_meanzero_logpostcond(c, j, delta_c, delta(i-1, c_index), yty.submat(c, c, c, c), 
                                             yty.submat(c, 0, c, c-1), yty.submat(0, 0, c-1, c-1),
                                             gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
        }
        
        double logmhratio = num - denom;
        
        double random_u = NA_REAL;
        if (logmhratio >= 0) {
          delta_c(j) = candidate_deltac;
          delta_accept(i-1, c_index) = true;
        } else {
          random_u = unif_rand();
          if (random_u < std::exp(logmhratio)) {
            delta_c(j) = candidate_deltac;
            delta_accept(i-1, c_index) = true;
          } else {
            delta_c(j) = delta(i-1, c_index);
            delta_accept(i-1, c_index) = false;
          }
        }
        delta(i, c_index) = delta_c(j);
      }
    }
    
    // sampling gamma
    NumericVector delta_i = delta(i, _);
    arma::colvec delta_arma(delta_i.begin(), delta_i.size(), false);
    
    NumericVector tau_i = tau(i, _);
    
    for (int c = 0; c < q; c++) { 
      if (c == 0) {
        // compute shape
        double shape = n/2.0 + a;
        
        //compute likelihood portion of rate
        double l2term = arma::accu(yty.submat(c, c, c, c)); // always have yty for c = 0
        
        // compute rate
        double rate = bg + 0.5 * l2term;
        
        // sample
        gamma(i,c) = 1.0/rgamma(1, shape, 1.0/rate)(0);
      } else {
        // extract delta_c and tau_c and format
        IntegerVector ones(c, 1);
        IntegerVector temp_pos_j = pos[c - 1]; // remember these are R indices
        IntegerVector pos_j = temp_pos_j - ones; // convert position to C++ indices
        NumericVector delta_c = delta_i[pos_j];
        arma::colvec deltac_arma(delta_c.begin(), delta_c.size(), false);
        NumericVector tau_c = tau_i[pos_j];
        
        // compute shape
        double shape = n/2.0 + a + c/alpha_d(i-1);
        
        // compute likelihood portion of rate
        double l2term;
        if (c > n - 1) {
          // direct calculation
          arma::vec ycm = y.col(c) - y.cols(0, c-1) * deltac_arma;
          l2term = dot(ycm, ycm);
        } else {
          l2term = arma::accu(yty.submat(c, c, c, c)  +
            deltac_arma.t() * yty.submat(0, 0, c-1, c-1) * deltac_arma -
            2.0 * yty.submat(c, 0, c, c - 1) * deltac_arma);
        }
        
        // compute delta prior portion of rate
        NumericVector delta_prior_term = tau_c * Rcpp::pow(Rcpp::abs(delta_c), alpha_d(i-1));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * Rcpp::sum(delta_prior_term);
        
        // sample
        gamma(i,c) = 1.0/rgamma(1, shape, 1.0/rate)(0);
      }
    }
    
    // sampling alpha.d
    double candidate_alpha_d = alpha_d(i-1) + alpha_d_scale * norm_rand();
    double num_alpha_d = log_post_cond_alpha_d_mh(candidate_alpha_d, delta_i, 
                                               tau_i, gamma(i,_), k1, k2);
    double denom_alpha_d = log_post_cond_alpha_d_mh(alpha_d(i-1), delta_i, 
                                                 tau_i, gamma(i,_), k1, k2);
    
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
  
  NumericVector delta_accept_rate = Rcpp::colMeans(delta_accept);
  double alpha_accept_rate = Rcpp::mean(alpha_d_accept);
  
  List acceptances = List::create(Named("delta_accept_rate") = delta_accept_rate,
                                  Named("alpha_accept_rate") = alpha_accept_rate);
  
  NumericVector res(timer);
  for (int i = 0; i < res.size(); i++) {
    res[i] = res[i]/1e9; // convert from nanoseconds to seconds
  }
  
  return List::create(Named("mcmc") = mcmc, 
                      Named("acceptances") = acceptances,
                      Named("timing") = res);
}

