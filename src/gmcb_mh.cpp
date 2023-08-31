#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include "gmcb_help_mh.h"

// [[Rcpp::depends(RcppArmadillo)]]

// 6-1-23: This is the working version that does not require initializing values for lambda and tau

using namespace Rcpp;

// debugging output
// [[Rcpp::export]]
List gmcb_mh_debug(const arma::mat &y, const arma::mat &x, const NumericVector &b_init, 
                   const NumericVector &d_init, const NumericVector &gamma_init, const double &alpha_b_init, 
                   const double &alpha_d_init, const NumericVector &gamma_prior, 
                   const NumericVector &alpha_prior, const NumericMatrix &lambda_prior,
                   const NumericMatrix &tau_prior, const int &iter, 
                   const NumericVector &b_scale, const NumericVector &d_scale,
                   const double &alpha_b_scale, const double &alpha_d_scale, 
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
  int p = x.n_cols;
  int pq = p*q;
  int lowertcount = q*(q-1)/2;
  int iter_1 = iter + 1;
  
  // to hold output
  NumericMatrix lambda(iter_1, pq);
  NumericMatrix b(iter_1, pq);
  NumericVector alpha_b(iter_1);
  NumericMatrix tau(iter_1, lowertcount);
  NumericMatrix delta(iter_1, lowertcount);
  NumericMatrix gamma(iter_1, q);
  NumericVector alpha_d(iter_1);
  
  // initialize values
  lambda(0, _) = rep(NA_REAL, pq); // R version does not specify a starting value since it is not necessary
  b(0, _) = b_init;
  alpha_b(0) = alpha_b_init;
  tau(0, _) = rep(NA_REAL, pq); // R version does not specify a starting value since it is not necessary
  delta(0, _) = d_init;
  gamma(0, _) = gamma_init;
  alpha_d(0) = alpha_d_init;
  
  // acceptance indicator for MH steps
  NumericMatrix b_accept(iter, pq);
  NumericVector alpha_b_accept(iter);
  NumericMatrix delta_accept(iter, lowertcount);
  NumericVector alpha_d_accept(iter);
  
  // for storing random variates used in debugging
  NumericMatrix lambda_variates(iter, pq);
  NumericMatrix tau_variates(iter, lowertcount);
  
  NumericMatrix b_previous_state(iter, pq);
  NumericMatrix b_proposal(iter, pq);
  NumericMatrix b_log_ratio(iter, pq);
  NumericMatrix b_u(iter, pq);
  
  NumericVector alpha_b_previous_state(iter);
  NumericVector alpha_b_proposal(iter);
  NumericVector alpha_b_log_ratio(iter);
  NumericVector alpha_b_u(iter);
  
  NumericMatrix delta_previous_state(iter, lowertcount);
  NumericMatrix delta_proposal(iter, lowertcount);
  NumericMatrix delta_log_ratio(iter, lowertcount);
  NumericMatrix delta_u(iter, lowertcount);
  
  NumericVector alpha_d_previous_state(iter);
  NumericVector alpha_d_proposal(iter);
  NumericVector alpha_d_log_ratio(iter);
  NumericVector alpha_d_u(iter);
  
  // computational aids
  IntegerVector zero_to_pm1 = Rcpp::seq_len(p) - Rcpp::rep(1, p);
  
  std::map<std::string,arma::mat> precomputed = mh_precompute(p, q, n, x, y);
  
  timer.step("pre-sampling preparation complete");
  
  // sampler
  for (int i = 1; i < iter_1; i++) { // start at 1, since initializing values in row 0
    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    // arrange gamma values for sampling other parameters
    NumericVector rep_gamma = rep_each(gamma(i-1, _), p);
    
    // sample lambda first
    NumericMatrix lambda_postcond_par(4, pq);
    lambda_postcond_par(0, _) = lambda_prior(0, _) + rep(1.0/alpha_b(i-1), pq);
    lambda_postcond_par(1, _) = lambda_prior(1, _) + rep(0.5, pq)*Rcpp::pow(Rcpp::abs(b(i-1, _)), alpha_b(i-1))/rep_gamma;
    lambda_postcond_par(2, _) = lambda_prior(2, _) + rep(1.0/alpha_b(i-1), pq);
    lambda_postcond_par(3, _) = lambda_prior(3, _) + rep(0.5, pq)*Rcpp::pow(Rcpp::abs(b(i-1, _)), alpha_b(i-1))/rep_gamma;
    
    NumericVector intermediate_weight1(pq);
    NumericVector intermediate_weight2(pq);
    
    for (int j = 0; j < pq; j++) {
      intermediate_weight1(j) = std::pow(lambda_prior(1,j), lambda_prior(0, j))/std::tgamma(lambda_prior(0, j)) * 
        std::tgamma(lambda_postcond_par(0,j))/std::pow(lambda_postcond_par(1,j), lambda_postcond_par(0,j));
      
      intermediate_weight2(j) = std::pow(lambda_prior(3,j), lambda_prior(2, j))/std::tgamma(lambda_prior(2, j)) * 
        std::tgamma(lambda_postcond_par(2, j))/std::pow(lambda_postcond_par(3, j), lambda_postcond_par(2, j));
    }
    
    NumericVector lambda_postcond_weight = 
      intermediate_weight1/(intermediate_weight1 + intermediate_weight2);
    
    // success is defined as drawing from the distribution associated with intermediate_weight1
    NumericVector lambda_sample_aug(pq);
    for (int j = 0; j < pq; j++) {
      lambda_sample_aug(j) = unif_rand();
    }
    
    lambda_variates(i-1, _) = lambda_sample_aug;
    
    for (int j = 0; j < pq; j++) {
      if (lambda_sample_aug(j) <= lambda_postcond_weight(j)) {  // one row
        lambda(i,j) = rgamma(1, lambda_postcond_par(0, j), 1.0/lambda_postcond_par(1, j))(0); // uses scale parameterization
      } else {
        lambda(i,j) = rgamma(1, lambda_postcond_par(2, j), 1.0/lambda_postcond_par(3, j))(0); // uses scale parameterization
      }
    }
    
    // sampling B
    NumericVector temp_b = b(i-1,_);
    
    for (int j = 0; j < pq; j++) {
      double candidatebj = b(i-1,j) + b_scale(j) * norm_rand(); 
      
      // temp_b changes after each of these calculations, although its value is correct at the end of loop
      double num;
      double denom;
      if (p < n && q < n) {
        num = b_logpostcond(j, precomputed["yty"], precomputed["xtx"], precomputed["ytx"], 
                            temp_b, candidatebj, delta(i-1,_), gamma(i-1,_),
                            lambda(i,j), alpha_b(i-1));
        
        denom = b_logpostcond(j, precomputed["yty"], precomputed["xtx"], precomputed["ytx"], 
                              temp_b, b(i-1,j), delta(i-1,_), gamma(i-1,_),
                              lambda(i,j), alpha_b(i-1));
      }  else {
        
        num = b_logpostcond_direct(j, x, y, temp_b, candidatebj, delta(i-1,_), gamma(i-1,_),
                                   lambda(i,j), alpha_b(i-1));
        
        denom = b_logpostcond_direct(j, x, y, temp_b, b(i-1,j), delta(i-1,_), gamma(i-1,_),
                                     lambda(i,j), alpha_b(i-1));
      }
      
      
      double logmhratio = num - denom;
      double random_u = NA_REAL;
      if (logmhratio >= 0) {
        temp_b(j) = candidatebj;
        b_accept(i-1,j) = true;
      } else {
        random_u = unif_rand();
        if (random_u < std::exp(logmhratio)) {
          temp_b(j) = candidatebj;
          b_accept(i-1,j) = true;
        } else {
          temp_b(j) = b(i-1,j);
          b_accept(i-1,j) = false;
        }
      }
      
      b_previous_state(i-1,j) = b(i-1,j);
      b_proposal(i-1,j) = candidatebj;
      b_log_ratio(i-1,j) = logmhratio;
      b_u(i-1,j) = random_u;
    }
    b(i,_) = temp_b;
    
    // sample alpha.b
    double candidate_alpha_b = alpha_b(i-1) + alpha_b_scale * norm_rand();  
    double num_alpha_b = log_post_cond_alpha_b_mh(candidate_alpha_b, b(i,_), rep_gamma, lambda(i,_), k1, k2);
    double denom_alpha_b = log_post_cond_alpha_b_mh(alpha_b(i-1), b(i,_), rep_gamma, lambda(i,_), k1, k2);
    double logmhratio_alpha_b = num_alpha_b - denom_alpha_b;
    double random_u_alpha_b = NA_REAL;
    if (logmhratio_alpha_b >= 0) {
      alpha_b(i) = candidate_alpha_b;
      alpha_b_accept(i-1) = true;
    } else {
      random_u_alpha_b = unif_rand();
      if (random_u_alpha_b < std::exp(logmhratio_alpha_b)) {
        alpha_b(i) = candidate_alpha_b;
        alpha_b_accept(i-1) = true;
      } else {
        alpha_b(i) = alpha_b(i-1);
        alpha_b_accept(i-1) = false;
      }
    }
    
    alpha_b_previous_state(i-1) = alpha_b(i-1);
    alpha_b_proposal(i-1) = candidate_alpha_b;
    alpha_b_log_ratio(i-1) = logmhratio_alpha_b;
    alpha_b_u(i-1) = random_u_alpha_b;
    
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
    // construct matrix version of B for computation
    NumericVector b_i = b(i, _);
    arma::vec b_arma(b_i.begin(), b_i.size(), false);
    arma::mat b_mat = reshape(b_arma, p, q);
    
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
        if (p < n) {
          
          if (c > n - p - 1) {
            // direct calculation
            
            num = delta_logpostcond_direct(c, j, delta_c, candidate_deltac, y, x, b_mat,
                                           gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
            
            denom = delta_logpostcond_direct(c, j, delta_c, delta(i-1, c_index), y, x, b_mat,
                                             gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
            
          } else {
            
            num = delta_logpostcond(c, j, delta_c, candidate_deltac, precomputed["yty"].submat(c, c, c, c), 
                                    precomputed["yty"].submat(c, 0, c, c-1), 
                                    precomputed["yty"].submat(0, 0, c-1, c-1),
                                    precomputed["ytx"].row(c), precomputed["ytx"].rows(0, c-1), 
                                    precomputed["xtx"], b_mat,
                                    gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
            
            denom = delta_logpostcond(c, j, delta_c, delta(i-1, c_index), precomputed["yty"].submat(c, c, c, c), 
                                      precomputed["yty"].submat(c, 0, c, c-1), 
                                      precomputed["yty"].submat(0, 0, c-1, c-1),
                                      precomputed["ytx"].row(c), precomputed["ytx"].rows(0, c-1), 
                                      precomputed["xtx"], b_mat,
                                      gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
            
          }
          
        } else {
          // direct calculation
          
          num = delta_logpostcond_direct(c, j, delta_c, candidate_deltac, y, x, b_mat,
                                         gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
          
          denom = delta_logpostcond_direct(c, j, delta_c, delta(i-1, c_index), y, x, b_mat,
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
    // construct matrix version of lambda for computation
    NumericVector lambda_i = lambda(i, _);
    
    NumericVector delta_i = delta(i, _);
    arma::colvec delta_arma(delta_i.begin(), delta_i.size(), false);
    
    NumericVector tau_i = tau(i, _);
    
    for (int c = 0; c < q; c++) { 
      NumericVector b_c = wrap(b_mat.col(c));
      IntegerVector lambda_b_c = c*p + zero_to_pm1;
      NumericVector lambda_c = lambda_i[lambda_b_c];
      
      if (c == 0) {
        // compute shape
        double shape = n/2.0 + a + p/alpha_b(i);
        
        //compute likelihood portion of rate
        double l2term;
        if (p < n) {
          l2term = arma::accu(precomputed["yty"].submat(c, c, c, c) - 2.0 * precomputed["ytx"].row(0) * b_mat.col(0) + 
            b_mat.col(0).t() * precomputed["xtx"] * b_mat.col(0));
        } else {
          arma::vec y1mxb1 = y.col(0) - x * b_mat.col(0);
          l2term = dot(y1mxb1, y1mxb1);
        }
        
        // compute B prior portion of rate
        double b_prior_term = Rcpp::sum(lambda_c * Rcpp::pow(Rcpp::abs(b_c), alpha_b(i)));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * b_prior_term;
        
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
        double shape = n/2.0 + a + p/alpha_b(i) + c/alpha_d(i-1);
        
        // compute likelihood portion of rate
        double l2term;
        if (p < n) {
          
          if (c > n - p - 1) {
            // direct calculation
            arma::vec ycm = y.col(c) - x * b_mat.col(c) - (y.cols(0, c-1) - x * b_mat.cols(0, c-1)) * deltac_arma;
            l2term = dot(ycm, ycm);
          } else {
            l2term = arma::accu(precomputed["yty"].submat(c, c, c, c) + b_mat.col(c).t() * precomputed["xtx"] * b_mat.col(c) +
              deltac_arma.t() * (precomputed["yty"].submat(0, 0, c-1, c-1) - 2.0 * precomputed["ytx"].rows(0, c-1) * b_mat.cols(0, c - 1) + 
              b_mat.cols(0, c - 1).t() * precomputed["xtx"] * b_mat.cols(0, c-1)) * deltac_arma -
              2.0 * precomputed["ytx"].row(c) * b_mat.col(c) - 2.0 * (precomputed["yty"].submat(c, 0, c, c - 1) - 
              precomputed["ytx"].row(c) * b_mat.cols(0, c - 1)) * deltac_arma +
              2.0 * deltac_arma.t() * (precomputed["ytx"].rows(0, c-1) - b_mat.cols(0, c - 1).t() * precomputed["xtx"]) * b_mat.col(c));
          }
          
        } else {
          // direct calculation
          arma::vec ycm = y.col(c) - x * b_mat.col(c) - (y.cols(0, c-1) - x * b_mat.cols(0, c-1)) * deltac_arma;
          l2term = dot(ycm, ycm);
        }
        
        // compute delta prior portion of rate
        NumericVector delta_prior_term = tau_c * Rcpp::pow(Rcpp::abs(delta_c), alpha_d(i-1));
        
        // compute B prior portion of rate
        double b_prior_term = Rcpp::sum(lambda_c * Rcpp::pow(Rcpp::abs(b_c), alpha_b(i)));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * Rcpp::sum(delta_prior_term) + 0.5 * b_prior_term;
        
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
  
  List mcmc = List::create(Named("b") = b, 
                           Named("delta") = delta, 
                           Named("gamma") = gamma, 
                           Named("lambda") = lambda, 
                           Named("tau") = tau, 
                           Named("alpha_b") = alpha_b, 
                           Named("alpha_d") = alpha_d);
  
  List debug_out = List::create(Named("lambda_variates") = lambda_variates, 
                                Named("tau_variates") = tau_variates,
                                Named("b_previous_state") = b_previous_state, 
                                Named("b_proposal") = b_proposal,
                                Named("b_log_ratio") = b_log_ratio, 
                                Named("b_u") = b_u, 
                                Named("alpha_b_previous_state") = alpha_b_previous_state,
                                Named("alpha_b_proposal") = alpha_b_proposal, 
                                Named("alpha_b_log_ratio") = alpha_b_log_ratio,
                                Named("alpha_b_u") = alpha_b_u,
                                Named("delta_previous_state") = delta_previous_state, 
                                Named("delta_proposal") = delta_proposal,
                                Named("delta_log_ratio") = delta_log_ratio, 
                                Named("delta_u") = delta_u,
                                Named("alpha_d_previous_state") = alpha_d_previous_state, 
                                Named("alpha_d_proposal") = alpha_d_proposal, 
                                Named("alpha_d_log_ratio") = alpha_d_log_ratio,
                                Named("alpha_d_u") = alpha_d_u);
  
  List acceptances = List::create(Named("b_accept_rate") = b_accept,
                                  Named("alpha_b_accept_rate") = alpha_b_accept,
                                  Named("delta_accept_rate") = delta_accept,
                                  Named("alpha_d_accept_rate") = alpha_d_accept);
  
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
List gmcb_mh_nodebug(const arma::mat &y, const arma::mat &x, const NumericVector &b_init, 
                     const NumericVector &d_init, const NumericVector &gamma_init, const double &alpha_b_init, 
                     const double &alpha_d_init, const NumericVector &gamma_prior, 
                     const NumericVector &alpha_prior, const NumericMatrix &lambda_prior,
                     const NumericMatrix &tau_prior, const int &iter, 
                     const NumericVector &b_scale, const NumericVector &d_scale,
                     const double &alpha_b_scale, const double &alpha_d_scale, 
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
  int p = x.n_cols;
  int pq = p*q;
  int lowertcount = q*(q-1)/2;
  int iter_1 = iter + 1;
  
  // to hold output
  NumericMatrix lambda(iter_1, pq);
  NumericMatrix b(iter_1, pq);
  NumericVector alpha_b(iter_1);
  NumericMatrix tau(iter_1, lowertcount);
  NumericMatrix delta(iter_1, lowertcount);
  NumericMatrix gamma(iter_1, q);
  NumericVector alpha_d(iter_1);
  
  // initialize values
  lambda(0, _) = rep(NA_REAL, pq); // R version does not specify a starting value since it is not necessary
  b(0, _) = b_init;
  alpha_b(0) = alpha_b_init;
  tau(0, _) = rep(NA_REAL, pq); // R version does not specify a starting value since it is not necessary
  delta(0, _) = d_init;
  gamma(0, _) = gamma_init;
  alpha_d(0) = alpha_d_init;
  
  // acceptance indicator for MH steps
  NumericMatrix b_accept(iter, pq);
  NumericVector alpha_b_accept(iter);
  NumericMatrix delta_accept(iter, lowertcount);
  NumericVector alpha_d_accept(iter);
  
  // computational aids
  IntegerVector zero_to_pm1 = Rcpp::seq_len(p) - Rcpp::rep(1, p);
  
  std::map<std::string,arma::mat> precomputed = mh_precompute(p, q, n, x, y);
  
  timer.step("pre-sampling preparation complete");

  // sampler
  for (int i = 1; i < iter_1; i++) { // start at 1, since initializing values in row 0
    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    // arrange gamma values for sampling other parameters
    NumericVector rep_gamma = rep_each(gamma(i-1, _), p);
    
    // sample lambda first
    NumericMatrix lambda_postcond_par(4, pq);
    lambda_postcond_par(0, _) = lambda_prior(0, _) + rep(1.0/alpha_b(i-1), pq);
    lambda_postcond_par(1, _) = lambda_prior(1, _) + rep(0.5, pq)*Rcpp::pow(Rcpp::abs(b(i-1, _)), alpha_b(i-1))/rep_gamma;
    lambda_postcond_par(2, _) = lambda_prior(2, _) + rep(1.0/alpha_b(i-1), pq);
    lambda_postcond_par(3, _) = lambda_prior(3, _) + rep(0.5, pq)*Rcpp::pow(Rcpp::abs(b(i-1, _)), alpha_b(i-1))/rep_gamma;
    
    NumericVector intermediate_weight1(pq);
    NumericVector intermediate_weight2(pq);
    
    for (int j = 0; j < pq; j++) {
      intermediate_weight1(j) = std::pow(lambda_prior(1,j), lambda_prior(0, j))/std::tgamma(lambda_prior(0, j)) * 
        std::tgamma(lambda_postcond_par(0,j))/std::pow(lambda_postcond_par(1,j), lambda_postcond_par(0,j));
      
      intermediate_weight2(j) = std::pow(lambda_prior(3,j), lambda_prior(2, j))/std::tgamma(lambda_prior(2, j)) * 
        std::tgamma(lambda_postcond_par(2, j))/std::pow(lambda_postcond_par(3, j), lambda_postcond_par(2, j));
    }
    
    NumericVector lambda_postcond_weight = 
      intermediate_weight1/(intermediate_weight1 + intermediate_weight2);
    
    // success is defined as drawing from the distribution associated with intermediate_weight1
    NumericVector lambda_sample_aug(pq);
    for (int j = 0; j < pq; j++) {
      lambda_sample_aug(j) = unif_rand();
    }
    
    for (int j = 0; j < pq; j++) {
      if (lambda_sample_aug(j) <= lambda_postcond_weight(j)) {  // one row
        lambda(i,j) = rgamma(1, lambda_postcond_par(0, j), 1.0/lambda_postcond_par(1, j))(0); // uses scale parameterization
      } else {
        lambda(i,j) = rgamma(1, lambda_postcond_par(2, j), 1.0/lambda_postcond_par(3, j))(0); // uses scale parameterization
      }
    }
    
    
    // sampling B
    NumericVector temp_b = b(i-1,_);
    
    for (int j = 0; j < pq; j++) {
      double candidatebj = b(i-1,j) + b_scale(j) * norm_rand();
      
      // temp_b changes after each of these calculations, although its value is correct at the end of loop
      double num;
      double denom;
      if (p < n && q < n) {
        num = b_logpostcond(j, precomputed["yty"], precomputed["xtx"], precomputed["ytx"], 
                            temp_b, candidatebj, delta(i-1,_), gamma(i-1,_),
                            lambda(i,j), alpha_b(i-1));
        
        denom = b_logpostcond(j, precomputed["yty"], precomputed["xtx"], precomputed["ytx"], 
                              temp_b, b(i-1,j), delta(i-1,_), gamma(i-1,_),
                              lambda(i,j), alpha_b(i-1));
      }  else {
        
        num = b_logpostcond_direct(j, x, y, temp_b, candidatebj, delta(i-1,_), gamma(i-1,_),
                                   lambda(i,j), alpha_b(i-1));
        
        denom = b_logpostcond_direct(j, x, y, temp_b, b(i-1,j), delta(i-1,_), gamma(i-1,_),
                                     lambda(i,j), alpha_b(i-1));
      }
      
      double logmhratio = num - denom;
      double random_u = NA_REAL;
      if (logmhratio >= 0) {
        temp_b(j) = candidatebj;
        b_accept(i-1,j) = true;
      } else {
        random_u = unif_rand();
        if (random_u < std::exp(logmhratio)) {
          temp_b(j) = candidatebj;
          b_accept(i-1,j) = true;
        } else {
          temp_b(j) = b(i-1,j);
          b_accept(i-1,j) = false;
        }
      }
    }
    b(i,_) = temp_b;
    
    // sample alpha.b
    double candidate_alpha_b = alpha_b(i-1) + alpha_b_scale * norm_rand();
    double num_alpha_b = log_post_cond_alpha_b_mh(candidate_alpha_b, b(i,_), rep_gamma, lambda(i,_), k1, k2);
    double denom_alpha_b = log_post_cond_alpha_b_mh(alpha_b(i-1), b(i,_), rep_gamma, lambda(i,_), k1, k2);
    double logmhratio_alpha_b = num_alpha_b - denom_alpha_b;
    double random_u_alpha_b = NA_REAL;
    if (logmhratio_alpha_b >= 0) {
      alpha_b(i) = candidate_alpha_b;
      alpha_b_accept(i-1) = true;
    } else {
      random_u_alpha_b = unif_rand();
      if (random_u_alpha_b < std::exp(logmhratio_alpha_b)) {
        alpha_b(i) = candidate_alpha_b;
        alpha_b_accept(i-1) = true;
      } else {
        alpha_b(i) = alpha_b(i-1);
        alpha_b_accept(i-1) = false;
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
    // construct matrix version of B for computation
    NumericVector b_i = b(i, _);
    arma::vec b_arma(b_i.begin(), b_i.size(), false);
    arma::mat b_mat = reshape(b_arma, p, q);
    
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
        if (p < n) {
          
          if (c > n - p - 1) {
            // direct calculation
            
            num = delta_logpostcond_direct(c, j, delta_c, candidate_deltac, y, x, b_mat,
                                           gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
            
            denom = delta_logpostcond_direct(c, j, delta_c, delta(i-1, c_index), y, x, b_mat,
                                             gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
            
          } else {
            
            num = delta_logpostcond(c, j, delta_c, candidate_deltac, precomputed["yty"].submat(c, c, c, c), 
                                    precomputed["yty"].submat(c, 0, c, c-1), 
                                    precomputed["yty"].submat(0, 0, c-1, c-1),
                                    precomputed["ytx"].row(c), precomputed["ytx"].rows(0, c-1), 
                                    precomputed["xtx"], b_mat,
                                    gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
            
            denom = delta_logpostcond(c, j, delta_c, delta(i-1, c_index), precomputed["yty"].submat(c, c, c, c), 
                                      precomputed["yty"].submat(c, 0, c, c-1), 
                                      precomputed["yty"].submat(0, 0, c-1, c-1),
                                      precomputed["ytx"].row(c), precomputed["ytx"].rows(0, c-1), 
                                      precomputed["xtx"], b_mat,
                                      gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
            
          }
          
        } else {
          // direct calculation
          
          num = delta_logpostcond_direct(c, j, delta_c, candidate_deltac, y, x, b_mat,
                                         gamma(i-1,c), alpha_d(i-1), tau(i, c_index));
          
          denom = delta_logpostcond_direct(c, j, delta_c, delta(i-1, c_index), y, x, b_mat,
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
    // construct matrix version of lambda for computation
    NumericVector lambda_i = lambda(i, _);
    
    NumericVector delta_i = delta(i, _);
    arma::colvec delta_arma(delta_i.begin(), delta_i.size(), false);
    
    NumericVector tau_i = tau(i, _);
    
    for (int c = 0; c < q; c++) { 
      NumericVector b_c = wrap(b_mat.col(c));
      IntegerVector lambda_b_c = c*p + zero_to_pm1;
      NumericVector lambda_c = lambda_i[lambda_b_c];
      if (c == 0) {
        // compute shape
        double shape = n/2.0 + a + p/alpha_b(i);
        
        //compute likelihood portion of rate
        double l2term;
        if (p < n) {
          l2term = arma::accu(precomputed["yty"].submat(c, c, c, c) - 2.0 * precomputed["ytx"].row(0) * b_mat.col(0) + 
            b_mat.col(0).t() * precomputed["xtx"] * b_mat.col(0));
        } else {
          arma::vec y1mxb1 = y.col(0) - x * b_mat.col(0);
          l2term = dot(y1mxb1, y1mxb1);
        }
        
        // compute B prior portion of rate
        double b_prior_term = Rcpp::sum(lambda_c * Rcpp::pow(Rcpp::abs(b_c), alpha_b(i)));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * b_prior_term;
        
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
        double shape = n/2.0 + a + p/alpha_b(i) + c/alpha_d(i-1);
        
        double l2term;
        if (p < n) {
          
          if (c > n - p - 1) {
            // direct calculation
            arma::vec ycm = y.col(c) - x * b_mat.col(c) - (y.cols(0, c-1) - x * b_mat.cols(0, c-1)) * deltac_arma;
            l2term = dot(ycm, ycm);
          } else {
            l2term = arma::accu(precomputed["yty"].submat(c, c, c, c) + b_mat.col(c).t() * precomputed["xtx"] * b_mat.col(c) +
              deltac_arma.t() * (precomputed["yty"].submat(0, 0, c-1, c-1) - 2.0 * precomputed["ytx"].rows(0, c-1) * b_mat.cols(0, c - 1) + 
              b_mat.cols(0, c - 1).t() * precomputed["xtx"] * b_mat.cols(0, c-1)) * deltac_arma -
              2.0 * precomputed["ytx"].row(c) * b_mat.col(c) - 2.0 * (precomputed["yty"].submat(c, 0, c, c - 1) - 
              precomputed["ytx"].row(c) * b_mat.cols(0, c - 1)) * deltac_arma +
              2.0 * deltac_arma.t() * (precomputed["ytx"].rows(0, c-1) - b_mat.cols(0, c - 1).t() * precomputed["xtx"]) * b_mat.col(c));
          }
          
        } else {
          // direct calculation
          arma::vec ycm = y.col(c) - x * b_mat.col(c) - (y.cols(0, c-1) - x * b_mat.cols(0, c-1)) * deltac_arma;
          l2term = dot(ycm, ycm);
        }
        
        // compute delta prior portion of rate
        NumericVector delta_prior_term = tau_c * Rcpp::pow(Rcpp::abs(delta_c), alpha_d(i-1));
        
        // compute B prior portion of rate
        double b_prior_term = Rcpp::sum(lambda_c * Rcpp::pow(Rcpp::abs(b_c), alpha_b(i)));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * Rcpp::sum(delta_prior_term) + 0.5 * b_prior_term;
        
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
  
  
  List mcmc = List::create(Named("b") = b, 
                           Named("delta") = delta, 
                           Named("gamma") = gamma, 
                           Named("lambda") = lambda, 
                           Named("tau") = tau, 
                           Named("alpha_b") = alpha_b, 
                           Named("alpha_d") = alpha_d);
  
  List acceptances = List::create(Named("b_accept_rate") = b_accept,
                                  Named("alpha_b_accept_rate") = alpha_b_accept,
                                  Named("delta_accept_rate") = delta_accept,
                                  Named("alpha_d_accept_rate") = alpha_d_accept);
  
  NumericVector res(timer);
  for (int i = 0; i < res.size(); i++) {
    res[i] = res[i]/1e9; // convert from nanoseconds to seconds
  }
  
  return List::create(Named("mcmc") = mcmc, 
                      Named("acceptances") = acceptances,
                      Named("timing") = res);
}

