#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>
#include "retstable.h"
#include "gmcb_help_smn.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

// with debugging output
// [[Rcpp::export]]
List gmcb_smn_debug(const arma::mat &y, const arma::mat &x, const arma::rowvec &b_init,
                    const arma::rowvec &d_init, const arma::rowvec &gamma_init, 
                    const double &alpha_b_init, 
                    const double &alpha_d_init, const arma::rowvec &lambda_init, 
                    const arma::rowvec &tau_init,
                    const arma::vec &gamma_prior, 
                    const arma::vec &alpha_prior, const arma::mat &lambda_prior,
                    const arma::mat &tau_prior, const int &iter, 
                    const double &alpha_b_scale, const double &alpha_d_scale, 
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
  int p = x.n_cols;
  int pq = p*q;
  int lowertcount = q*(q-1)/2;
  int iter_1 = iter + 1;
  
  // Rcout << "dimensions defined" << std::endl;
  
  // to hold output
  arma::mat lambda(iter_1, pq);
  arma::mat b(iter_1, pq);
  arma::vec alpha_b(iter_1);
  arma::mat tau(iter_1, lowertcount);
  arma::mat delta(iter_1, lowertcount);
  arma::mat gamma(iter_1, q);
  arma::vec alpha_d(iter_1);
  
  // Rcout << "output matrices created" << std::endl;
  
  // initialize values
  lambda.row(0) = lambda_init; 
  b.row(0) = b_init;
  alpha_b(0) = alpha_b_init;
  tau.row(0) = tau_init; 
  delta.row(0) = d_init;
  gamma.row(0) = gamma_init;
  alpha_d(0) = alpha_d_init;
  
  // Rcout << "initialization values saved" << std::endl;
  
  // acceptance indicator for MH steps
  arma::vec alpha_b_accept(iter);
  arma::vec alpha_d_accept(iter);
  
  // debugging output, produced based on size of p and q relative to n
  std::map<std::string,arma::mat> debug = debugmatrices(p, q, n, iter, pos);
  
  // computational aids
  arma::uvec zero_to_pm1(p); // for sampling gamma
  for (int j = 0; j < p; j++) {
    zero_to_pm1(j) = j;
  }
  
  // for sampling B
  std::map<std::string,arma::mat> bsamp = b_compute(p, q, n, x, y);
  
  timer.step("pre-sampling preparation complete");
  
  // sampler
  for (int i = 1; i < iter_1; i++) { // start at 1, since initializing values in row 0

    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    // arrange gamma values for sampling other parameters
    arma::rowvec rep_gamma_row = arma::repelem(gamma.row(i-1), 1, p); // repeat each element of the row vector p times
    arma::vec rep_gamma = arma::conv_to< arma::colvec >::from(rep_gamma_row);
    
    // sample omega, latent variables for B
    double half_alpha_b = alpha_b(i-1)/2.0;
    
    // NumericVector omega_tilt = Rcpp::pow(lambda(i-1,_)/(2*rep_gamma), 2/alpha_b(i-1)) * Rcpp::pow(b(i-1,_), 2.0);
    arma::vec omega(pq);
    for (int j = 0; j < pq; j++) {
      double omega_tilt = std::pow(lambda(i-1,j)/(2.0 * rep_gamma(j)), 2.0/alpha_b(i-1)) * std::pow(b(i-1,j), 2.0)/2.0;
      
      omega(j) = retstable_LD(half_alpha_b, std::pow(2.0, half_alpha_b), omega_tilt);
    }
    
    debug["omega_variates"].row(i-1) = arma::conv_to< arma::rowvec >::from(omega);
    
    // B sampling
    if (p >= n) {
      // preparation for sampling B following Bhattacharya et al. (2016)
      arma::vec delta_im1 = arma::conv_to< arma::colvec >::from(delta.row(i-1));
      arma::mat t_matrix = delta_to_matrix_inner_smn(delta_im1, q); // T
      arma::mat t_matrix_inv = arma::inv(trimatl(t_matrix)); // T^(-1)
      
      arma::vec gamma_im1 = arma::conv_to< arma::colvec >::from(gamma.row(i-1));
      arma::mat d_matrix_inv_half = arma::diagmat(sqrt(1/gamma_im1)); // D^(-1/2)
      
      arma::mat chol_prec = t_matrix.t() * d_matrix_inv_half;
      arma::mat t_kronecker_vt = kron(t_matrix, bsamp["x_vt"]);
      
      arma::vec lambda_im1 = arma::conv_to< arma::colvec >::from(lambda.row(i-1));
      arma::vec Delta_inv_diag = arma::pow((2*rep_gamma)/lambda_im1, 2/alpha_b(i-1)) / omega;
      arma::mat Delta_inv = arma::diagmat(Delta_inv_diag);
      
      arma::mat phi = kron(d_matrix_inv_half, bsamp["x_c"]);
      arma::vec alpha_bhat = vectorise(bsamp["uty"] * chol_prec);
      arma::mat t_kronecker_vt_kronecker_phi = kron(chol_prec, bsamp["vct"]);
      arma::mat b_w_lhs = t_kronecker_vt_kronecker_phi.t() * Delta_inv * t_kronecker_vt_kronecker_phi;
      
      // sampling B
      arma::vec eta_delta(n*q);
      for (int j = 0; j < (n*q); j++) {
        eta_delta(j) = norm_rand();
      }
      arma::vec eta_u_int_sd = arma::sqrt(Delta_inv_diag);
      
      arma::vec eta_u_rv(pq);
      for (int j = 0; j < pq; j++) {
        eta_u_rv(j) = norm_rand();
      }
      
      arma::vec eta_u_int = eta_u_int_sd % eta_u_rv;
      arma::vec eta_u = t_kronecker_vt * eta_u_int;
      arma::mat eta_nu = phi * eta_u + eta_delta;
      arma::mat inq(n*q, n*q, arma::fill::eye);
      arma::mat eta_w = solve(b_w_lhs + inq, alpha_bhat - eta_nu, arma::solve_opts::likely_sympd);
      
      arma::vec eta_vec = eta_u + t_kronecker_vt * Delta_inv* t_kronecker_vt_kronecker_phi * eta_w;
      arma::mat eta_mat = reshape(eta_vec, p, q);
      
      arma::vec sample_b = vectorise(bsamp["x_v"] * eta_mat * t_matrix_inv.t());
      b.row(i) = arma::conv_to< arma::rowvec >::from(sample_b);
      
      debug["b_u"].row(i-1) = arma::conv_to< arma::rowvec >::from(eta_u);
      debug["b_z"].row(i-1) = arma::conv_to< arma::rowvec >::from(eta_delta);
      debug["b_eta"].row(i-1) = arma::conv_to< arma::rowvec >::from(eta_vec);
      
    } else {
      // sampling B directly
      arma::vec delta_im1 = arma::conv_to< arma::colvec >::from(delta.row(i-1));
      arma::vec gamma_im1 = arma::conv_to< arma::colvec >::from(gamma.row(i-1));
      arma::mat sigmainv = sigmainv_calc_smn(delta_im1, gamma_im1, false);
      arma::mat thetainv = kron(sigmainv, bsamp["xtx"]);
      
      arma::vec om_diag = omega % arma::pow(arma::conv_to< arma::colvec >::from(lambda.row(i-1)) / (2 * rep_gamma), 1/half_alpha_b);
      arma::mat om_mat = diagmat(om_diag);
      arma::mat phi = thetainv + om_mat;

      arma::mat phi_chol = arma::chol(phi); // upper triangular form
      arma::vec postmean_v = solve(trimatl(phi_chol.t()), thetainv * bsamp["b_ols_vec"]);
      arma::vec postmean = solve(trimatu(phi_chol), postmean_v);
      
      arma::vec b_randomvariates(pq);
      for (int j = 0; j < pq; j++) {
        b_randomvariates(j) = norm_rand();
      }
      arma::vec b_w = solve(trimatu(phi_chol), b_randomvariates);
      b.row(i) = arma::conv_to< arma::rowvec >::from(postmean + b_w);
      debug["b_z"].row(i-1) = arma::conv_to< arma::rowvec >::from(b_randomvariates);
    }
    
    // sample lambda first
    arma::mat lambda_postcond_par(4, pq);
    NumericVector lambda_postcond_weight(pq);
    arma::rowvec lambda_sample_aug(pq);
    for (int j = 0; j < pq; j++) {
      lambda_postcond_par(0, j) = lambda_prior(0, j) + 1.0/alpha_b(i-1);
      lambda_postcond_par(1, j) = lambda_prior(1, j) + 0.5*std::pow(std::abs(b(i, j)), alpha_b(i-1))/rep_gamma(j);
      lambda_postcond_par(2, j) = lambda_prior(2, j) + 1.0/alpha_b(i-1);
      lambda_postcond_par(3, j) = lambda_prior(3, j) + 0.5*std::pow(std::abs(b(i, j)), alpha_b(i-1))/rep_gamma(j);
      
      double intermediate_weight1 = std::pow(lambda_prior(1,j), lambda_prior(0, j))/std::tgamma(lambda_prior(0, j)) * 
        std::tgamma(lambda_postcond_par(0,j))/std::pow(lambda_postcond_par(1,j), lambda_postcond_par(0,j));
      
      double intermediate_weight2 = std::pow(lambda_prior(3,j), lambda_prior(2, j))/std::tgamma(lambda_prior(2, j)) * 
        std::tgamma(lambda_postcond_par(2, j))/std::pow(lambda_postcond_par(3, j), lambda_postcond_par(2, j));
      
      lambda_postcond_weight(j) = intermediate_weight1/(intermediate_weight1 + intermediate_weight2);
      
      lambda_sample_aug(j) = unif_rand();
    }
    debug["lambda_variates"].row(i-1) = lambda_sample_aug;
    
    for (int j = 0; j < pq; j++) {
      if (lambda_sample_aug(j) <= lambda_postcond_weight(j)) {  // one row
        lambda(i,j) = rgamma(1, lambda_postcond_par(0, j), 1.0/lambda_postcond_par(1, j))(0); // uses scale parameterization
      } else {
        lambda(i,j) = rgamma(1, lambda_postcond_par(2, j), 1.0/lambda_postcond_par(3, j))(0); // uses scale parameterization
      }
    }
    
    // sample alpha.b
    double candidate_alpha_b = alpha_b(i-1) + alpha_b_scale * norm_rand();
    double num_alpha_b = log_post_cond_alpha_b_smn(candidate_alpha_b, arma::conv_to< arma::colvec >::from(b.row(i)), 
                                               rep_gamma, arma::conv_to< arma::colvec >::from(lambda.row(i)), k1, k2);
    double denom_alpha_b = log_post_cond_alpha_b_smn(alpha_b(i-1), arma::conv_to< arma::colvec >::from(b.row(i)), 
                                                 rep_gamma, arma::conv_to< arma::colvec >::from(lambda.row(i)), k1, k2);
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
    
    debug["alpha_b_previous_state"].row(i-1) = alpha_b(i-1);
    debug["alpha_b_proposal"].row(i-1) = candidate_alpha_b;
    debug["alpha_b_log_ratio"].row(i-1) = logmhratio_alpha_b;
    debug["alpha_b_u"].row(i-1) = random_u_alpha_b;
    
    // construct matrix version of B for computation
    arma::mat b_mat = reshape(b.row(i), p, q);
    
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
        
        arma::mat wc_gammac = 1/sqrt(gamma(i-1,c)) * (y.cols(0, c-1) - x * b_mat.cols(0, c-1));
        arma::vec zc_gammac = 1/sqrt(gamma(i-1,c)) * (y.col(c) - x * b_mat.col(c));
        
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
        arma::mat phiinv_mat = diagmat(phiinv);
        arma::mat w = solve(wc_gammac * phiinv_mat * wc_gammac.t() + i_n, zc_gammac - v, arma::solve_opts::likely_sympd);
        
        delta.submat(i, pos_j_c(0), i, pos_j_c(pos_j_c.n_elem - 1)) = 
          arma::conv_to< arma::rowvec >::from(u + phiinv_mat * wc_gammac.t() * w);
        
        // indexing for saving the debugging random variates
        arma::uvec npos1 = pos[n-1]; // to re-index, find the element of list pos that corresponds to n + 1 - this index is n in R, so n - 1 in C++
        arma::uvec npos1_c = npos1 - arma::ones<arma::uvec>(n); // indices are R indices, so subtract 1 to obtain the indices for C++
        debug["delta_u_gn"].submat(i-1, pos_j_c(0) - npos1_c(0), 
                         i-1, pos_j_c(pos_j_c.n_elem - 1) - npos1_c(0)) = arma::conv_to< arma::rowvec >::from(u);
        debug["delta_z_gn"].submat(i-1, (c - n) * n, i-1, (c - n) * n + (n-1)) = arma::conv_to< arma::rowvec >::from(z);
        
        
      } else {
        
        arma::vec gammacphi = gamma(i-1,c) * epsilon_c % arma::pow(tau_c / (2.0 * gamma(i-1,c)), 2.0/alpha_d(i-1));
        
        arma::mat gammacphi_mat = diagmat(gammacphi);
        
        arma::mat wctwc;
        arma::vec wctzc;
        if (p < n) {
          
          wctwc = bsamp["yty"].submat(0, 0, c-1, c-1) - bsamp["ytx"].rows(0, c-1) * b_mat.cols(0, c-1) -
            b_mat.cols(0, c-1).t() * bsamp["ytx"].rows(0, c-1).t() + b_mat.cols(0, c-1).t() * bsamp["xtx"] * b_mat.cols(0, c-1);
          
          wctzc = bsamp["yty"].submat(0, c, c-1, c) - bsamp["ytx"].rows(0, c-1) * b_mat.col(c) - 
            b_mat.cols(0, c-1).t() * bsamp["ytx"].row(c).t() + b_mat.cols(0, c-1).t() * bsamp["xtx"] * b_mat.col(c);
          
        } else {
          
          arma::mat wc = y.cols(0, c-1) - x * b_mat.cols(0, c-1);
          arma::vec zc = y.col(c) - x * b_mat.col(c);
          wctwc = wc.t() * wc;
          wctzc = wc.t() * zc;
          
        }
        
        arma::mat postcond_prec_mat = wctwc + gammacphi_mat;
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
      
      arma::uvec lambda_b_c = c*p + zero_to_pm1;
      arma::vec lambda_c = arma::conv_to< arma::colvec >::from(lambda.submat(i, lambda_b_c(0), i, lambda_b_c(lambda_b_c.n_elem - 1)));
      if (c == 0) {
        // compute shape
        double shape = n/2.0 + a + p/alpha_b(i);
        
        //compute likelihood portion of rate
        double l2term;
        if (p < n) {
          l2term = arma::accu(bsamp["yty"].submat(c, c, c, c) - 2.0 * bsamp["ytx"].row(0) * b_mat.col(0) + 
            b_mat.col(0).t() * bsamp["xtx"] * b_mat.col(0));
        } else {
          arma::vec y1mxb1 = y.col(0) - x * b_mat.col(0);
          l2term = dot(y1mxb1, y1mxb1);
        }
        
        // compute B prior portion of rate
        double b_prior_term = arma::sum(lambda_c % arma::pow(arma::abs(b_mat.col(c)), alpha_b(i)));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * b_prior_term;
        
        // sample
        gamma(i,c) = 1.0/rgamma(1, shape, 1.0/rate)(0);
      } else {
        // extract delta_c and tau_c and format
        arma::uvec temp_pos_j = pos[c-1]; // R indices
        arma::uvec pos_j = temp_pos_j - arma::ones<arma::uvec>(c); // convert to C++ indices
        arma::vec delta_c = arma::conv_to< arma::colvec >::from(delta.submat(i, pos_j(0), i, pos_j(pos_j.n_elem - 1)));
        
        arma::vec tau_c = arma::conv_to< arma::colvec >::from(tau.submat(i, pos_j(0), i, pos_j(pos_j.n_elem - 1)));
        
        // compute shape
        double shape = n/2.0 + a + p/alpha_b(i) + c/alpha_d(i-1);
        
        // compute likelihood portion of rate
        double l2term;
        if (p < n) {
          
          if (c > n - p - 1) {
            // direct calculation
            arma::vec ycm = y.col(c) - x * b_mat.col(c) - (y.cols(0, c-1) - x * b_mat.cols(0, c-1)) * delta_c;
            l2term = dot(ycm, ycm);
          } else {
            l2term = arma::accu(bsamp["yty"].submat(c, c, c, c) + b_mat.col(c).t() * bsamp["xtx"] * b_mat.col(c) +
              delta_c.t() * (bsamp["yty"].submat(0, 0, c-1, c-1) - 2.0 * bsamp["ytx"].rows(0, c-1) * b_mat.cols(0, c - 1) + 
              b_mat.cols(0, c - 1).t() * bsamp["xtx"] * b_mat.cols(0, c-1)) * delta_c -
              2.0 * bsamp["ytx"].row(c) * b_mat.col(c) - 2.0 * (bsamp["yty"].submat(c, 0, c, c - 1) - 
              bsamp["ytx"].row(c) * b_mat.cols(0, c - 1)) * delta_c +
              2.0 * delta_c.t() * (bsamp["ytx"].rows(0, c-1) - b_mat.cols(0, c - 1).t() * bsamp["xtx"]) * b_mat.col(c));
          }
          
        } else {
          // direct calculation
          arma::vec ycm = y.col(c) - x * b_mat.col(c) - (y.cols(0, c-1) - x * b_mat.cols(0, c-1)) * delta_c;
          l2term = dot(ycm, ycm);
        }
        
        // compute delta prior portion of rate
        double delta_prior_term = arma::sum(tau_c % arma::pow(arma::abs(delta_c), alpha_d(i-1)));
        
        // compute B prior portion of rate
        double b_prior_term = arma::sum(lambda_c % arma::pow(arma::abs(b_mat.col(c)), alpha_b(i)));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * delta_prior_term + 0.5 * b_prior_term;
        
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
  
  List mcmc = List::create(Named("b") = b, 
                           Named("delta") = delta, 
                           Named("gamma") = gamma, 
                           Named("lambda") = lambda, 
                           Named("tau") = tau, 
                           Named("alpha_b") = alpha_b, 
                           Named("alpha_d") = alpha_d);
  
  // outputting as list preserves the matrix structure 
  List debug_out;
  debug_out["omega_variates"] = debug["omega_variates"];
  debug_out["epsilon_variates"] = debug["epsilon_variates"];
  debug_out["lambda_variates"] = debug["lambda_variates"];
  debug_out["tau_variates"] = debug["tau_variates"];
  debug_out["alpha_b_previous_state"] = debug["alpha_b_previous_state"];
  debug_out["alpha_b_proposal"] = debug["alpha_b_proposal"];
  debug_out["alpha_b_log_ratio"] = debug["alpha_b_log_ratio"];
  debug_out["alpha_b_u"] = debug["alpha_b_u"];
  debug_out["alpha_d_previous_state"] = debug["alpha_d_previous_state"];
  debug_out["alpha_d_proposal"] = debug["alpha_d_proposal"];
  debug_out["alpha_d_log_ratio"] = debug["alpha_d_log_ratio"];
  debug_out["alpha_d_u"] = debug["alpha_d_u"];
  
  // // allow check of precomputed values
  // List precomputed;
  // precomputed["yty"] = yty;
  // precomputed["xtx"] = xtx;
  // precomputed["ytx"] = ytx;

  if (p >= n) {
    debug_out["b_u"] = debug["b_u"];
    debug_out["b_z"] = debug["b_z"];
    debug_out["b_eta"] = debug["b_eta"];
    
    // precomputed["x_u"] = bsamp["x_u"];
    // precomputed["x_v"] = bsamp["x_v"];
    // precomputed["x_vt"] = bsamp["x_vt"];
    // precomputed["x_c"] = bsamp["x_c"];
    // precomputed["uty"] = bsamp["uty"];
    // precomputed["vct"] = bsamp["vct"];
  } else {
    debug_out["b_z"] = debug["b_z"];
    
    // precomputed["xtxinv"] = bsamp["xtxinv"];
    // precomputed["b_ols_vec"] = bsamp["b_ols_vec"];
  }

  if (q > n) {
    debug_out["delta_u_gn"] = debug["delta_u_gn"];
    debug_out["delta_z_gn"] = debug["delta_z_gn"];
    debug_out["delta_z_len"] = debug["delta_z_len"];
  } else {
    debug_out["delta_z"] = debug["delta_z"];
  }
  
  List acceptances = List::create(Named("alpha_b_accept_rate") = arma::mean(alpha_b_accept),
                                  Named("alpha_d_accept_rate") = arma::mean(alpha_d_accept));
  
  NumericVector res(timer);
  for (int i = 0; i < res.size(); i++) {
    res[i] = res[i]/1e9; // convert from nanoseconds to seconds
  }
  
  return List::create(Named("mcmc") = mcmc, 
                      Named("debug_out") = debug_out,
                      Named("acceptances") = acceptances,
                      // Named("precomputed") = precomputed,
                      Named("timing") = res);
}

// no debug output
// [[Rcpp::export]]
List gmcb_smn_nodebug(const arma::mat &y, const arma::mat &x, const arma::rowvec &b_init,
                      const arma::rowvec &d_init, const arma::rowvec &gamma_init, 
                      const double &alpha_b_init, 
                      const double &alpha_d_init, const arma::rowvec &lambda_init, 
                      const arma::rowvec &tau_init,
                      const arma::vec &gamma_prior, 
                      const arma::vec &alpha_prior, const arma::mat &lambda_prior,
                      const arma::mat &tau_prior, const int &iter, 
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
  arma::mat lambda(iter_1, pq);
  arma::mat b(iter_1, pq);
  arma::vec alpha_b(iter_1);
  arma::mat tau(iter_1, lowertcount);
  arma::mat delta(iter_1, lowertcount);
  arma::mat gamma(iter_1, q);
  arma::vec alpha_d(iter_1);
  
  // initialize values
  lambda.row(0) = lambda_init; 
  b.row(0) = b_init;
  alpha_b(0) = alpha_b_init;
  tau.row(0) = tau_init; 
  delta.row(0) = d_init;
  gamma.row(0) = gamma_init;
  alpha_d(0) = alpha_d_init;
  
  // acceptance indicator for MH steps
  arma::vec alpha_b_accept(iter);
  arma::vec alpha_d_accept(iter);
  
  // computational aids
  arma::uvec zero_to_pm1(p); // for sampling gamma
  for (int j = 0; j < p; j++) {
    zero_to_pm1(j) = j;
  }
  
  // for sampling B
  std::map<std::string,arma::mat> bsamp = b_compute(p, q, n, x, y);
  
  timer.step("pre-sampling preparation complete");
  
  // sampler
  for (int i = 1; i < iter_1; i++) { // start at 1, since initializing values in row 0
    if (i % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
    // arrange gamma values for sampling other parameters
    arma::rowvec rep_gamma_row = repelem(gamma.row(i-1), 1, p);
    arma::vec rep_gamma = arma::conv_to< arma::colvec >::from(rep_gamma_row);
    
    // sample omega, latent variables for B
    double half_alpha_b = alpha_b(i-1)/2.0;
    
    // NumericVector omega_tilt = Rcpp::pow(lambda(i-1,_)/(2*rep_gamma), 2/alpha_b(i-1)) * Rcpp::pow(b(i-1,_), 2.0);
    arma::vec omega(pq);
    for (int j = 0; j < pq; j++) {
      double omega_tilt = std::pow(lambda(i-1,j)/(2.0 * rep_gamma(j)), 2.0/alpha_b(i-1)) * std::pow(b(i-1,j), 2.0)/2.0;
      
      omega(j) = retstable_LD(half_alpha_b, std::pow(2.0, half_alpha_b), omega_tilt);
    }
    
    // B sampling
    if (p >= n) {
      // preparation for sampling B following Bhattacharya et al. (2016)
      arma::vec delta_im1 = arma::conv_to< arma::colvec >::from(delta.row(i-1));
      arma::mat t_matrix = delta_to_matrix_inner_smn(delta_im1, q); // T
      arma::mat t_matrix_inv = arma::inv(trimatl(t_matrix)); // T^(-1)
      
      arma::vec gamma_im1 = arma::conv_to< arma::colvec >::from(gamma.row(i-1));
      arma::mat d_matrix_inv_half = arma::diagmat(sqrt(1/gamma_im1)); // D^(-1/2)
      
      arma::mat chol_prec = t_matrix.t() * d_matrix_inv_half;
      arma::mat t_kronecker_vt = kron(t_matrix, bsamp["x_vt"]);
      
      arma::vec lambda_im1 = arma::conv_to< arma::colvec >::from(lambda.row(i-1));
      arma::vec Delta_inv_diag = arma::pow((2*rep_gamma)/lambda_im1, 2/alpha_b(i-1)) / omega;
      arma::mat Delta_inv = arma::diagmat(Delta_inv_diag);
      
      arma::mat phi = kron(d_matrix_inv_half, bsamp["x_c"]);
      arma::vec alpha_bhat = vectorise(bsamp["uty"] * chol_prec);
      arma::mat t_kronecker_vt_kronecker_phi = kron(chol_prec, bsamp["vct"]);
      arma::mat b_w_lhs = t_kronecker_vt_kronecker_phi.t() * Delta_inv * t_kronecker_vt_kronecker_phi;
      
      // sampling B
      arma::vec eta_delta(n*q);
      for (int j = 0; j < (n*q); j++) {
        eta_delta(j) = norm_rand();
      }
      arma::vec eta_u_int_sd = arma::sqrt(Delta_inv_diag);
      
      arma::vec eta_u_rv(pq);
      for (int j = 0; j < pq; j++) {
        eta_u_rv(j) = norm_rand();
      }
      
      arma::vec eta_u_int = eta_u_int_sd % eta_u_rv;
      arma::vec eta_u = t_kronecker_vt * eta_u_int;
      arma::mat eta_nu = phi * eta_u + eta_delta;
      arma::mat inq(n*q, n*q, arma::fill::eye);
      arma::mat eta_w = solve(b_w_lhs + inq, alpha_bhat - eta_nu, arma::solve_opts::likely_sympd);
      
      arma::vec eta_vec = eta_u + t_kronecker_vt * Delta_inv* t_kronecker_vt_kronecker_phi * eta_w;
      arma::mat eta_mat = reshape(eta_vec, p, q);
      
      arma::vec sample_b = vectorise(bsamp["x_v"] * eta_mat * t_matrix_inv.t());
      b.row(i) = arma::conv_to< arma::rowvec >::from(sample_b);
      
    } else {
      // sampling B directly
      arma::vec delta_im1 = arma::conv_to< arma::colvec >::from(delta.row(i-1));
      arma::vec gamma_im1 = arma::conv_to< arma::colvec >::from(gamma.row(i-1));
      arma::mat sigmainv = sigmainv_calc_smn(delta_im1, gamma_im1, false);
      arma::mat thetainv = kron(sigmainv, bsamp["xtx"]);
      
      arma::vec om_diag = omega % arma::pow(arma::conv_to< arma::colvec >::from(lambda.row(i-1)) / (2 * rep_gamma), 1/half_alpha_b);
      arma::mat om_mat = diagmat(om_diag);
      arma::mat phi = thetainv + om_mat;
      
      arma::mat phi_chol = arma::chol(phi); // upper triangular form
      arma::vec postmean_v = solve(trimatl(phi_chol.t()), thetainv * bsamp["b_ols_vec"]);
      arma::vec postmean = solve(trimatu(phi_chol), postmean_v);
      
      arma::vec b_randomvariates(pq);
      for (int j = 0; j < pq; j++) {
        b_randomvariates(j) = norm_rand();
      }
      arma::vec b_w = solve(trimatu(phi_chol), b_randomvariates);
      b.row(i) = arma::conv_to< arma::rowvec >::from(postmean + b_w);
    }
    
    // sample lambda first
    arma::mat lambda_postcond_par(4, pq);
    NumericVector lambda_postcond_weight(pq);
    arma::rowvec lambda_sample_aug(pq);
    for (int j = 0; j < pq; j++) {
      lambda_postcond_par(0, j) = lambda_prior(0, j) + 1.0/alpha_b(i-1);
      lambda_postcond_par(1, j) = lambda_prior(1, j) + 0.5*std::pow(std::abs(b(i, j)), alpha_b(i-1))/rep_gamma(j);
      lambda_postcond_par(2, j) = lambda_prior(2, j) + 1.0/alpha_b(i-1);
      lambda_postcond_par(3, j) = lambda_prior(3, j) + 0.5*std::pow(std::abs(b(i, j)), alpha_b(i-1))/rep_gamma(j);
      
      double intermediate_weight1 = std::pow(lambda_prior(1,j), lambda_prior(0, j))/std::tgamma(lambda_prior(0, j)) * 
        std::tgamma(lambda_postcond_par(0,j))/std::pow(lambda_postcond_par(1,j), lambda_postcond_par(0,j));
      
      double intermediate_weight2 = std::pow(lambda_prior(3,j), lambda_prior(2, j))/std::tgamma(lambda_prior(2, j)) * 
        std::tgamma(lambda_postcond_par(2, j))/std::pow(lambda_postcond_par(3, j), lambda_postcond_par(2, j));
      
      lambda_postcond_weight(j) = intermediate_weight1/(intermediate_weight1 + intermediate_weight2);
      
      lambda_sample_aug(j) = unif_rand();
    }
    
    for (int j = 0; j < pq; j++) {
      if (lambda_sample_aug(j) <= lambda_postcond_weight(j)) {  // one row
        lambda(i,j) = rgamma(1, lambda_postcond_par(0, j), 1.0/lambda_postcond_par(1, j))(0); // uses scale parameterization
      } else {
        lambda(i,j) = rgamma(1, lambda_postcond_par(2, j), 1.0/lambda_postcond_par(3, j))(0); // uses scale parameterization
      }
    }
    
    // sample alpha.b
    double candidate_alpha_b = alpha_b(i-1) + alpha_b_scale * norm_rand();
    double num_alpha_b = log_post_cond_alpha_b_smn(candidate_alpha_b, arma::conv_to< arma::colvec >::from(b.row(i)), 
                                               rep_gamma, arma::conv_to< arma::colvec >::from(lambda.row(i)), k1, k2);
    double denom_alpha_b = log_post_cond_alpha_b_smn(alpha_b(i-1), arma::conv_to< arma::colvec >::from(b.row(i)), 
                                                 rep_gamma, arma::conv_to< arma::colvec >::from(lambda.row(i)), k1, k2);
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
    
    // construct matrix version of B for computation
    arma::mat b_mat = reshape(b.row(i), p, q);
    
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
        
        arma::mat wc_gammac = 1/sqrt(gamma(i-1,c)) * (y.cols(0, c-1) - x * b_mat.cols(0, c-1));
        arma::vec zc_gammac = 1/sqrt(gamma(i-1,c)) * (y.col(c) - x * b_mat.col(c));
        
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
        arma::mat phiinv_mat = diagmat(phiinv);
        arma::mat w = solve(wc_gammac * phiinv_mat * wc_gammac.t() + i_n, zc_gammac - v, arma::solve_opts::likely_sympd);
        
        delta.submat(i, pos_j_c(0), i, pos_j_c(pos_j_c.n_elem - 1)) = 
          arma::conv_to< arma::rowvec >::from(u + phiinv_mat * wc_gammac.t() * w);
        
      } else {
        
        arma::vec gammacphi = gamma(i-1,c) * epsilon_c % arma::pow(tau_c / (2.0 * gamma(i-1,c)), 2.0/alpha_d(i-1));
        
        arma::mat gammacphi_mat = diagmat(gammacphi);
        
        arma::mat wctwc;
        arma::vec wctzc;
        if (p < n) {
          
          wctwc = bsamp["yty"].submat(0, 0, c-1, c-1) - bsamp["ytx"].rows(0, c-1) * b_mat.cols(0, c-1) -
            b_mat.cols(0, c-1).t() * bsamp["ytx"].rows(0, c-1).t() + b_mat.cols(0, c-1).t() * bsamp["xtx"] * b_mat.cols(0, c-1);
          
          wctzc = bsamp["yty"].submat(0, c, c-1, c) - bsamp["ytx"].rows(0, c-1) * b_mat.col(c) - 
            b_mat.cols(0, c-1).t() * bsamp["ytx"].row(c).t() + b_mat.cols(0, c-1).t() * bsamp["xtx"] * b_mat.col(c);
          
        } else {
          
          arma::mat wc = y.cols(0, c-1) - x * b_mat.cols(0, c-1);
          arma::vec zc = y.col(c) - x * b_mat.col(c);
          wctwc = wc.t() * wc;
          wctzc = wc.t() * zc;
          
        }
        
        arma::mat postcond_prec_mat = wctwc + gammacphi_mat;
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
      arma::uvec lambda_b_c = c*p + zero_to_pm1;
      arma::vec lambda_c = arma::conv_to< arma::colvec >::from(lambda.submat(i, lambda_b_c(0), i, lambda_b_c(lambda_b_c.n_elem - 1)));
      
      if (c == 0) {
        // compute shape
        double shape = n/2.0 + a + p/alpha_b(i);
        
        //compute likelihood portion of rate
        double l2term;
        if (p < n) {
          l2term = arma::accu(bsamp["yty"].submat(c, c, c, c) - 2.0 * bsamp["ytx"].row(0) * b_mat.col(0) + 
            b_mat.col(0).t() * bsamp["xtx"] * b_mat.col(0));
        } else {
          arma::vec y1mxb1 = y.col(0) - x * b_mat.col(0);
          l2term = dot(y1mxb1, y1mxb1);
        }
        
        // compute B prior portion of rate
        double b_prior_term = arma::accu(lambda_c % arma::pow(arma::abs(b_mat.col(c)), alpha_b(i)));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * b_prior_term;
        
        // sample
        gamma(i,c) = 1.0/rgamma(1, shape, 1.0/rate)(0);
      } else {
        // extract delta_c and tau_c and format
        arma::uvec temp_pos_j = pos[c-1]; // R indices
        arma::uvec pos_j = temp_pos_j - arma::ones<arma::uvec>(c); // convert to C++ indices
        arma::vec delta_c = arma::conv_to< arma::colvec >::from(delta.submat(i, pos_j(0), i, pos_j(pos_j.n_elem - 1)));
        
        arma::vec tau_c = arma::conv_to< arma::colvec >::from(tau.submat(i, pos_j(0), i, pos_j(pos_j.n_elem - 1)));
        
        // compute shape
        double shape = n/2.0 + a + p/alpha_b(i) + c/alpha_d(i-1);
        
        // compute likelihood portion of rate
        double l2term;
        if (p < n) {
          
          if (c > n - p - 1) {
            // direct calculation
            arma::vec ycm = y.col(c) - x * b_mat.col(c) - (y.cols(0, c-1) - x * b_mat.cols(0, c-1)) * delta_c;
            l2term = dot(ycm, ycm);
          } else {
            l2term = arma::accu(bsamp["yty"].submat(c, c, c, c) + b_mat.col(c).t() * bsamp["xtx"] * b_mat.col(c) +
              delta_c.t() * (bsamp["yty"].submat(0, 0, c-1, c-1) - 2.0 * bsamp["ytx"].rows(0, c-1) * b_mat.cols(0, c - 1) + 
              b_mat.cols(0, c - 1).t() * bsamp["xtx"] * b_mat.cols(0, c-1)) * delta_c -
              2.0 * bsamp["ytx"].row(c) * b_mat.col(c) - 2.0 * (bsamp["yty"].submat(c, 0, c, c - 1) - 
              bsamp["ytx"].row(c) * b_mat.cols(0, c - 1)) * delta_c +
              2.0 * delta_c.t() * (bsamp["ytx"].rows(0, c-1) - b_mat.cols(0, c - 1).t() * bsamp["xtx"]) * b_mat.col(c));
          }
          
        } else {
          // direct calculation
          arma::vec ycm = y.col(c) - x * b_mat.col(c) - (y.cols(0, c-1) - x * b_mat.cols(0, c-1)) * delta_c;
          l2term = dot(ycm, ycm);
        }
        
        // compute delta prior portion of rate
        arma::vec delta_prior_term = tau_c % arma::pow(arma::abs(delta_c), alpha_d(i-1));
        
        // compute B prior portion of rate
        double b_prior_term = arma::accu(lambda_c % arma::pow(arma::abs(b_mat.col(c)), alpha_b(i)));
        
        // compute rate
        double rate = bg + 0.5 * l2term + 0.5 * arma::accu(delta_prior_term) + 0.5 * b_prior_term;
        
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
  
  List mcmc = List::create(Named("b") = b, 
                           Named("delta") = delta, 
                           Named("gamma") = gamma, 
                           Named("lambda") = lambda, 
                           Named("tau") = tau, 
                           Named("alpha_b") = alpha_b, 
                           Named("alpha_d") = alpha_d);
  
  List acceptances = List::create(Named("alpha_b_accept_rate") = arma::mean(alpha_b_accept),
                                  Named("alpha_d_accept_rate") = arma::mean(alpha_d_accept));
  
  NumericVector res(timer);
  for (int i = 0; i < res.size(); i++) {
    res[i] = res[i]/1e9; // convert from nanoseconds to seconds
  }
  
  return List::create(Named("mcmc") = mcmc,
                      Named("acceptances") = acceptances,
                      Named("timing") = res);
}

