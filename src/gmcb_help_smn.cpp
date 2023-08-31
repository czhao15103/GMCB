#include <RcppArmadillo.h>
using namespace Rcpp;

// Taken from Dirk Eddelbuettel's answer here: 
// https://stackoverflow.com/questions/28442582/reproducing-r-rep-with-the-times-argument-in-c-and-rcpp
// [[Rcpp::export]]
arma::vec rep_times_smn(arma::vec x, arma::vec y) {
  int n = y.size();
  arma::vec out(sum(y));
  int ind=0;
  for (int i=0; i < n; ++i) {
    int p = y[i];
    std::fill(out.begin()+ind, out.begin()+ind+p, x[i]);
    ind += p;
  }
  return out;
}

// log posterior conditional for alpha.b
double log_post_cond_alpha_b_smn(double alpha_b, arma::vec b, arma::vec rep_gamma,
                             arma::vec lambda, double k1, double k2) {
  if (alpha_b < k1 || alpha_b > k2) {
    return -std::numeric_limits<double>::infinity();
  } else {
    int pq = b.n_elem;
    double exponent_term = 0;
    for (int j = 0; j < pq; j++) {
      exponent_term += -0.5*lambda(j) * std::pow(std::abs(b(j)), alpha_b)/rep_gamma(j);
    }
    double pq_term = pq * (std::log(alpha_b) - 1.0/alpha_b*std::log(2.0) - std::lgamma(1.0/alpha_b));
    return pq_term + 1.0/alpha_b*sum(log(lambda) - log(rep_gamma)) + exponent_term;
  }
}

// log posterior conditional for alpha.d
double log_post_cond_alpha_d_smn(double alpha_d, arma::vec delta, arma::vec tau,
                             arma::vec gamma, double k1, double k2) {
  if (alpha_d < k1 || alpha_d > k2) {
    return -std::numeric_limits<double>::infinity();
  } else {
    int q = gamma.n_elem;
    int lowertcount = q*(q-1)/2;
    arma::vec reps = arma::linspace(1, q-1, q-1);
    arma::vec calc_gamma = rep_times_smn(gamma.subvec(1, q - 1), reps);
    
    double alpha_only_term = q*(q-1)/2*(std::log(alpha_d) - 1.0/alpha_d*std::log(2.0) - std::lgamma(1.0/alpha_d));
    double gamma_term = 0;
    for (int k = 0; k < q - 1; k++) {
      gamma_term += reps(k) * log(gamma(k + 1))/ alpha_d;
    }
    double tau_term = 1.0/alpha_d * arma::sum(arma::log(tau));
    double exp_term = 0;
    for (int k = 0; k < lowertcount; k++) {
      exp_term += -0.5* tau(k)/calc_gamma(k) * std::pow(std::abs(delta(k)), alpha_d);
    }
    
    return alpha_only_term - gamma_term + tau_term + exp_term;
  }
}

// convert a q(q-1)/2 vector into a unit lower triangular matrix
// [[Rcpp::export]]
arma::mat delta_to_matrix_inner_smn(arma::vec &delta, int q) {
  arma::mat out(q, q, arma::fill::eye); // identity matrix
  arma::uvec upper_indices = trimatu_ind( arma::size(out) , 1); // don't include indices of diagonal
  out.elem(upper_indices) = -delta;
  return out.t();
}


// calculate the precision or covariance matrix from delta and gamma
arma::mat sigmainv_calc_smn(arma::vec &delta, arma::vec &gamma, bool cov) {
  int q = gamma.n_elem;
  arma::mat unit_t = delta_to_matrix_inner_smn(delta, q); //T
  
  // construct D^(-1)
  arma::mat d = diagmat(gamma);
  arma::mat dinv = diagmat(1/gamma);
  
  arma::mat out(q, q);
  if (cov) { // return covariance matrix
    arma::mat u = arma::inv(trimatl(unit_t));
    out = u * d * u.t(); 
  } else { // return precision matrix
    out = unit_t.t() * dinv * unit_t;
  }
  return out;
}

std::map<std::string,arma::mat> debugmatrices(int p, int q, int n, int iter, const List &pos) {
  // always present when debugging
  arma::mat omega_variates(iter, p*q);
  arma::mat epsilon_variates(iter, q*(q-1)/2);
  
  arma::mat lambda_variates(iter, p*q);
  arma::mat tau_variates(iter, q*(q-1)/2);
  
  arma::mat alpha_b_previous_state(iter, 1);
  arma::mat alpha_b_proposal(iter, 1);
  arma::mat alpha_b_log_ratio(iter, 1);
  arma::mat alpha_b_u(iter, 1);
  
  arma::mat alpha_d_previous_state(iter, 1);
  arma::mat alpha_d_proposal(iter, 1);
  arma::mat alpha_d_log_ratio(iter, 1);
  arma::mat alpha_d_u(iter, 1);
  
  std::map<std::string, arma::mat>m;
  m.insert(std::pair<std::string, arma::mat>("omega_variates", omega_variates));
  m.insert(std::pair<std::string, arma::mat>("epsilon_variates", epsilon_variates));
  m.insert(std::pair<std::string, arma::mat>("lambda_variates", lambda_variates));
  m.insert(std::pair<std::string, arma::mat>("tau_variates", tau_variates));
  m.insert(std::pair<std::string, arma::mat>("alpha_b_previous_state", alpha_b_previous_state));
  m.insert(std::pair<std::string, arma::mat>("alpha_b_proposal", alpha_b_proposal));
  m.insert(std::pair<std::string, arma::mat>("alpha_b_log_ratio", alpha_b_log_ratio));
  m.insert(std::pair<std::string, arma::mat>("alpha_b_u", alpha_b_u));
  m.insert(std::pair<std::string, arma::mat>("alpha_d_previous_state", alpha_d_previous_state));
  m.insert(std::pair<std::string, arma::mat>("alpha_d_proposal", alpha_d_proposal));
  m.insert(std::pair<std::string, arma::mat>("alpha_d_log_ratio", alpha_d_log_ratio));
  m.insert(std::pair<std::string, arma::mat>("alpha_d_u", alpha_d_u));
  
  if (p >= n) {
    arma::mat b_u(iter, p*q);
    arma::mat b_z(iter, n*q);
    arma::mat b_eta(iter, p*q);
    
    m.insert(std::pair<std::string, arma::mat>("b_u", b_u));
    m.insert(std::pair<std::string, arma::mat>("b_z", b_z));
    m.insert(std::pair<std::string, arma::mat>("b_eta", b_eta));
  } else {
    arma::mat b_z(iter, p*q);
    m.insert(std::pair<std::string, arma::mat>("b_z", b_z));
  }
  
  if (q > n) { 
    IntegerVector qseq = Rcpp::seq_len(q - 1) + rep(1, q - 1);
    
    LogicalVector index_greater_n(q - 1);
    for (int j = 0; j < q - 1; j++) {
      index_greater_n(j) = qseq(j) > n;
    }
    List pos_q_greater_n = pos[index_greater_n];
    List pos_c_less_n = pos[!index_greater_n];
    
    double num_indices_new_pos = 0;
    for (int j = 0; j < pos_q_greater_n.length(); j++) {
      IntegerVector posc = pos_q_greater_n[j];
      
      num_indices_new_pos += posc.length();
    }
    
    arma::mat delta_u(iter, num_indices_new_pos);
    arma::mat delta_z(iter, n*pos_q_greater_n.length());
    arma::mat delta_z_len(iter, q*(q-1)/2 - num_indices_new_pos);
    
    m.insert(std::pair<std::string, arma::mat>("delta_u_gn", delta_u));
    m.insert(std::pair<std::string, arma::mat>("delta_z_gn", delta_z));
    m.insert(std::pair<std::string, arma::mat>("delta_z_len", delta_z_len));
  } else {
    arma::mat delta_z(iter, q*(q-1)/2);
    m.insert(std::pair<std::string, arma::mat>("delta_z", delta_z));
  }
  return(m);
}

// b sampling quantities that depend on whether we have p > n
std::map<std::string,arma::mat> b_compute(int p, int q, int n, const arma::mat &x, 
                                          const arma::mat &y) {
  // modify so that xtx, yty, ytx are also computed on a case by case basis
  
  std::map<std::string, arma::mat>m;
  
  if (p < n) {
    arma::mat xtx = x.t() * x;
    arma::mat ytx = y.t() * x;
    arma::mat xtxinv = arma::inv_sympd(xtx);
    arma::mat b_ols = arma::solve(xtx, ytx.t(), arma::solve_opts::likely_sympd);
    arma::vec b_ols_vec = arma::vectorise(b_ols);
    arma::mat b_ols_vec_asmat = arma::conv_to< arma::mat >::from(b_ols_vec);
    
    m.insert(std::pair<std::string, arma::mat>("xtx", xtx));
    m.insert(std::pair<std::string, arma::mat>("ytx", ytx));
    m.insert(std::pair<std::string, arma::mat>("xtxinv", xtxinv));
    m.insert(std::pair<std::string, arma::mat>("b_ols_vec", b_ols_vec_asmat));
    
    if (q <= n) {
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
    arma::mat x_u;
    arma::mat x_v;
    arma::vec x_d;
    svd(x_u, x_d, x_v, x);
    arma::mat x_vt = x_v.t();
    arma::mat x_psi = diagmat(x_d);
    arma::mat x_c = resize(x_psi, n, p);
    
    
    // for sampling arma::vec(eta)
    arma::mat uty = x_u.t() * y;
    arma::mat vct = x_v * x_c.t();
    
    m.insert(std::pair<std::string, arma::mat>("x_u", x_u));
    m.insert(std::pair<std::string, arma::mat>("x_v", x_v));
    m.insert(std::pair<std::string, arma::mat>("x_vt", x_vt));
    m.insert(std::pair<std::string, arma::mat>("x_c", x_c));
    m.insert(std::pair<std::string, arma::mat>("uty", uty));
    m.insert(std::pair<std::string, arma::mat>("vct", vct));
    
    return(m);
  }
}

