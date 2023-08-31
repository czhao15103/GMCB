#ifndef GMCB_HELP_SMN
#define GMCB_HELP_SMN

arma::vec rep_times_smn(arma::vec x, arma::vec y);

double log_post_cond_alpha_b_smn(double alpha_b, arma::vec b, arma::vec rep_gamma,
                             arma::vec lambda, double k1, double k2);

double log_post_cond_alpha_d_smn(double alpha_d, arma::vec delta, arma::vec tau,
                             arma::vec gamma, double k1, double k2);

arma::mat delta_to_matrix_inner_smn(arma::vec &delta, int q);

arma::mat sigmainv_calc_smn(arma::vec &delta, arma::vec &gamma, bool cov);

std::map<std::string,arma::mat> debugmatrices(int p, int q, int n, int iter, 
                                              const Rcpp::List &pos);

std::map<std::string,arma::mat> b_compute(int p, int q, int n, const arma::mat &x, 
                                          const arma::mat &y);
  
#endif