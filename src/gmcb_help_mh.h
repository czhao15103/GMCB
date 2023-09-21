#ifndef GMCB_HELP_MH
#define GMCB_HELP_MH

Rcpp::NumericVector rep_times_mh(Rcpp::NumericVector x, Rcpp::NumericVector y);

double log_post_cond_alpha_b_mh(double alpha_b, Rcpp::NumericVector b, 
                                Rcpp::NumericVector rep_gamma,
                             Rcpp::NumericVector lambda, double k1, double k2);

double log_post_cond_alpha_d_mh(double alpha_d, Rcpp::NumericVector delta, 
                                Rcpp::NumericVector tau,
                             Rcpp::NumericVector gamma, double k1, double k2);

arma::mat delta_to_matrix_inner_mh(Rcpp::NumericVector &delta, int q);
  
arma::mat sigmainv_calc_mh(Rcpp::NumericVector &delta, Rcpp::NumericVector &gamma,
                        bool cov);

double b_logpostcond(int index, const arma::mat &yty, const arma::mat &xtx, 
                     const arma::mat &ytx, Rcpp::NumericVector b,
                     double bij, Rcpp::NumericVector delta, Rcpp::NumericVector gamma,
                     double lambdaij, double alpha_b);

double b_logpostcond_direct(int index, const arma::mat &x, 
                            const arma::mat &y, Rcpp::NumericVector b,
                            double bij, Rcpp::NumericVector delta, Rcpp::NumericVector gamma,
                            double lambdaij, double alpha_b);

double delta_logpostcond(int c, int j, Rcpp::NumericVector deltac, double deltacj,
                         const arma::mat &yctyc, const arma::mat &yctycm1, 
                         const arma::mat &ycm1tycm1, const arma::mat &yctx, 
                         const arma::mat &ycm1tx, const arma::mat&xtx,
                         arma::mat &b,
                         double gammac, double alphad, double taucj);

double delta_logpostcond_direct(int c, int j, Rcpp::NumericVector deltac, double deltacj,
                                const arma::mat &y, const arma::mat &x, 
                                arma::mat &b,
                                double gammac, double alphad, double taucj);

double delta_meanzero_logpostcond(int c, int j, Rcpp::NumericVector deltac, double deltacj,
                                  const arma::mat &yctyc, const arma::mat &yctycm1, 
                                  const arma::mat &ycm1tycm1,
                                  double gammac, double alphad, double taucj);

double delta_meanzero_logpostcond_direct(int c, int j, Rcpp::NumericVector deltac, double deltacj,
                                         const arma::mat &y, 
                                         double gammac, double alphad, double taucj);

std::map<std::string,arma::mat> mh_precompute(int p, int q, int n, const arma::mat &x, 
                                              const arma::mat &y);
  
#endif