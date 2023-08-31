
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

// The following code is from the rbridge package, which is licensed under GPL-3 and modified
// for our use

arma::uvec SetDiff(arma::uvec & x, arma::uvec & y){
  
  for (size_t j = 0; j < y.n_elem; j++) {
    arma::uword q1 = arma::conv_to<arma::uword>::from(arma::find(x == y[j]));
    x.shed_row(q1);
  }
  return x;
}

arma::mat BetaIntial(arma::mat x, 
                     arma::colvec y, 
                     arma::colvec lambda) {
  arma::uword n_lam = lambda.n_elem;
  arma::mat U,V;
  arma::vec s;
  svd_econ(U,s,V,x);
  
  arma::mat rhs = trans(U)*y;
  arma::uword sx = s.n_elem;
  arma::vec div = repmat(pow(s,2),n_lam,1)+repelem(lambda, sx, 1);
  arma::vec a = repmat(s%rhs,n_lam,1)/div;
  arma::mat A;
  A.insert_cols(0, a);
  A.reshape(sx, n_lam);
  arma::mat coef = V*A; 
  return(coef);
}

// [[Rcpp::export]]
arma::sp_mat Bridge( const arma::mat & x,
                     const arma::colvec & y,
                     const double & q,
                     const arma::colvec & lambda,
                     const double & converge,
                     const double & eta) {
  //
  arma::uword p = x.n_cols;
  arma::uword n_lambda = lambda.n_elem;
  arma::sp_mat beta_mat(p, n_lambda);
  arma::mat beta_start = BetaIntial(x,y,lambda);
  
  
  
  
  for (arma::uword i = 0; i < n_lambda; i++){ // for look 1 begin
    
    // initial beta without intercept
    arma::vec beta_prev = beta_start.col(i);
    
    // initial converge
    double converge = 1e+10;
    
    arma::uword iteration = 0;
    // iteration until converge or two many iterations
    
    while((converge > eta) & (iteration <= 100)){
      iteration += 1;
      
      arma::uvec del = find(abs(beta_prev)<eta);
      
      
      double l_del   = del.n_elem;
      double l_beta_prev   = beta_prev.n_elem;
      
      
      // if all coefficient are small enough then stop the iteration
      if(l_del == l_beta_prev){
        beta_prev.rows(del).zeros();
        converge = 0;
        
        
      }else if(l_del>0){
        beta_prev.rows(del).zeros();
        
        
        // update design matrix x
        arma::mat x_new = x;
        arma::uvec ind_p = linspace<uvec>(0, x.n_cols-1, x.n_cols);
        arma::uvec set_diff_del = SetDiff(ind_p, del);
        double l_set_diff_del   = set_diff_del.n_elem;
        x_new = x_new.cols(set_diff_del);
        
        
        
        // calculate the diagonal matrix involving penalty terms
        arma::mat diag_mat = zeros(l_set_diff_del, l_set_diff_del);
        if(l_beta_prev-l_del==1){
          diag_mat = lambda(i)*q*pow(abs(beta_prev(set_diff_del)),q-2)/2;
        }else {
          diag_mat = diagmat(lambda(i)*q*pow(abs(beta_prev(set_diff_del)),q-2)/2);
        }
        
        
        // next beta
        arma::vec beta_curr = beta_prev;
        arma::mat  C_new    =  trans(x_new)*x_new;
        beta_curr(set_diff_del) = inv(C_new+diag_mat)*trans(x_new)*y;
        
        // new converge value
        converge  = accu(square(beta_curr-beta_prev));
        
        
        // next iteration
        beta_prev = beta_curr;
        
      }else{
        arma::mat x_new = x;
        arma::mat C_new = trans(x_new)*x_new;
        arma::mat diag_mat = diagmat(lambda(i)*q*(pow(abs(beta_prev),q-2)/2));
        arma::vec beta_curr = inv(C_new+diag_mat)*trans(x_new)*y;
        
        
        // new converge value
        converge  = accu(square(beta_curr-beta_prev));
        
        // next iteration
        beta_prev = beta_curr;
        
      }
    }
    
    beta_mat.col(i) = beta_prev;
  }
  
  return(beta_mat);
  
}// for look 1 end

// [[Rcpp::export]]
arma::mat Prediction_Grid(const arma::mat & x_test,
                          const arma::mat & x_train,
                          const arma::vec & y_train,
                          const arma::mat & grid_betas){
  // Function that returns predictions from a sequence of betas (coefficients)
  arma::uword n = x_test.n_rows;
  arma::uword len_grid = grid_betas.n_cols;
  arma::mat predictions = zeros(n, len_grid);
  arma::rowvec mu_x = mean(x_train);
  double mu_y = arma::mean(y_train);
  for(arma::uword i = 0; i < len_grid; i++){
    predictions.col(i) =  mu_y + x_test * grid_betas.col(i);
    predictions.col(i).each_row() -= mu_x * grid_betas.col(i);
  }
  return(predictions);
}
