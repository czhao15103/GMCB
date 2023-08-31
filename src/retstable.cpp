// The following code is from the R package copula, modified for use with gmcb_smn
// The code from copula is licensed under GNU GPL v3.

#include <RcppArmadillo.h>

/*
 * Fast and accurate evaluation of sinc(x) := sin(x)/x, including the limit x=0
 *
 * @param x any (double precision) number
 * @return sinc(x)
 * @author Martin Maechler (2010-04-28)
 */
double sinc_MM(double x) {
  double ax = fabs(x);
  if(ax < 0.006) {
    if(x == 0.) return 1;
    double x2 = x*x;
    if(ax < 2e-4)
      return 1. - x2/6.;
    else return 1. - x2/6.*(1 - x2/20.);
  }
  /*< else */
  return sin(x)/x;
}

/*
 * Evaluation of Zolotarev's function, see Devroye (2009), to the power 1-alpha.
 * The 3-arg. version allows more precision for  alpha ~=~ 1
 *
 * @param x argument
 * @param alpha parameter in (0,1]
 * @return sin(alpha*x)^alpha * sin((1-alpha)*x)^(1-alpha) / sin(x)
 * @author Martin Maechler (2010-04-28)
 */
/*
#define _A_3(_x, _alpha_, _I_alpha)				                       \
 pow(_I_alpha* sinc_MM(_I_alpha*_x), _I_alpha) *		            \
 pow(_alpha_ * sinc_MM(_alpha_ *_x), _alpha_ ) / sinc_MM(_x)  \
 */

double A_3(double x, double alpha, double I_alpha) {
  return pow(I_alpha * sinc_MM(I_alpha * x), I_alpha) * pow(alpha * sinc_MM(alpha *x), alpha ) / sinc_MM(x);
}

// double A_(double x, double alpha) {
//   double Ialpha = 1.-alpha;
//   return A_3(x, alpha, Ialpha);
// }

/*
 * Evaluation of B(x)/B(0), see Devroye (2009).
 *
 * @param x argument
 * @param alpha parameter in (0,1]
 * @return sinc(x) / (sinc(alpha*x)^alpha * sinc((1-alpha)*x)^(1-alpha))
 * @author Martin Maechler (2010-04-28)
 */
double BdB0(double x,double alpha) {
  double Ialpha = 1.-alpha;
  double den = pow(sinc_MM(alpha*x),alpha) * pow(sinc_MM(Ialpha*x),Ialpha);
  return sinc_MM(x) / den;
}

/*
 * Sample St ~ \tilde{S}(alpha, 1, (cos(alpha*pi/2)*V_0)^{1/alpha},
 *			 V_0*I_{alpha = 1}, h*I_{alpha != 1}; 1)
 * with Laplace-Stieltjes transform exp(-V_0((h+t)^alpha-h^alpha)),
 * see Nolan's book for the parametrization, via double rejection,
 * see Devroye (2009).
 *
 * @param St vector of random variates (result) // 
 * @param V0 vector of random variates V0
 * @param h parameter in [0,infinity)
 * @param alpha parameter in (0,1]
 * @param n length of St
 * @return none
 * @author Marius Hofert, Martin Maechler
 */

double retstable_LD(double alpha, double V0, double h)
{
  /*
   * alpha == 1 => St corresponds to a point mass at V0 with Laplace-Stieltjes
   * transform exp(-V0*t)
   */
  if(alpha == 1.) {
    return V0;
    // for(R_xlen_t i = 0; i < n; i++)
    //   St[i] = V0[i];
    // return;
  }
  
  // compute variables not depending on V0
  const double c1 = sqrt(M_PI_2);
  const double c2 = 2.+c1;
  double Ialpha = 1. - alpha, // "FIXME" : allow Ialpha to be argument
    b = Ialpha/alpha,
    h_a = pow(h, alpha);
  
  // for(R_xlen_t i = 0; i < n; i++) { // for each of the n required variates
  
  /*< set lambda for our parameterization */
  double lambda_alpha = h_a*V0; // < Marius Hofert: work directly with lambda^alpha (numerically more stable for small alpha)
  
  /*
   * Apply the algorithm of Devroye (2009) to draw from
   * \tilde{S}(alpha, 1, (cos(alpha*pi/2))^{1/alpha}, I_{alpha = 1},
   * 	     lambda*I_{alpha != 1};1) with Laplace-Stieltjes transform
   * exp(-((lambda+t)^alpha-lambda^alpha))
   */
  double gamma = lambda_alpha*alpha*Ialpha;
  double sgamma = sqrt(gamma);
  double c3 = c2* sgamma;
  double xi = (1. + M_SQRT2 * c3)/M_PI; /*< according to John Lau */
  double psi = c3*exp(-gamma*M_PI*M_PI/8.)/M_SQRT_PI;
  double w1 = c1*xi/sgamma;
  double w2 = 2.*M_SQRT_PI * psi;
  double w3 = xi*M_PI;
  double X, c, E;
  do {
    double U, z, Z;
    do {
      double V = unif_rand();
      if(gamma >= 1) {
        if(V < w1/(w1+w2)) U = fabs(norm_rand())/sgamma;
        else{
          double W_ = unif_rand();
          U = M_PI*(1.-W_*W_);
        }
      }
      else{
        double W_ = unif_rand();
        if(V < w3/(w2+w3)) U = M_PI*W_;
        else U = M_PI*(1.-W_*W_);
      }
      double W = unif_rand();
      double zeta = sqrt(BdB0(U,alpha));
      z = 1/(1-pow(1+alpha*zeta/sgamma,-1/alpha)); /*< Marius Hofert: numerically more stable for small alpha */
  /*< compute rho */
  double rho = M_PI*exp(-lambda_alpha*(1. - 1./(zeta*zeta))) /
    ((1.+c1)*sgamma/zeta + z);
  double d = 0.;
  if(U >= 0 && gamma >= 1) d += xi*exp(-gamma*U*U/2.);
  if(U > 0 && U < M_PI) d += psi/sqrt(M_PI-U);
  if(U >= 0 && U <= M_PI && gamma < 1) d += xi;
  rho *= d;
  Z = W*rho;
    } while( !(U < M_PI && Z <= 1.)); /* check rejection condition */
  
  double
    // a = pow(_A_3(U,alpha,Ialpha), 1./Ialpha),
    a = pow(A_3(U,alpha,Ialpha), 1./Ialpha),
      m = pow(b/a,alpha)*lambda_alpha,
      delta = sqrt(m*alpha/a),
      a1 = delta*c1,
      a3 = z/a,
      s = a1+delta+a3;
  double V_ = unif_rand(), N_ = 0., E_ = 0.; // -Wall
  if(V_ < a1/s) {
    N_ = norm_rand();
    X = m-delta*fabs(N_);
  } else {
    if(V_ < (a1+delta)/s) /*< according to John Lau */
  X = m+delta*unif_rand();
    else {
      E_ = exp_rand();
      X = m+delta+E_*a3;
    }
  }
  E = -log(Z);
  // /*< check rejection condition */
  // c = a*(X-m)+exp((1/alpha)*log(lambda_alpha)-b*log(m))*(pow(m/X,b)-1); /*< Marius Hofert: numerically more stable for small alpha */
  
  // Christina Zhao: modify the rejection condition to allow lambda = 0
  // with the original computation, c = nan when lambda = 0 because m = 0
  c = a*(X-m);
  c += (m!=0) ? exp((1/alpha)*log(lambda_alpha)-b*log(m))*(pow(m/X,b)-1) : 0.0;
  
  if(X < m) c -= N_*N_/2.;
  else if(X > m+delta) c -= E_;
  
  } while (!(X >= 0 && c <= E));
  /*
   * Transform variates from the distribution corresponding to the
   * Laplace-Stieltjes transform exp(-((lambda+t)^alpha-lambda^alpha))
   * to those of the distribution corresponding to the Laplace-Stieltjes
   * transform exp(-V_0((h+t)^alpha-h^alpha)).
   */
  // St[i] = exp(1/alpha*log(V0[i])-b*log(X)); /*< Marius Hofert: numerically more stable for small alpha */
  
  // } /*< end for */
  return exp(1/alpha*log(V0)-b*log(X));
}