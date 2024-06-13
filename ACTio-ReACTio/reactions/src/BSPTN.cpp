#if HAVE_CONFIG_H
#include <config.h>
#endif

#include <cfloat>
#include <cmath>

#include "Common.h"
#include "Quadrature.h"
#include "SpecialFunctions.h"
#include "BSPTN.h"
#include "BeyondLCDM.h"

//#include<omp.h>

#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <math.h>       /* pow */
#include <iostream>
#include <stdlib.h>
#include <functional>


using std::cref;
using std::bind;

// 1st ORDER
double F1b_k1;
double G1b_k1;
double F1b_k2;
double G1b_k2;
double F1b_k3;
double G1b_k3;

double F1b_1;
double G1b_1;

double F1b_2;
double G1b_2;

double F1b_3;
double G1b_3;


// 2nd order
double F2b_k12;
double G2b_k12;
double F2b_k13;
double G2b_k13;
double F2b_k23;
double G2b_k23;

double F2b_p2mp;
double G2b_p2mp;

double F2b_p3mp;
double G2b_p3mp;

double F2b_12a;
double G2b_12a;

double F2b_13a;
double G2b_13a;

double F2b_23a;
double G2b_23a;

// 3rd order

double F3b_12mp;
double G3b_12mp;

double F3b_13mp;
double G3b_13mp;

double F3b_21mp;
double G3b_21mp;

double F3b_23mp;
double G3b_23mp;

double F3b_31mp;
double G3b_31mp;

double F3b_32mp;
double G3b_32mp;

double F3b_1pp;
double G3b_1pp;

double F3b_2pp;
double G3b_2pp;

double F3b_3pp;
double G3b_3pp;

// 4th order
double F4b_12pp;
double G4b_12pp;
double F4b_13pp;
double G4b_13pp;
double F4b_23pp;
double G4b_23pp;



/* Euler and Continuity equations for tree level Bispectrum numerical kernels */

//alphai(k1,k2)
inline double alphai(double k1, double k2, double u1){
	return 1.+k2*u1/k1;
}

//beta(k1,k2)
inline double betai(double k1, double k2, double u1){
	return u1*(k1*k1+k2*k2+2.*k1*k2*u1)/(2.*k1*k2);
}

/* Parameters passed to system of Euler and continuity equations*/
// k (magnitude) and x (angular) values for the system of equations
// args hold wave vector magnitudes of sums
// beta and alpha are precomputed alpha and beta kernels.
// extpars array holds beyond LCDM parameters
// omega0 = Omega_{m,0}
// maxpars - maximum number of extended parameters - specified in SpecialFunctions.h

//ZW' edits starts
struct param_type3 {
  double kv[3]; // 0: |k1|, 1: |k2|, 2: |k3|
  double xv[3]; // 0: k2.k3, 1: k1.k2 , 2: k1.k3
  double args[3]; // 0: |k1+k2| , 1: |k1+k3| , 2: |k2+k3|
	double beta[3]; // beta arguments
 	double alpha[6]; // alpha arguments
  double omega0; // can be replaced by pars[maxbasepars] ....
	double extpars[maxpars];
	//ZW
	double omeganu;
	//ZW
};
//ZW' edits ends

int jacb (double a, const double G[], double *dfdy, double dfdt[], void *params)
{
	return GSL_SUCCESS;
}


// k1.k2=k1k2x2 , k1.k3 = k1k3x3, k2.k3=k2k3x1
int funcbtdgp(double a, const double G[], double F[], void *params)
{
	param_type3 p = *(param_type3 *)(params);

    /* Background quantities - LCDM  */
      double hub = HA2(a,p.omega0);
      double hub1= HA(a,p.omega0);
      double hubsqr = pow2(hub1);
      double acub = pow3(a);

	  /* DGP Poisson modifications up to 2nd order */
	  // ZW's edit starts
      double betadgp = 1./(3.*(mu(a,1.,p.omega0,p.omeganu, p.extpars,3)-1.));
	  // ZW's edit ends
      double betadgpcub = pow3(betadgp);
      double mua = 1.+1./(3.*betadgp);
      double gam2= -1./(hubsqr*24.*betadgpcub*p.extpars[0])*pow2(p.omega0/acub);
      double gamk[3];
      gamk[0] = (1.-pow2(p.xv[0]));
      gamk[1] = (1.-pow2(p.xv[1]));
      gamk[2] = (1.-pow2(p.xv[2]));

	/* 1st order */
	//1. F1/G1(k1)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hub)*G[1]-hub*G[0]*mua);

	//2. F1/G1(k3)
	F[2] = -G[3]/a;
	F[3] =1./a*(-(2.-hub)*G[3]-hub*G[2]*mua);

	//3. F1/G1(k2)
	F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hub)*G[5]-hub*G[4]*mua);

	/* 2nd order */

	//4. F2/G2(k2,k3) (P22)
	F[6] =1./a*(-(p.alpha[0]*G[5]*G[2]+p.alpha[1]*G[3]*G[4])/2.-G[7]);
	F[7] =1./a*(-(2.-hub)*G[7]-hub*G[6]*mua - gam2*gamk[0]*G[4]*G[2] - p.beta[0]*G[5]*G[3]);

	//5. F2/G2(k1,k3)
	F[8] =1./a*(-(p.alpha[2]*G[1]*G[2]+p.alpha[3]*G[3]*G[0])/2.-G[9]);
	F[9] =1./a*(-(2.-hub)*G[9]-hub*G[8]*mua - gam2*gamk[2]*G[2]*G[0] - p.beta[1]*G[3]*G[1]);

	//7. F2/G2(k1,k2)
	F[10] =1./a*(-(p.alpha[4]*G[5]*G[0]+p.alpha[5]*G[1]*G[4])/2.-G[11]) ;
	F[11] =1./a*(-(2.-hub)*G[11]-hub*G[10]*mua - gam2*gamk[1]*G[4]*G[0]-p.beta[2]*G[5]*G[1]);


	return GSL_SUCCESS;
}


// equations for fofr
int funcbtfofr(double a, const double G[], double F[], void *params)
{
	param_type3 p = *(param_type3 *)(params);

		/* Background quantities - LCDM  */
      double hub = HA2(a,p.omega0);
      double hub1= HA(a,p.omega0);
      double hubsqr = pow2(hub1);
      double acub = pow3(a);
			double ap6 = pow2(acub);
			double ap12 = pow2(ap6);
 			double h0 = myh0sqr;

		/* f(R) Poisson modifications up to 2nd order */
			double fofrp = p.extpars[0]/h0;
			double fofrp2 = pow2(fofrp);
			double om0a3 = pow2(p.omega0/acub);
			double term0 = pow2(3.*p.omega0-4.);
			double term0sqr = pow2(term0);
			double term1a = p.omega0-4.*acub*(p.omega0-1.);
			double term1b = pow2(term1a);
			double term1bsqr= pow2(term1b);
			double term1 = term1b*term1a/(2.*ap6*acub*fofrp*term0);
			double koa1 = pow2(p.kv[0]/a) ;
			double koa2 = pow2(p.kv[1]/a) ;
			double koa3 = pow2(p.kv[2]/a) ;
			double koa23 = pow2(p.args[2]/a) ;
			double koa13 = pow2(p.args[1]/a) ;
			double koa12 = pow2(p.args[0]/a) ;

      // f(R)
			// could optimise in terms of computations
			// 1st order
      double muak1 = 1. + koa1/(3.*(koa1 + term1));
			double muak2 = 1. + koa2/(3.*(koa2 + term1));
			double muak3 = 1. + koa3/(3.*(koa3 + term1));
			double muak23 = 1. + koa23/(3.*(koa23 + term1));
			double muak13 = 1. + koa13/(3.*(koa13 + term1));
			double muak12 = 1. + koa12/(3.*(koa12 + term1));

      double gam223 =  -(9.*koa23*om0a3*term1bsqr*term1a)/
						    				(48.*ap12*acub*fofrp2*hubsqr*term0sqr
			 			   					*(koa23+term1)
			 			   					*(koa2+term1)
			 	 		   					*(koa3+term1));


			double gam213 =  -(9.*koa13*om0a3*term1bsqr*term1a)/
						    				(48.*ap12*acub*fofrp2*hubsqr*term0sqr
			 			   					*(koa13+term1)
			 			   					*(koa1+term1)
			 	 		   					*(koa3+term1));


  		double gam212 =  -(9.*koa12*om0a3*term1bsqr*term1a)/
					    			   	(48.*ap12*acub*fofrp2*hubsqr*term0sqr
		 			   				  	*(koa12+term1)
		 			   				  	*(koa2+term1)
		 	 		   				  	*(koa1+term1));

	/* 1st order */
	//1. F1/G1(k1)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hub)*G[1]-hub*G[0]*muak1);

	//2. F1/G1(k3)
	F[2] = -G[3]/a;
	F[3] =1./a*(-(2.-hub)*G[3]-hub*G[2]*muak3);

	//3. F1/G1(k2)
	F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hub)*G[5]-hub*G[4]*muak2);

	/* 2nd order */

	//4. F2/G2(k2,k3) (P22)
	F[6] =1./a*(-(p.alpha[0]*G[5]*G[2]+p.alpha[1]*G[3]*G[4])/2.-G[7]);
	F[7] =1./a*(-(2.-hub)*G[7]-hub*G[6]*muak23 - gam223*G[4]*G[2] - p.beta[0]*G[5]*G[3]);


	//5. F2/G2(k1,k3)
	F[8] =1./a*(-(p.alpha[2]*G[1]*G[2]+p.alpha[3]*G[3]*G[0])/2.-G[9]);
	F[9] =1./a*(-(2.-hub)*G[9]-hub*G[8]*muak13 - gam213*G[2]*G[0] - p.beta[1]*G[3]*G[1]);

	//7. F2/G2(k1,k2)
	F[10] =1./a*(-(p.alpha[4]*G[5]*G[0]+p.alpha[5]*G[1]*G[4])/2.-G[11]) ;
	F[11] =1./a*(-(2.-hub)*G[11]-hub*G[10]*muak12 - gam212*G[4]*G[0]-p.beta[2]*G[5]*G[1]);


	return GSL_SUCCESS;
}

//ZW' edit starts
void BSPTN::initnb0_dgp(double pars[], double extpars[], double k[], double x[], double kargs[])
{
				// initial scale factor
				double a = 0.0001;
				double omega0 = pars[1]; // total matter fraction
				double A = pars[0]; // target scale factor
				//ZW's edit starts
				double omeganu = pars[2]; // massive neutrino fraction 
				//ZW's edit ends

        // Non-Eds ICs
			  double G[12] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.};

			/*Parameters passed to system of equations */
				struct param_type3 mypars;

			// set all parameters
				for(int i=0; i<3;i++){
					mypars.kv[i]=k[i];
					mypars.xv[i]=x[i];
					mypars.args[i]=kargs[i];
				}
			// calculate alpha and beta kernels
				mypars.alpha[0] = alphai(k[1],k[2],x[0]);
				mypars.alpha[1] = alphai(k[2],k[1],x[0]);
				mypars.beta[0] = betai(k[1],k[2],x[0]);

				mypars.alpha[2] = alphai(k[0],k[2],x[2]);
				mypars.alpha[3] = alphai(k[2],k[0],x[2]);
				mypars.beta[1] = betai(k[0],k[2],x[2]);

				mypars.alpha[4] = alphai(k[1],k[0],x[1]);
				mypars.alpha[5] = alphai(k[0],k[1],x[1]);
				mypars.beta[2] = betai(k[0],k[1],x[1]);

				mypars.omega0 = omega0;
				//ZW's edit starts
				mypars.omeganu = omeganu;
				//ZW's edit ends

				for (int i=0; i<maxpars;i++){
				  mypars.extpars[i] = 	extpars[i];
				}

				gsl_odeiv2_system sys = {funcbtdgp, jacb, 12, &mypars};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-6, 1e-6, 1e-6);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


		/*Allocation of array values */

		// For tree level initialisation

  		/*1st order */

			//F1(k;a), G1(k;a)
			F1b_k1 = G[0] ;
			G1b_k1 = G[1] ;

			// F1(k3;a), G1(k3;a)
			F1b_k3 = G[2];
			G1b_k3 = G[3];

			//F1(k2;a), G1(k2;a)
			F1b_k2 = G[4];
			G1b_k2 = G[5];

			/*2nd order*/

			//F2/G2(k3,k2) for tree,  F2/G2(k1-p,k2+p) for B222
			F2b_k23 =  G[6];
			G2b_k23 =  G[7];

      //F2/G2(k1,k3) for tree, F2/G2(-p,k2+p) for B222
      F2b_k13 =  G[8];
      G2b_k13 =  G[9];

      //F2/G2(k1,k2) for tree, F2/G2(p,k1-p)
      F2b_k12 =  G[10];
      G2b_k12 =  G[11];


			gsl_odeiv2_driver_free(d);
}
//ZW' edit ends

void BSPTN::initnb0_fofr(double pars[], double extpars[], double k[], double x[], double kargs[])
{
				double a = 0.0001;

				double A = pars[0];
				double omega0 = pars[1];

        // Non-Eds ICs
			  double G[12] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.};


				/*Parameters passed to system of equations */
					struct param_type3 mypars;

				// set all parameters
					for(int i=0; i<3;i++){
						mypars.kv[i]=k[i];
						mypars.xv[i]=x[i];
						mypars.args[i]=kargs[i];
					}

				// calculate alpha and beta kernels
					mypars.alpha[0] = alphai(k[1],k[2],x[0]);
					mypars.alpha[1] = alphai(k[2],k[1],x[0]);
					mypars.beta[0] = betai(k[1],k[2],x[0]);

					mypars.alpha[2] = alphai(k[0],k[2],x[2]);
					mypars.alpha[3] = alphai(k[2],k[0],x[2]);
					mypars.beta[1] = betai(k[0],k[2],x[2]);

					mypars.alpha[4] = alphai(k[1],k[0],x[1]);
					mypars.alpha[5] = alphai(k[0],k[1],x[1]);
					mypars.beta[2] = betai(k[0],k[1],x[1]);

					mypars.omega0 = omega0;

					for (int i=0; i<maxpars;i++){
					  mypars.extpars[i] = 	extpars[i];
					}

				gsl_odeiv2_system sys = {funcbtfofr, jacb, 12, &mypars};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
										   1e-6, 1e-6, 1e-6);

				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

				// For tree level initialisation

  		/*1st order */

			//F1(k;a), G1(k;a)
			F1b_k1 = G[0] ;
			G1b_k1 = G[1] ;

			// F1(k3;a), G1(k3;a)
			F1b_k3 = G[2];
			G1b_k3 = G[3];

			//F1(k2;a), G1(k2;a)
			F1b_k2 = G[4];
			G1b_k2 = G[5];

			/*2nd order*/

			//F2/G2(k3,k2) for tree
			F2b_k23 =  G[6];
			G2b_k23 =  G[7];

      //F2/G2(k1,k3) for tree
      F2b_k13 =  G[8];
      G2b_k13 =  G[9];

      //F2/G2(k1,k2) for tree, F2/G2(p,k1-p)
      F2b_k12 =  G[10];
      G2b_k12 =  G[11];


			gsl_odeiv2_driver_free(d);
}



/* Euler and Continuity equations for 1-loop Bispectrum numerical kernels */

/* Parameters passed to system of Euler and continuity equations*/
// k (magnitude) and x (angular) values for the system of equations
// args hold wave vector magnitudes of sums
// beta and alpha are precomputed alpha and beta kernels.
// extpars array holds beyond LCDM parameters
// omega0 = Omega_{m,0}
// maxpars - maximum number of extended parameters - specified in SpecialFunctions.h

struct param_type5 {
  double kv[8]; // 0: |p|, 1: |k1-p|, 2: |k2+p|, 3: |k2-p|, 4: |k3-p|, 5: |k1|, 6: |k2|, 7:|k3=-k1-k2|
	double xv[8]; // 0: k1-p.k2+p, 1: p.k1-p, 2: -p.k2+p , 3: p.k2-p, 4: p.k3-p, 5:  p.k1, 6: p.k2, 6: p.k3
  double args[25]; // see below
	double beta[38]; // beta arguments // see initnb1_dgp or initnb1_fr
	double alpha[76]; // alpha arguments  // see initnb1_dgp or initnb1_fr
	double gamk[30]; // gamma2,gamma3 k-dep pieces   // see initnb1_dgp or initnb1_fr
  double omega0;
	double extpars[maxpars];
};

//  p.arg1; // k2 + 2p // Note: arg1 = arg[0], arg2=arg[1] etc
//  p.arg2; // k1 + p
//  p.arg3; // -k3-p
//  p.arg4; // k2-p . k1
//  p.arg5; // k2-p . p+k1
//  p.arg6; // k1.k2 = x
//  p.arg7; // p. -k3-p
//  p.arg8; // k3-p. k1
//  p.arg9; // k3-p. k1+p
//  p.arg10; // k1.k3
//  p.arg11; // k2.k1+p
//  p.arg12; // k3-p.k2
//  p.arg13; // k2+p.k3-p
//  p.arg14; // k2.k3
//  p.arg15; // - p . k1+p
//  p.arg16; // k1-p.k2
//  p.arg17; // k2. -k3-p
//  p.arg18; // k2-p.k3
//  p.arg19; // k2-p.p+k3
//  p.arg20; // k1-p . k3
//  p.arg21; // k1-p.p+k3
//  p.arg22; // k1.k2+p
//  p.arg23; // k3.k1+p
//  p.arg24; // k2+p.k3
//  p.arg25; // k1.k3+p

int funcdgp(double a, const double G[], double F[], void *params)
{
	param_type5 p = *(param_type5 *)(params);

/* Background quantities - LCDM  */
  double hub = HA2(a,p.omega0);
  double hub1= HA(a,p.omega0);
  double hubsqr = pow2(hub1);
  double acub = pow3(a);

	/* DGP Poisson modifications up to 4th order */
  double betadgp = 1.+hub1/sqrt(p.extpars[0])*(1.+HA1(a,p.omega0)/(3.*hubsqr));
  double betadgpcub = pow3(betadgp);
  double om0a = pow2(p.omega0/acub);
  double mua = 1.+1./(3.*betadgp);
  double gam2 = -1./(hubsqr*24.*betadgpcub*p.extpars[0])*om0a;
  double gam3 = 1./(hubsqr*144.*betadgpcub*pow2(betadgp)*pow2(p.extpars[0]))*om0a*(p.omega0/acub);
  double gam4 = -1./(hubsqr*1728.*pow2(betadgpcub)*betadgp*pow3(p.extpars[0]))*pow2(om0a);

	// GR
	// double mua = 1.;
	// double gam2 = 1.;
	// double gam3 = 1.;
	// double gam4 = 1.;

	/* 1st order */
	//1. F1/G1(k1) or F1/G1(p)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hub)*G[1]-hub*G[0]*mua);

	//2. F1/G1(k3) or F1/G1(k2+p)
	F[2] = -G[3]/a;
	F[3] =1./a*(-(2.-hub)*G[3]-hub*G[2]*mua);

	//3. F1/G1(k2) or F1/G1(k1-p)
	F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hub)*G[5]-hub*G[4]*mua);

	/* 2nd order for B222 */

	//4. F2/G2(k2,k3)  or F2/G2(k1-p,k2+p)
	F[6] =1./a*(-(p.alpha[0]*G[5]*G[2]+p.alpha[1]*G[3]*G[4])/2.-G[7]);
	F[7] =1./a*(-(2.-hub)*G[7]-hub*G[6]*mua - gam2*p.gamk[0]*G[4]*G[2] - p.beta[0]*G[5]*G[3]);


	//5. F2/G2(-k1,k3)  or F2/G2(-p,k2+p)
	F[8] =1./a*(-(p.alpha[2]*G[1]*G[2]+p.alpha[3]*G[3]*G[0])/2.-G[9]);
	F[9] =1./a*(-(2.-hub)*G[9]-hub*G[8]*mua - gam2*p.gamk[2]*G[2]*G[0] - p.beta[1]*G[3]*G[1]);



	//7. F2/G2(k2,k1) or F2/G2(k1-p,p)
	F[10] =1./a*(-(p.alpha[4]*G[5]*G[0]+p.alpha[5]*G[1]*G[4])/2.-G[11]) ;
	F[11] =1./a*(-(2.-hub)*G[11]-hub*G[10]*mua - gam2*p.gamk[1]*G[4]*G[0]-p.beta[2]*G[5]*G[1]);


  /* 1st order for B321 - I */

  // F1/G1(k4) or F1/G1(k2-p)

  F[12] = -G[13]/a;
  F[13] =1./a*(-(2.-hub)*G[13]-hub*G[12]*mua);

  // F1/G1(k5) or F1/G1(k3-p)

  F[14] = -G[15]/a;
  F[15] =1./a*(-(2.-hub)*G[15]-hub*G[14]*mua);

  // F1/G1(k6) or F1/G1(k1)

  F[16] = -G[17]/a;
  F[17] =1./a*(-(2.-hub)*G[17]-hub*G[16]*mua);

  // F1/G1(k7) or F1/G1(k2)

  F[18] = -G[19]/a;
  F[19] =1./a*(-(2.-hub)*G[19]-hub*G[18]*mua);

  // F1/G1(k8) or F1/G1(k3)

  F[20] = -G[21]/a;
  F[21] =1./a*(-(2.-hub)*G[21]-hub*G[20]*mua);

  // 2nd order used in 3rd order equations of B321-I

  // F2/G2(k4,k1)  =  (k2-p,p)
  F[22] =1./a*(-(p.alpha[6]*G[13]*G[0]+p.alpha[7]*G[1]*G[12])/2.-G[23]) ;
  F[23] =1./a*(-(2.-hub)*G[23]-hub*G[22]*mua - gam2*p.gamk[3]*G[12]*G[0]-p.beta[3]*G[13]*G[1]);



  // F2/G2(k5,k1) = (k3-p,p)
  F[24] =1./a*(-(p.alpha[8]*G[15]*G[0]+p.alpha[9]*G[1]*G[14])/2.-G[25]) ;
  F[25] =1./a*(-(2.-hub)*G[25]-hub*G[24]*mua - gam2*p.gamk[4]*G[14]*G[0]-p.beta[4]*G[15]*G[1]);



  // F2/G2(k6,k1) = (k1,p)
  F[26] =1./a*(-(p.alpha[10]*G[17]*G[0]+p.alpha[11]*G[1]*G[16])/2.-G[27]) ;
  F[27] =1./a*(-(2.-hub)*G[27]-hub*G[26]*mua - gam2*p.gamk[5]*G[16]*G[0]-p.beta[5]*G[17]*G[1]);



  // F2/G2(k6,k4) = (k1,k2-p)
  F[28] =1./a*(-(p.alpha[12]*G[17]*G[12]+p.alpha[13]*G[13]*G[16])/2.-G[29]) ;
  F[29] =1./a*(-(2.-hub)*G[29]-hub*G[28]*mua - gam2*p.gamk[8]*G[16]*G[12]-p.beta[6]*G[17]*G[13]);

  // F3/G3(k6,k1,k4) = (k1,p,k2-p)
 F[30] = - 1./(3.*a)*(p.alpha[14]*G[26]*G[13]
           + p.alpha[15]*G[28]*G[1]
           + p.alpha[16]*G[22]*G[17]

           + p.alpha[17]*G[27]*G[12]
           + p.alpha[18]*G[29]*G[0]
           + p.alpha[19]*G[23]*G[16]

           +3.*G[31]) ;


 F[31] =1./(3.*a)*(-3.*(2.-hub)*G[31]-3.*hub*G[30]*mua

          -2.*p.beta[9]*G[17]*G[23]

          -2.*p.beta[8]*G[1]*G[29]

          -2.*p.beta[7]*G[13]*G[27]

          -2.*gam2*p.gamk[10]*G[16]*G[22]
          -2.*gam2*p.gamk[11]*G[0]*G[28]
          -2.*gam2*p.gamk[9]*G[12]*G[26]

          -gam3*(p.gamk[10]*p.gamk[3] + p.gamk[11]*p.gamk[8] + p.gamk[9]*p.gamk[5])*G[12]*G[16]*G[0]);



  // F2/G2(k5,k6) = (k3-p,k1)
  F[32] =1./a*(-(p.alpha[20]*G[15]*G[16]+p.alpha[21]*G[14]*G[17])/2.-G[33]) ;
  F[33] =1./a*(-(2.-hub)*G[33]-hub*G[32]*mua - gam2*p.gamk[12]*G[16]*G[14]-p.beta[10]*G[17]*G[15]);


// F3/G3(k6,k1,k5) = (k1,p,k3-p)
 F[34] = - 1./(3.*a)*(p.alpha[22]*G[26]*G[15]
           + p.alpha[2]*G[32]*G[1]
           + p.alpha[23]*G[24]*G[17]

           + p.alpha[24]*G[27]*G[14]
           + p.alpha[3]*G[33]*G[0]
           + p.alpha[25]*G[25]*G[16]

           +3.*G[35]) ;


 F[35] =1./(3.*a)*(-3.*(2.-hub)*G[35]-3.*hub*G[34]*mua

          -2.*p.beta[12]*G[17]*G[25]

          -2.*p.beta[1]*G[1]*G[33]

          -2.*p.beta[11]*G[15]*G[27]

          -2*gam2*p.gamk[14]*G[16]*G[24]
          -2*gam2*p.gamk[2]*G[0]*G[32]
          -2*gam2*p.gamk[13]*G[14]*G[26]

            -gam3*(p.gamk[14]*p.gamk[4] + p.gamk[2]*p.gamk[12] + p.gamk[13]*p.gamk[5])*G[14]*G[16]*G[0]);


 // F2/G2(k7,k1) = (k2,p)
 F[36] =1./a*(-(p.alpha[26]*G[19]*G[0]+p.alpha[27]*G[1]*G[18])/2.-G[37]) ;
 F[37] =1./a*(-(2.-hub)*G[37]-hub*G[36]*mua - gam2*p.gamk[6]*G[18]*G[0]-p.beta[13]*G[19]*G[1]);

 // F2/G2(k5,k7) = (k3-p,k2)
 F[38] =1./a*(-(p.alpha[28]*G[15]*G[18]+p.alpha[29]*G[19]*G[14])/2.-G[39]) ;
 F[39] =1./a*(-(2.-hub)*G[39]-hub*G[38]*mua - gam2*p.gamk[16]*G[14]*G[18]-p.beta[14]*G[15]*G[19]);


 // F3/G3(k7,k1,k5) = (k2,p,k3-p)
  F[40] = - 1./(3.*a)*(p.alpha[30]*G[36]*G[15]
            + p.alpha[31]*G[38]*G[1]
            + p.alpha[32]*G[24]*G[19]

            + p.alpha[33]*G[37]*G[14]
            + p.alpha[34]*G[39]*G[0]
            + p.alpha[35]*G[25]*G[18]

            +3.*G[41]) ;


  F[41] =1./(3.*a)*(-3.*(2.-hub)*G[41]-3.*hub*G[40]*mua

           -2.*p.beta[17]*G[19]*G[25]

           -2.*p.beta[16]*G[1]*G[39]

           -2.*p.beta[15]*G[15]*G[37]

           -2*gam2*p.gamk[18]*G[18]*G[24]
           -2*gam2*p.gamk[19]*G[0]*G[38]
           -2*gam2*p.gamk[17]*G[14]*G[36]

             -gam3*(p.gamk[18]*p.gamk[4] + p.gamk[19]*p.gamk[16] + p.gamk[17]*p.gamk[6])*G[14]*G[18]*G[0]);


 // F2/G2(k2,k7) = (k1-p,k2)
 F[42] =1./a*(-(p.alpha[36]*G[19]*G[2]+p.alpha[37]*G[3]*G[18])/2.-G[43]) ;
 F[43] =1./a*(-(2.-hub)*G[43]-hub*G[42]*mua - gam2*p.gamk[20]*G[18]*G[2]-p.beta[18]*G[19]*G[3]);


//F3/G3(k1-p,p,k2)
 F[44] = - 1./(3.*a)*(p.alpha[0]*G[36]*G[5]
           + p.alpha[15]*G[42]*G[1]
           + p.alpha[19]*G[10]*G[19]

           + p.alpha[1]*G[37]*G[4]
           + p.alpha[18]*G[43]*G[0]
           + p.alpha[16]*G[11]*G[18]

           +3.*G[45]) ;


 F[45] =1./(3.*a)*(-3.*(2.-hub)*G[45]-3.*hub*G[44]*mua

          -2.*p.beta[9]*G[19]*G[11]

          -2.*p.beta[8]*G[1]*G[43]

          -2.*p.beta[0]*G[5]*G[37]

          -2*gam2*p.gamk[0]*G[4]*G[36]
          -2*gam2*p.gamk[11]*G[0]*G[42]
          -2*gam2*p.gamk[10]*G[18]*G[10]

            -gam3*(p.gamk[0]*p.gamk[6] + p.gamk[11]*p.gamk[20] + p.gamk[10]*p.gamk[1])*G[4]*G[18]*G[0]);

 // F2/G2(k8,k1) = (k3,p)
 F[46] =1./a*(-(p.alpha[38]*G[21]*G[0]+p.alpha[39]*G[1]*G[20])/2.-G[47]) ;
 F[47] =1./a*(-(2.-hub)*G[47]-hub*G[46]*mua - gam2*p.gamk[7]*G[20]*G[0]-p.beta[19]*G[21]*G[1]);

 // F2/G2(k4,k8) = (k2-p,k3)
 F[48] =1./a*(-(p.alpha[40]*G[21]*G[12]+p.alpha[41]*G[13]*G[20])/2.-G[49]) ;
 F[49] =1./a*(-(2.-hub)*G[49]-hub*G[48]*mua - gam2*p.gamk[22]*G[20]*G[12]-p.beta[20]*G[21]*G[13]);




 // F3/G3(k8,k1,k4) = (k3,p,k2-p)
 F[50] = - 1./(3.*a)*(
             p.alpha[42]*G[46]*G[13]
           + p.alpha[31]*G[48]*G[1]
           + p.alpha[35]*G[22]*G[21]

           + p.alpha[43]*G[47]*G[12]
           + p.alpha[34]*G[49]*G[0]
           + p.alpha[32]*G[23]*G[20]

           +3.*G[51]) ;


 F[51] =1./(3.*a)*(-3.*(2.-hub)*G[51]-3.*hub*G[50]*mua

          -2.*p.beta[17]*G[21]*G[23]

          -2.*p.beta[16]*G[1]*G[49]

          -2.*p.beta[21]*G[13]*G[47]

          -2*gam2*p.gamk[18]*G[20]*G[22]
          -2*gam2*p.gamk[19]*G[0]*G[48]
          -2*gam2*p.gamk[23]*G[12]*G[46]

            -gam3*(p.gamk[18]*p.gamk[3] + p.gamk[19]*p.gamk[22] + p.gamk[23]*p.gamk[7])*G[12]*G[20]*G[0]);

 // F2/G2(k2,k8) = (k1-p,k3)
 F[52] =1./a*(-(p.alpha[44]*G[21]*G[2]+p.alpha[45]*G[3]*G[20])/2.-G[53]) ;
 F[53] =1./a*(-(2.-hub)*G[53]-hub*G[52]*mua - gam2*p.gamk[24]*G[20]*G[2]-p.beta[22]*G[21]*G[3]);


// F3/G3(k3,p,k1-p)
 F[54] = - 1./(3.*a)*(p.alpha[46]*G[46]*G[3]
           + p.alpha[2]*G[52]*G[1]
           + p.alpha[25]*G[10]*G[21]

           + p.alpha[47]*G[47]*G[2]
           + p.alpha[3]*G[53]*G[0]
           + p.alpha[23]*G[11]*G[20]

           +3.*G[55]) ;


 F[55] =1./(3.*a)*(-3.*(2.-hub)*G[55]-3.*hub*G[54]*mua

          -2.*p.beta[12]*G[21]*G[11]

          -2.*p.beta[1]*G[1]*G[53]

          -2.*p.beta[23]*G[3]*G[47]

          -2*gam2*p.gamk[14]*G[20]*G[10]
          -2*gam2*p.gamk[2]*G[0]*G[52]
          -2*gam2*p.gamk[25]*G[2]*G[46]

            -gam3*(p.gamk[14]*p.gamk[1] + p.gamk[2]*p.gamk[24] + p.gamk[25]*p.gamk[7])*G[2]*G[20]*G[0]);


// F2/G2(k6,k7) = F2/G2(k1,k2)
F[56] =1./a*(-(p.alpha[16]*G[17]*G[18]+p.alpha[19]*G[19]*G[16])/2.-G[57]) ;
F[57] =1./a*(-(2.-hub)*G[57]-hub*G[56]*mua - gam2*p.gamk[10]*G[18]*G[16]-p.beta[9]*G[17]*G[19]);


// F2/G2(k6,k8) = F2/G2(k1,k3)
F[58] =1./a*(-(p.alpha[23]*G[17]*G[20]+p.alpha[25]*G[21]*G[16])/2.-G[59]) ;
F[59] =1./a*(-(2.-hub)*G[59]-hub*G[58]*mua - gam2*p.gamk[14]*G[20]*G[16]-p.beta[12]*G[17]*G[21]);


// F2/G2(k7,k8) = F2/G2(k2,k3)
F[60] =1./a*(-(p.alpha[32]*G[19]*G[20]+p.alpha[35]*G[21]*G[18])/2.-G[61]) ;
F[61] =1./a*(-(2.-hub)*G[61]-hub*G[60]*mua - gam2*p.gamk[18]*G[20]*G[18]-p.beta[17]*G[19]*G[21]);


// F2/G2(-k3,p)
F[62] =1./a*(-(p.alpha[48]*G[21]*G[0]+p.alpha[49]*G[1]*G[20])/2.-G[63]) ;
F[63] =1./a*(-(2.-hub)*G[63]-hub*G[62]*mua - gam2*p.gamk[7]*G[20]*G[0]-p.beta[25]*G[21]*G[1]);


 // F2/G2(-k2,p)
 F[64] =1./a*(-(p.alpha[50]*G[19]*G[0]+p.alpha[51]*G[1]*G[18])/2.-G[65]) ;
 F[65] =1./a*(-(2.-hub)*G[65]-hub*G[64]*mua- gam2*p.gamk[6]*G[18]*G[0]-p.beta[26]*G[19]*G[1]);


 // F2/G2(k6,k1) = (-k1,p)
 F[66] =1./a*(-(p.alpha[52]*G[17]*G[0]+p.alpha[53]*G[1]*G[16])/2.-G[67]) ;
 F[67] =1./a*(-(2.-hub)*G[67]-hub*G[66]*mua - gam2*p.gamk[5]*G[16]*G[0]-p.beta[27]*G[17]*G[1]);


 // F3/G3(k6,k7,k1) = (-k1,-k2,p)
 F[68] = - 1./(3.*a)*(p.alpha[39]*G[56]*G[1]
                    + p.alpha[36]*G[66]*G[19]
                    + p.alpha[12]*G[64]*G[17]

           + p.alpha[38]*G[57]*G[0]
           + p.alpha[37]*G[67]*G[18]
           + p.alpha[13]*G[65]*G[16]

           +3.*G[69]) ;

 F[69] =1./(3.*a)*(-3.*(2.-hub)*G[69]-3.*hub*G[68]*mua

          -2.*p.beta[6]*G[65]*G[17]

          -2.*p.beta[18]*G[67]*G[19]

          -2.*p.beta[19]*G[57]*G[1]

          -2.*gam2*p.gamk[8]*G[16]*G[64]
          -2.*gam2*p.gamk[20]*G[18]*G[66]
          -2.*gam2*p.gamk[7]*G[0]*G[56]

            -gam3*(p.gamk[8]*p.gamk[6] + p.gamk[20]*p.gamk[5] + p.gamk[7]*p.gamk[10])*G[0]*G[18]*G[16]);



// F3/G3(k6,k7,k1) = (-k1,-k2,-p)
F[70] = - 1./(3.*a)*(p.alpha[49]*G[56]*G[1]
                  + p.alpha[54]*G[26]*G[19]
                  + p.alpha[55]*G[36]*G[17]

         + p.alpha[48]*G[57]*G[0]
         + p.alpha[56]*G[27]*G[18]
         + p.alpha[57]*G[37]*G[16]

         +3.*G[71]) ;


F[71] =1./(3.*a)*(-3.*(2.-hub)*G[71]-3.*hub*G[70]*mua

        -2.*p.beta[29]*G[37]*G[17]

        -2.*p.beta[28]*G[27]*G[19]

        -2.*p.beta[25]*G[57]*G[1]

        -2.*gam2*p.gamk[26]*G[16]*G[36]
        -2.*gam2*p.gamk[15]*G[18]*G[26]
        -2.*gam2*p.gamk[7]*G[0]*G[56]

          -gam3*(p.gamk[26]*p.gamk[6] + p.gamk[15]*p.gamk[5] + p.gamk[7]*p.gamk[10])*G[0]*G[18]*G[16]);


F[72] =1./a*(-G[73]) ; // set 1st term to 0 because x1 is not identically -1 (needed to keep evolution at 0 for GR)
F[73] =1./a*(-(2.-hub)*G[73]-hub*G[72]*mua);


// commented out the terms that are 0.

F[74] = - 1./(3.*a)*(p.alpha[5]*G[66]*G[1]
                  + p.alpha[31]*G[26]*G[1]
            //      + p.alpha[58]*G[72]*G[17]

         				+ p.alpha[4]*G[67]*G[0]
         				+ p.alpha[34]*G[27]*G[0]
         		//		+ p.alpha[59]*G[73]*G[16]

         				+3.*G[75]) ;


F[75] =1./(3.*a)*(-3.*(2.-hub)*G[75]-3.*hub*G[74]*mua

    //    -2.*p.beta[30]*G[73]*G[17]

        -2.*p.beta[16]*G[27]*G[1]

        -2.*p.beta[2]*G[67]*G[1]

      //  -2*gam2*0.*G[16]*G[72]
        -2*gam2*p.gamk[19]*G[0]*G[26]
        -2*gam2*p.gamk[1]*G[0]*G[66]

          -gam3*(p.gamk[19]*p.gamk[5] + p.gamk[1]*p.gamk[5])*G[0]*G[0]*G[16]);


// F3/G3(k7,k1,k1) = (-k2,p,-p)
F[76] = - 1./(3.*a)*(p.alpha[7]*G[64]*G[1]
                   + p.alpha[2]*G[36]*G[1]
        //           + p.alpha[60]*G[72]*G[19]

                  + p.alpha[6]*G[65]*G[0]
                  + p.alpha[3]*G[37]*G[0]
          //        + p.alpha[61]*G[73]*G[18]

                  +3.*G[77]) ;


F[77] =1./(3.*a)*(-3.*(2.-hub)*G[77]-3.*hub*G[76]*mua

        -2.*p.beta[3]*G[65]*G[1]

        -2.*p.beta[1]*G[37]*G[1]

    //    -2.*p.beta[31]*G[73]*G[19]

      //  -2.*gam2*(1.-pow2(XMIN))*G[18]*G[72]
        -2.*gam2*p.gamk[2]*G[0]*G[36]
        -2.*gam2*p.gamk[3]*G[0]*G[64]

        -gam3*(p.gamk[2]*p.gamk[6] + p.gamk[3]*p.gamk[6] + (1.-pow2(XMIN)))*G[0]*G[0]*G[18]);


// F4/G4(-k6,-k7,k1,-k1) = (-k1,-k2,p,-p)
F[78] = -1./(24.*a)*(
										6.*p.alpha[16]*G[76]*G[17] // a(-k1,-k2) F3(-k2,p,-p) G1(k1)
                   +6.*p.alpha[19]*G[74]*G[19] // a(-k2,-k1) F3(-k1,p,-p) G1(k2)
                   +6.*p.alpha[9]*G[70]*G[1] // a(p,k3-p) F3(-p,-k1,-k2) G1(p)
                   +6.*p.alpha[15]*G[68]*G[1] // a(-p,k3+p) F3(p,-k1,-k2) G1(p)

                    +6.*p.alpha[19]*G[77]*G[16] // a(-k2,-k1) G3(-k2,p,-p) F1(k1)
                    +6.*p.alpha[16]*G[75]*G[18] // a(-k2,-k1) G3(-k1,p,-p) F1(k2)
                    +6.*p.alpha[8]*G[71]*G[0] // a(p,k3-p) G3(-p,-k1,-k2) F1(p)
                    +6.*p.alpha[18]*G[69]*G[0] // a(p,k3-p) G3(-p,-k1,-k2) F1(p)

                //    +4.*p.alpha[62]*G[56]*G[73] // a(0,k3) F2(k1,k2) G2(p,-p)
                //    +4.*p.alpha[63]*G[72]*G[57] // a(k3,0) F2(p,-p) G2(k1,k2)

                    + 4.*p.alpha[17]*G[64]*G[27] // a(-k1-p,-k2+p,) F2(-k2,p) G2(-p,-k1)
                    + 4.*p.alpha[14]*G[26]*G[65] // a(-k2+p,-k1-p) F2(-k1,-p) G2(-k2,p)

                    + 4.*p.alpha[0]*G[36]*G[67] // a(-k1+p,-k2-p) F2(-k2,-p)G2(-k1,p)
                    + 4.*p.alpha[1]*G[66]*G[37] // a(-k2-p,-k1+p) F2(-k1,p) G2(-k2,-p)

                    +24.*G[79]);

F[79] = 1./a*(-(2.-hub)*G[79]-hub*G[78]*mua

          -1./24.*(
                  12.*p.beta[9]*G[77]*G[17] // b(-k1,-k2) G3(-k2,p,-p) G1(k1)
                  +12.*p.beta[9]*G[75]*G[19] // b(-k2,-k1) G3(-k1,p,-p) G1(k2)
                  +12.*p.beta[4]*G[71]*G[1] // b(p,k3-p) G3(-p,-k1,-k2) G1(p)
                  +12.*p.beta[8]*G[69]*G[1] // b(-p,k3+p) G3(p,-k1,-k2) G1(p)

              //    + 8.*p.beta[24]*G[57]*G[73] // b(k3,0) G2(k1,k2) G2(p,-p)
                  + 8.*p.beta[7]*G[65]*G[27]  // b(k2-p,k1+p) G2(-k2,p) G2(k1,p)
                  + 8.*p.beta[0]*G[67]*G[37] // b(k2+p,k1-p) G2(k2,p) G2(-k1,p)

//2nd order
            //      + 8.*gam2*(1.-pow2(XMIN))*G[56]*G[72]
                  + 8.*gam2*p.gamk[9]*G[64]*G[26]
                  + 8.*gam2*p.gamk[0]*G[66]*G[36]

//3rd order
                  + 12.*gam2*p.gamk[10]*G[76]*G[16]
                  + 12.*gam2*p.gamk[10]*G[74]*G[18]
                  + 12.*gam2*p.gamk[4]*G[70]*G[0]
                  + 12.*gam2*p.gamk[11]*G[68]*G[0]


									+4.*gam3*((p.gamk[5]*p.gamk[9] + p.gamk[10]*p.gamk[3] + p.gamk[8]*p.gamk[11])*G[16]*G[0]*G[64]
										 			 +(p.gamk[6]*p.gamk[9] + p.gamk[10]*p.gamk[19] + p.gamk[15]*p.gamk[4])*G[18]*G[0]*G[26]
												   +(p.gamk[0]*p.gamk[6] + p.gamk[10]*p.gamk[1] + p.gamk[11]*p.gamk[20])*G[18]*G[0]*G[66]
												 	 +(p.gamk[0]*p.gamk[5] + p.gamk[10]*p.gamk[2] + p.gamk[4]*p.gamk[26])*G[16]*G[0]*G[36]
												 	 +(p.gamk[10]+ p.gamk[10] + p.gamk[10])*G[16]*G[18]*G[72]
												   +(p.gamk[7]*p.gamk[4] + p.gamk[7]*p.gamk[11])*G[0]*G[0]*G[56])

// 4th order
									 +  gam4*4.*(p.gamk[5]*p.gamk[6]*p.gamk[0]
										 				 + p.gamk[5]*p.gamk[6]*p.gamk[9]

										 				   + p.gamk[10]*p.gamk[6]*p.gamk[2]
															 + p.gamk[4]*p.gamk[6]*p.gamk[26]
															 + p.gamk[10]*p.gamk[6]*p.gamk[3]
															 + p.gamk[11]*p.gamk[6]*p.gamk[8]
															 + p.gamk[10]*p.gamk[5]*p.gamk[19]
															 + p.gamk[4]*p.gamk[5]*p.gamk[15]
															 + p.gamk[10]*p.gamk[5]*p.gamk[1]
															 + p.gamk[11]*p.gamk[5]*p.gamk[20]
															 + p.gamk[4]*p.gamk[10]*p.gamk[7]
															 + p.gamk[11]*p.gamk[10]*p.gamk[7])
															   *G[16]*G[18]*G[0]*G[0]));

// F3/G3(k6,k8,k1) = (-k1,-k3,p)
F[80] = - 1./(3.*a)*(p.alpha[21]*G[62]*G[17] // a(k1,k3-p) F2(-k3,p) G1(k1)
                   + p.alpha[44]*G[66]*G[21] // a(k3,k1-p) F2(-k1,p) G1(k3)
                   + p.alpha[27]*G[58]*G[1] // a(p,k2) F2(k1,k3) G1(p)

          + p.alpha[20]*G[63]*G[16]
          + p.alpha[45]*G[67]*G[20]
          + p.alpha[26]*G[59]*G[0]

          +3.*G[81]) ;

F[81] =1./(3.*a)*(-3.*(2.-hub)*G[81]-3.*hub*G[80]*mua

         -2.*p.beta[10]*G[63]*G[17]

         -2.*p.beta[22]*G[67]*G[21]

         -2.*p.beta[13]*G[59]*G[1]

         -2*gam2*p.gamk[12]*G[16]*G[62]
         -2*gam2*p.gamk[24]*G[20]*G[66]
         -2*gam2*p.gamk[6]*G[0]*G[58]

           -gam3*(p.gamk[12]*p.gamk[7] + p.gamk[24]*p.gamk[5] + p.gamk[6]*p.gamk[14])*G[0]*G[20]*G[16]);


// F3/G3(k8,k1,k1) = (-k3,p,-p)
F[82] = - 1./(3.*a)*(p.alpha[9]*G[62]*G[1] // a(p,k3-p) F2(k3,-p) G1(p)
                 + p.alpha[15]*G[46]*G[1] // a(p,-k3-p) F2(k3,p) G1(p)
          //       + p.alpha[63]*G[72]*G[20] // a(-k3,0.) F2(p,-p) G1(k3)

        + p.alpha[8]*G[63]*G[0]
        + p.alpha[18]*G[47]*G[0]
    //    + p.alpha[62]*G[73]*G[21]

        +3.*G[83]) ;


F[83] =1./(3.*a)*(-3.*(2.-hub)*G[83]-3.*hub*G[82]*mua

       -2.*p.beta[4]*G[63]*G[1]

       -2.*p.beta[8]*G[47]*G[1]

       -2.*p.beta[24]*G[73]*G[21]


  //     -2*gam2*0.*G[20]*G[72]
       -2*gam2*p.gamk[11]*G[0]*G[46]
       -2*gam2*p.gamk[4]*G[0]*G[62]

         -gam3*(p.gamk[11]*p.gamk[7] + p.gamk[4]*p.gamk[7])*G[0]*G[0]*G[20]);


// F3/G3(k6,k8,k1) = (-k1,-k3,-p)
F[84] = - 1./(3.*a)*(p.alpha[64]*G[46]*G[17] // a(k1,k3+p) F2(k3,p) G1(k1)
                  + p.alpha[65]*G[26]*G[21] // a(k3,k1+p) F2(k1,p) G1(k3)
                  + p.alpha[51]*G[58]*G[1] // a(p,-k2) F2(k1,k3) G1(p)

         + p.alpha[66]*G[47]*G[16]
         + p.alpha[67]*G[27]*G[20]
         + p.alpha[50]*G[59]*G[0]

         +3.*G[85]) ;

F[85] =1./(3.*a)*(-3.*(2.-hub)*G[85]-3.*hub*G[84]*mua

        -2.*p.beta[32]*G[47]*G[17]

        -2.*p.beta[33]*G[27]*G[21]

        -2.*p.beta[26]*G[59]*G[1]

			 -2*gam2*p.gamk[29]*G[16]*G[46]
       -2*gam2*p.gamk[27]*G[20]*G[26]
			 -2*gam2*p.gamk[6]*G[0]*G[58]

         -gam3*(p.gamk[29]*p.gamk[7] + p.gamk[27]*p.gamk[5]  + p.gamk[6]*p.gamk[14])*G[0]*G[20]*G[16]);


// F4/G4(-k6,-k8,k1,-k1) = (-k1,-k3,p,-p)
F[86] = -1./(24.*a)*(6.*p.alpha[23]*G[82]*G[17] // a(k1,k3) F3(-k3,p,-p) G1(k1)
                    +6.*p.alpha[25]*G[74]*G[21] // a(k3,k1) F3(-k1,p,-p) G1(k3)
                    +6.*p.alpha[7]*G[84]*G[1] // a(p,k2-p) F3(-p,-k1,-k3) G1(p)
                    +6.*p.alpha[70]*G[80]*G[1] // a(-p,k2+p) F3(p,-k1,-k3) G1(p)

                   +6.*p.alpha[25]*G[83]*G[16]
                   +6.*p.alpha[23]*G[75]*G[20]
                   +6.*p.alpha[6]*G[85]*G[0]
                   +6.*p.alpha[71]*G[81]*G[0]


          //         + 4.*p.alpha[60]*G[58]*G[73] // a(0,k3) F2(k1,k3) G2(p,-p)
          //         + 4.*p.alpha[61]*G[72]*G[59] // a(k3,0) F2(p,-p) G2(k1,k3)

                   + 4.*p.alpha[24]*G[62]*G[27] // a(k1+p,k3-p) F2(-k3,p) G2(-p,-k1)
                   + 4.*p.alpha[22]*G[26]*G[63] // a(k3-p,k1+p) F2(-k1,-p) G2(-k3,p)

                   + 4.*p.alpha[69]*G[46]*G[67] // a(k1-p,k3+p) F2(-k3,-p)G2(-k1,p)
                   + 4.*p.alpha[68]*G[66]*G[47] // a(k3+p,k1-p) F2(-k1,p) G2(-k3,-p)

                   +24.*G[87]);

F[87] = 1./a*(-(2.-hub)*G[87]-hub*G[86]*mua

            -1./24.*(
                  12.*p.beta[12]*G[83]*G[17] // b(k1,k3) G3(-k3,p,-p) G1(k1)
                 +12.*p.beta[12]*G[75]*G[21] // b(k3,k1) G3(-k1,p,-p) G1(k3)
                 +12.*p.beta[3]*G[85]*G[1] // b(p,k2-p) G3(-p,-k1,-k3) G1(p)
                 +12.*p.beta[34]*G[81]*G[1]  // b(-p,k2+p) G3(p,-k1,-k3) G1(p)

                 + 8.*p.beta[31]*G[59]*G[73] // b(k2,0) G2(k1,k3) G2(p,-p)
                 + 8.*p.beta[11]*G[63]*G[27]  // b(k3-p,k1+p) G2(-k3,p) G2(k1,p)
                 + 8.*p.beta[35]*G[67]*G[47] // b(k3+p,k1-p) G2(k3,p) G2(-k1,p)

            //     + 8.*gam2*0.*G[58]*G[72]
                 + 8.*gam2*p.gamk[13]*G[62]*G[26]
                 + 8.*gam2*p.gamk[25]*G[66]*G[46]

                 + 12.*gam2*p.gamk[14]*G[82]*G[16]
                 + 12.*gam2*p.gamk[14]*G[74]*G[20]
                 + 12.*gam2*p.gamk[3]*G[84]*G[0]
                 + 12.*gam2*p.gamk[2]*G[80]*G[0]

								 + 4.*gam3*((p.gamk[13]*p.gamk[5] + p.gamk[14]*p.gamk[4] + p.gamk[2]*p.gamk[12])*G[16]*G[0]*G[62]
													+(p.gamk[13]*p.gamk[7] + p.gamk[14]*p.gamk[19] + p.gamk[27]*p.gamk[3])*G[20]*G[0]*G[26]
													+(p.gamk[7]*p.gamk[25] + p.gamk[14]*p.gamk[1] + p.gamk[2]*p.gamk[24])*G[20]*G[0]*G[66]
													+(p.gamk[5]*p.gamk[25] + p.gamk[14]*p.gamk[11] + p.gamk[3]*p.gamk[29])*G[16]*G[0]*G[46]
													+(p.gamk[14]+ p.gamk[14] + p.gamk[14])*G[16]*G[20]*G[72]
													+(p.gamk[6]*p.gamk[3] + p.gamk[6]*p.gamk[2])*G[0]*G[0]*G[58])


								 +  gam4*4.*(p.gamk[5]*p.gamk[7]*p.gamk[25]
									 		     + p.gamk[5]*p.gamk[7]*p.gamk[13]

									 				   + p.gamk[14]*p.gamk[7]*p.gamk[11]
														 + p.gamk[3]*p.gamk[7]*p.gamk[29]
														 + p.gamk[14]*p.gamk[7]*p.gamk[4]
														 + p.gamk[2]*p.gamk[7]*p.gamk[12]
														 + p.gamk[14]*p.gamk[5]*p.gamk[19]
														 + p.gamk[3]*p.gamk[5]*p.gamk[27]
														 + p.gamk[14]*p.gamk[5]*p.gamk[1]
														 + p.gamk[2]*p.gamk[5]*p.gamk[24]
														 + p.gamk[3]*p.gamk[14]*p.gamk[6]
														 + p.gamk[2]*p.gamk[14]*p.gamk[6])
														 *G[16]*G[20]*G[0]*G[0]));




// F3/G3(k6,k8,k1) = (-k2,-k3,-p)
F[88] = - 1./(3.*a)*(p.alpha[72]*G[46]*G[19] // a(k2,k3+p) F2(k3,p) G1(k2)
                 + p.alpha[74]*G[36]*G[21] // a(k3,k2+p) F2(k2,p) G1(k3)
                 + p.alpha[53]*G[60]*G[1] // a(p,-k1) F2(k2,k3) G1(p)


        + p.alpha[73]*G[47]*G[18]
        + p.alpha[75]*G[37]*G[20]
        + p.alpha[52]*G[61]*G[0]

        +3.*G[89]) ;


F[89] =1./(3.*a)*(-3.*(2.-hub)*G[89]-3.*hub*G[88]*mua

       -2.*p.beta[36]*G[47]*G[19]

       -2.*p.beta[37]*G[37]*G[21]

       -2.*p.beta[27]*G[61]*G[1]


      -2*gam2*p.gamk[21]*G[18]*G[46]
      -2*gam2*p.gamk[28]*G[20]*G[36]
      -2*gam2*p.gamk[5]*G[0]*G[60]

         -gam3*(p.gamk[21]*p.gamk[7] + p.gamk[28]*p.gamk[6] + p.gamk[5]*p.gamk[18])*G[0]*G[20]*G[18]);



// F3/G3(k6,k8,k1) = (-k2,-k3,p)

 F[90] = - 1./(3.*a)*(p.alpha[29]*G[62]*G[19] // a(k2,k3-p) F2(k3,-p) G1(k2)
                  + p.alpha[40]*G[64]*G[21] // a(k3,k2-p) F2(k2,-p) G1(k3)
                  + p.alpha[11]*G[60]*G[1] // a(p,k1) F2(k2,k3) G1(p)

         + p.alpha[28]*G[63]*G[18]
         + p.alpha[41]*G[65]*G[20]
         + p.alpha[10]*G[61]*G[0]

         +3.*G[91]) ;



 F[91] =1./(3.*a)*(-3.*(2.-hub)*G[91]-3.*hub*G[90]*mua

        -2.*p.beta[14]*G[63]*G[19]

        -2.*p.beta[20]*G[65]*G[21]

        -2.*p.beta[5]*G[61]*G[1]

        -2*gam2*p.gamk[16]*G[18]*G[62]
        -2*gam2*p.gamk[22]*G[20]*G[64]
        -2*gam2*p.gamk[5]*G[0]*G[60]

          -gam3*(p.gamk[16]*p.gamk[7] + p.gamk[22]*p.gamk[6] + p.gamk[5]*p.gamk[18])*G[0]*G[20]*G[18]);


// F4/G4(-k6,-k8,k1,-k1) = (-k2,-k3,p,-p)
F[92] = -1./(24.*a)*(6.*p.alpha[32]*G[82]*G[19] // a(k2,k3) F3(-k3,p,-p) G1(k2)
                   +6.*p.alpha[35]*G[76]*G[21] // a(k3,k2) F3(-k2,p,-p) G1(k3)
                   +6.*p.alpha[5]*G[88]*G[1] // a(p,k1-p) F3(-p,-k2,-k3) G1(p)
                   +6.*p.alpha[31]*G[90]*G[1] // a(-p,k1+p) F3(p,-k2,-k3) G1(p)

                  +6.*p.alpha[35]*G[83]*G[18]
                  +6.*p.alpha[32]*G[77]*G[20]
                  +6.*p.alpha[4]*G[89]*G[0]
                  +6.*p.alpha[34]*G[91]*G[0]


            //      + 4.*p.alpha[58]*G[60]*G[73] // a(0,k1) F2(k2,k3) G2(p,-p)
            //      + 4.*p.alpha[59]*G[72]*G[61] // a(k1,0) F2(p,-p) G2(k2,k3)
                  + 4.*p.alpha[33]*G[62]*G[37] // a(k2+p,k3-p) F2(-k3,p) G2(p,k2)
                  + 4.*p.alpha[30]*G[36]*G[63] // a(k3-p,k2+p) F2(k2,p) G2(-k3,p)

                  + 4.*p.alpha[42]*G[46]*G[65] // a(k2-p,k3+p) F2(-k3,-p)G2(-k2,p)
                  + 4.*p.alpha[43]*G[64]*G[47] // a(k3+p,k2-p) F2(-k2,p) G2(-k3,-p)


                  +24.*G[93]);

F[93] = 1./a*(-(2.-hub)*G[93]-hub*G[92]*mua

        -1./24.*(
                 12.*p.beta[17]*G[83]*G[19] // b(k2,k3) G3(-k3,p,-p) G1(k2)
                +12.*p.beta[17]*G[77]*G[21] // b(k3,k2) G3(-k2,p,-p) G1(k3)
                +12.*p.beta[2]*G[89]*G[1]// b(p,k1-p) G3(-p,-k2,-k3) G1(p)
                +12.*p.beta[16]*G[91]*G[1]  // b(-p,k1+p) G3(p,-k2,-k3) G1(p)

          //      + 8.*p.beta[30]*G[61]*G[73] // b(k1,0) G2(k2,k3) G2(p,-p)
                + 8.*p.beta[15]*G[63]*G[37]  // b(k3-p,k2+p) G2(-k3,p) G2(k2,p)
                + 8.*p.beta[21]*G[65]*G[47] // b(k3+p,k2-p) G2(k3,p) G2(-k2,p)

        //        + 8.*gam2*0.*G[60]*G[72]
                + 8.*gam2*p.gamk[17]*G[62]*G[36]
                + 8.*gam2*p.gamk[23]*G[64]*G[46]

                + 12.*gam2*p.gamk[18]*G[82]*G[18]
                + 12.*gam2*p.gamk[18]*G[76]*G[20]
                + 12.*gam2*p.gamk[1]*G[88]*G[0]
                + 12.*gam2*p.gamk[19]*G[90]*G[0]

								+4.*gam3*((p.gamk[6]*p.gamk[17] + p.gamk[18]*p.gamk[4] + p.gamk[16]*p.gamk[19])*G[18]*G[0]*G[62]
												 +(p.gamk[17]*p.gamk[7] + p.gamk[18]*p.gamk[2] + p.gamk[28]*p.gamk[1])*G[20]*G[0]*G[36]
												 +(p.gamk[7]*p.gamk[23] + p.gamk[18]*p.gamk[3] + p.gamk[19]*p.gamk[22])*G[20]*G[0]*G[64]
												 +(p.gamk[6]*p.gamk[23] + p.gamk[18]*p.gamk[11] + p.gamk[1]*p.gamk[21])*G[18]*G[0]*G[46]
												 +(p.gamk[18]+ p.gamk[18] + p.gamk[18])*G[18]*G[20]*G[72]
												 +(p.gamk[5]*p.gamk[1] + p.gamk[5]*p.gamk[19])*G[0]*G[0]*G[60])

							  +  gam4*4.*(p.gamk[6]*p.gamk[7]*p.gamk[23] + p.gamk[6]*p.gamk[7]*p.gamk[17]

								 				   + p.gamk[18]*p.gamk[7]*p.gamk[11] + p.gamk[1]*p.gamk[7]*p.gamk[21]
													 + p.gamk[18]*p.gamk[7]*p.gamk[4] + p.gamk[19]*p.gamk[7]*p.gamk[16]
													 + p.gamk[18]*p.gamk[6]*p.gamk[2] + p.gamk[1]*p.gamk[6]*p.gamk[2]
													 + p.gamk[18]*p.gamk[6]*p.gamk[3] + p.gamk[19]*p.gamk[6]*p.gamk[22]
													 + p.gamk[1]*p.gamk[18]*p.gamk[5] + p.gamk[19]*p.gamk[18]*p.gamk[5])
													 *G[18]*G[20]*G[0]*G[0]));



  	return GSL_SUCCESS;
  }


// Bispectrum B411, B321 AND B222 initialisations
// k array holds vector magnitudes
// x array holds angle between vectors
// kargs holds magnitudes of combinations of vectors
void BSPTN::initnb1_dgp(double pars[], double extpars[], double epars[], double k[], double x[], double kargs[])
{
				double a = 0.0001;
				double A = pars[0];
				double omega0 = pars[1];

				// Non-EdS initial conditions
       double G[94] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.,a,-a,a,-a,a,-a,a,-a,a,-a,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
                            0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};


			/*Parameters passed to system of equations */
          struct param_type5 mypars;
					// set all parameters

						for(int i=0; i<8;i++){
							mypars.kv[i]=k[i];
							mypars.xv[i]=x[i];
						}
						for(int i=0; i<25;i++){
							mypars.args[i]=kargs[i];
						}

					// precompute alpha and beta kernels
					mypars.alpha[0] = alphai(k[1],k[2],x[0]);
					mypars.alpha[1] = alphai(k[2],k[1],x[0]);
					mypars.alpha[2] = alphai(k[0],k[2],x[2]);
					mypars.alpha[3] = alphai(k[2],k[0],x[2]);
					mypars.alpha[4] = alphai(k[1],k[0],x[1]);
					mypars.alpha[5] = alphai(k[0],k[1],x[1]);
					mypars.alpha[6] = alphai(k[3],k[0],x[3]);
					mypars.alpha[7] = alphai(k[0],k[3],x[3]);
					mypars.alpha[8] = alphai(k[4],k[0],x[4]);
					mypars.alpha[9] = alphai(k[0],k[4],x[4]);
					mypars.alpha[10] = alphai(k[5],k[0],x[5]);
					mypars.alpha[11] = alphai(k[0],k[5],x[5]);
					mypars.alpha[12] = alphai(k[5],k[3],kargs[3]);
					mypars.alpha[13] = alphai(k[3],k[5],kargs[3]);
					mypars.alpha[14] = alphai(k[3],kargs[1],kargs[4]);
					mypars.alpha[15] = alphai(k[0],kargs[2],kargs[6]);
					mypars.alpha[16] = alphai(k[5],k[6],kargs[5]);
					mypars.alpha[17] = alphai(kargs[1],k[3],kargs[4]);
					mypars.alpha[18] = alphai(kargs[2],k[0],kargs[6]);
					mypars.alpha[19] = alphai(k[6],k[5],kargs[5]);
					mypars.alpha[20] = alphai(k[4],k[5],kargs[7]);
					mypars.alpha[21] = alphai(k[5],k[4],kargs[7]);
					mypars.alpha[22] = alphai(k[4],kargs[1],kargs[8]);
					mypars.alpha[23] = alphai(k[5],k[7],kargs[9]);
					mypars.alpha[24] = alphai(kargs[1],k[4],kargs[8]);
					mypars.alpha[25] = alphai(k[7],k[5],kargs[9]);
					mypars.alpha[26] = alphai(k[6],k[0],x[6]);
					mypars.alpha[27] = alphai(k[0],k[6],x[6]);
					mypars.alpha[28] = alphai(k[4],k[6],kargs[11]);
					mypars.alpha[29] = alphai(k[6],k[4],kargs[11]);
					mypars.alpha[30] = alphai(k[4],k[2],kargs[12]);
					mypars.alpha[31] = alphai(k[0],kargs[1],kargs[14]);
					mypars.alpha[32] = alphai(k[6],k[7],kargs[13]);
					mypars.alpha[33] = alphai(k[2],k[4],kargs[12]);
					mypars.alpha[34] = alphai(kargs[1],k[0],kargs[14]);
					mypars.alpha[35] = alphai(k[7],k[6],kargs[13]);
					mypars.alpha[36] = alphai(k[6],k[1],kargs[15]);
					mypars.alpha[37] = alphai(k[1],k[6],kargs[15]);
					mypars.alpha[38] = alphai(k[7],k[0],x[7]);
					mypars.alpha[39] = alphai(k[0],k[7],x[7]);
					mypars.alpha[40] = alphai(k[7],k[3],kargs[17]);
					mypars.alpha[41] = alphai(k[3],k[7],kargs[17]);
					mypars.alpha[42] = alphai(k[3],kargs[2],kargs[18]);
					mypars.alpha[43] = alphai(kargs[2],k[3],kargs[18]);
					mypars.alpha[44] = alphai(k[7],k[1],kargs[19]);
					mypars.alpha[45] = alphai(k[1],k[7],kargs[19]);
					mypars.alpha[46] = alphai(k[1],kargs[2],kargs[20]);
					mypars.alpha[47] = alphai(kargs[2],k[1],kargs[20]);
					mypars.alpha[48] = alphai(k[7],k[0],-x[7]);
					mypars.alpha[49] = alphai(k[0],k[7],-x[7]);
					mypars.alpha[50] = alphai(k[6],k[0],-x[6]);
					mypars.alpha[51] = alphai(k[0],k[6],-x[6]);
					mypars.alpha[52] = alphai(k[5],k[0],-x[5]);
					mypars.alpha[53] = alphai(k[0],k[5],-x[5]);
					mypars.alpha[54] = alphai(k[6],kargs[1],kargs[10]);
					mypars.alpha[55] = alphai(k[5],k[2],kargs[21]);
					mypars.alpha[56] = alphai(kargs[1],k[6],kargs[10]);
					mypars.alpha[57] = alphai(k[2],k[5],kargs[21]);
					mypars.alpha[58] = 1.;//alphai(k[5],1.+XMIN,0.);
					mypars.alpha[59] = 1.;//alphai(1.+XMIN,k[5],0.);
					mypars.alpha[60] = alphai(k[6],k[0]*sqrt(2.*(1.+XMIN)),0.);
					mypars.alpha[61] = alphai(k[0]*sqrt(2.*(1.+XMIN)),k[6],0.);
					mypars.alpha[62] = 1.;//alphai(1.-XMAX,k[7],0.);
					mypars.alpha[63] = 1.; //alphai(k[7],1.-XMAX,0.);
					mypars.alpha[64] = alphai(k[5],kargs[2], kargs[24]); // a(k[0],k[2]+p)
					mypars.alpha[65] = alphai(k[7],kargs[1],kargs[22]); // a(k[2],k1+p)
					mypars.alpha[66] = alphai(kargs[2],k[5],kargs[24]); // a(k[2]+p,k1)
					mypars.alpha[67] = alphai(kargs[1],k[7],kargs[22]); // a(k1+p,k[2])
					mypars.alpha[68] = alphai(kargs[2],k[1],kargs[20]);
					mypars.alpha[69] = alphai(k[1],kargs[2],kargs[20]);
					mypars.alpha[70] = alphai(k[0],k[2],x[2]);
					mypars.alpha[71] = alphai(k[2],k[0],x[2]);
					mypars.alpha[72] = alphai(k[6],kargs[2],-kargs[16]); // a(k[1],k[2]+p)
					mypars.alpha[73] = alphai(kargs[2],k[6],-kargs[16]);
					mypars.alpha[74] = alphai(k[7],k[2],kargs[23]); // a(k[2],k[1]+p)
					mypars.alpha[75] = alphai(k[2],k[7],kargs[23]);

          mypars.beta[0] = betai(k[1],k[2],x[0]);
					mypars.beta[1] = betai(k[0],k[2],x[2]);
					mypars.beta[2] = betai(k[0],k[1],x[1]);
					mypars.beta[3] = betai(k[0],k[3],x[3]);
					mypars.beta[4] = betai(k[0],k[4],x[4]);
					mypars.beta[5] = betai(k[0],k[5],x[5]);
					mypars.beta[6] = betai(k[3],k[5],kargs[3]);
					mypars.beta[7] = betai(k[3],kargs[1],kargs[4]);
					mypars.beta[8] = betai(k[0],kargs[2],kargs[6]);
					mypars.beta[9] = betai(k[6],k[5],kargs[5]);
					mypars.beta[10] = betai(k[4],k[5],kargs[7]);
					mypars.beta[11] = betai(k[4],kargs[1],kargs[8]);
					mypars.beta[12] = betai(k[7],k[5],kargs[9]);
					mypars.beta[13] = betai(k[0],k[6],x[6]);
					mypars.beta[14] = betai(k[4],k[6],kargs[11]);
					mypars.beta[15] = betai(k[4],k[2],kargs[12]);
					mypars.beta[16] = betai(k[0],kargs[1],kargs[14]);
					mypars.beta[17] = betai(k[7],k[6],kargs[13]);
					mypars.beta[18] = betai(k[1],k[6],kargs[15]);
					mypars.beta[19] = betai(k[0],k[7],x[7]);
					mypars.beta[20] = betai(k[3],k[7],kargs[17]);
					mypars.beta[21] = betai(k[3],kargs[2],kargs[18]);
					mypars.beta[22] = betai(k[1],k[7],kargs[19]);
					mypars.beta[23] = betai(k[1],kargs[2],kargs[20]);
					mypars.beta[24] = 0.;// betai(1.-XMAX,k[7],0.);
					mypars.beta[25] = betai(k[0],k[7],-x[7]);
					mypars.beta[26] = betai(k[0],k[6],-x[6]);
					mypars.beta[27] = betai(k[0],k[5],-x[5]);
					mypars.beta[28] = betai(k[6],kargs[1],kargs[10]);
					mypars.beta[29] = betai(k[5],k[2],kargs[21]);
					mypars.beta[31] = betai(k[0]*sqrt(2.*(1.+XMIN)),k[6],0.);
					mypars.beta[32] = betai(kargs[2],k[5],kargs[24]); // a(k[2]+p,k1)
					mypars.beta[33] = betai(k[7],kargs[1],kargs[22]); // a(k[2]+p,k[2])
					mypars.beta[34] = betai(k[2],k[0],x[2]);
					mypars.beta[35] = betai(kargs[2],k[1],kargs[20]);
					mypars.beta[36] = betai(k[6],kargs[2],-kargs[16]);
					mypars.beta[37] = betai(k[7],k[2],kargs[23]);

					// precompute gamma k-dependencies
					for(int i=0;i<8;i++){
						mypars.gamk[i]=1.-pow2(x[i]);
					}
					for(int i=8;i<30;i++){
						mypars.gamk[i] = 1.-pow2(kargs[i-5]);
					}

					// cosmo and beyond lcdm params
					mypars.omega0 = omega0;

					for (int i=0; i<maxpars;i++){
						mypars.extpars[i] = 	extpars[i];
					}


				gsl_odeiv2_system sys = {funcdgp, jacb, 94, &mypars};

			//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
				gsl_odeiv2_driver * d =
				gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
													epars[1], epars[2] ,epars[3]);
// smallest possible accuracy and initial step when comparing to analytic result
// must reduce initial step size for earlier times!
// rk8pd is optimal
				int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */

// For tree level initialisation

  /*1st order */

			//F1(k;a), G1(k;a)
			F1b_k1 = G[0] ;
			G1b_k1 = G[1] ;

			// F1(k3;a), G1(k3;a)
			F1b_k3 = G[2];
			G1b_k3 = G[3];

			//F1(k2;a), G1(k2;a)
			F1b_k2 = G[4];
			G1b_k2 = G[5];


// for B321 and B222 initialisaiton

      /* 1st order */

      //F1(k1;a), G1(k1;a)
			F1b_1 = G[16] ;
			G1b_1 = G[17] ;

			// F1(k2;a), G1(k2;a)
			F1b_2 = G[18];
			G1b_2 = G[19];

			//F1(k3;a), G1(k3;a)
			F1b_3 = G[20];
			G1b_3 = G[21];


			/*2nd order*/

			//F2/G2(k3,k2) for tree,  F2/G2(k1-p,k2+p) for B222
			F2b_k23 =  G[6];
			G2b_k23 =  G[7];

      //F2/G2(k1,k3) for tree, F2/G2(-p,k2+p) for B222
      F2b_k13 =  G[8];
      G2b_k13 =  G[9];

      //F2/G2(k1,k2) for tree, F2/G2(p,k1-p)
      F2b_k12 =  G[10];
      G2b_k12 =  G[11];


      // F2/G2(p,k2-p)
      F2b_p2mp = G[22];
      G2b_p2mp = G[23];

      // F2/G2(p,k3-p)
      F2b_p3mp = G[24];
      G2b_p3mp = G[25];

      // F3/G3(k1,p,k2-p)
      F3b_12mp = G[30];
      G3b_12mp = G[31];

      // F3/G3(k1,p,k3-p)
      F3b_13mp = G[34];
      G3b_13mp = G[35];

      // F3/G3(k2,p,k3-p)
      F3b_23mp = G[40];
      G3b_23mp = G[41];

      // F3/G3(k2,p,k1-p)
      F3b_21mp = G[44];
      G3b_21mp = G[45];


      // F3/G3(k3,p,k2-p)
      F3b_32mp = G[50];
      G3b_32mp = G[51];

      // F3/G3(k3,p,k1-p)
      F3b_31mp = G[54];
      G3b_31mp = G[55];

      // F3/G3(k1,k2) for B321-II
      F2b_12a = G[56];
      G2b_12a = G[57];

      // F3/G3(k1,k3) for B321-II
      F2b_13a = G[58];
      G2b_13a = G[59];

      // F3/G3(k2,k3) for B321-II
      F2b_23a = G[60];
      G2b_23a = G[61];

      //F3/G3(k1,p,-p)
      F3b_1pp = G[74];
      G3b_1pp = G[75];

      //F3/G3(k2,p,-p)
      F3b_2pp = G[76];
      G3b_2pp = G[77];

      //F3/G3(k3,p,-p)
      F3b_3pp = G[82];
      G3b_3pp = G[83];

      // F4/G4(-k1,-k2,p,-p)

      F4b_12pp = G[78];
      G4b_12pp = G[79];

      // F4/G4(-k1,-k3,p,-p)

      F4b_13pp = G[86];
      G4b_13pp = G[87];

      // F4/G4(-k2,-k3,p,-p)

      F4b_23pp = G[92];
      G4b_23pp = G[93];

			gsl_odeiv2_driver_free(d);
}



int funcfr(double a, const double G[], double F[], void *params)
{
	param_type5 p = *(param_type5 *)(params);

	/* Background quantities - LCDM  */
	double hub = HA2(a,p.omega0);
	double hub1= HA(a,p.omega0);
	double hubsqr = pow2(hub1);
	double asqr = pow2(a);
	double acub = asqr*a;

	/* fR Poisson modifications up to 4th order */
	double epsi = (p.omega0+4.*acub*(1.-p.omega0))/acub;
	double epsi2 = pow2(epsi);
	double epsi3 = epsi*epsi2;
	double epsi5 = epsi2*epsi3;
	double epsi6 = epsi3*epsi3;
	double epsi7 = epsi2*epsi5;
	double epsi9 = epsi2*epsi7;

	double p1p = p.extpars[0]/myh0sqr * pow2(3.*p.omega0-4.);
	double p1p2 = pow2(p1p);
	double p1p3 = p1p*p1p2;
	double p1p4 = p1p2*p1p2;

	double om0a2 = pow2(p.omega0/acub);
	double om0a3 = om0a2*p.omega0/acub;;
	double om0a4 = pow2(om0a2);

	double term1 = epsi3/2./p1p;

	double sqrsa[11];
	double pifs[11];
	double muak[11];

	//  mu and PI
	for(int i =0; i<11; i++){
		sqrsa[i] = p.gamk[i]/asqr;
		pifs[i] = sqrsa[i] + term1;
		muak[i] = 1.+ sqrsa[i]/ 3. / pifs[i];
		// redefine to reduce computations in gamma functions (sqrsa always divided by pifs)
		sqrsa[i] = p.gamk[i]/asqr / pifs[i];
	}

	double gam2 = -3./16./hubsqr*om0a2*epsi5/p1p2;
	double gam3a = -5./32./hubsqr*om0a3*epsi7/p1p3;
	double gam3b = -9./10.*gam3a*epsi3/p1p;

	double gam4a = -35./256/hubsqr*om0a4*epsi9/p1p4;
	double gam4b = gam4a*27./140.*epsi6/p1p2;
	double gam4c = gam4a*9./7.*epsi3/p1p;
	double gam4d = gam4b*6.;

  double gamk2[36];
  //gamma2, gamma3 and gamma4 scale dependencies  (sorry but need to put them by hand for optimisation)
	// arguments commented on right

	// optimise computations - lots of repetitions and scale dep parts can be pulled out!
	gamk2[21] = sqrsa[1]/(pifs[6]*pifs[10]); // k2.k3+p
	gamk2[28] = sqrsa[1]/(pifs[7]*pifs[2]);  //k3.k2+p
	gamk2[30] = sqrsa[1]/(pifs[0]*pifs[5]); // -p.k1

	gamk2[6] = sqrsa[2]/(pifs[0]*pifs[6]); // p.k2
	gamk2[12] = sqrsa[2]/(pifs[5]*pifs[4]); // k1.k3-p
	gamk2[24] = sqrsa[2]/(pifs[7]*pifs[1]); // k3.k1-p

	gamk2[27] = sqrsa[3]/(pifs[7]*pifs[9]);  // k3.k1+p
	gamk2[29] = sqrsa[3]/(pifs[5]*pifs[10]);  //k1.k3+p
	gamk2[31] = sqrsa[3]/(pifs[0]*pifs[6]); // -p.k2

	gamk2[15] = sqrsa[4]/(pifs[6]*pifs[9]); // k2.k1+p
	gamk2[26] = sqrsa[4]/(pifs[5]*pifs[2]); // k1.k2+p
	gamk2[32] = sqrsa[4]/(pifs[0]*pifs[7]); // -p.k3

	gamk2[1] = sqrsa[5]/(pifs[0]*pifs[1]); // p.k1-p
	gamk2[17] = sqrsa[5]/(pifs[2]*pifs[4]);// k2+p.k3-p
	gamk2[18] = sqrsa[5]/(pifs[6]*pifs[7]);// k2.k3
	gamk2[19] = sqrsa[5]/(pifs[0]*pifs[9]);// -p.k1+p
	gamk2[23] = sqrsa[5]/(pifs[10]*pifs[3]);  // k3+p.k2-p
	gamk2[33] = sqrsa[5]/(pifs[5]*term1); // k1,0

	gamk2[2] = sqrsa[6]/(pifs[0]*pifs[2]); // -p.k2+p
	gamk2[3] = sqrsa[6]/(pifs[0]*pifs[3]); // p. k2-p
	gamk2[13] = sqrsa[6]/(pifs[9]*pifs[4]);// k1+p.k3-p
	gamk2[14] = sqrsa[6]/(pifs[5]*pifs[7]); // k1.k3
	gamk2[25] = sqrsa[6]/(pifs[10]*pifs[3]);  // k3+p.k1-p
	gamk2[34] = sqrsa[6]/(pifs[6]*term1); // k2,0

	gamk2[0] = sqrsa[7]/(pifs[1]*pifs[2]); //k1-p.k2+p
	gamk2[4] = sqrsa[7]/(pifs[0]*pifs[4]); // p.k3-p
	gamk2[9] = sqrsa[7]/(pifs[9]*pifs[3]); // p+k1.k2-p
	gamk2[10] = sqrsa[7]/(pifs[5]*pifs[6]); // k1.k2
	gamk2[11] = sqrsa[7]/(pifs[0]*pifs[10]); // -p.k3+p
	gamk2[35] = sqrsa[7]/(pifs[7]*term1); // k3,0

	gamk2[5] = sqrsa[9]/(pifs[0]*pifs[5]); // p.k1
	gamk2[22] = sqrsa[9]/(pifs[7]*pifs[3]); // k3.k2-p
	gamk2[16] = sqrsa[9]/(pifs[6]*pifs[4]);// k2.k3-p

	gamk2[7] = sqrsa[10]/(pifs[0]*pifs[7]); // p.k3
	gamk2[8] = sqrsa[10]/(pifs[5]*pifs[3]); // k1.k2-p
	gamk2[20] = sqrsa[10]/(pifs[6]*pifs[1]); // k2.k1-p


	double gamk3[27];
	gamk3[13] = sqrsa[1]/(pifs[6]*pifs[7]*pifs[0]); // -k2,-k3,-p
	gamk3[10] = sqrsa[2]/(pifs[5]*pifs[7]*pifs[0]); //-k1,-k3,p
	gamk3[12] = sqrsa[3]/(pifs[5]*pifs[7]*pifs[0]); // -k1,-k3,-p
	gamk3[7] = sqrsa[4]/(pifs[5]*pifs[6]*pifs[0]); // -k1,-k2,-p


	gamk3[2] = sqrsa[5]/(pifs[6]*pifs[0]*pifs[4]); // k2,p,k3-p
	gamk3[4] = sqrsa[5]/(pifs[7]*pifs[0]*pifs[3]); // k3,p,k2-p
	gamk3[8] = sqrsa[5]/(pifs[5]*pifs[0]*pifs[0]); //-k1,p,-p
	gamk3[19] = sqrsa[5]/(pifs[10]*pifs[6]*pifs[0]); // -k3-p,-k2,p
	gamk3[20] = sqrsa[5]/(pifs[2]*pifs[7]*pifs[0]); // k2+p,-k3,p
	gamk3[25] = sqrsa[5]/(pifs[0]*pifs[0]*pifs[5]); // k1,p,-p
	gamk3[26] = sqrsa[5]/(term1*pifs[6]*pifs[7]); // 0.,-k2,-k3

	gamk3[1] = sqrsa[6]/(pifs[5]*pifs[0]*pifs[4]); // k1,p,k3-p
	gamk3[5] = sqrsa[6]/(pifs[7]*pifs[0]*pifs[1]); // k3,p,k1-p
	gamk3[17] = sqrsa[6]/(pifs[9]*pifs[7]*pifs[0]); // -k1-p,-k3,p
	gamk3[18] = sqrsa[6]/(pifs[10]*pifs[5]*pifs[0]); // -k3-p,-k1,p
	gamk3[23] = sqrsa[6]/(pifs[0]*pifs[0]*pifs[6]); // k2,p,-p
	gamk3[24] = sqrsa[6]/(term1*pifs[5]*pifs[7]); // 0.,-k1,-k3
	gamk3[9] = sqrsa[6]/(pifs[6]*pifs[0]*pifs[0]); // -k2,p,-p


	gamk3[0] = sqrsa[7]/(pifs[5]*pifs[0]*pifs[3]); // k1,p,k2-p
	gamk3[3] = sqrsa[7]/(pifs[1]*pifs[0]*pifs[6]); // k1-p,p,k2
	gamk3[11] = sqrsa[7]/(pifs[7]*pifs[0]*pifs[0]); // -k3,p,-p
	gamk3[15] = sqrsa[7]/(pifs[9]*pifs[6]*pifs[0]); // -k1-p,-k2,p
	gamk3[16] = sqrsa[7]/(pifs[2]*pifs[5]*pifs[0]);// -k2-p,-k1,p
	gamk3[21] = sqrsa[7]/(pifs[0]*pifs[0]*pifs[7]); // k3,p,-p
	gamk3[22] = sqrsa[7]/(term1*pifs[6]*pifs[5]); // 0.,-k1,-k2

	gamk3[14] = sqrsa[9]/(pifs[6]*pifs[7]*pifs[0]); // -k2,-k3,p
	gamk3[6] = sqrsa[10]/(pifs[5]*pifs[6]*pifs[0]); // -k1,-k2,p

	double gamk4[3];
	gamk4[0] = sqrsa[7]/(pifs[5]*pifs[6]*pifs[0]*pifs[0]); // -k1,-k2,p,-p
	gamk4[1] = sqrsa[6]/(pifs[5]*pifs[7]*pifs[0]*pifs[0]); // -k1,-k3,p,-p
	gamk4[2] = sqrsa[5]/(pifs[7]*pifs[6]*pifs[0]*pifs[0]); // -k2,-k3,p,-p


	/* 1st order */
	//1. F1/G1(p)
	F[0] = -G[1]/a;
	F[1] =1./a*(-(2.-hub)*G[1]-hub*G[0]*muak[0]);

	//2.  F1/G1(k2+p)
	F[2] = -G[3]/a;
	F[3] =1./a*(-(2.-hub)*G[3]-hub*G[2]*muak[2]);

	//3. F1/G1(k1-p)
	F[4] = -G[5]/a;
	F[5] =1./a*(-(2.-hub)*G[5]-hub*G[4]*muak[1]);

	/* 2nd order for B222 */

	//4. F2/G2(k1-p,k2+p)
	F[6] =1./a*(-(p.alpha[0]*G[5]*G[2]+p.alpha[1]*G[3]*G[4])/2.-G[7]);
	F[7] =1./a*(-(2.-hub)*G[7]-hub*G[6]*muak[7] - gam2*gamk2[0]*G[4]*G[2] - p.beta[0]*G[5]*G[3]);

	//5. F2/G2(-p,k2+p)
	F[8] =1./a*(-(p.alpha[2]*G[1]*G[2]+p.alpha[3]*G[3]*G[0])/2.-G[9]);
	F[9] =1./a*(-(2.-hub)*G[9]-hub*G[8]*muak[6] - gam2*gamk2[2]*G[2]*G[0] - p.beta[1]*G[3]*G[1]);

	//7. F2/G2(k2,k1) or F2/G2(k1-p,p)
	F[10] =1./a*(-(p.alpha[4]*G[5]*G[0]+p.alpha[5]*G[1]*G[4])/2.-G[11]) ;
	F[11] =1./a*(-(2.-hub)*G[11]-hub*G[10]*muak[5] - gam2*gamk2[1]*G[4]*G[0]-p.beta[2]*G[5]*G[1]);


  /* 1st order for B321 - I */

  //  F1/G1(k2-p)

  F[12] = -G[13]/a;
  F[13] =1./a*(-(2.-hub)*G[13]-hub*G[12]*muak[3]);

  // F1/G1(k3-p)

  F[14] = -G[15]/a;
  F[15] =1./a*(-(2.-hub)*G[15]-hub*G[14]*muak[4]);

  //  F1/G1(k1)

  F[16] = -G[17]/a;
  F[17] =1./a*(-(2.-hub)*G[17]-hub*G[16]*muak[5]);

  //  F1/G1(k2)

  F[18] = -G[19]/a;
  F[19] =1./a*(-(2.-hub)*G[19]-hub*G[18]*muak[6]);

  //  F1/G1(k3)

  F[20] = -G[21]/a;
  F[21] =1./a*(-(2.-hub)*G[21]-hub*G[20]*muak[7]);

  // 2nd order used in 3rd order equations of B321-I


  // F2/G2(k2-p,p)
  F[22] =1./a*(-(p.alpha[6]*G[13]*G[0]+p.alpha[7]*G[1]*G[12])/2.-G[23]) ;
  F[23] =1./a*(-(2.-hub)*G[23]-hub*G[22]*muak[6] - gam2*gamk2[3]*G[12]*G[0]-p.beta[3]*G[13]*G[1]);

  // F2/G2(k5,k1) = (k3-p,p)
  F[24] =1./a*(-(p.alpha[8]*G[15]*G[0]+p.alpha[9]*G[1]*G[14])/2.-G[25]) ;
  F[25] =1./a*(-(2.-hub)*G[25]-hub*G[24]*muak[7] - gam2*gamk2[4]*G[14]*G[0]-p.beta[4]*G[15]*G[1]);


  // F2/G2(k6,k1) = (k1,p)
  F[26] =1./a*(-(p.alpha[10]*G[17]*G[0]+p.alpha[11]*G[1]*G[16])/2.-G[27]) ;
  F[27] =1./a*(-(2.-hub)*G[27]-hub*G[26]*muak[9] - gam2*gamk2[5]*G[16]*G[0]-p.beta[5]*G[17]*G[1]);


  // F2/G2(k6,k4) = (k1,k2-p)
  F[28] =1./a*(-(p.alpha[12]*G[17]*G[12]+p.alpha[13]*G[13]*G[16])/2.-G[29]) ;
  F[29] =1./a*(-(2.-hub)*G[29]-hub*G[28]*muak[10] - gam2*gamk2[8]*G[16]*G[12]-p.beta[6]*G[17]*G[13]);




  // F3/G3(k6,k1,k4) = (k1,p,k2-p)
 F[30] = - 1./(3.*a)*(p.alpha[14]*G[26]*G[13]
           + p.alpha[15]*G[28]*G[1]
           + p.alpha[16]*G[22]*G[17]

           + p.alpha[17]*G[27]*G[12]
           + p.alpha[18]*G[29]*G[0]
           + p.alpha[19]*G[23]*G[16]

           +3.*G[31]) ;


 F[31] =1./(3.*a)*(-3.*(2.-hub)*G[31]-3.*hub*G[30]*muak[7]

          -2.*p.beta[9]*G[17]*G[23]

          -2.*p.beta[8]*G[1]*G[29]

          -2.*p.beta[7]*G[13]*G[27]

          -2.*gam2*gamk2[10]*G[16]*G[22] // gam2(k1,p + k2- p )
          -2.*gam2*gamk2[11]*G[0]*G[28] // gam2(p, k1 + k2-p )
          -2.*gam2*gamk2[9]*G[12]*G[26] // gam2(k2-p, k1 + p)

					-gamk3[0]*(3.*gam3a + gam3b*(1./pifs[10] + 1./pifs[9] + 1./pifs[6]))*G[12]*G[16]*G[0]);


  // F2/G2(k3-p,k1)
  F[32] =1./a*(-(p.alpha[20]*G[15]*G[16]+p.alpha[21]*G[14]*G[17])/2.-G[33]) ;
  F[33] =1./a*(-(2.-hub)*G[33]-hub*G[32]*muak[2] - gam2*gamk2[12]*G[16]*G[14]-p.beta[10]*G[17]*G[15]);


// F3/G3(k6,k1,k5) = (k1,p,k3-p)
 F[34] = - 1./(3.*a)*(p.alpha[22]*G[26]*G[15]
           + p.alpha[2]*G[32]*G[1]
           + p.alpha[23]*G[24]*G[17]

           + p.alpha[24]*G[27]*G[14]
           + p.alpha[3]*G[33]*G[0]
           + p.alpha[25]*G[25]*G[16]

           +3.*G[35]) ;


 F[35] =1./(3.*a)*(-3.*(2.-hub)*G[35]-3.*hub*G[34]*muak[6]

          -2.*p.beta[12]*G[17]*G[25]

          -2.*p.beta[1]*G[1]*G[33]

          -2.*p.beta[11]*G[15]*G[27]

          -2.*gam2*gamk2[14]*G[16]*G[24] // gam2(k1, p + k3-p)
          -2.*gam2*gamk2[2]*G[0]*G[32] // gam2(p, k3-p + k1)
          -2.*gam2*gamk2[13]*G[14]*G[26] // gam2( k3-p, p + k1)

				-gamk3[1]*(3.*gam3a + gam3b*(1./pifs[9] + 1./pifs[7] + 1./pifs[2]))*G[14]*G[16]*G[0]);


 // F2/G2(k7,k1) = (k2,p)
 F[36] =1./a*(-(p.alpha[26]*G[19]*G[0]+p.alpha[27]*G[1]*G[18])/2.-G[37]) ;
 F[37] =1./a*(-(2.-hub)*G[37]-hub*G[36]*muak[2] - gam2*gamk2[6]*G[18]*G[0]-p.beta[13]*G[19]*G[1]);


 // F2/G2(k5,k7) = (k3-p,k2)
 F[38] =1./a*(-(p.alpha[28]*G[15]*G[18]+p.alpha[29]*G[19]*G[14])/2.-G[39]) ;
 F[39] =1./a*(-(2.-hub)*G[39]-hub*G[38]*muak[9] - gam2*gamk2[16]*G[14]*G[18]-p.beta[14]*G[15]*G[19]);


 // F3/G3(k7,k1,k5) = (k2,p,k3-p)

  F[40] = - 1./(3.*a)*(p.alpha[30]*G[36]*G[15]
            + p.alpha[31]*G[38]*G[1]
            + p.alpha[32]*G[24]*G[19]

            + p.alpha[33]*G[37]*G[14]
            + p.alpha[34]*G[39]*G[0]
            + p.alpha[35]*G[25]*G[18]

            +3.*G[41]) ;



  F[41] =1./(3.*a)*(-3.*(2.-hub)*G[41]-3.*hub*G[40]*muak[5]

           -2.*p.beta[17]*G[19]*G[25]

           -2.*p.beta[16]*G[1]*G[39]

           -2.*p.beta[15]*G[15]*G[37]

           -2.*gam2*gamk2[18]*G[18]*G[24] // gam2(k2, p + k3-p)
           -2.*gam2*gamk2[19]*G[0]*G[38]  // gam2(p, k3-p + k2)
           -2.*gam2*gamk2[17]*G[14]*G[36]  // gam2(k3-p, k2 + p)

				-gamk3[2]*(3.*gam3a + gam3b*(1./pifs[2] + 1./pifs[7] + 1./pifs[9]))*G[14]*G[18]*G[0]);


 // F2/G2(k2,k7) = (k1-p,k2)
 F[42] =1./a*(-(p.alpha[36]*G[19]*G[2]+p.alpha[37]*G[3]*G[18])/2.-G[43]) ;
 F[43] =1./a*(-(2.-hub)*G[43]-hub*G[42]*muak[10] - gam2*gamk2[20]*G[18]*G[2]-p.beta[18]*G[19]*G[3]);


//F3/G3(k1-p,p,k2)

 F[44] = - 1./(3.*a)*(p.alpha[0]*G[36]*G[5]
           + p.alpha[15]*G[42]*G[1]
           + p.alpha[19]*G[10]*G[19]

           + p.alpha[1]*G[37]*G[4]
           + p.alpha[18]*G[43]*G[0]
           + p.alpha[16]*G[11]*G[18]

           +3.*G[45]) ;


 F[45] =1./(3.*a)*(-3.*(2.-hub)*G[45]-3.*hub*G[44]*muak[7]

          -2.*p.beta[9]*G[19]*G[11]

          -2.*p.beta[8]*G[1]*G[43]

          -2.*p.beta[0]*G[5]*G[37]


          -2*gam2*gamk2[0]*G[4]*G[36] // gam2(k1-p, k2, p  )
          -2*gam2*gamk2[11]*G[0]*G[42] // gam2(p, k1-p + k2 )
          -2*gam2*gamk2[10]*G[18]*G[10] // gam2(k2, k1-p + p  )

				-gamk3[3]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[10] + 1./pifs[2]))*G[4]*G[18]*G[0]);


 // F2/G2(k8,k1) = (k3,p)
 F[46] =1./a*(-(p.alpha[38]*G[21]*G[0]+p.alpha[39]*G[1]*G[20])/2.-G[47]) ;
 F[47] =1./a*(-(2.-hub)*G[47]-hub*G[46]*muak[10] - gam2*gamk2[7]*G[20]*G[0]-p.beta[19]*G[21]*G[1]);

 // F2/G2(k4,k8) = (k2-p,k3)
 F[48] =1./a*(-(p.alpha[40]*G[21]*G[12]+p.alpha[41]*G[13]*G[20])/2.-G[49]) ;
 F[49] =1./a*(-(2.-hub)*G[49]-hub*G[48]*muak[9] - gam2*gamk2[22]*G[20]*G[12]-p.beta[20]*G[21]*G[13]);



 // F3/G3(k8,k1,k4) = (k3,p,k2-p)
 F[50] = - 1./(3.*a)*(
             p.alpha[42]*G[46]*G[13]
           + p.alpha[31]*G[48]*G[1]
           + p.alpha[35]*G[22]*G[21]

           + p.alpha[43]*G[47]*G[12]
           + p.alpha[34]*G[49]*G[0]
           + p.alpha[32]*G[23]*G[20]

           +3.*G[51]) ;


 F[51] =1./(3.*a)*(-3.*(2.-hub)*G[51]-3.*hub*G[50]*muak[5]

          -2.*p.beta[17]*G[21]*G[23]

          -2.*p.beta[16]*G[1]*G[49]

          -2.*p.beta[21]*G[13]*G[47]

          -2.*gam2*gamk2[18]*G[20]*G[22] // gam2(k3, p + k2-p)
          -2.*gam2*gamk2[19]*G[0]*G[48] // gam2(p, k3 + k2-p)
          -2.*gam2*gamk2[23]*G[12]*G[46] // gam2(k2-p, k3 + p)

					-gamk3[4]*(3.*gam3a + gam3b*(1./pifs[10] + 1./pifs[6] + 1./pifs[9]))*G[12]*G[20]*G[0]);

 // F2/G2(k2,k8) = (k1-p,k3)

 F[52] =1./a*(-(p.alpha[44]*G[21]*G[2]+p.alpha[45]*G[3]*G[20])/2.-G[53]) ;
 F[53] =1./a*(-(2.-hub)*G[53]-hub*G[52]*muak[2] - gam2*gamk2[24]*G[20]*G[2]-p.beta[22]*G[21]*G[3]);


// F3/G3(k3,p,k1-p)
 F[54] = - 1./(3.*a)*(p.alpha[46]*G[46]*G[3]
           + p.alpha[2]*G[52]*G[1]
           + p.alpha[25]*G[10]*G[21]

           + p.alpha[47]*G[47]*G[2]
           + p.alpha[3]*G[53]*G[0]
           + p.alpha[23]*G[11]*G[20]

           +3.*G[55]) ;


 F[55] =1./(3.*a)*(-3.*(2.-hub)*G[55]-3.*hub*G[54]*muak[6]

          -2.*p.beta[12]*G[21]*G[11]

          -2.*p.beta[1]*G[1]*G[53]

          -2.*p.beta[23]*G[3]*G[47]

          -2*gam2*gamk2[14]*G[20]*G[10] // gam2(k3, p + k1-p)
          -2*gam2*gamk2[2]*G[0]*G[52] // gam2(p, k3 + k1-p)
          -2*gam2*gamk2[25]*G[2]*G[46] // gam2(k1-p, k3 + p)


				-gamk3[5]*(3.*gam3a + gam3b*(1./pifs[10] + 1./pifs[2] + 1./pifs[5]))*G[2]*G[20]*G[0]);



// F2/G2(k6,k7) = F2/G2(k1,k2)

F[56] =1./a*(-(p.alpha[16]*G[17]*G[18]+p.alpha[19]*G[19]*G[16])/2.-G[57]) ;
F[57] =1./a*(-(2.-hub)*G[57]-hub*G[56]*muak[7] - gam2*gamk2[10]*G[18]*G[16]-p.beta[9]*G[17]*G[19]);


// F2/G2(k6,k8) = F2/G2(k1,k3)

F[58] =1./a*(-(p.alpha[23]*G[17]*G[20]+p.alpha[25]*G[21]*G[16])/2.-G[59]) ;
F[59] =1./a*(-(2.-hub)*G[59]-hub*G[58]*muak[6] - gam2*gamk2[14]*G[20]*G[16]-p.beta[12]*G[17]*G[21]);


// F2/G2(k7,k8) = F2/G2(k2,k3)

F[60] =1./a*(-(p.alpha[32]*G[19]*G[20]+p.alpha[35]*G[21]*G[18])/2.-G[61]) ;
F[61] =1./a*(-(2.-hub)*G[61]-hub*G[60]*muak[5] - gam2*gamk2[18]*G[20]*G[18]-p.beta[17]*G[19]*G[21]);


// F2/G2(-k3,p)
F[62] =1./a*(-(p.alpha[48]*G[21]*G[0]+p.alpha[49]*G[1]*G[20])/2.-G[63]) ;
F[63] =1./a*(-(2.-hub)*G[63]-hub*G[62]*muak[4] - gam2*gamk2[32]*G[20]*G[0]-p.beta[25]*G[21]*G[1]);


 // F2/G2(-k2,p)
 F[64] =1./a*(-(p.alpha[50]*G[19]*G[0]+p.alpha[51]*G[1]*G[18])/2.-G[65]) ;
 F[65] =1./a*(-(2.-hub)*G[65]-hub*G[64]*muak[3]- gam2*gamk2[31]*G[18]*G[0]-p.beta[26]*G[19]*G[1]);


 // F2/G2(k6,k1) = (-k1,p)
 F[66] =1./a*(-(p.alpha[52]*G[17]*G[0]+p.alpha[53]*G[1]*G[16])/2.-G[67]) ;
 F[67] =1./a*(-(2.-hub)*G[67]-hub*G[66]*muak[1] - gam2*gamk2[30]*G[16]*G[0]-p.beta[27]*G[17]*G[1]);


 // F3/G3(k6,k7,k1) = (-k1,-k2,p)

 F[68] = - 1./(3.*a)*(p.alpha[39]*G[56]*G[1]
                    + p.alpha[36]*G[66]*G[19]
                    + p.alpha[12]*G[64]*G[17]

           + p.alpha[38]*G[57]*G[0]
           + p.alpha[37]*G[67]*G[18]
           + p.alpha[13]*G[65]*G[16]

           +3.*G[69]) ;



 F[69] =1./(3.*a)*(-3.*(2.-hub)*G[69]-3.*hub*G[68]*muak[10]

          -2.*p.beta[6]*G[65]*G[17]

          -2.*p.beta[18]*G[67]*G[19]

          -2.*p.beta[19]*G[57]*G[1]

          -2.*gam2*gamk2[8]*G[16]*G[64]
          -2.*gam2*gamk2[20]*G[18]*G[66]
          -2.*gam2*gamk2[7]*G[0]*G[56]

					-gamk3[6]*(3.*gam3a + gam3b*(1./pifs[7] + 1./pifs[1] + 1./pifs[3]))*G[0]*G[18]*G[16]);


// F3/G3(k6,k7,k1) = (-k1,-k2,-p)
F[70] = - 1./(3.*a)*(p.alpha[49]*G[56]*G[1]
                  + p.alpha[54]*G[26]*G[19]
                  + p.alpha[55]*G[36]*G[17]

         + p.alpha[48]*G[57]*G[0]
         + p.alpha[56]*G[27]*G[18]
         + p.alpha[57]*G[37]*G[16]

         +3.*G[71]) ;


F[71] =1./(3.*a)*(-3.*(2.-hub)*G[71]-3.*hub*G[70]*muak[4]

        -2.*p.beta[29]*G[37]*G[17]

        -2.*p.beta[28]*G[27]*G[19]

        -2.*p.beta[25]*G[57]*G[1]

        -2.*gam2*gamk2[26]*G[16]*G[36]
        -2.*gam2*gamk2[15]*G[18]*G[26]
        -2.*gam2*gamk2[7]*G[0]*G[56]

				-gamk3[7]*(3.*gam3a + gam3b*(1./pifs[7] + 1./pifs[9] + 1./pifs[2]))*G[0]*G[18]*G[16]);


F[72] =1./a*(-G[73]) ; // set 1st term to 0 because x1 is not identically -1 (needed to keep evolution at 0 for GR)
F[73] =1./a*(-(2.-hub)*G[73]-hub*G[72] - betai(p.kv[0],p.kv[0],XMIN)*G[1]*G[1]);


// F3/G3(k6,k1,k1) = (-k1,p,-p)
F[74] = - 1./(3.*a)*(p.alpha[5]*G[66]*G[1]
                  + p.alpha[31]*G[26]*G[1]
                  + p.alpha[58]*G[72]*G[17]

         				+ p.alpha[4]*G[67]*G[0]
         				+ p.alpha[34]*G[27]*G[0]
         				+ p.alpha[59]*G[73]*G[16]

         				+3.*G[75]) ;


F[75] =1./(3.*a)*(-3.*(2.-hub)*G[75]-3.*hub*G[74]*muak[5]

        -2.*p.beta[30]*G[73]*G[17]

        -2.*p.beta[16]*G[27]*G[1]

        -2.*p.beta[2]*G[67]*G[1]

        -2*gam2*gamk2[33]*G[16]*G[72]  //gam2 (k1,0)
        -2*gam2*gamk2[19]*G[0]*G[26]
        -2*gam2*gamk2[1]*G[0]*G[66]

				-gamk3[8]*(3.*gam3a + gam3b*(1./pifs[1] + 1./pifs[9] + 1./term1))*G[0]*G[0]*G[16]);



// F3/G3(k7,k1,k1) = (-k2,p,-p)
F[76] = - 1./(3.*a)*(p.alpha[7]*G[64]*G[1]
                   + p.alpha[2]*G[36]*G[1]
                   + p.alpha[60]*G[72]*G[19]

                  + p.alpha[6]*G[65]*G[0]
                  + p.alpha[3]*G[37]*G[0]
                  + p.alpha[61]*G[73]*G[18]

                  +3.*G[77]) ;


F[77] =1./(3.*a)*(-3.*(2.-hub)*G[77]-3.*hub*G[76]*muak[6]

        -2.*p.beta[3]*G[65]*G[1]

        -2.*p.beta[1]*G[37]*G[1]

        -2.*p.beta[31]*G[73]*G[19]

        -2.*gam2*gamk2[34]*G[18]*G[72]
        -2.*gam2*gamk2[2]*G[0]*G[36]
        -2.*gam2*gamk2[3]*G[0]*G[64]

				-gamk3[9]*(3.*gam3a + gam3b*(1./pifs[3] + 1./pifs[2] + 1./term1))*G[0]*G[0]*G[18]);



// F4/G4(-k6,-k7,k1,-k1) = (-k1,-k2,p,-p)
F[78] = -1./(24.*a)*(
										6.*p.alpha[16]*G[76]*G[17] // a(-k1,-k2) F3(-k2,p,-p) G1(k1)
                   +6.*p.alpha[19]*G[74]*G[19] // a(-k2,-k1) F3(-k1,p,-p) G1(k2)
                   +6.*p.alpha[9]*G[70]*G[1] // a(p,k3-p) F3(-p,-k1,-k2) G1(p)
                   +6.*p.alpha[15]*G[68]*G[1] // a(-p,k3+p) F3(p,-k1,-k2) G1(p)

                    +6.*p.alpha[19]*G[77]*G[16] // a(-k2,-k1) G3(-k2,p,-p) F1(k1)
                    +6.*p.alpha[16]*G[75]*G[18] // a(-k2,-k1) G3(-k1,p,-p) F1(k2)
                    +6.*p.alpha[8]*G[71]*G[0] // a(p,k3-p) G3(-p,-k1,-k2) F1(p)
                    +6.*p.alpha[18]*G[69]*G[0] // a(p,k3-p) G3(-p,-k1,-k2) F1(p)

                    +4.*p.alpha[62]*G[56]*G[73] // a(0,k3) F2(k1,k2) G2(p,-p)
                    +4.*p.alpha[63]*G[72]*G[57] // a(k3,0) F2(p,-p) G2(k1,k2)

                    + 4.*p.alpha[17]*G[64]*G[27] // a(-k1-p,-k2+p,) F2(-k2,p) G2(-p,-k1)
                    + 4.*p.alpha[14]*G[26]*G[65] // a(-k2+p,-k1-p) F2(-k1,-p) G2(-k2,p)

                    + 4.*p.alpha[0]*G[36]*G[67] // a(-k1+p,-k2-p) F2(-k2,-p)G2(-k1,p)
                    + 4.*p.alpha[1]*G[66]*G[37] // a(-k2-p,-k1+p) F2(-k1,p) G2(-k2,-p)

                    +24.*G[79]);

F[79] = 1./a*(-(2.-hub)*G[79]-hub*G[78]*muak[7]

          -1./24.*(
                  12.*p.beta[9]*G[77]*G[17] // b(-k1,-k2) G3(-k2,p,-p) G1(k1)
                  +12.*p.beta[9]*G[75]*G[19] // b(-k2,-k1) G3(-k1,p,-p) G1(k2)
                  +12.*p.beta[4]*G[71]*G[1] // b(p,k3-p) G3(-p,-k1,-k2) G1(p)
                  +12.*p.beta[8]*G[69]*G[1] // b(-p,k3+p) G3(p,-k1,-k2) G1(p)

                  + 8.*p.beta[24]*G[57]*G[73] // b(k3,0) G2(k1,k2) G2(p,-p)
                  + 8.*p.beta[7]*G[65]*G[27]  // b(k2-p,k1+p) G2(-k2,p) G2(k1,p)
                  + 8.*p.beta[0]*G[67]*G[37] // b(k2+p,k1-p) G2(k2,p) G2(-k1,p)

//2nd order
                  + 8.*gam2*gamk2[35]*G[56]*G[72] // gam2(k3,0)
                  + 8.*gam2*gamk2[9]*G[64]*G[26]  // gam2(k2-p,k1+p)
                  + 8.*gam2*gamk2[0]*G[66]*G[36]  // gam2(k1-p,k2+p)


                  + 12.*gam2*gamk2[10]*G[76]*G[16] // gam2(-k2+p-p,-k1)
                  + 12.*gam2*gamk2[10]*G[74]*G[18] // gam2(-k1+p-p,-k2)
                  + 12.*gam2*gamk2[4]*G[70]*G[0] // gam2(k3-p,p)
                  + 12.*gam2*gamk2[11]*G[68]*G[0] // gam2(k3+p,-p)

//3rd order

							+4.*(gamk3[0]*(3.*gam3a + gam3b*(1./pifs[9] + 1./pifs[6] + 1./pifs[10]))*G[16]*G[0]*G[64]
									   + gamk3[3]*(3.*gam3a + gam3b*(1./pifs[2] + 1./pifs[5] + 1./pifs[10]))*G[18]*G[0]*G[26]
										 + gamk3[15]*(3.*gam3a + gam3b*(1./pifs[3] + 1./pifs[5] + 1./pifs[4]))*G[18]*G[0]*G[66]
										 + gamk3[16]*(3.*gam3a + gam3b*(1./pifs[1] + 1./pifs[6] + 1./pifs[4]))*G[16]*G[0]*G[36]
									   + gamk3[21]*(3.*gam3a + gam3b*(1./term1 + 1./pifs[10] + 1./pifs[4]))*G[0]*G[0]*G[56]
								     + gamk3[22]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[6] + 1./pifs[7]))*G[16]*G[18]*G[72])


//4th order
							+ gamk4[0]*(24.*gam4a + 8.*gam4b*(1./(pifs[7]*term1) + 1./(pifs[1]*pifs[2]) + 1./(pifs[9]*pifs[3]))
																		+ 4.*gam4c*(1./pifs[7] + 1./term1 + 1./pifs[1] + 1./pifs[2] + 1./pifs[9] + 1./pifs[3])
																		+ 4./3.*gam4d*(1./(pifs[4]*pifs[7]) + 1./(pifs[4]*pifs[2]) + 1./(pifs[4]*pifs[9])
																			 			   + 1./(pifs[10]*pifs[7]) + 1./(pifs[10]*pifs[1]) + 1./(pifs[10]*pifs[3])
																							 + 1./(pifs[5]*term1) + 1./(pifs[5]*pifs[1]) +	1./(pifs[5]*pifs[9])
																							 + 1./(pifs[6]*term1) + 1./(pifs[6]*pifs[2]) + 1./(pifs[6]*pifs[3]))
																		+ 12./3.*gam4c*(1./pifs[10] + 1./pifs[4] + 1./pifs[5]  + 1./pifs[6]))*G[16]*G[18]*G[0]*G[0]));


// F3/G3(k6,k8,k1) = (-k1,-k3,p)
F[80] = - 1./(3.*a)*(p.alpha[21]*G[62]*G[17] // a(k1,k3-p) F2(-k3,p) G1(k1)
                   + p.alpha[44]*G[66]*G[21] // a(k3,k1-p) F2(-k1,p) G1(k3)
                   + p.alpha[27]*G[58]*G[1] // a(p,k2) F2(k1,k3) G1(p)

          + p.alpha[20]*G[63]*G[16]
          + p.alpha[45]*G[67]*G[20]
          + p.alpha[26]*G[59]*G[0]

          +3.*G[81]) ;


F[81] =1./(3.*a)*(-3.*(2.-hub)*G[81]-3.*hub*G[80]*muak[2]

         -2.*p.beta[10]*G[63]*G[17]

         -2.*p.beta[22]*G[67]*G[21]

         -2.*p.beta[13]*G[59]*G[1]

         -2.*gam2*gamk2[12]*G[16]*G[62]
         -2.*gam2*gamk2[24]*G[20]*G[66]
         -2.*gam2*gamk2[6]*G[0]*G[58]

				 -gamk3[10]*(3.*gam3a + gam3b*(1./pifs[6] + 1./pifs[2] + 1./pifs[4]))*G[0]*G[20]*G[16]);



// F3/G3(k8,k1,k1) = (-k3,p,-p)
F[82] = - 1./(3.*a)*(p.alpha[9]*G[62]*G[1] // a(p,k3-p) F2(k3,-p) G1(p)
                 + p.alpha[15]*G[46]*G[1] // a(p,-k3-p) F2(k3,p) G1(p)
                 + p.alpha[63]*G[72]*G[20] // a(-k3,0.) F2(p,-p) G1(k3)

        + p.alpha[8]*G[63]*G[0]
        + p.alpha[18]*G[47]*G[0]
        + p.alpha[62]*G[73]*G[21]

        +3.*G[83]) ;


F[83] =1./(3.*a)*(-3.*(2.-hub)*G[83]-3.*hub*G[82]*muak[7]

       -2.*p.beta[4]*G[63]*G[1]

       -2.*p.beta[8]*G[47]*G[1]

       -2.*p.beta[24]*G[73]*G[21]


       -2*gam2*gamk2[35]*G[20]*G[72]
       -2.*gam2*gamk2[11]*G[0]*G[46]
       -2.*gam2*gamk2[4]*G[0]*G[62]

			 -gamk3[11]*(3.*gam3a + gam3b*(1./pifs[4] + 1./pifs[10] + 1./term1))*G[0]*G[20]*G[0]);


// F3/G3(k6,k8,k1) = (-k1,-k3,-p)
F[84] = - 1./(3.*a)*(p.alpha[64]*G[46]*G[17] // a(k1,k3+p) F2(k3,p) G1(k1)
                  + p.alpha[65]*G[26]*G[21] // a(k3,k1+p) F2(k1,p) G1(k3)
                  + p.alpha[51]*G[58]*G[1] // a(p,-k2) F2(k1,k3) G1(p)

         + p.alpha[66]*G[47]*G[16]
         + p.alpha[67]*G[27]*G[20]
         + p.alpha[50]*G[59]*G[0]

         +3.*G[85]) ;

F[85] =1./(3.*a)*(-3.*(2.-hub)*G[85]-3.*hub*G[84]*muak[3]

        -2.*p.beta[32]*G[47]*G[17]

        -2.*p.beta[33]*G[27]*G[21]

        -2.*p.beta[26]*G[59]*G[1]

			 -2*gam2*gamk2[29]*G[16]*G[46]
       -2*gam2*gamk2[27]*G[20]*G[26]
			 -2*gam2*gamk2[6]*G[0]*G[58]

			 -gamk3[12]*(3.*gam3a + gam3b*(1./pifs[6] + 1./pifs[9] + 1./pifs[10]))*G[0]*G[20]*G[16]);


// F4/G4(-k6,-k8,k1,-k1) = (-k1,-k3,p,-p)
F[86] = -1./(24.*a)*(6.*p.alpha[23]*G[82]*G[17] // a(k1,k3) F3(-k3,p,-p) G1(k1)
                    +6.*p.alpha[25]*G[74]*G[21] // a(k3,k1) F3(-k1,p,-p) G1(k3)
                    +6.*p.alpha[7]*G[84]*G[1] // a(p,k2-p) F3(-p,-k1,-k3) G1(p)
                    +6.*p.alpha[70]*G[80]*G[1] // a(-p,k2+p) F3(p,-k1,-k3) G1(p)

                   +6.*p.alpha[25]*G[83]*G[16]
                   +6.*p.alpha[23]*G[75]*G[20]
                   +6.*p.alpha[6]*G[85]*G[0]
                   +6.*p.alpha[71]*G[81]*G[0]



                   + 4.*p.alpha[60]*G[58]*G[73] // a(0,k3) F2(k1,k3) G2(p,-p)
                   + 4.*p.alpha[61]*G[72]*G[59] // a(k3,0) F2(p,-p) G2(k1,k3)

                   + 4.*p.alpha[24]*G[62]*G[27] // a(k1+p,k3-p) F2(-k3,p) G2(-p,-k1)
                   + 4.*p.alpha[22]*G[26]*G[63] // a(k3-p,k1+p) F2(-k1,-p) G2(-k3,p)

                   + 4.*p.alpha[69]*G[46]*G[67] // a(k1-p,k3+p) F2(-k3,-p)G2(-k1,p)
                   + 4.*p.alpha[68]*G[66]*G[47] // a(k3+p,k1-p) F2(-k1,p) G2(-k3,-p)

                   +24.*G[87]);


F[87] = 1./a*(-(2.-hub)*G[87]-hub*G[86]*muak[6]

            -1./24.*(
                  12.*p.beta[12]*G[83]*G[17] // b(k1,k3) G3(-k3,p,-p) G1(k1)
                 +12.*p.beta[12]*G[75]*G[21] // b(k3,k1) G3(-k1,p,-p) G1(k3)
                 +12.*p.beta[3]*G[85]*G[1] // b(p,k2-p) G3(-p,-k1,-k3) G1(p)
                 +12.*p.beta[34]*G[81]*G[1]  // b(-p,k2+p) G3(p,-k1,-k3) G1(p)

                 + 8.*p.beta[31]*G[59]*G[73] // b(k2,0) G2(k1,k3) G2(p,-p)
                 + 8.*p.beta[11]*G[63]*G[27]  // b(k3-p,k1+p) G2(-k3,p) G2(k1,p)
                 + 8.*p.beta[35]*G[67]*G[47] // b(k3+p,k1-p) G2(k3,p) G2(-k1,p)

                 + 8.*gam2*gamk2[34]*G[58]*G[72]
                 + 8.*gam2*gamk2[13]*G[62]*G[26]
                 + 8.*gam2*gamk2[25]*G[66]*G[46]

                 + 12.*gam2*gamk2[14]*G[82]*G[16]
                 + 12.*gam2*gamk2[14]*G[74]*G[20]
                 + 12.*gam2*gamk2[3]*G[84]*G[0]
                 + 12.*gam2*gamk2[2]*G[80]*G[0]


	 							+4.*(gamk3[1]*(3.*gam3a + gam3b*(1./pifs[9] + 1./pifs[7] + 1./pifs[2]))*G[16]*G[0]*G[62]
	 									   + gamk3[5]*(3.*gam3a + gam3b*(1./pifs[10] + 1./pifs[5] + 1./pifs[2]))*G[20]*G[0]*G[26]
	 										 + gamk3[17]*(3.*gam3a + gam3b*(1./pifs[4] + 1./pifs[5] + 1./pifs[3]))*G[20]*G[0]*G[66]
	 										 + gamk3[18]*(3.*gam3a + gam3b*(1./pifs[1] + 1./pifs[7] + 1./pifs[3]))*G[16]*G[0]*G[46]
											 + gamk3[23]*(3.*gam3a + gam3b*(1./term1 + 1./pifs[2] + 1./pifs[3]))*G[0]*G[0]*G[58]
									     + gamk3[24]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[7] + 1./pifs[6]))*G[16]*G[20]*G[72])


							 + gamk4[1]*(24.*gam4a+ 8.*gam4b*(1./(pifs[6]*term1) + 1./(pifs[1]*pifs[10]) + 1./(pifs[9]*pifs[4]))
 																		 + 4.*gam4c*(1./pifs[6] + 1./term1 + 1./pifs[1] + 1./pifs[10] + 1./pifs[9] + 1./pifs[4])
																		 + 4./3.*gam4d*(1./(pifs[3]*pifs[6]) + 1./(pifs[3]*pifs[10]) + 1./(pifs[3]*pifs[9])
 																			 			   + 1./(pifs[2]*pifs[6]) + 1./(pifs[2]*pifs[1]) + 1./(pifs[2]*pifs[4])
 																							 + 1./(pifs[5]*term1) + 1./(pifs[5]*pifs[1]) +	1./(pifs[5]*pifs[9])
 																							 + 1./(pifs[7]*term1) + 1./(pifs[7]*pifs[10]) + 1./(pifs[7]*pifs[4]))
 																		 + 12./3.*gam4c*(1./pifs[2] + 1./pifs[3] + 1./pifs[5]  + 1./pifs[7]))*G[16]*G[20]*G[0]*G[0]));





// F3/G3(k6,k8,k1) = (-k2,-k3,-p)
F[88] = - 1./(3.*a)*(p.alpha[72]*G[46]*G[19] // a(k2,k3+p) F2(k3,p) G1(k2)
                 + p.alpha[74]*G[36]*G[21] // a(k3,k2+p) F2(k2,p) G1(k3)
                 + p.alpha[53]*G[60]*G[1] // a(p,-k1) F2(k2,k3) G1(p)


        + p.alpha[73]*G[47]*G[18]
        + p.alpha[75]*G[37]*G[20]
        + p.alpha[52]*G[61]*G[0]

        +3.*G[89]) ;


F[89] =1./(3.*a)*(-3.*(2.-hub)*G[89]-3.*hub*G[88]*muak[1]

       -2.*p.beta[36]*G[47]*G[19]

       -2.*p.beta[37]*G[37]*G[21]

       -2.*p.beta[27]*G[61]*G[1]


      -2.*gam2*gamk2[21]*G[18]*G[46]
      -2.*gam2*gamk2[28]*G[20]*G[36]
      -2.*gam2*gamk2[5]*G[0]*G[60]

			-gamk3[13]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[2] + 1./pifs[10]))*G[0]*G[20]*G[18]);


// F3/G3(k6,k8,k1) = (-k2,-k3,p)

 F[90] = - 1./(3.*a)*(p.alpha[29]*G[62]*G[19] // a(k2,k3-p) F2(k3,-p) G1(k2)
                  + p.alpha[40]*G[64]*G[21] // a(k3,k2-p) F2(k2,-p) G1(k3)
                  + p.alpha[11]*G[60]*G[1] // a(p,k1) F2(k2,k3) G1(p)

         + p.alpha[28]*G[63]*G[18]
         + p.alpha[41]*G[65]*G[20]
         + p.alpha[10]*G[61]*G[0]

         +3.*G[91]) ;

 F[91] =1./(3.*a)*(-3.*(2.-hub)*G[91]-3.*hub*G[90]*muak[9]

        -2.*p.beta[14]*G[63]*G[19]

        -2.*p.beta[20]*G[65]*G[21]

        -2.*p.beta[5]*G[61]*G[1]

        -2.*gam2*gamk2[16]*G[18]*G[62]
        -2.*gam2*gamk2[22]*G[20]*G[64]
        -2.*gam2*gamk2[5]*G[0]*G[60]

				-gamk3[14]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[3] + 1./pifs[4]))*G[0]*G[20]*G[18]);


// F4/G4(-k6,-k8,k1,-k1) = (-k2,-k3,p,-p)

F[92] = -1./(24.*a)*(6.*p.alpha[32]*G[82]*G[19] // a(k2,k3) F3(-k3,p,-p) G1(k2)
                   +6.*p.alpha[35]*G[76]*G[21] // a(k3,k2) F3(-k2,p,-p) G1(k3)
                   +6.*p.alpha[5]*G[88]*G[1] // a(p,k1-p) F3(-p,-k2,-k3) G1(p)
                   +6.*p.alpha[31]*G[90]*G[1] // a(-p,k1+p) F3(p,-k2,-k3) G1(p)

                  +6.*p.alpha[35]*G[83]*G[18]
                  +6.*p.alpha[32]*G[77]*G[20]
                  +6.*p.alpha[4]*G[89]*G[0]
                  +6.*p.alpha[34]*G[91]*G[0]


                  + 4.*p.alpha[58]*G[60]*G[73] // a(0,k1) F2(k2,k3) G2(p,-p)
                  + 4.*p.alpha[59]*G[72]*G[61] // a(k1,0) F2(p,-p) G2(k2,k3)
                  + 4.*p.alpha[33]*G[62]*G[37] // a(k2+p,k3-p) F2(-k3,p) G2(p,k2)
                  + 4.*p.alpha[30]*G[36]*G[63] // a(k3-p,k2+p) F2(k2,p) G2(-k3,p)

                  + 4.*p.alpha[42]*G[46]*G[65] // a(k2-p,k3+p) F2(-k3,-p)G2(-k2,p)
                  + 4.*p.alpha[43]*G[64]*G[47] // a(k3+p,k2-p) F2(-k2,p) G2(-k3,-p)


                  +24.*G[93]);

F[93] = 1./a*(-(2.-hub)*G[93]-hub*G[92]*muak[5]

        -1./24.*(
                 12.*p.beta[17]*G[83]*G[19] // b(k2,k3) G3(-k3,p,-p) G1(k2)
                +12.*p.beta[17]*G[77]*G[21] // b(k3,k2) G3(-k2,p,-p) G1(k3)
                +12.*p.beta[2]*G[89]*G[1]// b(p,k1-p) G3(-p,-k2,-k3) G1(p)
                +12.*p.beta[16]*G[91]*G[1]  // b(-p,k1+p) G3(p,-k2,-k3) G1(p)

                + 8.*p.beta[30]*G[61]*G[73] // b(k1,0) G2(k2,k3) G2(p,-p)
                + 8.*p.beta[15]*G[63]*G[37]  // b(k3-p,k2+p) G2(-k3,p) G2(k2,p)
                + 8.*p.beta[21]*G[65]*G[47] // b(k3+p,k2-p) G2(k3,p) G2(-k2,p)

                + 8.*gam2*gamk2[33]*G[60]*G[72]
                + 8.*gam2*gamk2[17]*G[62]*G[36]
                + 8.*gam2*gamk2[23]*G[64]*G[46]

                + 12.*gam2*gamk2[18]*G[82]*G[18]
                + 12.*gam2*gamk2[18]*G[76]*G[20]
                + 12.*gam2*gamk2[1]*G[88]*G[0]
                + 12.*gam2*gamk2[19]*G[90]*G[0]


						+4.*(gamk3[4]*(3.*gam3a + gam3b*(1./pifs[10] + 1./pifs[6] + 1./pifs[2]))*G[18]*G[0]*G[62]
							 + gamk3[2]*(3.*gam3a + gam3b*(1./pifs[2] + 1./pifs[7] + 1./pifs[2]))*G[20]*G[0]*G[36]
							 + gamk3[19]*(3.*gam3a + gam3b*(1./pifs[3] + 1./pifs[7] + 1./pifs[1]))*G[20]*G[0]*G[64]
							 + gamk3[20]*(3.*gam3a + gam3b*(1./pifs[4] + 1./pifs[6] + 1./pifs[1]))*G[18]*G[0]*G[46]
							 + gamk3[25]*(3.*gam3a + gam3b*(1./term1 + 1./pifs[9] + 1./pifs[1]))*G[0]*G[0]*G[60]
							 + gamk3[26]*(3.*gam3a + gam3b*(1./pifs[5] + 1./pifs[7] + 1./pifs[6]))*G[18]*G[20]*G[72])


			 + gamk4[2]*(24.*gam4a+ 8.*gam4b*(1./(pifs[5]*term1) + 1./(pifs[3]*pifs[10]) + 1./(pifs[2]*pifs[4]))
														 + 4.*gam4c*(1./pifs[5] + 1./term1 + 1./pifs[3] + 1./pifs[10] + 1./pifs[2] + 1./pifs[4])
														 + 4./3.*gam4d*(1./(pifs[1]*pifs[5]) + 1./(pifs[1]*pifs[10]) + 1./(pifs[1]*pifs[2])
																			 + 1./(pifs[9]*pifs[5]) + 1./(pifs[9]*pifs[3]) + 1./(pifs[9]*pifs[4])
																			 + 1./(pifs[6]*term1) + 1./(pifs[6]*pifs[3]) +	1./(pifs[6]*pifs[2])
																			 + 1./(pifs[7]*term1) + 1./(pifs[7]*pifs[10]) + 1./(pifs[7]*pifs[4]))
														 + 12./3.*gam4c*(1./pifs[9] + 1./pifs[1] + 1./pifs[6]  + 1./pifs[7]))*G[18]*G[20]*G[0]*G[0]));


  	return GSL_SUCCESS;
  }



	void BSPTN::initnb1_fr(double pars[], double extpars[], double epars[], double k[], double x[], double kargs[])
	{
					double a = 0.0001;
					double A = pars[0];
					double omega0 = pars[1];

	       double G[94] = { a,-a,a,-a,a,-a, 0.,0.,0.,0.,0.,0.,a,-a,a,-a,a,-a,a,-a,a,-a,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,
	                            0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

				/*Parameters passed to system of equations */
	          struct param_type5 mypars;
						// set all parameters

							for(int i=0; i<8;i++){
								mypars.kv[i]=k[i];
								mypars.xv[i]=x[i];
							}
							for(int i=0; i<25;i++){
								mypars.args[i]=kargs[i];
							}

						// precompute alpha and beta kernels
						mypars.alpha[0] = alphai(k[1],k[2],x[0]);
						mypars.alpha[1] = alphai(k[2],k[1],x[0]);
						mypars.alpha[2] = alphai(k[0],k[2],x[2]);
						mypars.alpha[3] = alphai(k[2],k[0],x[2]);
						mypars.alpha[4] = alphai(k[1],k[0],x[1]);
						mypars.alpha[5] = alphai(k[0],k[1],x[1]);
						mypars.alpha[6] = alphai(k[3],k[0],x[3]);
						mypars.alpha[7] = alphai(k[0],k[3],x[3]);
						mypars.alpha[8] = alphai(k[4],k[0],x[4]);
						mypars.alpha[9] = alphai(k[0],k[4],x[4]);
						mypars.alpha[10] = alphai(k[5],k[0],x[5]);
						mypars.alpha[11] = alphai(k[0],k[5],x[5]);
						mypars.alpha[12] = alphai(k[5],k[3],kargs[3]);
						mypars.alpha[13] = alphai(k[3],k[5],kargs[3]);
						mypars.alpha[14] = alphai(k[3],kargs[1],kargs[4]);
						mypars.alpha[15] = alphai(k[0],kargs[2],kargs[6]);
						mypars.alpha[16] = alphai(k[5],k[6],kargs[5]);
						mypars.alpha[17] = alphai(kargs[1],k[3],kargs[4]);
						mypars.alpha[18] = alphai(kargs[2],k[0],kargs[6]);
						mypars.alpha[19] = alphai(k[6],k[5],kargs[5]);
						mypars.alpha[20] = alphai(k[4],k[5],kargs[7]);
						mypars.alpha[21] = alphai(k[5],k[4],kargs[7]);
						mypars.alpha[22] = alphai(k[4],kargs[1],kargs[8]);
						mypars.alpha[23] = alphai(k[5],k[7],kargs[9]);
						mypars.alpha[24] = alphai(kargs[1],k[4],kargs[8]);
						mypars.alpha[25] = alphai(k[7],k[5],kargs[9]);
						mypars.alpha[26] = alphai(k[6],k[0],x[6]);
						mypars.alpha[27] = alphai(k[0],k[6],x[6]);
						mypars.alpha[28] = alphai(k[4],k[6],kargs[11]);
						mypars.alpha[29] = alphai(k[6],k[4],kargs[11]);
						mypars.alpha[30] = alphai(k[4],k[2],kargs[12]);
						mypars.alpha[31] = alphai(k[0],kargs[1],kargs[14]);
						mypars.alpha[32] = alphai(k[6],k[7],kargs[13]);
						mypars.alpha[33] = alphai(k[2],k[4],kargs[12]);
						mypars.alpha[34] = alphai(kargs[1],k[0],kargs[14]);
						mypars.alpha[35] = alphai(k[7],k[6],kargs[13]);
						mypars.alpha[36] = alphai(k[6],k[1],kargs[15]);
						mypars.alpha[37] = alphai(k[1],k[6],kargs[15]);
						mypars.alpha[38] = alphai(k[7],k[0],x[7]);
						mypars.alpha[39] = alphai(k[0],k[7],x[7]);
						mypars.alpha[40] = alphai(k[7],k[3],kargs[17]);
						mypars.alpha[41] = alphai(k[3],k[7],kargs[17]);
						mypars.alpha[42] = alphai(k[3],kargs[2],kargs[18]);
						mypars.alpha[43] = alphai(kargs[2],k[3],kargs[18]);
						mypars.alpha[44] = alphai(k[7],k[1],kargs[19]);
						mypars.alpha[45] = alphai(k[1],k[7],kargs[19]);
						mypars.alpha[46] = alphai(k[1],kargs[2],kargs[20]);
						mypars.alpha[47] = alphai(kargs[2],k[1],kargs[20]);
						mypars.alpha[48] = alphai(k[7],k[0],-x[7]);
						mypars.alpha[49] = alphai(k[0],k[7],-x[7]);
						mypars.alpha[50] = alphai(k[6],k[0],-x[6]);
						mypars.alpha[51] = alphai(k[0],k[6],-x[6]);
						mypars.alpha[52] = alphai(k[5],k[0],-x[5]);
						mypars.alpha[53] = alphai(k[0],k[5],-x[5]);
						mypars.alpha[54] = alphai(k[6],kargs[1],kargs[10]);
						mypars.alpha[55] = alphai(k[5],k[2],kargs[21]);
						mypars.alpha[56] = alphai(kargs[1],k[6],kargs[10]);
						mypars.alpha[57] = alphai(k[2],k[5],kargs[21]);
						mypars.alpha[58] = 1.;//alphai(k[5],1.+XMIN,0.);
						mypars.alpha[59] = 1.;//alphai(1.+XMIN,k[5],0.);
						mypars.alpha[60] = alphai(k[6],k[0]*sqrt(2.*(1.+XMIN)),0.);
						mypars.alpha[61] = alphai(k[0]*sqrt(2.*(1.+XMIN)),k[6],0.);
						mypars.alpha[62] = 1.;//alphai(1.-XMAX,k[7],0.);
						mypars.alpha[63] = 1.; //alphai(k[7],1.-XMAX,0.);
						mypars.alpha[64] = alphai(k[5],kargs[2], kargs[24]); // a(k[0],k[2]+p)
						mypars.alpha[65] = alphai(k[7],kargs[1],kargs[22]); // a(k[2],k1+p)
						mypars.alpha[66] = alphai(kargs[2],k[5],kargs[24]); // a(k[2]+p,k1)
						mypars.alpha[67] = alphai(kargs[1],k[7],kargs[22]); // a(k1+p,k[2])
						mypars.alpha[68] = alphai(kargs[2],k[1],kargs[20]);
						mypars.alpha[69] = alphai(k[1],kargs[2],kargs[20]);
						mypars.alpha[70] = alphai(k[0],k[2],x[2]);
						mypars.alpha[71] = alphai(k[2],k[0],x[2]);
						mypars.alpha[72] = alphai(k[6],kargs[2],-kargs[16]); // a(k[1],k[2]+p)
						mypars.alpha[73] = alphai(kargs[2],k[6],-kargs[16]);
						mypars.alpha[74] = alphai(k[7],k[2],kargs[23]); // a(k[2],k[1]+p)
						mypars.alpha[75] = alphai(k[2],k[7],kargs[23]);

	          mypars.beta[0] = betai(k[1],k[2],x[0]);
						mypars.beta[1] = betai(k[0],k[2],x[2]);
						mypars.beta[2] = betai(k[0],k[1],x[1]);
						mypars.beta[3] = betai(k[0],k[3],x[3]);
						mypars.beta[4] = betai(k[0],k[4],x[4]);
						mypars.beta[5] = betai(k[0],k[5],x[5]);
						mypars.beta[6] = betai(k[3],k[5],kargs[3]);
						mypars.beta[7] = betai(k[3],kargs[1],kargs[4]);
						mypars.beta[8] = betai(k[0],kargs[2],kargs[6]);
						mypars.beta[9] = betai(k[6],k[5],kargs[5]);
						mypars.beta[10] = betai(k[4],k[5],kargs[7]);
						mypars.beta[11] = betai(k[4],kargs[1],kargs[8]);
						mypars.beta[12] = betai(k[7],k[5],kargs[9]);
						mypars.beta[13] = betai(k[0],k[6],x[6]);
						mypars.beta[14] = betai(k[4],k[6],kargs[11]);
						mypars.beta[15] = betai(k[4],k[2],kargs[12]);
						mypars.beta[16] = betai(k[0],kargs[1],kargs[14]);
						mypars.beta[17] = betai(k[7],k[6],kargs[13]);
						mypars.beta[18] = betai(k[1],k[6],kargs[15]);
						mypars.beta[19] = betai(k[0],k[7],x[7]);
						mypars.beta[20] = betai(k[3],k[7],kargs[17]);
						mypars.beta[21] = betai(k[3],kargs[2],kargs[18]);
						mypars.beta[22] = betai(k[1],k[7],kargs[19]);
						mypars.beta[23] = betai(k[1],kargs[2],kargs[20]);
						mypars.beta[24] = 0.;// betai(1.-XMAX,k[7],0.);
						mypars.beta[25] = betai(k[0],k[7],-x[7]);
						mypars.beta[26] = betai(k[0],k[6],-x[6]);
						mypars.beta[27] = betai(k[0],k[5],-x[5]);
						mypars.beta[28] = betai(k[6],kargs[1],kargs[10]);
						mypars.beta[29] = betai(k[5],k[2],kargs[21]);
						mypars.beta[31] = betai(k[0]*sqrt(2.*(1.+XMIN)),k[6],0.);
						mypars.beta[32] = betai(kargs[2],k[5],kargs[24]); // a(k[2]+p,k1)
						mypars.beta[33] = betai(k[7],kargs[1],kargs[22]); // a(k[2]+p,k[2])
						mypars.beta[34] = betai(k[2],k[0],x[2]);
						mypars.beta[35] = betai(kargs[2],k[1],kargs[20]);
						mypars.beta[36] = betai(k[6],kargs[2],-kargs[16]);
						mypars.beta[37] = betai(k[7],k[2],kargs[23]);

						// precompute squares of k vectors used in gammas
						//  p.gamk[0]  // |kv[0]|^2
						//  p.gamk[1]  // |kv[1]|^2
						//  p.gamk[2]  // |kv[2]|^2
						//  p.gamk[3]  // |kv[3]|^2
						//  p.gamk[4]  // |kv[4]|^2
						//  p.gamk[5] // |kv[5]|^2
						//  p.gamk[6] // |kv[6]|^2
						//  p.gamk[7] // |kv[7]|^2
						//  p.gamk[8]  // |k2+2p|^2
						//  p.gamk[9] // |p+k1|^2
						//  p.gamk[10] // |k3+p|^2

						for(int i=0;i<11;i++){
							mypars.gamk[i]=kargs[25+i];
						}

						// cosmo and beyond lcdm params
						mypars.omega0 = omega0;

						for (int i=0; i<maxpars;i++){
							mypars.extpars[i] = extpars[i];
						}

	          gsl_odeiv2_system sys = {funcfr, jacb, 94, &mypars};

				//  this gives the system, the method, the initial interval step, the absolute error and the relative error.
					gsl_odeiv2_driver * d =
					gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
														epars[1], epars[2] ,epars[3]);

	// smallest possible accuracy and initial step when comparing to analytic result
	// must reduce initial step size for earlier times!
	// rk8pd is optimal
					int status1 = gsl_odeiv2_driver_apply (d, &a, A, G);


				/*Allocation of array values */


				// For tree level initialisation

	  		/*1st order */

				//F1(k;a), G1(k;a)
				F1b_k1 = G[0] ;
				G1b_k1 = G[1] ;

				// F1(k3;a), G1(k3;a)
				F1b_k3 = G[2];
				G1b_k3 = G[3];

				//F1(k2;a), G1(k2;a)
				F1b_k2 = G[4];
				G1b_k2 = G[5];


				// for B321 and B222 initialisaiton

	      /* 1st order */

	      //F1(k1;a), G1(k1;a)
				F1b_1 = G[16] ;
				G1b_1 = G[17] ;

				// F1(k2;a), G1(k2;a)
				F1b_2 = G[18];
				G1b_2 = G[19];

				//F1(k3;a), G1(k3;a)
				F1b_3 = G[20];
				G1b_3 = G[21];


				/*2nd order*/

				//F2/G2(k3,k2) for tree,  F2/G2(k1-p,k2+p) for B222
				F2b_k23 =  G[6];
				G2b_k23 =  G[7];

	      //F2/G2(k1,k3) for tree, F2/G2(-p,k2+p) for B222
	      F2b_k13 =  G[8];
	      G2b_k13 =  G[9];

	      //F2/G2(k1,k2) for tree, F2/G2(p,k1-p)
	      F2b_k12 =  G[10];
	      G2b_k12 =  G[11];


	      // F2/G2(p,k2-p)
	      F2b_p2mp = G[22];
	      G2b_p2mp = G[23];

	      // F2/G2(p,k3-p)
	      F2b_p3mp = G[24];
	      G2b_p3mp = G[25];

	      // F3/G3(k1,p,k2-p)
	      F3b_12mp = G[30];
	      G3b_12mp = G[31];

	      // F3/G3(k1,p,k3-p)
	      F3b_13mp = G[34];
	      G3b_13mp = G[35];

	      // F3/G3(k2,p,k3-p)
	      F3b_23mp = G[40];
	      G3b_23mp = G[41];

	      // F3/G3(k2,p,k1-p)
	      F3b_21mp = G[44];
	      G3b_21mp = G[45];


	      // F3/G3(k3,p,k2-p)
	      F3b_32mp = G[50];
	      G3b_32mp = G[51];

	      // F3/G3(k3,p,k1-p)
	      F3b_31mp = G[54];
	      G3b_31mp = G[55];

	      // F3/G3(k1,k2) for B321-II
	      F2b_12a = G[56];
	      G2b_12a = G[57];

	      // F3/G3(k1,k3) for B321-II
	      F2b_13a = G[58];
	      G2b_13a = G[59];

	      // F3/G3(k2,k3) for B321-II
	      F2b_23a = G[60];
	      G2b_23a = G[61];

	      //F3/G3(k1,p,-p)
	      F3b_1pp = G[74];
	      G3b_1pp = G[75];

	      //F3/G3(k2,p,-p)
	      F3b_2pp = G[76];
	      G3b_2pp = G[77];

	      //F3/G3(k3,p,-p)
	      F3b_3pp = G[82];
	      G3b_3pp = G[83];

	      // F4/G4(-k1,-k2,p,-p)

	      F4b_12pp = G[78];
	      G4b_12pp = G[79];

	      // F4/G4(-k1,-k3,p,-p)

	      F4b_13pp = G[86];
	      G4b_13pp = G[87];

	      // F4/G4(-k2,-k3,p,-p)

	      F4b_23pp = G[92];
	      G4b_23pp = G[93];

				gsl_odeiv2_driver_free(d);
	}
