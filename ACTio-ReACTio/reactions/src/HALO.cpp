#if HAVE_CONFIG_H
# include <config.h>
#endif

#include <cfloat>
#include <cmath>
#include <iostream>
#include <stdlib.h>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <math.h>       /* pow */


#include "Common.h"
#include "Quadrature.h"
#include "HALO.h"
#include "SCOL.h"
#include "Spline.h"
#include "SpecialFunctions.h"
#include "PowerSpectrum.h"
#include "SPT.h"
#include "BSPT.h"
#include "BSPTN.h"


#include <stdio.h>
#include <gsl/gsl_sf.h>       /* pow */
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_math.h>

#include <random>
#include <functional>
#include<omp.h>

using std::cref;
using std::bind;

extern const int ERROR_MESSAGE_LEN_HALO = 512;
char error_message_halo[ERROR_MESSAGE_LEN_HALO];

void react_error_halo(const char* msg) {
    const int truncate = ERROR_MESSAGE_LEN_HALO-10;

    #pragma omp critical
    {
        snprintf(error_message_halo, truncate, "HALO.cpp error: %s", msg);
        printf("%s\n", error_message_halo);
    }
}


/* Code to calculate halo model and reaction terms */

/* Based on work in halo model reaction papers:  */
/* 1812.05594, 2005.12184, 2105.12114   */


/* Makes use of spherical collapse ODE systems SCOL.cpp*/
/* See header file for argument values and descriptions */

/*Load the linear power spectra */
// For massless neutrinos or no neutrinos, user should set P_cb = P_nu = P_l (total matter spectrum)
HALO::HALO(const Cosmology& C_, const PowerSpectrum& P_l_, const PowerSpectrum& P_cb_, const PowerSpectrum& P_nu_, const PowerSpectrum& P_cbl_, real epsrel_)
: C(C_), P_l(P_l_), P_cb(P_cb_), P_nu(P_nu_), P_cbl(P_cbl_)
{
    epsrel = epsrel_;
}

// Window function (top hat)
static real W(real x) {
    if(x < 1e-5)
        return 1. - (1./30.)*x*x;
    else
        return 3./pow3(x) * (sin(x) - x*cos(x));
}

// Window function derivative wrt R (top hat)
static real Wd(double k, double R) {
  double kR = k*R;
  double skr = sin(kR);
  double ckr = cos(kR);
        return (3.*skr - 9.*(skr-kR*ckr)/pow2(kR))/kR/R ;
}


///// Modified Spherical Collapse  ////////

/* Initialise splines for delta_col , a_vir, delta_avir  and nu = delta_col/sigma as functions of mass index (10^mass index) */
Spline delta_col;
Spline a_vir;
Spline delta_avir;
Spline mysig;
Spline mysigp;

Spline linear_growth;

// sigma^2 (EQ 42 of 1812.05594) in terms of log_10(mass)
// --- LCDM @ z=0 value --
// include scale factor later
static real sigma_integrand(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * pow2(W(k*R)*Dl_spt/dnorm_spt);
}


// d(sigma_8^2)/dR
static real sigma8d_integrand(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * pow2(Dl_spt/dnorm_spt) * 2.*W(k*R)*Wd(k,R) ;
}


// Used when input power spectra is given at z and is already modified
static real sigma_integrand_modcamb(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * pow2(W(k*R));

}

// d(sigma_8^2)/dR
static real sigma8d_integrand_modcamb(const PowerSpectrum& P, real R, real k) {
  return k*k/(2.*M_PI*M_PI) * P(k) * 2.*W(k*R)*Wd(k,R) ;

}




// Initialise spherical collapse quantities for real cosmology
int HALO::scol_init(double pars[], double extpars[], bool modcamb, int model) const{
  SCOL scol;
  IOW iow;
  // number of points in mass loop (default is 30)
  int loop_N = (int)pars[5];

  int status = 0;

// arrays to store values
double lgmass[loop_N],sigar[loop_N],scolar0[loop_N],scolar1[loop_N],scolar2[loop_N];

// initialise linear growth spline for the theory of gravity/dark energy.
int loop_Nk = 300;
double kmin = 1e-5;
double kmax = 100.;
double kval_tab[loop_Nk],ling_tab[loop_Nk];

double omegacb = pars[1]-pars[2];

// initialise linear growth factor if we are using LCDM P(k,z=0) input so we can rescale
if(!modcamb){
for(int i = 0; i< loop_Nk; i++){
  kval_tab[i] =  kmin * exp(i*log(kmax/kmin)/(loop_Nk-1.));
  // numerical initialisation of linear growth in general theory of gravity/dark energy - see SpecialFunctions.cpp
  iow.initn_lin(pars, extpars, kval_tab[i], model);
  ling_tab[i] = F1_nk/dnorm_spt; // = D(z)/D_LCDM(z=0)
    }
}
// otherwise linear_growth should be set to unity (so functions like plinear_cosmosis still work fine)
else{
  for(int i = 0; i< loop_Nk; i++){
    kval_tab[i] =  kmin * exp(i*log(kmax/kmin)/(loop_Nk-1.));
    ling_tab[i] = 1.;
      }
}

    // spline linear growth factor
     linear_growth = LinearSpline(loop_Nk,kval_tab,ling_tab);


 double sig1,sig2;

  if (!modcamb) {
  // calculate sigma_8^2 and it's smoothing scale derivative @ R=8.
   sig1 = Integrate<ExpSub>(bind(sigma_integrand, cref(P_l), 8., std::placeholders::_1), 1e-4, 50., 1e-5,1e-5);
   sig2 = Integrate<ExpSub>(bind(sigma8d_integrand, cref(P_l), 8., std::placeholders::_1), 1e-4, 50., 1e-5,1e-5);
  }
  else{
     sig1 = Integrate<ExpSub>(bind(sigma_integrand_modcamb, cref(P_cbl), 8., std::placeholders::_1), 1e-4, 50., 1e-5,1e-5);
     sig2 = Integrate<ExpSub>(bind(sigma8d_integrand_modcamb, cref(P_cbl), 8., std::placeholders::_1), 1e-4, 50., 1e-5,1e-5);
  }

 if (!gsl_finite(sig1)) {
  react_error_halo("sigma_8 evaluated to non-numerical value");
}

if (!gsl_finite(sig2)) {
  react_error_halo("sigma_8 derivative evaluated to non-numerical value");
}


//ZW: comment out to improve speed
//printf("%s %e \n","LCDM sigma8: ", sqrt(sig1) );

//ZW: temporary no mass dependence
// Removed from loop if no mass dependence
   double myscolparams[3];

// initialise delta_c, a_vir, delta_vir
   scol.myscol(myscolparams, pars[0], omegacb, pars[2], 1.0, sig1, sig2, extpars, model);

   if(scol.error.errorno != 0) {
    status = scol.error.errorno;
   }

// log(mass) loop for calculating the collapse quantities
//#pragma omp parallel for
for(int i = 0; i< loop_N; i++){
  //ZW: temporary remove no mass dependence
  //  double myscolparams[3];

// log mass
   lgmass[i] = Mmin + i*(Mmax-Mmin)/(loop_N-1.); // linear spacing

// top hat radius
   double Rth = 0.1*pow((Gnewton*pow(10, lgmass[i]))/(5.*omegacb),ONE/THREE);

// initialise delta_c, a_vir, delta_vir
  //ZW: temporary remove no mass dependence
  //  scol.myscol(myscolparams, pars[0], omegacb, pars[2], Rth, sig1, sig2, extpars, model);
  //  if(scol.error.errorno != 0) {
  //   status = scol.error.errorno;
  //  }

// calculate sigma^2
  if (!modcamb) {
   sigar[i] = sqrt(Integrate<ExpSub>(bind(sigma_integrand, cref(P_l), Rth, std::placeholders::_1), 1e-4, 50., 1e-5,1e-5));
 }
 else{
   if (model!=9) {
   sigar[i] = sqrt(Integrate<ExpSub>(bind(sigma_integrand_modcamb, cref(P_cbl), Rth, std::placeholders::_1), 1e-4, 50., 1e-5,1e-5));
    }
    // if unscreened, we should use the cb spectrum of the modified gravity model
    else{
  sigar[i] = sqrt(Integrate<ExpSub>(bind(sigma_integrand_modcamb, cref(P_cb), Rth, std::placeholders::_1), 1e-4, 50., 1e-5,1e-5));
    }
   }



// store values in arrays
   scolar0[i] = myscolparams[0];
   scolar1[i] = myscolparams[1];
   scolar2[i] = myscolparams[2];

}


// Spline over arrays
// sigma^2
  mysig = CubicSpline(loop_N,lgmass,sigar);
// delta_c
  delta_col = CubicSpline(loop_N,lgmass,scolar0);
// a_virial
  a_vir = CubicSpline(loop_N,lgmass,scolar1);
// delta_virial
  delta_avir = CubicSpline(loop_N,lgmass,scolar2);

  return status;

}



///// Pseudo spherical collapse  ////////

Spline delta_colp;
Spline a_virp;
Spline delta_avirp;

// sigma^2 for pseudo cosmology --- modified
static real sigma_integrandp(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * pow2(linear_growth(k)*W(k*R));
}

// d(sigma_8^2)/dR
static real sigma8d_integrandp(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * pow2(linear_growth(k)) * 2.*W(k*R)*Wd(k,R) ;
}


// Used when input power spectra is given at z and is already modified
// sigma^2 for pseudo cosmology --- modified
static real sigma_integrandp_modcamb(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * pow2(W(k*R));
}
// d(sigma_8^2)/dR
static real sigma8d_integrandp_modcamb(const PowerSpectrum& P, real R, real k) {
    return k*k/(2.*M_PI*M_PI) * P(k) * 2.*W(k*R)*Wd(k,R) ;
}

// Initialise spherical collapse quantities for pseudo cosmology
int HALO::scol_initp(double pars[], bool modcamb) const{
  SCOL scol;

  int loop_N = 30;
  int status = 0;

    // arrays to store values
  double myscolparams[3],lgmass[loop_N],sigar[loop_N],scolar0[loop_N],scolar1[loop_N],scolar2[loop_N];

double sig1,sig2;
  if (!modcamb) {
  // calculate sigma_8^2 and it's smoothing scale derivative @ R=8.
   sig1 = Integrate<ExpSub>(bind(sigma_integrandp, cref(P_l), 8., std::placeholders::_1), 1e-4, 50., 1e-5,1e-5);
   sig2 = Integrate<ExpSub>(bind(sigma8d_integrandp, cref(P_l), 8., std::placeholders::_1), 1e-4, 50., 1e-5,1e-5);

  }
  else{
     sig1 = Integrate<ExpSub>(bind(sigma_integrandp_modcamb, cref(P_l), 8., std::placeholders::_1), 1e-4, 50., 1e-5,1e-5);
     sig2 = Integrate<ExpSub>(bind(sigma8d_integrandp_modcamb, cref(P_l), 8., std::placeholders::_1), 1e-4, 50., 1e-5,1e-5);
  }

  //ZW: comment out to improve speed
  //printf("%s %e \n", "Modified sigma8: ", sqrt(sig1) );

// feed back modified sigma8 for other uses (see pyreact or cosmosis modules)
    pars[6] = sqrt(sig1);

   // theory params for GR - to be safe
   double extpars[maxpars];
   extpars[0] = 1e-15;
   extpars[1] = 1.;
   extpars[2] = 1.;

  // set collapse model to GR
  int mymodel = 1;

// initialise spherical collapse quantities in GR (independent of R (or M))
  scol.myscol(myscolparams, pars[0], pars[1], 0, 1., sig1, sig2, extpars, mymodel);
  if(scol.error.errorno != 0) {
   status = scol.error.errorno;
  }

//#pragma omp parallel for
for(int i = 0; i< loop_N; i++){
   lgmass[i] = Mmin + i*(Mmax-Mmin)/(loop_N-1.); // linear spacing

   double Rth = 0.1*pow((Gnewton*pow(10, lgmass[i]))/(5.*pars[1]),ONE/THREE);

// calculate sigma^2
  if (!modcamb) {
    sigar[i] = sqrt(Integrate<ExpSub>(bind(sigma_integrandp, cref(P_l), Rth, std::placeholders::_1), 1e-4, 50., 1e-5,1e-5));
  }
  else{
    sigar[i] = sqrt(Integrate<ExpSub>(bind(sigma_integrandp_modcamb, cref(P_l), Rth, std::placeholders::_1), 1e-4, 50., 1e-5,1e-5));
  }
// store values in arrays
    scolar0[i] = myscolparams[0];
    scolar1[i] = myscolparams[1];
    scolar2[i] = myscolparams[2];
      }

// Spline over arrays
  // sigma^2
    mysigp = CubicSpline(loop_N,lgmass,sigar);
  // delta_c
    delta_colp = CubicSpline(loop_N,lgmass,scolar0);
  // a_virial
    a_virp = CubicSpline(loop_N,lgmass,scolar1);
  // delta_virial
    delta_avirp = CubicSpline(loop_N,lgmass,scolar2);

    return status;

}


/* HALO MODEL QUANTITIES */

//  nu and dnu/dlnM , a being the scale factor,  used in mass functions

// real
static double mynu(double lgmass, void * params){
  (void)(params);
  return  delta_col(lgmass)/mysig(lgmass);
}
static double nu_der(double lgmass){
  gsl_function F;
  double result, abserr;
  F.params = 0;
  F.function = &mynu;
  gsl_deriv_central(&F, lgmass, 1e-5, &result, &abserr);
return result;
}

// pseudo
static double mynup(double lgmass, void * params){
  (void)(params);
  return delta_colp(lgmass)/mysigp(lgmass);
}
static double nu_derp(double lgmass){
  gsl_function F;
  double result, abserr;
  F.params = 0;
  F.function = &mynup;
  gsl_deriv_central(&F, lgmass, 1e-5, &result, &abserr);
return result;
}



////////////////////////////////////////////////////////////////////////////////////////

/* Virial radius */
// See equation 36 of 1812.05594
// Note that we use Gnewton in units of [Mpc/M_solar  * (km/s)^2 ] -- see SCOL.h
// Correspondingly, we have rho_m = 3 x 10^4 / (8 x pi x Gnewton) x Omega_m [M_solar/h / (Mpc/h)^3]

// pars:
// 0 = scale factor of collapse , 1= Omega_{m,0}, 2 = Omega_{nu,0}
// real
double HALO::rvirial(double lgmass, double pars[]) const {
      //virial overdensity
      double omegacb = pars[1]-pars[2];
      double Deltavir = (1.+ delta_avir(lgmass))*pow3(pars[0]/a_vir(lgmass));
      return 0.1*pow((pow(10.,lgmass)*Gnewton)/(5.*omegacb*Deltavir), 1./3.);
    }

//pseudo
double HALO::rvirialp(double lgmass, double pars[]) const {
      //virial overdensity
      double Deltavir = (1.+ delta_avirp(lgmass))*pow3(pars[0]/a_virp(lgmass));
      return 0.1*pow((pow(10.,lgmass)*Gnewton)/(5.*pars[1]*Deltavir), 1./3.);
    }


////////////////////////////////////////////////////////////////////////////////////////


/* Virial concentration */
// see equation 46 of 1812.05594
// real
double HALO::cvirial(double lgmass, double acol) const {

      double myc0 = 9.;
      double alpha = 0.13;

      double lgmstar;

      const double upper = Mmax;
      const double lower = Mmin;


      std::mt19937 gen(2);

      std::uniform_real_distribution<> dis(lower, upper);

      double pos_pt = dis(gen);
      double neg_pt = dis(gen);

      while ( delta_col(lgmass)/mysig(pos_pt)-1.  < 0.){
          pos_pt = dis(gen);}

      int mycount = 1;
      while ( delta_col(lgmass)/mysig(neg_pt)-1. > 0.){
      neg_pt = dis(gen);
      mycount += 1;
       if (mycount>100) {
         fprintf(stderr, "Failed to converge in cvirial.\n");
         return g_de*myc0*pow(10.,-alpha*(lgmass-Mmin))*acol;
                       }
                    }

       const double about_zero_mag = 1e-3;
      for (;;)
      {
          double mid_pt = (pos_pt + neg_pt)/2.0;
          double f_mid_pt = delta_col(lgmass)/mysig(mid_pt)-1.;

          if (fabs(f_mid_pt)  < about_zero_mag){
              lgmstar = mid_pt;
              break;
                }

          if (f_mid_pt >= 0.){
              pos_pt = mid_pt;}
          else{
              neg_pt = mid_pt;}
      }

    return g_de*myc0*pow(10.,-alpha*(lgmass-lgmstar))*acol;
  }




//pseudo
  double HALO::cvirialp(double lgmass, double acol) const {

        double myc0 = 9.;
        double alpha = 0.13;
        double lgmstar;

        const double upper = Mmax;
        const double lower = Mmin;


        std::mt19937 gen(2);

        std::uniform_real_distribution<> dis(lower, upper);

        double pos_pt = dis(gen);
        double neg_pt = dis(gen);

        while ( delta_colp(lgmass)/mysigp(pos_pt)-1.  < 0.){
            pos_pt = dis(gen);}

        int mycount = 1;
        while ( delta_colp(lgmass)/mysigp(neg_pt)-1. > 0.){
        neg_pt = dis(gen);
        mycount += 1;
        if (mycount > 100 ) {
          fprintf(stderr, "Failed to converge in cvirialp.\n");
          return myc0*pow(10.,-alpha*(lgmass-Mmin))*acol;
                        }
                      }

         const double about_zero_mag = 1e-5;
        for (;;)
        {
            double mid_pt = (pos_pt + neg_pt)/2.0;
            double f_mid_pt = delta_colp(lgmass)/mysigp(mid_pt)-1.;

            if (fabs(f_mid_pt)  < about_zero_mag){
                lgmstar = mid_pt;
                break;
                  }

            if (f_mid_pt >= 0.){
                pos_pt = mid_pt;}
            else{
                neg_pt = mid_pt;}
        }

      return myc0*pow(10.,-alpha*(lgmass-lgmstar))*acol;

    }

////////////////////////////////////////////////////////////////////////////////////////

/* Virial number density */
// See equation 40 of 1812.05594
// real
  double HALO::nvirial(double lgmass, double omega0) const {
       void *params = nullptr;

       // Quantities as in https://arxiv.org/abs/1812.05594
       double q = 0.75;
       double p = 0.3;
       double aconst = 1./3.056 /1.012;

       // Quantities as in https://arxiv.org/abs/1909.02561
       // double q = 0.7665;
       // double p = 0.2488;
       // double aconst = 0.3292;


       double v2 = pow2(mynu(lgmass,params));
       double vfv = aconst*sqrt(2./M_PI * q * v2 ) * (1.+ pow(v2*q,-p)) * exp(-q*v2/2.);
       double myrho =  15.* omega0/(4.*M_PI*Gnewton*pow3(0.1));
       double Mvir = pow(10.,lgmass);

       double nuprime = 0.434294*nu_der(lgmass); // log10 (Mass) derivative of  nu, constant is conversion from ln to log10 derivative (=log_10(e))

       return myrho * vfv/Mvir *  1./ mynu(lgmass,params) * nuprime;    /// Halo mass function (ST for now)
  }

// pseudo
  double HALO::nvirialp(double lgmass, double omega0) const {
       void *params = nullptr;

       // Quantities as in https://arxiv.org/abs/1812.05594
       double q = 0.75;
       double p = 0.3;
       double aconst = 1./3.056 /1.012;

       // Quantities as in https://arxiv.org/abs/1909.02561
       // double q = 0.7665;
       // double p = 0.2488;
       // double aconst = 0.3292;


       double v2 = pow2(mynup(lgmass,params));
       double vfv = aconst*sqrt(2./M_PI * q * v2 ) * (1.+ pow(v2*q,-p)) * exp(-q*v2/2.);
       double myrho =  15.* omega0/(4.*M_PI*Gnewton*pow3(0.1));
       double Mvir = pow(10.,lgmass);

       double nuprime = 0.434294*nu_derp(lgmass); // log10 (Mass) derivative of  nu, constant is conversion from ln to log10 derivative (=log_10(e))

       return myrho * vfv/Mvir *  1./ mynup(lgmass,params) * nuprime;    /// Halo mass function (ST for now)

  }


  ////////////////////////////////////////////////////////////////////////////////////////

// STANDARD FT OF NFW PROFILE (SEE EQ. 81 OF 0206508 for example)
double HALO::halo_profileK2(double k, double Rvir, double mycvir) const {
     double rs = Rvir/mycvir;
     double krs = k*rs;
     double csp = 1.+mycvir;
     double term1 = 1./(log(csp) - mycvir/csp);

     return term1  * (sin(krs) * (gsl_sf_Si(csp*krs) - gsl_sf_Si(krs))
                    - sin(krs*mycvir)/(csp*krs)
                    + cos(krs) * (gsl_sf_Ci(csp*krs) - gsl_sf_Ci(krs)));
}



// 1 - halo term integrand
// pars:
// 0 = scale factor of collapse , 1= Omega_{m,0}, 2 = Omega_{nu,0}

static double halo_integrand(const Cosmology& C, const PowerSpectrum& P_l, double pars[], double k, double lgmvir){
        HALO halo(C, P_l,P_l,P_l,P_l, 1e-3); // doesn't matter which Pk we feed here

        // background matter density
        double omegacb = pars[1]-pars[2];
        double myrho = 15.*omegacb/(4.*M_PI*Gnewton*pow3(0.1));

        // Halo density profile components //

        // virial radius
        double Rvir = halo.rvirial(lgmvir, pars); // run first to initialise delta_collapse used in c_vir

        // virial concentration
        double mycvir = halo.cvirial(lgmvir,pars[0]);

        // virial mass
        double Mvir = pow(10,lgmvir);

        // halo density profile
        double hpk = fabs(halo.halo_profileK2(k,Rvir,mycvir));
        double halof2 = pow2(hpk);

        // virial number density
        double mynvir = halo.nvirial(lgmvir, omegacb);

        double prefac = pow2(Mvir/myrho);

        return prefac * mynvir * halof2;

}


// real one halo term

// pars:
// 0 = scale factor of collapse , 1= Omega_{m,0}, 2 = Omega_{nu,0}

double HALO::one_halo(double k, double pars[]) const {
    // constant converts from an integral over ln M to log10 M
    return  Integrate(bind(halo_integrand, cref(C), cref(P_l), pars, k, std::placeholders::_1),  1.01*Mmin , 0.99*Mmax  , epsrel)/0.434294;
  }



// 1 - halo term integrand for pseudo
// pars:
// 0 = scale factor of collapse , 1= Omega_{m,0}, 2 = Omega_{nu,0}
  static double halo_integrandp(const Cosmology& C, const PowerSpectrum& P_l, double pars[], double k, double lgmvir){
          HALO halo(C, P_l,P_l,P_l,P_l, 1e-3); // doesn't matter which Pk we feed here

          // background matter density
          double myrho = 15.*pars[1]/(4.*M_PI*Gnewton*pow3(0.1));

          // Halo density profile components //

          // virial radius
          double Rvir = halo.rvirialp(lgmvir, pars); // run first to initialise delta_collapse used in c_vir

          // virial concentration
          double mycvir = halo.cvirialp(lgmvir,pars[0]);

          // virial mass
          double Mvir = pow(10,lgmvir);

          // halo density profile
          double hpk = halo.halo_profileK2(k, Rvir,mycvir);
          double halof2 = pow2(hpk);


          // virial number density
          double mynvir = halo.nvirialp(lgmvir, pars[1]);

          double prefac = pow2(Mvir/myrho);

          return prefac * mynvir * halof2;

  }

// pseudo one halo term
  double HALO::one_halop(double k, double pars[]) const {
      // constant converts from an integral over ln M to log10 M
      return  Integrate(bind(halo_integrandp, cref(C), cref(P_l), pars, k, std::placeholders::_1), 1.01*Mmin ,0.99*Mmax , epsrel)/0.434294;
    }



/* Reaction with massive neutrinos */

/* Initialise reaction parameters */

// kstar, mathcalE and linear power spectrum spline
double kstar, bigE;

/// reaction parameter intiialisation for a single redshift
// pars:
// 0 = scale factor of collapse , 1= Omega_{m,0}, 2 = Omega_{nu,0}
// extpars - beyond LCDM parameters (see HALO.h for more details on arguments )

void HALO::react_init_nu(double pars[], double extpars[], bool modcamb, bool modg, int model) const{
  SPT spt(C, P_cb, epsrel);
  SPT spt2(C, P_l, epsrel);

  if ( fabs(pars[2])>1e-6 && !modcamb) {
    react_error_halo("Omega neutrino is not 0 but you have specified a LCDM transfer - set modcamb = true or set Omega_nu = 0");
  }

// adding in because of issues of mathcal{E} not going to 1 for some particular transfers without MG
  if (!modg) {
    kstar = 1.;
    bigE = 1.;
    return;
  }

  double pspt, psptp, p1h, p1hp, plm, plcb, plnu, psptcbv, argument;
// double psptpnosc,psptpcbv; // used for unscreened approximation - see below
  // fv
  double fv  = pars[2]/pars[1];
  double fv2 = pow(fv,2);
  double fvt = 1.-fv;
  double fvt2  = pow(fvt,2);

// See HALO.h
// Scale at which to calculate mathcal{E}. see 1812.05594 eq 54. This should just be sufficiently small so that the 1-halo term is constant.
  double ke = EHMR;
// Smallest scale at which to trust 1-loop PT (at all relevant redshifts). see  see 1812.05594 eq 55. The authors find 0.06h/Mpc to be the optimal value.
  double k0 = K0HMR;

/// mathcal E///
// See equation 5 of 2005.12184
  bigE = fvt2*one_halo(ke,pars)/one_halop(ke,pars);

  /// Calculate kstar///
  // See equation 6 of 2005.12184
// 1-halo terms
  p1h = one_halo(k0,pars); // real
  p1hp = one_halop(k0,pars); // pseudo

// spt terms
if(!modcamb){
  // linear terms
    plm = pow2(linear_growth(k0))*P_l(k0);
    plcb = plm;
    plnu = plm;

// Real PT spectrum
  pspt = spt.PLOOPn2(1, k0, pars, extpars, 1e-3, model) + p1h ;
// GR-1-loop spectrum with linear growth replaced by modified gravity growth (unscreened)
  psptp = spt.PLOOPn2(4, k0, pars, extpars, 1e-3, model);
}
else{
  plm = P_l(k0);
  plcb = P_cb(k0);
  plnu = P_nu(k0);

// Real PT spectrum
  pspt = spt.PLOOPn2_nu(1, k0, pars, extpars, 1e-3, model) + p1h ;

  // //  Unscreened approx
  // psptpnosc = spt.PLOOPn2_nu(2,  k0, pars, extpars, 1e-3, model);
  //
  // psptpcbv = sqrt(psptpnosc*plnu);
  // // pseudo SPT term
  // psptp = (fvt2 * psptpnosc + 2.*fv*fvt*psptpcbv + fv2*plnu);

  psptp = spt2.PLOOPn2_nu(2,  k0, pars, extpars, 1e-3, model);

  }
// cbv terms
  psptcbv = sqrt(pspt*plnu);

// reaction  @ k0 with 2halo term = 1 loop spt
  double rsptk0 = (fvt2 * pspt + 2.*fv*fvt*psptcbv + fv2*plnu)/(psptp + p1hp);


  double a1 = ((1.-bigE)*plcb*fvt2);

  double term1 =  fv2*plnu + (plm+p1hp)*rsptk0 -fvt2*(bigE*plcb + p1h);
  double term2 = 2.*sqrt(fv2*(plm+p1hp)*plnu*rsptk0);

  if (bigE>1.) {
    argument =  (term1 + term2)/a1; // 1st root
  }
  else{
    argument = (term1 - term2)/a1; // 2nd root
    }

// For very small modifications, numerics can sometimes generate a negative argument.
// To prevent 'nans' in the log taken below in such cases, we manually set kstar and bigE to their (approximately) GR values.
  if(argument<0.001 || argument>1.){
    kstar = 1.;
    bigE = 1.;
  }
  else{
  kstar = (-k0/log(argument));
  }

 // Make sure the root is actually a root
  double prefac = ((1.-bigE)*exp(-k0/kstar) + bigE);
  double ph = prefac*plcb + p1h; // real
  double php = plm + p1hp; // pseudo
  double hmrk0 = (fvt2*ph +  2.*fv*fvt*sqrt(ph*plnu) + fv2*plnu)/php;

  double dif = fabs(hmrk0 - rsptk0);
  if (dif>1e-4) {
    kstar = 1.;
    bigE = 1.;
        }

}

// reaction parameter initialisation for multiple redshifts --- useful for lensing --  much faster than running react_init multiple times
// Reaction using 1-loop splines (ploopr and ploopp) as a function of redshift -- these should be initialised with the ploop_init or ploop_init_nu functions in SPT.cpp
// pars:
// 0 = scale factor of collapse , 1= Omega_{m,0}, 2 = Omega_{nu,0}

void HALO::react_init_nu_multiz(double pars[], Spline ploopr, Spline ploopp, bool modcamb, bool modg) const{
  SPT spt(C, P_cb, epsrel);
  SPT spt2(C, P_l, epsrel);

  if ( fabs(pars[2])>1e-6 && !modcamb) {
    react_error_halo("Omega neutrino is not 0 but you have specified a LCDM transfer - set modcamb = true or set Omega_nu = 0");
  }

// adding in because of issues of mathcal{E} not going to 1 for some particular transfers without MG
  if (!modg) {
    kstar = 1.;
    bigE = 1.;
    return;
  }

  double pspt, psptp, p1h, p1hp, plm, plcb, plnu, psptcbv, argument;

  // fv
  double fv  = pars[2]/pars[1];
  double fv2 = pow(fv,2);
  double fvt = 1.-fv;
  double fvt2  = pow(fvt,2);

  // See HALO.h
// Scale at which to calculate mathcal{E}. see 1812.05594 eq 54. This should just be sufficiently small so that the 1-halo term is constant.
  double ke = EHMR;
// Smallest scale at which to trust 1-loop PT (at all relevant redshifts). see  see 1812.05594 eq 55. The authors find 0.06h/Mpc to be the optimal value.
  double k0 = K0HMR;

/// mathcal E///
// See equation 5 of 2005.12184
  bigE = fvt2*one_halo(ke,pars)/one_halop(ke,pars);

  /// Calculate kstar///
  // See equation 6 of 2005.12184
// 1-halo terms
  p1h = one_halo(k0,pars); // real
  p1hp = one_halop(k0,pars); // pseudo

// spt terms
if(!modcamb){
  // linear terms
    plm = pow2(linear_growth(k0))*P_l(k0);
    plcb = plm;
    plnu = plm;
}
else{
  plm = P_l(k0);
  plcb = P_cb(k0);
  plnu = P_nu(k0);
  }

  // Real PT spectrum
    pspt = ploopr(pars[0]) + p1h ;
  // GR-1-loop spectrum with linear growth replaced by modified gravity growth (unscreened)
    psptp = ploopp(pars[0]);

// cbv terms
  psptcbv = sqrt(pspt*plnu);

// reaction  @ k0 with 2halo term = 1 loop spt
  double rsptk0 = (fvt2 * pspt + 2.*fv*fvt*psptcbv + fv2*plnu)/(psptp + p1hp);

  double a1 = ((1.-bigE)*plcb*fvt2);

  double term1 =  fv2*plnu + (plm+p1hp)*rsptk0 -fvt2*(bigE*plcb + p1h);
  double term2 = 2.*sqrt(fv2*(plm+p1hp)*plnu*rsptk0);

  if (bigE>1.) {
    argument =  (term1 + term2)/a1; // 1st root
  }
  else{
    argument = (term1 - term2)/a1; // 2nd root
    }

// For very small modifications, numerics can sometimes generate a negative argument.
// To prevent 'nans' in the log taken below in such cases, we manually set kstar and bigE to their (approximately) GR values.
  if(argument<0.001 || argument>1.){
    kstar = 1.;
    bigE = 1.;
  }
  else{
  kstar = (-k0/log(argument));
  }

      // Make sure the root is actually a root
    double prefac = ((1.-bigE)*exp(-k0/kstar) + bigE);
    double ph = prefac*plcb + p1h; // real
    double php = plm + p1hp; // pseudo
    double hmrk0 = (fvt2*ph +  2.*fv*fvt*sqrt(ph*plnu) + fv2*plnu)/php;

    double dif = fabs(hmrk0 - rsptk0);
    if (dif>1e-4) {
      kstar = 1.;
      bigE = 1.;
          }

}


// Calculate the reaction
double HALO::reaction_nu(double k, double pars[], bool modcamb) const {
  if (k<1e-3) {
    return 1.;
  }

  if ( fabs(pars[2])>1e-6 && !modcamb) {
    react_error_halo("Omega neutrino is not 0 but you have specified a LCDM transfer - set modcamb = true or set Omega_nu = 0");
  }

  double fv  = pars[2]/pars[1];
  double fv2 = pow(fv,2);
  double fvt = 1.-fv;
  double fvt2  = pow(fvt,2);

  // linear terms
    double plm;
    double plcb;
    double plnu;
  if(!modcamb){
    // linear terms
      plm = pow2(linear_growth(k))*P_l(k);
      plcb = plm;
      plnu = 0.;
      }
  else{
    plm = P_l(k);
    plcb = P_cb(k);
    plnu = fabs(P_nu(k)); // take absolute value to avoid negative values in spline when P_nu is very small.
      }

// HM terms
// prefac for 1hr
  double prefac = ((1.-bigE)*exp(-k/kstar) + bigE);
  double ph = prefac*plcb + one_halo(k,pars); // real
  double php = plm + one_halop(k,pars); // pseudo

  // real linear spectrums
 return (fvt2*ph +  2.*fv*fvt*sqrt(ph*plnu) + fv2*plnu)/php;
}


// Initialise everything - for multiple redshifts, must initialise react_init2 separately
void HALO::initialise(double pars[], double extpars[], bool modcamb, bool modg, int model) const{
  IOW iow;
  if ( fabs(pars[2])>1e-6 && !modcamb) {
    react_error_halo("Omega neutrino is not 0 but you have specified a LCDM transfer - set modcamb = true or set Omega_nu = 0");
  }

  iow.initnorm(pars,extpars,model);   // LCDM (or wCDM) linear growth rates
  scol_init(pars, extpars, modcamb,model); // real spherical collapse quantities
  scol_initp(pars, modcamb); // pseudo spherical collapse quantities
  react_init_nu(pars, extpars, modcamb,modg,model); // kstar and mathcal E for reaction
  printf("%s %e %s %e \n", "E:", bigE, "kstar:", kstar);
}


// Initialise everything for multiple redshifts
void HALO::initialise_multiz(double pars[], double extpars[], Spline ploopr, Spline ploopp, bool modcamb, bool modg, int model) const{
  IOW iow;
  if ( fabs(pars[2])>1e-6 && !modcamb) {
    react_error_halo("Omega neutrino is not 0 but you have specified a LCDM transfer - set modcamb = true or set Omega_nu = 0");
  }
  iow.initnorm(pars, extpars, model);   // LCDM (or wCDM) linear growth rates
  scol_init(pars, extpars, modcamb,model); // real spherical collapse quantities
  scol_initp(pars,modcamb); // pseudo spherical collapse quantities
  react_init_nu_multiz(pars,ploopr, ploopp, modcamb, modg); // kstar and mathcal E for reaction for multiple redshifts
  printf("%s %e %s %e \n", "E:", bigE, "kstar:", kstar);
}


// Compute linear power spectrum to be handed to CosmoSIS -  linear_growth must be initialised with scol_initp (needed for pseudo halo model quantities)
double HALO::plinear_cosmosis(double k) const {
  // real linear spectrum
  return pow2(linear_growth(k))*P_l(k);
}

// Linear growth only
double HALO::Lin_Grow(double k) const{
    return linear_growth(k);
}



/* The 1-loop SPT reaction */

double HALO::reaction_spt(double k0, double pars[], double extpars[], bool modcamb, int model) const{
  SPT spt(C, P_cb, epsrel);
  SPT spt2(C, P_l, epsrel);

  double pspt, psptp, p1h, p1hp, plm, plcb, plnu, psptcbv;
  // double psptpnosc,psptpcbv; // used for unscreened approximation - see below
  // fv
  double fv  = pars[2]/pars[1];
  double fv2 = pow(fv,2);
  double fvt = 1.-fv;
  double fvt2  = pow(fvt,2);


// 1-halo terms
  p1h = one_halo(k0,pars); // real
  p1hp = one_halop(k0,pars); // pseudo

// spt terms
if(!modcamb){
  // linear terms
    plm = pow2(linear_growth(k0))*P_l(k0);
    plcb = plm;
    plnu = plm;

// Real PT spectrum
  pspt = spt.PLOOPn2(1, k0, pars, extpars, 1e-3,  model) + p1h ;
// GR-1-loop spectrum with linear growth replaced by modified gravity growth
  psptp = spt.PLOOPn2(4, k0, pars, extpars, 1e-3,  model);
}
else{
  plm = P_l(k0);
  plcb = P_cb(k0);
  plnu = fabs(P_nu(k0));

// Real PT spectrum
  pspt = spt.PLOOPn2_nu(1, k0, pars, extpars, 1e-3, model) + p1h ;

  // Unscreened approx
  // psptpnosc = spt.PLOOPn2_nu(2, k0, pars, extpars, 1e-3, model);
  //
  // psptpcbv = sqrt(psptpnosc*plnu);
  // // pseudo SPT term
  // psptp = (fvt2 * psptpnosc + 2.*fv*fvt*psptpcbv + fv2*plnu);

  // Pseudo
   psptp = spt2.PLOOPn2_nu(2, k0, pars, extpars, 1e-3, model);

  }

// cbv terms
  psptcbv = sqrt(pspt*plnu);

// reaction  @ k0 with 2halo term = 1 loop spt
  return  (fvt2 * pspt + 2.*fv*fvt*psptcbv + fv2*plnu)/(psptp + p1hp);

}



// Pseudo halofit prescription
////////////////////////////////////////////////
///////HALO-FIT FOR NON-SCALE DEP MODELS : takashi et al 2012 fits  1208.2701///////////
///////////////////////////////////////////////

// Uses linear growth initialised in scol_init : linear_growth

double wintcambs_pseudo[3];
static double wintcamb_pseudo(bool modcamb, const PowerSpectrum& P_L, double r){
 int n1 = 3000;
 double sum1 = 0.;
 double sum2 = 0.;
 double sum3 = 0.;
 double anorm = 1./(2.*pow2(M_PI));
 double t,k, pkterm,x,x2,w1,w2,w3,fac;
 for(int i = 1; i<=n1 ; i++){
   t = (i-0.5)/n1;
   k = -1. + 1./t;
   if (modcamb) {
    pkterm = P_L(k)*pow3(k)*anorm;
   }
   else{
   pkterm = pow2(linear_growth(k))*P_L(k)*pow3(k)*anorm;
    }
   x = k*r;
   x2 = pow2(x);
   w1 = exp(-x2);
   w2 = 2.*x2*w1;
   w3 = 4.*x2*(1.-x2)*w1;
   fac = pkterm/k/t/t;
   sum1 += w1*fac;
   sum2 += w2*fac;
   sum3 += w3*fac;
 }
sum1 /= n1;
sum2 /= n1;
sum3 /= n1;
wintcambs_pseudo[0] = sqrt(sum1);
wintcambs_pseudo[1] = -sum2/sum1;
wintcambs_pseudo[2] = -sum2*sum2/sum1/sum1 - sum3/sum1;
return 0.;
}



double phpars_pseudo[13];
void parscamb_pseudo(bool modcamb, const PowerSpectrum& P_L)
{
  double xlogr1 = -2.;
  double xlogr2 = 3.5;
  double rmid = 0.;
  double diff = 1.;
  for (;;){
    rmid = (xlogr2 + xlogr1)/2.;
    rmid = pow(10,rmid);
    wintcamb_pseudo(modcamb, cref(P_L),rmid);
    diff = wintcambs_pseudo[0]-1.;
    if (fabs(diff)<=0.0001) {
      phpars_pseudo[0] = 1./rmid;
      phpars_pseudo[11] = -3.-wintcambs_pseudo[1];
      phpars_pseudo[12] = -wintcambs_pseudo[2];
      break;
    }
   else if(diff>0.0001){
     xlogr1 = log10(rmid);
   }
   else if (diff < -0.0001){
     xlogr2 = log10(rmid);
   }
    if (xlogr2<-1.9999){
      phpars_pseudo[0] = 1000.;
      break;
    }
    else if(xlogr2>3.4999){
      phpars_pseudo[0] = 1000.;
      break;
    }
  }
}

// Initialises components
void HALO::phinit_pseudo(double pars[], bool modcamb)const{
    double scalef = pars[0];
    double omega0 = pars[1];
    parscamb_pseudo(modcamb, cref(P_l));

    double neff = phpars_pseudo[11];
    double curv = phpars_pseudo[12];

   double neff2 = pow2(neff);
   double neff3 = neff2*neff;
   double neff4 = pow2(neff2);


   double ane = 1.5222 + 2.8553*neff + 2.3706*neff2+ 0.9903*neff3 + 0.2250*neff4- 0.6038*curv;
   double bne = -0.5642 + 0.5864*neff + 0.5716*neff2- 1.5474*curv;
   double cne = 0.3698 + 2.0404*neff + 0.8161*neff2 + 0.5869*curv;

   double nue = 5.2105 + 3.6902*neff;
   double acub = pow3(scalef);
   double omgm = (omega0/acub)/(omega0/acub + (1.-omega0));
   double omgl = (1.-omega0)/(omega0/acub + (1.-omega0));

     phpars_pseudo[1] = pow(10.,ane); // an
     phpars_pseudo[2] = pow(10.,bne); // bn
     phpars_pseudo[3] = pow(10.,cne); //cn

     phpars_pseudo[4] = 0.1971 - 0.0843*neff + 0.8460*curv; //gan
     phpars_pseudo[5] = fabs(6.0835 + 1.3373*neff - 0.1959*neff2- 5.5274*curv); // alpha
     phpars_pseudo[6] = 2.0379 - 0.7354*neff + 0.3157*neff2 + 1.2490*neff3 + 0.3980*neff4- 0.1682*curv; //beta

     phpars_pseudo[7] = pow(10.,nue); // nun

     if (fabs(1.-omgm) > 0.01) {
       double f1a = pow(omgm,-0.0732);
       double f2a = pow(omgm,-0.1423);
       double f3a = pow(omgm,0.0725);

       double f1b = pow(omgm,-0.0307);
       double f2b = pow(omgm,-0.0585);
       double f3b = pow(omgm,0.0743);

        double frac = omgl/(1. - omgm);
       phpars_pseudo[8] = frac*f1b + (1.-frac)*f1a;
       phpars_pseudo[9] = frac*f2b + (1.-frac)*f2a;
       phpars_pseudo[10] = frac*f3b + (1.-frac)*f3a;

        }
      else{
        phpars_pseudo[8]=1.;
        phpars_pseudo[9]=1.;
        phpars_pseudo[10]=1.;
      }

  }

// pseudo spectrum from halofit
double HALO::PHALO_pseudo(double k, bool modcamb) const{
  double mypkl;
  if (modcamb) {
    mypkl = P_l(k);
  }
  else{
    mypkl = pow2(linear_growth(k))*P_l(k);
  }

  if (phpars_pseudo[0] == 1000. || k<=0.005) {
     return mypkl;
  }
  else{
  double kcub = pow3(k)/2./pow2(M_PI);
  double deltahp = phpars_pseudo[1]*pow(k/phpars_pseudo[0],phpars_pseudo[8]*3.)/(1.+phpars_pseudo[2]*pow(k/phpars_pseudo[0],phpars_pseudo[9])+pow(phpars_pseudo[3]*phpars_pseudo[10]*k/phpars_pseudo[0],3.-phpars_pseudo[4]));
  double deltah = deltahp/(1.+ phpars_pseudo[7]/pow2(k/phpars_pseudo[0]));
  double deltaplin = kcub*mypkl;
  double deltaq = deltaplin*pow(1.+ deltaplin,phpars_pseudo[6])/(1.+phpars_pseudo[5]*deltaplin)*exp(-(k/phpars_pseudo[0])/4.-pow2(k/phpars_pseudo[0])/8.);
  double deltanl = deltaq + deltah;

  return  deltanl/kcub;
}
}



/* Extra stuff */

/// Halo density profile in real space - NFW for now  but  can edit this easily
double HALO::halo_profileR(double Mvir, double Rvir,  double mycvir, double r) const {
  double rs = Rvir/mycvir;
  double rhos = Mvir/(4.*M_PI*pow3(rs))/(log(1.+mycvir) - mycvir/(1.+mycvir));
  // NFW
  return  rhos/(r/rs * pow2(1.+ r/rs));
}

// profile in real space
static double halo_profile_integrand(const Cosmology& C, const PowerSpectrum& P_l, double Mvir, double Rvir, double mycvir, double k, double r){
      HALO halo(C, P_l, P_l, P_l, P_l, 1e-3);
      return  sin(k*r)*r/k*halo.halo_profileR(Mvir, Rvir, mycvir,r);
  }

/// Halo density profile in k space

// normalise it
double profile_norm = 1.;

void profn_init(const Cosmology& C, const PowerSpectrum& P_l, double Mvir, double Rvir, double mycvir){
  HALO halo(C, P_l,  P_l, P_l, P_l, 1e-3);
  profile_norm =  4.*M_PI*Integrate(bind(halo_profile_integrand, cref(C), cref(P_l), Mvir, Rvir, mycvir, 1e-4, std::placeholders::_1), QMINp, Rvir , 1e-3);
}

// ONE CAN USE THIS FUNCTION IF YOU WANT TO EDIT NFW PROFILE USED IN COMPUTATIONS
double HALO::halo_profileK(double k, double Mvir, double Rvir, double mycvir) const {
     profn_init(cref(C), cref(P_l),  Mvir,  Rvir,  mycvir);
     return 4.*M_PI/profile_norm*Integrate(bind(halo_profile_integrand, cref(C), cref(P_l), Mvir, Rvir, mycvir, k, std::placeholders::_1), QMINp, Rvir , epsrel);

}
