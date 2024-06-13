#ifndef SPT_H
#define SPT_H

#include "Common.h"

/******************************************************************************
 * SPT
 *
 * 1-loop standard perturbation theory.  Indices (a,b) are used to indicate
 * either the density-density, density-velocity, or velocity-velocity spectra,
 *   P_{11} <--> P_{\delta\delta}
 *   P_{12} <--> P_{\delta\theta}
 *   P_{22} <--> P_{\theta\theta}
 * whereas the labels "22" or "13" refer to different terms in the perturbation
 * series.
 ******************************************************************************/

extern double phpars[13];

class SPT {
public:
    SPT(const Cosmology& C, const PowerSpectrum& P_L, real epsrel = 1e-3);


    // Powerspectra at 1-loop level and fitting formulae in real and redshift space

    // pars: base parameters. Currently used: 0: scale factor, 1: total matter fraction today, 2: total massive neutrino fraction today
    // extpars: extended parameters. extpars[0] = Omega_rc for nDGP, fr0 for f(R), w0 for CPL, extpars[1] = wa for CPL
    // model selects MG or DE model (1 = GR, MG: 2 = f(R), 3 = DGP, DE models: 4 = quintessence, 5 = CPL, 6 = HyP)


      /* REAL SPACE */

    	/* P_{ab}(k) 1-loop: analytical (EdS) LCDM and nDGP */
    	/* The division by Dl is my own normalisation of P_L - it requires z=0 always in spt2.cpp */
      // values of a :
      // 1: P_dd, 2: P_dt , 3: P_tt ; LCDM
    	// 4:  P_dd, 5: P_dt , 6: P_tt ; nDGP
      real PLOOP(real k, int a) const;

    	/* P_{ab}(k) 1-loop : numerical for arbitrary model of gravity for massless/no neutrinos- see  1606.02520  */
      // values of a :
      // 0: P_dd linear
      // 1: P_dd, 2: P_dt , 3: P_tt for mg @ 1-loop level
      // 4: P_dd pseudo @ 1-loop level - used for halo model reactions (see HALO.cpp)
      real PLOOPn2(int a, double k, double pars[], double extpars[], double err, int model = 1) const;

      /* P_{ab}(k) 1-loop : numerical for arbitrary model of gravity with massive neutrinos - kernel dependent */
      // 0: P_dd for mg linear
      // 1: P_dd for mg @ 1-loop
      // 2: P_dd pseudo @ 1-loop level - used for halo model reactions (see HALO.cpp)
      real PLOOPn2_nu(int a, double k, double pars[], double extpars[], double err, int model = 1) const;



      // initialise p_loop values over redshifts[] at k=k0 for k_star in reaction code (HALO.cpp) for massless neutrino cosmologies
      // ploopr and ploopp hold initialised values
      // redshifts holds the redshifts we want to spline over
      // noz is the number of redshifts
      void ploop_init(double ploopr[], double ploopp[], double redshifts[], int noz, double pars[], double extpars[], double k0, int model = 1);

      // initialise p_loop values over redshifts[] at k=k0 for k_star in reaction code (HALO.cpp) for massive neutrino cosmologies
      // first 6 arrays hold values of 1-loop term integrands to be integrated - see function in SPT.cpp for more details
      void ploop_init_nu(double pkz0[], double pkz1[][50], double pkz2[][50*32], double pkz0p[], double pkz1p[][50], double pkz2p[][50*32], double ploopr[], double ploopp[], double redshifts[], int noz, double pars[], double extpars[], double k0, int model = 1);

      // halo fit coefficient initialisation
      void phinit(double scalef, double omega0) const;

      // halofit model with P_L input given by (D_spt/dnorm_spt)^2 * P_L(z=0)
      // See SpecialFunctions for initialisation of D_spt and dnorm_spt
      double PHALO(double k) const;



    /*Redshift Space power spectra */

    /* Analytic expressions - for scale independent models - EdS approximation*/

      // bl is galaxy bias
      // u is the cosine of the angle between k and the line of sight

    	/* P(k,u) KAISER  */
    	// a =1 : LCDM
    	// a =2 : nDGP
      real PRSD(real k, real u, real bl, int a) const;


    	/* P(k,u) TNS model 1006.0699 */
      // sigma_v is the damping factor
      // u is the cosine of the angle between k and LOS
    	// a = 1 : LCDM
    	// a = 2 : nDGP
      real PTNS(real k, real u, real bl, real sigma_v, int a) const;

    	/*Multipoles*/

    	/*  LCDM and scale independent model Kaiser Multipoles with DFoG term (set sigma_v=0 for linear kaiser multipoles) */

    	// a = 1 : Monopole
    	// a = 2 : Quadrupole
    	// a = 3 : Hexdecapole
    	real KASM(real k, real bl, real sigma_v, int a) const ;

    	/*  LCDM  TNS  Multipoles */
    	// a = 1 : Monopole
    	// a = 2 : Quadrupole
    	// a = 3 : Hexdecapole

      //linear bias
    	real PTNSM(real k, real bl, real sigma_v, int a) const;

      // qbias - see Eq.23 of 1507.01592
      real PTNSMq(real k, double barr[], real sigma_v, int a) const;

      // lagrangian bias - see model of 1607.03149
      real PTNSMl(real k, double barr[], real sigma_v, int a) const;

    	/*  nDGP TNS  Multipoles  */
    	// a = 1 : Monopole
    	// a = 2 : Quadrupole
    	// a = 3 : Hexdecapole

      //linear bias
      real PTNSMnDGP(real k, real bl, real sigma_v, int a) const ;

      // qbias
      real PTNSMnDGPq(real k, double barr[], real sigma_v, int a) const;

      // lagrangian bias terms
      real Lag_bias(int a, real k, real bias[]) const;



      // Linear velocity dispersion
      real sigmav_init(double pars[], double extpars[], int model = 1) const;


      // modified gravity 1-loop, TNS and Kaiser model - numerically calculated
      double PRSD_mg(int a, int b, double pars[], double extpars[], double rsdpars[], double bias[],  double k, double err, int model = 1) const;




// individual loop terms (EdS - analytic )

	//LCDM 1-loop terms

    /* $P_{ab}^{(2,2)}(k)$ */
    real P22(real k, int a=1, int b=1) const;

    /* $P_{ab}^{(1,3)}(k)$ */
    real P13(real k, int a=1, int b =1) const;

    /* $P_{\delta\delta}^{(2,2)}(k)$ */
    real P22_dd(real k) const;

    /* $P_{\delta\theta}^{(2,2)}(k)$ */
    real P22_dt(real k) const;

    /* $P_{\theta\theta}^{(2,2)}(k)$ */
  	real P22_tt(real k) const;

    /* $P_{\delta\delta}^{(1,3)}(k)$ */
    real P13_dd(real k) const;

    /* $P_{\delta\theta}^{(1,3)}(k)$ */
    real P13_dt(real k) const;

    /* $P_{\theta\theta}^{(1,3)}(k)$ */
	  real P13_tt(real k) const ;


	//DGP 1-loop terms

 /* $P_{ab}^{(2,2)}(k)$ */
    real P22D(real k, int a) const;

    /* $P_{ab}^{(1,3)}(k)$ */
    real P13D(real k, int a) const;

    /* $P_{\delta\delta}^{(2,2)}(k)$ */
    real P22D_dd(real k) const;

    /* $P_{\delta\theta}^{(2,2)}(k)$ */
    real P22D_dt(real k) const;

    /* $P_{\theta\theta}^{(2,2)}(k)$ */
     real P22D_tt(real k) const;

    /* $P_{\delta\delta}^{(1,3)}(k)$ */
    real P13D_dd(real k) const;

    /* $P_{\delta\theta}^{(1,3)}(k)$ */
    real P13D_dt(real k) const;

    /* $P_{\theta\theta}^{(1,3)}(k)$ */
     real P13D_tt(real k) const ;


	/* A+B+DGP terms analytical*/
  // a = 1 : Monopole
  // a = 2 : Quadrupole
  // a = 3 : Hexdecapole
	real AB_lcdm(real k, real bl, real sigmav, int a) const;
  real AB_dgp(real k, real bl, real sigmav, int a) const;

  // for 2d spectra
  real AB_mu_lcdm(real k, real bl, real u, int a) const;
  real AB_mu_dgp(real k, real bl, real u, int a) const;



private:
    const Cosmology& C;
    const PowerSpectrum& P_L;
    real epsrel;

};

#endif // SPT_H
