#ifndef HALO_H
#define HALO_H
#include "Common.h"


// limits on mass exponents
// The maximum mass is chosen sufficiently high so that we get convergence of the 1-halo term
// The minimum mass is chosen so that we always find a solution for nu = delta_nl/sigma (M) = 1 used in the virial concentration (see  eq 46 of 1812.05594)
// The current values were chosen to satisfy these properties for a wide range of test cases in f(R) gravity where nu is mass dependent.
// Ideally you want as small a range as possible so that you are able to solve spherical collapse less often
// We sample pars[5] (or mass_loop parameter in cosmosis pipeline) times in this range - see below
const double Mmax = 20.;
const double Mmin = 5.;


const double EHMR = 0.01; // value at which to calculate \mathcal{E}
const double K0HMR = 0.06; // value at which to calculate k_star

// For massless neutrinos or no neutrinos, user should set P_cb = P_nu = P_l (total matter spectrum)
class HALO {
public:
  HALO(const Cosmology& C, const PowerSpectrum& P_l, const PowerSpectrum& P_cb,const PowerSpectrum& P_nu,const PowerSpectrum& P_cbl, real epsrel = 1e-4);

// Initialisation of spherical collapse quantities & halo model based functions
// pars: base parameters. Currently used:
// 0: scale factor
// 1: total matter fraction today
// 2: total massive neutrino fraction today
// 3,4 : currently unused
// 5: number of mass bins over which to spline spherical collapse quantities

// extpars: extended parameters. extpars[0] = Omega_rc for nDGP, fr0 for f(R), w0 for CPL, extpars[1] = wa for CPL
// model selects MG or DE model (1 = GR, MG: 2 = f(R), 3 = DGP, DE models: 4 = quintessence, 5 = CPL, 6 = HyP)

// modcamb
// true: we are using the real cosmology transfer function or power spectrum input
// false: we are using LCDM linear transfer function or power spectrum at z=0 input and need to rescale

// modg
// true: do 1-loop computations for k_star
// false: don't do 1-loop computations, sets \mathcal{E}=1 - LCDM value

// initialise spherical collapse quantities
      int scol_init(double pars[], double extpars[], bool modcamb = false , int model = 1) const;
      int scol_initp(double pars[], bool modcamb = false) const;

// halo model components
      double rvirial(double Mvir, double pars[]) const;
      double cvirial(double Mvir, double acol) const;
      double nvirial(double Mvir, double omega0) const;

// pseudo components
      double rvirialp(double Mvir, double pars[]) const;
      double cvirialp(double Mvir, double acol) const;
      double nvirialp(double Mvir, double omega0) const;

// standard FT of NFW
      double halo_profileK2(double k,double Rvir, double mycvir) const; // standard FT of NFW

// halo model power spectra terms
      double one_halo(double k, double pars[]) const;
      double one_halop(double k, double pars[]) const;


// reactions with massive neutrinos
      void react_init_nu(double pars[], double extpars[], bool modcamb = false, bool modg = true, int model = 1) const;

// for multiple redshifts - P_1-loop(k0,z) is initialised as a spline beforehand - see SPT.h : ploop_init or ploop_init_nu
      void react_init_nu_multiz(double pars[], Spline ploopr, Spline ploopp, bool modcamb=false , bool modg = true) const;

// reaction function
      double reaction_nu(double k, double pars[], bool modcamb = true) const;
// 1-loop SPT reaction
      double reaction_spt(double k0, double pars[], double extpars[], bool modcamb = false, int model = 1) const;


// Initialise everything for single redshift
      void initialise(double pars[], double extpars[], bool  modcamb = false, bool modg = true, int model = 1) const;

// Optimised initialisation for multiple redshifts
// requires that 1-loop spectra P_{real}(k0,z) & P_{pseudo}(k0,z) as a function of z is initialised
      void initialise_multiz(double pars[], double extpars[], Spline ploopr, Spline ploopp, bool  modcamb = false, bool modg = true, int model = 1) const;

// Linear spectrum for CosmoSIS
      double plinear_cosmosis(double k) const;

// Linear growth
      double Lin_Grow(double k) const;

      // halofit pseudo spectrum
      // requires running of scol_init to initialise linear growth factor spline D(k)
      double PHALO_pseudo(double k, bool modcamb = false) const;
      // initialiser for halofit quantities - pars is as in all other functions, only call once for all k but at fixed scale factor(a = pars[0])
      void phinit_pseudo(double pars[], bool modcamb = false)const;

// extras
    // density profile in real and fourier space (see expression in .cpp to edit to desired profile )
          double halo_profileR(double Mvir, double Rvir,  double mycvir, double r) const;
    // FT of above profile
          double halo_profileK(double k, double Mvir, double Rvir, double mycvir) const; // used for general profiles rho(r)


private:
    const Cosmology& C;
    const PowerSpectrum& P_l;
    const PowerSpectrum& P_cb;
    const PowerSpectrum& P_nu;
    const PowerSpectrum& P_cbl;
    real epsrel;


} ;



#endif // HALO
