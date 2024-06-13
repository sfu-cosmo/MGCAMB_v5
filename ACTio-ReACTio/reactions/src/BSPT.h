#ifndef BSPT_H
#define BSPT_H

#include "Common.h"

// splines needed for n_effective and linear power spectrum
extern Spline neffehu;
extern Spline neffnw;
extern Spline mylinearps;

class BSPT {
public:
    BSPT(const Cosmology& C, const PowerSpectrum& P_L, real epsrel = 1e-3);

    // Bispectra at tree and 1-loop levels  and fitting formulae in real space 

    // Arguments:
    // k1, k2 wave vector magnitudes, x = k1.k2 (cosine of angle between wave modes)
    // model : 1 GR, 2 fR, 3 nDGP
    // pars - base cosmological parameters : 0 : scale factor , 1 : Omega_{m,0} total matter fraction today
    // extpars - extended model parameters : 0: fr0 or Omega_rc
    // epars - error parameters : 0 : integral abs error, 1 : ODE initial step, 2: ODE abs error, 3: ODE rel error


    /* Analytic expressions (EdS approximation) 1808.01120  */

    // tree level analytic (LCDM + DGP only)
    real Btree(int model, double k1, double k2, real x) const;
    // 1-loop analytic (LCDM and DGP)
    real Bloop(int model, double k1, double k2, double x) const;
    // 1-loop terms except for B123
    real Bloopterms(int model, real k1, real k2, real k3, real x) const;

    /* Numerical expressions 1808.01120  */

    // tree level numerical
    real Btreen(double pars[], double extpars[], int model, double k1, double k2, double x) const;
    // 1-loop numerical
    real Bloopn(double pars[], double extpars[], double epars[], int model, double k1, double k2, double x) const;

    /* Fitting formulae */

    /* Gil-marin/Namikawa fitting formula  : 1805.10567  */
    real Bfit(double k1, double k2, double x) const;
    // Initialize the NW power spectrum, A^{nw,0l}, sigma8 and knl - used in fitting formula
    // internal = True initialises k_nl and sigma8 at the redshift set by inite (linear growth is DGP or GR - D_spt, see SpecialFunctions)
    // internal = False, k_nl and sigma8 are set by user to pars[2], pars[3] respectively.
    void mypspline(double pars[], bool internal = true) const;

private:
    const Cosmology& C;
    const PowerSpectrum& P_L;
    real epsrel;

};

#endif // BSPT_H
