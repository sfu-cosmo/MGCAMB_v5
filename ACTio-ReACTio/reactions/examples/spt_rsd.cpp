#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <sys/stat.h>
#include <cmath>


#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <time.h>
#include <cmath>

#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/SPT.h>
#include <Copter/GrowthFunction.h>
#include <Copter/SpecialFunctions.h>
#include <Copter/NoWigglePS.h>
#include <Copter/BeyondLCDM.h>

using namespace std;

/* Example code to output the 1-loop powerspectrum for modified gravity in real and redshift space*/

int main(int argc, char* argv[]) {

  // Which gravity or dark energy model?
  // 1: GR  2: f(R) 3: DGP 4: quintessence 5: CPL, 6: Dark scattering with CPL
  // 7 -9 : EFTofDE k->infinity limit,  w PPF, unscreened, superscreened respectively
  // 10-12: EFTofDE with PPF, unscreened and Error function Phenomenological model respectively.
    int mymodel = 3;

    // chosen redshift
      double myz = 0.;
    // chosen omega_matter (total)
      double omega0 = 0.281;
     // MG parameter (for f(R) this is |fr0| and for nDGP this is Omega_rc)
      double mgpar = 0.25;

      //output file name
      const char* output = "DGP_P024_z0.dat";

      // number of bins between kmin and kmax
      int Nk =30;
      double kmin = 0.001;
      double kmax = 0.3;

      // choose RSD model :
      // 0: Kaiser
      // 1: TNS with q-bias (see 1507.01592 for example)
      // 2: TNS with Lagrangian bias (incomplete!)
      // 3: SPT 1-loop RSD spectrum
      double rsd_model = 1;

      // RSD parameters
      double rsdpars[3];

      // sigma_v for TNS (rsd_model = 1,2)
      double mysigv = 5.;

      // Lagrangian bias params (See 1607.03149 for example) OR Q-bias params (see PRSD_mg in src/SPT.cpp).
      double bias[3];
      bias[0] = 1.; // b1
      bias[1] = 0.; // b2
      bias[2] = 0.; // N

      // absolute error in RSD integrals
      double err = 1e-2;

    // cosmo file for P_L(k)
    const char* cstr = "transfers/dgp";

    // Keep it z=0 to keep Copter's Growth @ 1
    real z = 0;
    // Relative error in magnitude integrations
    real epsrel = 1e-3;

    // Initialise all classes
    Cosmology C(cstr);
    LinearPS P_l(C, z);
    SPT spt(C, P_l, epsrel);
    IOW iow;

    // base parameter values
    double pars[10];
    pars[0] = 1./(1.+myz); // scale factor
    pars[1] =  omega0; // Omega_{m,0}
    pars[2] =  0.; // Omega_{nu,0}

    // ODE solver initial step
    pars[7] = 1e-4;
    // ODE solver relative error
    pars[8] = 1e-3;
    // ODE solver absolute error
    pars[9] = 1e-3;

    // extended model parameters
    double extpars[maxpars]; // Currently maxed out at 20 extra params
    extpars[0] = mgpar;

    // initialise LCDM or wCDM lin growth for PS normalisation
    iow.initnorm(pars,extpars,mymodel);

    if (rsd_model==3) {
      //linear vel dispersion for 1-loop (rsd_model = 3)
      rsdpars[0] = spt.sigmav_init(pars, extpars, mymodel); // calculates linear theory dispersion
    }
    else{
      rsdpars[0] = mysigv;
    }

    /* Output section */

    /* Open output file */
    FILE* fp = fopen(output, "w");

    real p1,p2,p3,p4,p5,p6,p7,p8,p9;


for(int i =0; i <Nk;  i ++) {

     real k = kmin * exp(i*log(kmax/kmin)/(Nk-1));

    /* for more quantities check out the SPT.h file in the src directory */
      p1 = spt.PRSD_mg(rsd_model, 1, pars, extpars, rsdpars, bias, k, err, mymodel); // P0 SPT
      p2 = spt.PRSD_mg(rsd_model, 2, pars, extpars, rsdpars, bias, k, err, mymodel); // P2 SPT
      p3 = spt.PRSD_mg(rsd_model, 3, pars, extpars, rsdpars, bias, k, err, mymodel); // P4 SPT

     printf("%e %e %e %e  \n", k, p1,p2,p3); // print to terminal
     fprintf(fp,"%e %e %e %e  \n", k, p1,p2,p3); // print to file

}


	/*close output file*/
    fclose(fp);
    return 0;
}
