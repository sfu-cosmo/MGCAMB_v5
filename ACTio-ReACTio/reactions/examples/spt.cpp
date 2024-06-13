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
vector<vector<double> > allDatamult;

/* Example code to output the 1-loop powerspectrum for modified gravity in real space : exact numerical and EdS approximation */

int main(int argc, char* argv[]) {

    // Which gravity or dark energy model?
    // 1: GR  2: f(R) 3: DGP 4: quintessence 5: CPL
    int mymodel = 3;

    // chosen redshift
       double myz = 1.;
    // chosen omega_matter (total)
       double omega0 = 0.281;
    // MG parameter (for f(R) this is |fr0| and for nDGP this is Omega_rc)
       double mgpar = 0.25;


    //output file name
    const char* output = "dgp_eds_vs_full.dat";

    // Number of k-modes you wish to output between kmin and kmax
      int Nk = 50;
      real kmin = 0.01;
      real kmax = 0.5;

      // error in 1-loop integrals
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
  double pars[7];
  pars[0] = 1./(1.+myz);
  pars[1] =  omega0;
  pars[2] =  0.;

  // extended model parameters
  double extpars[maxpars]; // Currently maxed out at 20 extra params
  extpars[0] = mgpar;

// normalise the growth + calculate all growth factors
  iow.inite(pars,extpars,mymodel);

/* Output section */

 /* Open output file */
 FILE* fp = fopen(output, "w");

real p1,p2,p3,p4,p5,p6,p7,p8,p9;

for(int i =0; i <Nk;  i ++) {

    real k = kmin * exp(i*log(kmax/kmin)/(Nk-1));

/* for more quantities check out the SPT.h file in the src directory */

  // exact linear and 1-loop matter, cross and velocity spectra
  p1 = spt.PLOOPn2(0, k, pars, extpars, err, mymodel); // linear spectrum
  p2 = spt.PLOOPn2(1, k, pars, extpars, err, mymodel); // 1-loop dd spectrum
  p3 = spt.PLOOPn2(2, k, pars, extpars, err, mymodel); // 1-loop dt spectrum
  p4 = spt.PLOOPn2(3, k, pars, extpars, err, mymodel); // 1-loop tt spectrum

  // EdS approximated equivalents
  p5 = pow2(F1_nk/dnorm_spt)*P_l(k);
  p6 = spt.PLOOP(k,4);
  p7 = spt.PLOOP(k,5);
  p8 = spt.PLOOP(k,6);

   printf("%e %e %e %e %e %e %e %e %e \n", k, p1,p2,p3,p4,p5,p6,p7,p8); // print to terminal
   fprintf(fp,"%e %e %e %e %e %e %e %e %e  \n", k, p1,p2,p3,p4,p5,p6,p7,p8 ); // print to file

  }


	/*close output file*/
    fclose(fp);
    return 0;
}
