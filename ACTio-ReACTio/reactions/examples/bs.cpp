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

#include <omp.h>

#include <Copter/HALO.h>
#include <Copter/Cosmology.h>
#include <Copter/LinearPS.h>
#include <Copter/SpecialFunctions.h>
#include <Copter/SPT.h>
#include <Copter/BSPT.h>
#include <Copter/BeyondLCDM.h>

using namespace std;

/* Example code to output the real space bispectrum for modified gravity */

int main() {


  // Which gravity or dark energy model?
  // 1: GR  2: f(R) 3: DGP
  int mymodel = 3;

  // Keep it z=0 to keep Copter's Growth @ 1
  real z = 0;
  // Relative error in magnitude integrations
  real epsrel = 1e-3;

  // cosmo file for P_L(k)
  const char* cstr = "transfers/wmap9b";

  // Initialise all classes
  Cosmology C(cstr);
  LinearPS P_l(C, z);
  BSPT bspt(C, P_l, epsrel);
  SPT spt(C, P_l, epsrel);
  IOW iow;

  // chosen redshift
   double myz = 0.541;
  // chosen omega_matter (total)
   double omega0 = 0.279;
  // MG parameter (for f(R) this is |fr0| and for nDGP this is Omega_rc)
   double mgpar = 0.25;


   // base parameter values
   double pars[7];
   pars[0] = 1./(1.+myz); // scale factor
   pars[1] =  omega0; // total matter fraction
   pars[2] =  0.0; // total massive neutrino fraction

   // extended model parameters
   double extpars[maxpars]; // Currently maxed out at 20 extra params
   extpars[0] = mgpar;


   // computation error parameters
   double epars[4];
   epars[0] = 1e-2; // 1-loop integral abs error
   epars[1] = 1e-3; // ODE initial step size
   epars[2] = 1e-3; // ODE abs error
   epars[3] = 1e-3; // ODE rel error


// initialise GR or DGP lin growth for PS normalisation
if (mymodel==3) {
  iow.inite2(pars,extpars,mymodel);
}
else{
  iow.inite2(pars,extpars,mymodel);
}

// initialise fitting function params (n_effective, k_nl, sigma_8)
bspt.mypspline(pars,true);

/* Output section */

//output file name
const char* output = "DGP_orc0p25_z0p54.dat";

/* Open output file */
FILE* fp = fopen(output, "w");

real b1,b2,b3,b4;

// Number of k-modes you wish to output between kmin and kmax
int Nk =10;
double kmin = 0.01;
double kmax = 0.5;

 for(int i =0; i < Nk;  i ++) {

  real k =  kmin * exp(i*log(kmax/kmin)/(Nk-1));

// bispectra for equilateral config

 double mu = -0.5; // cosine of angle between k and k1
 double k1 = k;


 b1 = bspt.Btreen(pars,extpars,mymodel,k,k1,mu); // numerically computed tree level
 b2 = bspt.Bloopn(pars,extpars,epars,mymodel,k,k1,mu); //numerically computed
 b3 = bspt.Bloop(3,k,k1,mu); //EdS approximated GR or DGP spectrum
 b4 = bspt.Bfit(k,k1,mu); // Gil-Marin/Namikawa et al fitting formula (GR and DGP only)


 printf("%e %e %e %e %e %e \n",k,b1,b2,b3,b4,b1/b2); // print to terminal
 fprintf(fp,"%e %e %e %e %e \n",k,b1,b2,b3,b4); // print to terminal


}

	/*close output file*/
    fclose(fp);
    return 0;
}
