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
#include <Copter/BeyondLCDM.h>

using namespace std;

using std::ifstream;
using std::string;
using std::istringstream;


/* Example code to output the halo model powerspectrum for modified gravity */

int main(int argc, char* argv[]) {

  // Which gravity or dark energy model?
  // 1: GR  2: Hu-Sawicki f(R) 3: DGP 4: quintessence 5: CPL, 6: Dark scattering with CPL
  // 7 -9 : EFTofDE k->infinity limit,  w PPF, unscreened, superscreened respectively
  // 10-12: EFTofDE with PPF, unscreened and Error function Phenomenological model respectively.
  int mymodel = 2;

  // chosen redshift
     double myz = 0.;
  // chosen omega_matter (total)
    double omegam0 = 0.3072;
  // chosen omega_nu (total)
    double omeganu0 = 0.;

  // Value of f_{R0}
    double fr0 = 1e-5;

  //output file name
  const char* output = "halo_ps.dat";

  // number of bins between kmin and kmax
  int Nk =100;
  double kmin = 0.01;
  double kmax = 10.;


  // perform 1-loop corrections? This prompts the calculation of k_star and \mathcal{E}.
  bool modg = true;

  // Is the transfer being fed to ReACT of the target cosmology?
  // If false, the transfer should be LCDM at z=0 and ReACT will rescale P_L using internally computed modified growth - see README.
  // If true, the transfer function should be that of the target cosmology at the target redshift (with MG or/and massive neutrinos)
  // Note that ReACT does not calculate growth factors for massive neutrino cosmologies and so the target transfer function should be supplied.
  bool mgcamb = false;

   // cosmo file for P_L(k)
   const char* cstr = "transfers/Matteo_fr";

   // Keep it z=0 to keep Copter's Growth @ 1
   real z = 0;
   // Relative error in magnitude integrations
   real epsrel = 1e-3;
   // Initialise all classes
   Cosmology C(cstr);
   LinearPS P_l(C, z);
   HALO halo(C, P_l,P_l,P_l,P_l, epsrel);
   SPT spt(C, P_l, epsrel);
   IOW iow;


   // extended model parameters
   double extpars[maxpars]; // Currently maxed out at 20 extra params
   extpars[0] = fr0; // fr0 for f(R)

   // Base parameters
    double pars[7];
    pars[0] = 1./(1.+myz); // scale factor
    pars[1] = omegam0;  // chosen omega_matter (total)
    pars[2] = omeganu0;  // chosen omega_matter (total)
    pars[5] = 50.; // number of mass bins for spherical collapse solver
    pars[6] = 0.0; // omega_neutrinos

//initialise spherical collapse quantities and reaction quantities
halo.initialise(pars,extpars,mgcamb,modg,mymodel);

/* Output section */

/* Open output file */
FILE* fp = fopen(output, "w");

real p1,p2,p3,p4,p5,p6;


 for(int i =0; i < Nk;  i ++) {

  real k = kmin * exp(i*log(kmax/kmin)/(Nk-1));

      p1 =  halo.one_halo(k, pars); // one halo term : real
      p2 =  halo.one_halop(k, pars); // one halo term : pseudo
      p3 =  halo.reaction_nu(k, pars, mgcamb); // halo model reaction

     printf("%d %e %e %e %e \n", i, k, p1,p2,p3); // print to terminal
     fprintf(fp,"%e %e %e %e \n", k, p1,p2, p3); // print to file

}

	/*close output file*/
    fclose(fp);
    return 0;
}
