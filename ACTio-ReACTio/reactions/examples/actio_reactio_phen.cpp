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

//using namespace std;
using std::ifstream;
using std::string;
using std::istringstream;

vector<vector<double> > mytrans;
vector<vector<double> > mytransl;
vector<vector<double> > mypk;


/* Example code to output the reaction and halo spectra for EFTofDE & Phenomenological parametrisations */

int main(int argc, char* argv[]) {

// Which gravity or dark energy model?
// 1: GR  2: f(R) 3: DGP 4: quintessence 5: CPL, 6: Dark scattering with CPL
// 7 -9 : EFTofDE k->infinity limit,  w PPF, unscreened, superscreened respectively
// 10-12: EFTofDE with PPF, unscreened and Error function Phenomenological model respectively.
int mymodel = 12;

// target redshift
double myz = 0.;

// neutrino mass
double mnu = 0.00;  // mv = 0.0ev

// Alpha-basis : Set their scale factor dependence, 1st and 2nd derivatives in reactions/src/BeyondLCDM.cpp alphai_eft, dalphai_eft and ddalphai_eft
// Default is linear scale factor dependence : alpha_i(a) = alphai0 a
double alphak0 = 1.;
double alphab0 = 0.5;
double alpham0 = 0.;
double alphat0 = 0.;
double m2 = 0.;
double c0 = 0.; // where c(a)->c0*c(a) and c(a) is solved using Eq. A16 assuming LCDM H(a).

// Non-linear parametrisation
double myp1 = 0.5; // Screening scale
double myp2 = 0.5; // its mass dependence
double myp3 = 1.0; // its environment dependence
double myp4 = 0.0; // Yukawa suppression scale

//output file name
const char* output = "KGB_Phen_z0.dat";

// number of bins between kmin and kmax
int Nk = 300;
double kmin = 0.01;
double kmax = 5.;

// perform 1-loop corrections? This prompts the calculation of k_star and \mathcal{E}.
bool modg = false;

// Is the transfer being fed to ReACT of the target cosmology?
// If false, the transfer should be LCDM at z=0 and ReACT will rescale P_L using internally computed modified growth - see README.
// If true, the transfer function should be that of the target cosmology at the target redshift (with MG or/and massive neutrinos)
// Note that ReACT does not calculate growth factors for massive neutrino cosmologies and so the target transfer function should be supplied.
bool mgcamb = false;

// Load modified transfer function at with all species at some redshift
// This is a LCDM z=0 transfer with the cosmology specified below
ifstream fin("transfers/EFT/transfer_out0.dat");
// Base cosmology associated with the transfer function
double h  = 0.6737;
double n_s = 0.96605; // spectral index
double Omega_b  = 0.0491989; //  baryon fraction
double Omega_c = 0.265373;

double Omega_nu = mnu/93.14/pow2(h);
double pscale = 0.05;
double As = 2.09681e-9;

double Omega_m = Omega_c + Omega_b + Omega_nu;

// Load in the transfer data
string line;
    while (getline(fin, line)) {      // for each line
            vector<double> lineData;           // create a new row
            double val;
            istringstream lineStream(line);
            while (lineStream >> val) {          // for each value in line
                    lineData.push_back(val);           // add to the current row
            }
            mytrans.push_back(lineData);         // add row to allData
    }


// Populate arrays and normalise to 1 at large scales for Cosmology class input
int Nkr = mytrans.size();
int* Nkt = &Nkr;
array Tm(*Nkt);
array Tcb(*Nkt);
array Tnu(*Nkt);
array ki(*Nkt);

          for(int i = 0; i< Nkr; i++){
                  ki[i] = mytrans[i][0];
                  Tm[i] = mytrans[i][6];
                  Tcb[i] = mytrans[i][7];
                  Tnu[i] = mytrans[i][6]; // should be adjusted if this column is all 0s ....
              }

// integration error
real epsrel = 1e-3;

// number of mass bins between 5<Log10[M]<18 for the spherical collapse solver to iterate over
double massb = 50.;

// store params for passing into React functions
double pars[7];
    pars[0] = 1./(1.+myz); //  scale factor
    pars[1] = Omega_m;
    pars[2] = Omega_nu ; //  extra param
    pars[5] = massb; // number of mass bins between 5<Log10[M]<18

// store beyond-LCDM parameters
double extpars[maxpars];
    // Linear and background
    extpars[0] = alphak0;
    extpars[1] = alphab0;
    extpars[2] = alpham0;
    extpars[3] = alphat0;
    extpars[4] = m2;
    extpars[19] = c0;
    // non-linear
    extpars[5] = myp1;
    extpars[6] = myp2;
    extpars[7] = myp3;
    extpars[8] = myp4;

    // Load cosmology classes
    Cosmology Cm(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tm);
    Cosmology Ccb(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tcb);
    Cosmology Cnu(h, n_s, Omega_m, Omega_b, As, pscale, ki, Tnu);

    // Get linear P(k) from input transfer
    LinearPS_as P_l(Cm, 0.);
    LinearPS_as P_cb(Ccb, 0.);
    LinearPS_as P_nu(Cnu, 0.);
    LinearPS_as P_cbl(Ccb, 0.);

    IOW iow;

    // Load halo class witth all linear P(k)
    HALO halo(Cm, P_l, P_cb, P_nu, P_cbl, epsrel);


//initialise spherical collapse quantities and reaction quantities
halo.initialise(pars,extpars,mgcamb,modg,mymodel);
// initialise halofit parameters
halo.phinit_pseudo(pars,mgcamb);

/* Output section */

/* Open output file */
FILE* fp = fopen(output, "w");

real p1,p2,p3,p4,p5,p6;


 for(int i =0; i <Nk;  i ++) {

      real k = kmin* exp(i*log(kmax/kmin)/(Nk-1));

      p1 = halo.plinear_cosmosis(k); // Linear spectrum
      p2 = halo.reaction_nu(k,pars, mgcamb); // halo model reaction
      p3 = halo.PHALO_pseudo(k,mgcamb); // halofit pseudo spectrum

      printf("%e %e %e %e %e  \n", k, p1,p2, p3,p3*p2); // output to terminal
      fprintf(fp,"%e %e %e %e %e \n", k, p1,p2,p3,p3*p2); // output to file : k , P_linear, R, P_pseudo, P_nl = RxP_pseudo

}

	/*close output file*/
    fclose(fp);
    return 0;
}
