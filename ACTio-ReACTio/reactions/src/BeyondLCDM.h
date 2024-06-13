#ifndef BEYONDLCDM_H
#define BEYONDLCDM_H

#include <vector>

// ZW's edit starts
// maximum number of beyond-LCDM theory parameters: can easily adjust here.
const int maxpars = 50;
void setMGflags(int N_flags, long int flag_arr[]);
void setMGparams(int N_params, double params_arr[]);
void printMGflags();
void printMGparams();
void set_reconstruction_arr(double omega0);
//ZW's edit ends

// Background functions
double bespokehub(double a, double omega0, double extpars[], int model);
double bespokehubd(double a, double omega0, double extpars[], int model);
double bespokehubdd(double a, double omega0, double extpars[], int model);

double HAg(double a, double omega0, double extpars[], int model = 1);
double HA1g(double a, double omega0, double extpars[], int model = 1);
double HA2g(double a, double omega0, double extpars[], int model = 1);
double HA2g2(double a, double omega0, double extpars[], int model = 1);


// Dark scattering friction term (1605.05623)
double myfricF(double a, double omega0, double extpars[], int model = 1);

// 1-loop PT Poisson modification (1606.02520)
//ZW's edit starts
struct dark_energy {
  double grhov_t;
  double gpresv_t;
};
double mu(double a, double k0, double omega0, double omegacb, double extpars[], int model = 1  );
double mgcamb_m(double a, double omega0, double omegacb);
double mgcamb_beta(double a);
void MGCAMB_DarkEnergy(double a, double omega0, double extpars[], dark_energy* de_terms);
//ZW's edit ends

double gamma2(double a, double omega0, double k0, double k1, double k2, double u1, double extpars[] , int model = 1 );
double gamma3(double a, double omega0, double k0, double k1, double k2, double k3, double u1,double u2, double u3, double extpars[], int model = 1 );

// Spherical collapse Poisson modification (1812.05594 - appendix)
//ZW's edit starts
double mymgF(double a, double yh, double yenv, double Rth, double omega0, double omegacb, double extpars[], double delta_initial, int model = 1);
//ZW's edit ends

// DEBUGGING
// void print_mu() ;

// Virial theorem contribution of DE
double WEFF(double a, double omega0, double extpars[], int model = 1 );


#endif // BEYONDLCDM_H
