# modified gravity parameters

from ctypes import c_int, c_double, c_bool
from .baseconfig import F2003Class, fortran_class

class BaseModGravityModel(F2003Class):
    """
    Abstract base class for modified gravity models
    """

    _fields_ = [
        ("MG_wrapped", c_bool, "Run the code in Fortran or Python? If true, uses Python"),
    ]


@fortran_class
class ModGravityModel(BaseModGravityModel):

    # TODO: write one line documentation about these params

    _fields_ = [

        ("MG_flag", c_int, "MG_flag"),
        ("GRtrans", c_double, "GRtrans"),
        ("pure_MG_flag", c_int, "pure_MG_flag"),
        ("alt_MG_flag", c_int, "alt_MG_flag"),
        ("QSA_flag", c_int, "QSA_flag"),
        ("CDM_flag", c_int, "CDM_flag"),
        ("muSigma_flag", c_int, "muSigma_flag"),
        ("mugamma_par", c_int, "mugamma_par"),
        ("B1", c_double, "B1"),
        ("lambda1_2", c_double, "lambda1_2"),
        ("B2", c_double, "B2"),
        ("lambda2_2", c_double, "lambda2_2"),
        ("ss", c_double, "ss"),
        ("E11", c_double, "E11"),
        ("E22", c_double, "E22"),
        ("ga", c_double, "ga"),
        ("nn", c_double, "nn"),
        ("musigma_par", c_int, "musigma_par"),
        ("mu0", c_double, "mu0"),
        ("sigma0", c_double, "sigma0"),
        ("QR_par", c_int, "QR_par"),
        ("MGQfix", c_double, "MGQfix"),
        ("MGRfix", c_double, "MGRfix"),
        ("Qnot", c_double, "Qnot"),
        ("Rnot", c_double, "Rnot"),
        ("sss", c_double, "sss"),
        ("Linder_gamma", c_double, "Linder_gamma"),
        ("B0", c_double, "B0"),
        ("beta_star", c_double, "beta_star"),
        ("a_star", c_double, "a_star"),
        ("xi_star", c_double, "xi_star"),
        ("beta0", c_double, "beta0"),
        ("xi0", c_double, "xi0"),
        ("DilS", c_double, "DilS"),
        ("DilR", c_double, "DilR"),
        ("F_R0", c_double, "F_R0"),
        ("FRn", c_double, "FRn"),
        ("DE_model", c_int, "DE_model"),
        ("w0DE", c_double, "w0DE"),
        ("waDE", c_double, "waDE"),
        ("MGDE_pert", c_bool, "MGDE_pert"),
		("MGCAMB_Mu_idx_1", c_double, "MGCAMB_Mu_idx(1)"),
		("MGCAMB_Mu_idx_2", c_double, "MGCAMB_Mu_idx(2)"),
		("MGCAMB_Mu_idx_3", c_double, "MGCAMB_Mu_idx(3)"),
		("MGCAMB_Mu_idx_4", c_double, "MGCAMB_Mu_idx(4)"),
		("MGCAMB_Mu_idx_5", c_double, "MGCAMB_Mu_idx(5)"),
		("MGCAMB_Mu_idx_6", c_double, "MGCAMB_Mu_idx(6)"),
		("MGCAMB_Mu_idx_7", c_double, "MGCAMB_Mu_idx(7)"),
		("MGCAMB_Mu_idx_8", c_double, "MGCAMB_Mu_idx(8)"),
		("MGCAMB_Mu_idx_9", c_double, "MGCAMB_Mu_idx(9)"),
		("MGCAMB_Mu_idx_10", c_double, "MGCAMB_Mu_idx(10)"),
		("MGCAMB_Mu_idx_11", c_double, "MGCAMB_Mu_idx(1)1"),
		("MGCAMB_Sigma_idx_1", c_double, "MGCAMB_Sigma_idx(1)"),
		("MGCAMB_Sigma_idx_2", c_double, "MGCAMB_Sigma_idx(2)"),
		("MGCAMB_Sigma_idx_3", c_double, "MGCAMB_Sigma_idx(3)"),
		("MGCAMB_Sigma_idx_4", c_double, "MGCAMB_Sigma_idx(4)"),
		("MGCAMB_Sigma_idx_5", c_double, "MGCAMB_Sigma_idx(5)"),
		("MGCAMB_Sigma_idx_6", c_double, "MGCAMB_Sigma_idx(6)"),
		("MGCAMB_Sigma_idx_7", c_double, "MGCAMB_Sigma_idx(7)"),
		("MGCAMB_Sigma_idx_8", c_double, "MGCAMB_Sigma_idx(8)"),
		("MGCAMB_Sigma_idx_9", c_double, "MGCAMB_Sigma_idx(9)"),
		("MGCAMB_Sigma_idx_10", c_double, "MGCAMB_Sigma_idx(10)"),
		("MGCAMB_Sigma_idx_11", c_double, "MGCAMB_Sigma_idx(11)"),
		("Funcofw_1", c_double, "Funcofw(1)"),
		("Funcofw_2", c_double, "Funcofw(2)"),
		("Funcofw_3", c_double, "Funcofw(3)"),
		("Funcofw_4", c_double, "Funcofw(4)"),
		("Funcofw_5", c_double, "Funcofw(5)"),
		("Funcofw_6", c_double, "Funcofw(6)"),
		("Funcofw_7", c_double, "Funcofw(7)"),
		("Funcofw_8", c_double, "Funcofw(8)"),
		("Funcofw_9", c_double, "Funcofw(9)"),
		("Funcofw_10", c_double, "Funcofw(10)"),
		("Funcofw_11", c_double, "Funcofw(11)")
        ]
    
    _fortran_class_module_ = 'classes'
    _fortran_class_name_ = 'TModGravityModel'