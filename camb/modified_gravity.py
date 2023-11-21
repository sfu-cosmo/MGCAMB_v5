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


    def set_params(self, MG_wrapped=True, MG_flag=0, GRtrans=.001, 
                   pure_MG_flag=1, alt_MG_flag=1, QSA_flag=1, CDM_flag=1, 
                   muSigma_flag=1, mugamma_par=1, B1=1.333, lambda1_2=1000,
                   B2=.5, lambda2_2=1000, ss=4, E11=1., E22=1.,
				   ga=.5, nn=2, musigma_par=1, mu0=0., sigma0=0, 
                   QR_par=1, MGQfix=1, MGRfix=1, Qnot=1., Rnot=1.,
                   sss=0, Linder_gamma=.545, B0=.001, beta_star=1.,
                   a_star=.5, xi_star=.001, beta0=0., xi0=.0001,
                   DilS=.24, DilR=1., F_R0=.0001, FRn=1., 
				   DE_model=0, w0DE=-1., waDE=0.0, MGDE_pert=False,
                   MGCAMB_Mu_idx_1=1., MGCAMB_Mu_idx_2=1., 
				   MGCAMB_Mu_idx_3=1., MGCAMB_Mu_idx_4=1.,
                   MGCAMB_Mu_idx_5=1., MGCAMB_Mu_idx_6=1.,
                   MGCAMB_Mu_idx_7=1., MGCAMB_Mu_idx_8 = 1.,
                   MGCAMB_Mu_idx_9=1., MGCAMB_Mu_idx_10=1.,
                   MGCAMB_Mu_idx_11=1.,MGCAMB_Sigma_idx_1=1.,
				   MGCAMB_Sigma_idx_2=1., MGCAMB_Sigma_idx_3=1.,
                   MGCAMB_Sigma_idx_4=1., MGCAMB_Sigma_idx_5=1.,
				   MGCAMB_Sigma_idx_6=1., MGCAMB_Sigma_idx_7=1.,
                   MGCAMB_Sigma_idx_8=1., MGCAMB_Sigma_idx_9=1.,
				   MGCAMB_Sigma_idx_10=1.,MGCAMB_Sigma_idx_11=1.,
                   Funcofw_1=.7, Funcofw_2=.7, Funcofw_3=.7, Funcofw_4=.7,
				   Funcofw_5=.7,Funcofw_6=.7,Funcofw_7=.7,
                   Funcofw_8=.7,Funcofw_9=.7,Funcofw_10=.7,Funcofw_11=.7):

        self.MG_wrapped = MG_wrapped
        self.MG_flag = MG_flag
        self.pure_MG_flag = pure_MG_flag
        self.alt_MG_flag = alt_MG_flag
        self.QSA_flag = QSA_flag
        self.muSigma_flag = muSigma_flag
        self.CDM_flag  = CDM_flag
        self.mugamma_par  = mugamma_par
        self.musigma_par  = musigma_par
        self.QR_par = QR_par
        self.DE_model = DE_model
        self.GRtrans = GRtrans                 
        self.B1 =  B1
        self.B2 =  B2
        self.lambda1_2 =  lambda1_2
        self.lambda2_2 = lambda2_2
        self.ss = ss
        self.B0 = B0
        self.E11 = E11
        self.E22 = E22
        self.MGQfix = MGQfix
        self.MGRfix = MGRfix
        self.Qnot = Qnot
        self.Rnot = Rnot
        self.sss = sss
        self.Linder_gamma = Linder_gamma    
        self.xi_star = xi_star 
        self.beta_star = beta_star
        self.a_star  = a_star
        self.beta0 = beta0
        self.xi0 = xi0
        self.DilR = DilR
        self.DilS = DilS
        self.F_R0 = F_R0
        self.FRn = FRn
        self.mu0 = mu0
        self.sigma0 = sigma0
        self.ga = ga
        self.nn = nn
        self.w0DE = w0DE    
        self.waDE = waDE      
        self.MGDE_pert = MGDE_pert
        self.MGCAMB_Mu_idx_1 =  MGCAMB_Mu_idx_1
        self.MGCAMB_Mu_idx_2 =  MGCAMB_Mu_idx_2
        self.MGCAMB_Mu_idx_3 =  MGCAMB_Mu_idx_3
        self.MGCAMB_Mu_idx_4 =  MGCAMB_Mu_idx_4
        self.MGCAMB_Mu_idx_5 =  MGCAMB_Mu_idx_5
        self.MGCAMB_Mu_idx_6 =  MGCAMB_Mu_idx_6
        self.MGCAMB_Mu_idx_7 =  MGCAMB_Mu_idx_7
        self.MGCAMB_Mu_idx_8 =  MGCAMB_Mu_idx_8
        self.MGCAMB_Mu_idx_9 =  MGCAMB_Mu_idx_9
        self.MGCAMB_Mu_idx_10 =  MGCAMB_Mu_idx_10
        self.MGCAMB_Mu_idx_11 =  MGCAMB_Mu_idx_11
        self.MGCAMB_Sigma_idx_1 =  MGCAMB_Sigma_idx_1
        self.MGCAMB_Sigma_idx_2 =  MGCAMB_Sigma_idx_2
        self.MGCAMB_Sigma_idx_3 =  MGCAMB_Sigma_idx_3
        self.MGCAMB_Sigma_idx_4 =  MGCAMB_Sigma_idx_4
        self.MGCAMB_Sigma_idx_5 =  MGCAMB_Sigma_idx_5
        self.MGCAMB_Sigma_idx_6 =  MGCAMB_Sigma_idx_6
        self.MGCAMB_Sigma_idx_7 =  MGCAMB_Sigma_idx_7
        self.MGCAMB_Sigma_idx_8 =  MGCAMB_Sigma_idx_8
        self.MGCAMB_Sigma_idx_9 =  MGCAMB_Sigma_idx_9
        self.MGCAMB_Sigma_idx_10 =  MGCAMB_Sigma_idx_10
        self.MGCAMB_Sigma_idx_11 =  MGCAMB_Sigma_idx_11
        self.Funcofw_1 =  Funcofw_1
        self.Funcofw_2 =  Funcofw_2
        self.Funcofw_3 =  Funcofw_3
        self.Funcofw_4 =  Funcofw_4
        self.Funcofw_5 =  Funcofw_5
        self.Funcofw_6 =  Funcofw_6
        self.Funcofw_7 =  Funcofw_7
        self.Funcofw_8 =  Funcofw_8
        self.Funcofw_9 =  Funcofw_9
        self.Funcofw_10 =  Funcofw_10
        self.Funcofw_11 =  Funcofw_11