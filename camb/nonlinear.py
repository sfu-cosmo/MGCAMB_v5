from .baseconfig import F2003Class, fortran_class
from ctypes import c_int, c_double

#ZW
from ctypes import byref
import numpy as np
import pyreact
import copy

Transfer_tot = c_int(7)
Transfer_nonu = c_int(8)

class NonLinearModel(F2003Class):
    """
    Abstract base class for non-linear correction models
    """
    _fields_ = [("Min_kh_nonlinear", c_double, "minimum k/h at which to apply non-linear corrections")]





halofit_original = 'original'
halofit_bird = 'bird'
halofit_peacock = 'peacock'
halofit_takahashi = 'takahashi'
halofit_mead = 'mead'
halofit_halomodel = 'halomodel'
halofit_casarini = 'casarini'
halofit_mead2015 = 'mead2015'
halofit_mead2016 = 'mead2016'
halofit_mead2020 = 'mead2020'
halofit_mead2020_feedback = 'mead2020_feedback'

halofit_default = halofit_mead2020

halofit_version_names = {halofit_original: 1,
                         halofit_bird: 2,
                         halofit_peacock: 3,
                         halofit_takahashi: 4,
                         halofit_mead: 5,
                         halofit_halomodel: 6,
                         halofit_casarini: 7,
                         halofit_mead2015: 8,
                         halofit_mead2016: 5,
                         halofit_mead2020: 9,
                         halofit_mead2020_feedback: 10}




@fortran_class
class Halofit(NonLinearModel):
    """
    Various specific approximate non-linear correction models based on HaloFit.
    """
    _fields_ = [
        ("halofit_version", c_int, {"names": halofit_version_names}),
        ("HMCode_A_baryon", c_double, "HMcode parameter A_baryon"),
        ("HMCode_eta_baryon", c_double, "HMcode parameter eta_baryon"),
        ("HMCode_logT_AGN", c_double, "HMcode parameter log10(T_AGN/K)"),
        #ZW
        ("p1", c_double, "p1"),
        ("p2", c_double, "p2"),
        ("p3", c_double, "p3"),
        ("p4", c_double, "p4")
    ]

    _fortran_class_module_ = 'NonLinear'
    _fortran_class_name_ = 'THalofit'

    def get_halofit_version(self):
        return self.halofit_version

    #ZW
    def get_react_function(self, CAMBdata = None, hubble_units=True, nz = 4, nk = 200, kh = None, z_lin = None,    calc_PK_lin = None):


        Params = CAMBdata.Params
        mymodel = "mgcamb_minimal"
        num_pars = 49
        extrapars = np.zeros(num_pars)     
        extrapars[0] = Params.w0DE
        extrapars[1] = Params.waDE
        extrapars[3] = self.p1  # Screening scale 
        extrapars[4] = self.p2  # Mass dependence 
        extrapars[5] = self.p3  # environment dependence 
        extrapars[6] = self.p4
        extrapars[7] = Params.GRtrans
        extrapars[8] = Params.B1
        extrapars[9] = Params.lambda1_2
        extrapars[10] = Params.ss
        extrapars[11] = Params.E11
        extrapars[12] = Params.Linder_gamma
        extrapars[13] = Params.B0
        extrapars[14] = Params.beta_star
        extrapars[15] = Params.a_star
        extrapars[16] = Params.xi_star
        extrapars[17] = Params.beta0
        extrapars[18] = Params.xi0
        extrapars[19] = Params.DilR
        extrapars[20] = Params.DilS
        extrapars[21] = Params.F_R0
        extrapars[22] = Params.FRn
        extrapars[23] = Params.mu0
        extrapars[24] = Params.ga
        extrapars[25] = Params.nn
        extrapars[26] = Params.MGCAMB_Mu_idx_1
        extrapars[27] = Params.MGCAMB_Mu_idx_2
        extrapars[28] = Params.MGCAMB_Mu_idx_3
        extrapars[29] = Params.MGCAMB_Mu_idx_4
        extrapars[30] = Params.MGCAMB_Mu_idx_5
        extrapars[31] = Params.MGCAMB_Mu_idx_6
        extrapars[32] = Params.MGCAMB_Mu_idx_7
        extrapars[33] = Params.MGCAMB_Mu_idx_8
        extrapars[34] = Params.MGCAMB_Mu_idx_9
        extrapars[35] = Params.MGCAMB_Mu_idx_10
        extrapars[36] = Params.MGCAMB_Mu_idx_11
        extrapars[37] = Params.Funcofw_1
        extrapars[38] = Params.Funcofw_2
        extrapars[39] = Params.Funcofw_3
        extrapars[40] = Params.Funcofw_4
        extrapars[41] = Params.Funcofw_5
        extrapars[42] = Params.Funcofw_6
        extrapars[43] = Params.Funcofw_7
        extrapars[44] = Params.Funcofw_8
        extrapars[45] = Params.Funcofw_9
        extrapars[46] = Params.Funcofw_10
        extrapars[47] = Params.Funcofw_11
        #pass h to react for for conversion in mu models
        h = Params.H0/100
        extrapars[48] = h

        #MG flags
        MG_flag = Params.MG_flag
        pure_MG_flag = Params.pure_MG_flag
        alt_MG_flag = Params.alt_MG_flag
        QSA_flag = Params.QSA_flag
        CDM_flag = Params.CDM_flag
        muSigma_flag = Params.muSigma_flag
        mugamma_par = Params.mugamma_par
        musigma_par = Params.musigma_par
        DE_model = Params.DE_model

        flagarr = [MG_flag, pure_MG_flag, alt_MG_flag, QSA_flag, CDM_flag, muSigma_flag, mugamma_par, musigma_par, DE_model]

        Omega_m = Params.omegam
        Omega_b = Params.ombh2/h**2
        Omega_nu = Params.omnuh2/h**2
        n_s = Params.InitPower.ns
        A_s = Params.InitPower.As


        PK_lin_tot, PK_lin_cb, PK_lcdm_cb = self.get_power_spectrum_asinput(CAMBdata=CAMBdata,
                                         hubble_units=hubble_units, nz=nz, nk=nk, calc_PK_lin=calc_PK_lin)

        #Now run ReACT to get the reaction and the modified gravity linear power spectrum
        react = pyreact.ReACT()
        z_lin = np.array(z_lin)

        z_react = z_lin[z_lin < 2.5]
        z_lin_ind = len(z_react)
        PK_lin_tot_cut = PK_lin_tot[:z_lin_ind, :]
        PK_lin_cb_cut = PK_lin_cb[:z_lin_ind, :]
        PK_lcdm_cb_cut = PK_lcdm_cb[:z_lin_ind, :]

        massloop = 30

        #setting wrappers
        react.setMGflags_wrapper(flagarr)
        react.setMGparams_wrapper(extrapars[7:len(extrapars)])

        #make sure to call reconstruction ahead
        react.get_reconstruction_arr(Omega_m)

        react, pofk_lin_MG_react, sigma_8, pseudo = react.compute_reaction_nu_ext(
                                        h, n_s, Omega_m, Omega_b, Omega_nu, A_s, 
                                        z_react, kh, PK_lin_tot_cut.flatten(), PK_lin_cb_cut.flatten(),
                                        kh, PK_lcdm_cb_cut.flatten(), 
                                        pscale = 0.05,
                                        model=mymodel, 
                                        extpars=extrapars,
                                        compute_pseudo=True,
                                        is_transfer=False, mass_loop=massloop,
                                        verbose=False)
        return react, pofk_lin_MG_react, z_react, pseudo


    def get_power_spectrum_asinput(self, CAMBdata, hubble_units=True, nz = 4, 
                           nk = 200, calc_PK_lin = None):
        r"""
        get different types of matter power spectrums as the input for react

        """

        PK_lin_tot_in = np.empty((nz, nk))
        PK_lin_cb_in = np.empty((nz, nk))
        PK_lcdm_cb_in = np.empty((nz, nk))


        #default MG 
        calc_PK_lin(byref(CAMBdata), PK_lin_tot_in, byref(Transfer_tot), byref(Transfer_tot), byref(hubble_units))
        calc_PK_lin(byref(CAMBdata), PK_lin_cb_in, byref(Transfer_nonu), byref(Transfer_nonu), byref(hubble_units))

        #shallow copy so int MG_flag could be changed without affecting the orginal value
        new_CAMBdata = copy.copy(CAMBdata)

        #lcdm case
        new_CAMBdata.Params.set_mgparams(MG_flag=0)

        new_CAMBdata.calc_power_spectra(new_CAMBdata.Params)

        calc_PK_lin(byref(new_CAMBdata), PK_lcdm_cb_in, byref(Transfer_nonu), byref(Transfer_nonu), byref(hubble_units))

        return PK_lin_tot_in, PK_lin_cb_in, PK_lcdm_cb_in

    #ZW
    def set_params(self, halofit_version=halofit_default, HMCode_A_baryon=3.13, HMCode_eta_baryon=0.603,
                   HMCode_logT_AGN=7.8, p1 = 0.5, p2 = 0.5 ,p3 = 0.5, p4 = 0.):
        """
        Set the halofit model for non-linear corrections.

        :param halofit_version: One of

            - original: `astro-ph/0207664 <https://arxiv.org/abs/astro-ph/0207664>`_
            - bird: `arXiv:1109.4416 <https://arxiv.org/abs/1109.4416>`_
            - peacock: `Peacock fit <http://www.roe.ac.uk/~jap/haloes/>`_
            - takahashi: `arXiv:1208.2701 <https://arxiv.org/abs/1208.2701>`_
            - mead: HMCode `arXiv:1602.02154 <https://arxiv.org/abs/1602.02154>`_
            - halomodel: basic halomodel
            - casarini: PKequal `arXiv:0810.0190 <https://arxiv.org/abs/0810.0190>`_, `arXiv:1601.07230 <https://arxiv.org/abs/1601.07230>`_
            - mead2015: original 2015 version of HMCode `arXiv:1505.07833 <https://arxiv.org/abs/1505.07833>`_
            - mead2016: Alias for 'mead'.
            - mead2020: 2020 version of HMcode `arXiv:2009.01858 <https://arxiv.org/abs/2009.01858>`_
            - mead2020_feedback: 2020 version of HMcode with baryonic feedback `arXiv:2009.01858 <https://arxiv.org/abs/2009.01858>`_
        :param HMCode_A_baryon: HMcode parameter A_baryon. Default 3.13. Used only in models mead2015 and mead2016 (and its alias mead).
        :param HMCode_eta_baryon: HMcode parameter eta_baryon. Default 0.603. Used only in mead2015 and mead2016 (and its alias mead).
        :param HMCode_logT_AGN: HMcode parameter logT_AGN. Default 7.8. Used only in model mead2020_feedback.
        """
        self.halofit_version = halofit_version
        self.HMCode_A_baryon = HMCode_A_baryon
        self.HMCode_eta_baryon = HMCode_eta_baryon
        self.HMCode_logT_AGN = HMCode_logT_AGN

        # nonlinear parameters in React
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4


@fortran_class
class SecondOrderPK(NonLinearModel):
    """
    Third-order Newtonian perturbation theory results for the non-linear correction.
    Only intended for use at very high redshift (z>10) where corrections are perturbative, it will not give
    sensible results at low redshift.

    See Appendix F of `astro-ph/0702600 <https://arxiv.org/abs/astro-ph/0702600>`_ for equations and references.

    Not intended for production use, it's mainly to serve as an example alternative non-linear model implementation.
    """

    _fortran_class_module_ = 'SecondOrderPK'
    _fortran_class_name_ = 'TSecondOrderPK'

    def set_params(self):
        pass

