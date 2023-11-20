module ModGravityInterface
    use classes
    implicit none
    public

    type, extends(TModGravModel) :: TBaseModGravModel

    logical :: MG_wrapped = .False.
    integer :: MG_flag = 0
    real(dl):: GRtrans = 0.001
    integer :: pure_MG_flag = 1
    integer :: alt_MG_flag = 1
    integer :: QSA_flag = 1
    integer :: CDM_flag = 1
    integer :: muSigma_flag = 1
    integer :: mugamma_par = 1
    real(dl):: B1 = 1.333d0
    real(dl):: lambda1_2  = 1000
    real(dl):: B2 = 0.5
    real(dl):: lambda2_2 = 1000
    real(dl):: ss = 4
    real(dl):: E11 = 1.0d0
    real(dl):: E22 = 1.0d0
    real(dl):: ga = 0.5
    real(dl):: nn = 2
    integer :: musigma_par = 1
    real(dl):: mu0 = 0.d0
    real(dl):: sigma0 = 0.d0
    integer :: QR_par = 1
    real(dl):: MGQfix = 1
    real(dl):: MGRfix = 1
    real(dl):: Qnot = 1.d0
    real(dl):: Rnot= 1.d0
    real(dl):: sss = 0
    real(dl):: Linder_gamma = 0.545
    real(dl):: B0 = 1.d-3
    real(dl):: beta_star  = 1.0
    real(dl):: a_star = 0.5 
    real(dl):: xi_star = 1.d-3
    real(dl):: beta0 = 1.d0
    real(dl):: xi0 = 1.d-4
    real(dl):: DilS = 0.24d0
    real(dl):: DilR = 1.d0
    real(dl):: F_R0 = 1.d-4
    real(dl):: FRn = 1.d0
    integer :: DE_model  = 0
    real(dl):: w0DE = -1.d0
    real(dl):: waDE = 0.d0
    logical :: MGDE_pert = .False.
    real(dl):: MGCAMB_Mu_idx_1 = 1.d0
    real(dl):: MGCAMB_Mu_idx_2 = 1.0d0
    real(dl):: MGCAMB_Mu_idx_3 = 1.0d0
    real(dl):: MGCAMB_Mu_idx_4 = 1.0d0
    real(dl):: MGCAMB_Mu_idx_5 = 1.0d0
    real(dl):: MGCAMB_Mu_idx_6 = 1.0d0
    real(dl):: MGCAMB_Mu_idx_7 = 1.0d0
    real(dl):: MGCAMB_Mu_idx_8 = 1.0d0
    real(dl):: MGCAMB_Mu_idx_9 = 1.0d0
    real(dl):: MGCAMB_Mu_idx_10 = 1.0d0
    real(dl):: MGCAMB_Mu_idx_11 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_1 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_2 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_3 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_4 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_5 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_6 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_7 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_8 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_9 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_10 = 1.0d0
    real(dl):: MGCAMB_Sigma_idx_11 = 1.0d0
    real(dl):: Funcofw_1 = 0.7d0
    real(dl):: Funcofw_2 = 0.7d0
    real(dl):: Funcofw_3 = 0.7d0
    real(dl):: Funcofw_4 = 0.7d0
    real(dl):: Funcofw_5= 0.7d0
    real(dl):: Funcofw_6= 0.7d0
    real(dl):: Funcofw_7= 0.7d0
    real(dl):: Funcofw_8= 0.7d0
    real(dl):: Funcofw_9= 0.7d0
    real(dl):: Funcofw_10= 0.7d0
    real(dl):: Funcofw_11= 0.7d0


    contains
    !procedure :: ReadParams => MGCAMB_read_in_MGparams
    procedure :: Init => TBaseModGravModel_Init
    procedure, nopass :: PythonClass => TBaseModGravModel_PythonClass
    procedure, nopass :: SelfPointer => TBaseModGravModel_SelfPointer
    !procedure :: Validate => TBaseModGravModel_Validate

    end type TBaseModGravModel

    contains

    function TBaseModGravModel_PythonClass()
        character(LEN=:), allocatable :: TBaseModGravModel_PythonClass
        TBaseModGravModel_PythonClass = 'BaseModGravModel'
    end function TBaseModGravModel_PythonClass



    subroutine TBaseModGravModel_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TBaseModGravModel), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P
    
    call c_f_pointer(cptr, PType)
    P => PType
    
    end subroutine TBaseModGravModel_SelfPointer


    subroutine TBaseModGravModel_Init(this, Params)
        use classes
        class(TBaseModGravModel) :: this
        class(TCAMBParameters), intent(in) :: Params
        
        
    end subroutine TBaseModGravModel_Init

end module ModGravityInterface