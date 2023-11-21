    module classes
    use precision
    use MpiUtils
    use interpolation, only : TSpline1D, TCubicSpline, TLogRegularCubicSpline
    use IniObjects
    implicit none

    Type MatterTransferData
        !Computed data
        integer   ::  num_q_trans   !    number of steps in k for transfer calculation
        real(dl), dimension (:), allocatable :: q_trans
        real(dl), dimension (:), allocatable ::  sigma_8
        real(dl), dimension (:), allocatable ::  sigma2_vdelta_8 !growth from sigma_{v delta}
        real(dl), dimension(:,:,:), allocatable :: TransferData
        !TransferData(entry,k_index,z_index) for entry=Tranfer_kh.. Transfer_tot
    contains
    procedure :: Free => MatterTransferData_Free
    end Type MatterTransferData

    Type MatterPowerData
        !everything is a function of k/h
        integer   ::  num_k, num_z
        real(dl), dimension(:), allocatable :: log_kh, redshifts
        !matpower is log(P_k)
        real(dl), dimension(:,:), allocatable :: matpower, ddmat
        !if NonLinear, nonlin_ratio =  sqrt(P_nonlinear/P_linear)
        !function of k and redshift NonLinearScaling(k_index,z_index)
        real(dl), dimension(:,:), allocatable :: nonlin_ratio
        !Sources
        real(dl), dimension(:), allocatable :: log_k
        real(dl), dimension(:,:), allocatable :: vvpower, ddvvpower
        real(dl), dimension(:,:), allocatable :: vdpower, ddvdpower

        real(dl), dimension(:,:), allocatable :: nonlin_ratio_vv
        real(dl), dimension(:,:), allocatable :: nonlin_ratio_vd

    end Type MatterPowerData

    !Classes that can be accessed from python and contain allocatable objects and arrays or other classes
    !In Python inherited from F2003Class (defined in baseconfig)
    !All python-accessible inherited classes must define SelfPointer, and use @fortran_class decorator in python
    !Allocatable objects can only contain instances of python-accessible classes so they can be identified from python.
    !Class could be abstract and SelfPointer deferred, but Fortran then doesn't allow called inherited "super()"
    !methods implemented in abstract classes
    Type TPythonInterfacedClass
    contains
    !SelfPointer msust be overriden for each class to be referenced from python.
    !Gets a class pointer from an untyped pointer.
    procedure, nopass :: SelfPointer
    !PythonClass gets string of python class name; not actually used internally
    procedure, nopass :: PythonClass
    !Replace is for copying over this type with an instance of the same type
    !It must be defined for any class used as a non-allocatable python structure field that can be assigned to
    procedure :: Replace
    end type TPythonInterfacedClass

    Type PythonClassPtr
        class(TPythonInterfacedClass), pointer :: Ref => null()
    end type

    Type PythonClassAllocatable
        class(TPythonInterfacedClass), allocatable :: P
    end Type PythonClassAllocatable


    type, extends(TPythonInterfacedClass) :: TCambComponent
    contains
    procedure :: ReadParams => TCambComponent_ReadParams
    procedure :: Validate =>  TCambComponent_Validate
    end type TCambComponent


    type, extends(TPythonInterfacedClass) :: TCAMBParameters
        !Actual type defined in model.f90
    end type TCAMBParameters

    Type, extends(TPythonInterfacedClass) :: TCAMBdata
        !Actual type defined in results.f90
    end type TCAMBdata

    type, extends(TCambComponent) :: TNonLinearModel
        real(dl) :: Min_kh_nonlinear  = 0.005_dl
    contains
    procedure :: Init => TNonLinearModel_init
    procedure :: GetNonLinRatios => TNonLinearModel_GetNonLinRatios
    procedure :: GetNonLinRatios_All => TNonLinearModel_GetNonLinRatios_All
    end type TNonLinearModel

    !> MGCAMB MOD START
    type, extends(TCambComponent) :: TModGravModel

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

        procedure :: Init => TModGravModel_Init
        !procedure :: ReadParams => MGCAMB_read_in_MGparams
        procedure, nopass :: PythonClass => TModGravModel_PythonClass
        procedure, nopass :: SelfPointer => TModGravModel_SelfPointer
        !procedure :: Validate => TModGravModel_Validate

    end type TModGravModel
    !< MGCAMB MOD END

    type, extends(TCambComponent) :: TInitialPower
    contains
    procedure :: ScalarPower => TInitialPower_ScalarPower
    procedure :: TensorPower => TInitialPower_TensorPower
    procedure :: Init => TInitialPower_Init
    procedure :: Effective_ns => TInitalPower_Effective_ns
    end type TInitialPower

    Type, extends(TCambComponent) :: TRecombinationModel
        real(dl) :: min_a_evolve_Tm = 1/(1+900.) !scale factor at which to start evolving Delta_TM
    contains
    procedure :: Init => TRecombinationModel_init
    procedure :: x_e => TRecombinationModel_xe !ionization fraction
    procedure :: xe_Tm => TRecombinationModel_xe_Tm !ionization fraction and baryon temperature
    procedure :: T_m => TRecombinationModel_tm !baryon temperature
    procedure :: T_s => TRecombinationModel_ts !Spin temperature
    procedure :: Version => TRecombinationModel_version
    procedure :: dDeltaxe_dtau => TRecombinationModel_dDeltaxe_dtau
    procedure :: get_Saha_z => TRecombinationModel_Get_Saha_z
    end type

    Type, extends(TCambComponent) :: TReionizationModel
        logical  :: Reionization = .true.
    contains
    procedure :: Init => TReionizationModel_Init
    procedure :: x_e => TReionizationModel_xe
    procedure :: get_timesteps => TReionizationModel_get_timesteps
    end Type TReionizationModel


    interface
    subroutine TClassDverk (this,n, fcn, x, y, xend, tol, ind, c, nw, w)
    use Precision
    import
    class(TCambComponent), target :: this
    integer n, ind
    real(dl) x, y(n), xend, tol, c(*), w(nw,9)
    external fcn
    end subroutine
    end interface

    contains

    subroutine TCambComponent_ReadParams(this, Ini)
    class(TCambComponent) :: this
    class(TIniFile), intent(in) :: Ini

    end subroutine TCambComponent_ReadParams

    subroutine  TCambComponent_Validate(this, OK)
    class(TCambComponent),intent(in) :: this
    logical, intent(inout) :: OK

    end subroutine TCambComponent_Validate


    function PythonClass()
    character(LEN=:), allocatable :: PythonClass

    PythonClass = ''
    error stop 'PythonClass Not implemented'

    end function PythonClass

    subroutine SelfPointer(cptr, P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TPythonInterfacedClass), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P

    call c_f_pointer(cptr, PType)
    P => PType
    error stop 'Class must define SelfPointer function returning pointer to actual type'

    end subroutine SelfPointer

    subroutine Replace(this, replace_with)
    class(TPythonInterfacedClass), target :: this
    class(TPythonInterfacedClass) :: replace_with

    error stop 'Assignment not implemented for this class'

    end subroutine Replace

    subroutine TNonLinearModel_Init(this, State)
    class(TNonLinearModel) :: this
    class(TCAMBdata), target :: State
    end subroutine TNonLinearModel_Init

    subroutine TNonLinearModel_GetNonLinRatios(this,State,CAMB_Pk)
    class(TNonLinearModel) :: this
    class(TCAMBdata) :: State
    type(MatterPowerData), target :: CAMB_Pk
    error stop 'GetNonLinRatios Not implemented'
    end subroutine TNonLinearModel_GetNonLinRatios

    subroutine TNonLinearModel_GetNonLinRatios_All(this,State,CAMB_Pk)
    class(TNonLinearModel) :: this
    class(TCAMBdata) :: State
    type(MatterPowerData), target :: CAMB_Pk
    error stop 'GetNonLinRatios_all  not supported (no non-linear velocities)'
    end subroutine TNonLinearModel_GetNonLinRatios_All

    !> MGCAMB MOD START
    subroutine TModGravModel_Init(this, Params)
        class(TModGravModel) :: this
        class(TCAMBParameters), intent(in) :: Params
    end subroutine TModGravModel_Init

    function TModGravModel_PythonClass()
        character(LEN=:), allocatable :: TModGravModel_PythonClass
        TModGravModel_PythonClass = 'ModGravModel'
    end function TModGravModel_PythonClass

    subroutine TModGravModel_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TModGravModel), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P
    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TModGravModel_SelfPointer
    !> MGCAMB MOD END

    function TInitialPower_ScalarPower(this, k)
    class(TInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TInitialPower_ScalarPower

    TInitialPower_ScalarPower = 0
    error stop 'ScalarPower not implemented'
    end function TInitialPower_ScalarPower

    function TInitialPower_TensorPower(this, k)
    class(TInitialPower) :: this
    real(dl), intent(in) ::k
    real(dl) TInitialPower_TensorPower

    TInitialPower_TensorPower = 0
    error stop 'TensorPower not implemented'
    end function TInitialPower_TensorPower


    subroutine TInitialPower_Init(this, Params)
    class(TInitialPower) :: this
    class(TCAMBParameters), intent(in) :: Params
    end subroutine TInitialPower_Init

    function TInitalPower_Effective_ns(this)
    class(TInitialPower) :: this
    real(dl) :: TInitalPower_Effective_ns

    TInitalPower_Effective_ns=1
    call MpiStop('Power spectrum does not support an effective n_s for halofit')

    end function TInitalPower_Effective_ns


    subroutine MatterTransferData_Free(this)
    class(MatterTransferData):: this
    integer st

    deallocate(this%q_trans, STAT = st)
    deallocate(this%TransferData, STAT = st)
    deallocate(this%sigma_8, STAT = st)
    deallocate(this%sigma2_vdelta_8, STAT = st)

    end subroutine MatterTransferData_Free


    function TRecombinationModel_tm(this,a)
    class(TRecombinationModel) :: this
    real(dl), intent(in) :: a
    real(dl) TRecombinationModel_tm

    call MpiStop('TRecombinationModel_tm not implemented')
    TRecombinationModel_tm=0

    end function TRecombinationModel_tm


    function TRecombinationModel_ts(this,a)
    class(TRecombinationModel) :: this
    !zrec(1) is zinitial-delta_z
    real(dl), intent(in) :: a
    real(dl) TRecombinationModel_ts

    call MpiStop('TRecombinationModel_ts not implemented')
    TRecombinationModel_ts=0

    end function TRecombinationModel_ts

    function TRecombinationModel_xe(this,a)
    class(TRecombinationModel) :: this
    real(dl), intent(in) :: a
    real(dl) TRecombinationModel_xe

    call MpiStop('TRecombinationModel_xe not implemented')
    TRecombinationModel_xe=0

    end function TRecombinationModel_xe

    subroutine TRecombinationModel_xe_Tm(this,a, xe, Tm)
    !Not required to implement, but may be able to optimize
    class(TRecombinationModel) :: this
    real(dl), intent(in) :: a
    real(dl), intent(out) :: xe, Tm

    xe = this%x_e(a)
    Tm = this%T_m(a)

    end subroutine TRecombinationModel_xe_Tm

    function TRecombinationModel_version(this) result(version)
    class(TRecombinationModel) :: this
    character(LEN=:), allocatable :: version

    version = ''

    end function TRecombinationModel_version

    subroutine TRecombinationModel_init(this,State, WantTSpin)
    class(TRecombinationModel), target :: this
    class(TCAMBdata), target :: State
    logical, intent(in), optional :: WantTSpin
    end subroutine TRecombinationModel_init

    function TRecombinationModel_dDeltaxe_dtau(this,a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb,adotoa)
    !d x_e/d tau
    class(TRecombinationModel) :: this
    real(dl) TRecombinationModel_dDeltaxe_dtau
    real(dl), intent(in):: a, Delta_xe,Delta_nH, Delta_Tm, hdot, kvb,adotoa

    call MpiStop('TRecombinationModel_dDeltaxe_dtau not implemented')
    TRecombinationModel_dDeltaxe_dtau=0

    end function TRecombinationModel_dDeltaxe_dtau

    real(dl) function TRecombinationModel_Get_Saha_z(this)
    class(TRecombinationModel) :: this
    call MpiStop('TRecombinationModel_Get_Saha_z not implemented')
    TRecombinationModel_Get_Saha_z = 0
    end function

    function TReionizationModel_xe(this, z, tau, xe_recomb)
    !a and time tau and redundant, both provided for convenience
    !xe_recomb is xe(tau_start) from recombination (typically very small, ~2e-4)
    !xe should map smoothly onto xe_recomb
    class(TReionizationModel) :: this
    real(dl), intent(in) :: z
    real(dl), intent(in), optional :: tau, xe_recomb
    real(dl) TReionizationModel_xe

    call MpiStop('TReionizationModel_xe not implemented')
    TReionizationModel_xe=0

    end function TReionizationModel_xe

    subroutine TReionizationModel_get_timesteps(this, n_steps, z_start, z_complete)
    !minimum number of time steps to use between tau_start and tau_complete
    !Scaled by AccuracyBoost later
    !steps may be set smaller than this anyway
    class(TReionizationModel) :: this
    integer, intent(out) :: n_steps
    real(dl), intent(out):: z_start, z_Complete
    call MpiStop('TReionizationModel_get_timesteps not implemented')
    n_steps=0
    z_start=0
    z_complete=0
    end subroutine TReionizationModel_get_timesteps

    subroutine TReionizationModel_Init(this, State)
    class(TReionizatioNModel) :: this
    class(TCAMBdata), target :: State
    end subroutine TReionizationModel_Init

    end module classes
