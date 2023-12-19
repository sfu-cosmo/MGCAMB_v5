module ModGravityInterface
    use classes
    use precision
    use IniObjects
    implicit none

    ! Base type for both dark energy and modified gravity


    type, extends(TDarkSector) :: TModGravityModel

        character(len=30) :: debug_root = "debug_"
        real(dl):: GRtrans = 0.001

        ! Whether or not we want additional DE perturbations on top of what MGCAMB already
        ! provides in terms of the phenomenological functions
        logical :: MGDE_pert = .False.

        ! 1. Background quantities
        real(dl) :: adotoa
        real(dl) :: Hdot
        real(dl) :: grho
        real(dl) :: gpres 
        real(dl) :: grhob_t    
        real(dl) :: grhoc_t
        real(dl) :: grhog_t 
        real(dl) :: grhor_t
        real(dl) :: grhov_t  
        real(dl) :: gpresv_t 
        real(dl) :: grhonu_t
        real(dl) :: gpresnu_t

        ! 2. Perturbation quantities
        real(dl) :: k
        real(dl) :: k2
        real(dl) :: dgrho
        real(dl) :: dgrhoc
        real(dl) :: dgq
        real(dl) :: dgqc 
        real(dl) :: pidot_sum
        real(dl) :: dgpi_w_sum
        real(dl) :: dgpi   
        real(dl) :: dgpi_diff
        real(dl) :: dgpidot
        real(dl) :: rhoDelta  
        real(dl) :: rhoDeltadot  
        real(dl) :: rhoDeltac
        real(dl) :: rhoDeltacdot 

        ! 3. MG functions
        real(dl) :: mu = 0.
        real(dl) :: mudot = 0.
        real(dl) :: gamma = 0.
        real(dl) :: gammadot = 0.
        real(dl) :: q = 0.
        real(dl) :: qdot = 0.
        real(dl) :: r = 0.
        real(dl) :: rdot = 0.
        real(dl) :: BigSigma = 0.
        real(dl) :: BigSigmadot = 0.
        real(dl) :: C_phi = 0.
        real(dl) :: C_phidot = 0.

        !> 4. Perturbations evolution variables
        real(dl) :: z
        real(dl) :: sigma
        real(dl) :: sigmadot
        real(dl) :: etak
        real(dl) :: etadot

        !> 5. ISW and lensing realted quantities
        real(dl) :: MG_alpha
        real(dl) :: MG_alphadot
        real(dl) :: MG_phi
        real(dl) :: MG_phidot
        real(dl) :: MG_psi
        real(dl) :: MG_psidot
        real(dl) :: MG_ISW
        real(dl) :: MG_lensing
        real(dl) :: source1
        real(dl) :: source3

        contains

        procedure :: Init => TModGravityModel_Init
        procedure, nopass :: PythonClass => TModGravityModel_PythonClass
        procedure, nopass :: SelfPointer => TModGravityModel_SelfPointer
        procedure :: ResetCache => TModGravityModel_ResetCache
        procedure :: MGCAMB_open_cache_files
        procedure :: MGCAMB_close_cache_files
        procedure :: MGCAMB_dump_cache

    end type TModGravityModel


    contains

    subroutine TModGravityModel_Init(this, State)
        class(TModGravityModel) :: this
        class(TCAMBdata), intent(in), target :: State
    end subroutine TModGravityModel_Init

    function TModGravityModel_PythonClass()
        character(LEN=:), allocatable :: TModGravityModel_PythonClass
        TModGravityModel_PythonClass = 'ModGravityModel'
    end function TModGravityModel_PythonClass

    subroutine TModGravityModel_SelfPointer(cptr,P)
    use iso_c_binding
    Type(c_ptr) :: cptr
    Type (TModGravityModel), pointer :: PType
    class (TPythonInterfacedClass), pointer :: P
    call c_f_pointer(cptr, PType)
    P => PType
    end subroutine TModGravityModel_SelfPointer

    ! TODO: this should be put somewhere more sensible, for now it's here
    !> Subroutine that sets the mgcamb_cache to zero
    subroutine TModGravityModel_ResetCache( this )

        class(TModGravityModel) :: this

        ! 1. Background quantities
        this%adotoa     = 0._dl
        this%Hdot       = 0._dl
        this%grho       = 0._dl
        this%gpres      = 0._dl
        this%grhob_t    = 0._dl
        this%grhoc_t    = 0._dl
        this%grhog_t    = 0._dl
        this%grhor_t    = 0._dl
        this%grhov_t    = 0._dl
        this%gpresv_t   = 0._dl
        this%grhonu_t   = 0._dl
        this%gpresnu_t  = 0._dl

        ! 2. Perturbation quantities
        this%k          = 0._dl
        this%k2         = 0._dl
        this%dgrho      = 0._dl
        this%dgq        = 0._dl
        this%pidot_sum  = 0._dl
        this%dgpi_w_sum = 0._dl
        this%dgpi       = 0._dl
        this%dgpi_diff  = 0._dl
        this%dgpidot    = 0._dl
        this%rhoDelta   = 0._dl
        this%rhoDeltadot= 0._dl

        ! 3. MG functions
        this%mu         = 0._dl
        this%mudot      = 0._dl
        this%gamma      = 0._dl
        this%gammadot   = 0._dl
        this%q          = 0._dl
        this%qdot       = 0._dl
        this%r          = 0._dl
        this%rdot       = 0._dl

        !> 4. Perturbations evolution variables
        this%z          = 0._dl
        this%sigma      = 0._dl
        this%sigmadot   = 0._dl
        this%etak       = 0._dl
        this%etadot     = 0._dl

        !> 5. ISW and lensing realted quantities
        this%MG_alpha   = 0._dl
        this%MG_alphadot= 0._dl
        this%MG_phi     = 0._dl
        this%MG_phidot  = 0._dl
        this%MG_psi     = 0._dl
        this%MG_psidot  = 0._dl
        this%MG_ISW     = 0._dl
        this%MG_lensing = 0._dl
        this%source1    = 0._dl
        this%source3    = 0._dl

    end subroutine

        ! ---------------------------------------------------------------------------------------------
    !> Subroutine that opens the MGCAMB cache files (for Debug)
    subroutine MGCAMB_open_cache_files( this )

        class(TModGravityModel) :: this

        ! 1. Open sources file
        open(unit=111, file=TRIM(this%debug_root) // 'MGCAMB_debug_sources.dat', status="new", &
            & action="write")
        write(111,*)  'k  ', 'a  ', 'MG_ISW  ', 'MG_Lensing  ', 'S_T  ', 'S_lensing'

        ! 2 Open MG functions file
        open(unit=222, file=TRIM(this%debug_root) // 'MGCAMB_debug_MG_fncs.dat', status="new",&
            & action="write")
        write(222,*)  'k  ', 'a  ', 'mu  ', 'gamma ', 'Q ', 'R ', 'Phi ', 'Psi ', 'dPhi ', 'dPsi '

        ! 3. Open Einstein solutions file
        open(unit=333, file=TRIM(this%debug_root) // 'MGCAMB_debug_EinsteinSol.dat', status="new",&
            & action="write")
        write(333,*) 'k  ', 'a  ', 'etak  ', 'z  ', 'sigma  ', 'etadot  ', 'sigmadot  '

        ! 4. Open Perturbation solution file
        open(unit=444, file=TRIM(this%debug_root) // 'MGCAMB_debug_PerturbSol.dat', status="new",&
        & action="write")
        write(444,*)  'k  ', 'a  ', 'dgrho  ', 'dgq  ', 'rhoDelta  ', 'dgpi  ', 'pidot_sum  ', 'dgpi_w_sum  '

        ! 5. Open Background file
        open(unit=555, file=TRIM(this%debug_root) // 'MGCAMB_debug_Background.dat', status="new",&
            & action="write")
        write(555,*)  'k  ', 'a  ', 'H  ', 'Hdot  ', 'grhov_t  ', 'gpresv_t  '

    end subroutine MGCAMB_open_cache_files


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that closes the MGCAMB cache files (for Debug)
    subroutine MGCAMB_close_cache_files( this )

        class(TModGravityModel) :: this

        close(111);close(222); close(333);close(444);close(555)

    end subroutine MGCAMB_close_cache_files


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints the MGCAMB cache on a file
    subroutine MGCAMB_dump_cache( this, a )

        class(TModGravityModel) :: this

        real(dl), intent(in) :: a   !< scale factor
        character(*), parameter :: cache_output_format = 'e18.8'

        ! 1. Write the sources
        write(111,'(14'//cache_output_format//')') this%k, a, this%MG_ISW, this%MG_Lensing,&
                                                    & this%source1, this%source3

        ! 2. Write the MG functions and the potentials
        write(222,'(14'//cache_output_format//')') this%k, a, this%mu, this%gamma, this%q, this%r, &
                                                & this%MG_phi, this%MG_psi, this%MG_phidot, this%MG_psidot

        ! 3. Write the Einstein equations solutions
        write(333,'(14'//cache_output_format//')') this%k, a, this%etak, this%z, this%sigma,&
                                                & this%etadot,this%sigmadot

        ! 4. Write the Perturbations Solutions
        write(444,'(14'//cache_output_format//')') this%k, a, this%dgrho, this%dgq, this%rhoDelta,&
                                                    & this%dgpi, this%pidot_sum, this%dgpi_w_sum

        !5. Write the background
        write(555,'(14'//cache_output_format//')') this%k, a, this%adotoa, this%Hdot, this%grhov_t,&
                                                    & this%gpresv_t

    end subroutine MGCAMB_dump_cache


end module ModGravityInterface