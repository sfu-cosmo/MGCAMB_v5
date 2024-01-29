module ModGravityInterface
    use classes, only : TCambComponent, TCAMBdata
    use precision
    !use IniObjects
    implicit none



    type, abstract, extends(TCambComponent) :: TModGravityModel

        ! MGCAMB debug flag.This will turn on printing of many things to aid debugging the code.
        logical :: DebugMGCAMB = .true.
        character(len=30) :: debug_root = "debug_"

        ! Choose at which time to turn on MG
        real(dl):: GRtrans = 0.001

        ! Whether or not we want additional DE perturbations on top of what MGCAMB already
        ! provides in terms of the phenomenological functions
        logical :: MGDE_pert = .false.

        ! type of modified growth model
        integer :: MG_flag = 0
        ! TODO: see if this is really needed
        !###### Part 2.4 - CDM-only coupling
        !#CDM_flag = 1: QSA models
        !#need to set other QSA models parameters as well
        integer :: CDM_flag = 0

        ! 1. Background quantities
        real(dl) :: Hdot
        !real(dl) :: grhov_t  
        !real(dl) :: gpresv_t 

        ! 2. Perturbation quantities
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
        real(dl) :: z = 0. ! this z is not the redshift
        real(dl) :: sigma = 0.
        real(dl) :: sigmadot = 0.
        real(dl) :: etadot = 0.

        !> 5. ISW and lensing realted quantities
        real(dl) :: MG_alpha = 0.
        real(dl) :: MG_alphadot = 0.
        real(dl) :: MG_phi = 0.
        real(dl) :: MG_phidot = 0.
        real(dl) :: MG_psi = 0.
        real(dl) :: MG_psidot = 0.
        real(dl) :: MG_ISW = 0.
        real(dl) :: MG_lensing = 0.
        real(dl) :: source1 = 0.
        real(dl) :: source3 = 0.

        contains

        ! subroutines that are really defined within the subtypes
        procedure(ComputeMGFunctionsInterface), deferred :: ComputeMGFunctions
        procedure(ComputesigmaInterface), deferred :: Computesigma
        procedure(ComputezInterface), deferred :: Computez
        procedure(ComputeISWInterface), deferred :: ComputeISW
        procedure(ComputeLensingInterface), deferred :: ComputeLensing

        ! other subroutines
        procedure :: Init => TModGravityModel_Init
        !procedure, nopass :: PythonClass => TModGravityModel_PythonClass
        !procedure, nopass :: SelfPointer => TModGravityModel_SelfPointer
        procedure :: ReadParams => TModGravityModel_ReadParams
        ! TODO: after the code works, see if you need to keep this (although it's useful for debugging)
        procedure :: PrintAttributes => TModGravityModel_PrintAttributes
        ! TODO: do we need these?
        procedure :: ResetCache => TModGravityModel_ResetCache
        procedure :: MGCAMB_open_cache_files
        procedure :: MGCAMB_close_cache_files
        procedure :: MGCAMB_dump_cache

    end type TModGravityModel

    ! these functions are actually defined within the subtypes
    abstract interface
        subroutine ComputeMGFunctionsInterface(this, a, k2, adotoa )
            use Precision
            import :: TModGravityModel
            class(TModGravityModel) :: this
            real(dl), intent(in) :: a
            real(dl), intent(in) :: k2
            real(dl), intent(in) :: adotoa
        end subroutine ComputeMGFunctionsInterface

        subroutine ComputesigmaInterface( this, a, k, k2, etak, adotoa, dgpi )
            use Precision
            import :: TModGravityModel
            class(TModGravityModel) :: this
            real(dl), intent(in) :: a
            real(dl), intent(in) :: k, k2 ! k^2
            real(dl), intent(in) :: etak
            real(dl), intent(in) :: adotoa
            real(dl), intent(in) :: dgpi
        end subroutine ComputesigmaInterface

        subroutine ComputezInterface(this, a, k, k2, adotoa, &
                                    & grhoc_t, grhob_t, grhor_t, grhog_t, grhonu_t, gpresnu_t, &
                                    & grhov_t, gpresv_t, &
                                    & dgq, dgpi, dgpi_w_sum, pidot_sum ) 

            use Precision
            import :: TModGravityModel
            class(TModGravityModel) :: this
            real(dl), intent(in) :: a
            real(dl), intent(in) :: k, k2
            real(dl), intent(in) :: adotoa
            real(dl), intent(in) :: grhoc_t, grhob_t, grhor_t, grhog_t, grhonu_t, gpresnu_t
            real(dl), intent(in) :: grhov_t, gpresv_t
            real(dl), intent(in) :: dgq, dgpi, dgpi_w_sum, pidot_sum 

        end subroutine ComputezInterface

        subroutine ComputeISWInterface( this, a, adotoa, k, k2, & 
                                      & grho, gpres, pidot_sum, dgq, dgpi, dgpi_diff )
            use Precision
            import :: TModGravityModel
            class(TModGravityModel) :: this
            real(dl), intent(in) :: a
            real(dl), intent(in) :: adotoa
            real(dl), intent(in) :: k, k2 
            real(dl), intent(in) :: grho, gpres
            real(dl), intent(in) :: pidot_sum
            real(dl), intent(in) :: dgq, dgpi, dgpi_diff
        end subroutine ComputeISWInterface

        subroutine ComputeLensingInterface(this, a)
            use Precision
            import :: TModGravityModel
            class(TModGravityModel) :: this
            real(dl), intent(in) :: a
        end subroutine ComputeLensingInterface

    end interface


    contains

    subroutine TModGravityModel_Init(this)!, State)
        class(TModGravityModel) :: this
        !class(TCAMBdata), intent(in), target :: State
    end subroutine TModGravityModel_Init

    ! function TModGravityModel_PythonClass()
    !     character(LEN=:), allocatable :: TModGravityModel_PythonClass
    !     TModGravityModel_PythonClass = 'ModGravityModel'
    ! end function TModGravityModel_PythonClass

    ! subroutine TModGravityModel_SelfPointer(cptr,P)
    ! use iso_c_binding
    ! Type(c_ptr) :: cptr
    ! Type (TModGravityModel), pointer :: PType
    ! class (TPythonInterfacedClass), pointer :: P
    ! call c_f_pointer(cptr, PType)
    ! P => PType
    ! end subroutine TModGravityModel_SelfPointer

    subroutine TModGravityModel_ReadParams( this, Ini )

        use IniObjects
        class(TModGravityModel) :: this
        class(TIniFile), intent(in) :: Ini

        this%DebugMGCAMB = Ini%Read_Logical( 'DebugMGCAMB', .false. )
        this%debug_root = Ini%Read_String_Default( 'output_root', 'debug_' )
        
        this%GRtrans = Ini%Read_Double( 'GRtrans', 0.001d0 )

        ! TODO: ultimately this flag will only be set internally in the Init
        ! of each child class
        this%MG_flag = Ini%Read_Int( 'MG_flag', 1 )
        this%MGDE_pert = Ini%Read_Logical( 'MGDE_pert', .false. )
        this%CDM_flag = Ini%Read_Int( 'CDM_flag', 1 )
        
    end subroutine TModGravityModel_ReadParams




    ! useful for debugging
    subroutine TModGravityModel_PrintAttributes(this)

        class(TModGravityModel), intent(in) :: this

        write (*,*) 'DebugMGCAMB = ', this%DebugMGCAMB
        write (*,*) 'debug_root = ', this%debug_root

        write (*,*) 'GRtrans = ', this%GRtrans
        write (*,*) 'MGDE_pert = ', this%MGDE_pert
        write (*,*) 'MG_flag = ', this%MG_flag
        write (*,*) 'CDM_flag = ', this%CDM_flag

    end subroutine TModGravityModel_PrintAttributes



    ! TODO: this should be put somewhere more sensible, for now it's here
    !> Subroutine that sets the mgcamb_cache to zero
    subroutine TModGravityModel_ResetCache( this )

        class(TModGravityModel) :: this

        ! 1. Background quantities
        !this%grhov_t    = 0._dl
        !this%gpresv_t   = 0._dl

        ! 2. Perturbation quantities
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
        open(unit=111, file=TRIM(this%debug_root) // 'MGCAMB_debug_sources.dat', status="replace", &
            & action="write")
        write(111,*)  'k  ', 'a  ', 'MG_ISW  ', 'MG_Lensing  ', 'S_T  ', 'S_lensing'

        ! 2 Open MG functions file
        open(unit=222, file=TRIM(this%debug_root) // 'MGCAMB_debug_MG_fncs.dat', status="replace",&
            & action="write")
        write(222,*)  'k  ', 'a  ', 'mu  ', 'gamma ', 'Q ', 'R ', 'Phi ', 'Psi ', 'dPhi ', 'dPsi '

        ! 3. Open Einstein solutions file
        open(unit=333, file=TRIM(this%debug_root) // 'MGCAMB_debug_EinsteinSol.dat', status="replace",&
            & action="write")
        write(333,*) 'k  ', 'a  ', 'etak  ', 'z  ', 'sigma  ', 'etadot  ', 'sigmadot  '

        ! 4. Open Perturbation solution file
        open(unit=444, file=TRIM(this%debug_root) // 'MGCAMB_debug_PerturbSol.dat', status="replace",&
        & action="write")
        write(444,*)  'k  ', 'a  ', 'dgrho  ', 'dgq  ', 'rhoDelta  ', 'dgpi  ', 'pidot_sum  ', 'dgpi_w_sum  '

        ! 5. Open Background file
        open(unit=555, file=TRIM(this%debug_root) // 'MGCAMB_debug_Background.dat', status="replace",&
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
    subroutine MGCAMB_dump_cache( this, a, k, etak, grhov_t, gpresv_t, dgrho, dgq, dgpi, adotoa, pidot_sum, dgpi_w_sum )

        class(TModGravityModel) :: this
        real(dl), intent(in) :: a   !< scale factor
        real(dl), intent(in) :: k
        real(dl), intent(in) :: etak
        real(dl), intent(in) :: grhov_t, gpresv_t
        real(dl), intent(in) :: dgrho, dgq, dgpi
        real(dl), intent(in) :: adotoa, pidot_sum, dgpi_w_sum
        character(*), parameter :: cache_output_format = 'e18.8'

        ! 1. Write the sources
        write(111,'(14'//cache_output_format//')') k, a, this%MG_ISW, this%MG_Lensing, this%source1, this%source3

        ! 2. Write the MG functions and the potentials
        write(222,'(14'//cache_output_format//')') k, a, this%mu, this%gamma, this%q, this%r, &
                                                & this%MG_phi, this%MG_psi, this%MG_phidot, this%MG_psidot

        ! 3. Write the Einstein equations solutions
        write(333,'(14'//cache_output_format//')') k, a, etak, this%z, this%sigma, this%etadot, this%sigmadot

        ! 4. Write the Perturbations Solutions
        write(444,'(14'//cache_output_format//')') k, a, dgrho, dgq, this%rhoDelta,&
                                                    & dgpi, pidot_sum, dgpi_w_sum

        !5. Write the background
        write(555,'(14'//cache_output_format//')') k, a, adotoa, this%Hdot, grhov_t,  gpresv_t

    end subroutine MGCAMB_dump_cache


! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints to screen the MGCAMB header.
    subroutine print_MGCAMB_header

        implicit none

        ! print the header:
        write(*,'(a)') "***************************************************************"
        write(*,'(a)') "     __  _________  ________   __  ______  "
        write(*,'(a)') "    /  \/  / ____/ / ___/ _ | /  |/  / _ ) "
        write(*,'(a)') "   / /\_/ / /_,-, / /__/ __ |/ /|_/ / _  | "
        write(*,'(a)') "  /_/  /_/_____/  \___/_/ |_/_/  /_/____/  "
        write(*,'(a)') "  "
        write(*,'(a)') "        Modified Growth with CAMB "
        write(*,'(a)') "  "
        write(*,'(a)') "***************************************************************"

    end subroutine print_MGCAMB_header



end module ModGravityInterface