module MGCAMB
    use precision
    use ModGravityInterface
    use splines

    implicit none

    ! TODO: when some function is not relevant to a parameterization/model
    ! create it and print an error code if it is called

    ! TODO: it's probably better to retain the MG flags somewhere to test what model it is

    ! =========== PURE MG MODELS ============

    ! ============ mu, gamma parameterisation ============
    type, extends(TModGravityModel) :: TMuGammaParameterization

        real(dl):: B1 = 1.333d0
        real(dl):: B2 = 0.5
        real(dl):: lambda1_2  = 1000
        real(dl):: lambda2_2 = 1000
        real(dl) :: ss

        contains

        ! deferred procedures
        procedure :: Computesigma => TMuGammaParameterization_Computesigma
        procedure :: ComputeMGFunctions => TMuGammaParameterization_ComputeMGFunctions
        procedure :: Computez => TMuGammaParameterization_Computez
        procedure :: ComputeISWAndLensing => TMuGammaParameterization_ComputeISWAndLensing

        procedure :: ReadParams => TMuGammaParameterization_ReadParams
        procedure, nopass :: PythonClass => TMuGammaParameterization_PythonClass
        procedure, nopass :: SelfPointer => TMuGammaParameterization_SelfPointer
        procedure :: Init => TMuGammaParameterization_Init
        procedure :: ComputeMu => TMuGammaParameterization_ComputeMu
        procedure :: ComputeMudot => TMuGammaParameterization_ComputeMudot
        procedure :: ComputeGamma => TMuGammaParameterization_ComputeGamma
        procedure :: ComputeGammadot => TMuGammaParameterization_ComputeGammadot
        procedure :: ComputeBigSigma => TMuGammaParameterization_ComputeBigSigma
        procedure :: ComputeBigSigmadot => TMuGammaParameterization_ComputeBigSigmadot
        procedure :: PrintAttributes => TMuGammaParameterization_PrintAttributes

    end type TMuGammaParameterization


    contains


    ! CACCA: devo finire di settare read params in both TmodGrav Model and here
    ! then set MGFLAG within the child class, make sure init works

    subroutine TMuGammaParameterization_ReadParams( this, Ini )
        use IniObjects
        class(TMuGammaParameterization) :: this
        class(TIniFile), intent(in) :: Ini

        call TModGravityModel_ReadParams( this, Ini )

        this%B1 = Ini%Read_Double( 'B1', 1.333d0 )
        this%B2 = Ini%Read_Double( 'B2', 0.5d0 )
        this%lambda1_2 = Ini%Read_Double( 'lambda1_2', 1000.d0 )
        this%lambda2_2 = Ini%Read_Double( 'lambda2_2', 1000.d0 )
        this%ss = Ini%Read_Double( 'ss', 4.d0 )

    end subroutine TMuGammaParameterization_ReadParams

    function TMuGammaParameterization_PythonClass()
        character(LEN=:), allocatable :: TMuGammaParameterization_PythonClass
        TMuGammaParameterization_PythonClass = 'MuGammaParameterization'
    end function TMuGammaParameterization_PythonClass

    subroutine TMuGammaParameterization_SelfPointer(cptr,P)
        use iso_c_binding
        use classes, only: TPythonInterfacedClass
        Type(c_ptr) :: cptr
        Type (TMuGammaParameterization), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P
        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TMuGammaParameterization_SelfPointer

    subroutine TMuGammaParameterization_Init(this)!, State)
        class(TMuGammaParameterization) :: this
        !class(TCAMBdata), intent(in), target :: State

        ! to keep track of models in different parts of the code
        this%MG_flag = 1

    end subroutine TMuGammaParameterization_Init


    !> mu(a,k) function
    subroutine TMuGammaParameterization_ComputeMu( this, a, k2, mu )
            
        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  ! scale factor
        real(dl), intent(in) :: k2 ! k^2
        
        real(dl), intent(out) :: mu

        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s

        LKA1 = this%lambda1_2 * k2 * a**this%ss

        mu = (1._dl + this%B1 * LKA1)/(1._dl + LKA1)  

    end subroutine TMuGammaParameterization_ComputeMu

    ! \dot{mu}(a,k) function
    subroutine TMuGammaParameterization_ComputeMudot( this, a, k2, adotoa, mudot )

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  !< scale factor
        real(dl), intent(in) :: k2 ! k^2
        real(dl), intent(in) :: adotoa

        real(dl), intent(out) :: mudot

        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: num
        real(dl) :: den

        LKA1 = this%lambda1_2 * k2 * a**this%ss
        num = (this%B1 - 1._dl) * adotoa * this%ss * LKA1
        den = (1._dl+LKA1)**2._dl

        mudot =  num / den

    end subroutine TMuGammaParameterization_ComputeMudot

    ! gamma(a,k) function
    subroutine TMuGammaParameterization_ComputeGamma( this, a, k2, gamma )
  
        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  !< scale factor
        real(dl), intent(in) :: k2 ! k^2

        real(dl), intent(out) :: gamma

        real(dl) :: LKA2 ! \lambda_2^2 k^2 a^s

        LKA2 = this%lambda2_2 * k2 * a**this%ss
        gamma = (1._dl + this%B2 * LKA2)/(1._dl + LKA2)

    end subroutine TMuGammaParameterization_ComputeGamma


    subroutine TMuGammaParameterization_ComputeGammadot( this, a, k2, adotoa, gammadot )
  
        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  !< scale factor
        real(dl), intent(in) :: k2 ! k^2
        real(dl), intent(in) :: adotoa

        real(dl), intent(out) :: gammadot

        real(dl) :: LKA2 ! \lambda_2^2 k^2 a^s
        real(dl) :: num
        real(dl) :: den

        LKA2 = this%lambda2_2 * k2 * a**this%ss
        num = (this%B2 -1._dl)* adotoa * this%ss * LKA2
        den = (1._dl + LKA2)**2._dl
        gammadot = num / den

    end subroutine TMuGammaParameterization_ComputeGammadot


    ! Sigma(a,k) function - derived from mu(a,k) and gamma(a,k)
    function TMuGammaParameterization_ComputeBigSigma( this, a, k2 ) result(BigSigma)

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  !< scale factor
        real(dl), intent(in) :: k2 ! k^2

        real(dl) :: mu, gamma, BigSigma

        call this%ComputeMu( a, k2, mu )
        call this%ComputeGamma( a, k2, gamma )

        BigSigma = mu * (1._dl + gamma)/2._dl   

    end function TMuGammaParameterization_ComputeBigSigma  


    ! dot{Sigma}(a,k) function - obtained from mu, gamma, mudot, gammadot
    ! postprocessing purposes only, not used in the calculation
    function TMuGammaParameterization_ComputeBigSigmadot( this, a, k2, adotoa ) result(BigSigmadot)

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  !< scale factor
        real(dl), intent(in) :: k2 ! k^2
        real(dl), intent(in) :: adotoa

        real(dl) :: mu, gamma, mudot, gammadot, BigSigmadot

        call this%ComputeMu( a, k2, mu )
        call this%ComputeGamma( a, k2, gamma )
        call this%ComputeMudot( a, k2, adotoa, mudot )
        call this%ComputeGammadot( a, k2, adotoa, gammadot )

        BigSigmadot = (mudot * (1._dl + gamma) + mu * gammadot)/2._dl
     
    end function TMuGammaParameterization_ComputeBigSigmadot



    !> this subroutine computes the MG functions at a time-step
    subroutine TMuGammaParameterization_ComputeMGFunctions( this, a, k2, adotoa, &
                                                            & mu, mudot, gamma, gammadot, &
                                                            & m, mdot, beta, betadot )
        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a   !< scale factor
        real(dl), intent(in) :: k2 ! k^2
        real(dl), intent(in) :: adotoa

        ! MG functions - relevant to this parameterisation
        real(dl), intent(out) :: mu, mudot, gamma, gammadot
        ! MG functions that are unused in this parameterisation - set to 0
        real(dl), intent(out) :: m, mdot, beta, betadot

        m = 0d0
        mdot = 0d0
        beta = 0d0
        betadot = 0d0

        call this%ComputeMu( a, k2, mu )
        call this%ComputeMuDot( a, k2, adotoa, mudot )
        call this%ComputeGamma( a, k2, gamma )
        call this%ComputeGammaDot( a, k2, adotoa, gammadot )

    end subroutine TMuGammaParameterization_ComputeMGFunctions


    !---------------------------------------------------------------------------
    !> this subroutine computes the shear sigma in MG or sigma^star in the notes
    subroutine TMuGammaParameterization_Computesigma( this, a, k, k2, etak, adotoa, dgpi, &
                                                    & mu, gamma, rhoDelta, &
                                                    & sigma, MG_alpha )

        ! TODO: design consideration: maybe mu, gamma, etc should just be computed
        ! inside each function that needs them

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a   !< scale factor
        real(dl), intent(in) :: k, k2 ! k^2
        real(dl), intent(in) :: etak, adotoa, dgpi
        real(dl), intent(in) :: mu, gamma, rhoDelta

        real(dl), intent(out) :: sigma, MG_alpha
        
        MG_alpha = ( etak / k + mu * ( gamma * rhoDelta + &
                            ( gamma -1._dl )*2._dl* dgpi )/(2._dl * k2)) / adotoa

        sigma = k * MG_alpha  

    end subroutine TMuGammaParameterization_Computesigma    

    !---------------------------------------------------------------------------
    !> this subroutine computes the perturbation Z in MG
    subroutine TMuGammaParameterization_Computez( this, a, k, k2, adotoa, Hdot, rhoDelta, &
                                                & grhoc_t, grhob_t, grhor_t, grhog_t, grhonu_t, gpresnu_t, &
                                                & grhov_t, gpresv_t, &
                                                & dgq, dgpi, dgpi_w_sum, pidot_sum, &
                                                & mu, gamma, mudot, gammadot, sigma, MG_alpha, &
                                                & MG_phi, MG_psi, MG_phidot, z, sigmadot, etadot ) 

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a   !< scale factor
        real(dl), intent(in) :: k, k2  ! k^2
        real(dl), intent(in) :: adotoa, Hdot, rhoDelta
        real(dl), intent(in) :: grhoc_t, grhob_t, grhor_t, grhog_t, grhonu_t, gpresnu_t
        real(dl), intent(in) :: grhov_t, gpresv_t
        real(dl), intent(in) :: dgq, dgpi, dgpi_w_sum, pidot_sum 
        real(dl), intent(in) :: mu, gamma, mudot, gammadot, sigma, MG_alpha

        real(dl), intent(out) :: MG_phi, MG_psi, MG_phidot, z, sigmadot, etadot

        real(dl) :: fmu
        real(dl) :: f1
        real(dl) :: term1
        real(dl) :: term2
        real(dl) :: term3
        real(dl) :: term4
        real(dl) :: term5
        real(dl) :: term6

        !> adding the massive neutrinos contibutions, but no DE parts
        if (this%MGDE_pert) then
            fmu = k2 + .5d0 * gamma * mu * (3._dl * (grhoc_t + grhob_t) &
                & + 4._dl * ( grhog_t + grhor_t ) + 3._dl * (grhonu_t + gpresnu_t ) &
                & + 3._dl * ( grhov_t + gpresv_t ))
        else
            fmu = k2 + .5d0 * gamma * mu * (3._dl * (grhoc_t + grhob_t) &
                & + 4._dl * ( grhog_t + grhor_t ) + 3._dl * (grhonu_t + gpresnu_t ) )
        end if

        !> adding massive neutrinos contributions
        f1 = k2 + 3._dl * ( adotoa**2 - Hdot )
        term1 = gamma * mu * f1 * dgq / k  

        !> adding massive neutrinos contribution, if w_DE /= -1 this has to be changed
        if (this%MGDE_pert) then
            term2 = k2 * MG_alpha * ( mu * gamma*( grhoc_t + grhob_t   &
                    & +(4._dl/3._dl) * ( grhog_t + grhor_t ) + ( grhonu_t + gpresnu_t )  &
                    & + ( grhov_t + gpresv_t )) -2._dl*( adotoa**2 - Hdot ))  
        else
            term2 = k2 * MG_alpha * ( mu * gamma *( grhoc_t + grhob_t   &
                & + (4._dl/3._dl) * ( grhog_t + grhor_t ) + ( grhonu_t + gpresnu_t ) ) &
                & - 2._dl * ( adotoa**2 - Hdot))
        end if

        term3 = ( mu * ( gamma -1._dl ) * adotoa - gamma * mudot &
                & - gammadot * mu ) * rhoDelta  

        ! typo corrected here
        term4 = 2._dl * mu * (gamma - 1._dl) * adotoa * dgpi_w_sum 

        ! separated from the previous term
        term5 = -2._dl * ((gamma-1._dl) * mudot - gammadot * mu) * dgpi 

        !> adding massive neutrinos contribution
        term6 = 2._dl * mu * (1._dl - gamma) * pidot_sum 

        !> calculate etadot
        etadot = (term1 + term2 + term3 + term4 + term5 + term6) / ( 2._dl * fmu) 

        !> finally calculate Z
        z = sigma - 3._dl * etadot/ k  

        !> Calculate the Newtonian potential  
        MG_psi = - mu * ( rhoDelta + 2._dl * dgpi)/(2._dl * k2) 

        !> calculate the curvature perturbation potential 
        MG_phi = gamma * MG_psi + mu * 1._dl * dgpi/ k2 

        MG_phidot = etadot - adotoa * (MG_psi - adotoa * MG_alpha) - Hdot * MG_alpha  

        sigmadot = k * (MG_psi - adotoa * MG_alpha) 

    end subroutine TMuGammaParameterization_Computez 


    !------------------------------------------------------------GaugeInterface_EvolveScal---------------
    !> this subroutine computes the ISW term in MG 
    subroutine TMuGammaParameterization_ComputeISWAndLensing( this, k, k2, a, adotoa, Hdot, &
                                                            & rhoDelta, mu, mudot, &
                                                            & grho, gpres, pidot_sum, &
                                                            & z, MG_phi, MG_psi, MG_alpha, MG_phidot, & 
                                                            & dgq, dgpi, dgpi_diff, &
                                                            & MG_alphadot, MG_psidot, MG_ISW, MG_lensing )

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: k, k2  ! k^2
        real(dl), intent(in) :: a, adotoa, Hdot
        real(dl), intent(in) :: rhoDelta, mu, mudot 
        real(dl), intent(in) :: grho, gpres
        real(dl), intent(in) :: pidot_sum
        real(dl), intent(in) :: z, MG_phi, MG_psi, MG_alpha, MG_phidot
        real(dl), intent(in) :: dgq, dgpi, dgpi_diff

        real(dl), intent(out) :: MG_alphadot, MG_psidot, MG_ISW, MG_lensing

        !local variables
        real(dl) :: term0, dgpidot, rhoDeltadot

        term0 = k2 + 3._dl* ( adotoa**2._dl - Hdot)

        !adding MG_rhoDeltadot
        rhoDeltadot = -term0 * dgq / k - ( grho + gpres ) * k * z - adotoa * rhoDelta - 2._dl * adotoa * dgpi 

        !adding dgpidot
        dgpidot = pidot_sum - (2._dl*dgpi + dgpi_diff )* adotoa  

        !! derived from 1st Possion eq
        MG_psidot = - 0.5d0 * mu/k2 * ( rhoDeltadot + 2._dl*dgpidot ) &
                            & - 0.5d0 * mudot/k2 * ( rhoDelta + 2._dl*dgpi ) 

        MG_ISW = MG_phidot + MG_psidot
        MG_alphadot = MG_psi - adotoa * MG_alpha
        MG_lensing = MG_phi + MG_psi

    end subroutine TMuGammaParameterization_ComputeISWAndLensing


    !---------------------------------------------------------------------------
    !> this subroutine prints the class attributes - useful for debugging
    subroutine TMuGammaParameterization_PrintAttributes( this )

        class(TMuGammaParameterization), intent(in) :: this 

        call TModGravityModel_PrintAttributes(this)

        write (*,*) 'mu-gamma parameterization, BZ parameterization'
        write (*,*) 'B1 = ', this%B1
        write (*,*) 'B2 = ', this%B2
        write (*,*) 'lambda1_2 = ', this%lambda1_2
        write (*,*) 'lambda2_2 = ', this%lambda2_2
        write (*,*) 'ss = ', this%ss

    end subroutine TMuGammaParameterization_PrintAttributes

end module MGCAMB

