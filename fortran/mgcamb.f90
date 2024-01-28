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
        procedure :: ComputeISW => TMuGammaParameterization_ComputeISW
        procedure :: ComputeLensing => TMuGammaParameterization_ComputeLensing

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
    subroutine TMuGammaParameterization_ComputeMu( this, a, k2 )
            
        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  ! scale factor
        real(dl), intent(in) :: k2 ! k^2
        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: mu

        LKA1 = this%lambda1_2 * k2 * a**this%ss

        this%mu = (1._dl + this%B1 * LKA1)/(1._dl + LKA1)  

    end subroutine TMuGammaParameterization_ComputeMu

    ! \dot{mu}(a,k) function
    subroutine TMuGammaParameterization_ComputeMudot( this, a, k2, adotoa )

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  !< scale factor
        real(dl), intent(in) :: k2 ! k^2
        real(dl), intent(in) :: adotoa

        ! local variables
        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: num
        real(dl) :: den

        LKA1 = this%lambda1_2 * k2 * a**this%ss
        num = (this%B1 - 1._dl) * adotoa * this%ss * LKA1
        den = (1._dl+LKA1)**2._dl

        this%mudot =  num / den

    end subroutine TMuGammaParameterization_ComputeMudot

    ! gamma(a,k) function
    subroutine TMuGammaParameterization_ComputeGamma( this, a, k2 )
  
        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  !< scale factor
        real(dl), intent(in) :: k2 ! k^2
        real(dl) :: LKA2 ! \lambda_2^2 k^2 a^s

        LKA2 = this%lambda2_2 * k2 * a**this%ss
        this%gamma = (1._dl + this%B2 * LKA2)/(1._dl + LKA2)

    end subroutine TMuGammaParameterization_ComputeGamma


    subroutine TMuGammaParameterization_ComputeGammadot( this, a, k2, adotoa )
  
        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  !< scale factor
        real(dl), intent(in) :: k2 ! k^2
        real(dl), intent(in) :: adotoa

        real(dl) :: LKA2 ! \lambda_2^2 k^2 a^s
        real(dl) :: num
        real(dl) :: den

        LKA2 = this%lambda2_2 * k2 * a**this%ss
        num = (this%B2 -1._dl)* adotoa * this%ss * LKA2
        den = (1._dl + LKA2)**2._dl
        this%gammadot = num / den

    end subroutine TMuGammaParameterization_ComputeGammadot


    ! Sigma(a,k) function - derived from mu(a,k) and gamma(a,k)
    ! TODO: this could just be made an attribute of `this`
    function TMuGammaParameterization_ComputeBigSigma( this, a, k2 ) result(BigSigma)

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  !< scale factor
        real(dl), intent(in) :: k2 ! k^2

        real(dl) :: mu, gamma, BigSigma

        call this%ComputeMu( a, k2 )
        call this%ComputeGamma( a, k2 )

        BigSigma = this%mu * (1._dl + this%gamma)/2._dl   

    end function TMuGammaParameterization_ComputeBigSigma  


    !dot{Sigma}(a,k) function - obtained from mu, gamma, mudot, gammadot - not used
    ! in the calculation
    function TMuGammaParameterization_ComputeBigSigmadot( this, a, k2, adotoa ) result(BigSigmadot)

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a  !< scale factor
        real(dl), intent(in) :: k2 ! k^2
        real(dl), intent(in) :: adotoa

        real(dl) :: mu, gamma, mudot, gammadot, BigSigmadot

        call this%ComputeMGFunctions( a, k2, adotoa )
        BigSigmadot = (this%mudot * (1._dl + this%gamma) + this%mu * this%gammadot)/2._dl
     
    end function TMuGammaParameterization_ComputeBigSigmadot



    !> this subroutine computes the MG functions at a time-step
    subroutine TMuGammaParameterization_ComputeMGFunctions( this, a, k2, adotoa )

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a   !< scale factor
        real(dl), intent(in) :: k2 ! k^2
        real(dl), intent(in) :: adotoa
    
        call this%ComputeMu( a, k2 )
        call this%ComputeMuDot( a, k2, adotoa )
        call this%ComputeGamma( a, k2 )
        call this%ComputeGammaDot( a, k2, adotoa )

    end subroutine TMuGammaParameterization_ComputeMGFunctions


    !---------------------------------------------------------------------------
    !> this subroutine computes the shear sigma in MG or sigma^star in the notes
    subroutine TMuGammaParameterization_Computesigma( this, a, k, k2, etak, adotoa, dgpi )

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a   !< scale factor
        real(dl), intent(in) :: k, k2 ! k^2
        real(dl), intent(in) :: etak
        real(dl), intent(in) :: adotoa
        real(dl), intent(in) :: dgpi
        
        ! first calculate MG_alpha 
        this%MG_alpha = ( etak / k + this%mu * ( this%gamma * this%rhoDelta + &
                            ( this%gamma -1._dl )*2._dl* dgpi )/(2._dl * k2)) / adotoa

        ! then calculate sigma 
        this%sigma = k * this%MG_alpha  

    end subroutine TMuGammaParameterization_Computesigma    

    !---------------------------------------------------------------------------
    !> this subroutine computes the perturbation Z in MG
    subroutine TMuGammaParameterization_Computez( this, a, k, k2, adotoa, &
                                                & grhoc_t, grhob_t, grhor_t, grhog_t, grhonu_t, gpresnu_t, &
                                                & dgq, dgpi, dgpi_w_sum, pidot_sum  ) 

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a   !< scale factor
        real(dl), intent(in) :: k, k2  ! k^2
        real(dl), intent(in) :: adotoa
        real(dl), intent(in) :: grhoc_t, grhob_t, grhor_t, grhog_t, grhonu_t, gpresnu_t
        real(dl), intent(in) :: dgq, dgpi, dgpi_w_sum, pidot_sum 

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
            fmu = k2 + .5d0 * this%gamma * this%mu * (3._dl * (grhoc_t + grhob_t) &
                & + 4._dl * (grhog_t + grhor_t ) + 3._dl * (grhonu_t + gpresnu_t ) &
                & + 3._dl * (this%grhov_t + this%gpresv_t ))
        else
            fmu = k2 + .5d0*this%gamma * this%mu * (3._dl * (grhoc_t + grhob_t) &
                & + 4._dl * (grhog_t + grhor_t ) + 3._dl * (grhonu_t + gpresnu_t ) )
        end if

        !> adding massive neutrinos contributions
        f1 = k2 + 3._dl * ( adotoa**2 - this%Hdot )
        term1 = this%gamma * this%mu* f1 * dgq / k  

        !> adding massive neutrinos contribution, if w_DE /= -1 this has to be changed
        if(this%MGDE_pert) then
            term2 = k2 * this%MG_alpha * (this%mu* this%gamma*( grhoc_t + grhob_t   &
                    & +(4._dl/3._dl) * ( grhog_t + grhor_t ) + ( grhonu_t + gpresnu_t )  &
                    & + (this%grhov_t + this%gpresv_t))- 2._dl*( adotoa**2 - this%Hdot ))  
        else
            term2 = k2 * this%MG_alpha * (this%mu * this%gamma *( grhoc_t + grhob_t   &
                & + (4._dl/3._dl) * ( grhog_t + grhor_t ) + ( grhonu_t + gpresnu_t ) ) &
                & - 2._dl * ( adotoa**2 - this%Hdot))
        end if

        term3 = ( this%mu * ( this%gamma -1._dl ) * adotoa - this%gamma * this%mudot &
                & - this%gammadot * this%mu ) * this%rhoDelta  

        ! typo corrected here
        term4 = 2._dl * this%mu * (this%gamma - 1._dl) * adotoa * dgpi_w_sum 

        ! separated from the previous term
        term5 = -2._dl * ((this%gamma-1._dl) * this%mudot -this%gammadot * this%mu) * dgpi 

        !> adding massive neutrinos contribution
        term6 = 2._dl * this%mu * (1._dl - this%gamma) * pidot_sum 

        !> calculate etadot
        this%etadot = (term1 + term2 + term3 + term4 + term5 + term6) / ( 2._dl * fmu) 

        !> finally calculate Z
        this%z = this%sigma - 3._dl * this%etadot/ k  

        !> Calculate the Newtonian potential  
        this%MG_psi = - this%mu * ( this%rhoDelta + 2._dl * dgpi)/(2._dl * k2) 

        !> calculate the curvature perturbation potential 
        this%MG_phi = this%gamma * this%MG_psi + this%mu * 1._dl * dgpi/ k2 

        this%MG_phidot = this%etadot - adotoa * (this%MG_psi - adotoa * this%MG_alpha) &
                            & - this%Hdot * this%MG_alpha  

        this%sigmadot = k * (this%MG_psi - adotoa * this%MG_alpha) 

    end subroutine TMuGammaParameterization_Computez 


    !------------------------------------------------------------GaugeInterface_EvolveScal---------------
    !> this subroutine computes the ISW term in MG 
    subroutine TMuGammaParameterization_ComputeISW( this, a, adotoa, k, k2, & 
                                                  & grho, gpres, pidot_sum, dgq, dgpi, dgpi_diff )

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a   !< scale factor
        real(dl), intent(in) :: adotoa
        real(dl), intent(in) :: k, k2  ! k^2
        real(dl), intent(in) :: grho, gpres
        real(dl), intent(in) :: pidot_sum
        real(dl), intent(in) :: dgq, dgpi, dgpi_diff

        !local variables
        real(dl) :: term0

        term0 = k2 + 3._dl* ( adotoa**2._dl - this%Hdot)

        !adding MG_rhoDeltadot
        this%rhoDeltadot = -term0 * dgq / k - ( grho + gpres) * k * this%z &
                            & - adotoa * this%rhoDelta - 2._dl * adotoa * dgpi 

        !adding dgpidot
        this%dgpidot = pidot_sum - (2._dl*dgpi + dgpi_diff )* adotoa  

        !! derived from 1st Possion eq
        this%MG_psidot = - 0.5d0 * this%mu/k2 * ( this%rhoDeltadot + 2._dl*this%dgpidot ) &
                            & - 0.5d0 * this%mudot/k2 * ( this%rhoDelta + 2._dl*dgpi ) 

        this%MG_ISW = this%MG_phidot+this%MG_psidot

        this%MG_alphadot = this%MG_psi - adotoa * this%MG_alpha

    end subroutine TMuGammaParameterization_ComputeISW

    !---------------------------------------------------------------------------
    !> this subroutine computes the lensing term in MG
    subroutine TMuGammaParameterization_ComputeLensing( this, a )

        class(TMuGammaParameterization) :: this
        real(dl), intent(in) :: a   !< scale factor

        this%MG_lensing = this%MG_phi + this%MG_psi

    end subroutine TMuGammaParameterization_ComputeLensing

    !---------------------------------------------------------------------------
    !> this subroutine prints the class attributes - useful for debugging
    subroutine TMuGammaParameterization_PrintAttributes( this )

        class(TMuGammaParameterization), intent(in) :: this 

        call TModGravityModel_PrintAttributes(this)

        write (*,*) 'B1 = ', this%B1
        write (*,*) 'B2 = ', this%B2
        write (*,*) 'lambda1_2 = ', this%lambda1_2
        write (*,*) 'lambda2_2 = ', this%lambda2_2
        write (*,*) 'ss = ', this%ss

    end subroutine TMuGammaParameterization_PrintAttributes

end module MGCAMB

