module MGCAMB
    use precision
    use splines

    ! new model selection flags
    integer :: MG_flag
    integer :: pure_MG_flag
    integer :: alt_MG_flag
    integer :: QSA_flag
    integer :: mugamma_par  
    integer :: muSigma_par  
    integer :: QR_par
    integer :: muSigma_flag
    integer :: CDM_flag

    ! DE model flag
    integer :: DE_model

    real(dl) :: GRtrans                     !< scale factor at which MG is switched on

    ! BZ parametrization (and QS f(R))
    real(dl) :: B1
    real(dl) :: B2
    real(dl) :: lambda1_2
    real(dl) :: lambda2_2
    real(dl) :: ss

    ! Planck Parametrization
    real(dl) :: E11
    real(dl) :: E22

    ! Q-R parametrization 1
    real(dl) :: MGQfix
    real(dl) :: MGRfix

    ! Q-R parametrization 2
    real(dl) :: Qnot
    real(dl) :: Rnot
    real(dl) :: sss

    ! Growth rate gamma
    real(dl) :: Linder_gamma

    ! Symmetron
    real(dl) :: beta_star
    real(dl) :: a_star
    real(dl) :: xi_star

    ! Dilaton
    real(dl) :: beta0
    real(dl) :: xi0
    real(dl) :: DilR
    real(dl) :: DilS

    ! Hu-Sawicki f(R) gravity
    real(dl) :: F_R0
    real(dl) :: FRn

    ! DES parametrization
    real(dl) :: mu0
    real(dl) :: sigma0


    ! effective Newton's constant  !! not sure
    real(dl) :: ga
    real(dl) :: nn

    ! DE model parameters

    real(dl) :: w0DE              !< w0 parameters for DE
    real(dl) :: waDE              !< waDE parameters for DE

	logical :: MGDE_const = .True.
    logical :: MGDE_pert = .False.

    !character(len=(10)) :: MGCAMB_version = 'v 4.0'


! =============MGXrecon=============
	integer, parameter  :: nnode=10 ! number of fine bins
	real(dl), parameter :: zstart=3.d0, ztanh=4.d0
	real(dl), parameter :: astart=1.d0/(1.d0+zstart), atanh=1.d0/(1.d0+ztanh)
	real(dl), parameter :: aend=1.d0

	real(dl) :: a_arr(2*nnode)
	real(dl) :: mu_arr(2*nnode) =  0.d0
	real(dl) :: dmu_arr(2*nnode), ddmu_arr(2*nnode), dddmu_arr(2*nnode)
	real(dl) :: sigma_arr(2*nnode) = 0.d0
	real(dl) :: gamma_arr(2*nnode), dgamma_arr(2*nnode), ddgamma_arr(2*nnode), dddgamma_arr(2*nnode)
	real(dl) :: X_arr(2*nnode) = 0.d0
	real(dl) :: dX_arr(2*nnode), ddX_arr(2*nnode), dddX_arr(2*nnode)
! =============MGXrecon=============

    ! define the type MGCAMB_par_cache
    type :: MGCAMB_parameter_cache
        real(dl) :: omegab
        real(dl) :: omegac
        real(dl) :: omegav
        real(dl) :: h0
        real(dl) :: h0_Mpc   

        character(len=30) :: output_root
    end type MGCAMB_parameter_cache

    type(MGCAMB_parameter_cache) :: mgcamb_par_cache

    ! define the tyoe MGCAMB_timestep_cache
    type :: MGCAMB_timestep_cache 

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
        real(dl) :: mu
        real(dl) :: mudot
        real(dl) :: gamma
        real(dl) :: gammadot
        real(dl) :: q
        real(dl) :: qdot
        real(dl) :: r
        real(dl) :: rdot
        real(dl) :: BigSigma
        real(dl) :: BigSigmadot
        real(dl) :: C_phi
        real(dl) :: C_phidot
    

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

    end type MGCAMB_timestep_cache


#ifdef DEBUG
    logical , parameter :: DebugMGCAMB = .true.              !< MGCAMB debug flag.This will turn on printing of many things to aid debugging the code.
#else
    logical , parameter :: DebugMGCAMB = .false.             !< MGCAMB debug flag.This will turn on printing of many things to aid debugging the code.
#endif


contains

! =============MGXrecon=============
	subroutine reconstruction_arr

		use precision
        !use Interpolation
        !use splines
		implicit none

		integer :: i
		real(dl),parameter :: d0lo=1.d32, d0hi=1.d32
			
		  do i=nnode, 2*nnode
		   gamma_arr(i) = 2.d0*sigma_arr(i)/mu_arr(i)-1.d0
		  end do
		
		  do i=1, nnode
		   a_arr(i)        = atanh*dble(i-1)/dble(nnode-1)
		   a_arr(nnode+i)  = astart+(1.d0-astart)*dble(i-1)/dble(nnode-1)
			if(i<nnode) then
			  mu_arr(i) = (mu_arr(nnode)-1.d0)/2.d0*(1.d0+tanh((a_arr(i)-atanh/2.d0)/0.04d0))+1.d0
					  gamma_arr(i) = (gamma_arr(nnode)-1.d0)/2.d0*(1.d0+tanh((a_arr(i)-atanh/2.d0)/0.04d0))+1.d0
					  !X_arr(i) = (X_arr(nnode)-1.d0)/2.d0*(1.d0+tanh((a_arr(i)-atanh/2.d0)/0.04d0))+1.d0
					  X_arr(i) = (X_arr(nnode)-mgcamb_par_cache%omegav)/2.d0*(1.d0+tanh((a_arr(i)-atanh/2.d0)/0.04d0))+mgcamb_par_cache%omegav !SP:X is CP%omegav at high z
			end if
		  end do

		!   call spline(a_arr,mu_arr,2*nnode,d0lo,d0hi,ddmu_arr)
		!   call spline_deriv(a_arr,mu_arr,ddmu_arr,dmu_arr,2*nnode)
		!   call spline(a_arr,dmu_arr,2*nnode,d0lo,d0hi,dddmu_arr)
		  
		!   call spline(a_arr,gamma_arr,2*nnode,d0lo,d0hi,ddgamma_arr)
		!   call spline_deriv(a_arr,gamma_arr,ddgamma_arr,dgamma_arr,2*nnode)
		!   call spline(a_arr,dgamma_arr,2*nnode,d0lo,d0hi,dddgamma_arr)
		  
		!   call spline(a_arr,X_arr,2*nnode,d0lo,d0hi,ddX_arr)
		!   call spline_deriv(a_arr,X_arr,ddX_arr,dX_arr,2*nnode)
		!   call spline(a_arr,dX_arr,2*nnode,d0lo,d0hi,dddX_arr)	  

		  call spline_def(a_arr,mu_arr,2*nnode,ddmu_arr)
		  call spline_deriv(a_arr,mu_arr,ddmu_arr,dmu_arr,2*nnode)
		  call spline_def(a_arr,dmu_arr,2*nnode,dddmu_arr)
		  
		  call spline_def(a_arr,gamma_arr,2*nnode,ddgamma_arr)
		  call spline_deriv(a_arr,gamma_arr,ddgamma_arr,dgamma_arr,2*nnode)
		  call spline_def(a_arr,dgamma_arr,2*nnode,dddgamma_arr)
		  
		  call spline_def(a_arr,X_arr,2*nnode,ddX_arr)
		  call spline_deriv(a_arr,X_arr,ddX_arr,dX_arr,2*nnode)
		  call spline_def(a_arr,dX_arr,2*nnode,dddX_arr)


	end subroutine reconstruction_arr
! =============MGXrecon=============

    !> this subroutine computes the MG functions at a time-step
    subroutine MGCAMB_compute_MG_functions( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

       

        ! Divide the cases here
        if (( MG_flag == 1 .and. pure_MG_flag /= 3 ) &   ! generic mu-gamma parametrization
            .or. MG_flag == 2 &
            .or. MG_flag == 3 .or. MG_flag == 6) then


            mg_cache%mu         = MGCAMB_Mu( a, mg_par_cache, mg_cache )
            mg_cache%mudot      = MGCAMB_MuDot( a, mg_par_cache, mg_cache )
            mg_cache%gamma      = MGCAMB_Gamma( a, mg_par_cache, mg_cache )
            mg_cache%gammadot   = MGCAMB_GammaDot( a, mg_par_cache, mg_cache )


        ! other EFT functions are zero
            mg_cache%q      =      0._dl
            mg_cache%qdot   =      0._dl
            mg_cache%r      =      0._dl
            mg_cache%rdot   =      0._dl
            mg_cache%C_phi   =      0._dl
            mg_cache%C_phidot=      0._dl
            mg_cache%BigSigma =    0._dl
            mg_cache%BigSigmadot = 0._dl

        else if (MG_flag == 4)  then ! only-CDM coupling

            if(CDM_flag == 1) then  !CDM QSA

                mg_cache%C_phi  = MGCAMB_C_phi( a, mg_par_cache, mg_cache )
                mg_cache%C_phidot  = MGCAMB_C_phidot( a, mg_par_cache, mg_cache )

            ! other MG functions are zero
                mg_cache%q           =      0._dl
                mg_cache%qdot        =      0._dl
                mg_cache%r           =      0._dl
                mg_cache%rdot        =      0._dl
                mg_cache%mu          =      0._dl
                mg_cache%mudot       =      0._dl
                mg_cache%gamma       =      0._dl
                mg_cache%gammadot    =      0._dl
                mg_cache%BigSigma    =      0._dl
                mg_cache%BigSigmadot =      0._dl
         
            end if            


        else if (MG_flag == 5)  then !direct mu-Sigma parametrization

            mg_cache%mu         = MGCAMB_Mu( a, mg_par_cache, mg_cache )
            mg_cache%mudot      = MGCAMB_MuDot( a, mg_par_cache, mg_cache )
            mg_cache%BigSigma  = MGCAMB_BigSigma( a, mg_par_cache, mg_cache )
            mg_cache%BigSigmadot   = MGCAMB_BigSigmadot( a, mg_par_cache, mg_cache )
            mg_cache%gamma       =  MGCAMB_Gamma( a, mg_par_cache, mg_cache )
            mg_cache%gammadot    =  MGCAMB_Gammadot( a, mg_par_cache, mg_cache )

            ! other MG functions are zero
            mg_cache%q           =      0._dl
            mg_cache%qdot        =      0._dl
            mg_cache%r           =      0._dl
            mg_cache%rdot        =      0._dl
            mg_cache%C_phi        =      0._dl
            mg_cache%C_phidot     =      0._dl
    

        else if (  MG_flag == 1 .and. pure_MG_flag == 3  ) then ! the Q,R parametrization

            mg_cache%q      = MGCAMB_Q( a, mg_par_cache, mg_cache )
            mg_cache%qdot   = MGCAMB_Qdot( a, mg_par_cache, mg_cache )
            mg_cache%r      = MGCAMB_R( a, mg_par_cache, mg_cache )
            mg_cache%rdot   = MGCAMB_Rdot( a, mg_par_cache, mg_cache )

            ! other MG functions are zero
            mg_cache%mu            = 0._dl
            mg_cache%mudot         = 0._dl
            mg_cache%gamma         = 0._dl
            mg_cache%gammadot      = 0._dl
            mg_cache%BigSigma      = 0._dl
            mg_cache%BigSigmadot   = 0._dl
            mg_cache%C_phi           = 0._dl
            mg_cache%C_phidot        = 0._dl
            
        end if

    end subroutine MGCAMB_compute_MG_functions


    !---------------------------------------------------------------------------
    !> this subroutine computes the shear sigma in MG or sigma^star in the notes
    subroutine MGCAMB_compute_sigma( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        if (( MG_flag == 1 .and. pure_MG_flag /= 3 ) & ! mu-gamma 
            .or. MG_flag == 2 &
            .or. MG_flag == 3 .or. MG_flag == 6) then

            ! first calculate MG_alpha 
            mg_cache%MG_alpha = ( mg_cache%etak/mg_cache%k + mg_cache%mu * ( mg_cache%gamma*mg_cache%rhoDelta+ &
                                ( mg_cache%gamma- 1._dl )*2._dl* mg_cache%dgpi)/(2._dl*mg_cache%k2)) / mg_cache%adotoa

            ! then calculate sigma 
            mg_cache%sigma = mg_cache%k * mg_cache%MG_alpha  

        else if (  MG_flag == 1 .and. pure_MG_flag == 3  ) then ! Q-R parametrization

            mg_cache%MG_phi      = - mg_cache%rhoDelta * mg_cache%q/(2._dl*mg_cache%k2)  
            mg_cache%sigma       = (mg_cache%etak - mg_cache%k * mg_cache%MG_phi)/mg_cache%adotoa  
            mg_cache%MG_alpha    = mg_cache%sigma/mg_cache%k 


        else if ( MG_flag == 4) then ! only-CDM coupling

            if (CDM_flag == 1) then ! CDM QSA

            
                mg_cache%MG_alpha = (mg_cache%etak/mg_cache%k + (mg_cache%rhoDelta + (mg_cache%C_phi)  &
                    & * mg_cache%rhoDeltac)/(2._dl*mg_cache%k2))/mg_cache%adotoa

                mg_cache%sigma = mg_cache%k * mg_cache%MG_alpha

            end if 

        else if (MG_flag == 5) then ! direct mu-Sigma parametrization

            mg_cache%MG_alpha = (mg_cache%etak/mg_cache%k + ((2._dl*mg_cache%BigSigma - mg_cache%mu)*mg_cache%rhoDelta &
                    +(mg_cache%BigSigma - mg_cache%mu)*2._dl* mg_cache%dgpi)/(2._dl*mg_cache%k2))/mg_cache%adotoa

            mg_cache%sigma = mg_cache%k * mg_cache%MG_alpha

        end if

    end subroutine MGCAMB_compute_sigma    

    !---------------------------------------------------------------------------
    !> this subroutine computes the perturbation Z in MG
    subroutine MGCAMB_compute_z( a, mg_par_cache, mg_cache ) 
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        !> other parameters
        real(dl) :: fmu
        real(dl) :: f1
        real(dl) :: fQ
        real(dl) :: fs
        real(dl) :: term1
        real(dl) :: term2
        real(dl) :: term3
        real(dl) :: term4
        real(dl) :: term5
        real(dl) :: term6
        real(dl) :: k2alpha
      
        real(dl) :: m,beta,betadot 

        m = MGCAMB_M(a, mg_par_cache, mg_cache)
        beta = MGCAMB_beta(a, mg_par_cache, mg_cache)
        betadot = MGCAMB_betadot(a, mg_par_cache, mg_cache)

        if (( MG_flag == 1 .and. pure_MG_flag /= 3 ) & 
            .or. MG_flag == 2 &
            .or. MG_flag == 3 .or. MG_flag == 6) then

            !> adding the massive neutrinos contibutions, but no DE parts
            !fmu = mg_cache%k2+0.5d0*mg_cache%gamma*mg_cache%mu*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
            !    & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) +3._dl * (mg_cache%grhonu_t + mg_cache%gpresnu_t )) 
            if(MG_flag == 1 .and. MGDE_pert) then
                fmu = mg_cache%k2+0.5d0*mg_cache%gamma*mg_cache%mu*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                    & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) +3._dl * (mg_cache%grhonu_t + mg_cache%gpresnu_t ) &
                    & + 3._dl * (mg_cache%grhov_t + mg_cache%gpresv_t ))
            else
                fmu = mg_cache%k2+0.5d0*mg_cache%gamma*mg_cache%mu*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                    & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) +3._dl * (mg_cache%grhonu_t + mg_cache%gpresnu_t ))
            end if
            !> adding massive neutrinos contributions

            f1 = mg_cache%k2+3._dl*( mg_cache%adotoa**2 - mg_cache%Hdot )  
            !f1 = mg_cache%k2+0.5d0*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
            !    & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) + 3._dl*(mg_cache%grhonu_t+mg_cache%gpresnu_t) &
            !    & + 3._dl*(mg_cache%grhov_t+mg_cache%gpresv_t))

            term1 = mg_cache%gamma*mg_cache%mu* f1 * mg_cache%dgq/mg_cache%k  

            !> adding massive neutrinos contribution, if w_DE /= -1 this has to be changed

            !term2 = mg_cache%k2*mg_cache%MG_alpha* (mg_cache%mu* mg_cache%gamma*( mg_cache%grhoc_t+mg_cache%grhob_t   &
            !        & +(4._dl/3._dl)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t) ) &
            !        & - 2._dl*(mg_cache%adotoa**2 - mg_cache%Hdot))  
            if(MG_flag == 1 .and. MGDE_pert) then
                term2 = mg_cache%k2*mg_cache%MG_alpha* (mg_cache%mu* mg_cache%gamma*( mg_cache%grhoc_t+mg_cache%grhob_t   &
                        & +(4._dl/3._dl)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t)  &
                        & + (mg_cache%grhov_t + mg_cache%gpresv_t))- 2._dl*(mg_cache%adotoa**2 - mg_cache%Hdot))  
            else
                term2 = mg_cache%k2*mg_cache%MG_alpha* (mg_cache%mu* mg_cache%gamma*( mg_cache%grhoc_t+mg_cache%grhob_t   &
                    & +(4._dl/3._dl)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t) ) &
                    & - 2._dl*(mg_cache%adotoa**2 - mg_cache%Hdot))
            end if

            term3= (mg_cache%mu * ( mg_cache%gamma -1._dl)* mg_cache%adotoa - mg_cache%gamma*mg_cache%mudot &
                    & - mg_cache%gammadot*mg_cache%mu )*mg_cache%rhoDelta  

            ! typo corrected here
            term4 = 2._dl*mg_cache%mu*(mg_cache%gamma - 1._dl)*mg_cache%adotoa*mg_cache%dgpi_w_sum 

            ! separated from the previous term
            term5 = -2._dl*((mg_cache%gamma-1._dl)*mg_cache%mudot -mg_cache%gammadot*mg_cache%mu)*mg_cache%dgpi 

            !> adding massive neutrinos contribution
            term6= 2._dl * mg_cache%mu*(1._dl - mg_cache%gamma)* mg_cache%pidot_sum 

            !> calculate etadot
            mg_cache%etadot = (term1 + term2 + term3 + term4 + term5 + term6)/( 2._dl * fmu) 

            !> finally calculate Z
            mg_cache%z = mg_cache%sigma - 3._dl * mg_cache%etadot/mg_cache%k  

            !> Calculate the Newtonian potential  
            mg_cache%MG_psi = - mg_cache%mu * ( mg_cache%rhoDelta + 2._dl* mg_cache%dgpi)/(2._dl*mg_cache%k2) 

            !> calculate the curvature perturbation potential 
            mg_cache%MG_phi = mg_cache%gamma * mg_cache%MG_psi + mg_cache%mu* 1._dl*mg_cache%dgpi/mg_cache%k2 

            mg_cache%MG_phidot = mg_cache%etadot - mg_cache%adotoa * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha) &
                                & - mg_cache%Hdot * mg_cache%MG_alpha  

        else if (  MG_flag == 1 .and. pure_MG_flag == 3  ) then  

            ! adding massive neutrinos contributions
            fQ = mg_cache%k2 + 0.5d0*mg_cache%q * (3._dl*(mg_cache%grhob_t+mg_cache%grhoc_t)+&
                & 4._dl*(mg_cache%grhor_t+mg_cache%grhog_t)+3._dl*(mg_cache%grhonu_t + mg_cache%gpresnu_t)) 

            ! fixed for w_DE /= -1
            !f1=mg_cache%k2+3._dl*( mg_cache%adotoa**2 - mg_cache%Hdot )  

            f1 = mg_cache%k2+0.5d0*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) + 3._dl*(mg_cache%grhonu_t+mg_cache%gpresnu_t) &
                & + 3._dl*(mg_cache%grhov_t+mg_cache%gpresv_t))

            k2alpha= mg_cache%k * mg_cache%sigma

            term1 = mg_cache%q * f1 * mg_cache%dgq/mg_cache%k

            term2 = k2alpha * ((mg_cache%q - 1._dl) * ( mg_cache%grhob_t+mg_cache%grhoc_t+(4._dl/3._dl) &
                    & *(mg_cache%grhor_t+mg_cache%grhog_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t) &
                    & ) -mg_cache%grhov_t - mg_cache%gpresv_t)

            !term2 = k2alpha * ((mg_cache%q) * ( mg_cache%grhob_t+mg_cache%grhoc_t+(4._dl/3._dl) &
            !    & *(mg_cache%grhor_t+mg_cache%grhog_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t)) &
            !    & - 2._dl *(mg_cache%adotoa**2 - mg_cache%Hdot))

            term3 = -( mg_cache%qdot + (mg_cache%r-1._dl) * mg_cache%q * mg_cache%adotoa ) * mg_cache%rhoDelta

            mg_cache%etadot = (term1 + term2 + term3)/( 2._dl * fQ )

            mg_cache%z = mg_cache%sigma - 3._dl * mg_cache%etadot/mg_cache%k

            !calculating also ISW related quantities
            mg_cache%MG_psi     = mg_cache%r * mg_cache%MG_phi - mg_cache%q * 1._dl * mg_cache%dgpi/mg_cache%k2
            mg_cache%MG_phidot  = mg_cache%etadot - mg_cache%adotoa * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha) &
                                & - mg_cache%Hdot * mg_cache%MG_alpha


        else if ( MG_flag == 4 )   then 

            if (CDM_flag == 1) then 

                fs = mg_cache%k2+0.5d0*(3._dl*(mg_cache%grhob_t) + &
                    & 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) + 3._dl * (mg_cache%grhonu_t + mg_cache%gpresnu_t) + &
                    & 3._dl*((1._dl + mg_cache%C_phi)*mg_cache%grhoc_t))

                f1 = mg_cache%k2+3._dl*( mg_cache%adotoa**2 - mg_cache%Hdot )

                term1 = f1*mg_cache%dgq/mg_cache%k

                term2 = -mg_cache%rhoDeltac*mg_cache%C_phidot

                term3 = mg_cache%k2*mg_cache%MG_alpha*((mg_cache%C_phi)*mg_cache%grhoc_t - &
                      & (mg_cache%grhov_t + mg_cache%gpresv_t))

                term4 = - (1._dl+mg_cache%C_phi)* (-1._dl)*(beta*betadot+3._dl*mg_cache%adotoa*beta**2) &
                    *(mg_cache%grhoc_t*mg_cache%dgrhoc - 3._dl*mg_cache%grhoc_t**2*mg_cache%MG_alpha*mg_cache%adotoa) &
                    /(mg_cache%k2 + m**2*a**2)    


                mg_cache%etadot = (term1 + term2 + term3 + term4)/( 2._dl * fs)

                mg_cache%z = mg_cache%sigma - 3._dl * mg_cache%etadot/mg_cache%k

                
                mg_cache%MG_psi = -(mg_cache%rhoDelta + (mg_cache%C_phi)*mg_cache%rhoDeltac + 2._dl* mg_cache%dgpi) &
                     & /(2._dl*mg_cache%k2)


                mg_cache%MG_phi = - mg_cache%MG_psi-(mg_cache%rhoDelta + (mg_cache%C_phi)*mg_cache%rhoDeltac+mg_cache%dgpi) &
                        & /mg_cache%k2


                mg_cache%MG_phidot = mg_cache%etadot - mg_cache%adotoa * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha) &
                                    & - mg_cache%Hdot * mg_cache%MG_alpha 
          
            end if      

        else if ( MG_flag == 5 )  then

			if(muSigma_flag ==1 .and. MGDE_pert) then
				fs = mg_cache%k2+0.5d0*(2._dl*mg_cache%BigSigma - mg_cache%mu)*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
					& + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) +3._dl * (mg_cache%grhonu_t + mg_cache%gpresnu_t ) &
					& + 3._dl * (mg_cache%grhov_t + mg_cache%gpresv_t )) 
			else
                fs = mg_cache%k2+0.5d0*(2._dl*mg_cache%BigSigma - mg_cache%mu)*(3._dl*(mg_cache%grhoc_t+mg_cache%grhob_t) &
                    & + 4._dl*(mg_cache%grhog_t+mg_cache%grhor_t) +3._dl * (mg_cache%grhonu_t + mg_cache%gpresnu_t )) 
            end if

            f1 = mg_cache%k2+3._dl*( mg_cache%adotoa**2 - mg_cache%Hdot )  

            term1 = (2._dl*mg_cache%BigSigma - mg_cache%mu)*f1*mg_cache%dgq/mg_cache%k

            term2 = mg_cache%rhoDelta*(2._dl*mg_cache%adotoa*(mg_cache%BigSigma - mg_cache%mu) &
                     - (2._dl*mg_cache%BigSigmadot - mg_cache%mudot))

			if(muSigma_flag ==1 .and. MGDE_pert) then
				term3 = mg_cache%k2*mg_cache%MG_alpha*((2._dl*mg_cache%BigSigma - mg_cache%mu)*( mg_cache%grhoc_t+mg_cache%grhob_t &
						& +(4._dl/3._dl)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t) &
						& + (mg_cache%grhov_t + mg_cache%gpresv_t))- 2._dl*(mg_cache%adotoa**2 - mg_cache%Hdot))
			else
                term3 = mg_cache%k2*mg_cache%MG_alpha*((2._dl*mg_cache%BigSigma - mg_cache%mu)*( mg_cache%grhoc_t+mg_cache%grhob_t &
                        & +(4._dl/3._dl)*(mg_cache%grhog_t+mg_cache%grhor_t) + (mg_cache%grhonu_t + mg_cache%gpresnu_t)) &
                        & - 2._dl*(mg_cache%adotoa**2 - mg_cache%Hdot))
            end if

            term4 =  - 2._dl*(mg_cache%BigSigma - mg_cache%mu)*mg_cache%pidot_sum

            term5 = 2._dl*mg_cache%adotoa*(mg_cache%BigSigma - mg_cache%mu)*(mg_cache%dgpi_w_sum + mg_cache%dgpi) &
                    & - 2._dl*(mg_cache%BigSigmadot-mg_cache%mudot)*mg_cache%dgpi

            mg_cache%etadot = (term1 + term2 + term3 + term4 + term5)/( 2._dl * fs)

            mg_cache%z = mg_cache%sigma - 3._dl * mg_cache%etadot/mg_cache%k

            mg_cache%MG_psi = - mg_cache%mu * ( mg_cache%rhoDelta + 2._dl* mg_cache%dgpi)/(2._dl*mg_cache%k2)

            mg_cache%MG_phi = - mg_cache%MG_psi - mg_cache%BigSigma*(mg_cache%rhoDelta+mg_cache%dgpi)/mg_cache%k2
                  

            mg_cache%MG_phidot = mg_cache%etadot - mg_cache%adotoa * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha) &
                                & - mg_cache%Hdot * mg_cache%MG_alpha             

        end if
        

        ! calculate sigmadot
        mg_cache%sigmadot = mg_cache%k * (mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha) 

    end subroutine MGCAMB_compute_z     

    !------------------------------------------------------------GaugeInterface_EvolveScal---------------
    !> this subroutine computes the ISW term in MG 
    subroutine MGCAMB_compute_ISW( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        !local variables
        real(dl) :: term0
     
        real(dl) :: m,beta,betadot 

        m = MGCAMB_M(a, mg_par_cache, mg_cache)
        beta = MGCAMB_beta(a, mg_par_cache, mg_cache)
        betadot = MGCAMB_betadot(a, mg_par_cache, mg_cache)

        term0 = mg_cache%k2 + 3._dl* (mg_cache%adotoa**2._dl - mg_cache%Hdot)

        !adding MG_rhoDeltadot
        mg_cache%rhoDeltadot = -term0 * mg_cache%dgq/mg_cache%k - (mg_cache%grho + mg_cache%gpres)* mg_cache%k*mg_cache%z &
                            & - mg_cache%adotoa * mg_cache%rhoDelta - 2._dl * mg_cache%adotoa * mg_cache%dgpi 

        !adding dgpidot
        mg_cache%dgpidot = mg_cache%pidot_sum - (2._dl*mg_cache%dgpi+ mg_cache%dgpi_diff )*mg_cache%adotoa  

        if (( MG_flag == 1 .and. pure_MG_flag /= 3 ) & ! all the mu, gamma parametrizations
            .or. MG_flag == 2 &
            .or. MG_flag == 3 .or. MG_flag == 6) then

            !! derived from 1st Possion eq
            mg_cache%MG_psidot = - 0.5d0*mg_cache%mu/mg_cache%k2*(mg_cache%rhoDeltadot+2._dl*mg_cache%dgpidot) &
                                & - 0.5d0*mg_cache%mudot/mg_cache%k2*(mg_cache%rhoDelta+2._dl*mg_cache%dgpi) 

        else if (  MG_flag == 1 .and. pure_MG_flag == 3  ) then

            mg_cache%MG_psidot = mg_cache%R * mg_cache%MG_phidot + mg_cache%Rdot * mg_cache%MG_phi - &
                            & mg_cache%Qdot*mg_cache%dgpi/mg_cache%k2 - mg_cache%Q * mg_cache%dgpidot /mg_cache%k2

        else if( MG_flag == 4) then ! only-CDM coupling

            if(CDM_flag == 1) then ! CDM QSA

                mg_cache%rhoDeltadot = -term0 * mg_cache%dgq/mg_cache%k - (mg_cache%grho + mg_cache%gpres)* mg_cache%k*mg_cache%z &
                                & - mg_cache%adotoa * mg_cache%rhoDelta - 2._dl * mg_cache%adotoa * mg_cache%dgpi - &
                                (beta*betadot+3._dl*mg_cache%adotoa*beta**2) &
                                *(mg_cache%grhoc_t*mg_cache%dgrhoc - 3._dl*mg_cache%grhoc_t**2*mg_cache%MG_alpha*mg_cache%adotoa) &
                                /(mg_cache%k2 + m**2*a**2)

                mg_cache%rhoDeltacdot = -term0 * mg_cache%dgqc/mg_cache%k- mg_cache%grhoc_t*mg_cache%k*mg_cache%z & 
                             - mg_cache%adotoa * mg_cache%rhoDeltac - (beta*betadot+3._dl*mg_cache%adotoa*beta**2) &
                    *(mg_cache%grhoc_t*mg_cache%dgrhoc - 3._dl*mg_cache%grhoc_t**2*mg_cache%MG_alpha*mg_cache%adotoa) &
                    /(mg_cache%k2 + m**2*a**2)


                mg_cache%MG_psidot = - 0.5d0/mg_cache%k2*(mg_cache%rhoDeltadot + (mg_cache%C_phi)*mg_cache%rhoDeltacdot &
                             & + mg_cache%C_phidot*mg_cache%rhoDeltac + 2._dl*mg_cache%dgpidot)   


            end if 

                  
        else if( MG_flag == 5) then 

            mg_cache%MG_psidot = - 0.5d0*mg_cache%mu/mg_cache%k2*(mg_cache%rhoDeltadot+2._dl*mg_cache%dgpidot) &
                                & - 0.5d0*mg_cache%mudot/mg_cache%k2*(mg_cache%rhoDelta+2._dl*mg_cache%dgpi)

        end if
   

        mg_cache%MG_ISW = mg_cache%MG_phidot+mg_cache%MG_psidot

        mg_cache%MG_alphadot = mg_cache%MG_psi - mg_cache%adotoa * mg_cache%MG_alpha

    end subroutine MGCAMB_compute_ISW

    !---------------------------------------------------------------------------
    !> this subroutine computes the lensing term in MG
    subroutine MGCAMB_compute_lensing( a, mg_par_cache, mg_cache )
        use precision
        implicit none

        real(dl) :: a   !< scale factor
        type(MGCAMB_timestep_cache), intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        mg_cache%MG_lensing = mg_cache%MG_phi + mg_cache%MG_psi

    end subroutine MGCAMB_compute_lensing

    !-----------------------------------------------
    !> mu(a,k) function
    function MGCAMB_Mu( a, mg_par_cache, mg_cache )
        !use ModelParams
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Mu                                       !< MG mu function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! local variables
        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: t1, t2, t1dot, t2dot
        real(dl) :: omm, ommdot

        real(dl) :: omegaDE_t

        ! beta, m parametrization
        real(dl) :: beta, m

        !> pure MG models 
        if ( MG_flag == 1 .and. pure_MG_flag /= 3 ) then  ! generic mu-gamma parametrization

            if ( pure_MG_flag == 1 ) then ! mu-gamma

                if ( mugamma_par == 1 ) then ! BZ parametrization 
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                    MGCAMB_Mu = (1._dl + B1 * LKA1)/(1._dl + LKA1)  

                else if ( mugamma_par == 2 ) then ! Planck parametrization

                    ! changing the following
                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2

                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    MGCAMB_Mu = 1._dl + E11*omegaDE_t

                else if ( mugamma_par == 3 ) then ! effective Newton constant
                    MGCAMB_Mu = 1._dl+ga*(1._dl)**nn - ga*(1._dl)**(2._dl*nn)

                else if ( mugamma_par == 4 ) then 
                    MGCAMB_Mu = 1._dl

                end if

            else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                if ( muSigma_par == 1 ) then ! DES parametrization

                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    !MGCAMB_Mu = 1._dl + mu0 * omegaDE_t/mg_par_cache%omegav

                    ! this is being changed
                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    MGCAMB_Mu = 1._dl + mu0 * omegaDE_t/mg_par_cache%omegav


                else if ( muSigma_par == 2 ) then
                    MGCAMB_Mu = 1._dl

                end if

            end if

        !> alternative MG
        else if ( MG_flag == 2 ) then

            if (alt_MG_flag == 1) then !(Linder Gamma)
                omm=(mg_par_cache%omegab+mg_par_cache%omegac)/((mg_par_cache%omegab+mg_par_cache%omegac) &
                & + (1-mg_par_cache%omegab-mg_par_cache%omegac)*a**3)
                ommdot=-3._dl*omm**2*a**3*mg_cache%adotoa*(1-mg_par_cache%omegab-mg_par_cache%omegac) &
                & /(mg_par_cache%omegab+mg_par_cache%omegac)

                MGCAMB_Mu=2._dl/3._dl*omm**(Linder_gamma-1._dl)*&
                (omm**Linder_gamma+2-3._dl*Linder_gamma+3._dl*(Linder_gamma-0.5d0)*omm)

            else if ( alt_MG_flag == 2 ) then
                MGCAMB_Mu = 1._dl
            end if


        !> QSA models
        else if ( MG_flag == 3 ) then

            if ( QSA_flag == 1 ) then ! f(R)
                LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                LKA2 = lambda2_2 * mg_cache%k2 * a**ss
                MGCAMB_Mu = (1._dl + B1 * LKA1)/(1._dl + LKA1)
                MGCAMB_Mu = MGCAMB_Mu/(1._dl - 1.4d-8 * lambda1_2 * a**3)

            else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                      QSA_flag == 3 .or. &
                      QSA_flag == 4 ) then
                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                t1      = (2._dl*beta**2._dl)*mg_cache%k2
                t2      = (m**2._dl)*a**2._dl

                MGCAMB_Mu = (mg_cache%k2 + t1 + t2)/(mg_cache%k2 + t2)
                


            else if ( QSA_flag == 5 )  then
                MGCAMB_Mu = 1._dl

            end if


        else if (MG_flag == 5) then  !direct mu-Sigma parametrization

            if(muSigma_flag == 1) then !pure MG models

                if ( pure_MG_flag == 1 ) then ! mu-gamma

                    if ( mugamma_par == 1 ) then ! BZ parametrization 
                        LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                        LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                        MGCAMB_Mu = (1._dl + B1 * LKA1)/(1._dl + LKA1)  

                    else if ( mugamma_par == 2 ) then ! Planck parametrization

                        ! changing the following
                        !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2

                        omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                        MGCAMB_Mu = 1._dl + E11*omegaDE_t

                    else if ( mugamma_par == 3 ) then ! effective Newton constant
                        MGCAMB_Mu = 1._dl+ga*(1._dl)**nn - ga*(1._dl)**(2._dl*nn)

                    else if ( mugamma_par == 4 ) then 
                        MGCAMB_Mu = 1._dl

                        end if

                else if ( pure_MG_flag == 2 ) then !mu-Sigma

                    if ( muSigma_par == 1 ) then ! DES parametrization

                        !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                        !MGCAMB_Mu = 1._dl + mu0 * omegaDE_t/mg_par_cache%omegav

                        ! this is being changed
                        omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                        MGCAMB_Mu = 1._dl + mu0 * omegaDE_t/mg_par_cache%omegav


                    else if ( muSigma_par == 2 ) then
                        MGCAMB_Mu = 1._dl

                    end if

                end if

            else if(muSigma_flag == 2) then  !alternative MG models

                if (alt_MG_flag == 1) then !(Linder Gamma)
                    omm=(mg_par_cache%omegab+mg_par_cache%omegac)/((mg_par_cache%omegab+mg_par_cache%omegac) &
                    & + (1-mg_par_cache%omegab-mg_par_cache%omegac)*a**3)
                    ommdot=-3._dl*omm**2*a**3*mg_cache%adotoa*(1-mg_par_cache%omegab-mg_par_cache%omegac) &
                    & /(mg_par_cache%omegab+mg_par_cache%omegac)

                    MGCAMB_Mu=2._dl/3._dl*omm**(Linder_gamma-1._dl)*&
                    (omm**Linder_gamma+2-3._dl*Linder_gamma+3._dl*(Linder_gamma-0.5d0)*omm)
                
                else if ( alt_MG_flag == 2 ) then
                    MGCAMB_Mu = 1._dl
                end if
            
            else if(muSigma_flag == 3) then  ! all-matter QSA models

                if ( QSA_flag == 1 ) then ! f(R)
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss
                    MGCAMB_Mu = (1._dl + B1 * LKA1)/(1._dl + LKA1)
                    MGCAMB_Mu = MGCAMB_Mu/(1._dl - 1.4d-8 * lambda1_2 * a**3)

                else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                        QSA_flag == 3 .or. &
                        QSA_flag == 4 ) then
                    beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                    m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                    t1      = (2._dl*beta**2._dl)*mg_cache%k2
                    t2      = (m**2._dl)*a**2._dl

                    MGCAMB_Mu = (mg_cache%k2 + t1 + t2)/(mg_cache%k2 + t2)
                    


                else if ( QSA_flag == 5 )  then
                    MGCAMB_Mu = 1._dl

                end if
            
            else if(muSigma_flag == 4) then ! reconstruction

                call splint1(a_arr,mu_arr,ddmu_arr,2*nnode,a,MGCAMB_Mu)

            else if(muSigma_flag == 5) then 

                write(*,*) 'Please write your own mu function for another model'  
                stop

            end if 

 ! =============MGXrecon=============           
        else if (MG_flag == 6) then	 !reconstruction
		     call splint1(a_arr,mu_arr,ddmu_arr,2*nnode,a,MGCAMB_Mu)
 ! =============MGXrecon=============

        end if
       
    end function MGCAMB_Mu        
           
    !-----------------------------------------------
    !> \dot{mu}(a,k) function
    function MGCAMB_Mudot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Mudot                                    !< MG mudot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! local variables
        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: t1,t2,t1dot,t2dot
        real(dl) :: omm, ommdot


        ! mapping beta,m into mu,gamma
        real(dl) :: beta, betadot, m, mdot
        real(dl) :: mu

        real(dl) :: omegaDEdot

        !> pure MG models
        if ( MG_flag == 1 .and. pure_MG_flag /= 3 ) then

            if ( pure_MG_flag == 1 ) then ! mu-gamma
                if ( mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                    MGCAMB_Mudot = ((B1 - 1._dl) * mg_cache%adotoa * ss * LKA1) / ((1._dl+LKA1)**2._dl)

                else if ( mugamma_par == 2 ) then ! Planck parametrization

                    ! changingh the following quantity
                    !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                    !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2

                    omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                            & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                    MGCAMB_Mudot = E11*omegaDEdot

                else if ( mugamma_par == 3 ) then  ! Newton's constants
                    MGCAMB_Mudot = mg_cache%adotoa*a*ga*nn*(-1._dl+2._dl*(1._dl-a)**nn)*(1._dl-a)**(nn-1._dl)

                else if ( mugamma_par == 4 ) then
                    MGCAMB_Mudot = 0._dl

                end if

            else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                if ( muSigma_par == 1 ) then ! DES parametrization
                    ! changing the following
                    !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                    !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                    MGCAMB_Mudot =  mu0 * omegaDEdot/mg_par_cache%omegav

                else if ( muSigma_par == 2 ) then
                    MGCAMB_Mudot = 0._dl

                end if

            end if

        !> alternative MG
        else if ( MG_flag == 2 ) then

            if (alt_MG_flag == 1) then !(Linder Gamma)
                mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )

                omm=(mg_par_cache%omegab+mg_par_cache%omegac)/((mg_par_cache%omegab+mg_par_cache%omegac) &
                    & +(1-mg_par_cache%omegab-mg_par_cache%omegac)*a**3)
                ommdot=-3._dl*omm**2*a**3*mg_cache%adotoa*(1-mg_par_cache%omegab-mg_par_cache%omegac) &
                    & /(mg_par_cache%omegab+mg_par_cache%omegac)

                MGCAMB_Mudot = mu/omm*(Linder_gamma-1._dl)*ommdot+&
                    2._dl/3._dl*omm**(Linder_gamma-1._dl)*ommdot*&
                    (Linder_gamma*omm**(Linder_gamma-1._dl)+3._dl*(Linder_gamma-0.5d0))

            else if ( alt_MG_flag == 2 ) then
                MGCAMB_Mudot = 0._dl
            end if


        !> QSA models
        else if ( MG_flag == 3 ) then

            if ( QSA_flag == 1 ) then ! f(R)
                LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                MGCAMB_Mudot = ((B1 - 1._dl) * mg_cache%adotoa * ss * LKA1) / ((1._dl+LKA1)**2._dl)
                mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )
                MGCAMB_Mudot = MGCAMB_Mudot/(1._dl - 1.4d-8 * lambda1_2 * a**3) + 3._dl * &
                                mu* mg_cache%adotoa *a**3 *(1.4d-8 * lambda1_2 ) &
                                /(1._dl - 1.4d-8 * lambda1_2 * a**3)

            else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                    QSA_flag == 3 .or. &
                    QSA_flag == 4 ) then

                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
                mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

                t1 = (2._dl*beta**2._dl)*mg_cache%k2
                t2 = (m**2._dl)*a**2._dl
                t1dot = 4._dl*beta*betadot*mg_cache%k2
                t2dot = (2._dl*a**2._dl)*(m*mdot+ (m**2._dl)*mg_cache%adotoa)

                MGCAMB_Mudot = (t1dot*(mg_cache%k2 + t2) - t1*t2dot)/((mg_cache%k2 + t2)**2._dl)


            else if ( QSA_flag == 5 )  then
                MGCAMB_Mudot = 0._dl

            end if


        else if( MG_flag == 5) then ! direct mu-Sigma parametrization

            if(muSigma_flag == 1) then !pure MG models

                if ( pure_MG_flag == 1 ) then ! mu-gamma
                    if ( mugamma_par == 1 ) then ! BZ parametrization
                        LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                        LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                        MGCAMB_Mudot = ((B1 - 1._dl) * mg_cache%adotoa * ss * LKA1) / ((1._dl+LKA1)**2._dl)

                    else if ( mugamma_par == 2 ) then ! Planck parametrization

                        ! changingh the following quantity
                        !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                        !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2

                        omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                        MGCAMB_Mudot = E11*omegaDEdot

                    else if ( mugamma_par == 3 ) then  ! Newton's constants
                        MGCAMB_Mudot = mg_cache%adotoa*a*ga*nn*(-1._dl+2._dl*(1._dl-a)**nn)*(1._dl-a)**(nn-1._dl)

                    else if ( mugamma_par == 4 ) then
                        MGCAMB_Mudot = 0._dl

                    end if

                else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                    if ( muSigma_par == 1 ) then ! DES parametrization
                        ! changing the following
                        !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                        !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                        omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                    & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                        MGCAMB_Mudot =  mu0 * omegaDEdot/mg_par_cache%omegav

                    else if ( muSigma_par == 2 ) then
                        MGCAMB_Mudot = 0._dl

                    end if

                end if
 
            else if(muSigma_flag == 2) then !alternative MG models

                if (alt_MG_flag == 1) then !(Linder Gamma)
                    mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )

                    omm=(mg_par_cache%omegab+mg_par_cache%omegac)/((mg_par_cache%omegab+mg_par_cache%omegac) &
                        & +(1-mg_par_cache%omegab-mg_par_cache%omegac)*a**3)
                    ommdot=-3._dl*omm**2*a**3*mg_cache%adotoa*(1-mg_par_cache%omegab-mg_par_cache%omegac) &
                        & /(mg_par_cache%omegab+mg_par_cache%omegac)

                    MGCAMB_Mudot = mu/omm*(Linder_gamma-1._dl)*ommdot+&
                        2._dl/3._dl*omm**(Linder_gamma-1._dl)*ommdot*&
                        (Linder_gamma*omm**(Linder_gamma-1._dl)+3._dl*(Linder_gamma-0.5d0))

                else if ( alt_MG_flag == 2 ) then
                    MGCAMB_Mudot = 0._dl
                end if 

            else if(muSigma_flag == 3) then  ! all-matter QSA models

                if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                      QSA_flag == 3 .or. &
                      QSA_flag == 4 ) then

                    beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                    m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                    betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
                    mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

                    t1 = (2._dl*beta**2._dl)*mg_cache%k2
                    t2 = (m**2._dl)*a**2._dl
                    t1dot = 4._dl*beta*betadot*mg_cache%k2
                    t2dot = (2._dl*a**2._dl)*(m*mdot+ (m**2._dl)*mg_cache%adotoa)

                    MGCAMB_Mudot = (t1dot*(mg_cache%k2 + t2) - t1*t2dot)/((mg_cache%k2 + t2)**2._dl)

                else
                    write(*,*) 'Please refer to params_MG.ini and choose the correct QSA flag'
                    stop  

                end if 

            else if(muSigma_flag == 4) then !reconstruction

                call splint1(a_arr,dmu_arr,dddmu_arr,2*nnode,a,MGCAMB_Mudot)

                MGCAMB_Mudot = MGCAMB_Mudot*mg_cache%adotoa *a

            else if(muSigma_flag == 5) then
        
                write(*,*) 'Please write your own mudot function for another model'  
                stop
			end if

! =============MGXrecon=============
		else if (MG_flag == 6) then !reconstruction
			call splint1(a_arr,dmu_arr,dddmu_arr,2*nnode,a,MGCAMB_Mudot)

			MGCAMB_Mudot = MGCAMB_Mudot*mg_cache%adotoa *a
            
        end if
! =============MGXrecon=============

    end function MGCAMB_Mudot    

    !-----------------------------------------------
    ! gamma(a,k) function
    function MGCAMB_Gamma( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Gamma                                    !< MG gamma function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: t1,t2, t1dot, t2dot

        real(dl) :: beta, m
        real(dl) :: omegaDE_t

        real(dl) :: sigma_t 
        real(dl) :: mu_t

        !> pure MG models
        if ( MG_flag == 1 .and. pure_MG_flag /= 3 ) then

            if ( pure_MG_flag == 1 ) then ! mu-gamma
                if ( mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                    MGCAMB_Gamma = (1._dl + B2 * LKA2)/(1._dl +LKA2)

                else if ( mugamma_par == 2 ) then ! Planck parametrization
                    ! changing the following
                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    MGCAMB_Gamma = 1._dl+E22*omegaDE_t

                else if ( mugamma_par == 3 ) then
                    MGCAMB_Gamma = 1._dl

                else if ( mugamma_par == 4 ) then
                    MGCAMB_Gamma = 1._dl
                end if

            else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                if ( muSigma_par == 1 ) then ! DES parametrization
                    ! changing the following
                    !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    sigma_t = 1._dl + sigma0 * omegaDE_t / mg_par_cache%omegav
                    mu_t    = 1._dl + mu0 * omegaDE_t / mg_par_cache%omegav
                    MGCAMB_Gamma = 2._dl * sigma_t / mu_t - 1._dl

                else if ( muSigma_par == 2 ) then
                    MGCAMB_Gamma = 1._dl

                end if

            end if

        !> alternative MG
        else if ( MG_flag == 2 ) then

            if (alt_MG_flag == 1) then !(Linder Gamma)
                MGCAMB_Gamma = 1._dl

            else if ( alt_MG_flag == 2 ) then
                MGCAMB_Gamma = 1._dl
            end if


        !> QSA models
        else if ( MG_flag == 3) then

            if ( QSA_flag == 1 ) then ! f(R)
                LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                MGCAMB_Gamma = (1._dl + B2 * LKA2)/(1._dl +LKA2)

            else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                    QSA_flag == 3 .or. &
                    QSA_flag == 4 ) then

                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )

                t1 = (2._dl*beta**2._dl)*mg_cache%k2
                t2 = (m**2._dl)*a**2._dl

                MGCAMB_Gamma = (mg_cache%k2 - t1 + t2)/(mg_cache%k2 + t1 + t2)


            else if ( QSA_flag == 5 )  then
                MGCAMB_Gamma = 1._dl

            end if

        else if (MG_flag == 5) then 

            if(muSigma_flag == 1) then 

                if ( pure_MG_flag == 1 ) then ! mu-gamma
                    if ( mugamma_par == 1 ) then ! BZ parametrization
                        LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                        LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                        MGCAMB_Gamma = (1._dl + B2 * LKA2)/(1._dl +LKA2)

                    else if ( mugamma_par == 2 ) then ! Planck parametrization
                        ! changing the following
                        !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                        omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                        MGCAMB_Gamma = 1._dl+E22*omegaDE_t

                    else if ( mugamma_par == 3 ) then
                        MGCAMB_Gamma = 1._dl

                    else if ( mugamma_par == 4 ) then
                        MGCAMB_Gamma = 1._dl
                    end if

                else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                    if ( muSigma_par == 1 ) then ! DES parametrization
                        ! changing the following
                        !omegaDE_t = mg_cache%grhov_t / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                        omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                        sigma_t = 1._dl + sigma0 * omegaDE_t / mg_par_cache%omegav
                        mu_t    = 1._dl + mu0 * omegaDE_t / mg_par_cache%omegav
                        MGCAMB_Gamma = 2._dl * sigma_t / mu_t - 1._dl

                    else if ( muSigma_par == 2 ) then
                        MGCAMB_Gamma = 1._dl

                    end if

                end if

            else if(muSigma_flag == 2) then  

                if (alt_MG_flag == 1) then !(Linder Gamma)
                    MGCAMB_Gamma = 1._dl

                else if ( alt_MG_flag == 2 ) then
                    MGCAMB_Gamma = 1._dl
                end if  

            else if(muSigma_flag == 3) then

                if ( QSA_flag == 1 ) then ! f(R)
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                    MGCAMB_Gamma = (1._dl + B2 * LKA2)/(1._dl +LKA2)

                else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                        QSA_flag == 3 .or. &
                        QSA_flag == 4 ) then

                    beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                    m       = MGCAMB_M( a, mg_par_cache, mg_cache )

                    t1 = (2._dl*beta**2._dl)*mg_cache%k2
                    t2 = (m**2._dl)*a**2._dl

                    MGCAMB_Gamma = (mg_cache%k2 - t1 + t2)/(mg_cache%k2 + t1 + t2)


                else if ( QSA_flag == 5 )  then
                    MGCAMB_Gamma = 1._dl

                end if

            else if(muSigma_flag == 4) then
                call splint1(a_arr,gamma_arr,ddgamma_arr,2*nnode,a,MGCAMB_Gamma)
            
			end if

 ! =============MGXrecon=============
		else if (MG_flag == 6) then

			call splint1(a_arr,gamma_arr,ddgamma_arr,2*nnode,a,MGCAMB_Gamma)

        end if
! =============MGXrecon=============


    end function MGCAMB_Gamma


    !-----------------------------------------------
    ! \dot{gamma}(a,k) function
    function MGCAMB_Gammadot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Gammadot                                 !< MG gammadot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters


        real(dl) :: LKA1 ! \lambda_1^2 k^2 a^s
        real(dl) :: LKA2 ! \lambda_1^2 k^2 a^s
        real(dl) :: t1,t2,t1dot,t2dot

        real(dl) :: beta, betadot, m, mdot
        real(dl) :: omegaDE_t, omegaDEdot
        real(dl) :: sigma_t, sigmadot_t  
        real(dl) :: mu_t, mudot_t



        !> pure MG models
        if ( MG_flag == 1 .and. pure_MG_flag /= 3 ) then

            if ( pure_MG_flag == 1 ) then ! mu-gamma
                if ( mugamma_par == 1 ) then ! BZ parametrization
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                    MGCAMB_Gammadot = ((B2 -1._dl)*mg_cache%adotoa * ss* LKA2)/((1._dl+LKA2)**2._dl)

                else if ( mugamma_par == 2 ) then ! Planck parametrization
                    ! changing the following
                    !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                    !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                    omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                    MGCAMB_Gammadot = E22*omegaDEdot

                else if ( mugamma_par == 3 ) then
                    MGCAMB_Gammadot = 0._dl

                else if ( mugamma_par == 4 ) then
                    MGCAMB_Gammadot = 0._dl

                end if

            else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                if ( muSigma_par == 1 ) then ! DES parametrization

                ! changing the following
                omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                omegaDEdot  =-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t
                sigma_t     = 1._dl + sigma0 * omegaDE_t / mg_par_cache%omegav
                sigmadot_t  = sigma0 * omegaDEdot / mg_par_cache%omegav
                mu_t        = 1._dl + mu0 * omegaDE_t / mg_par_cache%omegav
                mudot_t     = mu0 * omegaDEdot / mg_par_cache%omegav
                MGCAMB_Gammadot = 2._dl * sigmadot_t / mu_t - 2._dl *sigma_t*mudot_t/mu_t**2

                else if ( muSigma_par == 2 ) then
                    MGCAMB_Gammadot = 0._dl

                end if

            end if

        !> alternative MG
        else if ( MG_flag == 2 ) then

            if (alt_MG_flag == 1) then !(Linder Gamma)
                MGCAMB_Gammadot = 0._dl

            else if ( alt_MG_flag == 2 ) then
                MGCAMB_Gammadot = 0._dl
            end if


        !> QSA models
        else if ( MG_flag == 3) then

            if ( QSA_flag == 1 ) then ! f(R)
                LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                MGCAMB_Gammadot = ((B2 -1._dl)*mg_cache%adotoa * ss* LKA2)/((1._dl+LKA2)**2._dl)

            else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                    QSA_flag == 3  .or. &
                    QSA_flag == 4 ) then
                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
                mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

                t1      = (2._dl*beta**2._dl)*mg_cache%k2
                t2      = (m**2._dl)*a**2._dl
                t1dot   = 4._dl*beta*betadot*mg_cache%k2
                t2dot   = (2._dl*a**2._dl)*(m*mdot + (m**2._dl) *mg_cache%adotoa)

                MGCAMB_Gammadot = 2._dl*(t1*t2dot-t1dot*(mg_cache%k2 + t2))/((mg_cache%k2 + t1 + t2)**2._dl)

            else if ( QSA_flag == 5 )  then
                MGCAMB_Gammadot = 0._dl

            end if
        
        else if(MG_flag == 5) then 
        
            if(muSigma_flag == 1) then 

                if ( pure_MG_flag == 1 ) then ! mu-gamma
                    if ( mugamma_par == 1 ) then ! BZ parametrization
                        LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                        LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                        MGCAMB_Gammadot = ((B2 -1._dl)*mg_cache%adotoa * ss* LKA2)/((1._dl+LKA2)**2._dl)

                    else if ( mugamma_par == 2 ) then ! Planck parametrization
                        ! changing the following
                        !omegaDEdot = - 3._dl * mg_cache%adotoa * (mg_cache%grhov_t + mg_cache%gpresv_t) &
                        !            & / a**2 / 3._dl / mg_par_cache%h0_Mpc**2
                        omegaDEdot=-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                    & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t

                        MGCAMB_Gammadot = E22*omegaDEdot

                    else if ( mugamma_par == 3 ) then
                        MGCAMB_Gammadot = 0._dl

                    else if ( mugamma_par == 4 ) then
                        MGCAMB_Gammadot = 0._dl

                    end if

                else if ( pure_MG_flag == 2 ) then ! mu-Sigma

                    if ( muSigma_par == 1 ) then ! DES parametrization

                    ! changing the following
                    omegaDE_t = mg_cache%grhov_t / 3._dl / mg_cache%adotoa**2
                    omegaDEdot  =-(mg_cache%grhov_t+3._dl*mg_cache%gpresv_t)/3._dl/mg_cache%adotoa &
                                    & - 2._dl*mg_cache%Hdot/3._dl/mg_cache%adotoa**3*mg_cache%grhov_t
                    sigma_t     = 1._dl + sigma0 * omegaDE_t / mg_par_cache%omegav
                    sigmadot_t  = sigma0 * omegaDEdot / mg_par_cache%omegav
                    mu_t        = 1._dl + mu0 * omegaDE_t / mg_par_cache%omegav
                    mudot_t     = mu0 * omegaDEdot / mg_par_cache%omegav
                    MGCAMB_Gammadot = 2._dl * sigmadot_t / mu_t - 2._dl *sigma_t*mudot_t/mu_t**2

                    else if ( muSigma_par == 2 ) then
                        MGCAMB_Gammadot = 0._dl

                    end if

                end if   

            else if(muSigma_flag == 2)  then
            
                if (alt_MG_flag == 1) then !(Linder Gamma)
                    MGCAMB_Gammadot = 0._dl

                else if ( alt_MG_flag == 2 ) then
                    MGCAMB_Gammadot = 0._dl
                end if
            
            else if(muSigma_flag == 3) then 

                if ( QSA_flag == 1 ) then ! f(R)
                    LKA1 = lambda1_2 * mg_cache%k2 * a**ss
                    LKA2 = lambda2_2 * mg_cache%k2 * a**ss

                    MGCAMB_Gammadot = ((B2 -1._dl)*mg_cache%adotoa * ss* LKA2)/((1._dl+LKA2)**2._dl)

                else if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                        QSA_flag == 3  .or. &
                        QSA_flag == 4 ) then
                    beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                    m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                    betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
                    mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

                    t1      = (2._dl*beta**2._dl)*mg_cache%k2
                    t2      = (m**2._dl)*a**2._dl
                    t1dot   = 4._dl*beta*betadot*mg_cache%k2
                    t2dot   = (2._dl*a**2._dl)*(m*mdot + (m**2._dl) *mg_cache%adotoa)

                    MGCAMB_Gammadot = 2._dl*(t1*t2dot-t1dot*(mg_cache%k2 + t2))/((mg_cache%k2 + t1 + t2)**2._dl)

                else if ( QSA_flag == 5 )  then
                    MGCAMB_Gammadot = 0._dl

                end if
            
            else if(muSigma_flag == 4) then 
                call splint1(a_arr,dgamma_arr,dddgamma_arr,2*nnode,a,MGCAMB_Gammadot)
                MGCAMB_Gammadot = MGCAMB_Gammadot*mg_cache%adotoa *a
            end if

! =============MGXrecon=============
		else if (MG_flag == 6) then
			call splint1(a_arr,dgamma_arr,dddgamma_arr,2*nnode,a,MGCAMB_Gammadot)
			MGCAMB_Gammadot = MGCAMB_Gammadot*mg_cache%adotoa *a

        end if

! =============MGXrecon=============

    end function MGCAMB_Gammadot

    ! Sigma(a,k) function
    function MGCAMB_BigSigma( a, mg_par_cache, mg_cache ) 

        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_BigSigma                                 !< MG BigSigma function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: mu, gamma

        if (MG_flag == 5) then
            mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )
            gamma = MGCAMB_Gamma( a, mg_par_cache, mg_cache )

            MGCAMB_BigSigma = mu*(1._dl+gamma)/2._dl
  
        end if       

    end function MGCAMB_BigSigma  


    !dot{Sigma}(a,k) function
    function MGCAMB_BigSigmadot( a, mg_par_cache, mg_cache )

        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_BigSigmadot                              !< MG BigSigmadot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: mu, gamma, mudot, gammadot


        if (MG_flag == 5) then
            mu = MGCAMB_Mu( a, mg_par_cache, mg_cache )
            gamma = MGCAMB_Gamma( a, mg_par_cache, mg_cache )
            mudot = MGCAMB_Mudot( a, mg_par_cache, mg_cache )
            gammadot = MGCAMB_Gammadot( a, mg_par_cache, mg_cache )

            MGCAMB_BigSigmadot = (mudot*(1._dl+gamma) + mu*gammadot)/2._dl

        end if 
     
    end function MGCAMB_BigSigmadot

    !C_phi(a,k) function
    function MGCAMB_C_phi( a, mg_par_cache, mg_cache ) 

        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_C_phi                                    !< MG C_phi function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! local variables
        real(dl) :: t1, t2

        real(dl) :: beta, m

        if ( MG_flag == 4 ) then ! only-CDM coupling

            if (CDM_flag == 1) then ! CDM QSA

                if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                      QSA_flag == 3 .or. &
                      QSA_flag == 4 ) then     


                beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                t1      = (beta**2._dl)*mg_cache%grhoc_t
                t2      = (m**2._dl)*a**2._dl
                              
                MGCAMB_C_phi  = t1/(t2 + mg_cache%k2)


                else 
                     write(*,*) 'Please refer to params_MG.ini and choose the correct QSA flag'
                     stop 

                end if 

            end if 

        end if 

    end function MGCAMB_C_phi  

    !dot{C_phi}(a,k) function
    function MGCAMB_C_phidot( a, mg_par_cache, mg_cache )

        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_C_phidot                                 !< MG C_phidot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! local variables
        real(dl) :: t1,t2,t1dot,t2dot

        real(dl) :: beta, betadot, m, mdot

     
        if ( MG_flag == 4 ) then !only-CDM coupling

            if (CDM_flag == 1) then  !CDM QSA

                if ( QSA_flag == 2 .or. &  ! beta, m parametrization
                     QSA_flag == 3 .or. &
                     QSA_flag == 4 ) then     

                    beta    = MGCAMB_Beta( a, mg_par_cache, mg_cache )
                    m       = MGCAMB_M( a, mg_par_cache, mg_cache )
                    betadot = MGCAMB_Betadot( a, mg_par_cache, mg_cache )
                    mdot    = MGCAMB_Mdot( a, mg_par_cache, mg_cache )

                    t1      = (beta**2._dl)*mg_cache%grhoc_t
                    t2      = (m**2._dl)*a**2._dl
                    t1dot   = 2._dl*beta*betadot*mg_cache%grhoc_t - beta**2*mg_cache%grhoc_t*mg_cache%adotoa  
                    t2dot   = (2._dl*a**2._dl)*(m*mdot+ (m**2._dl)*mg_cache%adotoa)

                    MGCAMB_C_phidot =  (t1dot*(mg_cache%k2 + t2) - t1*t2dot)/((mg_cache%k2 + t2)**2._dl)

                else 
                    write(*,*) 'Please refer to params_MG.ini and choose the correct QSA flag'
                    stop
                end if 

            end if 

        end if 

    end function MGCAMB_C_phidot 


!----------------------------------------------------------------------------------------------
!> MGCAMB (beta, m) parametrization, QSA for scalar-tensor models

    !-----------------------------------------------
    !> m(a) function
    function MGCAMB_M( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_M                                        !< MG m function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: FRm0

        ! SYMMETRON
        if( QSA_flag ==  2 ) then
            MGCAMB_M = (mg_par_cache%H0/3.0D05) / (xi_star) * sqrt(1._dl-(a_star/a)**3._dl)

        ! DILATON: based on 1206.3568
        else if ( QSA_flag ==  3 ) then
            MGCAMB_M = (mg_par_cache%H0/3.0D05) /(xi0) * a**(- DilR)

        ! Hu-Sawicki f(R) model: m, beta parametrization as in 1305.5647
        else if ( QSA_flag ==  4 )then
            FRm0 = (mg_par_cache%h0/3.0D05)*sqrt((4._dl*mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac) &
                    & /((FRn+1._dl)*F_R0))!note factor of c here
            MGCAMB_M = FRm0 * ((4._dl * mg_par_cache%omegav + (mg_par_cache%omegab + mg_par_cache%omegac)*a**(-3._dl)) &
                    & /(4._dl * mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac))**(FRn/2._dl+1._dl)


        end if

    end function MGCAMB_M

    !-----------------------------------------------
    !> \dot{m}(a) function
    function MGCAMB_Mdot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Mdot                                     !< MG mdot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: FRm0
        real(dl) :: m

        m = MGCAMB_M( a, mg_par_cache, mg_cache )

        ! SYMMETRON
        if( QSA_flag ==  2 ) then
            MGCAMB_Mdot = 1.5d0*(mg_par_cache%H0/3.0D05)/(xi_star)*((a_star/a)**3._dl*mg_cache%adotoa)/&
                        & (sqrt(1._dl-(a_star/a)**3._dl))

        ! DILATON
        else if ( QSA_flag ==  3 ) then
            MGCAMB_Mdot = - DilR * m * mg_cache%adotoa


        ! Hu-Sawicki f(R) model
        else if ( QSA_flag ==  4 )then

            FRm0 = (mg_par_cache%h0/3.0D05)*sqrt((4._dl*mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac)/ &
                    & ((FRn+1._dl)*F_R0))
            MGCAMB_Mdot = m / (4._dl * mg_par_cache%omegav + (mg_par_cache%omegab + mg_par_cache%omegac)*a**(-3._dl)) &
                    & * (-3._dl*FRn/2._dl-3._dl)*((mg_par_cache%omegab + mg_par_cache%omegac)* a**(-3._dl)*mg_cache%adotoa)!/(4._dl * mg_par_cache%omegav + mg_par_cache%omegab + mg_par_cache%omegac)) ! complete this

        end if

    end function MGCAMB_Mdot

    !-----------------------------------------------
    !> beta(a) function
    function MGCAMB_Beta( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Beta                                     !< MG beta function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        ! SYMMETRON
        if( QSA_flag == 2 ) then
            MGCAMB_Beta =  beta_star * sqrt(1._dl-(a_star/a)**3._dl)

        ! DILATON
        else if ( QSA_flag == 3 ) then
            MGCAMB_Beta = beta0 * exp((DilS)/(2._dl* DilR - 3._dl)*(a**(2._dl* DilR - 3._dl)-1._dl))

        ! Hu-Sawicki f(R) model
        else if ( QSA_flag == 4 )then
            MGCAMB_Beta = beta0

        end if

    end function MGCAMB_Beta

    !-----------------------------------------------
    !> \dot{beta}(a) function
    function MGCAMB_Betadot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Betadot                                  !< MG betadot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        real(dl) :: beta

        beta = MGCAMB_Beta( a, mg_par_cache, mg_cache )

        ! SYMMETRON
        if( QSA_flag == 2 ) then
            MGCAMB_Betadot = 1.5d0 * (beta_star * (a_star/a)**3._dl * mg_cache%adotoa) /( sqrt(1._dl-(a_star/a)**3._dl))

        ! DILATON
        else if ( QSA_flag == 3 ) then
            MGCAMB_Betadot = beta * (DilS * a**(2._dl* DilR - 3._dl) *  mg_cache%adotoa)

        ! Hu-Sawicki f(R) model
        else if ( QSA_flag == 4 )then
            MGCAMB_Betadot = 0._dl


        end if

    end function MGCAMB_Betadot

    !----------------------------------------------------------------------------------------------
    !> MGCAMB (Q,R) parametrization, QSA for scalar-tensor models

    !-----------------------------------------------
    !> Q(a,k) function
    function MGCAMB_Q( a, mg_par_cache, mg_cache )
    implicit none
    real(dl) :: a                                               !< scale factor
    real(dl) :: MGCAMB_Q                                        !< MG Q function
    type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
    type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        if ( QR_par == 1) then
            MGCAMB_Q = MGQfix

        else if ( QR_par == 2 ) then
            MGCAMB_Q = 1._dl + (Qnot - 1._dl)* a**sss

        end if

    end function MGCAMB_Q

    !-----------------------------------------------
    !> \dot{Q}(a,k) function
    function MGCAMB_Qdot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Qdot                                     !< MG Qdot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        if ( QR_par == 1 ) then
            MGCAMB_Qdot = 0._dl

        else if ( QR_par == 2 ) then
            MGCAMB_Qdot = (Qnot - 1._dl)*mg_cache%adotoa* sss* a**(sss)
        end if

    end function MGCAMB_Qdot

    !-----------------------------------------------
    ! R(a,k) function
    function MGCAMB_R( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_R                                        !< MG R function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters


        if ( QR_par == 1 ) then
            MGCAMB_R=MGRfix

        else if ( QR_par == 2 ) then
            MGCAMB_R = 1._dl + (Rnot - 1._dl)* a**sss

        end if

    end function MGCAMB_R

    !-----------------------------------------------
    ! \dot{R}(a,k) function
    function MGCAMB_Rdot( a, mg_par_cache, mg_cache )
        implicit none
        real(dl) :: a                                               !< scale factor
        real(dl) :: MGCAMB_Rdot                                     !< MG Rdot function
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache        !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in) :: mg_par_cache    !< cache containing the parameters

        if ( QR_par == 1 ) then
            MGCAMB_Rdot = 0._dl

        else if ( QR_par == 2 ) then
            MGCAMB_Rdot = (Rnot - 1._dl)*mg_cache%adotoa* sss* a**(sss)

        end if

    end function MGCAMB_Rdot

    ! ---------------------------------------------------------------------------------------------
    !> Modifying the background

    subroutine MGCAMB_DarkEnergy( a, mg_par_cache, mg_cache ) 
        use precision
        implicit none

        real(dl) :: a   !< scale factor 
        type(MGCAMB_timestep_cache),  intent(inout) :: mg_cache      !< cache containing the time-dependent quantities
        type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

        real(dl) :: wnow

		real(dl) :: z, X, dXda, Y

        if ( DE_model == 0 ) then 
            mg_cache%grhov_t = 3._dl*mg_par_cache%h0_Mpc**2 * mg_par_cache%omegav *a**2
            mg_cache%gpresv_t = - mg_cache%grhov_t

        else if ( DE_model == 1 ) then
			if (a > 1e-10) then
            	mg_cache%grhov_t = 3._dl*mg_par_cache%h0_Mpc**2*mg_par_cache%omegav*a**(-1._dl-3._dl*w0DE)
			else 
				mg_cache%grhov_t = 0._dl
			end if          
            mg_cache%gpresv_t = mg_cache%grhov_t * w0DE
            
        else if (DE_model == 2 ) then
            wnow = w0DE+(1._dl-a)*waDE
			if (a > 1e-10) then
            	mg_cache%grhov_t = 3._dl*mg_par_cache%h0_Mpc**2*mg_par_cache%omegav* &
	                         a**(-1.d0-3.d0*w0DE-3.d0*waDE)*Exp(3.d0*waDE*(a-1.d0))
			else
				mg_cache%grhov_t = 0._dl
			end if
            mg_cache%gpresv_t = mg_cache%grhov_t * wnow

        else if ( DE_model == 3 ) then

			call splint1(a_arr,X_arr,ddX_arr,2*nnode,a,X)
			call splint1(a_arr,dX_arr,dddX_arr,2*nnode,a,dXda)

			Y = -X - dXda*a/3.d0 ! pg 5 of 1807.03772
			mg_cache%grhov_t  = X*3.d0*mg_par_cache%h0_Mpc**2.d0*a**2.d0
			mg_cache%gpresv_t = Y*3.d0*mg_par_cache%h0_Mpc**2.d0*a**2.d0

        else if ( DE_model == 4 ) then
            write(*,*) 'This will contain the reconstruction of rho_DE(a)'
            write(*,*) 'Not implemented yet'
            stop
        else
            write(*,*) 'choose a DE model'
            stop
        end if

    end subroutine MGCAMB_DarkEnergy

	subroutine MGCAMB_DE_perturb
        use precision
        implicit none

        if ( DE_model == 0 ) then
			MGDE_const = .True.

		else if ( DE_model == 1 ) then
			if(w0DE == -1._dl) then
				MGDE_const = .True.
			else
				MGDE_const = .False.
			end if
		else if (DE_model == 2) then
			if(w0DE == -1._dl .and. waDE == 0._dl) then
				MGDE_const = .True.
			else
				MGDE_const = .False.
			end if
		end if

	end subroutine MGCAMB_DE_perturb

	subroutine MGCAMB_DE_EoS(a, w)
		use precision
		implicit none

        real(dl) :: a   !< scale factor
		real(dl) :: w

        w = -1._dl
        
		if ( DE_model == 0 ) then
			w = -1._dl
		else if(DE_model == 1) then
			w = w0DE

		else if(DE_model == 2) then
			w = w0DE+(1._dl-a)*waDE
		end if

	end subroutine MGCAMB_DE_EoS
    
	!>this subroutine reads MG parameters into the module
	subroutine MGCAMB_read_in_MGparams(CP)
		use precision
		use constants
		use model
        implicit none

		type(CAMBparams) :: CP

        MG_flag = CP%ModGravity%MG_flag
        pure_MG_flag = CP%ModGravity%pure_MG_flag
        alt_MG_flag = CP%ModGravity%alt_MG_flag
        QSA_flag = CP%ModGravity%QSA_flag
        mugamma_par  = CP%ModGravity%mugamma_par
        muSigma_par  = CP%ModGravity%muSigma_par
        QR_par = CP%ModGravity%QR_par
        muSigma_flag = CP%ModGravity%muSigma_flag
        CDM_flag  = CP%ModGravity%CDM_flag

        DE_model = CP%ModGravity%DE_model

        GRtrans = CP%ModGravity%GRtrans                

        ! BZ parametrization (and QS f(R))
        B1 =  CP%ModGravity%B1
        B2 =  CP%ModGravity%B2
        lambda1_2 =  CP%ModGravity%lambda1_2
        lambda2_2 = CP%ModGravity%lambda2_2
        ss = CP%ModGravity%ss

        ! Planck Parametrization
        E11 = CP%ModGravity%E11
        E22 = CP%ModGravity%E22

        ! Q-R parametrization 1
        MGQfix = CP%ModGravity%MGQfix
        MGRfix = CP%ModGravity%MGRfix

        ! Q-R parametrization 2
        Qnot = CP%ModGravity%Qnot
        Rnot = CP%ModGravity%Rnot
        sss = CP%ModGravity%sss

        ! Growth rate gamma
        Linder_gamma = CP%ModGravity%Linder_gamma

        ! Symmetron
        beta_star = CP%ModGravity%beta_star
        a_star  = CP%ModGravity%a_star
        xi_star = CP%ModGravity%xi_star

        ! Dilaton
        beta0 = CP%ModGravity%beta0
        xi0 = CP%ModGravity%xi0
        DilR = CP%ModGravity%DilR
        DilS = CP%ModGravity%DilS

        ! Hu-Sawicki f(R) gravity
        F_R0 = CP%ModGravity%F_R0
        FRn = CP%ModGravity%FRn

        ! DES parametrization
        mu0 = CP%ModGravity%mu0
        sigma0 = CP%ModGravity%sigma0

        ! effective Newton's constant 
        ga = CP%ModGravity%ga
        nn = CP%ModGravity%nn

        ! DE model parameters
        w0DE = CP%ModGravity%w0DE             !< w0 parameters for DE
        waDE = CP%ModGravity%waDE             !< waDE parameters for DE

        !DE pertubations
        MGDE_pert = CP%ModGravity%MGDE_pert

        !reconstruction model parameters
		mu_arr(10) =  CP%ModGravity%MGCAMB_Mu_idx_1
		mu_arr(11) =  CP%ModGravity%MGCAMB_Mu_idx_2
		mu_arr(12) =  CP%ModGravity%MGCAMB_Mu_idx_3
		mu_arr(13) =  CP%ModGravity%MGCAMB_Mu_idx_4
		mu_arr(14) =  CP%ModGravity%MGCAMB_Mu_idx_5
		mu_arr(15) =  CP%ModGravity%MGCAMB_Mu_idx_6
		mu_arr(16) =  CP%ModGravity%MGCAMB_Mu_idx_7
		mu_arr(17) =  CP%ModGravity%MGCAMB_Mu_idx_8
		mu_arr(18) =  CP%ModGravity%MGCAMB_Mu_idx_9
		mu_arr(19) =  CP%ModGravity%MGCAMB_Mu_idx_10
		mu_arr(20) =  CP%ModGravity%MGCAMB_Mu_idx_11
		Sigma_arr(10) =  CP%ModGravity%MGCAMB_Sigma_idx_1
		Sigma_arr(11) =  CP%ModGravity%MGCAMB_Sigma_idx_2
		Sigma_arr(12) =  CP%ModGravity%MGCAMB_Sigma_idx_3
		Sigma_arr(13) =  CP%ModGravity%MGCAMB_Sigma_idx_4
		Sigma_arr(14) =  CP%ModGravity%MGCAMB_Sigma_idx_5
		Sigma_arr(15) =  CP%ModGravity%MGCAMB_Sigma_idx_6
		Sigma_arr(16) =  CP%ModGravity%MGCAMB_Sigma_idx_7
		Sigma_arr(17) =  CP%ModGravity%MGCAMB_Sigma_idx_8
		Sigma_arr(18) =  CP%ModGravity%MGCAMB_Sigma_idx_9
		Sigma_arr(19) =  CP%ModGravity%MGCAMB_Sigma_idx_10
		Sigma_arr(20) =  CP%ModGravity%MGCAMB_Sigma_idx_11
		X_arr(10) =  CP%ModGravity%Funcofw_1
		X_arr(11) =  CP%ModGravity%Funcofw_2
		X_arr(12) =  CP%ModGravity%Funcofw_3
		X_arr(13) =  CP%ModGravity%Funcofw_4
		X_arr(14) =  CP%ModGravity%Funcofw_5
		X_arr(15) =  CP%ModGravity%Funcofw_6
		X_arr(16) =  CP%ModGravity%Funcofw_7
		X_arr(17) =  CP%ModGravity%Funcofw_8
		X_arr(18) =  CP%ModGravity%Funcofw_9
		X_arr(19) =  CP%ModGravity%Funcofw_10
		X_arr(20) =  CP%ModGravity%Funcofw_11

		mgcamb_par_cache%omegab = CP%ombh2/(CP%H0/100)**2
		mgcamb_par_cache%omegac = CP%omch2/(CP%H0/100)**2
		mgcamb_par_cache%h0     = CP%H0
		mgcamb_par_cache%h0_Mpc = CP%H0 * (1.d3/c) 

		X_arr(2*nnode) = mgcamb_par_cache%omegav

        if ( MG_flag == 0 ) then
            ! copy the MGCAMB values for w0DE and waDE to the DarkEnergy class
            CP%DarkEnergy%is_cosmological_constant = .true.
            CP%DarkEnergy%w_lam = -1_dl
            CP%DarkEnergy%wa = 0._dl
            CP%DarkEnergy%cs2_lam = 1_dl
            CP%DarkEnergy%grho_de = 0.
            CP%DarkEnergy%we = -1._dl

		else if ( MG_flag /= 0 ) then
			if ( MG_flag == 1 ) then
				if( pure_MG_flag /= 1 .and. pure_MG_flag /= 2 .and. pure_MG_flag /= 3) then
					stop 'Choose pure_MG_flag properly!'
				end if
				if ( DE_model /= 0 .and. DE_model /=1 .and. DE_model /= 2 .and. DE_model /= 3) then
					stop 'Choose DE_model properly!'
				end if

			else if ( MG_flag == 2 ) then
				if ( DE_model /= 0 ) then
					stop 'Choose DE_model properly!'
				end if

			else if ( MG_flag == 3) then
				if ( QSA_flag > 4 .or. QSA_flag < 1 ) then
					stop 'Choose QSA_flag properly!'
				end if
				if ( QSA_flag ==  1 ) then
					B1 = 4._dl/3._dl
					lambda1_2= CP%ModGravity%B0 ! it is considered as the B0 parameter here
					lambda1_2 = (lambda1_2*(299792458.d-3)**2)/(2._dl*mgcamb_par_cache%H0**2)
					B2 = 0.5d0
					lambda2_2 = B1* lambda1_2
					ss = 4._dl
				else if( QSA_flag ==  2 ) then
					GRtrans = a_star
				else if ( QSA_flag ==  4 ) then
					beta0 = 1._dl/sqrt(6._dl)
				end if

				if ( DE_model /= 0 ) then
					stop 'Choose DE_model properly!'
				end if

			else if ( MG_flag == 4) then
				if(CDM_flag == 1) then
					if ( QSA_flag > 4 .or. QSA_flag < 2 ) then
						stop 'Choose QSA_flag properly!'
					end if
					if( QSA_flag ==  2 ) then
						GRtrans = a_star
					else if ( QSA_flag ==  4 ) then
						beta0 = 1._dl/sqrt(6._dl)
					end if
				else
					stop 'Please choose CDM_flag properly'
				end if

				if ( DE_model /= 0 ) then
					stop 'Choose DE_model properly!'
				end if

			else if ( MG_flag == 5) then
				if( muSigma_flag== 1) then
					if( pure_MG_flag /= 1 .and. pure_MG_flag /= 2) then
						stop 'Choose pure_MG_flag properly!'
					end if
					if ( DE_model /= 0 .and. DE_model /=1 .and. DE_model /= 2 .and. DE_model /= 3) then
						stop 'Choose DE_model properly!'
					end if

				else if(muSigma_flag==2) then
					if ( DE_model /= 0 ) then
						stop 'Choose DE_model properly!'
					end if

				else if(muSigma_flag==3) then
					if ( QSA_flag > 4 .or. QSA_flag < 1 ) then
						 stop 'Choose QSA_flag properly!'
					end if
					if ( QSA_flag ==  1 ) then
						B1 = 4._dl/3._dl
						lambda1_2= CP%ModGravity%B0 ! it is considered as the B0 parameter here
						lambda1_2 = (lambda1_2*(299792458.d-3)**2)/(2._dl*mgcamb_par_cache%H0**2)
						B2 = 0.5d0
						lambda2_2 = B1* lambda1_2
						ss = 4._dl
					else if( QSA_flag ==  2 ) then
						GRtrans = a_star
					else if ( QSA_flag ==  4 ) then
						beta0 = 1._dl/sqrt(6._dl)
					end if

					if ( DE_model /= 0 ) then
						stop 'Choose DE_model properly!'
					end if

				else if(muSigma_flag==4) then
					if ( DE_model /=1 .and. DE_model /= 2 .and. DE_model /= 3) then
						stop 'Choose DE_model properly!'
					end if
				else
					stop 'Please write your own model'
				end if

			else if ( MG_flag == 6 ) then
				if ( DE_model /=1 .and. DE_model /= 2 .and. DE_model /= 3) then
					stop 'Choose DE_model properly!'
				end if
			else
				stop 'Choose MG_flag properly!'
			end if
		end if
    end subroutine MGCAMB_read_in_MGparams


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that reads the MGCAMB model parameters
    subroutine MGCAMB_read_model_params( mg_par_cache, Ini )
        use IniObjects
	    type(TIniFile) :: Ini   
        Type(MGCAMB_parameter_cache), intent(in)   :: mg_par_cache  !< cache containing the parameters

		integer :: i
		character(len=(10)) :: cTemp

        ! 1. MG_flag
        MG_flag = Ini%Read_Int('MG_flag', 0) 
 
        if ( MG_flag /= 0 ) then
            call print_MGCAMB_header  !! call the printing subroutine down below
            write(*,*)
            write(*,*) 'MG_flag:', MG_flag

            write(*,*) 'Debug:', DebugMGCAMB


            ! read GRtrans
            GRtrans = Ini%Read_Double('GRtrans',0.01_dl)
            write(*,*) '    GRtrans:', GRtrans

            ! 1. pure MG models
            if ( MG_flag == 1 ) then

                pure_MG_flag = Ini%Read_Int('pure_MG_flag', 1)

                if ( pure_MG_flag == 1 ) then ! mu-gamma
                    write(*,*) '    MGCAMB: mu-gamma parametrization'
                    mugamma_par = Ini%Read_Int('mugamma_par' , 1)

                    if ( mugamma_par == 1 ) then
                        write(*,*) '        BZ parametrization'
                        B1= Ini%Read_Double('B1',0._dl)
                        B2= Ini%Read_Double('B2',0._dl)
                        lambda1_2= Ini%Read_Double('lambda1_2',0._dl)
                        lambda2_2= Ini%Read_Double('lambda2_2',0._dl)
                        ss= Ini%Read_Double('ss',0._dl)
                    else if ( mugamma_par == 2 ) then
                        write(*,*) '        Planck parametrization'
                        E11     = Ini%Read_Double('E11', 0._dl)
                        E22     = Ini%Read_Double('E22', 0._dl)
                        write(*,*) 'E11, E22', E11, E22
                    else if ( mugamma_par == 3 ) then
                        write(*,*) '        Effective Newton constant'
                        ga      = Ini%Read_Double('ga', 0._dl)
                        nn      = Ini%Read_Double('nn', 0._dl)
                        write(*,*) 'ga, nn:', ga, nn
                    else
                        write(*,*) ' write your own mu-gamma parametrization in mgcamb.f90'
                        stop
                    end if


                else if ( pure_MG_flag == 2 ) then ! mu-Sigma
                    write(*,*) '    MGCAMB: mu-Sigma parametrization'
                    muSigma_par = Ini%Read_Int('musigma_par', 1)
                    if ( muSigma_par == 1 ) then
                        write(*,*) '        DES parametrization'
                        mu0     = Ini%Read_Double('mu0', 0._dl)
                        sigma0  = Ini%Read_Double('sigma0', 0._dl)
                        write(*,*) 'mu0, sigma0:', mu0, sigma0
                    else if ( muSigma_par == 2 ) then
                        write(*,*) 'write you own mu-sigma parametrization in mgcamb.f90'
                        stop
                    else
                        write(*,*) 'Please choose a model in params_MG.ini'
                        stop
                    end if

                else if ( pure_MG_flag == 3 ) then ! Q-R
                    write(*,*) '    MGCAMB: Q-R parametrization'
                    QR_par = Ini%Read_Int('QR_par', 1)
                    if ( QR_par == 1 ) then
                        MGQfix=Ini%Read_Double('MGQfix', 0._dl)
                        MGRfix=Ini%Read_Double('MGRfix', 0._dl)
                    else if ( QR_par == 2 ) then
                        Qnot=Ini%Read_Double('Qnot', 0._dl)
                        Rnot=Ini%Read_Double('Rnot', 0._dl)
                        sss=Ini%Read_Double('sss', 0._dl)
                    else if ( QR_par == 3 ) then
                        write(*,*) 'write your own QR parametrization in mgcamb.f90'
                        stop
                    else
                        write(*,*) 'Please choose a model in params_MG.ini'
                        stop
                    end if

                end if

                ! Checking DE Model
                DE_model = Ini%Read_Int('DE_model', 0)

                write(*,*) 'DE_model:', DE_model

                MGDE_pert = Ini%Read_Logical('MGDE_pert',.false.)

                if ( DE_model == 1 ) then
                    w0DE = Ini%Read_Double('w0DE', -1._dl)
                    ! Assign value to DarkEnergy code component
                    CP%DarkEnergy%w_lam = w0DE
                    CP%DarkEnergy%wa = 0._dl
                else if ( DE_model == 2 ) then
                    w0DE = Ini%Read_Double('w0DE', -1._dl)
                    waDE = Ini%Read_Double('waDE', 0._dl)
                    ! Assign value to DarkEnergy code component
                    CP%DarkEnergy%w_lam = w0DE
                    CP%DarkEnergy%wa = waDE
                else if ( DE_model == 3 ) then
                    write(*,*) 'This is a reconstruction of w_DE(a)'

					do i = 1, nnode+1
						write( cTemp,'(i2)' ) i
						X_arr(i+nnode-1) =  Ini%Read_Double('Funcofw('//trim(adjustl(cTemp))//')', 1._dl)
					end do
					!X_arr(2*nnode) = mgcamb_par_cache%omegav
                else if ( DE_model == 4 ) then
                    write(*,*) 'This will contain the reconstruction of rho_DE(a)'
                    write(*,*) 'Not implemented yet'
                    stop
                else if ( DE_model /= 0 ) then
                    write(*,*) 'Please choose a DE model'
                    stop
                end if



            else if ( MG_flag == 2 ) then
                alt_MG_flag = Ini%Read_Int('alt_MG_flag', 1)
                if ( alt_MG_flag == 1 ) then
                    write(*,*) '    MGCAMB: Linder Gamma'
                    Linder_gamma = Ini%Read_Double('Linder_gamma', 0._dl)
                else if ( alt_MG_flag == 2 ) then
                    write(*,*) 'Please write your alternative MG model in mgcamb.f90'
                    stop
                else
                    write(*,*) 'Please choose a model in params_MG.ini'
                    stop
                end if

                ! Checking DE Model
                DE_model = Ini%Read_Int('DE_model', 0)

                if ( DE_model /= 0 ) then
                    write(*,*) 'alternative MG models supported only with cosmological constant!'
                end if


            else if ( MG_flag == 3 ) then 

                write(*,*) 'standard QSA with mu-gamma parametrization'

                QSA_flag = Ini%Read_Int('QSA_flag', 1)

                    if ( QSA_flag ==  1 ) then
                        write(*,*) '        QSA f(R)'
                        B1 = 4._dl/3._dl
                        lambda1_2= Ini%Read_Double('B0',0._dl) ! it is considered as the B0 parameter here
                        lambda1_2 = (lambda1_2*(299792458.d-3)**2)/(2._dl*mg_par_cache%H0**2)
                        B2 = 0.5d0
                        lambda2_2 = B1* lambda1_2
                        ss = 4._dl

                    else if ( QSA_flag ==  2 ) then
                        write(*,*) '        QSA Symmetron'
                        beta_star = Ini%Read_Double('beta_star', 0._dl)
                        xi_star = Ini%Read_Double ('xi_star', 0._dl)
                        a_star = Ini%Read_Double('a_star', 0._dl)
                        GRtrans = a_star

                    else if ( QSA_flag ==  3 ) then
                        write(*,*) '        QSA Dilaton'
                        ! GENERALIZED DILATON
                        beta0 = Ini%Read_Double('beta0', 0._dl)
                        xi0 = Ini%Read_Double('xi0', 0._dl)
                        DilR = Ini%Read_Double('DilR', 0._dl)
                        DilS = Ini%Read_Double('DilS', 0._dl)

                    else if ( QSA_flag ==  4 ) then
                        write(*,*) '        QSA Hu-Sawicki f(R)'
                        F_R0 = Ini%Read_Double('F_R0', 0._dl)
                        FRn = Ini%Read_Double('FRn', 0._dl)
                        beta0 = 1._dl/sqrt(6._dl)
                    else if ( QSA_flag ==  5 ) then
                        write(*,*) 'Please write your QSA model in mgcamb.f90'
                        stop

                    end if


                    ! Checking DE Model
                    DE_model = Ini%Read_Int('DE_model', 0)

                    if ( DE_model /= 0 ) then
                        write(*,*) 'QSA models supported only with cosmological constant!'
                    end if
 

            else if (MG_flag == 4) then

                CDM_flag = Ini%Read_Int('CDM_flag', 1)
                QSA_flag = Ini%Read_Int('QSA_flag', 1)

                write(*,*) 'only-CDM coupling QSA models'

                if(CDM_flag == 1) then 
                    write(*,*) 'CDM QSA'
                else
                    write(*,*) 'Please choose flag properly'
                    stop
                end if

                if(CDM_flag == 1) then 

                    if ( QSA_flag ==  2 ) then
                        write(*,*) '        QSA Symmetron'
                        beta_star = Ini%Read_Double('beta_star', 0._dl)
                        xi_star = Ini%Read_Double ('xi_star', 0._dl)
                        a_star = Ini%Read_Double('a_star', 0._dl)
                        GRtrans = a_star

                    else if ( QSA_flag ==  3 ) then
                        write(*,*) '        QSA Dilaton'
                        ! GENERALIZED DILATON
                        beta0 = Ini%Read_Double('beta0', 0._dl)
                        xi0 = Ini%Read_Double('xi0', 0._dl)
                        DilR = Ini%Read_Double('DilR', 0._dl)
                        DilS = Ini%Read_Double('DilS', 0._dl)

                    else if ( QSA_flag ==  4 ) then
                        write(*,*) '        QSA Hu-Sawicki f(R)'
                        F_R0 = Ini%Read_Double('F_R0', 0._dl)
                        FRn = Ini%Read_Double('FRn', 0._dl)
                        beta0 = 1._dl/sqrt(6._dl)

                    else
                        write(*,*) 'please choose QSA_flag properly'
                        stop

                    end if

                    DE_model = Ini%Read_Int('DE_model', 0)

                    if ( DE_model /= 0 ) then
                        write(*,*) 'QSA models supported only with cosmological constant!'
                    end if                    

                end if 

            else if (MG_flag == 5) then

                write(*,*) 'direct mu-Sigma parametrization'

                muSigma_flag = Ini%Read_Int('muSigma_flag', 1)
            
                if (muSigma_flag == 1) then

                    pure_MG_flag = Ini%Read_Int('pure_MG_flag', 1)

                    if ( pure_MG_flag == 1 ) then ! mu-gamma
                        write(*,*) '    MGCAMB: mu-gamma parametrization'
                        mugamma_par = Ini%Read_Int('mugamma_par' , 1)
                        if ( mugamma_par == 1 ) then
                            write(*,*) '        BZ parametrization'
                            B1= Ini%Read_Double('B1',0._dl)
                            B2= Ini%Read_Double('B2',0._dl)
                            lambda1_2= Ini%Read_Double('lambda1_2',0._dl)
                            lambda2_2= Ini%Read_Double('lambda2_2',0._dl)
                            ss= Ini%Read_Double('ss',0._dl)
                        else if ( mugamma_par == 2 ) then
                            write(*,*) '        Planck parametrization'
                            E11     = Ini%Read_Double('E11', 0._dl)
                            E22     = Ini%Read_Double('E22', 0._dl)
                            write(*,*) 'E11, E22', E11, E22
                        else if ( mugamma_par == 3 ) then
                            write(*,*) '        Effective Newton constant'
                            ga      = Ini%Read_Double('ga', 0._dl)
                            nn      = Ini%Read_Double('nn', 0._dl)
                            write(*,*) 'ga, nn:', ga, nn
                        else
                            write(*,*) ' write your own mu-gamma parametrization in mgcamb.f90'
                            stop
                        end if


                    else if ( pure_MG_flag == 2 ) then ! mu-Sigma
                        write(*,*) '    MGCAMB: mu-Sigma parametrization'
                        muSigma_par = Ini%Read_Int('musigma_par', 1)
                        if ( muSigma_par == 1 ) then
                            write(*,*) '        DES parametrization'
                            mu0     = Ini%Read_Double('mu0', 0._dl)
                            sigma0  = Ini%Read_Double('sigma0', 0._dl)
                            write(*,*) 'mu0, sigma0:', mu0, sigma0
                        else if ( muSigma_par == 2 ) then
                            write(*,*) 'write you own mu-sigma parametrization in mgcamb.f90'
                            stop
                        else
                            write(*,*) 'Please choose a model in params_MG.ini'
                            stop
                        end if 
                    end if

                    ! Checking DE Model
                    DE_model = Ini%Read_Int('DE_model', 0)

                    write(*,*) 'DE_model:', DE_model

                    MGDE_pert = Ini%Read_Logical('MGDE_pert',.false.)

                    if ( DE_model == 1 ) then
                        w0DE = Ini%Read_Double('w0DE', -1._dl)
                    else if ( DE_model == 2 ) then
                        w0DE = Ini%Read_Double('w0DE', -1._dl)
                        waDE = Ini%Read_Double('waDE', 0._dl)
                    else if ( DE_model == 3 ) then
                        write(*,*) 'This is a reconstruction of w_DE(a)'

                        do i = 1, nnode+1
                            write( cTemp,'(i2)' ) i
                            X_arr(i+nnode-1) =  Ini%Read_Double('Funcofw('//trim(adjustl(cTemp))//')', 1._dl)
                        end do
                        !X_arr(2*nnode) = mgcamb_par_cache%omegav

                    else if ( DE_model == 4 ) then
                        write(*,*) 'This will contain the reconstruction of rho_DE(a)'
                        write(*,*) 'Not implemented yet'
                        stop
                    else if ( DE_model /= 0 ) then
                        write(*,*) 'Please choose a DE model'
                        stop
                    end if

                else if(muSigma_flag == 2) then 

                    alt_MG_flag = Ini%Read_Int('alt_MG_flag', 1)
                    if ( alt_MG_flag == 1 ) then
                        write(*,*) '    MGCAMB: Linder Gamma'
                        Linder_gamma = Ini%Read_Double('Linder_gamma', 0._dl)
                    else if ( alt_MG_flag == 2 ) then
                        write(*,*) 'Please write your alternative MG model in mgcamb.f90'
                        stop
                    else
                        write(*,*) 'Please choose a model in params_MG.ini'
                        stop
                    end if

                    ! Checking DE Model
                    DE_model = Ini%Read_Int('DE_model', 0)

                    if ( DE_model /= 0 ) then
                        write(*,*) 'alternative MG models supported only with cosmological constant!'
                    end if
                
                else if (muSigma_flag == 3) then 

                    write(*,*) 'standard QSA for all-matter case'

                    QSA_flag = Ini%Read_Int('QSA_flag', 1)

                    if ( QSA_flag ==  1 ) then
                        write(*,*) '        QSA f(R)'
                        B1 = 4._dl/3._dl
                        lambda1_2= Ini%Read_Double('B0',0._dl) ! it is considered as the B0 parameter here
                        lambda1_2 = (lambda1_2*(299792458.d-3)**2)/(2._dl*mg_par_cache%H0**2)
                        B2 = 0.5d0
                        lambda2_2 = B1* lambda1_2
                        ss = 4._dl

                    else if ( QSA_flag ==  2 ) then
                        write(*,*) '        QSA Symmetron'
                        beta_star = Ini%Read_Double('beta_star', 0._dl)
                        xi_star = Ini%Read_Double ('xi_star', 0._dl)
                        a_star = Ini%Read_Double('a_star', 0._dl)
                        GRtrans = a_star

                    else if ( QSA_flag ==  3 ) then
                        write(*,*) '        QSA Dilaton'
                        ! GENERALIZED DILATON
                        beta0 = Ini%Read_Double('beta0', 0._dl)
                        xi0 = Ini%Read_Double('xi0', 0._dl)
                        DilR = Ini%Read_Double('DilR', 0._dl)
                        DilS = Ini%Read_Double('DilS', 0._dl)

                    else if ( QSA_flag ==  4 ) then
                        write(*,*) '        QSA Hu-Sawicki f(R)'
                        F_R0 = Ini%Read_Double('F_R0', 0._dl)
                        FRn = Ini%Read_Double('FRn', 0._dl)
                        beta0 = 1._dl/sqrt(6._dl)

                    else
                        write(*,*) 'please choose QSA_flag properly'
                        stop

                    end if

                    DE_model = Ini%Read_Int('DE_model', 0)

                    if ( DE_model /= 0 ) then
                        write(*,*) 'QSA models supported only with cosmological constant!'
                    end if  

                else if(muSigma_flag == 4) then 

                    write(*,*) 'Reconstruction'
                    do i = 1, nnode+1
                        write( cTemp,'(i2)' ) i
                        mu_arr(i+nnode-1) =  Ini%Read_Double('MGCAMB_Mu_idx('//trim(adjustl(cTemp))//')', 1._dl)
                        Sigma_arr(i+nnode-1) =  Ini%Read_Double('MGCAMB_Sigma_idx('//trim(adjustl(cTemp))//')', 1._dl)
                    end do

                    DE_model = Ini%Read_Int('DE_model', 0)

                    write(*,*) 'DE_model:', DE_model

                    if ( DE_model == 1 ) then
                        w0DE = Ini%Read_Double('w0DE', -1._dl)
                    else if ( DE_model == 2 ) then
                        w0DE = Ini%Read_Double('w0DE', -1._dl)
                        waDE = Ini%Read_Double('waDE', 0._dl)
                    else if ( DE_model == 3 ) then
                        write(*,*) 'This is a reconstruction of w_DE(a)'
                        do i = 1, nnode+1
                            write( cTemp,'(i2)' ) i
                            X_arr(i+nnode-1) =  Ini%Read_Double('Funcofw('//trim(adjustl(cTemp))//')', 1._dl)
                        end do
                        !(2*nnode) = mgcamb_par_cache%omegav
                
                    else
                        write(*,*) 'Please choose a DE model'
                        stop
                    end if             

                else 
                    write(*,*) 'Please write your own mu-Sigma model'
                    stop 
                end if

			else if ( MG_flag == 6) then

				write(*,*) 'Reconstruction'
				do i = 1, nnode+1
					write( cTemp,'(i2)' ) i
					mu_arr(i+nnode-1) =  Ini%Read_Double('MGCAMB_Mu_idx('//trim(adjustl(cTemp))//')', 1._dl)
					Sigma_arr(i+nnode-1) =  Ini%Read_Double('MGCAMB_Sigma_idx('//trim(adjustl(cTemp))//')', 1._dl)
				end do

				DE_model = Ini%Read_Int('DE_model', 0)
				write(*,*) 'DE_model:', DE_model

				if ( DE_model == 1 ) then
					w0DE = Ini%Read_Double('w0DE', -1._dl)
				else if ( DE_model == 2 ) then
					w0DE = Ini%Read_Double('w0DE', -1._dl)
					waDE = Ini%Read_Double('waDE', 0._dl)
				else if ( DE_model == 3 ) then
					write(*,*) 'This is a reconstruction of w_DE(a)'
					do i = 1, nnode+1
						write( cTemp,'(i2)' ) i
						X_arr(i+nnode-1) =  Ini%Read_Double('Funcofw('//trim(adjustl(cTemp))//')', 1._dl)
					end do
					!X_arr(2*nnode) = mgcamb_par_cache%omegav
				
				else 
					write(*,*) 'Please choose a DE model'
					stop
				end if

            else 
                write(*,*) 'Please choose a model'
                stop
            end if

        end if 


    end subroutine MGCAMB_read_model_params

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

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that sets the mgcamb_cache to zero
    subroutine MGCAMB_timestep_cache_nullify( mg_cache )
        use precision
        implicit none

        type(MGCAMB_timestep_cache),  intent(inout) :: mg_cache      !< cache containing the time-dependent quantities

        ! 1. Background quantities
        mg_cache%adotoa     = 0._dl
        mg_cache%Hdot       = 0._dl
        mg_cache%grho       = 0._dl
        mg_cache%gpres      = 0._dl
        mg_cache%grhob_t    = 0._dl
        mg_cache%grhoc_t    = 0._dl
        mg_cache%grhog_t    = 0._dl
        mg_cache%grhor_t    = 0._dl
        mg_cache%grhov_t    = 0._dl
        mg_cache%gpresv_t   = 0._dl
        mg_cache%grhonu_t   = 0._dl
        mg_cache%gpresnu_t  = 0._dl

        ! 2. Perturbation quantities
        mg_cache%k          = 0._dl
        mg_cache%k2         = 0._dl
        mg_cache%dgrho      = 0._dl
        mg_cache%dgq        = 0._dl
        mg_cache%pidot_sum  = 0._dl
        mg_cache%dgpi_w_sum = 0._dl
        mg_cache%dgpi       = 0._dl
        mg_cache%dgpi_diff  = 0._dl
        mg_cache%dgpidot    = 0._dl
        mg_cache%rhoDelta   = 0._dl
        mg_cache%rhoDeltadot= 0._dl

        ! 3. MG functions
        mg_cache%mu         = 0._dl
        mg_cache%mudot      = 0._dl
        mg_cache%gamma      = 0._dl
        mg_cache%gammadot   = 0._dl
        mg_cache%q          = 0._dl
        mg_cache%qdot       = 0._dl
        mg_cache%r          = 0._dl
        mg_cache%rdot       = 0._dl

        !> 4. Perturbations evolution variables
        mg_cache%z          = 0._dl
        mg_cache%sigma      = 0._dl
        mg_cache%sigmadot   = 0._dl
        mg_cache%etak       = 0._dl
        mg_cache%etadot     = 0._dl

        !> 5. ISW and lensing realted quantities
        mg_cache%MG_alpha   = 0._dl
        mg_cache%MG_alphadot= 0._dl
        mg_cache%MG_phi     = 0._dl
        mg_cache%MG_phidot  = 0._dl
        mg_cache%MG_psi     = 0._dl
        mg_cache%MG_psidot  = 0._dl
        mg_cache%MG_ISW     = 0._dl
        mg_cache%MG_lensing = 0._dl
        mg_cache%source1    = 0._dl
        mg_cache%source3    = 0._dl

    end subroutine

    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that opens the MGCAMB cache files (for Debug)
    subroutine MGCAMB_open_cache_files
        use precision
        implicit none

        ! 1. Open sources file
        open(unit=111, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_sources.dat', status="new", &
            & action="write")
        write(111,*)  'k  ', 'a  ', 'MG_ISW  ', 'MG_Lensing  ', 'S_T  ', 'S_lensing'

        ! 2 Open MG functions file
        open(unit=222, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_MG_fncs.dat', status="new",&
            & action="write")
        write(222,*)  'k  ', 'a  ', 'mu  ', 'gamma ', 'Q ', 'R ', 'Phi ', 'Psi ', 'dPhi ', 'dPsi '

        ! 3. Open Einstein solutions file
        open(unit=333, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_EinsteinSol.dat', status="new",&
            & action="write")
        write(333,*) 'k  ', 'a  ', 'etak  ', 'z  ', 'sigma  ', 'etadot  ', 'sigmadot  '

        ! 4. Open Perturbation solution file
        open(unit=444, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_PerturbSol.dat', status="new",&
        & action="write")
        write(444,*)  'k  ', 'a  ', 'dgrho  ', 'dgq  ', 'rhoDelta  ', 'dgpi  ', 'pidot_sum  ', 'dgpi_w_sum  '

        ! 5. Open Background file
        open(unit=555, file=TRIM(mgcamb_par_cache%output_root) // 'MGCAMB_debug_Background.dat', status="new",&
            & action="write")
        write(555,*)  'k  ', 'a  ', 'H  ', 'Hdot  ', 'grhov_t  ', 'gpresv_t  '

    end subroutine MGCAMB_open_cache_files


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that closes the MGCAMB cache files (for Debug)
    subroutine MGCAMB_close_cache_files
        use precision
        implicit none

        close(111);close(222); close(333);close(444);close(555)

    end subroutine MGCAMB_close_cache_files


    ! ---------------------------------------------------------------------------------------------
    !> Subroutine that prints the MGCAMB cache on a file
    subroutine MGCAMB_dump_cache( a, mg_cache )
        use precision
        implicit none

        real(dl), intent(in) :: a   !< scale factor
        type(MGCAMB_timestep_cache),  intent(in) :: mg_cache      !< cache containing the time-dependent quantities
        character(*), parameter :: cache_output_format = 'e18.8'


        ! 1. Write the sources
        write(111,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%MG_ISW, mg_cache%MG_Lensing,&
                                                    & mg_cache%source1, mg_cache%source3

        ! 2. Write the MG functions and the potentials
        write(222,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%mu, mg_cache%gamma, mg_cache%q, mg_cache%r, &
                                                & mg_cache%MG_phi, mg_cache%MG_psi, mg_cache%MG_phidot, mg_cache%MG_psidot

        ! 3. Write the Einstein equations solutions
        write(333,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%etak, mg_cache%z, mg_cache%sigma,&
                                                & mg_cache%etadot,mg_cache%sigmadot

        ! 4. Write the Perturbations Solutions
        write(444,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%dgrho, mg_cache%dgq, mg_cache%rhoDelta,&
                                                    & mg_cache%dgpi, mg_cache%pidot_sum, mg_cache%dgpi_w_sum

        !5. Write the background
        write(555,'(14'//cache_output_format//')') mg_cache%k, a, mg_cache%adotoa, mg_cache%Hdot, mg_cache%grhov_t,&
                                                    & mg_cache%gpresv_t



    end subroutine MGCAMB_dump_cache



end module MGCAMB

