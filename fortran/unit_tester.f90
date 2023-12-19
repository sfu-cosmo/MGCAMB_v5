program tester
    use CAMB
    use MGCAMB
    implicit none

    logical :: ReadParamFile = .false.

    type(CAMBparams) CPcz !defined in model.f90
    type(CAMBdata) res !contains computed quantities and functions for results (results.f90)
    
    Type(TIniFile) :: Ini
    logical bad


    !Need to set default classes for initial power spectrum, recombination etc.
    call CAMB_SetDefParams(CPcz)


    !ReadParamFile = CAMB_ReadParamFile( CPcz, "test1_params.ini", 30 )
    !write(*,*) "ReadParamFile = ", ReadParamFile

    CPcz%ombh2  = .0222_dl
    CPcz%omk  = 0._dl
    CPcz%H0      = 67._dl
    select type(InitPower=>CPcz%InitPower)
        class is (TInitialPowerLaw)
            InitPower%As = 2.1e-9
            InitPower%ns  = 1
            InitPower%r = 1
            !we don't use r here since we generate the Cls separately
            !so set to 1, and then put in the ratio after calculating the Cls
            write(*,*) "cacca", InitPower%r
    end select
    
    select type(ModGravity=>CPcz%ModGravity)
        class is (TMuGammaParameterization)
            !ModGravity%MG_flag=1
            call ModGravity%ReadParams( Ini )
            write (*,*) "debug", ModGravity%B1
        class default
            write(*,*) "I'm here"
    end select

    write(*,*) "set initial power spectrum"

    call Ini%Open( "test1_params.ini", bad, .false.)

    write(*,*) "reopened ini file", bad
    
    select type(ModGravity=>CPcz%ModGravity)
        class is (TMuGammaParameterization)
            !ModGravity%MG_flag=1
            call ModGravity%ReadParams( Ini )
            write (*,*) "debug"
            write (*,*) ModGravity%B1
        class default
            write(*,*) "I'm here"
    end select
    call Ini%Close()

    !write (*,*) CPcz%ModGravity%B1


    CPcz%WantScalars = .true.
    CPcz%WantTensors = .true.

    CPcz%Max_l=2500
    CPcz%Max_eta_k=6000
    CPcz%Max_l_tensor=200
    CPcz%Max_eta_k_tensor=500

    !call CAMB_GetResults(res, CPcz)


    
end program Tester
    
    
    