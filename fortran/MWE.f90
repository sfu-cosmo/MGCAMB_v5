program MWE_ReadParams
    use IniObjects
    use ModGravityInterface
    use MGCAMB

    implicit none

    type(TMuGammaParameterization) :: mgParam
    type(TIniFile) :: IniFile
    logical :: error

    ! Initialize INI file object
    call IniFile%Open('test1_params.ini', error)

    ! Check for errors during file opening
    if (error) then
        print *, "Error opening INI file"
        stop
    end if

    ! Read parameters from the INI file
    call mgParam%ReadParams(IniFile)

    ! Optionally print the read values for verification
    call mgParam%PrintAttributes()

    ! Close the INI file
    call IniFile%Close()

end program MWE_ReadParams
