program pentrc_nu_profile
    !-----------------------------------------------------------------------
    ! Regression test: PENTRC collisionality profiles when the kinetic
    ! profile ends before psi_n = 1.
    !
    ! Writes an analytic pedestal profile ending at psi_n = 0.995, reads it
    ! with read_kin and requires nu_i, nu_e > 0 and agreement of the splined
    ! nu_i, nu_e with the direct formula evaluated from the splined n, T
    ! to within 1 % on psi_n in [0.1, 0.99].
    !
    ! Build and run from install/ after "make pentrc":
    !   make -C install pentrc-nutest
    ! or see the command in regression/test/README_pentrc_nu.txt.
    !-----------------------------------------------------------------------
    use params, only : r8,mp,me
    use spline_mod, only : spline_eval
    use inputs, only : read_kin,kin
    implicit none

    integer, parameter :: np = 200, zi = 1, zimp = 6, mi = 2, mimp = 12
    character(512) :: file
    integer :: i, unit_kin, nbad
    real(r8) :: psi, x, ni, ti, te, zeff, zpitch, lnl, nui, nue, err_i, err_e
    real(r8) :: maxerr_i, maxerr_e, minnu
    real(r8), parameter :: ev = 1.602e-19_r8

    file = "nu_test.kin"
    open(newunit=unit_kin, file=trim(file), status="replace")
    write(unit_kin,'(a)') "psi ni ne ti te wexb"
    do i = 1, np
        psi = 0.995_r8*real(i,r8)/real(np,r8)
        x = 0.5_r8*(1.0_r8 - tanh((psi-0.96_r8)/0.015_r8))
        write(unit_kin,'(6es18.9)') psi, 6.0e19_r8*(0.2_r8+0.8_r8*x), &
            6.0e19_r8*(0.2_r8+0.8_r8*x), 5000.0_r8*(0.1_r8+0.9_r8*x), &
            4000.0_r8*(0.1_r8+0.9_r8*x), -3.0e3_r8*(1.0_r8-psi)
    end do
    close(unit_kin)

    call read_kin(file,zi,zimp,mi,mimp,1.0_r8,1.0_r8,1.0_r8,1.0_r8,.false.)

    nbad = 0
    maxerr_i = 0.0_r8
    maxerr_e = 0.0_r8
    minnu = huge(1.0_r8)
    do i = 0, 890
        psi = 0.1_r8 + 0.001_r8*real(i,r8)
        call spline_eval(kin,psi,0)
        ni = kin%f(1)
        ti = kin%f(3)/ev
        te = kin%f(4)/ev
        zeff = zimp-(kin%f(1)/kin%f(2))*zi*(zimp-zi)
        zpitch = 1.0_r8+(1.0_r8+mimp)/(2.0_r8*mimp)*zimp*(zeff-1.0_r8)/(zimp-zeff)
        lnl = 17.3_r8-0.5_r8*log(kin%f(2)/1.0e20_r8)+1.5_r8*log(kin%f(4)/1.602e-16_r8)
        nui = (zpitch/3.5e17_r8)*ni*lnl/(sqrt(1.0_r8*mi)*ti**1.5_r8)
        nue = (zpitch/3.5e17_r8)*kin%f(2)*lnl/(sqrt(me/mp)*te**1.5_r8)
        err_i = abs(kin%f(7)/nui-1.0_r8)
        err_e = abs(kin%f(8)/nue-1.0_r8)
        maxerr_i = max(maxerr_i,err_i)
        maxerr_e = max(maxerr_e,err_e)
        minnu = min(minnu,kin%f(7),kin%f(8))
        if(kin%f(7) <= 0.0_r8 .or. kin%f(8) <= 0.0_r8) nbad = nbad+1
    end do

    print '(a,es10.3)', "min(nu_i,nu_e)       = ", minnu
    print '(a,es10.3)', "max rel err nu_i     = ", maxerr_i
    print '(a,es10.3)', "max rel err nu_e     = ", maxerr_e
    if(nbad > 0 .or. maxerr_i > 0.01_r8 .or. maxerr_e > 0.01_r8) then
        print '(a,i0,a)', "FAIL: ", nbad, " non-positive nu points or error > 1 %"
        error stop 1
    end if
    print '(a)', "PASS"
end program pentrc_nu_profile
