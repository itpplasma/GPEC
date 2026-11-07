! Usage: pentrc_zero_omega_e <kinetic profile file>
! read_kin must keep the controls exact: wefac = 0 with wpfac = 1 gives
! omega_E = 0 at every knot, and wpfac = 1 leaves wefac*omega_E unchanged.
program pentrc_zero_omega_e
    use params, only: r8
    use inputs, only: read_kin, kin, chi1
    use spline_mod, only: spline_dealloc
    implicit none
    character(512) :: file
    real(r8), allocatable :: welec_full(:)

    call get_command_argument(1, file)
    chi1 = 1.0_r8

    call read_kin(file, 1, 6, 2, 12, 1.0_r8, 1.0_r8, 1.0_r8, 1.0_r8, .false.)
    welec_full = kin%fs(:, 5)
    if (all(welec_full == 0)) error stop 'fixture has no rotation'

    call spline_dealloc(kin)
    call read_kin(file, 1, 6, 2, 12, 1.0_r8, 1.0_r8, 0.0_r8, 1.0_r8, .false.)
    print '(a,es12.4)', 'max |omega_E| for wefac = 0: ', maxval(abs(kin%fs(:, 5)))
    if (any(kin%fs(:, 5) /= 0)) error stop 'wefac = 0 does not give omega_E = 0'

    call spline_dealloc(kin)
    call read_kin(file, 1, 6, 2, 12, 1.0_r8, 1.0_r8, 0.5_r8, 1.0_r8, .false.)
    if (any(kin%fs(:, 5) /= 0.5_r8*welec_full)) &
        error stop 'wpfac = 1 changes omega_E'
    print '(a)', 'PASS: exact omega_E controls'
end program pentrc_zero_omega_e
