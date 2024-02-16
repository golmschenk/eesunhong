module eesunhong_geo_par_interface
    use stdlib_kinds, only : dp, int32
    use iso_c_binding, only : c_double, c_int
    implicit none
contains
    subroutine geo_par_c_wrapper(qn, qe, hjd, alpha, delta, tfix) &
            bind(c, name = 'geo_par')
        implicit none
        real(c_double), intent(out) :: qn
        real(c_double), intent(out) :: qe
        real(c_double), intent(in) :: hjd
        real(c_double), intent(in) :: alpha
        real(c_double), intent(in) :: delta
        real(c_double), intent(in) :: tfix
        real(dp) :: f_qn
        real(dp) :: f_qe
        real(dp) :: f_hjd
        real(dp) :: f_alpha
        real(dp) :: f_delta
        real(dp) :: f_tfix

        f_hjd = hjd
        f_alpha = alpha
        f_delta = delta
        f_tfix = tfix

        call geo_par(f_qn, f_qe, f_hjd, f_alpha, f_delta, f_tfix)

        qn = f_qn
        qe = f_qe
    end subroutine geo_par_c_wrapper
end module eesunhong_geo_par_interface