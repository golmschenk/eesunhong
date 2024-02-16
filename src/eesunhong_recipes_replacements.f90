module eesunhong_recipes_replacements
    use stdlib_kinds, only : dp, int32
    use iso_c_binding, only : c_double, c_int
    use stdlib_sorting, only : sort_index, int_size
    implicit none
contains
    subroutine sort_light_curve_data_by_time(number_of_data_points, time, magnification, sig, iclr, iclrind)
        implicit none
        integer(int32), intent(in) :: number_of_data_points
        real(dp), intent(inout) :: time(:)
        real(dp), intent(inout) :: magnification(:)
        real(dp), intent(inout) :: sig(:)
        integer(int32), intent(inout) :: iclr(:)
        integer(int32), intent(inout) :: iclrind(:)
        integer(int_size), allocatable :: index_array(:)

        allocate(index_array(number_of_data_points))

        call sort_index(time(:number_of_data_points), index_array)
        magnification(:number_of_data_points) = magnification(index_array)
        sig(:number_of_data_points) = sig(index_array)
        iclr(:number_of_data_points) = iclr(index_array)
        iclrind(:number_of_data_points) = iclrind(index_array)
    end subroutine sort_light_curve_data_by_time

    subroutine sort_light_curve_data_by_time_c_wrapper(number_of_data_points, time, magnification, sig, iclr, iclrind) &
            bind(c, name = 'sort_light_curve_data_by_time')
        implicit none
        integer(c_int), intent(in) :: number_of_data_points
        real(c_double), intent(inout) :: time(number_of_data_points)
        real(c_double), intent(inout) :: magnification(number_of_data_points)
        real(c_double), intent(inout) :: sig(number_of_data_points)
        integer(c_int), intent(inout) :: iclr(number_of_data_points)
        integer(c_int), intent(inout) :: iclrind(number_of_data_points)
        real(dp) :: f_time(number_of_data_points)
        real(dp) :: f_magnification(number_of_data_points)
        real(dp) :: f_sig(number_of_data_points)
        integer(int32) :: f_iclr(number_of_data_points)
        integer(int32) :: f_iclrind(number_of_data_points)

        f_time = time
        f_magnification = magnification
        f_sig = sig
        f_iclr = iclr
        f_iclrind = iclrind

        call sort_light_curve_data_by_time(int(number_of_data_points, int32), f_time, f_magnification, f_sig, f_iclr, f_iclrind)

        time = f_time
        magnification = f_magnification
        sig = f_sig
        iclr = f_iclr
        iclrind = f_iclrind
    end subroutine sort_light_curve_data_by_time_c_wrapper

    function brent_wrapper_with_additional_lens_arguments(FUNC, X1, X2, y1, y2, TOL, ssx, ssy, Ustar2, sep, eps1) result(result)
        use stdlib_kinds, only : dp, int32
        use stdlib_sorting, only : sort_index, int_size
        use root_module, only : root_scalar
        implicit none
        interface
            function FUNC(A_, y1_, y2_, ssx_, ssy_, Ustar2_, sep_, eps1_) result(result)
                use stdlib_kinds, only : dp
                implicit none
                real(dp), intent(in) :: A_
                real(dp), intent(in) :: y1_
                real(dp), intent(in) :: y2_
                real(dp), intent(in) :: ssx_
                real(dp), intent(in) :: ssy_
                real(dp), intent(in) :: Ustar2_
                real(dp), intent(in) :: sep_
                real(dp), intent(in) :: eps1_
                real(dp) :: result
            end function FUNC
        end interface
        real(dp), intent(in) :: X1
        real(dp), intent(in) :: X2
        real(dp), intent(in) :: y1
        real(dp), intent(in) :: y2
        real(dp), intent(in) :: TOL
        real(dp), intent(in) :: ssx
        real(dp), intent(in) :: ssy
        real(dp), intent(in) :: Ustar2
        real(dp), intent(in) :: sep
        real(dp), intent(in) :: eps1
        real(dp) :: xzero
        real(dp) :: fzero
        integer(int32) :: iflag
        real(dp) :: result
        real(dp) :: rtol = 1.0e-9_dp

        call root_scalar('brent', partial, X1, X2, xzero, fzero, iflag, rtol=TOL)
        result = xzero
    contains
        function partial(A) result(result)
            use stdlib_kinds, only : dp
            implicit none
            real(dp), intent(in) :: A
            real(dp) :: result
            result = FUNC(A, y1, y2, ssx, ssy, Ustar2, sep, eps1)
        end function partial
    end function brent_wrapper_with_additional_lens_arguments
end module eesunhong_recipes_replacements
