module eesunhong_real_complex_conversion
    implicit none
contains
    subroutine from_2d_real_to_complex(real_array, complex_array)
        use stdlib_error, only : check
        use stdlib_kinds, only : dp, int32
        implicit none
        real(dp), dimension(:, :), intent(in) :: real_array
        complex(dp), dimension(:), intent(out) :: complex_array
        integer(int32) :: real_array_shape(2)
        integer(int32) :: complex_array_shape(1)
        integer(int32) :: index

        real_array_shape = shape(real_array)
        complex_array_shape = shape(complex_array)
        call check(real_array_shape(1) == 2, msg = 'Expected an a real array with a first dimension of 2.')
        call check(real_array_shape(2) == complex_array_shape(1), &
                msg = 'Expected the second dimension of the real array to match the first of the complex array.')

        do index = lbound(complex_array, dim=1), ubound(complex_array, dim=1)
            complex_array(index) = cmplx(real_array(1, index), real_array(2, index))
        end do
    end subroutine from_2d_real_to_complex

    subroutine from_complex_to_2d_real(complex_array, real_array)
        use stdlib_error, only : check
        use stdlib_kinds, only : dp, int32
        implicit none
        complex(dp), dimension(:), intent(in) :: complex_array
        real(dp), dimension(:, :), intent(out) :: real_array
        integer(int32) :: complex_array_shape(1)
        integer(int32) :: real_array_shape(2)
        integer(int32) :: index

        complex_array_shape = shape(complex_array)
        real_array_shape = shape(real_array)
        call check(real_array_shape(1) == 2, msg = 'Expected an a real array with a first dimension of 2.')
        call check(real_array_shape(2) == complex_array_shape(1), &
                msg = 'Expected the second dimension of the real array to match the first of the complex array.')

        do index = lbound(complex_array, dim=1), ubound(complex_array, dim=1)
            real_array(1, index) = real(complex_array(index))
            real_array(2, index) = aimag(complex_array(index))
        end do
    end subroutine from_complex_to_2d_real
end module eesunhong_real_complex_conversion
