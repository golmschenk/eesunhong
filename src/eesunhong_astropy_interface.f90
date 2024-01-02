module eesunhong_astropy_interface
    implicit none
contains
    function outer_multiply(outer_factor0, outer_factor1) result(outer_product)
        use, intrinsic :: iso_c_binding, only : c_float
        use stdlib_kinds, only : sp
        implicit none
        real(sp), intent(in) :: outer_factor0
        real(sp), intent(in) :: outer_factor1
        real(sp) :: outer_product

        interface
            real(c_float) function multiply(factor0, factor1) bind(c, name = 'multiply')
                import :: c_float
                real(c_float), value, intent(in) :: factor0
                real(c_float), value, intent(in) :: factor1
            end function
        end interface

        outer_product = multiply(outer_factor0, outer_factor1)
    end function
end module
