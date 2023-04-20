module grackle_recipes_wrappers
    implicit none
contains
    subroutine sort5_c_wrapper(N, RA, RB, RC, IC, id, WKSP, IWKSP, JWKSP) &
            bind(c, name = 'sort5')
        use stdlib_kinds, only : dp, int32
        use iso_c_binding, only : c_float, c_int
        implicit none
        integer(c_int), intent(in) :: N
        real(c_float), intent(inout) :: RA(N)
        real(c_float), intent(inout) :: RB(N)
        real(c_float), intent(inout) :: RC(N)
        integer(c_int), intent(inout) :: IC(N)
        integer(c_int), intent(inout) :: id(N)
        real(c_float), intent(in) :: WKSP(N)
        integer(c_int), intent(in) :: IWKSP(N)
        integer(c_int), intent(in) :: JWKSP(N)
        real(dp) :: f_RA(N)
        real(dp) :: f_RB(N)
        real(dp) :: f_RC(N)
        integer(int32) :: f_IC(N)
        integer(int32) :: f_id(N)
        real(dp) :: f_WKSP(N)
        integer(int32) :: f_IWKSP(N)
        integer(int32) :: f_JWKSP(N)

        f_RA = RA
        f_RB = RB
        f_RC = RC
        f_IC = IC
        f_id = id

        call sort5(int(N, int32), f_RA, f_RB, f_RC, f_IC, f_id, real(f_WKSP, dp), int(f_IWKSP, int32), int(f_JWKSP, int32))

        RA = f_RA
        RB = f_RB
        RC = f_RC
        IC = f_IC
        id = f_id

    end subroutine sort5_c_wrapper
end module grackle_recipes_wrappers
