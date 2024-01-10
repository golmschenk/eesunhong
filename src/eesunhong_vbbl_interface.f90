module eesunhong_vbbl_interface
    use, intrinsic :: iso_c_binding, only : c_float, c_ptr, c_double, c_char
    implicit none
    interface
        function create_vbbl_interface() bind(c, name = 'create_vbbl') result(vbbl)
            import :: c_ptr
            type(c_ptr) :: vbbl
        end function create_vbbl_interface

        subroutine destroy_vbbl_interface(vbbl) bind(c, name = 'destroy_vbbl')
            import :: c_ptr
            type(c_ptr), value :: vbbl
        end subroutine destroy_vbbl_interface

        subroutine set_object_coordinates_for_vbbl_interface(vbbl, coordinates_file, directory_for_satellite_tables) &
                bind(c, name = 'set_object_coordinates_for_vbbl')
            import :: c_ptr, c_char
            type(c_ptr), value :: vbbl
            character(kind=c_char) :: coordinates_file(*)
            character(kind=c_char) :: directory_for_satellite_tables(*)
        end subroutine set_object_coordinates_for_vbbl_interface

        subroutine compute_parallax_for_vbbl_interface(vbbl, t, t0, Et) &
                bind(c, name = 'compute_parallax_for_vbbl')
            import :: c_ptr, c_double
            type(c_ptr), value :: vbbl
            real(c_double), value :: t
            real(c_double), value :: t0
            real(c_double) :: Et(2)
        end subroutine compute_parallax_for_vbbl_interface
    end interface
contains
    function create_vbbl() result(vbbl)
        type(c_ptr) :: vbbl
        vbbl = create_vbbl_interface()
    end function create_vbbl

    subroutine destroy_vbbl(vbbl)
        type(c_ptr), intent(in) :: vbbl

        call destroy_vbbl_interface(vbbl)
    end subroutine destroy_vbbl

    subroutine set_object_coordinates_for_vbbl(vbbl, coordinates_file, directory_for_satellite_tables)
        type(c_ptr), intent(in) :: vbbl
        character(len=*), intent(in) :: coordinates_file
        character(len=*), intent(in) :: directory_for_satellite_tables

        call set_object_coordinates_for_vbbl_interface(vbbl, coordinates_file, directory_for_satellite_tables)
    end subroutine set_object_coordinates_for_vbbl

    subroutine compute_parallax_for_vbbl(vbbl, t, t0, Et)
        type(c_ptr), intent(in) :: vbbl
        real(kind=c_double), intent(in) :: t
        real(kind=c_double), intent(in) :: t0
        real(kind=c_double), intent(out) :: Et(2)

        call compute_parallax_for_vbbl_interface(vbbl, t, t0, Et)
    end subroutine compute_parallax_for_vbbl
end module eesunhong_vbbl_interface