module eesunhong_vbbl_interface
    use, intrinsic :: iso_c_binding, only : c_int, c_float, c_ptr, c_double, c_char, c_null_char
    use stdlib_kinds, only : dp
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

        subroutine set_parallax_system_for_vbbl_interface(vbbl, value) bind(c, name = 'set_parallax_system_for_vbbl')
            import :: c_ptr, c_int
            type(c_ptr), value :: vbbl
            integer(c_int), value :: value
        end subroutine set_parallax_system_for_vbbl_interface

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

    subroutine set_parallax_system_for_vbbl(vbbl, value)
        type(c_ptr), intent(in) :: vbbl
        integer(c_int), intent(in) :: value

        call set_parallax_system_for_vbbl_interface(vbbl, value)
    end subroutine set_parallax_system_for_vbbl

    subroutine set_object_coordinates_for_vbbl(vbbl, coordinates_file, directory_for_satellite_tables)
        type(c_ptr), intent(in) :: vbbl
        character(len=*), intent(in) :: coordinates_file
        character(len=*), intent(in) :: directory_for_satellite_tables

        call set_object_coordinates_for_vbbl_interface(vbbl, coordinates_file // c_null_char, &
                directory_for_satellite_tables // c_null_char)
    end subroutine set_object_coordinates_for_vbbl

    subroutine compute_parallax_for_vbbl(vbbl, t, t0, Et)
        type(c_ptr), intent(in) :: vbbl
        real(kind=c_double), intent(in) :: t
        real(kind=c_double), intent(in) :: t0
        real(kind=c_double), intent(out) :: Et(2)

        call compute_parallax_for_vbbl_interface(vbbl, t, t0, Et)
    end subroutine compute_parallax_for_vbbl

    subroutine create_vbbl_coordinates_file(a, b, c, d, e, f)
        implicit none
        real(dp), intent(in) :: a, b, c, d, e, f
        character(len=50) :: str_a, str_b, str_c, str_d, str_e, str_f
        character(len=256) :: line_to_write
        integer :: iounit = 132

        write(str_a, "(F0.4)") a
        write(str_b, "(F0.4)") b
        write(str_c, "(F0.4)") c
        write(str_d, "(F0.4)") d
        write(str_e, "(F0.4)") e
        write(str_f, "(F0.4)") f

        line_to_write = trim(adjustl(str_a)) // ":" // trim(adjustl(str_b)) // ":" // trim(adjustl(str_c)) // &
                " " // trim(adjustl(str_d)) // ":" // trim(adjustl(str_e)) // ":" // trim(adjustl(str_f))

        open(newunit=iounit, file="tmp_coordinates.txt", status="replace", action="write")
        write(iounit, '(A)') line_to_write
        close(iounit)
    end subroutine create_vbbl_coordinates_file

    subroutine delete_vbbl_coordinates_file()
        implicit none
        integer :: iounit = 132
        integer :: stat = -1

        open(unit=iounit, iostat=stat, file='tmp_coordinates.txt', status='old')
        if (stat == 0) close(iounit, status='delete')
    end subroutine delete_vbbl_coordinates_file
end module eesunhong_vbbl_interface
