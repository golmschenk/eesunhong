import errno
import os
from ctypes import CDLL, POINTER, c_int, c_double, c_void_p, c_char_p
from pathlib import Path
import platform

module_parent_directory = Path(__file__).parent.parent
library_directory = module_parent_directory
if library_directory.name == 'src' and library_directory.joinpath('../build').exists():
    library_directory = library_directory.joinpath('../build')
winmode = None

platform_system_name = platform.system()
if platform_system_name == 'Darwin':
    eesunhong_fortran_library_path = library_directory.joinpath('libeesunhong_fortran_library.dylib')
elif platform_system_name == 'Linux':
    eesunhong_fortran_library_path = library_directory.joinpath('libeesunhong_fortran_library.so')
elif platform_system_name == 'Windows':
    winmode = 0
    os.add_dll_directory(f'{library_directory.absolute()}')
    os.add_dll_directory(f'{library_directory.absolute()}/lib')
    os.add_dll_directory(f'{library_directory.absolute()}/eesunhong.libs')
    eesunhong_fortran_library_path = library_directory.joinpath('eesunhong_fortran_library.dll')
    if not eesunhong_fortran_library_path.exists():
        eesunhong_fortran_library_path = library_directory.joinpath('libeesunhong_fortran_library.dll')
else:
    raise ValueError(f'Platform system {platform_system_name} is not supported.')

if not eesunhong_fortran_library_path.exists():
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), str(eesunhong_fortran_library_path.absolute()))
eesunhong_fortran_library = CDLL(str(eesunhong_fortran_library_path), winmode=winmode)

sort_light_curve_data_by_time = eesunhong_fortran_library.sort_light_curve_data_by_time
sort_light_curve_data_by_time.argtypes = [POINTER(c_int),
                                          POINTER(c_double),
                                          POINTER(c_double),
                                          POINTER(c_double),
                                          POINTER(c_int),
                                          POINTER(c_int)]

compute_parallax_using_geo_par = eesunhong_fortran_library.geo_par
compute_parallax_using_geo_par.argtypes = [POINTER(c_double),
                                           POINTER(c_double),
                                           POINTER(c_double),
                                           POINTER(c_double),
                                           POINTER(c_double),
                                           POINTER(c_double)]

create_vbbl = eesunhong_fortran_library.create_vbbl
create_vbbl.restype = c_void_p

destroy_vbbl = eesunhong_fortran_library.destroy_vbbl
destroy_vbbl.argtypes = [c_void_p]

set_object_coordinates_for_vbbl = eesunhong_fortran_library.set_object_coordinates_for_vbbl
set_object_coordinates_for_vbbl.argtypes = [c_void_p, c_char_p, c_char_p]

compute_parallax_for_vbbl = eesunhong_fortran_library.compute_parallax_for_vbbl
compute_parallax_for_vbbl.argtypes = [c_void_p, c_double, c_double, POINTER(c_double)]

set_parallax_system_for_vbbl = eesunhong_fortran_library.set_parallax_system_for_vbbl
set_parallax_system_for_vbbl.argtypes = [c_void_p, c_int]
