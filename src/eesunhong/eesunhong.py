import errno
import os
from ctypes import CDLL, POINTER, c_int, c_float
from pathlib import Path
import platform

module_parent_directory = Path(__file__).parent.parent
library_directory = module_parent_directory
if library_directory.name == 'src' and library_directory.joinpath('../build').exists():
    library_directory = library_directory.joinpath('../build')
winmode = None

platform_system_name = platform.system()
if platform_system_name == 'Darwin':
    library_path = library_directory.joinpath('libeesunhong_fortran_library.dylib')
elif platform_system_name == 'Linux':
    library_path = library_directory.joinpath('libeesunhong_fortran_library.so')
elif platform_system_name == 'Windows':
    winmode = 0
    os.add_dll_directory(f'{library_directory.absolute()}')
    os.add_dll_directory(f'{library_directory.absolute()}/lib')
    os.add_dll_directory(f'{library_directory.absolute()}/eesunhong.libs')
    library_path = library_directory.joinpath('eesunhong_fortran_library.dll')
    if not library_path.exists():
        library_path = library_directory.joinpath('libeesunhong_fortran_library.dll')
else:
    raise ValueError(f'Platform system {platform_system_name} is not supported.')

if not library_path.exists():
    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), str(library_path.absolute()))
fortran_library = CDLL(str(library_path), winmode=winmode)
sort_light_curve_data_by_time = fortran_library.sort_light_curve_data_by_time
sort_light_curve_data_by_time.argtypes = [POINTER(c_int),
                                          POINTER(c_float),
                                          POINTER(c_float),
                                          POINTER(c_float),
                                          POINTER(c_int),
                                          POINTER(c_int)]
