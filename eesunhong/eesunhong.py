from ctypes import CDLL, POINTER, c_int, c_float
from pathlib import Path

module_parent_directory = Path(__file__).parent.parent
if module_parent_directory.joinpath('build').exists():
    library_directory = module_parent_directory.joinpath('build/lib')
else:
    library_directory = module_parent_directory.joinpath('lib')
library_path = library_directory.joinpath('libeesunhong_fortran_library.dylib')
if not library_path.exists():
    library_path = library_directory.joinpath('libeesunhong_fortran_library.so')
if not library_path.exists():
    library_path = library_directory.joinpath('libeesunhong_fortran_library.dll')
if not library_path.exists():
    raise FileNotFoundError

fortran_library = CDLL(str(library_path))
sort_light_curve_data_by_time = fortran_library.sort_light_curve_data_by_time
sort_light_curve_data_by_time.argtypes = [POINTER(c_int),
                                          POINTER(c_float),
                                          POINTER(c_float),
                                          POINTER(c_float),
                                          POINTER(c_int),
                                          POINTER(c_int)]
