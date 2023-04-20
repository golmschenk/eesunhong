from ctypes import CDLL, POINTER, c_int, c_float
from pathlib import Path

library_path = Path('./build/libeesunhong.dylib')
if not library_path.exists():
    library_path = Path('./build/libeesunhong.so')
if not library_path.exists():
    library_path = Path('./build/libeesunhong.dll')
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
