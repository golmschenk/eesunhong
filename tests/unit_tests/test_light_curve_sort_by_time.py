from ctypes import POINTER, c_double, c_int

import numpy as np
from eesunhong.eesunhong import sort_light_curve_data_by_time


def test_light_curve_sort_by_time():
    time = np.array([9004.0, 9000.0, 9001.0, 9003.0, 9002.0], dtype=c_double)
    magnification = np.array([0.0, 10.0, 20.0, 30.0, 40.0], dtype=c_double)
    sig = np.array([2, 3, 5, 1, 7], dtype=c_double)
    iclr = np.array([6, 7, 8, 9, 10], dtype=c_int)
    iclrind = np.array([15, 14, 13, 12, 11], dtype=c_int)

    sort_light_curve_data_by_time(c_int(len(time)),
                                  time.ctypes.data_as(POINTER(c_double)),
                                  magnification.ctypes.data_as(POINTER(c_double)),
                                  sig.ctypes.data_as(POINTER(c_double)),
                                  iclr.ctypes.data_as(POINTER(c_int)),
                                  iclrind.ctypes.data_as(POINTER(c_int)))

    expected_time = np.array([9000.0, 9001.0, 9002.0, 9003.0, 9004.0], dtype=c_double)
    expected_magnification = np.array([10.0, 20.0, 40.0, 30.0, 0.0], dtype=c_double)
    expected_sig = np.array([3, 5, 7, 1, 2], dtype=c_double)
    expected_iclr = np.array([7, 8, 10, 9, 6], dtype=c_int)
    expected_iclrind = np.array([14, 13, 11, 12, 15], dtype=c_int)

    assert np.allclose(time, expected_time)
    assert np.allclose(magnification, expected_magnification)
    assert np.allclose(sig, expected_sig)
    assert np.allclose(iclr, expected_iclr)
    assert np.allclose(iclrind, expected_iclrind)