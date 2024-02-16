from ctypes import POINTER, c_double, c_int
from pathlib import Path

import numpy as np
import pytest

from eesunhong.eesunhong import compute_parallax_using_geo_par, create_vbbl, destroy_vbbl, \
    set_object_coordinates_for_vbbl, compute_parallax_for_vbbl, set_parallax_system_for_vbbl


def test_parallax_computation_using_geo_par():
    qn_case0 = np.array(np.nan, dtype=c_double)
    qe_case0 = np.array(np.nan, dtype=c_double)
    t_case0 = np.array(9033.2148899999447, dtype=c_double)
    alpha_case0 = np.array(268.43250000000000, dtype=c_double)
    delta_case0 = np.array(-32.589311111111115, dtype=c_double)
    tfix_case0 = np.array(9100.8972159999994, dtype=c_double)

    compute_parallax_using_geo_par(qn_case0.ctypes.data_as(POINTER(c_double)),
                                   qe_case0.ctypes.data_as(POINTER(c_double)),
                                   t_case0.ctypes.data_as(POINTER(c_double)),
                                   alpha_case0.ctypes.data_as(POINTER(c_double)),
                                   delta_case0.ctypes.data_as(POINTER(c_double)),
                                   tfix_case0.ctypes.data_as(POINTER(c_double)))

    expected_qn_case0 = np.array(5.1377649825789830e-2, dtype=c_double)
    expected_qe_case0 = np.array(0.52891219662449951, dtype=c_double)

    assert qe_case0 == pytest.approx(expected_qe_case0)
    assert qn_case0 == pytest.approx(expected_qn_case0)

    qn_case1 = np.array(np.nan, dtype=c_double)
    qe_case1 = np.array(np.nan, dtype=c_double)
    t_case1 = np.array(9090.8258499999993, dtype=c_double)
    alpha_case1 = np.array(268.43250000000000, dtype=c_double)
    delta_case1 = np.array(-32.589311111111115, dtype=c_double)
    tfix_case1 = np.array(9100.8972159999994, dtype=c_double)

    compute_parallax_using_geo_par(qn_case1.ctypes.data_as(POINTER(c_double)),
                                   qe_case1.ctypes.data_as(POINTER(c_double)),
                                   t_case1.ctypes.data_as(POINTER(c_double)),
                                   alpha_case1.ctypes.data_as(POINTER(c_double)),
                                   delta_case1.ctypes.data_as(POINTER(c_double)),
                                   tfix_case1.ctypes.data_as(POINTER(c_double)))

    expected_qn_case1 = np.array(0.00048496270844741756, dtype=c_double)
    expected_qe_case1 = np.array(0.01419318, dtype=c_double)

    assert qn_case1 == pytest.approx(expected_qn_case1)
    assert qe_case1 == pytest.approx(expected_qe_case1)


def test_parallax_computation_using_vbbl():
    vbbl = create_vbbl()
    set_parallax_system_for_vbbl(vbbl, c_int(1))
    coordinates_file_path = Path(__file__).parent.joinpath('test_parallax_computation_resources/coordinates.txt'
                                                           ).absolute()
    set_object_coordinates_for_vbbl(vbbl, str(coordinates_file_path).encode('utf-8'), '.'.encode('utf-8'))

    t_case0 = np.array(9033.2148899999447, dtype=c_double)
    t0_case0 = np.array(9100.8972159999994, dtype=c_double)
    et_case0 = np.array([np.nan, np.nan], dtype=c_double)

    compute_parallax_for_vbbl(vbbl,
                              t_case0,
                              t0_case0,
                              et_case0.ctypes.data_as(POINTER(c_double)))

    expected_et_case0 = np.array([5.1377649825789830e-2, 0.52891219662449951], dtype=c_double)

    assert np.allclose(et_case0, expected_et_case0, rtol=0.005)

    t_case1 = np.array(9090.8258499999993, dtype=c_double)
    t0_case1 = np.array(9100.8972159999994, dtype=c_double)
    et_case1 = np.array([np.nan, np.nan], dtype=c_double)

    compute_parallax_for_vbbl(vbbl,
                              t_case1,
                              t0_case1,
                              et_case1.ctypes.data_as(POINTER(c_double)))

    expected_et_case1 = np.array([0.00048496270844741756, 0.01419318], dtype=c_double)

    assert np.allclose(et_case1, expected_et_case1, rtol=0.05)

    destroy_vbbl(vbbl)


def test_geo_par_and_vbbl_calculate_same_parallax():
    time = 9033.2148899999447
    time_fix = 9100.8972159999994

    qn_geo_par_case = np.array(np.nan, dtype=c_double)
    qe_geo_par_case = np.array(np.nan, dtype=c_double)
    t_geo_par_case = np.array(time, dtype=c_double)
    alpha_geo_par_case = np.array(268.43250000000000, dtype=c_double)
    delta_geo_par_case = np.array(-32.589311111111115, dtype=c_double)
    tfix_geo_par_case = np.array(time_fix, dtype=c_double)

    compute_parallax_using_geo_par(qn_geo_par_case.ctypes.data_as(POINTER(c_double)),
                                   qe_geo_par_case.ctypes.data_as(POINTER(c_double)),
                                   t_geo_par_case.ctypes.data_as(POINTER(c_double)),
                                   alpha_geo_par_case.ctypes.data_as(POINTER(c_double)),
                                   delta_geo_par_case.ctypes.data_as(POINTER(c_double)),
                                   tfix_geo_par_case.ctypes.data_as(POINTER(c_double)))

    q_array_geo_par_case = np.array([qn_geo_par_case, qe_geo_par_case])

    vbbl = create_vbbl()
    set_parallax_system_for_vbbl(vbbl, c_int(1))
    coordinates_file_path = Path(__file__).parent.joinpath('test_parallax_computation_resources/coordinates.txt'
                                                           ).absolute()
    set_object_coordinates_for_vbbl(vbbl, str(coordinates_file_path).encode('utf-8'), '.'.encode('utf-8'))

    t_vbbl_case = np.array(time, dtype=c_double)
    t0_vbbl_case = np.array(time_fix, dtype=c_double)
    et_vbbl_case = np.array([np.nan, np.nan], dtype=c_double)

    compute_parallax_for_vbbl(vbbl,
                              t_vbbl_case,
                              t0_vbbl_case,
                              et_vbbl_case.ctypes.data_as(POINTER(c_double)))

    destroy_vbbl(vbbl)

    assert np.allclose(q_array_geo_par_case, et_vbbl_case, rtol=0.005)
