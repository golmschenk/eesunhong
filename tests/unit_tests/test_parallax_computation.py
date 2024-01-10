from ctypes import POINTER, c_double

import numpy as np
import pytest

from eesunhong.eesunhong import compute_parallax_using_geo_par


def test_parallax_computation():
    #  t   9033.2148899999447      a   268.43250000000000      d  -32.589311111111115      f   9100.8972159999994
    #  qn   5.1377649825789830E-002 qe  0.52891219662449951

    #  t   9100.8258499999993      a   268.43250000000000      d  -32.589311111111115      f   9100.8972159999994
    #  qn   8.6430833360821114E-009 qe   7.3889854740073971E-007

    qn0 = np.array(np.nan, dtype=c_double)
    qe0 = np.array(np.nan, dtype=c_double)
    t0 = np.array(9033.2148899999447, dtype=c_double)
    alpha0 = np.array(268.43250000000000, dtype=c_double)
    delta0 = np.array(-32.589311111111115, dtype=c_double)
    tfix0 = np.array(9100.8972159999994, dtype=c_double)

    compute_parallax_using_geo_par(qn0.ctypes.data_as(POINTER(c_double)),
                                   qe0.ctypes.data_as(POINTER(c_double)),
                                   t0.ctypes.data_as(POINTER(c_double)),
                                   alpha0.ctypes.data_as(POINTER(c_double)),
                                   delta0.ctypes.data_as(POINTER(c_double)),
                                   tfix0.ctypes.data_as(POINTER(c_double)))

    expected_qn0 = np.array(5.1377649825789830e-2, dtype=c_double)
    expected_qe0 = np.array(0.52891219662449951, dtype=c_double)

    assert qe0 == pytest.approx(expected_qe0)
    assert qn0 == pytest.approx(expected_qn0)

    qn1 = np.array(np.nan, dtype=c_double)
    qe1 = np.array(np.nan, dtype=c_double)
    t1 = np.array(9100.8258499999993, dtype=c_double)
    alpha1 = np.array(268.43250000000000, dtype=c_double)
    delta1 = np.array(-32.589311111111115, dtype=c_double)
    tfix1 = np.array(9100.8972159999994, dtype=c_double)

    compute_parallax_using_geo_par(qn1.ctypes.data_as(POINTER(c_double)),
                                   qe1.ctypes.data_as(POINTER(c_double)),
                                   t1.ctypes.data_as(POINTER(c_double)),
                                   alpha1.ctypes.data_as(POINTER(c_double)),
                                   delta1.ctypes.data_as(POINTER(c_double)),
                                   tfix1.ctypes.data_as(POINTER(c_double)))

    expected_qn1 = np.array(8.6430833360821114e-9, dtype=c_double)
    expected_qe1 = np.array(7.3889854740073971E-7, dtype=c_double)

    assert qe1 == pytest.approx(expected_qe1)
    assert qn1 == pytest.approx(expected_qn1)