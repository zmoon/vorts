"""
Tests transforms for Vortons/Tracers
"""
import numpy as np

import vorts


def test_add():
    t0 = vorts.Tracers([1], [0])
    np.testing.assert_equal((t0 + (0, 1)).xy, np.c_[1, 1])


def test_scale():
    t0 = vorts.Tracers([2], [1])
    np.testing.assert_equal((1.5 * t0).xy, np.c_[3, 1.5])


def test_rotate():
    t0 = vorts.Tracers([1], [0])
    np.testing.assert_almost_equal(t0.rotate(90).xy, np.c_[0, 1])
