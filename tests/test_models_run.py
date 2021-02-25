"""
Test that the models run for a short time without error
"""
import pytest

import vorts


NT = 10  # number of time steps to run for


def test_model_f_default_runs():
    m = vorts.Model_f(nt=NT)
    assert m.int_scheme_name == "RK4"
    m.run()


def test_model_f_other_integ_runs():
    m = vorts.Model_f(nt=NT, int_scheme_name="FT")
    m.run()


def test_model_py_default_runs():
    m = vorts.Model_py(nt=NT)  # too many annoying Numba warnings currently
    assert m.int_scheme_name == "RK4"
    m.run()


@pytest.mark.parametrize(
    "int_scheme_name",
    [n for n in vorts.Model_py._allowed_int_scheme_names if n != "RK4"]
)
def test_model_py_other_integ_runs(int_scheme_name):
    m = vorts.Model_py(nt=NT, int_scheme_name=int_scheme_name)
    m.run()
