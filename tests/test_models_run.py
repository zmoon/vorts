"""
Test that the models run without error
"""
import vorts


def test_model_f_default_runs():
    m = vorts.Model_f(tracers=vorts.Tracers.randu(n=10))
    m.run()


def test_model_py_default_runs():
    m = vorts.Model_py(int_scheme_name="scipy_DOP853")  # too many annoying Numba warnings currently
    m.run()
