"""
Test that the models run without error
"""
import vorts


def test_model_f_runs():
    m = vorts.Model_f(tracers=vorts.Tracers.randu(n=10))
    m.run()


def test_model_py_runs():
    m = vorts.Model_py()
    m.run()
