"""
Integration routines in Python.

The integration functions `vorts.py.integ.integrate_manual` and `vorts.py.integ.integrate_scipy`
are available here (in the `vorts.py` namespace).
"""
from .integ import (  # noqa: F401
    MANUAL_STEPPERS,
    SCIPY_METHODS,
    calc_C,
    integrate_manual,
    integrate_scipy,
)
