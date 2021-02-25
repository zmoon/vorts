"""
Test that the plotting functions work without error
"""
import pytest

import vorts

# Run case to plot
m = vorts.Model_py(tracers=vorts.Tracers.grid(7, 7)).run()


@pytest.mark.parametrize(
    "which",
    ("vortons", "tracers", "poincare")
)
def test_plot_works_from_model(which):
    m.plot(which)
