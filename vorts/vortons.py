"""
Vorton class and Vortons container class (TODO).
"""
import warnings

import numpy as np

class Vorton():
    def __init__(self, G, x, y, nt):
        """Create vorton.

        x, y inputs can be either a single point (IC)
        or an array
        so that same class can be used for catching/organizing the Fortran model output
        """
#        self.x = xi
#        self.y = yi
        self.G = G

        try:
            self.xhist = np.zeros(nt+1)
            self.yhist = np.zeros(nt+1)
            self.xhist[0] = x
            self.yhist[0] = y

        except ValueError:  # ValueError: setting array with seq; TypeError: ints and floats have no len()
            self.xhist = x
            self.yhist = y

            if self.xhist.size != nt+1:
                warnings.warn(
                    f"Fortran model output should be size {nt+1:d} but is {self.xhist.size:d}"
                )
