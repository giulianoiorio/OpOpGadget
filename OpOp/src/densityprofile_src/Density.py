from __future__ import  division, print_function

import numpy as np

class Density:

    def __init__(self):
        a=1

    def dens(self,x):
        x=np.asarray(x)
        return self._evaluatedens(x)

    def sdens(self,x):
        x=np.asarray(x)
        return self._evaluatesdens(x)






