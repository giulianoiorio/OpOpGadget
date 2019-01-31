from __future__ import  division, print_function

from ..densityprofile_src import Density
import numpy as np

class King(Density.Density):

    def __init__(self,c=None,rt=None,rc=None):
        """
        :param:
        :return:
        """

        if c is not None: self.c=rt
        elif (rt is not None) and (rc is not None): self.c=rt/rc
        else: raise ValueError('Neither c or (rt,rc) defined')


    def _evaluatedens(self,x):
        """
        Formula 27 in King,92
        :param x:
        :return:
        """

        z=self._king_z(x)
        z0=self._king_z(0.)
        #d=self._king_functional(z)
        d0=self._king_functional(z0)

        d=self._king_functional(z)

        return d/d0

    def _evaluatesdens(self,x):
        """
        Formula 14 in Kint,92
        :param x:
        :return:
        """
        z=self._king_z(x)
        c=self.c
        amp=1/(1+c*c)

        ret=amp*(1/z-1)**2

        return ret


    def _king_z(self,x):

        c=self.c
        num=1+x*x
        den=1+c*c

        return np.sqrt(num/den)


    def _king_functional(self,z):
        den=z*z
        num=(1/z)*np.arccos(z) - np.sqrt((1-z*z))

        return num/den