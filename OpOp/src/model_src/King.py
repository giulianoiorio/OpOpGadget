import numpy as np
from scipy.interpolate import UnivariateSpline
from ..model_src import Model
from ..densityprofile_src import Density
from numpy.ctypeslib import ndpointer
import ctypes as ct
import os
from scipy.interpolate import interp1d
from astropy.constants import G as conG


#Bisogna mettere integrazione

class King(Model.Model):


    def __self__(self,c=None,rt=None,rc=None,Mmax=1,G='kpc km2 / (M_sun s2)', denorm=True,**kwargs):


        if c is not None: self.c=rt
        elif (rt is not None) and (rc is not None): self.c=rt/rc
        else: raise ValueError('Neither c or (rt,rc) defined')

        if denorm==True: self._set_denorm(self.Mmax)
        else:
            self.Mc=1
            self.dc=1
            self.pc=1

    def _evaluatedens(self,x):
        """
        Formula 29 in King,92
        :param x:
        :return:
        """

        z=self._king_z(x)
        z0=self._king_z(0.)
        #d=self._king_functional(z)
        d0=self._king_functional(z0)

        d=np.where(z<=1,self._king_functional(z),0.)

        return self.dc*d/d0

    def _evaluatesdens(self,x):
        """
        Formula 14 in Kint,92
        :param x:
        :return:
        """

        c=self.c

        b=1/np.sqrt((1+c*c))
        a=1/np.sqrt((1+x*x))

        ret=np.where(x<c,(a/b-1)**2,0.)

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

    def _set_denorm(self,Mmax):
        self.Mc=Mmax/self.mass_arr[-1]
        self.dc=self.Mc/(4*np.pi)
        self.pc=self.G*self.Mc


