from __future__ import  division, print_function

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

class King(GeneralModel.GeneralModel):

    def __init__(self, rc, Ms, c=None,rt=None, rtnorm=False, R=None, rini=3e-5,
                 rfin=300, kind='log', n=512, G='kpc km2 / (M_sun s2)', denorm=True, r_physic=False, normalise_tmodel=True, use_c=False,
                 **kwargs):
        """
        King density model (King, 1962)



        with rs=rcpar/c
        and dc depends on c and rho_crit

        :param Mc: Mass at the radius where the  density reach c_par times the critical density of the Univese.
        :param rcore: inner core radius in physical unit
        :param c: Concentration, if None it is calculated with the c-Mvir relation in Munos-Cuartas+10
        :param rt: Truncation radius (in physical unit o in unit of rcpar if rtnorm is True), if None it is equal to 2*rmax(no truncation).
                    if rtnorm is True it is
        :param rtnorm: if true the truncation radius is in unit of rcpar
        :param c_par: Concentration parameter, Mc is the mass at the radius where the density is c_par times the critical density
                of the universe. If c is None c_par is fixed to 200. [200]
        :param z: redshift (used to calculate c) [0]
        :param h: hubble constant (h=H/100, H in unit of  km/s/Mpc ) [0.67]
        :param R: if not None, use this list of normalized radius on rc (if r_physics=False) or non normalised radius uf if r_physics=True
        #Generate grid
        :param rini:  First radius (normalised on rc if r_physics=False)
        :param rfin: Last radius (normalised on rc if r_physics=False)
        :param kind: use a lin or log grid
        :param n: number of points to use to evaluate the density.
        :param G: Value of the gravitational constant G, it can be a number of a string.
                    If G=1, the physical value of the potential will be Phi/G.
                    If string it must follow the rule of the unity of the module.astropy constants.
                    E.g. to have G in unit of kpc3/Msun s2, the input string is 'kpc3 / (M_sun s2)'
                    See http://astrofrog-debug.readthedocs.org/en/latest/constants/
        :param denorm: If True, the output value of mass, dens and pot will be denormalized using Mmax and G.
        :param r_physic: If False all the R, rini and rfin in input ar normalised on rc.
        :param normalise_tmodel: If True and rt=True, normalise the model such that the density is equal for r<<rt to the same
        model without trunctation
        :param use_c: To calculate pot and mass with a C-cyle, WARNING it creates more noisy results
        :param kwargs:
        """


    def __self__(self,rc,Mmax,c=None,rt=None,G='kpc km2 / (M_sun s2)', denorm=True,**kwargs):


        if c is not None: self.c=c
        elif rt is not None: self.c=rt/rc
        else: raise ValueError('Neither c nor rt defined')
        if denorm==True: self._set_denorm(self.Mmax)
        else:
            self.Mc=1
            self.dc=1
            self.pc=1


        if R is None:
            if r_physic:  # rescale to normalised radius
                rini = rini / self.rc
                rfin = rfin / self.rc
            if kind == 'log':
                R = np.logspace(np.log10(rini + 0.01), np.log10(rfin + 0.01), n) - 0.01  # To avoid log(0)
            elif kind == 'lin':
                R = np.linspace(rini, rfin, n)
            self.rini = rini
            self.rfin = rfin
        else:
            if r_physic:
                R = np.asarray(R) / self.rc  # rescale
            else:
                R = np.asarray(R)
            self.rini=R[0]
            self.rfin=R[-1]

	
        
        super(King, self).__init__(R=R, dens=self._adens, rc=self.rc, Ms=Mmax, rs=None, G=G, use_c=use_c, denorm=denorm)


    def _adens(self,x):
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

        return self.dc*d/d0

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

    def _set_denorm(self,Mmax):
        self.Mc=Mmax/self.mass_arr[-1]
        self.dc=self.Mc/(4*np.pi)
        self.pc=self.G*self.Mc


