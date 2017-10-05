from __future__ import  division, print_function

import numpy as np
from astropy.constants import G as conG
from ..model_src import Model, GeneralModel
from scipy.special import gamma, gammaincc


class Einasto(GeneralModel.GeneralModel):

    def __init__(self,nslope,rc,Mmax,R=None,rini=3e-5,rfin=300,kind='log',n=512,G='kpc km2 / (M_sun s2)',denorm=True,use_c=False):
        """
        Einasto profile (not truncated)
        :param nslope: density exponent, it define the change of slope as function of r (r^(1/n)).
        :param rc: Scale radius
        :param Mmax: Physical Value of the Mass at Rmax (the last point of the R grid). The physical unity of dens and pot and mass
               will depends on the unity of Mmax
        :param R: if not None, use this list of normalized radius on rc.
        #Generate grid
        :param rini:  First radius normalized on rc
        :param rfin: Last normalized radius on rc
        :param kind: use a lin or log grid
        :param n: number of points to use to evaluate the density.
        :param G: Value of the gravitational constant G, it can be a number of a string.
                    If G=1, the physical value of the potential will be Phi/G.
                    If string it must follow the rule of the unity of the module.astropy constants.
                    E.g. to have G in unit of kpc3/Msun s2, the input string is 'kpc3 / (M_sun s2)'
                    See http://astrofrog-debug.readthedocs.org/en/latest/constants/
        :param denorm: If True, the output value of mass, dens and pot will be denormalized using Mmax and G.
        """

        if R is None:
            if kind=='log': R=np.logspace(np.log10(rini+0.01),np.log10(rfin+0.01),n)-0.01 #To avoid log(0)
            elif kind=='lin': R=np.linspace(rini,rfin,n)
        else:
            R=np.asarray(R)

        self.rc=rc
        self.Mmax=Mmax
        self.n=nslope
        self.dn=self._dn()
        self.ninv=1/self.n

        super(Einasto,self).__init__(R=R,dens=self._adens,rc=self.rc,Mmax=Mmax,G=G,use_c=use_c,denorm=denorm)

    def _adens(self,x):
        """
        rho(r)=rc * Exp(-s*dn) where s=(r/rc)^(1/n)
        :param x: normalises  radius
        :return:
        """
        s=x**self.ninv

        sn=s*self.dn

        return np.exp(-sn)


    def _dn(self):
        """
        Retana-Montenegro+12
        :return:
        """
        n=self.n

        c0=0.333333333333333333
        c1=0.006584362139917695
        c2=0.000801271583164587
        c3=0.000033805660996637

        return 3.*n -c0 +c1/n +c2/(n*n) +c3/(n*n*n)








#Einasto teorico va rivisto
class Einastoteo(Model.Model):

    def __init__(self, n, rc, Mmax, G='kpc km2 / (M_sun s2)', denorm=True):
        """
        Einasto density profile: the slope of the density law depends on r
        as r^(1/n). All the functional form for the dens, mass and pot have been
        taken from 'E. Retana-Montenegro et al., 2012, A&A, 540, A70'
        :param n: density exponente, it define the change of slope as function of r (r^(1/n)).
        :param rc: Scale radius
        :param Mmax: Total mass of the Halo
        :param G: Value of the gravitational constant G, it can be a number of a string.
                    If G=1, the physical value of the potential will be Phi/G.
                    If string it must follow the rule of the unity of the module.astropy constants.
                    E.g. to have G in unit of kpc3/Msun s2, the input string is 'kpc3 / (M_sun s2)'
                    See http://astrofrog-debug.readthedocs.org/en/latest/constants/
        :param denorm: If True, the output value of mass, dens and pot will be denormalized using Mmax and G.
        """
        self.rc=rc
        self.Mmax=Mmax
        self.n=n
        self.dn=self._dn()
        self.ninv=1/self.n
        self.gamma3n=gamma(3*self.n)
        self.h=(self.rc)/(self.dn**self.n)


        self.use_c=False
        self._use_nparray=True


        #G
        if isinstance(G,float) or isinstance(G,int): self.G=G
        else:
            GG=conG.to(G)
            self.G=GG.value

        #Denorm
        if denorm==True: self._set_denorm()
        else:
            self.Mc=1
            self.dc=1
            self.pc=1

    def _evaluatedens(self,r):
        """
        rho(r)=rc * Exp(-s^(1/n)) where s=r/h and h=rc/dn^n
        :param r:  radius
        :return:
        """
        s=r/self.h
        sn=s**self.ninv

        return self.dc*np.exp(-sn)

    def _evaluatemass(self,r):
        """
        M(r)=Mc * (1- gamma(3n, s^(1/n))/gamma(3n) )where s=r/h and h=rc/dn^n
        :param r:  radius
        :return:
        """
        s=r/self.h
        sn=s**self.ninv


        a=gammaincc(3*self.n,sn)  #La gamma inc in scipy è normalizzata per la gamma complete

        return self.Mc*(1 - a)

    def _evaluatepot(self,r):
        """
        psi(r)=psic * (1/s) * (1- ( gamma(3n, s^(1/n))/gamma(3n) ) + s*( gamma(2n, s^(1/n))/gamma(3n) )  ) where s=r/h and h=rc/dn^n
        :param r:  radius
        :return:
        """
        s =  self.h / r
        sn=s**self.ninv
        gamma3n=self.gamma3n
        gamma2n=gamma(2*self.n)

        a=(1/s)
        b=gammaincc(3*self.n, sn)  #La gamma inc in scipy è normalizzata per la gamma complete
        c=gammaincc(2*self.n, sn)*(gamma2n/gamma3n)

        return a*(1-b+s*c)

    def _evaluatesdens(self,r):

        raise NotImplementedError('Surface dens not yet implemented for Einasto')

    def _set_denorm(self):

        self.Mc=self.Mmax

        h3=self.h*self.h*self.h
        n=self.n
        gamma3=self.gamma3n
        cost=4*np.pi*h3*n*gamma3
        self.dc=self.Mc/cost


        self.pc=(self.G*self.Mc)/(self.h)

    def _dn(self):
        """
        Retana-Montenegro+12
        :return:
        """
        n=self.n

        c0=0.333333333333333333
        c1=0.006584362139917695
        c2=0.000801271583164587
        c3=0.000033805660996637

        return 3.*n -c0 +c1/n +c2/(n*n) +c3/(n*n*n)