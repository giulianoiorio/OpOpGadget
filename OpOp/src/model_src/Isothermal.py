from ..model_src import Model
from astropy.constants import G as conG
from fermi.tsintegrator import Tsintegrator1D
import numpy as np

class Isothermal(Model.Model):

    def __init__(self,rc,  G='kpc km2 / (M_sun s2)', Mmax=None, rmax=None, Vinf=None, d0=None):
        """
        PseudoIsothermal Model:

        d=d0/(sqrt(1+r^2/rc^2))

        Given that the Mass diverge at infinty it is possibile to initialize the
        model in different ways:
        Mmax  (rmax): The model will have a mass of Mmax at radius rmax, if rmax is not supplied
                      it is equal to 10*rc
        d0: Central density. The unity depends on the combinations of gthe value of G and rc. By default
                      is Msun/kpc^3.
        Vinf: asymptotic circular velocity in km/s.

        The routine gives priority firstly to Mmax, then to d0 and finally to Vinf.

        :param rc:
        :param G: Value of the gravitational constant G, it can be a number of a string.
                    If G=1, the physical value of the potential will be Phi/G.
                    If string it must follow the rule of the unity of the module.astropy constants.
                    E.g. to have G in unit of kpc3/Msun s2, the input string is 'kpc3 / (M_sun s2)'
                    See http://astrofrog-debug.readthedocs.org/en/latest/constants/
        :param Mmax: Total Mass in Msun at radius rmax.
        :param rmax: Radius to cut the density distribution, by default is equal to 10*rc.
        :param d0: Central density. The unity depends on the combinations of gthe value of G and rc. By default
                   is Msun/kpc^3. Available only if Mmax is None.
        :param Vinf: asymptotic circular velocity in km/s. Available only if Mmax and d0 are None.
        :return:
        """

        if rmax is None: self.rmax=10*rc
        else: self.rmax=rmax

        self.rc=rc
        self.use_c=False
        self._use_nparray=True
        if isinstance(G,float) or isinstance(G,int): self.G=G
        else:
            GG=conG.to(G)
            self.G=GG.value

        self.Mc=1
        self.dc=1
        self.pc=1
        if Mmax is not None:
            totmass=self._evaluatemass(self.rmax)
            self.Mc=Mmax/totmass
            self.dc=self.Mc/(4*np.pi)
            self.pc=self.G*self.Mc
        elif d0 is not None:
            self.dc=d0
            self.Mc=(4*np.pi)*self.dc
            self.pc=self.G*self.Mc
        elif Vinf is not None:
            self.pc=(Vinf*Vinf)/(self.rc*self.rc)
            self.Mc=self.pc/self.G
            self.dc=self.Mc/(4*np.pi)

    def _evaluatesdens(self,r,rcut=None):
        """
        This is the results of a in integration between R and Rcut. I
        :param r:
        :param rcut: If None Rcut is equals to Rmax. If equals to 'inf' the funtcion returns
                    the exact analytic solutions for a infinite isothermal sphere.
        :return:
        """

        x=r/self.rc

        d1=2*self.rc/(np.sqrt(1+x*x))

        if rcut is None:
            rcut=self.rmax
            num=rcut*rcut-r*r
            den=1+x*x
            d2=np.arctan( (1/self.rc)*(np.sqrt(num/den))  )
        elif rcut=='inf':
            d2=np.pi/2.
        else:
            num=rcut*rcut-r*r
            den=1+x*x
            d2=np.arctan( (1/self.rc)*(np.sqrt(num/den))  )

        return self.dc*d1*d2

    def _evaluatedens(self,r):

        x=r/self.rc

        return self.dc/(1+x*x)

    def _evaluatemass(self,r):

        x=r/self.rc
        rc3=self.rc**3

        return self.Mc*(rc3*(x-np.arctan(x)))

    def _evaluatepot(self,r,rcut=None):

        if rcut is None: rcut=self.rmax

        x=r/self.rc
        rc3=self.rc**3
        rc2=self.rc*self.rc

        p1=(rc3*(a-np.arctan(x)))/r
        p2=0.5*rc2 * (np.log(rc2+rcut*rcut) - np.log(rc2*r*r) )

        return self.pc*(p1+p2)

