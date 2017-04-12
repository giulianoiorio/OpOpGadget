import numpy as np
from astropy.constants import G as conG
from ..model_src import Model
from scipy.special import gamma,gammaincc

class Sersic(Model.Model):

    def __init__(self,m,rc,Mtot,G='kpc km2 / (M_sun s2)',denorm=True):
        """
        Analytic Sersic model:

        Surface density: S(R)=S(0) * Exp(-(R/rc)^(1/m))
        for the 3D deprojection we use the analytic approximation of Lima Neto, Gerbal and Marquez, 1999.
        Density: d(r)=d(0) * Exp(-(r/rc)^(1/m)) * (r/rc)^(-p) where p=1 - 0.6097*(1/m) + 0.05463*(1/m^2)

        :param m: Sersic exponent
        :param rc: Sersic scale length
        :param Mtot:  Sersic total mass
        :param G: Value of the gravitational constant G, it can be a number of a string.
                    If G=1, the physical value of the potential will be Phi/G.
                    If string it must follow the rule of the unity of the module.astropy constants.
                    E.g. to have G in unit of kpc3/Msun s2, the input string is 'kpc3 / (M_sun s2)'
                    See http://astrofrog-debug.readthedocs.org/en/latest/constants/
        :param denorm: If True, the output value of mass, dens and pot will be denormalized using Mmax and G.
        :return:
        """

        self.rc=rc
        self.Mmax=Mtot
        self.m=m
        self.nu=1/m
        self.p=1-0.6097*self.nu+0.05463*self.nu*self.nu
        self.use_c=False
        self._use_nparray=True
        self._analytic_radius=True

        if isinstance(G,float) or isinstance(G,int): self.G=G
        else:
            GG=conG.to(G)
            self.G=GG.value

        #set denorm
        if denorm==True: self._set_denorm(self.Mmax)
        else:
            self.Mc=1
            self.sdc=1
            self.dc=1
            self.pc=1


    def _set_denorm(self,Mmax):
        rc=self.rc
        dmgamma= self.m*gamma(self.m*(3-self.p)) #Mtot=4*pi*d0*rc^3 * m* *Gamma(m*(3-p))
        self.dc=Mmax/(4*np.pi*rc*rc*rc*dmgamma)
        sdmgamma= self.m*gamma(2*self.m)  #Mtot=2*pi*S0*rc^2 * m * Gamma(2*m)
        self.sdc=Mmax/(2*np.pi*rc*rc*sdmgamma)
        self.pc=(self.G*Mmax)
        self.Mc = Mmax



    def _evaluatedens(self,R):
        """
        Formula=d(r)=d(0) * Exp(-(r/rc)^(1/m)) * (r/rc)^(-p) where p=1 - 0.6097*(1/m) + 0.05463*(1/m^2)
        :param R:
        :return:
        """
        R=np.asarray(R)
        a=(R/self.rc)**(-self.p)
        b=np.exp( -(R/self.rc)**self.nu )

        return self.dc*a*b

    def _evaluatesdens(self,R):
        """
        formula: S(R)=S(0) * Exp(-(R/rc)^(1/m))
        :param R:
        :return:
        """
        R=np.asarray(R)
        b=np.exp( -(R/self.rc)**self.nu )

        return self.sdc*b

    def _evaluatemass(self,R):
        """
        formula: M(r)= Mtot * (1- Gamma(m(3-p),r/rc**(1/m))/Gamma(m(3-p))
        :param R:
        :return:
        """
        R=np.asarray(R)
        num=gammaincc(self.m*(3-self.p),(R/self.rc)**self.nu) #Non mi serve dividere per Gamma(2m), perchè in scipy la gammaincc è gia normalizzta
                                          #su gamma

        return self.Mc*(1-num)

        #return self.Mmax*(num/den)

    def _evaluatepot(self,R):
        """
        Pot= G*( M(r)/r + (Mtot/rc)*gamma(m(p-2),R^1/m)/gamma(m(p-3)   )
        ma M(r)= Mtot * (1- Gamma(m(3-p),r**(1/m))/Gamma(m(3-p))
        quindi

        Pot= (GMtot/rc)*((rc/r)*((1- Gamma(m(3-p),r/rc**(1/m))/Gamma(m(3-p)) + gamma(m(2-p),R^1/m)/gamma(m(3-p) )

        :param R:
        :return:
        """
        m=self.m
        p=self.p
        rc=self.rc
        x=m*(3-p)
        y=m*(2-p)

        R=np.asarray(R)

        a=(1/R)*(1-gammaincc(x,(R/rc)**(1/m)))

        b1=gammaincc(y,(R/rc)**(1/m))*gamma(y)
        b2=gamma(x)

        b=(1/rc)*b1/b2


        return self.pc*( a+b)

