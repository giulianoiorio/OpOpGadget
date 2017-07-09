from ..model_src import GeneralModel
import numpy as np
from astropy.constants import G as conG




class NFWc(GeneralModel.GeneralModel):

    def __init__(self, Mc, rcore, mode='h', c=None, rt=None, rtnorm=False, c_par=200, z=0, h=0.67, R=None, rini=3e-5,
                 rfin=300, kind='log', n=512, G='kpc km2 / (M_sun s2)', denorm=True, r_physic=False, normalise_tmodel=True, use_c=False,
                 **kwargs):
        """
        NFW density models with an inner core
        mode=s:
            d=dc/( (rcore/rc + r/rs)*(1+r/rs)^2 )*np.exp(-r*r/rt*rt)
        mode=h:
            d=dc/( ((rcore^2/rc^2 + r^2/rs^2)*0.5)*(1+r/rs)^2 )*np.exp(-r*r/rt*rt)
        mode=e:
            d=dc/( rcore/rs*exp(-r/rcore) + (r/rs)*(1+r/rs)^2  )*np.exp(-r*r/rt*rt)
        mode=a:
            d=dc/( ((rcore^2/rc^2 + r^2/rs^2)*0.5)*(1+r^2/rs^2)  )*np.exp(-r*r/rt*rt)

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

        self.Mc = Mc
        self.rcore = rcore

        self.mode=mode

        # cosmo. par:
        self.z = z  # redshift
        self.h = h  # hubble constant a z0
        self.c_par = c_par  # concentration parameter

        if c is None:
            # If c is None, use the c-Mvir, cosmological relation in Munos-Cuartas+10,
            # but in this cas the concentration parameter is fixed to 200.
            if c_par != 200:
                print(
                    'Warning, you are using a cvir-M200 relation, but your conc. parameter is %.3f, the new conc is fixed to 200' % c_par)
                self.c_par = 200

            self.c = self.cvir()

        else:
            self.c = c

        # Other
        self.rhoc = self.rho_critic()  # Rho critic
        self.denc = self.rhoc * self.c_par  # value of the density at Rcp (R200)
        self.rcp = self.rcpar()  # Radius where M=Mc
        self.rscale = self.rcpar() / self.c  # scale radius
        self.rho_s = self.rho0()

        #Scale ratio
        self.eps= self.rcore  / self.rscale



        if R is None:
            if r_physic:  # rescale to normalised radius
                rini = rini / self.rscale
                rfin = rfin / self.rscale
            if kind == 'log':
                R = np.logspace(np.log10(rini + 0.01), np.log10(rfin + 0.01), n) - 0.01  # To avoid log(0)
            elif kind == 'lin':
                R = np.linspace(rini, rfin, n)
            self.rini = rini
            self.rfin = rfin
        else:
            if r_physic:
                R = np.asarray(R) / self.rscale  # rescale
            else:
                R = np.asarray(R)
            self.rini=R[0]
            self.rfin=R[-1]



        if rt is None:
            self.rt=2*R[-1]
        elif rtnorm:
            if rt < 1.1: print('Warning the truncation radius is lower that 1.1 the virial radius of the halo')
            self.rt = self.rcp * rt
        else:
            if rt < 1.1 * self.rcp: print(
                'Warning the truncation radius is lower that 1.1 the virial radius of the halo')
            self.rt = rt



        # call Tbeta-Model


        if (rt is not None) and (normalise_tmodel==True):
            rtscale=0.1*self.rt
            super(NFWc, self).__init__(R=R, dens=self._adens_not, rc=self.rscale, Ms=self.Mc, rs=self.rcp, G=G, use_c=use_c,
                                        denorm=denorm)
            #print('Mass NFWc no trunc', self.mass(0.1 * self.rt))
            mass_scale=self.mass(rtscale)
            super(NFWc, self).__init__(R=R, dens=self._adens, rc=self.rscale, Ms=mass_scale, rs=rtscale, G=G, use_c=use_c,
                                        denorm=denorm)
            #print('Mass NFWc fin', self.mass(0.1 * self.rt))
        else:
            super(NFWc, self).__init__(R=R, dens=self._adens, rc=self.rscale, Ms=self.Mc, rs=self.rcp, G=G, use_c=use_c,
                                       denorm=denorm)

            #print('Mass NFWc_trunc', self.mass(0.1 * self.rt))

    def _adens(self,x):

        eps=self.eps

        if self.mode=='s':

            a=(eps+x)
            b=(1+x)*(1+x)

            dens=1/(a*b)

        elif self.mode=='h':

            a=(eps*eps+x*x)**(0.5)
            b = (1 + x) * (1 + x)

            dens=1/(a*b)

        elif self.mode=='e':

            y=x/eps

            a=eps*np.exp(-y)
            b=x*(1+x)*(1+x)

            dens=1/(a+b)

        elif self.mode=='a':

            a = (eps * eps + x * x) ** (0.5)
            b = (1+x*x)

            dens=1 / (a * b)

        else:
            raise NotImplementedError('mode %s not implementend in model NFWc'%self.mode)

        etrunc=self.rscale / self.rt

        return dens*np.exp(-x*x*etrunc*etrunc)


    def _adens_not(self,x):

        eps=self.eps

        if self.mode=='s':

            a=(eps+x)
            b=(1+x)*(1+x)

            dens=1/(a*b)

        elif self.mode=='h':

            a=(eps*eps+x*x)**(0.5)
            b = (1 + x) * (1 + x)

            dens=1/(a*b)

        elif self.mode=='e':

            y=x/eps

            a=eps*np.exp(-y)
            b=x*(1+x)*(1+x)

            dens=1/(a+b)

        elif self.mode=='a':

            a = (eps * eps + x * x) ** (0.5)
            b = (1+x*x)

            dens=1 / (a * b)

        else:
            raise NotImplementedError('mode %s not implementend in model NFWc'%self.mode)


        return dens



    def delta_c(self):

        c = self.c

        a = np.log(1 + c)
        b = c / (1 + c)

        return 1 / (a + b)

    def rho0(self):

        c = self.c
        gc = self.delta_c()

        return (self.denc * gc * c ** 3) / 3

    def rcpar(self):  # r200 ic c_par=200, is the radius at which M=Mc

        cost = 3. / 4.
        a = (np.pi * self.denc)

        return (cost * self.Mc / a) ** (1. / 3.)

    def cvir(self):
        """
        From Munos-Cuartas+10
        :return:
        """
        z = self.z
        h = self.h
        mvir = self.Mc

        x = np.log10(mvir * h)  # In the function M is in unit of Msun/h

        # fit constant
        w = 0.029
        m = 0.097
        alpha = -110.001
        beta = 2469.720
        gamma = 16.885

        a = w * z - m
        b = (alpha / (z + gamma)) + (beta / (z + gamma) ** 2)

        log10c = a * x + b

        return 10 ** log10c

    def rho_critic(self):
        """
        Critic density of the Universe given h=H/100
        :return:
        """
        cost = (3 / (8 * np.pi))
        G = conG.to('kpc km2 / (M_sun s2)')
        H = self.h * 100 * 1e-3  # in km/s Mpc

        num = cost * H * H

        return num / G.value

    def __str__(self):

        if self.mode=='s': mode='soft'
        elif self.mode=='h': mode='hard'
        elif self.mode=='e': mode='exp'
        elif self.mode=='a': mode='analog'

        h = ''
        h += '\nModel: NFWc'
        h += '\ncore mode:%s'%mode
        h += '\nconcentration parameter: %i' % self.c_par
        h += '\nCosm parameters,  z:%.3f   h:%.3f  rhoc:%.3f' % (self.z, self.h, self.rhoc)
        h += '\nc: %.3f' % self.c
        h += '\nM%i: %.3e' % (self.c_par, self.Ms)
        h += '\nr%i: %.3f' % (self.c_par, self.rcp)
        h += '\nrs: %.3f' % self.rscale
        h += '\nrcore: %.3f (physical)  %.3f (normalised)'%(self.rcore,self.rcore/self.rscale)
        h += '\nrt: %.3f (physical)  %.3f (normalised)' % (self.rt, self.rt / self.rs)
        h += '\nrini: %.3f (physical)  %.3f (normalised)' % (self.rini * self.rscale, self.rini)
        h += '\nrfin: %.3f (physical)  %.3f (normalised)' % (self.rfin * self.rscale, self.rfin)
        h += '\nTotal Mass: %.3e at last radius: %.3f' % (self.mass(self.rfin * self.rscale), self.rfin * self.rscale)
        h += '\nuse_c set to %s' % str(self.use_c)
        h += '\nuse_nparray set to %s' % str(self._use_nparray)
        h += '\n'

        return h