from __future__ import  division, print_function

from ..model_src import GeneralModel
import numpy as np
from astropy.constants import G as conG


class Tbetamodel(GeneralModel.GeneralModel):

    def __init__(self,rc,Ms=1,rs=None,gamma=1,beta=3,rt=None,R=None,rini=3e-5,rfin=300,kind='log',n=512,G='kpc km2 / (M_sun s2)',denorm=True,r_physic=False,use_c=False,**kwargs):
        """
        Truncated double power law model:
        dens=dens0 * (r/rc)^(-gamma) * (1+r/rc)^(- (beta-gamma)) * Exp[ -(r/rt)^2]
        It simpy call the class general model with the density law above evaluating it on a grid of radius normalized to rc. This grid
        can be supplied by the user directly or can be generate with the keyword rini,rfin,kind,n.

        :param rc: Scale length
        :param Ms: Physical Value of the Mass at rs.  if rs is None Ms is the value Mmax at last radius R[-1].
                Note: The physical unity of dens and pot and mass
               will depends on the unity of Mmax
        :param rs: Radius at which the mass is equal to Ms in physical unit, if rs is none rs is equalt o the last radius
                    of R.
        :param gamma: first power law exponent
        :param beta: secondo powe law exponent
        :param rt:  Truncation radius (in physical unit), if None it is equal to 3*rmax(no truncation)
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
        :param use_c: To calculate pot and mass with a C-cyle, WARNING it creates more noisy results
        :return:
        """


        self.rc = rc

        if 'Mmax' in kwargs:
            print('Warning keyword Mmax is deprecated for TbetaModel, use instead Ms',flush=True)
            self.Ms=kwargs['Mmax']
        else:
            self.Ms=Ms



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

        if rt is None: self.rt=2*R[-1]*self.rc
        else: self.rt=rt
        if rs is None: self.rs=R[-1]*self.rc
        else: self.rs=rs

        self.gamma=gamma
        self.beta=beta


        super(Tbetamodel,self).__init__(R=R,dens=self._adens,rc=self.rc,Ms=self.Ms,rs=self.rs,G=G,use_c=use_c,denorm=denorm)



    def _adens(self,x):

        y=self.rc/self.rt
        alpha_inn=self.gamma
        alpha_out=self.beta-self.gamma

        dens= ( x**alpha_inn ) * (  (1+x)**alpha_out   )

        return (1./dens)*np.exp(-x*x*y*y)

    def __str__(self):

        h=''
        h+='\nModel: TbetaModel'
        h+='\ngamma: %.2f'%self.gamma
        h+='\nbeta: %.2f'%self.beta
        h+='\nrc: %.3f'%self.rc
        h+='\nrt: %.3f (physical)  %.3f (normalised)'%(self.rt,self.rt/self.rc)
        h+='\nrini: %.3f (physical)  %.3f (normalised)'%(self.rini*self.rc,self.rini)
        h+='\nrfin: %.3f (physical)  %.3f (normalised)'%(self.rfin*self.rc,self.rfin)
        h+='\nMass: %.3e at scale radius rs: %.3f'%(self.Ms,self.rs)
        h+='\nTotal Mass: %.3e at last radius: %.3f'%(self.mass(self.rfin*self.rc),self.rfin*self.rc)
        h+='\nuse_c set to %s'%str(self.use_c)
        h+='\nuse_nparray set to %s' % str(self._use_nparray)
        h+='\n'

        return h


class NFW(Tbetamodel):

    def __init__(self,Mc,c=None,  rt=None, rtnorm=False, c_par=200, z=0, h=0.67, R=None,rini=3e-5,rfin=300,kind='log',n=512,G='kpc km2 / (M_sun s2)',denorm=True,r_physic=False,use_c=False,**kwargs):
        """
        NFW density models:
        d=dc/( (r/rs)*(1+r/rs) )*np.exp(-r*r/rt*rt)

        with rs=rcpar/c
        and dc depends on c and rho_crit

        :param Mc: Mass at the radius where the  density reach c_par times the critical density of the Univese.
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
        :param use_c: To calculate pot and mass with a C-cyle, WARNING it creates more noisy results
        :param kwargs:
        """


        self.Mc=Mc

        #cosmo. par:
        self.z=z #redshift
        self.h=h #hubble constant a z0
        self.c_par=c_par #concentration parameter

        if c is None:
            #If c is None, use the c-Mvir, cosmological relation in Munos-Cuartas+10,
            #but in this cas the concentration parameter is fixed to 200.
            if c_par!=200:
                print('Warning, you are using a cvir-M200 relation, but your conc. parameter is %.3f, the new conc is fixed to 200'%c_par)
                self.c_par=200

            self.c=self.cvir()

        else:
            self.c=c

        #Other
        self.rhoc=self.rho_critic() #Rho critic
        self.denc=self.rhoc*self.c_par #value of the density at Rcp (R200)
        self.rcp=self.rcpar() #Radius where M=Mc
        self.rscale=self.rcpar()/self.c #scale radius
        self.rho_s=self.rho0()



        if rt is None: pass
        elif rtnorm:
            #if rt<1.1: print('Warning the truncation radius is lower that 1.1 the virial radius of the halo')
            rt=self.rcp*rt
        else:
            #if rt<1.1*self.rcp: print('Warning the truncation radius is lower that 1.1 the virial radius of the halo')
            pass


        #call Tbeta-Model
        super(NFW,self).__init__(rc=self.rscale,Ms=self.Mc,rs=self.rcp,gamma=1,beta=3,rt=rt,R=R,rini=rini,rfin=rfin,kind=kind,n=n,G=G,denorm=denorm,r_physic=r_physic,use_c=use_c,**kwargs)


    def delta_c(self):

        c=self.c

        a=np.log(1+c)
        b=c/(1+c)

        return 1/(a+b)

    def rho0(self):

        c=self.c
        gc=self.delta_c()

        return (self.denc*gc*c**3)/3

    def rcpar(self): #r200 ic c_par=200, is the radius at which M=Mc

        cost=3./4.
        a=(np.pi*self.denc)

        return (cost*self.Mc/a)**(1./3.)

    def cvir(self):
        """
        From Munos-Cuartas+10
        :return:
        """
        z=self.z
        h=self.h
        mvir=self.Mc


        x=np.log10(mvir*h) #In the function M is in unit of Msun/h

        #fit constant
        w=0.029
        m=0.097
        alpha=-110.001
        beta=2469.720
        gamma=16.885

        a=w*z-m
        b= ( alpha/(z+gamma) ) + ( beta/(z+gamma)**2 )

        log10c=a*x+b

        return 10**log10c

    def rho_critic(self):
        """
        Critic density of the Universe given h=H/100
        :return:
        """
        cost= (3/(8*np.pi))
        G=conG.to('kpc km2 / (M_sun s2)')
        H=self.h*100*1e-3 #in km/s Mpc

        num=cost*H*H

        return num/G.value


    def __str__(self):

        h=''
        h+='\nModel: NFW'
        h+='\nconcentration parameter: %i'%self.c_par
        h+='\nCosm parameters,  z:%.3f   h:%.3f  rhoc:%.3f'%(self.z,self.h,self.rhoc)
        h+='\nc: %.3f'%self.c
        h+='\nM%i: %.3e'%(self.c_par,self.Ms)
        h+='\nr%i: %.3f'%(self.c_par,self.rcp)
        h+='\nrs: %.3f'%self.rscale
        h+='\nrt: %.3f (physical)  %.3f (normalised)'%(self.rt,self.rt/self.rs)
        h+='\nrini: %.3f (physical)  %.3f (normalised)'%(self.rini*self.rscale,self.rini)
        h+='\nrfin: %.3f (physical)  %.3f (normalised)'%(self.rfin*self.rscale,self.rfin)
        h+='\nTotal Mass: %.3e at last radius: %.3f'%(self.mass(self.rfin*self.rscale),self.rfin*self.rscale)
        h+='\nuse_c set to %s'%str(self.use_c)
        h+='\nuse_nparray set to %s' % str(self._use_nparray)
        h+='\n'

        return h

class Hernquist(Tbetamodel):

    def __init__(self,rc,Ms=1,rs=None,rt=None,R=None,rini=3e-5,rfin=300,kind='log',n=512,G='kpc km2 / (M_sun s2)',denorm=True,r_physic=False,use_c=False,**kwargs):
        """
        Hernquist:
        dens=dens0 * (r/rc)^(-1) * (1+r/rc)^(- (3)) * Exp[ -(r/rt)^2]
        It simpy call the class general model with the density law above evaluating it on a grid of radius normalized to rc. This grid
        can be supplied by the user directly or can be generate with the keyword rini,rfin,kind,n.

        :param rc: Scale length
        :param Ms: Physical Value of the Mass at rs.  if rs is None Ms is the value Mmax at last radius R[-1].
                Note: The physical unity of dens and pot and mass
               will depends on the unity of Mmax
        :param rs: Radius at which the mass is equal to Ms in physical unit, if rs is none rs is equalt o the last radius
                    of R.
        :param gamma: first power law exponent
        :param beta: secondo powe law exponent
        :param rt:  Truncation radius (in physical unit), if None it is equal to 3*rmax(no truncation)
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
        :param use_c: To calculate pot and mass with a C-cyle, WARNING it creates more noisy results
        :return:
        """


        self.rc = rc

        if 'Mmax' in kwargs:
            print('Warning keyword Mmax is deprecated for TbetaModel, use Ms instead',flush=True)
            self.Ms=kwargs['Mmax']
        else:
            self.Ms=Ms



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

        if rt is None: self.rt=2*R[-1]*self.rc
        else: self.rt=rt
        if rs is None: self.rs=R[-1]*self.rc
        else: self.rs=rs


        self.rscale=rc

        super(Hernquist, self).__init__(rc=self.rscale,Ms=self.Mc,rs=self.rcp,gamma=1,beta=4,rt=rt,R=R,rini=rini,rfin=rfin,kind=kind,n=n,G=G,denorm=denorm,r_physic=r_physic,use_c=use_c,**kwargs)


    def __str__(self):

        h=''
        h+='\nModel: Henquist'
        h+='\ngamma: %.2f'%self.gamma
        h+='\nbeta: %.2f'%self.beta
        h+='\nrscale: %.3f'%self.rc
        h+='\nrt: %.3f (physical)  %.3f (normalised)'%(self.rt,self.rt/self.rc)
        h+='\nrini: %.3f (physical)  %.3f (normalised)'%(self.rini*self.rc,self.rini)
        h+='\nrfin: %.3f (physical)  %.3f (normalised)'%(self.rfin*self.rc,self.rfin)
        h+='\nMass: %.3e at scale radius rs: %.3f'%(self.Ms,self.rs)
        h+='\nTotal Mass: %.3e at last radius: %.3f'%(self.mass(self.rfin*self.rc),self.rfin*self.rc)
        h+='\nuse_c set to %s'%str(self.use_c)
        h+='\nuse_nparray set to %s' % str(self._use_nparray)
        h+='\n'

        return h