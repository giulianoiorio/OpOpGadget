from scipy.integrate import quad
from math import sqrt
from astropy.constants import G as conG
import fermi.tsintegrator as ft
from scipy.interpolate import UnivariateSpline
from .jsolver_c_ext import CJsolver
from ..model_src.Model import Model, Multimodel
import numpy as np

class Jsolver():

    def __init__(self,dprof,sprof,mprof,amodel='isotropic',G='kpc km2 / (M_sun s2)'):
        """
        funnzioni , dprof, sprof, mprof
        """
        self.amodel=amodel



        self.dprof=dprof
        self.sprof=sprof
        self.mprof=mprof

        if isinstance(G,float) or isinstance(G,int): self.G=G
        else:
            GG=conG.to(G)
            self.G=GG.value

        if amodel=='isotropic':
            self.kernel=self._kerneliso
        else:
            raise NotImplementedError('Amodel %s not implemented'%amodel)


        self.cost=sqrt(2.*self.G)

    @classmethod
    def model(cls,*args,tracer=0,amodel='isotropic',G='kpc km2 / (M_sun s2)'):



        #Check
        for mod in args:
            if isinstance(mod,Model): pass
            else: raise ValueError('Model is not an instance of Class Model')

        #Define
        if len(args)>1:
            dprof = args[tracer].dens
            sprof = args[tracer].sdens
            mmod=Multimodel(*args)
            mprof=mmod.mass
        else:
            dprof = args[0].dens
            sprof = args[0].sdens
            mprof = args[0].mass

        return cls(dprof=dprof,sprof=sprof,mprof=mprof,amodel=amodel,G=G)

    def _kproj_iso(self,r,R):
        """
        Mamon and Lokas, 2005, A16-isotropic
        :param r: spherical radius
        :param R: cylindircal radius
        :return: kproj
        """
        return np.sqrt(1-(R*R)/(r*r))

    def _intkernel(self,r):
        """

        :param r: spherical radius
        :return:
        """

        return self.dprof(r)*self.mprof(r)/r

    def _kerneliso(self,r,R):

        return self._kproj_iso(r,R)*self._intkernel(r)

    def integrate(self,R,Rcut):
        return np.sqrt((2*self.G)*(1/self.sprof(R))*quad(self.kernel,R,Rcut,args=(R,))[0])

    def integratec(self,R=None,rini=0.01,rfin=10,n=512,kind='log',cpar=True):

        ret=CJsolver._integratec(self.G,self.dprof,self.sprof,self.mprof,R=R,rini=rini,rfin=rfin,n=n,kind=kind,cpar=cpar)

        return ret

    def integrateT(self,R,Rcut,N=20,hstep=3.):
        f=ft.Tsintegrator1D(N,hstep=hstep)
        a=f.integrate(self._kerneliso,R,Rcut,extra_args=(R,))


        return sqrt((2*self.G)*(1/self.sprof(R))*a)

    def vdisp(self,R,Rcut=None,use_c=True,mode='F',n=512,kind='log',nt=20):
        """
        Calculate the velocity dispersion for a series of point R.
        :param R: float, int or array. NB Do not use value <0
        :param Rcut: Cutting radius for the integration, by default is the last R*5.
        :param use_c: If True use the c for-loop
        :param mode: N-Use the numerical trapezoidal integration F-Use the fast tanhsinh numerical integration S-Use the slow adaptive quad integration
        :param n: if mode=N, the point to use to generate the grid of density and mass values
        :param kind: if mode=N, the mode (lin or log) to generate the grid
        :param nt: if mode=F, the functional evaluation of the tanh-sinh method
        :return:
        """



        if isinstance(R,float) or isinstance(R,type(1)):

            if Rcut is None: Rcut=R*5
            if mode=='N' or mode=='F':
                ret=self.integrateT(R,Rcut,N=nt)
            else:
                ret=self.integrate(R,Rcut)


        else:
            R=np.array(R,dtype=np.float64)
            if Rcut is None: Rcut=R[-1]*5

            if mode=='N':
                rr,ret_tmp=self.integratec(rini=R[0],rfin=Rcut,n=n,kind=kind)
                f=UnivariateSpline(rr,ret_tmp,s=0)
                ret=f(R)

            elif use_c:

                if mode=='F': f=self.integrateT
                elif mode=='S': f=self.integrate

                ret=CJsolver._evaluate(f,R,Rcut)

            else:

                 if mode=='F': int=np.vectorize(self.integrateT,excluded=['Rcut','N'])
                 elif mode=='S': int=np.vectorize(self.integrate,excluded=['Rcut','N'])
                 ret=int(R,Rcut)

        return ret
