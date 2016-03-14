import numpy as np
from scipy.interpolate import UnivariateSpline
from ..model_src import Model
from ..densityprofile_src import Density
from numpy.ctypeslib import ndpointer
import ctypes as ct
import os
from scipy.interpolate import interp1d
from astropy.constants import G as conG

class GeneralModel(Model.Model):

    def __init__(self,R,dens,rc=1,Mmax=1, G='kpc km2 / (M_sun s2)', denorm=True, use_c=False):
        """
        The purpose of the general model is to start from a density law R-dens to build a galaxy model.
        Attenzione per come è creato il modello assume sempre che
        per R>rmax la densita sia 0, la massa resti costante al suo valore massimo e il potenziale vada
        come M/r. Per modelli che raggiungono la massa massima all infinito questo potrebbe essere un problema,
        quindi si dovrebbero usare modelli con massa finita o troncarli e campionarli fino a quanto la massa non raggiunge
        il suo valore max. Per modelli non troncati è meglio utilizzare modelli analitici se possibile.
        Anche nel calcolo del potenziale Rinf è settato uguale all ultimo punto di R, poichè cmq per R>Rmax
        dens=0 e l integrale int_Rmax^inf dens r dr=0 sempre.
        :param R: list of radii, it needs to  be in the form  r/rc
        :param dens: list of dens at radii R. It can be also a function or a lambda function that depends
                     only on the variable R=r/rc
        :param rc: Scale length of the model, the R in input will be multiplyed by rc before start all the calculation
        :param Mmax: Physical Value of the Mass at Rmax (the last point of the R grid). The physical unity of dens and pot and mass
               will depends on the unity of Mmax
        :param G: Value of the gravitational constant G, it can be a number of a string.
                    If G=1, the physical value of the potential will be Phi/G.
                    If string it must follow the rule of the unity of the module.astropy constants.
                    E.g. to have G in unit of kpc3/Msun s2, the input string is 'kpc3 / (M_sun s2)'
                    See http://astrofrog-debug.readthedocs.org/en/latest/constants/
        :param denorm: If True, the output value of mass, dens and pot will be denormalized using Mmax and G.
        :param use_c: To calculate pot and mass with a C-cyle, WARNING it creates more noisy results
        """

        self.rc=rc
        self.Mmax=Mmax
        if isinstance(G,float) or isinstance(G,int): self.G=G
        else:
            GG=conG.to(G)
            self.G=GG.value


        if isinstance(dens,list) or isinstance(dens,tuple) or isinstance(dens,np.ndarray):  self.dens_arr=np.array(dens,dtype=float,order='C')
        elif isinstance(dens,Density.Density): self.dens_arr=dens.dens(R)
        else: self.dens_arr=dens(R)

        self.R=np.array(R,dtype=float,order='C')*self.rc
        self.mass_arr=np.empty_like(self.dens_arr,dtype=float,order='C')
        self.pot_arr=np.empty_like(self.dens_arr,dtype=float,order='C')
        self.use_c=use_c
        self._use_nparray=False

        self._dens=UnivariateSpline(self.R,self.dens_arr, k=1, s=0, ext=1) #for R>rmax, dens=0

        if self.use_c==True:
            #add to path to use relative path
            dll_name='model_c_ext/GeneralModel.so'
            dllabspath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + dll_name
            lib = ct.CDLL(dllabspath)
            #add to path to use relativ path
            mass_func=lib.evalmass
            mass_func.restype=None
            mass_func.argtypes=[ndpointer(ct.c_double, flags="C_CONTIGUOUS"), ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ct.c_int,ndpointer(ct.c_double, flags="C_CONTIGUOUS")]
            mass_func(self.R,self.dens_arr,len(self.dens_arr),self.mass_arr)
            self._mass_int=UnivariateSpline(self.R,self.mass_arr, k=1, s=0, ext=3) #ext=3, const for R>Rmax non ci osno piu particelle e la massa rimane uguale



            pot_func=lib.evalpot
            pot_func.restype=None
            pot_func.argtypes=[ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ct.c_int,ndpointer(ct.c_double, flags="C_CONTIGUOUS")]
            pot_func(self.R,self.dens_arr,self.mass_arr,len(self.dens_arr),self.pot_arr)
            self._pot_int=UnivariateSpline(self.R,self.pot_arr, k=1, s=0, ext=1)



        else:
            self._dm2=UnivariateSpline(self.R,self.R*self.R*self.dens_arr, k=2, s=0,ext=1)
            self._dm=UnivariateSpline(self.R,self.R*self.dens_arr, k=1, s=0,ext=1)


            #Evaluate mass and pot on the R grid in input
            #mass
            func=np.vectorize(self._dm2.integral)
            self.mass_arr=func(0,self.R)

            #pot
            a=(1/self.R)*self.mass_arr
            func=np.vectorize(self._dm.integral)
            b=func(self.R,self.R[-1])
            self.pot_arr=a+b

        if denorm==True: self._set_denorm(self.Mmax)
        else:
            self.Mc=1
            self.dc=1
            self.pc=1

    def _evaluatedens(self,R):
        return self.dc*self._dens(R)

    def _evaluatemass(self,R):
        return self.Mc*self._dm2.integral(0,R)

    def _evaluatemassc(self,R):

        return self.Mc*self._mass_int(R)

    def _evaluatepot(self,R):
            """
            NB specific potential=-Phi
            :param R:
            :return:
            """

            a=(1/R)*self._dm2.integral(0,R)
            b=self._dm.integral(R,self.R[-1])

            return self.pc*(a+b)

    def _evaluatepotc(self,R):
            """
            A differenza di evaluatepot, in questo caso per R>Rmax non
            abbiamo valori, per cui dividiamo in due parti per R<Rmax
            usiamo l 'array calcolato, mentre per R>Rmax, usiamo semplicemente
            Phi=M/R dato che a questi raggi non c è piu materia e tutta la massa
            totale sottesa è l ultimo valore della griglia.
            NB specific potential=-Phi
            :param R:
            :return:
            """

            ret_arr=np.where(R<=self.R[-1],self._pot_int(R),self.mass_arr[-1]/R)
            return self.pc*ret_arr

    def _evaluateradius(self,x,x_type='mass'):

        if x_type=='mass': ret_func=interp1d(self.mass_arr,self.R, kind=linear)
        if x_type=='pot': ret_func=interp1d(self.pot_arr,self.R, kind=linear) #we use this beacuse Univariate spline can have problem if some value on x are equals

        return ret_func(x)

    def _set_denorm(self,Mmax):
        self.Mc=Mmax/self.mass_arr[-1]
        self.dc=self.Mc/(4*np.pi)
        self.pc=self.G*self.Mc
