from ..model_src import GeneralModel
from ..model_src import Model
import numpy as np
from astropy.constants import G as conG


class Plummer2D(Model.Model):

    def __init__(self,rc,Mmax,R=None,rini=3e-5,rfin=300,kind='log',n=512,G='kpc km2 / (M_sun s2)',denorm=True):


        self.rc=rc
        self.Mmax=Mmax
        self.use_c=False
        self.use_nparray=True

        if isinstance(G,float) or isinstance(G,int): self.G=G
        else:
            GG=conG.to(G)
            self.G=GG.value


        #Select grid for mass e pot inversion
        if R is None:
            if kind=='log': self.R=np.logspace(np.log10(rini+0.01),np.log10(rfin+0.01),n)-0.01 #To avoid log(0)
            elif kind=='lin': self.R=np.linspace(rini,rfin,n)
        else:
            self.R=np.asarray(R)

        #set denorm
        if denorm==True: self._set_denorm(self.Mmax)
        else:
            self.Mc=1
            self.sdc=1
            self.dc=1
            self.pc=1

        # mass_arr e pot_arr per mass e pot inversion
        #Nota che cosi mass_arr e pot_arr sono denormalizzati nel caso
        self.mass_arr = self._evaluatemass(self.R*self.rc) #moltiplicato per rc, perchè vogliamo che evaluate radius isa funzione del raggio normalizzato
        self.pot_arr = self._evaluatepot(self.R*self.rc)



    def _set_denorm(self,Mmax):
        self.Mc=Mmax
        self.sdc=3*self.Mc/(2*np.pi*self.rc*self.rc)
        self.dc=4*self.Mc/(np.pi*np.pi*self.rc*self.rc*self.rc)
        self.pc=4*self.G*self.Mc/(np.pi*self.rc)

    def _evaluatedens(self,R):

        y=R/self.rc

        dd= (1 + ( y*y ) )

        return self.dc*(dd)**(-3)

    def _evaluatesdens(self,R):

        y = R / self.rc

        dd = (1 + (y * y))

        return self.sdc*dd**(-2.5)

    def _evaluatemass(self,R):

        y=R/self.rc
        y2=y*y

        cost=2/np.pi
        a=( y*(y2-1) ) / ( (1+y2)*(1+y2) )
        b=np.arctan(y)


        return self.Mc*cost*(a+b)

    def _evaluatepot(self,R):

        y=R/self.rc
        y2 = y * y
        den= 1/ ( (1+y2)*(1+y2) )

        a=( y*(y2-1) ) * den
        b=np.arctan(y)
        c=den

        res=(1/(2*y))*(a+b)+den

        return self.pc*res

    def _evaluateradius(self,x,x_type='mass'):
        """

        :param x:  normalised radius
        :param x_type:
        :return:
        """
        if x_type=='mass': ret_func=interp1d(self.mass_arr,self.R, kind=linear)
        if x_type=='pot': ret_func=interp1d(self.pot_arr,self.R, kind=linear) #we use this beacuse Univariate spline can have problem if some value on x are equals

        return ret_func(x)



