from ..model_src import GeneralModel
import numpy as np


class Tbetamodel(GeneralModel.GeneralModel):

    def __init__(self,rc,rt,Mmax,gamma=0,beta=3,R=None,rini=3e-5,rfin=300,kind='log',n=512,G='kpc km2 / (M_sun s2)',denorm=True,use_c=False):
        """
        Truncated double power law model:
        dens=dens0 * (r/rc)^(-gamma) * (1+r/rc)^(-beta) * Exp[ -(r/rt)^2]
        It simpy call the class general model with the density law above evaluating it on a grid of radius normalized to rc. This grid
        can be supplied by the user directly or can be generate with the keyword rini,rfin,kind,n.

        :param rc: Scale length
        :param rt:  Truncation radius
        :param Mmax: Physical Value of the Mass at Rmax (the last point of the R grid). The physical unity of dens and pot and mass
               will depends on the unity of Mmax
        :param gamma: first power law exponent
        :param beta: secondo powe law exponent
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
        :param use_c: To calculate pot and mass with a C-cyle, WARNING it creates more noisy results
        :return:
        """

        if R is None:
            if kind=='log': R=np.logspace(np.log10(rini+0.01),np.log10(rfin+0.01),n)-0.01 #To avoid log(0)
            elif kind=='lin': R=np.linspace(rini,rfin,n)
        else:
            R=np.asarray(R)

        self.rt=rt
        self.rc=rc
        self.gamma=gamma
        self.beta=beta
        super(Tbetamodel,self).__init__(R,self._adens,self.rc,Mmax,G,use_c=use_c,denorm=denorm)

    def _adens(self,x):

        y=self.rc/self.rt

        dens= ( x**self.gamma ) * (  (1+x)**self.beta   )

        return (1./dens)*np.exp(-x*x*y*y)
