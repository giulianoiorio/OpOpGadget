from ..model_src import GeneralModel
import numpy as np






class TruncatedPlummer(GeneralModel.GeneralModel):
    def __init__(self, rc, Mmax, rt=None, R=None, rini=3e-5, rfin=300, kind='log', n=512, G='kpc km2 / (M_sun s2)',
                 denorm=True, r_physic=False, use_c=False):
        """
        Truncated double power law model:
        dens=dens0 * (r/rc)^(-gamma) * (1+r/rc)^(- (beta-gamma)) * Exp[ -(r/rt)^2]
        It simpy call the class general model with the density law above evaluating it on a grid of radius normalized to rc. This grid
        can be supplied by the user directly or can be generate with the keyword rini,rfin,kind,n.

        :param rc: Scale length
        :param Mmax: Physical Value of the Mass at Rmax (the last point of the R grid). The physical unity of dens and pot and mass
               will depends on the unity of Mmax
        :param rt:  Truncation radius (in physical unit), if None it is equal to 2*rmax(no truncation)
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
        self.Mmax =Mmax

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
        super(TruncatedPlummer, self).__init__(R=R, dens=self._adens, rc=self.rc, Ms=self.Mmax, G=G, use_c=use_c,
                                               denorm=denorm)

    def _adens(self, x):

        y = self.rc / self.rt

        dens=(1 + x*x)**(-2.5)

        return dens * np.exp(-x * x * y * y)

    def __str__(self):

        h=''
        h+='\nModel: Truncated Plummer'
        h+='\nrc: %.3f'%self.rc
        h+='\nrt: %.3f (physical)  %.3f (normalised)'%(self.rt,self.rt/self.rc)
        h+='\nrini: %.3f (physical)  %.3f (normalised)'%(self.rini*self.rc,self.rini)
        h+='\nrfin: %.3f (physical)  %.3f (normalised)'%(self.rfin*self.rc,self.rfin)
        h+='\nMass: %.3e at scale radius rs: %.3f'%(self.Ms,self.rs)
        h+='\nInput Mass: %.3e'%(self.Mmax)
        h+='\nTotal Mass: %.3e at last radius: %.3f'%(self.mass(self.rfin*self.rc),self.rfin*self.rc)
        h+='\nuse_c set to %s'%str(self.use_c)
        h+='\nuse_nparray set to %s' % str(self._use_nparray)

        return h