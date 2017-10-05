from __future__ import  division, print_function

import numpy as np
from ..model_src.Model import Model
from scipy.interpolate import interp1d

#Only 1D
class grid:
    """
    Class to define a grid and calculate some value on them
    """


    def __init__(self, N=512,  type='log', galaxy_model=None, ext_pot_model=None, setvol=True, setsup=True, efact1=1,efact2=1, **kwargs):
        """
        Purpose: initialize the grid
        :param N: Number of points in the grid
        :param type: lin or log for linearly or logaritmically spaced grid
        :param galaxy_model: Model to calculate mass, density and potential on the grid
        :param key:
                -min and max: build the grid from a min to a max
                -bins: array with the bin edges
                -data: array with some radii, this is used only to set the min and the max
                      to the min and the max of the data.
        """


        #Case 1, min, max
        if 'min' in kwargs and 'max' in kwargs:
            self.min=kwargs['min']
            self.max=kwargs['max']
            self.N=N

            if type=='log': self.gedge=np.logspace(np.log10(self.min),np.log10(self.max),self.N,dtype=float)
            elif type=='lin': self.gedge=np.linspace(self.min,self.max,self.N,dtype=float)
            else: raise ValueError("%s not supported" % str(type))



        #Case 2, user supplied bins
        elif 'bins' in kwargs:
            get=kwargs['bins']
            #Make bin_edges as np array
            if isinstance(get,tuple): self.gedge=np.array(get,dtype=float)
            elif isinstance(get,list): self.gedge=np.array(get,dtype=float)
            elif isinstance(get,np.ndarray): self.gedge=get
            else: ValueError('Bins must be a tuple a list or a numpy array')

            self.gedge=np.sort(self.gedge)
            self.min=self.gedge.min()
            self.max=self.gedge.max()
            self.N=len(self.gedge-1)

        #Case 3, user supplied data (only to calc min and max)
        elif 'data' in kwargs:
            self.min=kwargs['data'].min()
            self.max=kwargs['data'].max()
            self.N=N

            if type=='log': self.gedge=np.logspace(np.log10(self.min),np.log10(self.max),self.N,dtype=float)
            elif type=='lin': self.gedge=np.linspace(self.min,self.max,N,dtype=float)
            else: raise ValueError("%s not supported" % str(type))




        self.gx=self.bin_rad(self.gedge)
        self.gdr=self.bin_dr(self.gedge, 1)


        if (isinstance(galaxy_model,Model)==False) and (galaxy_model is not None) : raise ValueError('Model unsupported')


        #Check ext_pot_model
        if ext_pot_model is None: ext_pot_model_tmp=ext_pot_model
        elif (isinstance(ext_pot_model,Model)):
                ext_pot_model_tmp=[ext_pot_model]
        elif (isinstance(ext_pot_model,list)) or (isinstance(ext_pot_model,tuple)) or (isinstance(ext_pot_model,np.ndarray)):
                badlist=[]
                i=0
                ext_pot_model_tmp=[]
                for pm in ext_pot_model:
                    if isinstance(pm,Model)==False: badlist.append(i+1)
                    ext_pot_model_tmp.append(pm)
                    i+=1
                if len(badlist)>0: raise ValueError('Models',badlist,'in ext_pot_model unsupported')
        else:
            raise ValueError('Ext Pot Model unsupported')

        ext_pot_model=ext_pot_model_tmp



        self.gmodel=galaxy_model
        if self.gmodel != None:
            self.dgrid=self.set_dens_grid(galaxy_model)
            self.mgrid=self.set_mass_grid(galaxy_model)
            self.pgrid=self.set_pot_grid(galaxy_model)
            self.tpgrid=self.set_totpot_grid(galaxy_model,ext_pot_model)

        else:
            self.dgrid=0.
            self.mgrid=0.
            self.pgrid=0.
            self.totpgrid=0.

        if setvol: self.setvol(efact1=efact1,efact2=efact2)
        if setsup: self.setsup(efact1=efact1)

    #Set
    def set_dens_grid(self,model):

        return model.dens(self.gx)

    def set_mass_grid(self,model):

        return model.mass(self.gx)

    def set_pot_grid(self, model):
        return model.pot(self.gx)

    def set_totpot_grid(self,model,ext_pot_model):

        totpot=model.pot(self.gx)

        if ext_pot_model is not None:
            for ep in ext_pot_model:

                totpot+=ep.pot(self.gx)


        return totpot

    #Eval
    def eval_rad(self, x, ivar='mass', normalize_ivar=1):
        """
        Eval dens interpolating the grid
        :param g:
        :return:
        """
        clist={'mass':self.mgrid, 'dens':self.dgrid, 'pot':self.pgrid, 'rad':self.gx}

        if (clist[ivar] is None) or (ivar not in clist): raise Exception("Grid does not set")
        else: func=interp1d(clist[ivar]/normalize_ivar,self.gx,kind='linear')
        return func(x)

    def eval_dens(self,x,ivar='rad'):
        """
        Eval dens interpolating the grid
        :param g:
        :return:
        """
        clist={'mass':self.mgrid, 'dens':self.dgrid, 'pot':self.pgrid, 'rad':self.gx}

        if (clist[ivar] is None) or (ivar not in clist): raise Exception("Grid does not set")
        else:
            func=interp1d(clist[ivar],self.dgrid,kind='linear')
        return func(x)

    def eval_mass(self,x,ivar='rad'):
        """
        Eval dens interpolating the grid
        :param g:
        :return:
        """

        clist={'mass':self.mgrid, 'dens':self.dgrid, 'pot':self.pgrid, 'rad':self.gx}

        if (clist[ivar] is None) or (ivar not in clist): raise Exception("Grid does not set")
        else:
            func=interp1d(clist[ivar],self.mgrid,kind='linear')
        return func(x)

    def eval_pot(self,x,ivar='rad'):
        """
        Eval dens interpolating the grid
        :param g:
        :return:
        """

        clist={'mass':self.mgrid, 'dens':self.dgrid, 'pot':self.pgrid, 'rad':self.gx}

        if (clist[ivar] is None) or (ivar not in clist): raise Exception("Grid does not set")
        else:
            func=interp1d(clist[ivar],self.pgrid,kind='linear')
        return func(x)

    def eval_tpot(self,x,ivar='rad'):
        """
        Eval dens interpolating the grid
        :param g:
        :return:
        """

        clist={'mass':self.mgrid, 'dens':self.dgrid, 'pot':self.pgrid, 'rad':self.gx}

        if (clist[ivar] is None) or (ivar not in clist): raise Exception("Grid does not set")
        else:
            func=interp1d(clist[ivar],self.tpgrid,kind='linear')
        return func(x)

    @staticmethod
    def bin_dr(bin_edges,pow):

            #Make bin_edges as np array
            if isinstance(bin_edges,np.ndarray): pass
            elif isinstance(bin_edges,tuple): bin_edges=np.array(bin_edges,dtype=float)
            elif isinstance(bin_edges,list): bin_edges=np.array(bin_edges,dtype=float)
            else: ValueError('Bins must be a tuple a list or a numpy array')

            return (bin_edges[1:]**pow-bin_edges[:-1]**pow)

    @staticmethod
    def bin_rad(bin_edges):

            #Make bin_edges as np array
            if isinstance(bin_edges,np.ndarray): pass
            elif isinstance(bin_edges,tuple): bin_edges=np.array(bin_edges,dtype=float)
            elif isinstance(bin_edges,list): bin_edges=np.array(bin_edges,dtype=float)
            else: ValueError('Bins must be a tuple a list or a numpy array')

            return (bin_edges[1:]+bin_edges[:-1])*0.5

    def setvol(self,efact1=1,efact2=1, fact=None):
        '''
        Set the volume inside a bin interval. The volume is for a general ellipsoide with axes: a,b,c
        We suppose that the value of the bin_edges are referred to a, so the volume is
        V=4/3pi a * a(b/a) * a(c/a)= 4/3 pi * (b/a)*(c/a)  *a **3= fact*a**3.
        IN the case of a sphere (default) b/a and c/a are equal to 1.
        :param efact1: b/a
        :param efact2: c/a
        :key fact: User defined factor to multiply bin_edges**3(i+1) - bin_edges**3(i)
        '''
        if fact is None: fact=(4./3.)*np.pi*efact1*efact2

        self.g_vol=fact*self.bin_dr(self.gedge, 3)

    def setsup(self,efact1=1, fact=None):
        '''
        Set the area inside a bin interval. The volume is for a general ellittic curve with axes: a,b
        We suppose that the value of the bin_edges are referred to a, so the Area is
        A= pi a * a(b/a) = pi * (b/a)  *a **2= fact*a**2
        IN the case of a sphere (default) b/a and c/a are equal to 1.
        :param efact1: b/a
        :key fact: User defined factor to multiply bin_edges**2(i+1) - bin_edges**2(i)
        '''

        if fact is None: fact=np.pi*efact1

        self.g_sup=fact*self.bin_dr(self.gedge, 2)

    def set_model(self,galaxy_model):

        if isinstance(galaxy_model,Model):
            self.dgrid=self.set_dens_grid(galaxy_model)
            self.mgrid=self.set_mass_grid(galaxy_model)
            self.pgrid=self.set_pot_grid(galaxy_model)
        else: raise ValueError('Galaxy model not from model_src')

    def scale_rad(self,rc=1):
        self.gedge=self.gedge*rc

        self.gx=self.bin_rad(self.gedge)
        self.gdr=self.bin_dr(self.gedge, 1)

        self.gmodel=galaxy_model
        if self.gmodel != None:
            self.dgrid=self.set_dens_grid(galaxy_model)
            self.mgrid=self.set_mass_grid(galaxy_model)
            self.pgrid=self.set_pot_grid(galaxy_model)
        else:
            self.dgrid=0.
            self.mgrid=0.
            self.pgrid=0.

        self.setvol(efact1=efact1,efact2=efact2)
        self.setsup(efact1=efact1)



