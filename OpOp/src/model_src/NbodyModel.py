from __future__ import  division, print_function #for compatibility with Python2
#internal import
from ..particle_src.particle import Header,Particles
from ..model_src.Model import Model
from ..grid_src.grid import grid
from ..df_src.spherical import df_isotropic
from ..analysis_src import analysis
#external import
import numpy as np
import ctypes as ct
from numpy.ctypeslib import ndpointer
import math as mt
import os
import  time

#@TODO: parallelise the generation process
#@TODO: generalise for flattened and/or non isotropic systems
class NbodyModel():
    """Class to generate a spherical and isotropic N-body realisation from a list of  dynamic models (instances of the class Model)"""
    def __init__(self,components,dff=df_isotropic, Ngrid=512, xmin=2E-3, xmax=200, kind='log', r_physic=False):
        """
        Class to Generate a Nbody models from dynamic models from the class Model.
        :param components: Array of the dynamic components, each element need to be a dictionary  {'type': , 'model': , 'npart': }:
                            -type: integer in the range (0,5),  type of particles to be used for the component, the numbers follows the Gadget-2 convention (0-Gas, 1-DM Halo, 2-Disk, 3-Bulge, 4-Star, 5-Bndry);
                            -model: instance of the class Model;
                            -npart: Number of particles to generate.
        :param dff: density distribution calculator. It is a function that takes in input the density and total potential grid and gives in output the
                    used total potential grid, the distribution function grid, and the distribution function interpolated from the df grid.
        !Grid params
        These  parameters set the properties of grids where the density, mass and total potential are computer for the each dynamic component.
        :param Ngrid: Number of sampling points, the standard is 512.
        :param xmin: Minimum value of the grid in units of r/rc if r_physic is False, otherwise is in units of kpc.
        :param xmax: Maximum value of the grid in units of r/rc if r_physic is False, otherwise is in units of kpc.
        :param r_physic: If False xmin and xmax in input ar normalised on rc (see Model Class).
        :param kind: linear or log spaced grid.

        !Method:
        The only method that users are allows to call is: generate

        !Important attributes:
        :var p: Instance of the Particle class, it  contains only zeros in the case the method generate has still not been called.
        :var components: Array with one dictionary for each dynamic models used. The keys in the dictionaries are (the values signed with
        [input] are the value given in input to the class):
            :key type [input]: integer from 0 to 5, it indicates the type  assigned to the the particles generated for the given component, the numbers follows the Gadget convention (0-Gas, 1-DM Halo, 2-Disk, 3-Bulge, 4-Star, 5-Bndry)
            :key model [input]: The dynamic model, an instance of the class Model.
            :key npart [input]: Number of particles.
            :key id: it is  the id number assigned to the considered dynamic component.
            :key grid: grid used to sample  the potential and  the density, it is an instance of the class grid.
            :key pindex: it stores the first and last index  of the particles of the considered component.
        """
        #Check components
        #Single component
        if isinstance(components,dict):
            self._check_component(components)
            self.components=(components,)
        #Multiple components
        elif hasattr(components, '__iter__'):
            i=0
            for c in components:
                self._check_component(c,i=i)
                i+=1
            components=sorted(components, key= lambda k: k['type']  )
            self.components=tuple(components)

        #Assign attributes
        self.Ngrid=Ngrid
        self.xmin=xmin
        self.xmax=xmax
        self.kind=kind
        self.r_physic=r_physic
        self.df=dff
        self._set_header()
        self._set_particles()
        for c in self.components: self._set_grid(c)

    def generate(self, use_c=True, po=None,vo=None,mq=70,set_vel=True):
        """
        Generate cartesian position and velocities for the particles of the dynamic components and assign them to the attribute p.
        The position are randomly picked from the cumulative mass distribution, while
        the velocity component are taken sampling the distribution function. If po and vo are None,
        the generated particles are centred in (0,0,0) with noi systemic motion.
        :param use_c: If True, use the fast c-implementation.
        :param po:  If not None, Move the COM of the generated particles to this cartesian position (3D tuple).
        :param vo:  If not None, change  the velocity of the COM of the generated particles to this cartesian vector (3D tuple).
        :param mq:  Fraction of mass used to estiamte the COM of the generated particles (e.g. mq=70, only stars within the radius containing 70% of the mass are used).
        :param set_vel:  If True generate velocities, if False generate only positions.
        :return: updated p instance of the Particles Class.
        """
        print('***Generate ICS: Start***')
        for c in self.components:
            print('-Component id:%i type:%i Npart=%i'%(c['id'],c['type'],c['npart']),flush=True)
            t1=time.time()
            print('     Generate Positions:',end=' ',flush=True)
            self._set_position(c)
            print('     Done',flush=True)
            if set_vel:
                print('     Generate Velocities:',end=' ',flush=True)
                self._set_vel(c,use_c=use_c)
                print('     Done',flush=True)
                print('     Done in %.3f'%(time.time()-t1),flush=True)

        if (po is not None)    or (vo is not None):
            if po is None: po=(0,0,0)
            if vo is None: vo=(0,0,0)
            print('Calculate COM and move it to (%.1f,%.1f,%.1f) with V (%.1f,%.1f,%.1f):'%(po[0],po[1],po[2],vo[0],vo[1],vo[2]),end=' ',flush=True)
            #Find the CM and set it at the position po with velocity vo
            analysis.Analysis(self.p,safe=False, mq=mq, auto_centre=True,po=po,vo=vo) #Autrocentre and if po or vo are not None change COM and COM velocity.
            print('     Done',flush=True)

        return self.p

    def _set_header(self):
        """
        Set the header for the class particles
        :return: Just a check value of True.
        """

        self.h=Header()
        self._memindex=[] #It takes into account the index in which we have to store the particles of the dynamic component,
        #If the Types of the particles are different for each component this is equal to the information obtained by Nall.

        i=0
        for c in self.components:
            mindex_ini=self.h.header['Ntot']
            self.h.header['Npart'][0][int(c['type'])]+=int(c['npart'])
            self.h.header['Massarr'][0][int(c['type'])]=c['model'].Mmax/c['npart']
            self.h.header['Ntot']+=int(c['npart'])
            mindex_ini_fin=self.h.header['Ntot']
            c['pindex']=[mindex_ini,mindex_ini_fin]
            c['id']=i
            i+=1
        self.h.header['Nall']=self.h.header['Npart'][0]

        return True

    def _set_particles(self):
        """
        Create an instance of the particle object
        :return: Just a check value of True.
        """

        self.p=Particles(h=self.h)

        return True

    def _check_component(self,component,i=0):
        """
        Check if the dictionary of the given component is in the right format.
        :param component:  input component.
        :param i:  Component index.
        :return: Just a check value of True.
        """
        fail=[]
        if 'type' not in component: fail.append('type')
        if 'model' not in component: fail.append('model')
        elif isinstance(component['model'],Model)==False: fail.append('model')
        if 'npart' not in component: fail.append('npart')
        if len(fail)>0: raise ValueError('Missed keyword in component', i ,':',fail)

        return True

    def _set_grid(self,component):
        """
        Set grid for a given component.
        :param component: set the grid for this component
        :return: Just a check value of 0.
        """

        rc=component['model'].rc

        ext_list=[]
        for c in self.components:
            if c['id']!=component['id']: ext_list.append(c['model'])

        if len(ext_list)==0: ext_list=None
        if self.r_physic: c_grid=grid(N=self.Ngrid, galaxy_model=component['model'], ext_pot_model=ext_list, type=self.kind, min=self.xmin, max=self.xmax )
        else: c_grid=grid(N=self.Ngrid, galaxy_model=component['model'], ext_pot_model=ext_list, type=self.kind, min=rc*self.xmin, max=rc*self.xmax )

        component['grid']=c_grid

        return 0

    def _set_position(self,component):
        """
        Set positions from the cumulative mass distribution of the given component.
        :param component: generate positions  for this component.
        :return: Just a check value of True.
        """

        g=component['grid']
        maxmass=np.max(g.mgrid)
        minmass=np.min(g.mgrid)
        Nrand=component['npart']
        mass_per_part=component['model'].Mmax/Nrand

        #Random generation in spherical coordinates
        u=np.random.uniform(minmass,maxmass,size=int(Nrand)) #for radius
        radius=g.eval_rad(u)
        theta=np.random.uniform(-1,1,size=int(Nrand)) #for theta angle (for a spherical simmetric distribution we have to sample uniformly in cos(theta))
        phi=np.random.uniform(0,2*np.pi,size=int(Nrand)) #uniform sample  for the phi angle from 0 to 360 degree

        #Transform in cartesian coordinates
        x=radius*(np.sqrt(1-theta*theta)*np.cos(phi)) #r (sin(teta) cos(phi))
        y=radius*(np.sqrt(1-theta*theta)*np.sin(phi)) #r (sin(teta) sin(phi))
        z=radius*theta #r (cos(teta))

        idxmax=component['pindex'][1]
        idxmin=component['pindex'][0]

        self.p.Radius[idxmin:idxmax]=radius #radius
        self.p.Pos[idxmin:idxmax,0]=x #x-coordinate
        self.p.Pos[idxmin:idxmax,1]=y #y-coordinate
        self.p.Pos[idxmin:idxmax,2]=z #z-coordinate
        self.p.Mass[idxmin:idxmax]=mass_per_part #mass of each particle
        self.p.Pot[idxmin:idxmax]=g.eval_tpot(radius) #potential at (x,yz)

        return True

    def _set_vel(self,component,use_c=True):
        """
        Set velocity components for the particles of the given component.
        :param component: generate velocities for this component.
        :param use_c: Fast implementation in c.
        :return: Just a check value of True.
        """

        g=component['grid']
        idxmax=component['pindex'][1]
        idxmin=component['pindex'][0]


        pot_grid,df_grid,df_func=self.df(g.dgrid,g.tpgrid,use_c=use_c)


        if use_c==True:
            vx,vy,vz,v=self._v_extract_c(pot_grid,df_grid,df_func,self.p.Pot[idxmin:idxmax])
        else:
            f=np.vectorize(self._v_extract,otypes=[np.float,np.float,np.float,np.float])
            vx,vy,vz,v=f(self.p.Pot[idxmin:idxmax],df_func)



        self.p.Vel[idxmin:idxmax,0]=vx
        self.p.Vel[idxmin:idxmax,1]=vy
        self.p.Vel[idxmin:idxmax,2]=vz
        self.p.Vel_tot[idxmin:idxmax]=v
        self.p.Energy[idxmin:idxmax]=self.p.Pot[idxmin:idxmax] - 0.5*v*v

        return True

    def _v_extract(self,pot,df_func):
        """
        Generate velocity using acceptance-rejection algorithm.
        It generates a spherical symmetric  velocity distribution.
        :param pot: Total potential at the position of each particle.
        :param df_func: distribution function, function of energy.
        :return: Just a check value of True.
        """

        #Initialize variables to enter in the cycle
        v=2
        ch=0

        while (v>= 1 or ch==0):
            vx,vy,vz=np.random.uniform(-1,1,size=3)
            v=vx*vx+vy*vy+vz*vz
            e=pot*(1-v)

            if v<1:
                u=np.random.random()
                umax=df_func(e)/df_func(pot)
                if u<=umax: ch=1
                else: ch=0

        norm=np.sqrt(2*pot)

        vx=vx*norm
        vy=vy*norm
        vz=vz*norm
        v=mt.sqrt(v)*norm

        return vx,vy,vz,v

    #We have to use static method to avoid problems with the load of the C extension.
    @staticmethod
    def _v_extract_c(pot_grid,df_grid,df_func,particle_pot):
        """
        Generate velocity using acceptance-rejection algorithm.
        It generates a spherical symmetric  velocity distribution.
        :param pot_grid: array that stores the potential on the grid.
        :param df_grid: array that stores the value of the distribution function on the grid.
        :param df_func: function that returns the distribution function (function of the energy).
        :param particle_pot: potential at the position of each particle.
        :return: 3 components of the velocity (vx, vy, vz) and the module of the total velocity.
        """
        ngrid=len(pot_grid)
        N=len(particle_pot)
        pot=np.zeros(N,order='C',dtype=float)
        potgpos=np.zeros(N,order='C',dtype='i4')
        dfmax=np.zeros(N,order='C',dtype=float)
        potgrid=np.zeros(ngrid,order='C',dtype=float)
        dfgrid=np.zeros(ngrid,order='C',dtype=float)
        vx=np.zeros(N,order='C')
        vy=np.zeros(N,order='C')
        vz=np.zeros(N,order='C')
        v=np.zeros(N,order='C')


        indx=np.searchsorted(pot_grid[::-1],particle_pot)
        pot[:]=particle_pot
        potgpos[:]=ngrid-indx[:]
        dfmax[:]=df_func(particle_pot)
        potgrid[:]=pot_grid
        dfgrid[:]=df_grid

        #Load C extension
        dll_name='model_c_ext/GenerateModel.so'
        dllabspath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + dll_name
        lib = ct.CDLL(dllabspath)
        v_gen=lib.v_gen
        v_gen.restype=None
        v_gen.argtypes=[ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_int, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ct.c_int,ct.c_int,ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS")]
        v_gen(pot,potgpos,dfmax,potgrid,dfgrid,N,ngrid,vx,vy,vz,v)

        return vx,vy,vz,v
