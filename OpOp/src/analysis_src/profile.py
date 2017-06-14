import numpy as np
from astropy.constants import G as conG
from scipy.interpolate import UnivariateSpline
from scipy.linalg import eigh
import copy

from ..particle_src.particle import  Particles
from ..io_src.LoadSnap import load_snap
from .analysis import Analysis
from ..grid_src.grid import grid




class Profile:

    def __init__(self,particles=None, type=None, center=False, mq=98,  ngrid=512, xmin=None, xmax=None, iter=False, safe=False, kind='log',**kwargs):

        #Check input
        if 'filename' in kwargs: part=load_snap(kwargs['filename'])
        elif isinstance(particles,Particles): part=particles
        else: raise IOError('Incorrect particles or filename')

        #center
        if center==True:
            a=Analysis(part,safe=safe)
            a.center(mq=mq, iter=iter,single=True)
            p=a.p
        else:
            if safe: p=copy.deepcopy(part)
            else: p=part
            p.setrad()
            p.setvelt()


        #extract arrays
        if type is None:
            self.pos=p.Pos[:]
            self.vel=p.Vel[:]
            self.rad=p.Radius[:]
            self.radcyl=None#np.sqrt(p.Pos[:,ax1]**2+p.Pos[:,ax2]**2)
            self.velpro=None#p.Vel[:,ax3]
            self.vel_tot=p.Vel_tot[:]
            self.pmass=p.Mass[:]

        elif isinstance(type,int):
            if p.header['Nall'][type]!=0: idx_type=p.Type==type
            else:   raise ValueError('Type %i not present in particles'%type)

            self.pos=p.Pos[idx_type]
            self.vel=p.Vel[idx_type]
            self.rad=p.Radius[idx_type]
            self.vel_tot=p.Vel_tot[idx_type]
            self.pmass=p.Mass[idx_type]
            self.radcyl=None#np.sqrt(p.Pos[idx_type,ax1]**2+p.Pos[idx_type,ax2]**2)
            self.velpro=None#p.Vel[idx_type,ax3]

        else: raise ValueError('type need to be None or an integer')



        #define grid
        if 'Ngrid' in kwargs: self.Ngrid=kwargs['Ngrid']
        else: self.Ngrid=ngrid
        if xmin is None: self.xmin=np.min(self.rad)
        else: self.xmin=xmin
        if xmax is None: self.xmax=np.max(self.rad)
        else: self.xmax=xmax
        self.kind=kind
        self.grid=grid(N=self.Ngrid, type=self.kind, min=self.xmin, max=self.xmax )
        self.grid.setvol()
        self.grid.setsup()


        self.massbin=None#np.histogram(self.rad,bins=self.grid.gedge,weights=self.mass)[0]
        self.masscum=None#self.massbin.cumsum()
        self.cdens=None#self.massbin/self.grid.g_vol


        self.paxdens=None
        self.paxvdisp2d=None
        self.paxvdisp3d=None
        self.massbinsup=None#np.histogram(self.radcyl,bins=self.grid.gedge,weights=self.mass)[0]
        self.masscumsup=None#self.massbinsup.cumsum()
        self.csupdens=None#self.massbinsup/self.grid.g_sup

        self.cvdisp2d=None#np.zeros_like(self.dens)
        self.cvdisp3d=None#np.zeros_like(self.dens)


        #for i in range(len(self.grid.gedge)-1):
            #cond=(self.radcyl>self.grid.gedge[i])&(self.radcyl<=self.grid.gedge[i+1])
            #print(self.grid.gedge[i],self.grid.gedge[i+1],np.max(self.radcyl[cond]))
            #self.vdisp[i]=np.std(self.velpro[cond])

    def dens(self,ret=True,func=True,s=0):
        """
        ret: If True return an array with rad e dens
        func: if True and ret True return also a spline of the dens
        s: smoothing of the spline, see scipy Univariate spline
        """

        if self.cdens is  None:

            if self.massbin is  None:
                self.massbin=np.histogram(self.rad,bins=self.grid.gedge,weights=self.pmass)[0]


            self.cdens=self.massbin/self.grid.g_vol


        if ret==True:
            retarray=np.zeros((len(self.cdens),2))
            retarray[:,0]=self.grid.gx
            retarray[:,1]=self.cdens
            if func==True:
                rfunc=UnivariateSpline(retarray[:,0],retarray[:,1],k=2,s=s)
                return retarray,rfunc
            else:
                return retarray

    def supdens(self,pax='z',ret=True,func=True,s=0):
        """
        pax: Projection axis (x,y or z)
        ret: If True return an array with rad e dens
        func: if True and ret True return also a spline of the dens
        s: smoothing of the spline, see scipy Univariate spline
        """

        if (self.csupdens is not None) and self.paxdens==pax:

            if ret==True:
                retarray=np.zeros((len(self.csupdens),2))
                retarray[:,0]=self.grid.gx
                retarray[:,1]=self.csupdens
                if func==True:
                    rfunc=UnivariateSpline(retarray[:,0],retarray[:,1],k=2,s=s)
                    return retarray,rfunc
                else:
                    return retarray

        else:
            #self.pax=pax
            if pax=='z':
                ax1=0
                ax2=1
                ax3=2
            elif pax=='y':
                ax1=0
                ax2=2
                ax3=1
            elif pax=='x':
                ax1=1
                ax2=2
                ax3=0

            self.radcyl=np.sqrt(self.pos[:,ax1]**2+self.pos[:,ax2]**2)

            if (self.massbinsup is  None) or (self.paxdens!=pax): self.massbinsup=np.histogram(self.radcyl,bins=self.grid.gedge,weights=self.pmass)[0]

            self.csupdens=self.massbinsup/self.grid.g_sup
            self.paxdens=pax

            if ret==True:
                retarray=np.zeros((len(self.csupdens),2))
                retarray[:,0]=self.grid.gx
                retarray[:,1]=self.csupdens
                if func==True:
                    rfunc=UnivariateSpline(retarray[:,0],retarray[:,1],k=2,s=s)
                    return retarray,rfunc
                else:
                    return retarray




            #self.velpro=self.vel[:,ax3]

    def mass(self,ret=True,func=True,s=0):

        if self.masscum is None:
            if self.massbin is None: self.massbin=np.histogram(self.rad,bins=self.grid.gedge,weights=self.pmass)[0]
            self.masscum=self.massbin.cumsum()

        if ret==True:
            retarray=np.zeros((len(self.masscum),2))
            retarray[:,0]=self.grid.gx
            retarray[:,1]=self.masscum
            if func==True:
                rfunc=UnivariateSpline(retarray[:,0],retarray[:,1],k=2,s=s)
                return retarray,rfunc
            else:
                return retarray

    def vdisp2d(self,pax='z',ret=True,func=True,s=0):

        if (self.cvdisp2d is None) or (self.paxvdisp2d!=pax):

            if pax=='z':
                ax1=0
                ax2=1
                ax3=2
            elif pax=='y':
                ax1=0
                ax2=2
                ax3=1
            elif pax=='x':
                ax1=1
                ax2=2
                ax3=0
            self.paxvdisp2d=pax
            self.radcyl=np.sqrt(self.pos[:,ax1]**2+self.pos[:,ax2]**2)
            self.cvdisp2d=np.zeros(len(self.grid.gedge)-1)
            for i in range(len(self.grid.gedge)-1):
                cond=(self.radcyl>self.grid.gedge[i])&(self.radcyl<=self.grid.gedge[i+1])
                self.cvdisp2d[i]=np.std(self.vel[cond,ax3])

            if ret==True:
                retarray=np.zeros((len(self.cvdisp2d),2))
                retarray[:,0]=self.grid.gx
                retarray[:,1]=self.cvdisp2d
                if func==True:
                    rfunc=UnivariateSpline(retarray[:,0],retarray[:,1],k=2,s=0)
                    return retarray,rfunc
                else:
                    return retarray

    def vdisp3d(self,ax='z',ret=True,func=True,s=0):

        if (self.cvdisp3d is None) or (self.paxvdisp3d!=ax):

            if ax=='z': ax3=2
            elif ax=='y': ax3=1
            elif ax=='x': ax3=0

            self.paxvdisp3d=ax


            self.cvdisp3d=np.zeros(len(self.grid.gedge)-1)
            for i in range(len(self.grid.gedge)-1):
                cond=(self.rad>self.grid.gedge[i])&(self.rad<=self.grid.gedge[i+1])
                self.cvdisp3d[i]=np.std(self.vel[cond,ax3])

        if ret==True:
            retarray=np.zeros((len(self.cvdisp3d),2))
            retarray[:,0]=self.grid.gx
            retarray[:,1]=self.cvdisp3d
            if func==True:
                rfunc=UnivariateSpline(retarray[:,0],retarray[:,1],k=2,s=0)
                return retarray,rfunc
            else:
                return retarray


