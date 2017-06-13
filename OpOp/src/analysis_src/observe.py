import copy
from ..io_src.LoadSnap import load_snap
from ..particle_src.particle import  Particles
from ..analysis_src.analysis import Analysis
from roteasy import align_frame
from math import cos,sin,sqrt
import numpy as np

angle_to_rad=np.pi/180.
Ks=1/4.74047 #from km/s  to  kpc mas  yr^-1

#TODO Automatically set all the quantities (mul, mub.. ecc)
#TODO Special Particle class Observation
class Observe():

    def __init__(self,particles,psun=(8.5,0,0),vsun=(-10,5.25,7.17),vrot=223.37875188,type=2,align_mode='auto',align_vec=(), com='iter', mq=50,**kwargs):
        """

        :param particles:
        :param psun:
        :param vsun:
        :param vrot:
        :param type:
        :param align_mode: auto, cartesian or sky
        :param align_vec: if auto align to the object centre of mass it is not used, if cartesian (xg,yg,zg) if sky (l,b,dist from Sun)
        """
        #Check input
        if 'filename' in kwargs: p=load_snap(kwargs['filename'])
        elif isinstance(particles,Particles): p=particles
        else: raise IOError('Incorrect particles or filename')


        self.psun=np.array(psun)
        self.dsung=np.sqrt(np.sum(self.psun*self.psun))
        self.vsun=np.array((vsun[0],vsun[1]+vrot,vsun[2]))
        self.type=type
        idx=p.Type==type
        self.oPos=p.Pos[idx]
        self.oVel=p.Vel[idx]
        self.align_mode=align_mode

        #set self.align_vec and self.dist_obj
        self._set_align(p,align_mode,align_vec,com=com,mq=mq)

        #rot Pos
        self.sPos, self.obsPos=self._rotate_pos()
        self.sVel, self.obsVel=self._rotate_vel()

    def _set_align(self,p,align_mode,align_vec,com,mq):
        #set align_vec=(x,y,z) of the object wrt to the Sun position
        #set dist_obj=distance of the object wrt to the Sun position
        if align_mode.lower()=='auto':
            an=Analysis(particles=p, safe=False, auto_centre=False)

            if com!='iter' or len(p.Id)<=10: pcom,vcom=an.com(mq,type=self.type)
            else: pcom,vcom=an.comit(fac=0.975,limit_mass=mq,maxiter=500,type=self.type)

            #Pcom and Vcom are the position and velocity of the centre of mass of the system
            #Move to sun
            print('pcom',pcom)
            print('vcom',vcom)
            self.align_vec = np.array(  (self.psun[0]-pcom[0], pcom[1]-self.psun[1], pcom[2]-self.psun[2])  )
            self.dist_obj=np.sqrt(np.sum(self.align_vec*self.align_vec))

        elif align_mode.lower()=='cartesian':
            self.align_vec = np.array(  (self.psun[0]-align_vec[0], align_vec[1]-self.psun[1], align_vec[2]-self.psun[2])  )
            self.dist_obj = np.sqrt(np.sum(self.align_vec * self.align_vec))

        elif align_mode.lower()=='sky':
            l=align_vec[0]
            b=align_vec[1]
            d=align_vec[2]
            self.dist_obj=d
            x=d*cos(b*angle_to_rad)*cos(l*angle_to_rad)
            y=d*cos(b*angle_to_rad)*sin(l*angle_to_rad)
            z=d*sin(b*angle_to_rad)
            self.align_vec=np.array(x,y,z)

        return 0

    def _rotate_pos(self):

        pos_sun=align_frame(self.oPos, pos_vec=self.psun, ax='x', cartesian=True, reference='l', xoff=self.dsung,
                    change_reference='x') #Cordd centered on the Sun

        pos_obj_tmp=align_frame(pos_sun,pos_vec=self.align_vec,ax='z',cartesian=True,reference='r',zoff=self.dist_obj)

        pos_obj=np.zeros_like(pos_obj_tmp)
        pos_obj[:,0]=pos_obj_tmp[:,1]
        pos_obj[:,1]=-pos_obj_tmp[:,0]
        pos_obj[:,2]=pos_obj_tmp[:,2]
        #x-l, y-b, z-r, centred on the object com

        return pos_sun, pos_obj

    def _rotate_vel(self):

        vel_sun = align_frame(self.oVel, pos_vec=self.psun, ax='x', cartesian=True, reference='l', xoff=self.vsun[0],
                            yoff=self.vsun[1], zoff=self.vsun[2], change_reference='x')

        vel_obj_tmp = align_frame(vel_sun,pos_vec=self.align_vec,ax='z',cartesian=True,reference='r')

        vel_obj=np.zeros_like(vel_obj_tmp)
        vel_obj[:,0]=vel_obj_tmp[:,1]
        vel_obj[:,1]=-vel_obj_tmp[:,0]
        vel_obj[:,2]=vel_obj_tmp[:,2]

        return vel_sun, vel_obj

    def _vobs_core(self,pos_sun,vel_sun):
        """
        Return the proper motion and ther radial velocity of an object with coordinate pos_sun_arr wrt to the
        position of the sun and velocities vel_sun_arr wrt to the Sun motion
        :param pos_sun_arr: position in Sun coordinates (Xs,Ys,Zs)
        :param vel_sun_arr:  (velocity in sun.coordiantes )  (VXs,VYs,VZs)
        :return:
        """

        distance=sqrt(pos_sun[0]*pos_sun[0]+pos_sun[1]*pos_sun[1]+pos_sun[2]*pos_sun[2])
        vel_obs = align_frame(vel_sun, pos_vec=pos_sun, ax='z', cartesian=True, reference='r')


        mul = vel_obs[1] * Ks /  distance
        mub = -vel_obs[0] *Ks  /  distance

        return (mul,mub,vel_obs[2])


    def vobs(self,pos_sun_arr,vel_sun_arr):

        if isinstance(pos_sun_arr,float) or isinstance(pos_sun_arr,int):
            return np.array(self._vobs_core(pos_sun_arr,vel_sun_arr))
        elif len(pos_sun_arr)==1:
            return np.array(self._vobs_core(pos_sun_arr, vel_sun_arr))


        vobs_tmp=[]

        for i in range(len(pos_sun_arr)):

            poss=pos_sun_arr[i]
            voss=vel_sun_arr[i]
            vobs_tmp.append(self._vobs_core(poss,voss))

        return np.array(vobs_tmp)