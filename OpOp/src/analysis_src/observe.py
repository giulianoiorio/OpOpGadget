from __future__ import  division, print_function
import copy
from ..io_src.LoadSnap import load_snap
from ..particle_src.particle import  Particles
from ..analysis_src.analysis import Analysis
from ..particle_src.sky_particle import Sky_Particles, Sky_Particle
from roteasy import align_frame
from math import cos,sin,sqrt
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord

angle_to_rad=np.pi/180.
rad_to_angle=180./np.pi
Ks=1/4.74047 #from km/s  to  kpc mas  yr^-1
Ksi=4.74047

#TODO Set xi, nu, ra, dec, rad (with respect to the COM)
#TODO Give in input after the observation the feature of the COM
#TODO Automatically set all the quantities (mul, mub.. ecc)
class Observe():

    def __init__(self,particles,type=None,**kwargs):
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

        if type is None:
            type_cond=('type>1',) #take only observed stars
        elif isinstance(type,int):
            type_cond=('type=%i'%type,)
        elif isinstance(type,list) or isinstance(type,tuple) or  isinstance(type,np.ndarray):
            type_cond=[]
            for ty in type:
                type_cond.append('type=%i'%ty)
        else:
            raise ValueError('Invalid type')


        self.header=p.header
        self.p=particles.extract(*type_cond,mode='or')
        #centre is  a Sky particle
        self.centre=Sky_Particle(type=9,mass=np.sum(self.p.Mass[:]),id=0)


    def observe(self,psun=(8.0,0,0),vsun=(-10,5.25,7.17),vrot=220,align_mode='auto',align_vec=(), com='iter', mq=50,**kwargs):

        self.set_sun_position(psun=psun,vsun=vsun,vrot=vrot) #Sun Position
        align_pos, align_vel=self.set_align(align_mode, align_vec, com=com, mq=mq) #align

        #rotate
        pos_sun, pos_obs, pcentre_sun, pcentre_obs = self._rotate_pos()
        vel_sun, vel_obs, vcentre_sun, vcentre_obs = self._rotate_vel(vcentre=align_vel)


        c = self._make_centre(pcentre_sun,vcentre_sun)
        c.Pos[:]=align_pos
        c.Vel[:]=align_vel
        c.Pcord = "Cen. on GC, X-toward Sun, Y-toward Sun motion"
        c.Vcord = "Cen. on GC, X-toward Sun, Y-toward Sun motion"
        c.setRadius()
        c.setVel()
        s = self._make_sky_particles(pos_sun,pos_obs,vel_sun,vel_obs)
        self.centre=c

        return s, c

    def set_align(self,align_mode,align_vec,com,mq):
        #set align_vec=(x,y,z) of the object wrt to the Sun position
        #set dist_obj=distance of the object wrt to the Sun position
        if align_mode.lower()=='auto':
            an=Analysis(particles=self.p, safe=False, auto_centre=False)

            if com!='iter' or len(self.p.Id)<=10: pcom,vcom=an.com(mq)
            else: pcom,vcom=an.comit(fac=0.975,limit_mass=mq,maxiter=500)

            align_pos=pcom
            align_vel=vcom

            #Pcom and Vcom are the position and velocity of the centre of mass of the system
            #Move to sun

            self.align_vec = np.array(  (self.psun[0]-pcom[0], pcom[1]-self.psun[1], pcom[2]-self.psun[2])  )
            self.dist_obj=np.sqrt(np.sum(self.align_vec*self.align_vec))

        elif align_mode.lower()=='cartesian':
            self.align_vec = np.array(  (self.psun[0]-align_vec[0], align_vec[1]-self.psun[1], align_vec[2]-self.psun[2])  )
            self.dist_obj = np.sqrt(np.sum(self.align_vec * self.align_vec))

            align_pos=align_vec
            align_vel=np.array([0,0,0])

        elif align_mode.lower()=='sky':
            l=align_vec[0]
            b=align_vec[1]
            d=align_vec[2]
            self.dist_obj=d
            x=d*cos(b*angle_to_rad)*cos(l*angle_to_rad)
            y=d*cos(b*angle_to_rad)*sin(l*angle_to_rad)
            z=d*sin(b*angle_to_rad)
            self.align_vec=np.array([x,y,z])

            align_pos=np.array([x,y,z])
            align_vel=np.array([0,0,0])

        return align_pos,align_vel

    def _rotate_pos(self):

        pos_sun=align_frame(self.p.Pos, pos_vec=self.psun, ax='x', cartesian=True, reference='l', xoff=self.dsung,
                    change_reference='x') #Cordd centered on the Sun

        pos_obj_tmp=align_frame(pos_sun,pos_vec=self.align_vec,ax='z',cartesian=True,reference='r')

        pos_obj=np.zeros_like(pos_obj_tmp)
        pos_obj[:,0]=pos_obj_tmp[:,1]
        pos_obj[:,1]=-pos_obj_tmp[:,0]
        pos_obj[:,2]=pos_obj_tmp[:,2]
        #x-l, y-b, z-r, centred on the object com


        align_vec_gal=[[self.psun[0]-self.align_vec[0], self.align_vec[1]-self.psun[1], self.align_vec[2]-self.psun[2]],]
        pcentre_sun=align_frame(align_vec_gal, pos_vec=self.psun, ax='x', cartesian=True, reference='l', xoff=self.dsung,
                    change_reference='x') #Cordd centered on the Sun
        pcentre_obs=align_frame(pcentre_sun,pos_vec=self.align_vec,ax='z',cartesian=True,reference='r')

        return pos_sun, pos_obj, pcentre_sun, pcentre_obs

    def _rotate_vel(self,vcentre=None):

        vel_sun = align_frame(self.p.Vel, pos_vec=self.psun, ax='x', cartesian=True, reference='l', xoff=self.vsun[0],
                            yoff=self.vsun[1], zoff=self.vsun[2], change_reference='x')

        vel_obj_tmp = align_frame(vel_sun,pos_vec=self.align_vec,ax='z',cartesian=True,reference='r')

        vel_obj=np.zeros_like(vel_obj_tmp)
        vel_obj[:,0]=vel_obj_tmp[:,1]
        vel_obj[:,1]=-vel_obj_tmp[:,0]
        vel_obj[:,2]=vel_obj_tmp[:,2]

        if vcentre is not None:
            vcentre_sun=align_frame([vcentre,], pos_vec=self.psun, ax='x', cartesian=True, reference='l', xoff=self.vsun[0],
                            yoff=self.vsun[1], zoff=self.vsun[2], change_reference='x')
            vcentre_obs = align_frame(vcentre_sun, pos_vec=self.align_vec, ax='z', cartesian=True, reference='r')
        else:
            vcentre_sun=np.array([0,0,0])
            vcentre_obs=np.array([0,0,0])

        return vel_sun, vel_obj, vcentre_sun, vcentre_obs

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
            raise ValueError()



        vobs_tmp=[]

        for i in range(len(pos_sun_arr)):

            poss=pos_sun_arr[i]
            voss=vel_sun_arr[i]
            vobs_tmp.append(self._vobs_core(poss,voss))

        return np.array(vobs_tmp)

    def set_sun_position(self,psun=(8.0,0,0),vsun=(-10,5.25,7.17),vrot=220):

        self.psun=np.array(psun)
        self.dsung=np.sqrt(np.sum(self.psun*self.psun))
        self.vsun=np.array((vsun[0],vsun[1]+vrot,vsun[2]))

    def _sun_to_galactic(self,pos_sun):

        dist=np.sqrt(np.sum(pos_sun*pos_sun,axis=1))
        b=np.arcsin(pos_sun[:,2]/dist)*rad_to_angle
        l=np.arctan2(pos_sun[:,1],pos_sun[:,0])*rad_to_angle

        return l,b

    #TODO documentation
    def _make_sky_particles(self,pos_sun,pos_obs,vel_sun,vel_obs):

        s = Sky_Particles(N=self.p.n)

        s.Id[:]=self.p.Id
        s.Type[:]=self.p.Type
        s.Mass[:]=self.p.Mass
        s.distance[:] = np.sqrt(np.sum(pos_sun * pos_sun, axis=1))

        s.Pos[:, :] = pos_obs
        s.Vel[:, :] = vel_obs

        vvlos = self.vobs(pos_sun, vel_sun)

        s.Vlos[:]=vvlos[:,2]
        s.mul[:]=vvlos[:,0]
        s.mub[:]=vvlos[:,1]
        s.Vl[:]=s.mul * s.distance * Ksi
        s.Vb[:]=s.mub * s.distance * Ksi
        s.l[:],s.b[:]=self._sun_to_galactic(pos_sun)
        gc=SkyCoord(l=s.l * u.degree, b=s.b * u.degree, frame='galactic')
        gc=gc.fk5
        s.ra[:]=np.array(gc.ra.value)
        s.dec[:]=np.array(gc.dec.value)
        s.xi[:] = np.arctan(s.Pos[:,0]/s.distance[:])*rad_to_angle*3600
        s.eta[:] = np.arctan(s.Pos[:,1]/s.distance[:])*rad_to_angle*3600
        s.header=self.header

        return s

    def _make_centre(self,pcentre_sun,vcentre_sun):

        centre = Sky_Particle(type=9,mass=np.sum(self.p.Mass), id=0)

        dsun = np.sqrt(np.sum(pcentre_sun[0]*pcentre_sun[0]))
        b = np.arcsin(pcentre_sun[:, 2] / dsun)*rad_to_angle
        l = np.arctan2(pcentre_sun[:, 1], pcentre_sun[:, 0])*rad_to_angle
        if l < 0: l = 360 + l
        centre.l=l
        centre.b=b
        gc=SkyCoord(l=l * u.degree, b=b * u.degree, frame='galactic')
        gc=gc.fk5
        centre.ra=gc.ra.value
        centre.dec=gc.dec.value
        centre.distance=dsun



        centre_vvlos = self.vobs(np.array(pcentre_sun), np.array(vcentre_sun))
        centre_mul, centre_mb, centre_Vlos = centre_vvlos[0]
        centre.mul=centre_mul
        centre.mub=centre_mb
        centre.Vlos=centre_Vlos

        centre.Vl= centre_mul * dsun * Ksi
        centre.Vb= centre_mb * dsun * Ksi

        centre.Pos=pcentre_sun[0]
        centre.Vel=vcentre_sun[0]


        return  centre


    #TODO documentation
    #TODO minimisation alg
    def locate_sun(self, dsun=8.0,zsun=0, N=1000, l_obs=None, b_obs=None, d_obs=None, mul_obs=None, mub_obs=None, vlos_obs=None, com='iter', mq=50, set_sun=True):

        if ( l_obs is None ) and (b_obs is None) and (d_obs is None) and (mul_obs is None) and (mub_obs is None) and (vlos_obs is None):
            print('Warning all parameters set to None in locate_sun')
            return None

        if l_obs is None: l_obs=np.nan
        if b_obs is None: b_obs=np.nan
        if d_obs is None: d_obs=np.nan
        if mul_obs is None: mul_obs=np.nan
        if mub_obs is None: mub_obs=np.nan
        if vlos_obs is None: vlos_obs=np.nan


        an = Analysis(particles=self.p, safe=False, auto_centre=False)

        if com != 'iter' or len(self.p.Id) <= 10:
            pcom, vcom = an.com(mq)
        else:
            pcom, vcom = an.comit(fac=0.975, limit_mass=mq, maxiter=500)

        pcom=((pcom[0],),(pcom[1],),(pcom[2],))
        vcom = ((vcom[0],), (vcom[1],), (vcom[2],))

        x=np.linspace(0,dsun,N)
        y=np.sqrt(dsun*dsun - x*x)

        llist=np.zeros(N,dtype=float)
        blist = np.zeros(N, dtype=float)
        dlist = np.zeros(N, dtype=float)
        vloslist = np.zeros(N, dtype=float)
        mllist = np.zeros(N, dtype=float)
        mblist = np.zeros(N, dtype=float)


        for i in range(N):


            psun=(x[i], y[i], zsun)

            #l and b of com
            pcom_sun=align_frame(pcom, pos_vec=psun, ax='x', cartesian=True, reference='l', xoff=dsun,
                        change_reference='x',unpacked=True)
            dcom_sun=sqrt(pcom_sun[0,0]*pcom_sun[0,0]+pcom_sun[0,1]*pcom_sun[0,1]+pcom_sun[0,2]*pcom_sun[0,2])

            l, b = self._sun_to_galactic(pcom_sun)
            if l < 0: l += 360


            llist[i]=(l[0] - l_obs ) / l_obs
            blist[i]=(b[0] - b_obs) / b_obs
            dlist[i]=(dcom_sun - d_obs) / d_obs

            #Vlos mul and mub of com
            vel_sun = align_frame(vcom, pos_vec=psun, ax='x', cartesian=True, reference='l',
                                  xoff=self.vsun[0],
                                  yoff=self.vsun[1], zoff=self.vsun[2], change_reference='x',unpacked=True)

            vel_obs = align_frame( vel_sun, pos_vec=(
                                   l, b), ax='z',cartesian=False, spherical=True,reference='r')



            mul = vel_obs[0,1] * Ks / dcom_sun
            mub = -vel_obs[0,0] * Ks / dcom_sun
            vlos= vel_obs[0,2]

            mllist[i]=(mul -mul_obs ) /mul_obs
            mblist[i]=(mub -mub_obs ) /mub_obs
            vloslist[i]=(vlos - vlos_obs ) / vlos_obs


        final_diff=np.abs(np.vstack((llist,blist,dlist,mllist,mblist,vloslist)).T)

        prun=np.nanmax(final_diff,axis=1)
        idx_max=np.argmin(prun)
        xbest=x[idx_max]
        ybest=y[idx_max]
        diff_best=final_diff[idx_max]

        if set_sun:
            self.set_sun_position(psun=(xbest,ybest,zsun))

        return x[idx_max], y[idx_max], zsun, diff_best


