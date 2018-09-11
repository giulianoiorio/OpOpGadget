from __future__ import  division, print_function

from .particle import Particles, Particle
from math import sqrt, cos, sin,  asin,atan2
from ..utility_src.utility import make_fits
import numpy as np
from roteasy import rotate_frame, align_frame
import astropy.units as u
from astropy.coordinates import SkyCoord

angle_to_rad=np.pi/180.
rad_to_angle=180./np.pi
Ks=4.74047 #from  kpc mas  yr^-1 to km/s

#TODO Add  class method to define Sky_Particle starting from Sun coordinate and Galactic coordinate e.g. Sky_Partilces.sun()
#TODO Documentation!
class Sky_Particle(Particle):

    def __init__(self,l=0,b=0,distance=1,mul=0,mub=0,vlos=0,id=0,type=2,mass=0,Pos=None,Vel=None,xi=None,eta=None,centre_loc=None):
        """

        :param l:
        :param b:
        :param distance:
        :param mul:
        :param mub:
        :param vlos:
        :param id:
        :param type:
        :param mass:
        :param centre_loc: if not None, align the los for Pos e Vel  to this object. Centre loc should be in coordintate l,b
        """

        self.distance=distance #from the Sun
        self.l=l
        self.b=b
        self.mul=mul
        self.mub=mub
        self.Vlos=vlos
        gc=SkyCoord(l=l * u.degree, b=b * u.degree, distance=distance * u.kpc , pm_l_cosb=mul * u.mas/u.yr, pm_b=mub * u.mas/u.yr, frame='galactic')
        gc=gc.fk5
        self.ra=gc.ra.value
        self.dec=gc.dec.value
        self.mura=gc.pm_ra_cosdec.value
        self.mudec=gc.pm_dec.value
        self.Vl= mul * Ks * distance
        self.Vb= mub * Ks *distance


        #self.vlos=self.Vlos
        #self.vl=self.Vl
        #self.vb=self.Vb

        #set Pos
        if Pos is None:
            pos=self._calculate_pos_sky(centre_loc=centre_loc)
            self.Pos=pos
        else:
            self.Pos=Pos

        #set Vel
        if Vel is None:
            vel=self._calculate_vel_sky(centre_loc=centre_loc)
            self.Vel=vel
        else:
            self.Vel=Vel

        super(Sky_Particle,self).__init__(id=id,type=type,pos=self.Pos,vel=self.Vel,mass=mass)
        if xi is not None:
            self.xi=np.arctan(self.Pos[0]/self.distance)*rad_to_angle*3600 #arcsec
        else:
            self.xi=xi
            
        if eta is not None:
            self.eta=np.arctan(self.Pos[1]/self.distance)*rad_to_angle*3600 #arcsec
        else:
            self.eta=eta

        if centre_loc is None:
            self.Pcord = "Cen. on Sun, X-toward Gal centre, Y-toward Sun motion"
            self.Vcord = "Cen. on Sun, X-toward Gal centre, Y-toward Sun motion and, Vel relative to Sun motion"
        else:
            self.Pcord= "Cen. on Sun, Z-toward obejct in l=%3.1f deg b=%3.1f deg, X-toward l, Y-toward b"%(centre_loc[0],centre_loc[1])
            self.Vcord= "Cen. on Sun, Z-toward obejct in l=%3.1f deg b=%3.1f deg, X-toward l, Y-toward b, Vel wrt Sun motion"%(centre_loc[0],centre_loc[1])


    def _calculate_pos_sky(self,centre_loc=None):
        """ Calculate the vector pos
        If centre_loc is None, the pos is calculated with respect to the Sun, otherwise it is
        calculated with respect to the point on the sky with l, b and dist set in centre_loc.

        :param centre_loc: None or centre of reference on the sky with (l,b,dist)
        :return: the pos vector
        """


        psun=self._from_sky_to_sun(self.l,self.b,self.distance)

        if centre_loc is None: #Set to the Sun
            return psun

        elif len(centre_loc)==2:


            xs,ys,zs=psun


            pos_obj= align_frame(cord=((xs,),(ys,),(zs,)), pos_vec=(centre_loc[0],centre_loc[1]), ax='z', unpacked=True, cartesian=False,spherical=True, reference='r')[0]

            xn=pos_obj[1]
            yn=-pos_obj[0]
            zn=pos_obj[2]

            return (xn,yn,zn)

        else:
            raise ValueError('Invalide centre_loc format is should be (l,b)')

    def _calculate_vel_sky(self,centre_loc=None):

        vel = ((self.Vl,), (self.Vb,), (self.Vlos,))

        if centre_loc is None: #Calculate with respect to the sun
            frot=-(90+self.l) #to align x with x #around y
            srot=-90 #around x
            trot=self.b #around y

            veln = rotate_frame(cord=vel,angles=(frot,srot,trot),axes='yxy',unpacked=True)

        elif len(centre_loc)==2:

            #Calculate with respect to a reference point on the sky

            lcent, bcent = centre_loc

            ld=lcent - self.l
            bd=bcent - self.b

            veln = rotate_frame(cord=vel,angles=(ld,-bd),axes='yx',unpacked=True)

        else:
            raise ValueError('Invalide centre_loc format is should be (l,b)')

        return veln[0]

    def _from_sky_to_sun(self,l,b,distance):
        #l, b in degree
        l*=angle_to_rad
        b*=angle_to_rad
        R=distance*cos(b)
        z=distance*sin(b)
        x=R*cos(l)
        y=R*sin(l)

        return (x,y,z)


    def _rad_sun(self,pos):

        return sqrt(pos[0]*pos[0]+pos[1]*pos[1]+pos[2]*pos[2])

    def _sky_galactic(self,pos):

        x,y,z=pos
        dist=self._rad_sun(pos)
        b=asin(z/dist)*rad_to_angle
        l=atan2(y,x)*rad_to_angle

        return l,b


    def print_sky(self):
        print("Sky coord".center(50, "*"))
        st = "Coo system: Galactic"
        print(st.center(50, " "))
        print('l: ', self.l, 'b: ', self.b, 'distance: ', self.distance)
        print('Vlos:', self.Vlos, 'mul (Vl): %.5f (%.3f)'%(self.mul,self.Vl), 'mub (Vb): %.5f (%.3f)'%(self.mub,self.Vb))

    def __str__(self):
        mess = ""

        line = '\n' + "General".center(35, "*") + '\n' + 'Id: ' + str(self.Id) + '|' + 'Type: ' + str(
            self.Type) + '| ' + 'Mass:' + str(self.Mass) + '\n'
        mess += line

        st = "Coo system: Galactic"
        line = "Sky".center(50, '*') + '\n' + st.center(50, " ") + '\n' + 'l: ' + str(self.l) + ' b: ' + str(
            self.b) + '\n'
        line += 'Vlos: %.3f'%(self.Vlos) +  ' mul(Vl): %.5f(%.3f)'%(self.mul,self.Vl) + ' mub(Vb): %.5f(%.3f)'%(self.mub,self.Vb)  + ' mura: %.5f'%(self.mura) +  ' mudec: %.5f'%(self.mudec) + '\n'
        mess += line


        st = "Coo system: " + self.Pcord
        line = "Position".center(50, '*') + '\n' + st.center(50, " ") + '\n' + 'X: ' + str(self.Pos[0]) + ' Y: ' + str(
            self.Pos[1]) + ' Z: ' + str(self.Pos[2]) + '\n'
        rad = sqrt((self.Pos[0]) ** 2 + (self.Pos[1]) ** 2 + (self.Pos[2]) ** 2)
        line += 'Radius: ' + str(self.Radius) + ' Radius_calc: ' + str(rad) + '\n'
        mess += line

        st = "Coo system: " + self.Vcord
        line = "Velocity".center(50, "*") + '\n' + st.center(50, " ") + '\n' + 'Vx: ' + str(
            self.Vel[0]) + ' Vy: ' + str(self.Vel[1]) + ' Vz: ' + str(self.Vel[2]) + '\n'
        vel = sqrt((self.Vel[0]) ** 2 + (self.Vel[1]) ** 2 + (self.Vel[2]) ** 2)
        line += 'Veltot: ' + str(self.Vel_tot) + ' Veltot_calc: ' + str(vel)
        mess += line

        return mess

#TODO Documentation!
class Sky_Particles(Particles):

    def __init__(self,p=None,h=None,N=1):

        if ((p is None) and (h is None)):
            self.n=N
            self._initialize_vars(self.n) #Create from scratch
            self._make_header()


        elif p is None:
            if (isinstance(h,Header)): self.header=h.header
            else: self.header=h

            self.n=np.sum(self.header['Nall'])
            self._initialize_vars(self.n)
            #id
            self._makeid()
            #make type and mass
            self._maketypemass()


        else:
            if not (isinstance(p, Sky_Particle) or ((isinstance(p, (list, tuple, np.ndarray))) and (isinstance(p[0], Sky_Particle)))): raise ValueError('Incorrect particle format')
            if isinstance(p, Sky_Particle): p=[p]
            self.n = (len(p))
            self._initialize_vars(self.n) #Create from scratch
            self._fill_from_particle(p)
            self._make_header()

        self.par_dic = {'id': self.Id, 'type':self.Type, 'mass':self.Mass, 'l':self.l, 'b':self.b, 'vlos':self.Vlos, 'vl':self.Vl, 'vb':self.Vb, 'mul':self.mul, 'mub':self.mub, 'mura':self.mura, 'mudec':self.mudec, 'distance':self.distance , 'rad':self.Radius, 'velt':self.Vel_tot ,'x': self.Pos[:,0], 'y': self.Pos[:,1], 'z': self.Pos[:,2], 'vx': self.Vel[:,0], 'vy': self.Vel[:,1], 'vz': self.Vel[:,2], 'pos': self.Pos, 'vel':self.Vel}

    def _initialize_vars(self,n):

        super(Sky_Particles,self)._initialize_vars(n)
        self.xi= np.zeros(n,dtype=float)
        self.eta= np.zeros(n,dtype=float)
        self.l  = np.zeros(n, dtype=float)
        self.b  = np.zeros(n, dtype=float)
        self.ra = np.zeros(n, dtype=float)
        self.dec = np.zeros(n, dtype=float)
        self.mura = np.zeros(n, dtype=float)
        self.mudec = np.zeros(n, dtype=float)
        self.Vlos = np.zeros(n, dtype=float)
        self.Vb  = np.zeros(n, dtype=float)
        self.Vl  = np.zeros(n, dtype=float)
        self.mul  = np.zeros(n, dtype=float)
        self.mub  = np.zeros(n, dtype=float)
        self.distance  = np.zeros(n, dtype=float)

    def _fill_from_particle(self,p):
        for i in range(self.n):
                #Sky
                self.xi[i] = p[i].xi
                self.eta[i] = p[i].eta
                self.l[i]=p[i].l
                self.b[i] = p[i].b
                self.ra[i] = p[i].ra
                self.dec[i] = p[i].dec
                self.mura[i] = p[i].mura
                self.mudec[i] = p[i].mudec
                self.Vlos[i] = p[i].Vlos
                self.Vl[i] = p[i].Vl
                self.Vb[i] = p[i].Vb
                self.mul[i] = p[i].mul
                self.mub[i] = p[i].mub
                self.distance[i] = p[i].distance

                #Pos
                self.Pos[i] = p[i].Pos  # (Cartesian) Coordinates of the particles (list of tuple) (C-type float)
                self.Vel[i] = p[i].Vel  # (Cartesian) componente of the particles velocities (list of tuple) (C-type float)
                self.Id[i] = p[i].Id  # Particle identification number (C-type int)
                self.Mass[i] = p[i].Mass  # Particle mass (C-type float)
                self.Type[i] = p[i].Type  # Particle Type (C-type int): 0-Gas, 1-Halo, 2-Disk, 3-Bulge, 4-Stars, 5-Bndry

                # The following blocks exist in the snapshot only if enabled in the makefile
                self.Pot[i] = p[i].Pot  # Gravitational potential of the particles (C-type float)
                self.Acce[i] = p[i].Acce  # Acceleration of particles (C-type float)
                self.Tstp[i] = p[i].Tstp  # Time at the start of simulation  (C-type float)

                # The following blocks exist only per SPH particles (Tyep-0)
                self.U[i] = p[i].U  # SPH Particles internal energy. (C-type float)
                self.Rho[i] = p[i].Rho  # SPH Particle density (C-type float)
                # The following blocks exist in the snapshot only if enabled in the makefile for SPH particles
                self.Hsml[i] = p[i].Hsml  # Sph smoothing length
                self.Endt[i] = p[i].Endt  # Rate of change of entropy (if enabled in the makefile)

                # The following variables are not imported from file and are used only in internal method of the class
                self.Radius[i] = p[i].Radius  # Distance of the particle from the axes-origin.
                self.Vel_tot[i] = p[i].Vel_tot  # Particle total velocity
                self.Energy[i] = p[i].Energy  # Particle total energy (in the sense of epsilon in BT)
                self.Pcord[i] = p[i].Pcord
                self.Vcord[i] = p[i].Vcord

    def order(self,key='Id'):

        allowed_order_keys=('Id','Mass','Type','Pot','Acce','U','Rho','Radius','Vel_tot','Energy','Vlos','Vl','Vb')
        if key not in allowed_order_keys: raise ValueError('key: %s. Not supported order key',key)

        if key=='Id':           sort_idx=np.argsort(self.Id)
        elif key=='Radius':     sort_idx= np.argsort(self.Radius)
        elif key=='Vel_tot':    sort_idx= np.argsort(self.Vel_tot)
        elif key=='Energy':     sort_idx= np.argsort(self.Energy)
        elif key=='Mass':       sort_idx= np.argsort(self.Mass)
        elif key=='Type':       sort_idx= np.argsort(self.Type)
        elif key=='Pot':        sort_idx= np.argsort(self.Pot)
        elif key=='Acce':       sort_idx= np.argsort(self.Acce)
        elif key=='U':          sort_idx= np.argsort(self.U)
        elif key=='Rho':        sort_idx= np.argsort(self.Rho)
        elif key=='Vlos':       sort_idx= np.argsort(self.Vlos)
        elif key=='Vl':       sort_idx= np.argsort(self.Vl)
        elif key=='Vb':       sort_idx= np.argsort(self.Vb)

        self.xi = self.xi[sort_idx]
        self.eta = self.eta[sort_idx]
        self.l = self.l[sort_idx]
        self.b = self.b[sort_idx]
        self.ra = self.ra[sort_idx]
        self.dec = self.dec[sort_idx]
        self.distance = self.distance[sort_idx]
        self.Vlos = self.Vlos[sort_idx]
        self.Vl = self.Vl[sort_idx]
        self.Vb = self.Vb[sort_idx]
        self.mul = self.mul[sort_idx]
        self.mub = self.mub[sort_idx]
        self.mura = self.mura[sort_idx]
        self.mudec = self.mudec[sort_idx]

        self.Pos=self.Pos[sort_idx]
        self.Vel=self.Vel[sort_idx]
        self.Id=self.Id[sort_idx]
        self.Mass=self.Mass[sort_idx]
        self.Type=self.Type[sort_idx]
        self.Pot= self.Pot[sort_idx]
        self.Acce=self.Acce[sort_idx]
        self.Tstp=self.Tstp[sort_idx]
        self.U=self.U[sort_idx]
        self.Rho=self.Rho[sort_idx]
        self.Hsml=self.Hsml[sort_idx]
        self.Endt=self.Endt[sort_idx]
        self.Radius=self.Radius[sort_idx]
        self.Vel_tot=self.Vel_tot[sort_idx]
        self.Pcord=self.Pcord[sort_idx]
        self.Vcord=self.Vcord[sort_idx]
        self.order_var=key

        return sort_idx

    def _getitem(self,id):
        """Get a single element

        :param id: int
        :return:  a Sky_Particle element with the properties of the ith particle
        """


        part=Sky_Particle(l=self.l[id],b=self.b[id],distance=self.distance[id],mul=self.mul[id],mub=self.mub[id],vlos=self.Vlos[id],type=self.Type[id],mass=self.Mass[id])

        return part

    def _getslice(self, idslice):
        """Return a slice of Sky_Particles

        :param idslice:  slice object
        :return: a list of Sky_Particle objects sliced from the original one
        """
        if idslice.start is None: start = 0
        else: start = idslice.start

        if idslice.stop is None: stop = self.n
        else: stop = idslice.stop

        if idslice.step is None: step = 1
        else: step = idslice.step

        plist = []

        for i in range(start, stop, step):
            plist.append(Sky_Particle(l=self.l[i], b=self.b[i], distance=self.distance[i], mul=self.mul[i],
                                  mub=self.mub[i], vlos=self.Vlos[i], type=self.Type[i], Pos=self.Pos[i], Vel=self.Vel[i], xi=self.xi[i], eta=self.eta[i], mass=self.Mass[i]))

        return np.array(plist)


    def _getboolean(self,bollist):
        """Return a list of Sky_Particle objects where bollist is True

        :param bollist: a boolean array with the dimension of Sky_Particles
        :return:  a list of Sky_Particle objects
        """

        l=self.l[bollist]
        b=self.b[bollist]
        distance=self.distance[bollist]
        mul=self.mul[bollist]
        mub=self.mub[bollist]
        vlos=self.Vlos[bollist]
        mass=self.Mass[bollist]
        type=self.Type[bollist]
        id=self.Id[bollist]
        pos=self.Pos[bollist]
        vel=self.Vel[bollist]
        xi=self.xi[bollist]
        eta=self.eta[bollist]

        plist=[]

        for i in range(len(id)):
            plist.append(Sky_Particle(l=l[i],b=b[i],distance=distance[i],mul=mul[i],mub=mub
                                      [i],vlos=vlos[i],id=id[i],type=type[i],mass=mass[i], Pos=pos[i], Vel=vel[i], xi=xi[i], eta=eta[i]) )

        return plist

    def extract_sample(self, Nsample=None,  rad_max=None, position_err=None, velocity_err=None, err_distibution='uniform', rad_deg=False, save_txt=None, save_fits=None):
        """

        :param Nsample:
        :param rad_max:
        :param position_err:
        :param velocity_err:
        :param err_distibution:
        :param save_txt:
        :param save_fits:
        :return:
        """

        Noriginal=len(self.l)

        
        
        if rad_max is None:
            Nrad=Noriginal
            id_rad=self.Id[:]
            idx_rad=np.arange(0,Nrad)
            
        
        else:
            idx_rad=rad<=rad_max
            Nrad=np.sum(idx_rad)
            id_rad=self.Id[idx_rad]
        
        if rad_deg:
            rad = np.sqrt(self.xi[idx_rad]*self.xi[idx_rad]+self.eta[idx_rad]*self.eta[idx_rad])/3600.
        else:
            rad = np.sqrt(self.Pos[idx_rad, 0] * self.Pos[idx_rad, 0] + self.Pos[idx_rad, 1] * self.Pos[idx_rad, 1])
                    

        ret_arr=np.zeros(shape=(Nrad,13))
        ret_arr[:,0]=self.l[idx_rad]
        ret_arr[:,1]=self.b[idx_rad]
        ret_arr[:,2]=self.ra[idx_rad]
        ret_arr[:,3]=self.dec[idx_rad]
        ret_arr[:,4]=self.xi[idx_rad]
        ret_arr[:,5]=self.eta[idx_rad]
        ret_arr[:,6]=rad
        ret_arr[:,7]=self.mul[idx_rad]
        ret_arr[:,8]=self.mub[idx_rad]
        ret_arr[:,9]=self.Vlos[idx_rad]
        ret_arr[:,11]=self.distance[idx_rad]
        ret_arr[:,12]=id_rad
        
        if Nsample is None: 
            Nsample=Nrad
            idx_extract=np.random.choice(Nrad, Nsample, replace=False)  
            Nfinal=len(idx_extract)
            ret_arr=ret_arr[idx_extract,:]

        if velocity_err is None:
            pass
        elif isinstance(velocity_err,int) or isinstance(velocity_err,float):
            ret_arr[:,9]=np.random.normal(ret_arr[:,9],velocity_err)
            ret_arr[:,10]=velocity_err
        elif len(velocity_err,2):
            if err_distibution[0].lower()=='u':
                ret_arr[:,10]=np.random.uniform(velocity_err[0],velocity_err[1],size=Nfinal)
            elif err_distibution[0].lower()=='g':
                ret_arr[:,10]=np.random.normal(velocity_err[0],velocity_err[1],size=Nfinal)
            else:
                raise NotImplementedError('err distrbution not implemented')

            ret_arr[:,9]=np.random,normal(ret_arr[:,9],ret_arr[:,10])

        else:

            raise NotImplementedError('Vel error needs to be None, a float a int, a tule ora a list')

        if save_txt is not None:

            np.savetxt(save_txt,ret_arr, fmt='%.4f',header='0-l[deg], 1-b[deg], 2-ra[deg], 3-dec[deg], 4-xi[asec], 5-eta[asec], 6-rc[kpc], 7-mul[mas/yr] 8-mub[mas/yr] 9-Vlos[km/s] 10-eVlos[km/s] 11-dist[kpc] 12-id')

        if save_fits is not None:

            dic={}
            idl=('l','b','ra','dec','xi','eta','rc','mul','mub','Vlos','eVlos','dist','id')
            tpl=('D','D','D','D','D','D','D','D','D','D','D','D','J')
            for i in range(ret_arr.shape[1]):
                dic[idl[i]]=(ret_arr[:,i],tpl[i])

            make_fits(dic, outname=save_fits)


        return ret_arr