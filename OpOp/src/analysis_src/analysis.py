import numpy as np
from astropy.constants import G as conG
from scipy.interpolate import UnivariateSpline
from scipy.linalg import eigh
import copy

from ..particle_src.particle import  Particles
from ..particle_src.sky_particle import Sky_Particles
from ..io_src.LoadSnap import load_snap
from ..utility_src.utility import nparray_check, Continue_check, stat_error
from ..grid_src.grid import grid
from .profile import Profile


class Analysis:

    def __init__(self, particles=None, safe=True, auto_centre=False, iter=False, **kwargs):


        #Check input
        if 'filename' in kwargs: self.p=load_snap(kwargs['filename'])
        elif isinstance(particles,Particles):
            if safe==True: self.p=copy.deepcopy(particles)
            else: self.p=particles
        else: raise IOError('Incorrect particles or filename')

        self.issky = False

        #Assure that rad and total vel are stored
        self.p.setrad()
        self.p.setvelt()
        self._set_pindex()
        self.pcom=None
        self.pvcom=None

        if 'single' in kwargs: single=kwargs['single']
        else: single=True

        if auto_centre==True:
            if 'po' in kwargs:po=kwargs['po']
            else: po=(0,0,0)
            if 'vo' in kwargs:vo=kwargs['vo']
            else: vo=(0,0,0)
            if 'mq' in kwargs:
                mq=kwargs['mq']
                self.center(mq=mq,single=single,iter=iter,po=po,vo=vo)
            else:
                if iter: self.center(mq=10,iter=True,single=single,po=po,vo=vo)
                else: self.center(mq=70,single=single,iter=False,po=po,vo=vo)

        if isinstance(particles,Sky_Particles):
            self.issky=True

    def mass(self,rad,type=None):

        if type is None:
            rad_array=self.p.Radius[:]
            mas_array=self.p.Mass[:]
        else:
            type= nparray_check(type)
            rad_array=self._make_array(self.p.Radius,type)
            mas_array=self._make_array(self.p.Mass,type)

        idx=rad_array<=rad

        mass=np.sum(mas_array[idx])

        return mass

    def massup(self,rad,pax='z',type=None):


        if pax == 'obs':

            if self.issky:
                radcyl = np.sqrt(self.p.Pos[:, 0] ** 2 + self.p.Pos[:, 1] ** 2)  # set R
            else:
                raise AttributeError('obs profile is available only for object of the sky particles class')

        else:

            if pax == 'z':
                ax1 = 0
                ax2 = 1
                ax3 = 2
            elif pax == 'y':
                ax1 = 0
                ax2 = 2
                ax3 = 1
            elif pax == 'x':
                ax1 = 1
                ax2 = 2
                ax3 = 0


            radcyl = np.sqrt(self.p.Pos[:, ax1] ** 2 + self.p.Pos[:, ax2] ** 2)

        if type is None:
            rad_array=radcyl
            mas_array=self.p.Mass[:]
        else:
            type= nparray_check(type)
            rad_array=self._make_array(radcyl,type)
            mas_array=self._make_array(self.p.Mass,type)

        idx = rad_array <= rad

        mass = np.sum(mas_array[idx])

        return mass

    def vdisp(self,rad_max=None,pax='z',type=None):

        if pax == 'obs':

            if self.issky:
                radcyl = np.sqrt(self.p.Pos[:, 0] ** 2 + self.p.Pos[:, 1] ** 2)  # set R
                vel_proj = self.p.Vlos[:]  # set Vproj
            else:
                raise AttributeError('obs profile is available only for object of the sky particles class')

        else:

            if pax == 'z':
                ax1 = 0
                ax2 = 1
                ax3 = 2
            elif pax == 'y':
                ax1 = 0
                ax2 = 2
                ax3 = 1
            elif pax == 'x':
                ax1 = 1
                ax2 = 2
                ax3 = 0

            radcyl = np.sqrt(self.p.Pos[:, ax1] ** 2 + self.p.Pos[:, ax2] ** 2)  # set R
            vel_proj = self.p.Vel[:, ax3]  # set Vproj

        if type is not None:
            type= nparray_check(type)
            radcyl=self._make_array(radcyl,type)
            vel_proj=self._make_array(vel_proj,type)


        if rad_max is not None:
            idx=radcyl<=rad_max
            vd=np.std(vel_proj[idx])
        else:
            vd=np.std(vel_proj)


        return vd

    def binned_dispersion(self,bins,range=None,pax='z',Nperbin=None,type=None, bins_kind='lin', velocity_err = None, err_distibution = 'uniform', nboot=1000):
        """

        :param bins:
        :param range:
        :param Nperbins:
        :return:
        """

        #Set projection
        if pax == 'obs':

            if self.issky:
                radcyl   = np.sqrt(self.p.Pos[:, 0] ** 2 + self.p.Pos[:, 1] ** 2)  # set R
                vel_proj = self.p.Vlos[:]  # set Vproj
            else:
                raise AttributeError('obs profile is available only for object of the sky particles class')

        else:

            if pax == 'z':
                ax1 = 0
                ax2 = 1
                ax3 = 2
            elif pax == 'y':
                ax1 = 0
                ax2 = 2
                ax3 = 1
            elif pax == 'x':
                ax1 = 1
                ax2 = 2
                ax3 = 0

            radcyl   = np.sqrt(self.p.Pos[:, ax1] ** 2 + self.p.Pos[:, ax2] ** 2)  # set R
            vel_proj = self.p.Vel[:,ax3]


        #Set type
        if type is not None:
            type= nparray_check(type)
            radcyl=self._make_array(radcyl,type)
            vel_proj=self._make_array(vel_proj,type)


        #Set bins
        if isinstance(bins,int) or isinstance(bins,float):

            #Set range
            if range is None:
                rmin=np.min(radcyl)
                rmax=np.max(radcyl)
            elif isinstance(range,list) or isinstance(range,tuple) or isinstance(range,np.ndarray):
                if len(range)==2:
                    rmin,rmax=range
                else:
                    raise ValueError('range needs to be None or a list/tuple/numpy array with len 2')
            else:
                raise ValueError('range needs to be None or a list/tuple/numpy array with len 2')


            if      bins_kind =='lin': bins=np.linspace(rmin,rmax,int(bins))
            elif    bins_kind =='log': bins=np.logspace(np.log10(rmin),np.log10(rmax),int(nbins))
            else:   raise ValueError('bins kind needs to be lin or log')

        #Nperbins
        if isinstance(Nperbin, int) or isinstance(Nperbin, float):

            Nperbin=np.array([int(Nperbin),]*(len(bins)-1))

        elif len(Nperbin)!=(len(bins)-1):

            raise ValueError('Nperbins needs to be an integer or a iterable with len=bins-1')

        #calculate

        rad_mean, rad_err, vmean, vmean_err, vdisp, vdisp_err=self._vdisp_bins_extractor(radcyl,vel_proj,bins,Nperbin,velocity_err,err_distibution,nboot)


        return rad_mean, rad_err, vmean, vmean_err, vdisp, vdisp_err

    def _vdisp_bins_extractor(self,rad,vproj,bins,Nperbin=None, velocity_err = None, err_distibution = 'uniform',nboot=1000):
        """

        :param rad:
        :param vproj:
        :param bins:
        :param Nperbin:
        :param velocity_err:
        :param err_distibution:
        :param nboot:
        :return:
        """
        rad_mean=np.zeros(shape=(len(bins)-1))
        rad_err=np.zeros_like(rad_mean)
        vmean=np.zeros_like(rad_mean)
        vmean_err=np.zeros_like(rad_mean)
        vdisp=np.zeros_like(rad_mean)
        vdisp_err=np.zeros_like(rad_mean)

        for i in range(len(bins)-1):

            rmin=bins[i]
            rmax=bins[i+1]

            idx=(rad>=rmin)&(rad<=rmax)
            rad_tmp=rad[idx]
            vproj_tmp=vproj[idx]

            #extract number
            if Nperbin is not None:
                Nthisbin=int(Nperbin[i])
                #print(i,Nthisbin,len(rad_tmp))
                idx_extract = np.random.choice(len(rad_tmp), Nthisbin, replace=False)
                rad_tmp=rad_tmp[idx_extract]
                vproj_tmp=vproj_tmp[idx_extract]

            vproj_tmp_err=np.ones_like(vproj_tmp)

            #set err
            if velocity_err is None:
                pass
            elif isinstance(velocity_err, int) or isinstance(velocity_err, float):
                vproj_tmp = np.random.normal(vproj_tmp, velocity_err)
                vproj_tmp_err = velocity_err
            elif len(velocity_err, 2):
                if err_distibution[0].lower() == 'u':
                    vproj_tmp_err = np.random.uniform(velocity_err[0], velocity_err[1], size=Nfinal)
                elif err_distibution[0].lower() == 'g':
                    vproj_tmp_err = np.random.normal(velocity_err[0], velocity_err[1], size=Nfinal)
                else:
                    raise NotImplementedError('err distrbution not implemented')

                vproj_tmp = np.random, normal(ret_arr[:, 9], ret_arr[:, 10])

            else:

                raise NotImplementedError('Vel error needs to be None, a float a int, a tule ora a list')


            #calculate
            mean_tmp, meanerr_tmp, std_tmp, stderr_tmp = stat_error(vproj_tmp, vproj_tmp_err, err_type='bootstrap',
                                                            bootstrap_error=True, n=nboot,
                                                            vsys_fix=None)

            rad_mean[i]  =  np.mean(rad_tmp)
            rad_err[i]   =  np.std(rad_tmp)
            vmean[i]     =  mean_tmp
            vmean_err[i] =  meanerr_tmp
            vdisp[i]     =  std_tmp
            vdisp_err[i] =  stderr_tmp

        return rad_mean, rad_err, vmean, vmean_err, vdisp, vdisp_err

    def qmass(self,q,safe_mode=True,type=None, rad_max=None):
        """
        Calculate the radius in which the mass is a q % fraction of the total mass.
        :param q: q is the mass fraction, it can range from 0 to 100
        :param safe_mode:
        :param type:
        :return: An array with the value of the radius where mass is a q % fraction of the total mass.
                 and the total value of the mass for the type choosen.
        """

        if type is None:
            rad_array=self.p.Radius[:]
            mas_array=self.p.Mass[:]
        else:
            type= nparray_check(type)
            rad_array=self._make_array(self.p.Radius,type)
            mas_array=self._make_array(self.p.Mass,type)

        if rad_max is None:
            ret_value,_=self.qradius_ext(rad_array,mas_array,q,safe_mode=safe_mode)
        else:
            idx=rad_array<=rad_max
            ret_value,_=self.qradius_ext(rad_array[idx],mas_array[idx],q,safe_mode=safe_mode)


        return  ret_value

    def qmassup(self,q,pax='z',safe_mode=True,type=None,rad_max=None):
        """
        Calculate the radius in which the mass is a q % fraction of the total mass in 2D.
        :param q: q is the mass fraction, it can range from 0 to 100
        :param safe_mode:
        :param type:
        :return: An array with the value of the radius where mass is a q % fraction of the total mass.
                 and the total value of the mass for the type choosen.
        """

        if pax == 'obs':

            if self.issky:
                radcyl = np.sqrt(self.p.Pos[:, 0] ** 2 + self.p.Pos[:, 1] ** 2)  # set R
            else:
                raise AttributeError('obs profile is available only for object of the sky particles class')

        else:
            if pax == 'z':
                ax1 = 0
                ax2 = 1
            elif pax == 'y':
                ax1 = 0
                ax2 = 2
            elif pax == 'x':
                ax1 = 1
                ax2 = 2

            radcyl = np.sqrt(self.p.Pos[:, ax1] ** 2 + self.p.Pos[:, ax2] ** 2)

        if type is None:
            rad_array=radcyl
            mas_array=self.p.Mass[:]
        else:
            type= nparray_check(type)
            rad_array=self._make_array(radcyl,type)
            mas_array=self._make_array(self.p.Mass,type)


        if rad_max is None:
            ret_value,_=self.qradius_ext(rad_array,mas_array,q,safe_mode=safe_mode)
        else:
            idx=rad_array<=rad_max
            ret_value,_=self.qradius_ext(rad_array[idx],mas_array[idx],q,safe_mode=safe_mode)

        return  ret_value

    def com(self,mq=None,minrad=None,maxrad=None,type=None):
        """
        Calculate the position of the center of mass and its velocity
        :param mq [0-100]: Fraction of the mass where to calc the com. If different from None, the com will be calculated from 0 to rq, where rq is the mass that contain the fraction mq
                            mq will have the proprity with respect of minrad and maxrad, that will be used only in the case mq=None.
        :param minrad: Min radius to use to calculate the com. If None use the minimum radius of the particles.
        :param maxrad: Max radius to use to calculate the com. If None use maximum radius of the particles.
        :param type: Particle type to use to calculate the mean. If None use all the component together
        :return: com: array with the xyz position of the com, vcom: array with vx vy vz of the com.
        """
        if type is None:
            rad_array=self.p.Radius[:]
            vel_array=self.p.Vel[:]
            pos_array=self.p.Pos[:]
            mas_array=self.p.Mass[:]
        else:
            type= nparray_check(type)
            rad_array=self._make_array(self.p.Radius,type)
            vel_array=self._make_array(self.p.Vel,type)
            pos_array=self._make_array(self.p.Pos,type)
            mas_array=self._make_array(self.p.Mass,type)

        if mq is None:
            if minrad is None: minrad=np.min(rad_array)
            if maxrad is None: maxrad=np.max(rad_array)
        else:
            minrad=np.min(rad_array)
            maxrad,_=self.qradius_ext(rad_array,mas_array,mq)


        idx_bound=  (rad_array>=minrad) & (rad_array<=maxrad)

        ret_com= self.massw_mean(pos_array[idx_bound], mas_array[idx_bound])
        ret_vom= self.massw_mean(vel_array[idx_bound], mas_array[idx_bound])


        return ret_com, ret_vom

    def comit(self,fac=0.975,limit_mass=10,maxiter=500,type=None):
        """
        Calculate iteratively the COM following Power et al, 2003.
        :param fac:  Each step the maximum radius will decrease of a factor equal to fac.
        :param limit_mass: The iteration will stop when the mass encircled on the current maxradius is less than limit_mass*Mtot/100
        :param maxiter: The iteration will stop when the number of steps is greater than maxiter.
        :param type: use only this type of particle to calculate com.
        :return: com, vcom
        """

        iterstep=0
        #pcom=(0,0,0)
        #vcom=(0,0,0)

        limit_mass=limit_mass/100.

        if type is None:
            p=np.copy(self.p.Pos[:])
            v=np.copy(self.p.Vel[:])
            rad_array=np.copy(self.p.Radius[:])
            mas_array=np.copy(self.p.Mass[:])
            mass_t=np.sum(mas_array)
        else:
            type= nparray_check(type)
            v=np.copy(self._make_array(self.p.Vel,type))
            p=np.copy(self._make_array(self.p.Pos,type))
            rad_array=np.copy(self._make_array(self.p.Radius,type))
            mas_array=np.copy(self._make_array(self.p.Mass,type))
            mass_t=np.sum(mas_array)
        
        pcom=self.massw_mean(p, mas_array)
        vcom=self.massw_mean(v, mas_array)
        p=p - pcom
        v=v - vcom
        rad_array=np.sqrt(np.sum(p*p,axis=1))


        maxrad=np.max(rad_array)*fac
        idx_bound=  (rad_array>=0) & (rad_array<=maxrad)
        mass=np.sum(mas_array[idx_bound])
        mbody=mass/mass_t

        while (mbody>limit_mass) and (iterstep<=maxiter):



            tpos= self.massw_mean(p[idx_bound], mas_array[idx_bound])
            tvel= self.massw_mean(v[idx_bound], mas_array[idx_bound])

            p=p - tpos
            v=v - tvel
            rad_array=np.sqrt(np.sum(p*p,axis=1))
            pcom=np.array(tpos)+pcom
            vcom=np.array(tvel)+vcom

            maxrad=maxrad*fac
            idx_bound=  (rad_array>=0) & (rad_array<=maxrad)
            mass=np.sum(mas_array[idx_bound])
            mbody=mass/mass_t

            iterstep+=1




        del p
        del v
        del mas_array
        del rad_array


        return pcom,vcom

    def tdyn(self,mq=100,type=None,G='(kpc3)/(M_sun s2)'):
        """
        Calculate the dynamical time of a stystem as Tdyn=0.5*pi*Sqrt(Rh^3/(G*Mh)).
        where Rh is the radius that contains the fraction Mh=h*M of the mass.
        This is the time for a particle at distance r to reach r=0 in a homogeneus spherical system
        with density rho. We do not have an homogenues sphere, but we use a medium value rho_mh
        equals to Mh/(4/3 * pi * Rh^3).
        :param mq: Fraction of the mass to use, it can ranges from 0 to 100
        :param type: Type of particle to use, it need to be an array. If None use all
        :param G: Value of the gravitational constant G, it can be a number of a string.
                    If G=1, the physical value of the potential will be Phi/G.
                    If string it must follow the rule of the unity of the module.astropy constants.
                    E.g. to have G in unit of kpc3/Msun s2, the input string is 'kpc3 / (M_sun s2)'
                    See http://astrofrog-debug.readthedocs.org/en/latest/constants/

        :return: Dynamical tyme. The units wil depends on the units used in G, If used the
                 G in default the time will be in unit of Gyr.
        """

        if type is None:
            rad_array=self.p.Radius[:]
            mas_array=self.p.Mass[:]
        else:
            type= nparray_check(type)
            rad_array=self._make_array(self.p.Radius,type)
            mas_array=self._make_array(self.p.Mass,type)


        if isinstance(G,float) or isinstance(G,int): G=G
        else:
            GG=conG.to(G)
            G=GG.value

        rq,_=self.qradius_ext(rad_array,mas_array,mq)

        
        mass_tot_rq=np.sum(self.p.Mass[self.p.Radius<=rq])

        mass_phy=G*(mq/100)*mass_tot_rq
        tdyn=0.5*np.pi*rq*np.sqrt(rq/mass_phy)


        return tdyn

    def softening_scale(self,mq=70,rq=None,auto=True,r=None,dens=None,mass=None,kernel='Gadget',type=None):
        """
        Calculate the optimal softening scale following Dehnen, 2012 eps=cost*a(dens)*N^-0.2. The output will be in unit
        of r. If Auto==True, r and dens will be not considered.
        :param mq: Mass fraction where calcualte the softening_scale.
        :param auto: If True calculate the r-dens gride using the grid and Profile class
                wit 512 points from 0.001*rq to 10*rq.
        :param r: Array with the sampling radii.
        :param dens: Array with the density at the sampling radii. Its unity need to be the same of mass/r^3
        :param mass: Total mass of the system, the method will calculate in automatic the fraction mq/100*mass
        :param kernel: Kernel to use. Different kernel have different constant C. The implemented kernels are:
                        -spline: generic cubic spline (as in Dehnen, 2012)
                        -Gadget: to calculate the softening_scale using the spline kernel of Gadget2
        :return: the softening scale.
        """
        opt_dict={'Gadget':0.698352, 'spline':0.977693}
        

        if rq is None: rq=self.qmass(mq,type=type)


        if mass is None:
            if type is None:
                mass=np.sum(self.p.Mass[:])
            else:
                type_tmp= nparray_check(type)
                masarr=self._make_array(self.p.Mass,type_tmp)
                mass=np.sum(masarr)


        if auto==True:
            prof=Profile(self.p,Ngrid=512,xmin=0.01*rq,xmax=10*rq,kind='lin',type=type)
            arr=prof.dens()[0]
            r=arr[:,0]
            dens=arr[:,1]



        dens_spline=UnivariateSpline(r, dens, k=1,s=0,ext=1)
        der_dens=dens_spline.derivative()




        derdens=der_dens(r)
        ap=UnivariateSpline(r, r*r*dens*derdens*derdens, k=1,s=0,ext=1)
        bp=UnivariateSpline(r, r*r*dens*dens, k=1,s=0,ext=1)

        B=bp.integral(0,rq)
        A=ap.integral(0,rq)/(mass*mq/100.)
        C=(B/(A))**(1/5)

        cost=opt_dict[kernel]
        N=len(self.p.Id)**(1/5)

        return C*cost/N

    @staticmethod
    def qradius_ext(radius_array,mass_array,q,safe_mode=True):
        """

        :param radius_array:
        :param mass_array:
        :param q:
        :param safe_mode:
        :return:
        """
        if safe_mode==True:

            #check array
            radius_array=nparray_check(radius_array)
            mass_array=nparray_check(mass_array)

            #check q
            qlist=nparray_check(q)
            check_q=qlist>100
            if np.any(check_q): raise ValueError('The percent massfraction must be in the range[0-100]')

            #check mass
            if np.any(mass_array<=0):  Continue_check('Warning one or more particles have Mass<=0')


        idx_sort=np.argsort(radius_array)
        radsort=radius_array[idx_sort]
        massort=mass_array[idx_sort]
        cum_mass=massort.cumsum()
        cum_mass=100*cum_mass/cum_mass[-1]

        qmass_func=UnivariateSpline(cum_mass,radsort,k=1,s=0,ext=1)
        radq=qmass_func(q)

        return radq,qmass_func

    @staticmethod
    def massw_mean(x_array,mass_array):
        """
        Mass weighted mean on the array x_array
        :param x_array: array where to calculate the wmean
        :param mass_array: mass array
        :return: massweighted mean
        """

        ret_arr=(mass_array*x_array.T)/(np.sum(mass_array))

        return np.sum(ret_arr,axis=1)

    def center(self,mq, iter=False, single=True, type=None,  fac=0.975, maxiter=500, po=(0,0,0),vo=(0,0,0), ):
        """
        Translate the system in the center of mass, Usually is the first thing to done after the call
        of the Class. The function has two modes, if single==False, the com and vcom will be calculated
        over all the particle specified in type, and all the particles will be moved in the com and vcom.
        If single==True, a separate com and vcom will be calculated for each particle type and only the
        related particle will be moved in such com and vcom.
        The second method seems to give better results.
        :param mq: if iter=False, Mass fraction of particles used to calculate the com.
                if iter=True, limit_mass of the iteration (see comit)
        :param iter: Calculate com iteratively?
        :param single: If true enable the single method, see above.
        :param type:  IF single==False, use only this type of particle to calculate com.
        :param fac: if iter=True, Each step the maximum radius will decrease of a factor equal to fac.
        :param maxiter: if iter=True, The iteration will stop when the number of steps is greater than maxiter.
        :param po: Move the particle COM at position po
        :param vo: Set the velocity of COM at vo
        """

        po=np.array(po)
        vo=np.array(vo)

        if single==False:

            if iter: com,vcom=self.comit(fac=fac,limit_mass=mq,maxiter=maxiter,type=type)
            else:    com,vcom=self.com(mq=mq, type=type)


            self.p.Pos=self.p.Pos - com + po
            self.p.Vel=self.p.Vel - vcom + vo
            self.p.setrad()
            self.p.setvelt()
            self.pcom = com
            self.pvcom = vcom

        else:
            pcom = []
            pvcom = []
            i=0
            for i in range(6):

                co=self.p.header['Nall'][i]

                if co!=0:

                    if iter: com,vcom=self.comit(fac=fac,limit_mass=mq,maxiter=maxiter,type=[i])
                    else:    com,vcom=self.com(mq=mq, type=[i])


                    idxmin=self.pindex[i][0]
                    idxmax=self.pindex[i][1]
                    self.p.Pos[idxmin:idxmax]=self.p.Pos[idxmin:idxmax] - com + po
                    self.p.Vel[idxmin:idxmax]=self.p.Vel[idxmin:idxmax] - vcom + vo
                    self.p.setrad()
                    self.p.setvelt()
                    pcom.append(com)
                    pvcom.append(vcom)

                i+=1
            self.pcom=np.array(pcom)
            self.pvcom=np.array(pvcom)

    def inertia_tensor(self,eig=True,mq=None,minrad=None,maxrad=None,type=None):
        """
        Calculate the inertia tensor and eventually the three eigenvector.
        :param eig: If True calculate the eigenvector of the symmetric inertia tensor matrix.
        :param mq [0-100]: Fraction of the mass where to calc the angolar momentum. If different from None, the com will be calculated from 0 to rq, where rq is the mass that contain the fraction mq
                            mq will have the proprity with respect of minrad and maxrad, that will be used only in the case mq=None.
        :param minrad: Min radius to use to calculate the com. If None use the minimum radius of the particles.
        :param maxrad: Max radius to use to calculate the com. If None use maximum radius of the particles.
        :param type: Particle type to use to calculate the mean. If None use all the component together
        :return: com: array with the xyz position of the com, vcom: array with vx vy vz of the com.
        """

        if type is None:
            rad_array=self.p.Radius[:]
            pos_array=self.p.Pos[:]
            mas_array=self.p.Mass[:]
        else:
            type= nparray_check(type)
            rad_array=self._make_array(self.p.Radius,type)
            pos_array=self._make_array(self.p.Pos,type)
            mas_array=self._make_array(self.p.Mass,type)

        if mq is None:
            if minrad is None: minrad=np.min(rad_array)
            if maxrad is None: maxrad=np.max(rad_array)
        else:
            minrad=np.min(rad_array)
            maxrad,_=self.qradius_ext(rad_array,mas_array,mq)

        idx_bound=  (rad_array>=minrad) & (rad_array<=maxrad)

        return self._set_tensor(pos_array[idx_bound],mas_array[idx_bound],eig=eig)

    def vinertia_tensor(self,eig=True,mq=None,minrad=None,maxrad=None,type=None):
        """
        Calculate the kinetic tensor  and eventually the three eigenvector.
        :param eig: If True calculate the eigenvector of the symmetric inertia tensor matrix.
        :param mq [0-100]: Fraction of the mass where to calc the angolar momentum. If different from None, the com will be calculated from 0 to rq, where rq is the mass that contain the fraction mq
                            mq will have the proprity with respect of minrad and maxrad, that will be used only in the case mq=None.
        :param minrad: Min radius to use to calculate the com. If None use the minimum radius of the particles.
        :param maxrad: Max radius to use to calculate the com. If None use maximum radius of the particles.
        :param type: Particle type to use to calculate the mean. If None use all the component together
        :return: com: array with the xyz position of the com, vcom: array with vx vy vz of the com.
        """

        if type is None:
            rad_array=self.p.Radius[:]
            vel_array=self.p.Vel[:]
            mas_array=self.p.Mass[:]
        else:
            type= nparray_check(type)
            rad_array=self._make_array(self.p.Radius,type)
            vel_array=self._make_array(self.p.Vel,type)
            mas_array=self._make_array(self.p.Mass,type)

        if mq is None:
            if minrad is None: minrad=np.min(rad_array)
            if maxrad is None: maxrad=np.max(rad_array)
        else:
            minrad=np.min(rad_array)
            maxrad,_=self.qradius_ext(rad_array,mas_array,mq)

        idx_bound=  (rad_array>=minrad) & (rad_array<=maxrad)

        return self._set_tensor(vel_array[idx_bound],mas_array[idx_bound],eig=eig)

    def _set_pindex(self):
        pindex=[[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]]
        cum_count=0
        i=0
        for npar_t in self.p.header['Nall']:
            if npar_t!=0:
                pindex[i]=[cum_count,cum_count+npar_t]
                cum_count+=npar_t
            i+=1
        self.pindex=pindex

    def _make_array(self,original_array,type):
        """
        Extract from the original array only the particles belonging to the tye in the array type
        :param original_array: Array where to extract value
        :param type: Type of particle to extract, it need to be an array
        :return: A sub-sample of the original_array with only
        """

        #Define array dimension

        ndim=0

        for i in type: ndim+=self.p.header['Nall'][i]

        original_shape=original_array.shape



        if len(original_shape)==1:
                ret_arr=np.zeros(shape=(ndim))
        elif len(original_shape)==2:
                axis2=original_shape[1]
                ret_arr=np.zeros(shape=(ndim,axis2))
        else: raise ValueError('Original array dimension need to have 1 or 2 axis')



        cumsum=0
        for i in np.sort(type):
            idxmin=self.pindex[i][0]
            idxmax=self.pindex[i][1]
            npart_type=self.p.header['Nall'][i]
            ret_arr[cumsum:cumsum+npart_type]=original_array[idxmin:idxmax]
            cumsum+=npart_type

        return ret_arr

    def _set_tensor(self,value_array,weight_array,eig=True):

        mat=weight_array*value_array.T
        mat2=mat.dot(value_array)

        if eig:
            evalue, evec = eigh(mat2, eigvals_only = False)
            evalue = evalue /np.max(evalue)
            eversor = evec.T / (np.sqrt(np.sum(evec*evec, axis=1)))
            eversor = eversor.T
            #eversor=evec
            return mat2, evalue, eversor

        else:
            return mat2

    def angolar_momentum(self,mq=None,minrad=None,maxrad=None,type=None):
        """
        Calculate the three components of the angular momentum.
        :param mq [0-100]: Fraction of the mass where to calc the angolar momentum. If different from None, the com will be calculated from 0 to rq, where rq is the mass that contain the fraction mq
                            mq will have the proprity with respect of minrad and maxrad, that will be used only in the case mq=None.
        :param minrad: Min radius to use to calculate the com. If None use the minimum radius of the particles.
        :param maxrad: Max radius to use to calculate the com. If None use maximum radius of the particles.
        :param type: Particle type to use to calculate the mean. If None use all the component together
        :return: com: array with the xyz position of the com, vcom: array with vx vy vz of the com.
        """

        if type is None:
            rad_array=self.p.Radius[:]
            vel_array=self.p.Vel[:]
            pos_array=self.p.Pos[:]
            mas_array=self.p.Mass[:]
        else:
            type= nparray_check(type)
            rad_array=self._make_array(self.p.Radius,type)
            vel_array=self._make_array(self.p.Vel,type)
            pos_array=self._make_array(self.p.Pos,type)
            mas_array=self._make_array(self.p.Mass,type)

        if mq is None:
            if minrad is None: minrad=np.min(rad_array)
            if maxrad is None: maxrad=np.max(rad_array)
        else:
            minrad=np.min(rad_array)
            maxrad,_=self.qradius_ext(rad_array,mas_array,mq)


        idx_bound=  (rad_array>=minrad) & (rad_array<=maxrad)

        return self._set_angolar_momentum(pos_array[idx_bound],vel_array[idx_bound],mas_array[idx_bound])

    def _set_angolar_momentum(self,pos_array,vel_array,mass_array):

        lx=mass_array*(pos_array[:,1]*vel_array[:,2]-pos_array[:,2]*vel_array[:,1])
        ly=mass_array*(pos_array[:,2]*vel_array[:,0]-pos_array[:,0]*vel_array[:,2])
        lz=mass_array*(pos_array[:,0]*vel_array[:,1]-pos_array[:,1]*vel_array[:,0])


        return np.sum(lx),np.sum(ly),np.sum(lz)



