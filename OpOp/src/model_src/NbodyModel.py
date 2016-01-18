from ..particle_src.particle import Header,Particles
from ..model_src.Model import Model
from ..grid_src.grid import grid
from ..df_src.spherical import df_isotropic
from ..analysis_src import analysis

import numpy as np
import ctypes as ct
from numpy.ctypeslib import ndpointer
import math as mt
import os
import  time


class NbodyModel():

    def __init__(self,components,dff=df_isotropic, Ngrid=512, xmin=2E-3, xmax=200, kind='log'):
        """
        Class to Generate a Nbody models from dynamic models from the class Model.
        :param components: Array of the dynamic components, each element need to be an array with the following keys:
                            -type: integer from 0 to 5, it is the type of particles to be used for the component, the numbers follows the Gadget convention (0-Gas, 1-DM Halo, 2-Disk, 3-Bulge, 4-Star, 5-Bndry)
                            -model: An instance of the class Model
                            -npart: Number of particles to use
                            At the end of the instance of the class the dictonary will have other keys as 'grid' that is the grid used to sample pot and dens
                            'pindex', that stores the first and last index where to find particles of the component in the particle instance made and
                            'id', that is a the id number of the dynamical model.
        :param dff: density distribution calculator, need to be a function that get in input the density and totpot grid and gives in output the
                    used tpotgrid, the df grid, and df interpolated function.

        !Grid params
        These are the parameters to set the grids where evaluate dens, mass and tot pot for the components
        :param Ngrid: Number of sampling points, the standard is 512
        :param xmin: Minimum value of the grid in unit of R/Rc.
        :param xmax: Maximum value of the grid i unit of R/Rc.
        :param kind: Kind of spaced sampling: linear or log.
        The Class produce a grid instance for every model in the component array from a minimum  equals to xmin*model.rc to a maximum equals to xmax*model.rc

        !Method:
        The only method to be used from the user is:
        -method generate:
            Use generate to generate the initial conditions in the phase space for the  dynamical model in input
            :param use_c: If use_c=True, the velocities will be generated with a fast C exstension. This
                          is strongly recommended for components with Npart > 1e5, in fact in this case the
                          pure Python implementation can be very slow (>10 min for 1e6, >100 min for 1e7). With use_c=True should be possibile
                          to generate up to 1e8 particles in the reasonable time (<5 min).
                          The time for pure Python version scale more or less like (N/1e4)*6 s
                          For the C versione like (N/1e6)*3 s.
                          This time are calibrated for a single component, for system with more than 1 dynamical component the times could
                          considerably larger, exspecially in the pure Python version.
            :return p: Return the updated p instance of the Particles Class (it fills all the field, also Pot, Radius, Vel_tot ).


        !Important variables:
        :var p: Instance of the Particle class, it can contains only zeros in the case the method generate has still not been called
        :var components: Array with one dictionary for each dynamical models used. The key in the dictionary will be (the values signed with
        [input] are the value given in input to the class):
            :key type [input]: integer from 0 to 5, it is the type of particles to be used for the component, the numbers follows the Gadget convention (0-Gas, 1-DM Halo, 2-Disk, 3-Bulge, 4-Star, 5-Bndry)
            :key model [input]: The dynamical mode, an instance of the class Model
            :key npart [input]: Number of particles.
            :key id: that is a the id number of the dynamical model.
            :key grid: grid used to sample pot and dens, it is an instance of the class grid
            :key pindex: 'pindex', that stores the first and last index where to find particles of the component in the particle instance made
        """

        self.Ngrid=Ngrid
        self.xmin=xmin
        self.xmax=xmax
        self.kind=kind

        #Check components
        if isinstance(components,dict):
            self._check_component(components)
            self.components=(components,)

        elif hasattr(components, '__iter__'):
            i=0
            for c in components:
                self._check_component(c,i=i)
                i+=1
            components=sorted(components, key= lambda k: k['type']  )
            self.components=tuple(components)



        self.df=dff
        self._set_header()
        self._set_particles()
        for c in self.components: self._set_grid(c)

    def generate(self, use_c=True, po=None,vo=None,mq=70):
        """
        Generate Position and velocities.
        The position are randomly picken from the cumulative mass distribution, while
        the velocity component are taken sampling the distribution function.
        The output system will be spherically symmetric and isotropic.
        :param use_c: If use_c use a fast implementation to generate velocity
        :return:
        """
        print('***Generate ICS: Start***')

        for c in self.components:
            print('-Component id:%i type:%i Npart=%i'%(c['id'],c['type'],c['npart']),flush=True)
            t1=time.time()
            print('     Generate Positions:',end=' ',flush=True)
            self._set_position(c)
            print('     Done',flush=True)
            print('     Generate Velocities:',end=' ',flush=True)
            self._set_vel(c,use_c=use_c)
            print('     Done',flush=True)
            print('     Done in %.3f'%(time.time()-t1),flush=True)

        if (po is not None)    or (vo is not None):
            if po is None: po=(0,0,0)
            if vo is None: vo=(0,0,0)
            print('     Calculate COM and move it to (%.1f,%.1f,%.1f) with V (%.1f,%.1f,%.1f):'%(po[0],po[1],po[2],vo[0],vo[1],vo[2]),end=' ',flush=True)
            #Find the CM and set it at the position po with velocity vo
            analysis.Analysis(self.p,safe=False, mq=mq, auto_centre=True,po=po,vo=vo)
            print('     Done',flush=True)

        return self.p

    def _set_header(self):
        """
        Set header for the class particles
        :return:
        """

        self.h=Header()
        self._memindex=[] #serve per tenere conto in che indici vanno posizionate le varie componenti,
                          #questa info è uguale a quella ottenuta da Nall nel caso le componenti inserite siano
                          #di Type diversi, ma è essenziale nel caso ci siano diverse componenti con lo stesso Type
                          #che Nall non distingue.

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

    def _set_particles(self):
        """
        Create an instance of the particle object
        :return:
        """

        self.p=Particles(h=self.h)

    def _check_component(self,component,i=0):

        fail=[]

        if 'type' not in component: fail.append('type')
        if 'model' not in component: fail.append('model')
        elif isinstance(component['model'],Model)==False: fail.append('model')
        if 'npart' not in component: fail.append('npart')

        if len(fail)>0: raise ValueError('Missed keyword in component', i ,':',fail)

    def _set_grid(self,component):
        """
        Set grid for the component c
        :param component:
        :return:
        """

        rc=component['model'].rc

        ext_list=[]
        for c in self.components:
            if c['id']!=component['id']: ext_list.append(c['model'])

        if len(ext_list)==0: ext_list=None


        c_grid=grid(N=self.Ngrid, galaxy_model=component['model'], ext_pot_model=ext_list, type=self.kind, min=rc*self.xmin, max=rc*self.xmax )
        component['grid']=c_grid

    def _set_position(self,component):
        """
        Set position from the cumulative distribution
        :param component:
        :return:
        """

        g=component['grid']
        maxmass=np.max(g.mgrid)
        minmass=np.min(g.mgrid)
        Nrand=component['npart']
        mass_per_part=component['model'].Mmax/Nrand

        #Random generation in spherical coordinate
        u=np.random.uniform(minmass,maxmass,size=int(Nrand)) #for Rad
        radius=g.eval_rad(u)
        theta=np.random.uniform(-1,1,size=int(Nrand)) #for theta angle (for a spheril simmetry need to ve uniform in cos(thteta))
        phi=np.random.uniform(0,2*np.pi,size=int(Nrand)) #pick random for phi angle from 0 to 360 degree

        #Transform in cartesian coordinates
        x=radius*(np.sqrt(1-theta*theta)*np.cos(phi)) #R (sin(teta) cos(phi))
        y=radius*(np.sqrt(1-theta*theta)*np.sin(phi)) #R (sin(teta) sin(phi))
        z=radius*theta #R (cos(teta))



        idxmax=component['pindex'][1]
        idxmin=component['pindex'][0]

        self.p.Radius[idxmin:idxmax]=radius
        self.p.Pos[idxmin:idxmax,0]=x
        self.p.Pos[idxmin:idxmax,1]=y
        self.p.Pos[idxmin:idxmax,2]=z
        self.p.Mass[idxmin:idxmax]=mass_per_part
        self.p.Pot[idxmin:idxmax]=g.eval_tpot(radius)

        return 0

    def _set_vel(self,component,use_c=True):
        """
        Set velocity for the component
        :param component:
        :param use_c: Fast implementation in c
        :return:
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

    def _v_extract(self,pot,df_func):
        """
        Generate velocity using acceptance-rejection algorithm.
        It generate a sphrical symmetric velocity distribution.
        :param pot: Potential of the particle
        :param df: distribution function, function of e
        :return:
        """
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

    @staticmethod
    def _v_extract_c(pot_grid,df_grid,df_func,particle_pot):
        """
        Generate velocity using acceptance-rejection algorithm.
        It generate a sphrical symmetric velocity distribution.
        :param pot: Potential of the particle
        :param df: distribution function, function of e
        :return:
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





        dll_name='model_c_ext/GenerateModel.so'
        dllabspath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + dll_name
        lib = ct.CDLL(dllabspath)

        v_gen=lib.v_gen
        v_gen.restype=None
        v_gen.argtypes=[ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_int, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ct.c_int,ct.c_int,ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS")]

        v_gen(pot,potgpos,dfmax,potgrid,dfgrid,N,ngrid,vx,vy,vz,v)



        return vx,vy,vz,v









if __name__=='main':
    R,mass_t,dens_t,pot_t=np.loadtxt('stellarcomp.txt',unpack=True,comments='#')








    mms=GeneralModel(R,dens_t,rc=0.6,G='(kpc km2)/(M_sun s2)',Mmax=1e7)
    mmdm=GeneralModel(R,dens_t,rc=1,G='(kpc km2)/(M_sun s2)',Mmax=1e8)



    s={'type':2,'model':mms, 'npart':int(1e3)}
    dm={'type':1, 'model':mmdm,'npart':int(1e4)}
    N=int(1e4)

    a=NbodyModel([dm],df_isotropic)
    print(a.components)

    outtab=np.zeros(shape=(20,8))
    for i in range(20):
        t1=time.time()
        a.generate()
        t2=time.time()
        #g=a.components[0]['grid']
        #print(np.std(a.p.Radius[:N]))
        #print(np.std(a.p.Radius[N:]))
        print(np.std(a.p.Vel[:N,0]),np.std(a.p.Vel[:N,1]),np.std(a.p.Vel[:N,2]),np.std(a.p.Vel_tot[:N]),np.mean(a.p.Vel[:N,0]))
        #print(np.std(a.p.Vel[N:,0]),np.std(a.p.Vel[N:,1]),np.std(a.p.Vel[N:,2]),np.std(a.p.Vel_tot[N:]),np.mean(a.p.Vel[N:,0]))
        outtab[i,0]=np.std(a.p.Vel[:N,0])
        outtab[i,1]=np.std(a.p.Vel[:N,1])
        outtab[i,2]=np.std(a.p.Vel[:N,2])
        outtab[i,3]=t2-t1

        t1=time.time()
        a.generate(use_c=False)
        t2=time.time()
        #N=int(1e5)
        #print(np.std(a.p.Radius[:N]))
        #print(np.std(a.p.Radius[N:]))
        print(np.std(a.p.Vel[:N,0]),np.std(a.p.Vel[:N,1]),np.std(a.p.Vel[:N,2]),np.std(a.p.Vel_tot[:N]),np.mean(a.p.Vel[:N,0]))
        outtab[i,4]=np.std(a.p.Vel[:N,0])
        outtab[i,5]=np.std(a.p.Vel[:N,1])
        outtab[i,6]=np.std(a.p.Vel[:N,2])
        outtab[i,7]=t2-t1
    np.savetxt('out_gen_1e4.txt',outtab)




    #print(np.std(a.p.Vel[N:,0]),np.std(a.p.Vel[N:,1]),np.std(a.p.Vel[N:,2]),np.std(a.p.Vel_tot[N:]),np.mean(a.p.Vel[N:,0]))

    #a.generate(use_c=False)
    #print(np.std(a.p.Vel[:N,0]),np.std(a.p.Vel[:N,1]),np.std(a.p.Vel[:N,2]),np.std(a.p.Vel_tot[:N]),np.mean(a.p.Vel[:N,0]))
    #print(np.std(a.p.Vel[N:,0]),np.std(a.p.Vel[N:,1]),np.std(a.p.Vel[N:,2]),np.std(a.p.Vel_tot[N:]),np.mean(a.p.Vel[N:,0]))






    '''
            k=max_indx;
            while(e<=grid_pot[k])
            {
                k=k+1;
            }
    '''

    '''
    i=0
    while v>1 or ch==0:
        vx,vy,vz=np.random.uniform(-1,1,size=3)
        v=mt.sqrt(vx*vx+vy*vy+vz*vz)
        e=pot*(1-v*v)
        print(v)
        if v<1:
            u=np.random.random()
            df=bf(e)/bf(pot)
            print('pot',pot,'v',v,'u',u,'e',e,'df',df,'df(e)',bf(e),'dfmax',bf(pot))
            if u<=df: ch=1
            else: ch=0
        i+=1
    print(v,i)
    '''

    '''
    def v_extract(pot,df):
        """
        Generate velocity using trial-error algorithm.
        It generate a sphrical symmetric velocity distribution.
        :param pot: Potential of the particle
        :param df: distribution function
        :return:
        """
        v=2
        ch=0

        while v>1 or ch==0:
            vx,vy,vz=np.random.uniform(-1,1,size=3)
            v=mt.sqrt(vx*vx+vy*vy+vz*vz)
            e=pot*(1-v*v)

            if v<=1:
                u=np.random.random()
                umax=df(e)/df(pot)
                if u<=umax: ch=1
                else: ch=0

        norm=np.sqrt(2*pot)

        vx=vx*norm
        vy=vy*norm
        vz=vz*norm
        v=v*norm


        return vx,vy,vz,v

    f=np.vectorize(v_extract,otypes=[np.float,np.float,np.float,np.float])



    t1=time.time()
    #vx,vy,vz,v=f(pot,bf)
    t2=time.time()

    #print(out)

    #val,bins=np.histogram(vx,bins=100)
    #plt.step(0.5*(bins[:-1]+bins[1:]),val)

    #plt.show()
    print(np.std(vx),np.std(vy),np.std(vz),np.std(v),np.mean(vx))
    print('tempo',t2-t1)
    '''
    '''
    #print(a.components)
    print(a.components[0]['grid'].dgrid)
    print(a.components[0]['grid'].tpgrid)

    bx,by,bf=df_isotropic(a.components[0]['grid'].dgrid,a.components[0]['grid'].tpgrid)

    plt.plot(bx,np.log10(bf(bx)/bf(10)))
    plt.ylim(-5,5)
    plt.show()
    '''

    '''
    mod1=a.components[0]['model']

    mod2=a.components[1]['model']




    rad=a.p.Radius[:10]
    pot=a.p.Pot[:10]
    print('R PotDm PotStar Pottot Calc ')
    i=0
    for r in rad:
        print('%.3f  %.3f  %.3f  %.3f %.3f'%(r,mod1.pot(r),mod2.pot(r),mod2.pot(r)+mod1.pot(r),pot[i]))
        i+=1


    rad=a.p.Radius[int(1e5):int(1e5)+10]
    pot=a.p.Pot[int(1e5):int(1e5)+10]
    print('R PotDm PotStar Pottot Calc ')
    i=0
    for r in rad:
        print('%.3f  %.3f  %.3f  %.3f %.3f'%(r,mod1.pot(r),mod2.pot(r),mod2.pot(r)+mod1.pot(r),pot[i]))
        i+=1
    '''

    '''
    R=np.linspace(0.001,100,1000)
    #print(mod1.R,mod2.R)
    plt.plot(R,mod1.pot(R)/mod2.pot(R))
    plt.show()
    '''

    '''
    print(a.components)
    print(np.max(a.p.Radius[:int(1e5)]))
    print(np.max(a.p.Radius[int(1e5):]))
    x=a.p.Pos[:int(1e5),0]
    y=a.p.Pos[:int(1e5),1]
    z=a.p.Pos[:int(1e5),2]
    plt.scatter(x,z,c='black',alpha=0.01)

    x=a.p.Pos[int(1e5):,0]
    y=a.p.Pos[int(1e5):,1]
    z=a.p.Pos[int(1e5):,2]
    plt.scatter(x,z,c='red',alpha=0.5)

    plt.axis('equal')

    plt.savefig('halo.pdf')
    '''

    '''
    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')

    x=a.p.Pos[:int(1e5),0]
    y=a.p.Pos[:int(1e5),1]
    z=a.p.Pos[:int(1e5),2]
    ax.scatter(x, y, z, c='black', marker='o',alpha=0.05)


    x=a.p.Pos[int(1e5):,0]
    y=a.p.Pos[int(1e5):,1]
    z=a.p.Pos[int(1e5):,2]
    ax.scatter(x, y, z, c='red', marker='o')

    ax.set_xlim(-4,4)
    ax.set_ylim(-4,4)
    ax.set_zlim(-4,4)

    plt.show()
    '''

    '''
    #prova per settare i raggi
    for i in range(10):
        npart=1e6
        mpart=1/npart

        mnorm=g.mgrid/np.max(g.mgrid)

        max_m=np.max(g.mgrid)
        min_m=np.min(g.mgrid)

        #u random tra 0 e 1
        #random tra [min_m, max_m) : (max-min)*u + min

        masscum_teo=mm.mass(g.gx)/np.max(mm.mass(g.gx))
        mass_idx=np.random.uniform(size=int(npart))

        pos=g.eval_rad(mass_idx)


        #mass_idx=np.random.randint(0,len(g.gx),int(npart))
        #pos=g.gx[mass_idx]

        #print(np.max(mass_idx),g.eval_rad(np.max(mass_idx)))
        #print(np.min(mass_idx),g.eval_rad(np.min(mass_idx)))

        masshist,bine=np.histogram(pos,bins=g.gedge)
        masscum=np.cumsum(masshist)
        masscum_mod_1=masscum*mpart
        if i==0: plt.plot(g.gx,np.log10(np.abs((masscum_teo-masscum_mod_1)/masscum_teo)),c='red',label='Npart 1e6')
        else: plt.plot(g.gx,np.log10(np.abs((masscum_teo-masscum_mod_1)/masscum_teo)),c='red')


        #Dens
        npart=1e3
        mpart=1/npart

        mnorm=g.mgrid/np.max(g.mgrid)

        max_m=np.max(g.mgrid)
        min_m=np.min(g.mgrid)

        #u random tra 0 e 1
        #random tra [min_m, max_m) : (max-min)*u + min


        mass_idx=np.random.uniform(size=int(npart))


        pos=g.eval_rad(mass_idx)

        print(np.max(mass_idx),g.eval_rad(np.max(mass_idx)))
        print(np.min(mass_idx),g.eval_rad(np.min(mass_idx)))

        masshist,bine=np.histogram(pos,bins=g.gedge)
        masscum=np.cumsum(masshist)

        #plt.step(g.gx,masshist*mpart)
        #plt.xlim(0,7)
        #plt.show()

        masscum_teo=mm.mass(g.gx)/np.max(mm.mass(g.gx))

        masscum_mod_1=masscum*mpart






        #plt.plot(g.gx,masscum*mpart)
        #plt.plot(g.gx,masscum_teo)
        if i==0: plt.plot(g.gx,np.log10(np.abs((masscum_teo-masscum_mod_1)/masscum_teo)),c='blue', ls='dashed', label='Npart 1e3')
        else: plt.plot(g.gx,np.log10(np.abs((masscum_teo-masscum_mod_1)/masscum_teo)),c='blue', ls='dashed')

    plt.xlim(-0.2,14)
    plt.xlabel('R/Rc',fontsize=18)
    plt.ylabel('$\log (\epsilon_M)$',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.suptitle('Cumulative mass',fontsize=20)
    plt.legend()
    plt.savefig('cumulative_mass_test.pdf')
    #plt.show()
    '''
    '''
    for i in range(10):
        npart=1e3
        mpart=1/npart
        mnorm=g.mgrid/np.max(g.mgrid)
        max_m=np.max(g.mgrid)
        min_m=np.min(g.mgrid)
        mass_idx=np.random.uniform(size=int(npart))
        pos=g.eval_rad(mass_idx)
        masshist,bine=np.histogram(pos,bins=g.gedge)
        dentm=masshist/g.g_vol




        dentm_f=UnivariateSpline(g.gx,dentm,s=0,k=1)
        dentt_f=UnivariateSpline(R,dens_t,s=0,k=1)


        dens_mod1=dentm_f(g.gx)/dentm_f(1)
        dens_teo=dentt_f(g.gx)/dentt_f(1)
        diff=np.abs(dens_mod1-dens_teo)/dens_teo
        if i==0: plt.plot(g.gx, np.log10(diff),c='blue', ls='dashed', label='Npart 1e3' )
        else: plt.plot(g.gx, np.log10(diff),c='blue', ls='dashed' )


        npart=1e6
        mpart=1/npart
        mnorm=g.mgrid/np.max(g.mgrid)
        max_m=np.max(g.mgrid)
        min_m=np.min(g.mgrid)
        mass_idx=np.random.uniform(size=int(npart))
        pos=g.eval_rad(mass_idx)
        masshist,bine=np.histogram(pos,bins=g.gedge)
        dentm=masshist/g.g_vol




        dentm_f=UnivariateSpline(g.gx,dentm,s=0,k=1)
        dentt_f=UnivariateSpline(R,dens_t,s=0,k=1)


        dens_mod1=dentm_f(g.gx)/dentm_f(1)

        dens_teo=dentt_f(g.gx)/dentt_f(1)
        diff=np.abs(dens_mod1-dens_teo)/dens_teo


        if i==0: plt.plot(g.gx, np.log10(diff),c='red', label='Npart 1e6' )
        else:  plt.plot(g.gx, np.log10(diff),c='red')





    plt.xlabel('R/Rc',fontsize=18)
    plt.ylabel(r'$\log (\epsilon_\rho) $',fontsize=20)
    plt.tick_params(axis='both', which='major', labelsize=16)
    plt.suptitle('Density profile',fontsize=20)
    plt.xlim(0,10)
    plt.ylim(-4,0.5)
    plt.legend(loc='lower right')
    plt.savefig('density_profile_test.pdf')
    #plt.show()
    '''
