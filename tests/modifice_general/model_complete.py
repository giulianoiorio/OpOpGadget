from OpOp.Model import Plummer2D, Plummer, Sersic, Isothermal, Tbetamodel
from OpOp.Model import NbodyModel
from OpOp.analysis import Analysis, Profile
from OpOp.io import write_snap
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
label_size =20
mpl.rcParams.update({'figure.autolayout':True})
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 
mpl.rcParams['mathtext.default']='regular'
from astropy.constants import G as conG

def iso_mass(r,rc,rm,Mm):
	"""
	Per ottenere la massa a un certo raggio M, conoscendo Rc e la massa ad un altro raggio
	Mm
	"""
	num=r-rc*np.arctan(r/rc)
	den=rm-rc*np.arctan(rm/rc)
	
	return Mm *num/den

def iso_vinf(rc,rm,Mm,G='kpc km2 / (M_sun s2)'):
    
	GG=conG.to(G)
	G=GG.value
	
	num=G*Mm
	den=rm-rc*np.arctan(rm/rc)
	
	return np.sqrt(num/den)

#parametri Sculptor da tesi Giuseppina
dist=79
#Star
rca=12.7  #primi d'angolo
rc_star=dist*np.tan(rca/60*(np.pi)/180)
mstar=4.15e6
#halo
#NFW-core
rc_halo=0.3258 
rt_halo=15
gamma_halo=0
beta_halo=3
mm_halo=10**8.7161#Mass at rm_halo
rmax=30 #rmax halo
rmin=0.001
rini=rmin/rc_halo
rfin=rmax/rc_halo




ms=Plummer(Mmax=mstar,rc=rc_star)
mh=Tbetamodel(rc=rc_halo,rt=rt_halo,Mmax=mm_halo,gamma=gamma_halo,beta=beta_halo,rini=rini,rfin=rfin)


s={'type':2,'model':ms, 'npart':int(5e4)}
dm={'type':1, 'model':mh,'npart':int(1e6)}
a=NbodyModel([s,dm],xmax=rmax/rc_halo)
p=a.generate(use_c=True, po=None,vo=None,mq=70,set_vel=True)
aa=Analysis(p,safe=False,auto_centre=True)
#write_snap(p,'Sculptor_plummer_beta3.dat',end='<',enable_mass=False,safe_write=True)

print(p.Pos[0])
print(p.Vel[0])

'''
prof=Profile(p,xmin=0.001,xmax=10,ngrid=512,kind='lin',type=2)



arr=prof.vdisp2d(pax='x',ret=True,func=True,s=None)[0]
r=arr[:,0]
vd=arr[:,1]
plt.plot(r,vd,label='Vdx 2D',c='red')

arr=prof.vdisp2d(pax='y',ret=True,func=True,s=None)[0]
r=arr[:,0]
vd=arr[:,1]
plt.plot(r,vd,label='Vdy 2D',c='orange')

arr=prof.vdisp2d(pax='z',ret=True,func=True,s=None)[0]
r=arr[:,0]
vd=arr[:,1]
plt.plot(r,vd,label='Vdz 2D',c='magenta')

plt.legend(loc='best')
plt.xlim(0.01,2)
#plt.xscale('log')

plt.xlabel('r/R [kpc]',fontsize=20)
plt.ylabel('$\sigma \  [km/s]$',fontsize=20)

plt.show()
'''







