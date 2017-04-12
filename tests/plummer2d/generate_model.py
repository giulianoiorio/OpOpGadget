from OpOp.Model import Plummer2D, Plummer, Sersic
from OpOp.Model import NbodyModel
import numpy as np
from OpOp.analysis import Analysis
import matplotlib.pyplot as plt

#parametri Sculptor da tesi Giuseppina
dist=79
rca=12.7  #primi d'angolo
rc=dist*np.tan(rca/60*(np.pi)/180)
print(rc)
mtot=4.15e6

#mms=GeneralModel(R,dens_t,rc=0.6,G='(kpc km2)/(M_sun s2)',Mmax=1e7)
#mmdm=GeneralModel(R,dens_t,rc=5,G='(kpc km2)/(M_sun s2)',Mmax=1e8)
mms=Plummer(Mmax=mtot,rc=rca/60)


s={'type':2,'model':mms, 'npart':int(1e5)}
#dm={'type':1, 'model':mmdm,'npart':int(1e6)}
#N=int(1e6)

a=NbodyModel([s,])
p=a.generate(use_c=True, po=None,vo=None,mq=70,set_vel=True)



#plt.scatter(p.Pos[:,0],p.Pos[:,1])
#plt.xlim(-5,5)
#plt.ylim(-5,5)
#plt.show()