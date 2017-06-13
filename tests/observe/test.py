from OpOp.particle import Particle,Particles
from OpOp.analysis import Observe
import matplotlib.pyplot as plt
import numpy as np

#Set a system with two particles

pp=[]

pcom=(5.41,-9.77,-85.39)
vcom=(20.98,202.21,-103.13)
for i in range(1000):
    pos=np.random.normal(loc=pcom,scale=1)
    vel=np.random.normal(loc=vcom,scale=1)
    pp.append(Particle(id=i,type=2,pos=pos,vel=vel,mass=1))



p=Particles(p=pp)


#plt.scatter(p.Pos[:,0],p.Pos[:,2])
#plt.show()

#o=Observe(particles=p,com='nor',mq=100)
o=Observe(particles=p)
print(o.align_vec)
print(o.dist_obj)

print('Median')
a,b,c=np.median(o.obsVel,axis=0)
print(a / (4.74047*86), b / (4.74047*86) )
print('#')


vo=o.vobs(o.sPos,o.sVel)

print(np.median(vo,axis=0))
plt.scatter(vo[:,0],vo[:,1])
plt.show()