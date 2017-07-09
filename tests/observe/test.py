from OpOp.particle import Particle,Particles, Sky_Particles
from OpOp.io import load_snap
from OpOp.analysis import Observe
import matplotlib.pyplot as plt
import numpy as np

#Set a system with two particles

p0=load_snap(filename='out00.bin')
p=load_snap(filename='out10.bin')


#o=Observe(particles=p,com='nor',mq=100)
o=Observe(particles=p)
s=o.observe()

o0=Observe(particles=p0)
s0=o0.observe()



#a=o.locate_sun(dsun=8.0,l_obs=287.534078654, b_obs=-83.1564986815,vlos_obs=111.4, mul_obs=-0.0919737066883, mub_obs=0.00639040515297, d_obs=86 , set_sun=False,com='iter',mq=50)
#print(a)

fig=plt.figure(figsize=(10,5))
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)


ax1.scatter(s.Pos[:,0],s.Pos[:,1],s=0.1)
ax1.set_title('T %.2f Gyr'%s.header['Time'])
ax1.set_xlim(-10,10)
ax1.set_ylim(-10,10)

ax2.scatter(s0.Pos[:,0],s0.Pos[:,1],s=0.1)
ax2.set_title('T %.2f Gyr'%s0.header['Time'])
ax2.set_xlim(-10,10)
ax2.set_ylim(-10,10)

plt.show()

'''
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
'''