from OpOp.particle import Particle,Particles, Sky_Particles
from OpOp.io import load_snap
from OpOp.analysis import Observe
import matplotlib.pyplot as plt
import numpy as np

#Set a system with two particles

p=load_snap(filename='out10.bin')


#o=Observe(particles=p,com='nor',mq=100)
o=Observe(particles=p)
s=o.observe()
print(s.header)

a=o.locate_sun(dsun=8.5,l_obs=287.534078654, b_obs=-83.1564986815,vlos_obs=111.4, mul_obs=-0.0919737066883, mub_obs=0.00639040515297, d_obs=86 , set_sun=False,com='iter',mq=50)
print(a)


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