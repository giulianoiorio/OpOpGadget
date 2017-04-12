from OpOp.Model import Plummer2D, Plummer
import numpy as np
import matplotlib.pyplot as plt

p=Plummer(1,1e6)
p2=Plummer2D(1,1e6)

R=np.logspace(np.log10(0.0001),np.log10(5),1000)

densp=p.dens(R)
sdensp=p.sdens(R)
densp2=p2.dens(R)
sdensp2=p2.sdens(R)


plt.plot(R,densp/p.dens(0.1),label='Plummer')
plt.plot(R,densp2/p2.dens(0.1),label='Plummer2d')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()