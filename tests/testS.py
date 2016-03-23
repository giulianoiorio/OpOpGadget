from OpOp.Model import Sersic,Plummer, GeneralModel, Tbetamodel, Isothermal
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma,gammaincc
from OpOp.jsolver import Jsolver
import time

r=np.linspace(0.01,30,100)

rc=0.5
i1=Isothermal(rc,Mmax=1e9)
i2=Isothermal(rc,Mmax=1e9,rmax=200)

j1=Jsolver(i1.dens,i1.sdens,i1.mass)

plt.ylim(0,25)
plt.plot(r,j1.vdisp(r))
#plt.plot(r,i1.sdens(r))
#plt.plot(r,i2.sdens(r,100))
plt.show()