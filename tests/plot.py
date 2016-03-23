import numpy as np
import matplotlib.pyplot as plt
from OpOp.Model import Sersic,Plummer, GeneralModel, Tbetamodel,Multimodel
from OpOp.jsolver import Jsolver


data=np.loadtxt('data.txt')
samples=np.loadtxt('samples.txt')

rc=np.percentile(samples[:,0],[16,50,84])
rt=np.percentile(samples[:,1],[16,50,84])
mtot=10**np.percentile(samples[:,2],[16,50,84])
print(mtot)
t1=Tbetamodel(rc=rc[0],rt=rt[0],Mmax=mtot[0],gamma=1,beta=3,n=512)
t2=Tbetamodel(rc=rc[1],rt=rt[1],Mmax=mtot[1],gamma=1,beta=3,n=512)
t3=Tbetamodel(rc=rc[2],rt=rt[2],Mmax=mtot[2],gamma=1,beta=3,n=512)

mtots=5e7
s=Sersic(0.71,0.581,mtots)
m1=Multimodel(t1,s)
m2=Multimodel(t2,s)
m3=Multimodel(t3,s)

j1=Jsolver(s.dens,s.sdens,m1.mass)
j2=Jsolver(s.dens,s.sdens,m2.mass)
j3=Jsolver(s.dens,s.sdens,m3.mass)

r=np.linspace(0.,2.5,100)
plt.errorbar(data[:,0],data[:,1],data[:,2],fmt='o')
plt.plot(r,j1.vdisp(r))
plt.plot(r,j2.vdisp(r))
plt.plot(r,j3.vdisp(r))
plt.xlim(0,2.5)
plt.show()