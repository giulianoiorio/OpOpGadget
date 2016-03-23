from OpOp.Model import Sersic,Plummer, GeneralModel, Tbetamodel,Multimodel, Tbetamodel2
import numpy as np
import matplotlib.pyplot as plt
from OpOp.jsolver import Jsolver
import emcee
import corner

import random
import time

'''
mtots=5e7
s=Sersic(0.71,0.581,mtots)

rc=30
rt=100
rt2=120
mtot=1e11
gamma=1
beta=3
n=512

r=np.linspace(0.001,10,50)

t1=time.time()
t=Tbetamodel(rc=rc,rt=rt,Mmax=mtot,gamma=gamma,beta=beta,n=n)
m=Multimodel(t,s)
j=Jsolver(s.dens,s.sdens,m.mass)
#plt.plot(r,t.mass(r))
plt.plot(r,j.vdisp(r))
print(time.time()-t1)

t1=time.time()
tc=Tbetamodel(rc=rc,rt=rt2,Mmax=mtot,gamma=gamma,beta=beta,n=n,use_c=True)
mc=Multimodel(tc,s)
jc=Jsolver(s.dens,s.sdens,mc.mass)
#plt.plot(r,tc.mass(r))
plt.scatter(r,jc.vdisp(r))
print(time.time()-t1)

plt.show()
'''


def plot_data(data,samples,stellarmod,gamma=0,beta=3,n=512):

    rc=np.percentile(samples[:,0],[16,50,84])
    rt=np.percentile(samples[:,1],[16,50,84])
    mtot=10**np.percentile(samples[:,2],[16,50,84])


    t1=Tbetamodel(rc=rc[0],rt=rt[0],Mmax=mtot[0],gamma=gamma,beta=beta,n=n)
    t2=Tbetamodel(rc=rc[1],rt=rt[1],Mmax=mtot[1],gamma=gamma,beta=beta,n=n)
    t3=Tbetamodel(rc=rc[2],rt=rt[2],Mmax=mtot[2],gamma=gamma,beta=beta,n=n)

    m1=Multimodel(t1,stellarmod)
    m2=Multimodel(t2,stellarmod)
    m3=Multimodel(t3,stellarmod)

    j1=Jsolver(s.dens,s.sdens,m1.mass)
    j2=Jsolver(s.dens,s.sdens,m2.mass)
    j3=Jsolver(s.dens,s.sdens,m3.mass)

    r=np.linspace(0.,2.5,100)
    plt.errorbar(data[:,0],data[:,1],data[:,2],fmt='o')
    plt.plot(r,j1.vdisp(r))
    plt.plot(r,j2.vdisp(r))
    plt.plot(r,j3.vdisp(r))
    plt.xlim(0,2.5)
    plt.savefig=('plot.pdf')

t=Tbetamodel(rc=3.5,rt=6,Mmax=3.05e9,gamma=0,beta=3,n=512)
t2=Tbetamodel(rc=28.4,rt=36.4,Mmax=1.89e10,gamma=1,beta=3,n=512)
mtots=5e7
s=Sersic(0.71,0.581,mtots)
m=Multimodel(t,s)
m2=Multimodel(t2,s)

r=np.logspace(np.log10(0.01),np.log10(2),20)
m=Multimodel(t,s)

j=Jsolver(s.dens,s.sdens,m.mass)
a=j.vdisp(r,mode='F')
out=np.zeros((len(a),3))
out[:,0]=r
noise=np.random.normal(scale=0.5,size=len(a))
out[:,1]=a+noise
out[:,2]=noise
plt.errorbar(out[:,0],out[:,1],out[:,2],fmt='o')

np.savetxt('data.txt',out)



def lnlike(theta,x,y,yerr,stellarmod):

    rc,rt,Mtot=theta
    t=Tbetamodel(rc=rc,rt=rt,Mmax=Mtot,gamma=0,beta=3,n=512,use_c=True)
    totmodel=Multimodel(stellarmod,t)
    j=Jsolver(stellarmod.dens,stellarmod.sdens,totmodel.mass)
    teodisp=j.vdisp(x,mode='F',use_c=False)
    return -np.sum( ((y-teodisp)/(yerr))**2)

def lnprior(theta,x):
    rc,rt,Mtot=theta
    if (0.0001 < rc < rt) and (x[-1]< rt < 50) and ( mtots < Mtot < 1e3*mtots ):
        return 0.0
    else:
        return -np.inf

def lnprob(theta,x,y,yerr,stellarmod):

    lp=lnprior(theta,x)

    if np.isfinite(lp):
        r=lnlike(theta,x,y,yerr,stellarmod)
        if np.isfinite(r):
            return lp + r
        else:
            return -np.inf
    else:
        return -np.inf



ndim,nwalkers=3,100


mm=0.5*(r[-1]+r[1])
posrc=np.random.uniform(mm,mm*5,nwalkers)
posrt=np.random.uniform(posrc,2*posrc,nwalkers)
posmtot=10**(np.random.uniform(np.log10(mtots),np.log10(mtots*1e3),nwalkers))
pos=np.vstack((posrc,posrt,posmtot)).T





sampler=emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(out[:,0],out[:,1],out[:,2],s),threads=3)

t1=time.time()
sampler.run_mcmc(pos,100)
print(time.time()-t1)

samples=sampler.chain[:,20:,:].reshape((-1,ndim))
samples[:,2]=np.log10(samples[:,2])

np.savetxt('samples.txt',samples)

del sampler

fig=plt.figure()
ax1=fig.add_subplot(3,1,1)
ax1.plot(samples[:,0])
ax1.set_ylabel('Rc')
ax2=fig.add_subplot(3,1,2)
ax2.plot(samples[:,1])
ax2.set_ylabel('Rt')
ax3=fig.add_subplot(3,1,3)
ax3.plot(samples[:,2])
ax3.set_ylabel('Mtot')
fig.savefig('chain.pdf')
fig.clf()

fig=corner.corner(samples,labels=["Rc","Rt","Mtot"])
fig.savefig('prova.pdf')

print('Rc',np.percentile(samples[:,0],[16,50,84]))
print('Rt',np.percentile(samples[:,1],[16,50,84]))
print('Mtot',10**np.percentile(samples[:,2],[16,50,84]))

plot_data(out,samples,s,1,3,512)
