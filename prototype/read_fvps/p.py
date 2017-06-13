from OpOp.io import load_snap
from OpOp.analysis import Analysis, Profile
import matplotlib.pyplot as plt
import matplotlib as mpl
import struct
import numpy as np
'''
t=np.loadtxt('data.txt')
mass=t[:,0]
id=t[:,4]
pos=t[:,1:4]
vel=t[:,5:8]
ipar=np.array([0,0,len(t),2,1,0])
rpar=np.array([0,5,1,0,0.5,0,0,0,1])
print('pos',pos)
print(ipar)
print(rpar)

pos_c=np.array(pos, dtype=np.dtype('=f'), order='C')
vel_c=np.array(vel, dtype=np.dtype('=f'), order='C')
mass_c=np.array(mass, dtype=np.dtype('=f'), order='C')
id_c=np.array(id, dtype=np.dtype('=f'), order='C')
ipar_c=np.array(ipar,dtype=np.dtype('=i'))
rpar_c=np.array(rpar,dtype=np.dtype('=f'))

filename='tol.bin'

wr.write_file_c(filename, ipar_c, rpar_c, pos_c, vel_c, id_c, mass_c)
'''
filename='out00.bin'

stream=open(filename,'rb')
#struct.unpack("i", stream.read(4))
ac=struct.unpack("i", stream.read(4))[0]
bl=int(ac/4)
print(ac)
a= struct.unpack("i"*bl, stream.read(ac))  # read the initial block checek (C-type int)
ac=struct.unpack("i", stream.read(4))
print(ac)
print(a)

ac=struct.unpack("i", stream.read(4))[0]
bl=int(ac/4)
print(ac)
a= struct.unpack("f"*bl, stream.read(ac))  # read the initial block checek (C-type int)
struct.unpack("i", stream.read(4))
print(ac)
print(a)


ac=struct.unpack("i", stream.read(4))[0]
bl=int(ac/4)
print(ac)
a= struct.unpack("f"*bl, stream.read(ac))  # read the initial block checek (C-type int)
struct.unpack("i", stream.read(4))
print(ac)
print(a)

ac=struct.unpack("i", stream.read(4))[0]
bl=int(ac/4)
print(ac)
a= struct.unpack("f"*bl, stream.read(ac))  # read the initial block checek (C-type int)
struct.unpack("i", stream.read(4))
print(ac)
print(a)
#print(p.Pos)
#print(p.Id)
#print(p.Type)
#print(p.Radius)
#print(p.header)
#print(p[3])
#stream.close()
#h=load_header(filename)
#print(h)

p=load_snap(filename)
#idxh=p.Type[:]==1
#x=p.Pos[idxh,0]
#y=p.Pos[idxh,1]
#plt.scatter(x,y,s=0.0001)
#idxs=p.Type[:]==2
#x=p.Pos[idxs,0]
#y=p.Pos[idxs,1]
#plt.scatter(x,y,s=10,c='red',zorder=10000)
#plt.show()

h=p.header
print(h)
print(p[0])
#print(np.sum(idxh),np.sum(idxs))

#a0=Analysis(p,safe=True,auto_centre=True,iter=True,single=False)
#prof0=Profile(a0.p,xmin=0.001,xmax=12,ngrid=50,kind='lin',type=2)
#arr=prof0.vdisp2d(pax='x')[0]
#r=arr[:,0]
#vd=arr[:,1]
#plt.plot(r,vd,lw=3,label='T=0 Gy')
#plt.xscale('log')
#plt.yscale('log')
#plt.show()


#print(a0.com(mq=90,type=2))
#print(a0.com(mq=90,type=1))
#print(a0.com(mq=90))



"""
dic,ida,posa,vela,massa=cr.read_file_c('tol.bin')

print(dic)
print(np.asarray(ida))
print(np.asarray(massa))
print(np.asarray(posa))
print(np.asarray(vela))
"""
