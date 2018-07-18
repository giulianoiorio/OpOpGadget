import ctypes as ct
import numpy as np
from numpy.ctypeslib import ndpointer
import os

#pot=np.ascontiguousarray(pot,dtype=float)
#dens=np.ascontiguousarray(dens, dtype=float)

N=1000000
pot=np.zeros(N,order='C')
potgpos=np.zeros(N,order='C')
df=np.zeros(N,order='C')
potgrid=np.zeros(N,order='C')
dfgrid=np.zeros(N,order='C')
vx=np.zeros(N,order='C')
vy=np.zeros(N,order='C')
vz=np.zeros(N,order='C')
v=np.zeros(N,order='C')

dll_name='GenerateModel.so'
dllabspath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + dll_name
lib = ct.CDLL(dllabspath)

df_func=lib.v_gen
df_func.restype=None
df_func.argtypes=[ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ct.c_int,ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS")]
df_func(pot,potgpos,df,potgrid,dfgrid,N,vx,vy,vz,v)
print(vx)
print(vy)
print(vz)
print(np.min(v),np.max(v))