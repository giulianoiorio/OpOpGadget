import numpy as np
#from cython.parallel cimport prange
from cython import boundscheck,wraparound,cdivision
from math import sqrt


cdef extern from "math.h" nogil:
    double c_sqrt "sqrt" (double)

@cdivision(True)
cdef double _kproj_iso(double r,double R) nogil:
    """
    Mamon and Lokas, 2005, A16-isotropic
    :param r: spherical radius
    :param R: cylindircal radius
    :return: kproj
    """

    return c_sqrt(1-(R*R)/(r*r))

@boundscheck(False)
@wraparound(False)
@cdivision(True)
cdef void _integrateiso( double[:] ret, int rlen ,double[:] r, double[:] dens, double[:] sdens, double[:] mass,int cpar):
    """

    :param r:
    :param dens:
    :param mass:
    :return:
    """

    cdef:
        unsigned int i, k, n=rlen


    if cpar==1:
        for k in range(n-1,nogil=True):


            for i in range(k,n):


                if (i==k): ret[k] +=  0.5*(r[i+1]-r[i])*( (_kproj_iso(r[i],r[k])*dens[i]*mass[i]/r[i]) + (_kproj_iso(r[i+1],r[k])*dens[i+1]*mass[i+1]/r[i+1]) )
                else: ret[k] +=  0.5*(r[i]-r[i-1])*( (_kproj_iso(r[i],r[k])*dens[i]*mass[i]/r[i]) + (_kproj_iso(r[i-1],r[k])*dens[i-1]*mass[i-1]/r[i-1]) )


            ret[k]=c_sqrt(ret[k]/sdens[k])
    else:
        for k in range(n-1):

            for i in range(k,n):


                if (i==k): ret[k] +=  0.5*(r[i+1]-r[i])*( (_kproj_iso(r[i],r[k])*dens[i]*mass[i]/r[i]) + (_kproj_iso(r[i+1],r[k])*dens[i+1]*mass[i+1]/r[i+1]) )
                else: ret[k] +=  0.5*(r[i]-r[i-1])*( (_kproj_iso(r[i],r[k])*dens[i]*mass[i]/r[i]) + (_kproj_iso(r[i-1],r[k])*dens[i-1]*mass[i-1]/r[i-1]) )


            ret[k]=c_sqrt(ret[k]/sdens[k])

    ret[-1]=ret[-2]





def _integratec(G,dprof,sprof,mprof,R=None,rini=0,rfin=10,n=512,kind='log',cpar=True):

    cost=sqrt(2*G)

    if R is None:
        if kind=='log': R=np.logspace(np.log10(rini+0.001),np.log10(rfin+0.001),n)-0.001
        elif kind=='lin': R=np.linspace(rini,rfin,n)

    dens=dprof(R)
    mass=mprof(R)
    sdens=sprof(R)

    ret=np.zeros(len(R),dtype=np.float64,order='C')

    if cpar==True: _integrateiso(ret,len(ret) ,R,dens,sdens,mass,1)
    else: _integrateiso(ret,len(ret) ,R,dens,sdens,mass,0)

    ret=np.array(ret)*cost

    return R,ret

@boundscheck(False)
@wraparound(False)
def _evaluate(func,double[:] R,double Rcut):

    cdef:
        int i,n=len(R)
        double[:] ret=np.zeros(n,dtype=np.float64)

    for i in range(n):
        ret[i]=func(R[i],Rcut)

    return np.array(ret)

