import numpy as np
cimport numpy as np
cimport cython
from libc.stdio cimport FILE, fopen, fread, stdout, fclose


cdef read_block_data(FILE *buf, int Nobj):

    cdef float[:, ::1] pos_block = np.empty((Nobj,3),  order='C', dtype=np.dtype('=f'))
    cdef float[:, ::1] vel_block = np.empty((Nobj,3),  order='C', dtype=np.dtype('=f'))
    cdef float[::1] id_block = np.empty((Nobj,),  order='C', dtype=np.dtype('=f'))
    cdef float[::1] mass_block = np.empty((Nobj,),  order='C', dtype=np.dtype('=f'))
    cdef int   ds=sizeof(float)



    for i from 0 <= i < Nobj:
        #mass
        fread(&mass_block[i], ds,1, buf)
        #pos
        fread(&pos_block[i,0], ds,1, buf)
        fread(&pos_block[i,1], ds,1, buf)
        fread(&pos_block[i,2], ds,1, buf)
        #id
        fread(&id_block[i], ds,1, buf)
        #pvel
        fread(&vel_block[i,0], ds,1, buf)
        fread(&vel_block[i,1], ds,1, buf)
        fread(&vel_block[i,2], ds,1, buf)

    return id_block, pos_block, vel_block, mass_block

cdef read_par_data(FILE *buf, int[:] ipar, float[:] rpar):

    cdef int di=sizeof(int)
    cdef int ds=sizeof(float)

    fread(&ipar[0], di,6, buf)
    fread(&rpar[0], ds,9, buf)

    return 0

cpdef read_file_c(str filename):

    cdef FILE *fo
    cdef float[:] rpar = np.empty((9,),  order='C',dtype=np.dtype('=f'))
    cdef int[:] ipar = np.empty((6,),  order='C',dtype=np.dtype('=i'))
    cdef int Ntot

    dic={}

    fo=fopen(filename.encode(), "rb")

    #header
    read_par_data(fo, ipar, rpar)
    dic['ng']=ipar[0]
    dic['nh']=ipar[1]
    dic['ns']=ipar[2]
    dic['comp']=ipar[3]
    dic['ngal']=ipar[4]
    dic['iout']=ipar[5]
    dic['mg']=rpar[0]
    dic['ms']=rpar[1]
    dic['mh']=rpar[2]
    dic['time']=rpar[3]
    dic['tdyn']=rpar[4]
    dic['e0']=rpar[5]
    dic['g0']=rpar[6]
    dic['b0']=rpar[7]
    dic['vir0']=rpar[8]
    #Particle
    Ntot=ipar[0]+ipar[1]+ipar[2]
    ida,posa,vela,massa=read_block_data(fo, Ntot)


    fclose(fo)

    return dic, ida, posa, vela,massa

cpdef read_header(str filename):

    cdef FILE *fo
    cdef float[:] rpar = np.empty((9,),  order='C',dtype=np.dtype('=f'))
    cdef int[:] ipar = np.empty((6,),  order='C',dtype=np.dtype('=i'))
    dic={}

    fo=fopen(filename.encode(), "rb")
    #header
    read_par_data(fo, ipar, rpar)
    dic['ng']=ipar[0]
    dic['nh']=ipar[1]
    dic['ns']=ipar[2]
    dic['comp']=ipar[3]
    dic['ngal']=ipar[4]
    dic['iout']=ipar[5]
    dic['mg']=rpar[0]
    dic['ms']=rpar[1]
    dic['mh']=rpar[2]
    dic['time']=rpar[3]
    dic['tdyn']=rpar[4]
    dic['e0']=rpar[5]
    dic['g0']=rpar[6]
    dic['b0']=rpar[7]
    dic['vir0']=rpar[8]

    return dic