import numpy as np
cimport numpy as np
cimport cython
from libc.stdio cimport FILE, fopen, fread, stdout, fclose, fprintf, stderr
from libc.stdlib cimport exit, EXIT_FAILURE



cdef read_block_data(FILE *buf, int Nobj):

    cdef float[:, ::1] pos_block = np.empty((Nobj,3),  order='C', dtype=np.dtype('=f'))
    cdef float[:, ::1] vel_block = np.empty((Nobj,3),  order='C', dtype=np.dtype('=f'))
    cdef float[::1] id_block = np.empty((Nobj,),  order='C', dtype=np.dtype('=f'))
    cdef float[::1] mass_block = np.empty((Nobj,),  order='C', dtype=np.dtype('=f'))
    cdef int   ds=sizeof(float)
    cdef int di=sizeof(int)
    cdef int check, check2
    cdef int block_l


    for i from 0 <= i < Nobj:
        #check block (32byte)
        fread(&check,di,1,buf)

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

        #check block
        fread(&check2,di,1,buf)

        #control
        if check!=check2:
            errmess='\nError: Control blocks do is not match (particle)... exiting\n'
            fprintf(stderr, errmess.encode())
            exit(EXIT_FAILURE)

    return id_block, pos_block, vel_block, mass_block

cdef read_par_data(FILE *buf, int[:] ipar, float[:] rpar):

    cdef int di=sizeof(int)
    cdef int ds=sizeof(float)
    cdef int check, check2
    cdef int block_l

    #fortran unformatted have a 4-byte int at begginig and at the end of each written block
    #containg the bit of the block

    #ipar block
    fread(&check,di,1,buf)
    block_l=check/di
    fread(&ipar[0], di,block_l, buf)
    fread(&check2,di,1,buf)

    if check!=check2:
        errmess='\nError: Control blocks do is not match (ipar)... exiting\n'
        fprintf(stderr, errmess.encode())
        exit(EXIT_FAILURE)

    #rpar block
    fread(&check,di,1,buf)
    block_l=check/di
    fread(&rpar[0], ds,block_l, buf)
    fread(&check2,di,1,buf)

    if check!=check2:
        errmess='\nError: Control blocks do is not match (rpar)... exiting\n'
        fprintf(stderr, errmess.encode())
        exit(EXIT_FAILURE)

    return 0

cpdef read_file_c(str filename):

    cdef FILE *fo
    cdef float[:] rpar = np.empty((10,),  order='C',dtype=np.dtype('=f'))
    cdef int[:] ipar = np.empty((10,),  order='C',dtype=np.dtype('=i'))
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
    dic['mh']=rpar[1]
    dic['ms']=rpar[2]
    dic['time']=rpar[3]
    dic['tdyn']=rpar[4]
    dic['e0']=rpar[5]
    dic['g0']=rpar[6]
    dic['b0']=rpar[7]
    dic['vir0']=rpar[8]
    #Particle
    Ntot=dic['ng']+dic['nh']+dic['ns']
    ida,posa,vela,massa=read_block_data(fo, Ntot)


    fclose(fo)

    return dic, ida, posa, vela,massa

cpdef read_header(str filename):

    cdef FILE *fo
    cdef float[:] rpar = np.empty((10,),  order='C',dtype=np.dtype('=f'))
    cdef int[:] ipar = np.empty((10,),  order='C',dtype=np.dtype('=i'))
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
    dic['mh']=rpar[1]
    dic['ms']=rpar[2]
    dic['time']=rpar[3]
    dic['tdyn']=rpar[4]
    dic['e0']=rpar[5]
    dic['g0']=rpar[6]
    dic['b0']=rpar[7]
    dic['vir0']=rpar[8]

    return dic