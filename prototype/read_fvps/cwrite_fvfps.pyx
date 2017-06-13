from libc.stdio cimport FILE, fopen, fwrite,  fclose


cdef write_block_data(FILE *buf, int Nobj, float[:, ::1] pos_block, float[:, ::1] vel_block, float[::1] id_block, float[::1] mass_block):

    cdef int   ds=sizeof(float)



    for i from 0 <= i < Nobj:
        #mass
        fwrite(&mass_block[i], ds,1, buf)
        #pos
        fwrite(&pos_block[i,0], ds,1, buf)
        fwrite(&pos_block[i,1], ds,1, buf)
        fwrite(&pos_block[i,2], ds,1, buf)
        #id
        fwrite(&id_block[i], ds,1, buf)
        #pvel
        fwrite(&vel_block[i,0], ds,1, buf)
        fwrite(&vel_block[i,1], ds,1, buf)
        fwrite(&vel_block[i,2], ds,1, buf)

    return 0

cdef write_par_data(FILE *buf, int[:] ipar, float[:] rpar):

    cdef int di=sizeof(int)
    cdef int ds=sizeof(float)

    fwrite(&ipar[0], di,6, buf)
    fwrite(&rpar[0], ds,9, buf)

    return 0

cpdef write_file_c(str filename, int[:] ipar, float[:] rpar, float[:, ::1] pos_block, float[:, ::1] vel_block, float[::1] id_block, float[::1] mass_block):

    cdef FILE *fo
    cdef int Ntot

    dic={}

    fo=fopen(filename.encode(), "wb")

    Ntot=ipar[0]+ipar[1]+ipar[2]

    if Ntot!=len(id_block): raise ValueError('Ntot in header and in the array do not match')


    #header
    write_par_data(fo, ipar, rpar)
    write_block_data(fo, Ntot, pos_block, vel_block, id_block, mass_block)

    fclose(fo)

    return 0

