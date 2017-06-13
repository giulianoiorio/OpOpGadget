from libc.stdio cimport FILE, fopen, fwrite,  fclose


cdef write_block_data(FILE *buf, int Nobj, float[:, ::1] pos_block, float[:, ::1] vel_block, float[::1] id_block, float[::1] mass_block):
    #In fvfps both each particle block has 32 byte
    cdef int   ds=sizeof(float)
    cdef int di=sizeof(int)
    cdef int check=32


    for i from 0 <= i < Nobj:
        #Start check_block
        fwrite(&check, di,1, buf)

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

        #End check_block
        fwrite(&check, di,1, buf)

    return 0

cdef write_par_data(FILE *buf, int[:] ipar, float[:] rpar):
    #In fvfps both ipar and rpar are block of 40 byte
    #ipar deve essere lungo 10 con
    #0-ngas 1-nhalo 2-nstar 3-ncomp 4-ngal 6-iout 7-0 8-0 9-0
    #rpar deve essere lungo 10 con
    #0-mgas 1-mhalo 2-mstar 3-time 4-tdyn 5-emu 6-gam 7-beta 8-vir0 9-0


    cdef int di=sizeof(int)
    cdef int ds=sizeof(float)
    cdef int check=40

    #Start check_block
    fwrite(&check, di,1, buf)
    fwrite(&ipar[0], di,10, buf)
    #End check_block
    fwrite(&check, di,1, buf)

    #Start check_block
    fwrite(&check, di,1, buf)
    fwrite(&rpar[0], ds,10, buf)
    #End check_block
    fwrite(&check, di,1, buf)

    return 0

cpdef write_file_c(str filename, int[:] ipar, float[:] rpar, float[:, ::1] pos_block, float[:, ::1] vel_block, float[::1] id_block, float[::1] mass_block):
    #ipar deve essere lungo 10 con
    #0-ngas 1-nhalo 2-nstar 3-ncomp 4-ngal 6-iout 7-0 8-0 9-0
    #rpar deve essere lungo 10 con
    #0-mgas 1-mhalo 2-mstar 3-time 4-tdyn 5-emu 6-gam 7-beta 8-vir0 9-0
    #the other block need to be in float (4 byte)

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

