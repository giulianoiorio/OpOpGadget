import struct
import os.path
import numpy as np

from .io_c_ext import cread_fvfps as cr
from ..particle_src.particle import Header,Particles
from ..utility_src import utility


#func
def _load_tstp(particles,fopen,nini, nfin, end='<'):
    """
    To load tstp
    fopen: open file
    :return:
    """
    f=fopen
    p=particles
    block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
    for ii in range(nini,nfin): p.Tstp[ii] = float(struct.unpack(end+"f", f.read(4))[0])  #read the 3-d pos and append it as a list of int in the Pos vector
    block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
    print('TSTP check',block_check_end,block_check_start)
    # block control check.
    if block_check_start != block_check_end:
        utility.Continue_check("Warning: Control TSTP Block failed")
        return 1
    else:
        return 0

def _load_pot(particles,fopen, nini, nfin, end='<'):
    """
    To load tstp
    fopen: open file
    :return:
    """
    f=fopen
    p=particles
    block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
    for ii in range(nini,nfin): p.Pot[ii] = float(struct.unpack(end+"f", f.read(4))[0])  #read the 3-d pos and append it as a list of int in the Pos vector
    block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
    print('Pot check',block_check_end,block_check_start)
    # block control check.
    if block_check_start != block_check_end:
        utility.Continue_check("Warning: Control Potential Block failed")
        return 1
    else:
        return 0

def _load_acce(particles, fopen, nini, nfin, end='<'):
        """
        To load acc
        fopen: open file
        :return:
        """
        f = fopen
        p = particles
        block_check_start = struct.unpack(end + "i", f.read(4))[0]  # read the initial block checek (C-type int)
        for ii in range(nini,nfin): p.Acce[ii] = struct.unpack(end + "fff",f.read(12))  # read the 3-d pos and append it as a list of float in the Pos vector
        block_check_end = struct.unpack(end + "i", f.read(4))[0]  # read the final block checek (C-type int)
        print('Acce check', block_check_end, block_check_start)
        # block control check.
        if block_check_start != block_check_end:
            utility.Continue_check("Warning: Control Acce Block failed")
            return 1

        return 0

def _load_pos(particles, fopen, nini, nfin, end='<'):
    """
    To load acc
    fopen: open file
    :return:
    """
    f = fopen
    p = particles

    #Load particle position (Pos)
    block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
    for ii in range(nini,nfin):
        p.Pos[ii] = struct.unpack(end+"fff", f.read(12))  #read the 3-d pos and append it as a list of float in the Pos vector
        p.Pcord[ii] = "Cartesian (X,Y,Z)" #Standard gadget coordinate system
    block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
    print('POS check',block_check_end,block_check_start)

    # block control check.
    if block_check_start != block_check_end:
        utility.Continue_check("Warning: Control Pos Block failed")
        return 1

    return 0

def _load_vel(particles, fopen, nini, nfin, end='<'):
    """
    To load acc
    fopen: open file
    :return:
    """
    f = fopen
    p = particles

    block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
    for ii in range(nini,nfin):
        p.Vel[ii] = struct.unpack(end+"fff", f.read(12))  #read the 3-d pos and append it as a list of float in the Pos vector
        p.Vcord[ii] = "Cartesian (Vx,Vy,Vz)" #Standard gadget coordinate system
    block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
    print('Vel check',block_check_end,block_check_start)

    # block control check.
    if block_check_start != block_check_end:
        utility.Continue_check("Warning: Control Vel Block failed")
        return 1

    return 0

def _load_id(particles, fopen, nini, nfin, end='<'):
    """
    To load acc
    fopen: open file
    :return:
    """
    f = fopen
    p = particles

    #Load particle Id (Id)
    block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
    for ii in range(nini,nfin): p.Id[ii] = int(struct.unpack(end+"i", f.read(4))[0])  #read the 3-d pos and append it as a list of int in the Pos vector
    block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
    print('Id check',block_check_end,block_check_start)

    # block control check.
    if block_check_start != block_check_end:
        utility.Continue_check("Warning: Control Id Block failed")
        return 1

    return 0

def _load_mass(particles, fopen, nini, end='<',findex=0):
    """
    To load acc
    fopen: open file
    :return:
    """
    f = fopen
    p = particles

    ntot_withmass = 0
    for k in range(6):
        if p.header['Massarr'][0][k] == 0: ntot_withmass = ntot_withmass + p.header["Npart"][findex][k]
    print("Total number of particles with mass block: %i" % (ntot_withmass))

    if ntot_withmass > 0:  # If exist at least one type of particle with mass block
        block_check_start = struct.unpack(end + "i", f.read(4))[0]  # read the initial block checek (C-type int)
        pcount = nini
        for k in range(6):  # repeat for each particle types
            if p.header["Massarr"][0][k] > 0:  # case 1- laod the mass from the header
                for n in range(pcount, p.header["Npart"][findex][k] + pcount):
                    p.Type[n] = k
                    p.Mass[n] = p.header["Massarr"][findex][k]
            else:  # case 2- laod the mass from the mass block
                for n in range(pcount, p.header["Npart"][findex][k] + pcount):
                    p.Type[n] = k
                    p.Mass[n] = float(struct.unpack(end + "f", f.read(4))[0])
            pcount = p.header["Npart"][findex][k] + nini
        block_check_end = struct.unpack(end + "i", f.read(4))[0]  # read the final block checek (C-type int)
        # block control check.
        if block_check_start != block_check_end:
            utility.Continue_check("Warning: Control Id Block failed")
            return 1
    else:  # case 1- laod the mass from the header
        pcount = nini
        for k in range(6):  # repeat for each particle types
            for n in range(pcount, p.header["Npart"][findex][k] + pcount):
                p.Type[n] = k
                p.Mass[n] = p.header["Massarr"][0][k]
            pcount = p.header["Npart"][findex][k] + nini
        print('Mass check')


    return 0


# def
dict_eblock = {'pot': _load_pot, 'acce': _load_acce, 'tstp': _load_tstp}
dict_eblock_ord = ('pot', 'acce', 'tstp')

#ok funge
def _load_single(filename,end='<',order_key='Id',verbose=False,extra_block=(),**kwargs):
    """
    Load data from snapshot written following the Gadget 1-type binary convention from a single file
    :param filename: Filename to read the data.
    :param end: type of binary writing- < littlendian, >bigendian, = native
    :param order_key: Key to order in ascending order the data. (It can be, Id, Radius, Vel_tot, Type, Mass, Energy, Potential, Acceleration)
                      If None, the data will not be ordered
    :return: Particles object with the data loaded by filename
    """

    #Check format
    wformat={'=':'Native', '<':'little-endian', '>':'big-endian'}
    if end not in wformat:
        print('Warning: incorrect  format type, it will be set to default value (little-endian)')
        end='<'
    f = open(filename, "rb") ##Open file in read binary mode ("rb")
    print ("\nReading header in file %s......" % (filename))

    h=Header()

    #Read header

    block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block check (C-type int)
    #Now read the binary file following the format=1 in Gadget.
    h.header['Npart'][0] = list( struct.unpack(end+"iiiiii", f.read(24)) )
    h.header['Massarr'][0] = list( struct.unpack(end+"dddddd", f.read(48)) )
    h.header['Time'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['Redshift'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['FlagSfr'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['FlagFeedback'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['Nall'] =  list(struct.unpack(end+"iiiiii", f.read(24)))
    h.header['FlagCooling'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['NumFiles'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['BoxSize'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['Omega0'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['OmegaLambda'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['HubbleParam'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['FlagAge'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['FlagMetals'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['NallHW'] =  list(struct.unpack(end+"iiiiii", f.read(24)))
    h.header['flag_entr_ics'] = struct.unpack(end+"i", f.read(4))[0]

    #Now read the last unused byte (the total dimension of the header is 256.)
    f.read(256 - 24 - 48 - 8 - 8 - 4 - 4 - 24 - 4 - 4 - 8 - 8 - 8 - 8 - 4 - 4 - 24-4 )
    block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the end block check (C-type int)
    print('header check',block_check_end,block_check_start)

    h.header['Ntot'] = sum(h.header['Npart'][0]) #Store the total number of particle in the header, we sum in Npart instead Nall
                                                           #because is safer. Indeed, in some Ics file the Nall file is left to 0.
    h.header['filename'] = filename #Store the name of the read file in the header
    print ("header data loaded from file %s \n" % (filename))


    #Create the particle object
    p=Particles(h=h)
    #Particle
    print ("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print ("Reading particle data from %s......" % (filename))
    warning_count=0 #Initialization of the counter for the warning in the file reading


    warning_count += _load_pos(particles=p, fopen=f,end=end,nini=0,nfin=p.n)
    warning_count += _load_vel(particles=p, fopen=f,end=end,nini=0,nfin=p.n)
    warning_count += _load_id(particles=p, fopen=f,end=end,nini=0,nfin=p.n)
    warning_count += _load_mass(particles=p, fopen=f, nini=0, end='<',findex=0)




    #Load particle Mass and Type
    #WARNING: The block Type does not exist in the Gadget output and initial condition, but it can be
    #found by the Id and the information of Npart in the array, in fact the first Npart[0] in the
    #file will be always gas particale the following Npart[1] will be always halo...ecc. The mass
    #is a bit more problematic: a block Mass for the particle k exists only if the related Massarr[k]
    #in the header is set to 0.

    """
    #check particle with mass: the following routine calculates the number of particles that have not
    # mass in the header, but have a mass block in the file.
    ntot_withmass = 0
    for k in range(6):
        if p.header['Massarr'][0][k] == 0: ntot_withmass = ntot_withmass + p.header["Npart"][0][k]

    if ntot_withmass > 0: #If exist at least one type of particle with mass block
        block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
        pcount=0
        for k in range(6): #repeat for each particle types
            if p.header["Massarr"][0][k] > 0: #case 1- laod the mass from the header
                for i in range(pcount, p.header["Npart"][0][k]+ pcount):
                    p.Type[i]=k
                    p.Mass[i]= p.header["Massarr"][0][k]
            else: #case 2- laod the mass from the mass block
                for i in range(pcount, p.header["Npart"][0][k]+ pcount):
                    p.Type[i]=k
                    p.Mass[i]= float( struct.unpack("<f", f.read(4))[0] )
            pcount= p.header["Npart"][0][k]
        block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
        print('Mass check',block_check_end,block_check_start)
        #block control check.
        if block_check_start != block_check_end:
            utility.Continue_check("Warning: Control Mass Block failed")
            warning_count+=1

    else: #case 1- laod the mass from the header
        pcount=0
        for k in range(6): #repeat for each particle types
            for i in range(pcount, p.header["Npart"][0][k]+ pcount):
                    p.Type[i]=k
                    p.Mass[i]= p.header["Massarr"][0][k]
            pcount= p.header["Npart"][0][k]
    """

    #Load SPH paramenters only for 0-Type particle
    if p.header["Npart"][0][0] > 0:
        #Load internal energy per unity mass U
        block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
        for i in range(0, p.header["Npart"][0][0]): p.U[i]=float((struct.unpack(end+"f", f.read(4))[0]))
        block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
        #block control check.
        if block_check_start != block_check_end:
            utility.Continue_check("Warning: Control U Block failed")
            warning_count+=1



    #Optional block
    for eblock in dict_eblock_ord:
        if eblock in extra_block:
            func_fill = dict_eblock[eblock]
            warning_count += func_fill(particles=p, fopen=f,end=end, nini=0, nfin=p.n)
        else:
            pass


    f.close()
    p.setrad()
    p.setvelt()
    p.set_pardic()


    if order_key is not None:
        print('Sorting by %s'% order_key)
        p.order(key=order_key)
        print('Sorted')

    return p

#ok funge
def _load_multi(filename,end='<',order_key='Id',verbose=False,extra_block=(),**kwargs):
    """
    Load data from snapshot written following the Gadget 1-type binary convention and from multiple subfile
    :param filename: Filename to read the data.  Write only
                     the name of the file without the .n
    :param end: type of binary writing- < littlendian, >bigendian, = native
    :param order_key: Key to order in ascending order the data. (It can be, Id, Radius, Vel_tot, Type, Mass, Energy, Potential, Acceleration)
                      If None, the data will not be ordered
    :return: Particles object with the data loaded by filename
    """

    #Check format
    wformat={'=':'Native', '<':'little-endian', '>':'big-endian'}
    if end not in wformat:
        print('Warning: incorrect write format type, it will be set to default value (little-endian)')
        end='<'
    buff = filename + ".0" #Name of the first file.
    f = open(buff, "rb")
    print ("\nReading header in file " + buff)

    h=Header()
    #Read header
    block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block check (C-type int)
    #Now read the binary file following the format=1 in Gadget.
    h.header['Npart'][0] = list( struct.unpack(end+"iiiiii", f.read(24)) )
    h.header['Massarr'][0] = list( struct.unpack(end+"dddddd", f.read(48)) )
    h.header['Time'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['Redshift'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['FlagSfr'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['FlagFeedback'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['Nall'] =  list(struct.unpack(end+"iiiiii", f.read(24)))
    h.header['FlagCooling'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['NumFiles'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['BoxSize'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['Omega0'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['OmegaLambda'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['HubbleParam'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['FlagAge'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['FlagMetals'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['NallHW'] =  list(struct.unpack(end+"iiiiii", f.read(24)))
    h.header['flag_entr_ics'] = struct.unpack(end+"i", f.read(4))[0]

    #Now read the last unused byte (the total dimension of the header is 256.)
    f.read(256 - 24 - 48 - 8 - 8 - 4 - 4 - 24 - 4 - 4 - 8 - 8 - 8 - 8 - 4 - 4 - 24-4 )
    block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the end block check (C-type int)
    #Now check if the initial and final block check is equal:
    if block_check_start != block_check_end: utility.Continue_check("Warning: Control Header Block failed at file" + "." + str(i))
    f.close()

    #Now add the information of Nall from the other files
    for i in range(1, h.header["NumFiles"]):        #Now add the information of Nall from the other files
        buff = filename + "." + str(i) #name of the current reading file
        f = open(buff, "rb")

        print ("Reading header in file " + buff)

        block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block check (C-type int)
        h.header["Npart"].append(list(struct.unpack(end+"iiiiii", f.read(24)))) #update the Npart information
        f.read(256 - 24) #Skip all the other byte that are equal among the files
        block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the final block check (C-type int)
        #Now check if the initial and final block check is equal:
        if block_check_start != block_check_end: utility.Continue_check("Warning: Control Header Block failed at file" + "." + str(i))
        f.close()

    for i in range(h.header["NumFiles"]): h.header['Ntot'] += sum(h.header['Npart'][i]) #Store the total number of particle in the header, we sum in Npart instead Nall
                                                                                                    #because is safer. Indeed, in some Ics file the Nall file is left to 0.
    h.header['filename'] = filename #Store the name of the read file in the header
    print ("header data loaded from file %s. \n" % (filename))
    warning_count_tot=0 #Initialization of the global counter for the warning for all the file

    p=Particles(h=h)

    print ("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print ("Starting reading particles data data from %i files (%s.n):" % (p.header["NumFiles"], p.header["filename"]))

    for i in range(p.header["NumFiles"]): #repeat for each of the n- files
        #open file filename.i in binary mode
        buff = filename + "." + str(i)
        f = open(buff, "rb")
        warning_count=0 #Initialization of the counter for the warning in the current file reading
        #skipping header and related control block (4 byte int + 256 byte header + 4 byte int)
        f.read(256 + 8)
        #Update the particle counter
        if i==0:
            start_count= 0
            end_count=sum(p.header["Npart"][0])
        else:
            start_count= end_count
            end_count+=sum(p.header["Npart"][i])

        if i==0: print ("\nReading file %s...... (particles from %i to %i)" % (buff,start_count,end_count))
        else: print ("Reading file %s...... (particles from %i to %i)" % (buff,start_count,end_count))

        warning_count += _load_pos(particles=p, fopen=f, end=end, nini=start_count, nfin=end_count)
        warning_count += _load_vel(particles=p, fopen=f, end=end, nini=start_count, nfin=end_count)
        warning_count += _load_id(particles=p, fopen=f, end=end, nini=start_count, nfin=end_count)
        warning_count += _load_mass(particles=p, fopen=f, nini=start_count, end='<', findex=i)

        """
        #Load particle position (Pos)
        block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
        for n in range(start_count, end_count):
            p.Pos[n] = struct.unpack(end+"fff", f.read(12))  #read the 3-d pos and append it as a list of float in the Pos vector
            p.Pcord[n] = "Cartesian (X,Y,Z)" #Standard gadget coordinate system
        block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
        #block control check.
        if block_check_start != block_check_end:
            utility.Continue_check("Warning: Control Position Block failed in file " + buff)
            warning_count+=1

        #Load particle velocity (Vel)
        block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
        for n in range(start_count, end_count):
            p.Vel[n] = struct.unpack(end+"fff", f.read(12))  #read the 3-d pos and append it as a list of float in the Pos vector
            p.Vcord[n] = "Cartesian (Vx,Vy,Vz)" #Standard gadget coordinate system
        block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
        #block control check.
        if block_check_start != block_check_end:
            utility.Continue_check("Warning: Control Velocity Block in file " + buff)
            warning_count+=1

        #Load particle Id (Id)
        block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
        for n in range(start_count, end_count): p.Id[n] = int(struct.unpack(end+"i", f.read(4))[0])  #read the 3-d pos and append it as a list of int in the Pos vector
        block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
        #block control check.
        if block_check_start != block_check_end:
            utility.Continue_check("Warning: Control Id Block in file " + buff)
            warning_count+=1


        #Load particle Mass and Type
        #WARNING: The block Type does not exist in the Gadget output and initial condition, but it can be
        #found by the Id and the information of Npart in the array, in fact the first Npart[0] in the
        #file will be always gas particale the following Npart[1] will be always halo...ecc. The mass
        #is a bit more problematic: a block Mass for the particle k exists only if the related Massarr[k]
        #in the header is set to 0.

        #check particle with mass: the following routine calculates the number of particles that have not
        #  mass in the header, but have a mass block in the file.
        ntot_withmass = 0
        for k in range(6):
            if p.header['Massarr'][0][k] == 0: ntot_withmass = ntot_withmass + p.header["Npart"][i][k]
        print ("Total number of particles with mass block: %i" % (ntot_withmass))

        if ntot_withmass > 0: #If exist at least one type of particle with mass block
            block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
            pcount=start_count
            for k in range(6): #repeat for each particle types
                if p.header["Massarr"][0][k] > 0: #case 1- laod the mass from the header
                    for n in range(pcount, p.header["Npart"][i][k]+ pcount):
                        p.Type[n]=k
                        p.Mass[n]= p.header["Massarr"][i][k]
                else: #case 2- laod the mass from the mass block
                    for n in range(pcount, p.header["Npart"][i][k]+ pcount):
                        p.Type[n]=k
                        p.Mass[n]= float( struct.unpack(end+"f", f.read(4))[0] )
                pcount= p.header["Npart"][i][k] + start_count
            block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
            #block control check.
            if block_check_start != block_check_end:
                utility.Continue_check("Warning: Control Mass in file " + buff)
                warning_count+=1
        else: #case 1- laod the mass from the header
            pcount=start_count
            for k in range(6): #repeat for each particle types
                for n in range(pcount, p.header["Npart"][i][k]+ pcount):
                        p.Type[n]=k
                        p.Mass[n]= p.header["Massarr"][0][k]
                pcount= p.header["Npart"][i][k] + start_count
        """

        #Load SPH paramenters only for 0-Type particle
        if p.header["Npart"][i][0] > 0:
            #Load internal energy per unity mass U
            block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
            for n in range(start_count, p.header["Npart"][i][0]+start_count):
                p.U[n]=float((struct.unpack(end+"f", f.read(4))[0]))
            block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
            #block control check.
            if block_check_start != block_check_end:
                utility.Continue_check("Warning: Control U Block in file " + buff)
                warning_count+=1

            #Load density RHO
            block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
            for n in range(start_count, h.header["Npart"][i][0]+start_count):
                p.Rho[n]=float((struct.unpack(end+"f", f.read(4))[0]))
            block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
            #block control check.
            if block_check_start != block_check_end:
                utility.Continue_check("Warning: Control Rho Block in file " + buff)
                warning_count+=1

        # Optional block
        for eblock in dict_eblock_ord:
            if eblock in extra_block:
                func_fill = dict_eblock[eblock]
                warning_count += func_fill(particles=p, fopen=f, end=end, nini=start_count, nfin=end_count)
            else:
                pass


        f.close()
        print ("Data load from "+buff+ " with " + str(warning_count) + " warnings")
        warning_count_tot+=warning_count

    print ("\nGlobal particles data load from " + filename +".n"  + " with " + str(warning_count_tot) + " warnings")
    print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")



    p.setrad()
    p.setvelt()
    p.set_pardic()


    if order_key is not None:
        print('Sorting by %s' % order_key)
        p.order(key=order_key)
        print('Sorted')

    return p


def _load_header_fvfps(dic,**kwargs):
    """
    Write header from rpar and ipar
    :param dic: dic with quantities from rpar and ipar, see cread_fvfps
    :param mscale: rescaling of the mass
    :param tscale: rescaline of the time
    :return:
    """
    if 'mscale' in kwargs: mscale=kwargs['mscale']
    else: mscale=1e10

    #Standard t scale in fvfps is lscale/vscale in s, but  in OpOp is in Gyr, so by default tscale=4.72e-3 yr
    if 'tscale' in kwargs: tscale=kwargs['tscale']
    else: tscale=4.718106379419217e-3 #yr


    if dic['ng']==0: massgas=0
    else: massgas=  ( dic['mg']/dic['ng'] )*mscale

    #print('ng',dic['ng'],'masstot',dic['mg'],'mpart',massgas)

    if dic['nh']==0: masshalo=0
    else: masshalo= ( dic['mh']/dic['nh'] ) *mscale

    #print('nh', dic['nh'], 'masstot', dic['mh'], 'mpart', masshalo)

    if dic['ns']==0: massstar=0
    else: massstar= ( dic['ms']/dic['ns'] ) *mscale

    #print('ns', dic['ns'], 'masstot', dic['ms'], 'mpart', massstar)

    h = Header()
    h.header['Npart'][0] = [dic['ng'],dic['nh'],dic['ns'],0,0,0]
    h.header['Massarr'][0] = [massgas,masshalo,massstar,0,0,0]
    h.header['Time'] = dic['time'] * tscale
    h.header['Tdyn'] = dic['tdyn'] * tscale
    h.header['Redshift'] = 0
    h.header['FlagSfr'] = 0
    h.header['FlagFeedback'] = 0
    h.header['Nall'] =  [dic['ng'],dic['nh'],dic['ns'],0,0,0]
    h.header['FlagCooling'] = 0
    h.header['NumFiles'] = 1
    h.header['BoxSize'] = 0.
    h.header['Omega0'] = 0.
    h.header['OmegaLambda'] = 0.
    h.header['HubbleParam'] = 0.
    h.header['FlagAge'] = 0
    h.header['FlagMetals'] = 0
    h.header['NallHW'] =  [0,0,0,0,0,0]
    h.header['flag_entr_ics'] = 0

    return h


def load_snap_fvfps(filename,order_key='Id',**kwargs):

    print ("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print ("Reading particle data from %s......" % (filename))
    dic, id, pos, vel, mass=cr.read_file_c(filename)
    print("Done")
    print ("\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

    #Standard length scale in fvfps is kpc as the one of OpOp, so by default lscale=1
    if 'lscale' in kwargs: lscale=kwargs['lscale']
    else: lscale=1

    #Standard mass scale in fvfps is 10^10 msun, but  in OpOp is 1msun, so by default mscale=1e10
    if 'mscale' in kwargs: mscale=kwargs['mscale']
    else: mscale=1e10

    #Standard mass scale in fvfps is Sqrt(G*mscale/lscale) in m/s, but  in OpOp is km/s, so by default vscale=207.38 km/s
    if 'vscale' in kwargs: vscale=kwargs['vscale']
    else: vscale=207.3844607923656


    #Standard t scale in fvfps is lscale/vscale in s, but  in OpOp is in Gyr, so by default tscale=4.72e-3 yr
    if 'tscale' in kwargs: tscale=kwargs['tscale']
    else: tscale=4.718106379419217e-3 #yr



    #part header
    h=_load_header_fvfps(dic,mscale=mscale,tscale=tscale)
    h.header['filename'] = filename


    #Create the particle object
    p=Particles(h=h)
    #Particle
    p.Pos=np.array(pos)*lscale
    p.Vel = np.array(vel)*vscale
    p.Id = np.array(id, dtype=int)
    p.Mass = np.array(mass)*mscale

    print('Sorting by Id')
    p.order(key='Id')
    print('Sorted')
    p._maketype()
    p.setrad()
    p.setvelt()
    p.set_pardic()

    if (order_key is not None) and (order_key!='Id'):
        print('Sorting by %s'% order_key)
        p.order(key=order_key)
        print('Sorted')




    return p

#ok funge
def load_snap(filename,end='<',order_key='Id',extra_block=(),kind='fvfps',**kwargs):
    """
    Load data from snapshot written following the Gadget 1-type binary convention
    :param filename: Filename to read the data. If reading multiple file, write only
                     the name of the file without the .n
    :param end: type of binary writing- < littlendian, >bigendian, = native
    :param order_key: Key to order in ascending order the data. (It can be, Id, Radius, Vel_tot, Type, Mass, Energy, Potential, Acceleration)
                      If None, the data will not be ordered
    :return: Particles object with the data loaded by filename
    """

    if kind.lower()=='gadget':
        if (os.path.isfile(filename)): particles=_load_single(filename,end=end,order_key=order_key,extra_block=extra_block,**kwargs)
        elif (os.path.isfile(filename+'.0')): particles=_load_multi(filename,end=end,order_key=order_key,extra_block=extra_block,**kwargs)
        else: raise IOError('File %s not found'%filename)
    elif kind.lower()=='fvfps':
        particles=load_snap_fvfps(filename,order_key=order_key,**kwargs)
    else:
        raise NotImplementedError('load file from %s not implemented' % kind)


    return particles

def _load_header_single(filename,end='<',**kwargs):
    """
    Load header from snapshot written following the Gadget 1-type binary convention from single file
    :param filename: Filename to read the data.
    :param end: type of binary writing- < littlendian, >bigendian, = native
    :return: Header object with the header loaded by filename
    """
    #Check format
    wformat={'=':'Native', '<':'little-endian', '>':'big-endian'}
    if end not in wformat:
        print('Warning: incorrect  format type, it will be set to default value (little-endian)')
        end='<'
    f = open(filename, "rb") ##Open file in read binary mode ("rb")
    print ("\nReading header in file %s......" % (filename))

    h=Header()

    #Read header

    block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block check (C-type int)
    #Now read the binary file following the format=1 in Gadget.
    h.header['Npart'][0] = list( struct.unpack(end+"iiiiii", f.read(24)) )
    h.header['Massarr'][0] = list( struct.unpack(end+"dddddd", f.read(48)) )
    h.header['Time'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['Redshift'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['FlagSfr'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['FlagFeedback'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['Nall'] =  list(struct.unpack(end+"iiiiii", f.read(24)))
    h.header['FlagCooling'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['NumFiles'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['BoxSize'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['Omega0'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['OmegaLambda'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['HubbleParam'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['FlagAge'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['FlagMetals'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['NallHW'] =  list(struct.unpack(end+"iiiiii", f.read(24)))
    h.header['flag_entr_ics'] = struct.unpack(end+"i", f.read(4))[0]

    #Now read the last unused byte (the total dimension of the header is 256.)
    f.read(256 - 24 - 48 - 8 - 8 - 4 - 4 - 24 - 4 - 4 - 8 - 8 - 8 - 8 - 4 - 4 - 24-4 )
    block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the end block check (C-type int)
    print('header check',block_check_end,block_check_start)

    h.header['Ntot'] = sum(h.header['Npart'][0]) #Store the total number of particle in the header, we sum in Npart instead Nall
                                                           #because is safer. Indeed, in some Ics file the Nall file is left to 0.
    h.header['filename'] = filename #Store the name of the read file in the header
    print ("header data loaded from file %s \n" % (filename))
    f.close()

    return h

def _load_header_multi(filename,end='<',**kwargs):
    """
    Load header from snapshot written following the Gadget 1-type binary convention from multiple subfile
    :param filename: Filename to read the data. Write only
                     the name of the file without the .n
    :param end: type of binary writing- < littlendian, >bigendian, = native
    :return: Header object with the header loaded by filename
    """

    #Check format
    wformat={'=':'Native', '<':'little-endian', '>':'big-endian'}
    if end not in wformat:
        print('Warning: incorrect write format type, it will be set to default value (little-endian)')
        end='<'
    buff = filename + ".0" #Name of the first file.
    f = open(buff, "rb")
    print ("\nReading header in file " + buff)

    h=Header()
    #Read header
    block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block check (C-type int)
    #Now read the binary file following the format=1 in Gadget.
    h.header['Npart'][0] = list( struct.unpack(end+"iiiiii", f.read(24)) )
    h.header['Massarr'][0] = list( struct.unpack(end+"dddddd", f.read(48)) )
    h.header['Time'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['Redshift'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['FlagSfr'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['FlagFeedback'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['Nall'] =  list(struct.unpack(end+"iiiiii", f.read(24)))
    h.header['FlagCooling'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['NumFiles'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['BoxSize'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['Omega0'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['OmegaLambda'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['HubbleParam'] = struct.unpack(end+"d", f.read(8))[0]
    h.header['FlagAge'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['FlagMetals'] = struct.unpack(end+"i", f.read(4))[0]
    h.header['NallHW'] =  list(struct.unpack(end+"iiiiii", f.read(24)))
    h.header['flag_entr_ics'] = struct.unpack(end+"i", f.read(4))[0]

    #Now read the last unused byte (the total dimension of the header is 256.)
    f.read(256 - 24 - 48 - 8 - 8 - 4 - 4 - 24 - 4 - 4 - 8 - 8 - 8 - 8 - 4 - 4 - 24-4 )
    block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the end block check (C-type int)
    #Now check if the initial and final block check is equal:
    if block_check_start != block_check_end: utility.Continue_check("Warning: Control Header Block failed at file" + "." + str(i))
    f.close()

    #Now add the information of Nall from the other files
    for i in range(1, h.header["NumFiles"]):        #Now add the information of Nall from the other files
        buff = filename + "." + str(i) #name of the current reading file
        f = open(buff, "rb")

        print ("Reading header in file " + buff)

        block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block check (C-type int)
        h.header["Npart"].append(list(struct.unpack(end+"iiiiii", f.read(24)))) #update the Npart information
        f.read(256 - 24) #Skip all the other byte that are equal among the files
        block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the final block check (C-type int)
        #Now check if the initial and final block check is equal:
        if block_check_start != block_check_end: utility.Continue_check("Warning: Control Header Block failed at file" + "." + str(i))
        f.close()

    for i in range(h.header["NumFiles"]): h.header['Ntot'] += sum(h.header['Npart'][i]) #Store the total number of particle in the header, we sum in Npart instead Nall
                                                                                                    #because is safer. Indeed, in some Ics file the Nall file is left to 0.
    h.header['filename'] = filename #Store the name of the read file in the header
    print ("header data loaded from file %s. \n" % (filename))
    warning_count_tot=0 #Initialization of the global counter for the warning for all the file
    f.close()

    return h

def load_header(filename,end='<', kind='fvfps',**kwargs):
    """
    Load header from snapshot written following the Gadget 1-type binary convention
    :param filename: Filename to read the data. If reading multiple file, write only
                     the name of the file without the .n
    :param end: type of binary writing- < littlendian, >bigendian, = native
    :return: Header object with the header loaded by filename
    """
    if kind.lower() == 'gadget':

        if (os.path.isfile(filename)): header_obj=_load_header_single(filename,end=end,**kwargs)
        elif (os.path.isfile(filename+'.0')): header_obj=_load_header_multi(filename,end=end,**kwargs)
        else: raise IOError('File %s not found'%filename)

    elif kind.lower() == 'fvfps':

        dic=cr.read_header(filename)
        header_obj = _load_header_fvfps(dic,**kwargs)
        header_obj.header['filename']=filename
    else:
        raise NotImplementedError('load file from %s not implemented' % kind)

    return header_obj



