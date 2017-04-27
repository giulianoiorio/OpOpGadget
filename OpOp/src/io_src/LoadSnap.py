import struct
import os.path

from ..particle_src.particle import Header,Particles
from ..utility_src import utility


#ok funge
def _load_single(filename,end='<',order_key='Id',verbose=False):
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

    #Load particle position (Pos)
    block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
    for i in range(0,p.n):
        p.Pos[i] = struct.unpack(end+"fff", f.read(12))  #read the 3-d pos and append it as a list of float in the Pos vector
        p.Pcord[i] = "Cartesian (X,Y,Z)" #Standard gadget coordinate system
    block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
    print('POS check',block_check_end,block_check_start)
    #block control check.
    if block_check_start != block_check_end:
        utility.Continue_check("Warning: Control Position Block failed")
        warning_count+=1


    #Load particle velocity (Vel)
    block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
    for i in range(0,p.n):
        p.Vel[i] = struct.unpack(end+"fff", f.read(12))  #read the 3-d pos and append it as a list of float in the Pos vector
        p.Vcord[i] = "Cartesian (Vx,Vy,Vz)" #Standard gadget coordinate system
    block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
    print('Vel check',block_check_end,block_check_start)
    #block control check.
    if block_check_start != block_check_end:
        utility.Continue_check("Warning: Control Velocity Block failed")
        warning_count+=1

    #Load particle Id (Id)
    block_check_start =  struct.unpack(end+"i", f.read(4))[0] #read the initial block checek (C-type int)
    for i in range(0,p.n): p.Id[i] = int(struct.unpack(end+"i", f.read(4))[0])  #read the 3-d pos and append it as a list of int in the Pos vector
    block_check_end =  struct.unpack(end+"i", f.read(4))[0] #read the final block checek (C-type int)
    print('Id check',block_check_end,block_check_start)
    # block control check.
    if block_check_start != block_check_end:
        utility.Continue_check("Warning: Control Id Block failed")
        warning_count+=1


    #Load particle Mass and Type
    #WARNING: The block Type does not exist in the Gadget output and initial condition, but it can be
    #found by the Id and the information of Npart in the array, in fact the first Npart[0] in the
    #file will be always gas particale the following Npart[1] will be always halo...ecc. The mass
    #is a bit more problematic: a block Mass for the particle k exists only if the related Massarr[k]
    #in the header is set to 0.


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
    f.close()

    p.setrad()
    p.setvelt()

    if order_key is not None:
        print('Sorting by %s'% order_key)
        p.order(key=order_key)
        print('Sorted')

    return p

#ok funge
def _load_multi(filename,end='<',order_key='Id'):
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

        f.close()
        print ("Data load from "+buff+ " with " + str(warning_count) + " warnings")
        warning_count_tot+=warning_count

    print ("\nGlobal particles data load from " + filename +".n"  + " with " + str(warning_count_tot) + " warnings")
    print ("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n")

    p.setrad()
    p.setvelt()

    if order_key is not None:
        print('Sorting by %s' % order_key)
        p.order(key=order_key)
        print('Sorted')

    return p

#ok funge
def load_snap(filename,end='<',order_key='Id'):
    """
    Load data from snapshot written following the Gadget 1-type binary convention
    :param filename: Filename to read the data. If reading multiple file, write only
                     the name of the file without the .n
    :param end: type of binary writing- < littlendian, >bigendian, = native
    :param order_key: Key to order in ascending order the data. (It can be, Id, Radius, Vel_tot, Type, Mass, Energy, Potential, Acceleration)
                      If None, the data will not be ordered
    :return: Particles object with the data loaded by filename
    """

    if (os.path.isfile(filename)): particles=_load_single(filename,end=end,order_key=order_key)
    elif (os.path.isfile(filename+'.0')): particles=_load_multi(filename,end=end,order_key=order_key)
    else: raise IOError('File %s not found'%filename)

    return particles

def _load_header_single(filename,end='<'):
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

def _load_header_multi(filename,end='<'):
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

def load_header(filename,end='<'):
    """
    Load header from snapshot written following the Gadget 1-type binary convention
    :param filename: Filename to read the data. If reading multiple file, write only
                     the name of the file without the .n
    :param end: type of binary writing- < littlendian, >bigendian, = native
    :return: Header object with the header loaded by filename
    """

    if (os.path.isfile(filename)): header_obj=_load_header_single(filename,end=end)
    elif (os.path.isfile(filename+'.0')): header_obj=_load_header_multi(filename,end=end)
    else: raise IOError('File %s not found'%filename)

    return header_obj



