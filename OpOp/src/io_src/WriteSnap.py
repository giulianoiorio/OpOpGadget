import struct


from ..utility_src import utility
from ..particle_src.particle import Particles


#ok funge
def write_snap(particles,filename,end='<',enable_mass=False,safe_write=True,verbose=True):
    """
    Write a Particles object in a file following the Gadget 1-type convention
    :param particles: Particles list object
    :param filename: Filename to write the Conditions
    :param end: type of binary writing- < littlendian, >bigendian, = native
    :param enable_mass:
    :param safe_write: If True, check the existence of the filename and eventually ask if
            overwrite or change name. If False, always overwrite.
    :return: Number of byte written
    """

    #Check input particle
    if isinstance(particles,Particles)==False: raise IOError('Incorrect particles or filename')


    #Check file
    if safe_write==True:
        s=utility.file_check(filename) #check if the output file already exist
        if s=='a': stream=open(filename,"ab")
        elif s=='': stream=open(filename,"wb")
        else: stream=open(s,"wb")
    else:
        stream=open(filename,"wb")

    #Check format
    wformat={'=':'Native', '<':'little-endian', '>':'big-endian'}
    if end not in wformat:
        print('Warning: incorrect write format type, it will be set to default value (little-endian)')
        end='<'

    cbyte=0 #Initialize byte counter

    #header
    print ("\nWriting header in file " + stream.name)

    block_check_start = struct.pack(end+"i", 256) #write the initial block check (C-type int)
    stream.write(block_check_start)
    cbyte+=4 #integer start block check

    #####
    #rattoppo da sistemare
    hpart=[0,0,0,0,0,0]

    for j in range(particles.header['NumFiles']):
        for i in range(6): hpart[i]=hpart[i]+particles.header['Npart'][j][i]

    for i in range(6):
        buf=struct.pack(end+"i",hpart[i])
        stream.write(buf)

    for i in range(6):
        buf=struct.pack(end+"d",particles.header['Massarr'][0][i])
        stream.write(buf)
    ###

    buf=struct.pack(end+"d",particles.header['Time'])
    stream.write(buf)

    buf=struct.pack(end+"d",particles.header['Redshift'])
    stream.write(buf)

    buf=struct.pack(end+"i",particles.header['FlagSfr'])
    stream.write(buf)

    buf=struct.pack(end+"i",particles.header['FlagFeedback'])
    stream.write(buf)

    for i in range(6):
        buf=struct.pack(end+"i",particles.header['Nall'][i])
        stream.write(buf)

    buf=struct.pack(end+"i",particles.header['FlagCooling'])
    stream.write(buf)


    #buf=struct.pack(end+"i",particles.header['NumFiles'])
    buf=struct.pack(end+"i",1) #per ora salviamo sempre e solo un file
    stream.write(buf)

    buf=struct.pack(end+"d",particles.header['BoxSize'])
    stream.write(buf)

    buf=struct.pack(end+"d",particles.header['Omega0'])
    stream.write(buf)

    buf=struct.pack(end+"d",particles.header['OmegaLambda'])
    stream.write(buf)

    buf=struct.pack(end+"d",particles.header['HubbleParam'])
    stream.write(buf)

    buf=struct.pack(end+"i",particles.header['FlagAge'])
    stream.write(buf)

    buf=struct.pack(end+"i",particles.header['FlagMetals'])
    stream.write(buf)

    for i in range(6):
        buf=struct.pack(end+"i",particles.header['NallHW'][i])
        stream.write(buf)

    buf=struct.pack(end+"i",particles.header['flag_entr_ics'])
    stream.write(buf)


    #Now write the last unused byte (the total dimension of the header is 256.)
    last_bys=256 - 24 - 48 - 8 - 8 - 4 - 4 - 24 - 4 - 4 - 8 - 8 - 8 - 8 - 4 - 4 - 24 - 4
    for i in range(last_bys):
        buf=struct.pack("<c",b'a')
        stream.write(buf)


    cbyte+=256#header fixed bytes

    block_check_end = struct.pack(end+"i", 256) #write the final block check (C-type int)
    stream.write(block_check_end)
    cbyte+=4 #integer end block check

    print ("header written in file %s." % (stream.name))

    #Particles

    print("\nWriting particle data in file " + stream.name)
    N=particles.n
    l=len(str(N))
    prog_in=str(0).zfill(l)+"/"+str(N)     #show the progress in a correct format

    block_byte=12*N #lenght of the Position and vel block: 3 float (12 by) per particle per total number of particles

    #Positions
    print("Writing Position Block.... Particle:",prog_in, end="\b"*len(prog_in))  #Initialize progres

    block_check_start = struct.pack(end+"i",block_byte)
    stream.write(block_check_start)

    cbyte+=4 #integer start block check

    for i in range(N):
        fmt=end+"fff"
        buf=struct.pack(fmt,particles.Pos[i][0], particles.Pos[i][1], particles.Pos[i][2])
        stream.write(buf)

        #output check
        prog=str(i).zfill(l)
        print(prog, end="\b"*l,flush=True)   #print progress

    print(str(N)+"/"+str(N),end="")        #print the completed progress

    cbyte+=block_byte #block bytes

    block_check_end = struct.pack(end+"i",block_byte)
    stream.write(block_check_end)

    cbyte+=4 #integer end block check

    print(".....Position Block Written","("+str(block_byte/1000),"KB)")

    #Velocities
    print("Writing Velocity Block.... Particle:",prog_in, end="\b"*len(prog_in))  #Initialize progres

    block_check_start = struct.pack(end+"i",block_byte)
    stream.write(block_check_start)
    cbyte+=4 #integer start block check

    for i in range(N):
        fmt=end+"fff"
        buf=struct.pack(fmt,particles.Vel[i][0], particles.Vel[i][1], particles.Vel[i][2])
        stream.write(buf)

        #output check
        prog=str(i).zfill(l)
        print(prog, end="\b"*l,flush=True)   #print progress
    print(str(N)+"/"+str(N),end="")        #print the completed progress

    cbyte+=block_byte #block bytes

    block_check_end = struct.pack(end+"i",block_byte)
    stream.write(block_check_end)
    cbyte+=4 #integer end block check

    print(".....Velocity Block Written","("+str(block_byte/1000),"KB)")

    #Id
    print("Writing Id Block.... Particle:",prog_in, end="\b"*len(prog_in))  #Initialize progres

    block_byte=4*N #lenght of the Id block: 1 integer (4 by) per particle per total number of particles

    block_check_start = struct.pack(end+"i",block_byte)
    stream.write(block_check_start)
    cbyte+=4 #integer start block check

    for i in range(N):
        buf=struct.pack(end+"i", particles.Id[i])
        stream.write(buf)

        #output check
        prog=str(i).zfill(l)
        print(prog, end="\b"*l,flush=True)   #print progress
    print(str(N)+"/"+str(N),end="")        #print the completed progress

    cbyte+=block_byte #block bytes

    block_check_end = struct.pack(end+"i",block_byte)
    stream.write(block_check_end)
    cbyte+=4 #integer end block check

    print(".....Id Block Written","("+str(block_byte/1000),"KB)")

    #Mass
    if enable_mass==True:

        print("Writing Mass Block.... Particle:",prog_in, end="\b"*len(prog_in))  #Initialize progres

        block_byte=4*N #lenght of the Mass block: 1 float (4 by) per particle per total number of particles

        block_check_start = struct.pack(end+"i",block_byte)
        stream.write(block_check_start)
        cbyte+=4 #integer start block check

        for i in range(N):
            buf=struct.pack(end+"f", particles.Mass[i])
            stream.write(buf)

            #output check
            prog=str(i).zfill(l)
            print(prog, end="\b"*l,flush=True)   #print progress
        print(str(N)+"/"+str(N),end="")        #print the completed progress

        cbyte+=block_byte #block bytes

        block_check_end = struct.pack(end+"i",block_byte)
        stream.write(block_check_end)

        cbyte+=4 #integer end block check

        print(".....Mass Block Written","("+str(block_byte/1000),"KB)")


    else:
        print("Mass block skipped")

    if particles.header['Nall'][0]>0:
        #U
        Ngas=particles.header['Nall'][0]

        l=len(str(Ngas))
        prog_in=str(0).zfill(l)+"/"+str(Ngas)     #show the progress in a correct format

        print("Writing U Block.... Particle:",prog_in, end="\b"*len(prog_in))  #Initialize progres


        block_byte=4*Ngas #lenght of the U block: 1 float (4 by) per particle per total number of gas particles

        block_check_start = struct.pack(end+"i",block_byte)
        stream.write(block_check_start)
        cbyte+=4 #integer start block check

        for i in range(Ngas):
            buf=struct.pack(end+"f",particles.U[i])
            stream.write(buf)

            #output check
            prog=str(i).zfill(l)
            print(prog, end="\b"*l,flush=True)   #print progress
        print(str(Ngas)+"/"+str(Ngas),end="")        #print the completed progress

        cbyte+=block_byte #block bytes

        block_check_end = struct.pack(end+"i",block_byte)
        stream.write(block_check_end)
        cbyte+=4 #integer end block check

        print(".....U Block Written","("+str(block_byte/1000),"KB)")

    else: print("No gas particles.... U block skipped")

    stream.close()

    print("particle written in file",stream.name,"with format",wformat[end],"("+str(cbyte/1000)+" KB)")

    return cbyte