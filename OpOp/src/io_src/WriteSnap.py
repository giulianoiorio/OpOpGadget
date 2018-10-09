from __future__ import  division, print_function
import struct
import numpy as np
from .io_c_ext  import cwrite_fvfps as cw

from ..utility_src import utility
from ..particle_src.particle import Particles


#ok funge
def write_snap_gadget(particles,filename,end='<',enable_mass=False,safe_write=True,verbose=True):
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
    if verbose: print ("\nWriting header in file " + stream.name)

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

    if verbose: print ("header written in file %s." % (stream.name))

    #Particles

    print("\nWriting particle data in file " + stream.name)
    N=particles.n
    l=len(str(N))
    prog_in=str(0).zfill(l)+"/"+str(N)     #show the progress in a correct format

    block_byte=12*N #lenght of the Position and vel block: 3 float (12 by) per particle per total number of particles

    #Positions
    if verbose: print("Writing Position Block.... Particle:",prog_in, end="\b"*len(prog_in))  #Initialize progres

    block_check_start = struct.pack(end+"i",block_byte)
    stream.write(block_check_start)

    cbyte+=4 #integer start block check

    for i in range(N):
        fmt=end+"fff"
        buf=struct.pack(fmt,particles.Pos[i][0], particles.Pos[i][1], particles.Pos[i][2])
        stream.write(buf)

        #output check
        prog=str(i).zfill(l)
        if verbose: print(prog, end="\b"*l,flush=True)   #print progress

    if verbose: print(str(N)+"/"+str(N),end="")        #print the completed progress

    cbyte+=block_byte #block bytes

    block_check_end = struct.pack(end+"i",block_byte)
    stream.write(block_check_end)

    cbyte+=4 #integer end block check

    if verbose: print(".....Position Block Written","("+str(block_byte/1000),"KB)")

    #Velocities
    if verbose: print("Writing Velocity Block.... Particle:",prog_in, end="\b"*len(prog_in))  #Initialize progres

    block_check_start = struct.pack(end+"i",block_byte)
    stream.write(block_check_start)
    cbyte+=4 #integer start block check

    for i in range(N):
        fmt=end+"fff"
        buf=struct.pack(fmt,particles.Vel[i][0], particles.Vel[i][1], particles.Vel[i][2])
        stream.write(buf)

        #output check
        prog=str(i).zfill(l)
        if verbose: print(prog, end="\b"*l,flush=True)   #print progress
    if verbose: print(str(N)+"/"+str(N),end="")        #print the completed progress

    cbyte+=block_byte #block bytes

    block_check_end = struct.pack(end+"i",block_byte)
    stream.write(block_check_end)
    cbyte+=4 #integer end block check

    if verbose: print(".....Velocity Block Written","("+str(block_byte/1000),"KB)")

    #Id
    if verbose: print("Writing Id Block.... Particle:",prog_in, end="\b"*len(prog_in))  #Initialize progres

    block_byte=4*N #lenght of the Id block: 1 integer (4 by) per particle per total number of particles

    block_check_start = struct.pack(end+"i",block_byte)
    stream.write(block_check_start)
    cbyte+=4 #integer start block check

    for i in range(N):
        buf=struct.pack(end+"i", particles.Id[i])
        stream.write(buf)

        #output check
        prog=str(i).zfill(l)
        if verbose: print(prog, end="\b"*l,flush=True)   #print progress
    if verbose: print(str(N)+"/"+str(N),end="")        #print the completed progress

    cbyte+=block_byte #block bytes

    block_check_end = struct.pack(end+"i",block_byte)
    stream.write(block_check_end)
    cbyte+=4 #integer end block check

    if verbose: print(".....Id Block Written","("+str(block_byte/1000),"KB)")

    #Mass
    if enable_mass==True:

        if verbose: print("Writing Mass Block.... Particle:",prog_in, end="\b"*len(prog_in))  #Initialize progres

        block_byte=4*N #lenght of the Mass block: 1 float (4 by) per particle per total number of particles

        block_check_start = struct.pack(end+"i",block_byte)
        stream.write(block_check_start)
        cbyte+=4 #integer start block check

        for i in range(N):
            buf=struct.pack(end+"f", particles.Mass[i])
            stream.write(buf)

            #output check
            prog=str(i).zfill(l)
            if verbose: print(prog, end="\b"*l,flush=True)   #print progress
        if verbose: print(str(N)+"/"+str(N),end="")        #print the completed progress

        cbyte+=block_byte #block bytes

        block_check_end = struct.pack(end+"i",block_byte)
        stream.write(block_check_end)

        cbyte+=4 #integer end block check

        if verbose: print(".....Mass Block Written","("+str(block_byte/1000),"KB)")


    else:
        if verbose: print("Mass block skipped")

    if particles.header['Nall'][0]>0:
        #U
        Ngas=particles.header['Nall'][0]

        l=len(str(Ngas))
        prog_in=str(0).zfill(l)+"/"+str(Ngas)     #show the progress in a correct format

        if verbose: print("Writing U Block.... Particle:",prog_in, end="\b"*len(prog_in))  #Initialize progres


        block_byte=4*Ngas #lenght of the U block: 1 float (4 by) per particle per total number of gas particles

        block_check_start = struct.pack(end+"i",block_byte)
        stream.write(block_check_start)
        cbyte+=4 #integer start block check

        for i in range(Ngas):
            buf=struct.pack(end+"f",particles.U[i])
            stream.write(buf)

            #output check
            prog=str(i).zfill(l)
            if verbose: print(prog, end="\b"*l,flush=True)   #print progress
        if verbose: print(str(Ngas)+"/"+str(Ngas),end="")        #print the completed progress

        cbyte+=block_byte #block bytes

        block_check_end = struct.pack(end+"i",block_byte)
        stream.write(block_check_end)
        cbyte+=4 #integer end block check

        if verbose: print(".....U Block Written","("+str(block_byte/1000),"KB)")

    elif verbose:
        print("No gas particles.... U block skipped")

    else:
        pass

    stream.close()

    print("particle written in file",stream.name,"with format",wformat[end],"("+str(cbyte/1000)+" KB)")

    return cbyte

def write_snap_fvfps(particles, filename, **kwargs):

    #tdyn in Gy, if not in kwargs tdyn=1 Myr

    #Scale
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

    #PAR
    if 'end' in kwargs: end = kwargs['end']
    else: end = '<'

    if 'tdyn' in kwargs: tdyn=kwargs['tdyn']/tscale
    else: tdyn=0.001/tscale #1 Myr


    pos=particles.Pos
    vel=particles.Vel
    id=particles.Id
    mass=particles.Mass
    leng=len(id)

    h=particles.header

    comp_count=0
    for item in h['Nall']:
        if item!=0: comp_count+=1

    #par array
    ipar=np.zeros((10,),  order='C', dtype=np.dtype(end+'i'))
    ipar[0]=h['Nall'][0]
    ipar[1]=h['Nall'][1]
    ipar[2]=h['Nall'][2]
    ipar[3]=comp_count
    ipar[4]=1


    rpar = np.zeros((10,), order='C', dtype=np.dtype(end + 'f'))
    rpar[0]=(h['Nall'][0]*h['Massarr'][0][0])/mscale
    rpar[1]=(h['Nall'][1]*h['Massarr'][0][1])/mscale
    rpar[2]=(h['Nall'][2]*h['Massarr'][0][2])/mscale



    #particle array
    posa=np.empty_like(pos,order='C', dtype=np.dtype(end + 'f'))
    if 'poscm' in kwargs:  pos[:,:]=pos[:,:]+np.array(pcom)
    posa[:,:]=pos[:,:]/lscale
    vela=np.empty_like(vel,order='C', dtype=np.dtype(end + 'f'))
    if 'velcm' in kwargs:  vel[:,:]=vel[:,:]+np.array(vcom)
    vela[:,:]=vel[:,:]/vscale
    massa=np.empty_like(mass,order='C', dtype=np.dtype(end + 'f'))
    massa[:]=mass[:]/mscale
    ida=np.empty_like(id,order='C', dtype=np.dtype(end + 'f'))
    ida[:]=id[:]

    cw.write_file_c(filename, ipar,  rpar, posa, vela, ida, massa)

def write_snap(particles, filename, kind='fvfps',**kwargs):


    if kind.lower()=='gadget':

        #param
        if 'end' in kwargs: end=kwargs['end']
        else: end='<'
        if 'enable_mass' in kwargs: enable_mass=kwargs['enable_mass']
        else: enable_mass=False
        if 'safe_write' in kwargs: safe_write=kwargs['safe_write']
        else: safe_write=True
        if 'verbose' in kwargs: verbose=kwargs['verbose']
        else: verbose=True

        write_snap_gadget(particles, filename, end=end, enable_mass=enable_mass, safe_write=safe_write, verbose=verbose)

    elif kind.lower()=='fvfps':

        write_snap_fvfps(particles, filename, **kwargs)

    else:

        raise NotImplementedError('kind %s not implemented in write snap'%kind)