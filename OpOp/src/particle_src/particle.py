__author__ = 'Giuliano'
#__all__=["Header","Particle","Particles"]

import numpy as np
import math as mt
from math import *
from ..utility_src.utility import find_symbol

class Header:
    '''
    The following object has a only single variable: a dictionary that stores the information of the header block.
    '''
    def __init__(self):
        self.header = {
            'Npart': [[0,0,0,0,0,0]], #Number of particles of each type in the current file.
            'Massarr': [[0,0,0,0,0,0]], #Mass of each type of particle. !!If 0 the mass must be read in the related block of particle_parameters.
            'Time':0, #Time of output.
            'Redshift': 0, #!!only set in cosmology.
            'FlagSfr': 0, #Flag for star formation !!actually unused
            'FlagFeedback': 0, #Flag for feedback !!actually unused
            'Nall': [0,0,0,0,0,0],  #Total number of particles of each type in the simulation. !!They can be divided in more than one file.
            'FlagCooling': 0, #Flag for cooling !!actually unused.
            'NumFiles': 1, #Number of files for each snapshot.
            'BoxSize': 0, #Box size !!used only in periodic boundary conditions.
            'Omega0': 0, #Matter density at z=0. !!Relevant only on cosmological integrations.
            'OmegaLambda': 0, #Vacuum energy density at z=0 !!Relevant only on cosmological integrations.
            'HubbleParam': 0.7, #Hublle constant in units of 100 km/s/Mpc !Relevant only on cosmological integrations.
            'FlagAge': 0, #Creation time of Stars !!Actually unused.
            'FlagMetals': 0, #Flag for metallicy values !!Currently unused.
            'NallHW': [0,0,0,0,0,0], #For simulation that use more than 2^32 particles, in our case set to 0.
            'flag_entr_ics': 0, #Flag for gas initial conditions: 0-ublock ----> Thermal energy, 1-ublock---> Entropy

             #Variables not in the standard Gadget header
            'Ntot': 0, #Total number of particle in the simulation.
            'filename': 'No file read' #Name of the last file read (not original Gadget parameter)

            }

    '''%%%%%%%%%%%%%%%%%%%'''
    ''' METHOD DEFINITION '''
    '''%%%%%%%%%%%%%%%%%%%'''

    '''////////////////////////////////////////////////'''
    '''Method: print
    '''
    def __str__(self):
        """
        :return: Print the values in the header.
        """
        check=0
        mess=""

        messt=""
        line="Npart: "+str(self.header["Npart"]) + " "*5
        messt+=line
        line="Massarr: "+str(self.header["Massarr"]) + " "*5
        messt+=line
        line="Ntot: "+str(self.header["Ntot"]) + " "*5
        messt+=line
        if len(messt)>check: check=len(messt)
        mess+=messt

        messt=""
        line="\nNall: "+str(self.header["Nall"]) + " "*5
        messt+=line
        line="NallHW: "+str(self.header["NallHW"]) + " "*5
        messt+=line
        line="NumFiles: "+str(self.header["NumFiles"]) + " "*5
        messt+=line
        if len(messt)>check: check=len(messt)
        mess+=messt

        messt=""
        line="\nTime: "+str(self.header["Time"]) + " "*5
        messt+=line
        line="Redshift: "+str(self.header["Redshift"]) + " "*5
        messt+=line
        line="BoxSize: "+str(self.header["BoxSize"]) + " "*5
        messt+=line
        line="HubbleParam: "+str(self.header["HubbleParam"]) + " "*5
        messt+=line
        line="Omega0: "+str(self.header["Omega0"]) + " "*5
        messt+=line
        line="OmegaLambda: "+str(self.header["OmegaLambda"]) + " "*5
        messt+=line
        if len(messt)>check: check=len(messt)
        mess+=messt

        messt=""
        line="\nFlagSfr: "+str(self.header["FlagSfr"]) + " "*5
        messt+=line
        line="FlagFeedback: "+str(self.header["FlagFeedback"]) + " "*5
        messt+=line
        line="FlagCooling: "+str(self.header["FlagCooling"]) + " "*5
        messt+=line
        line="FlagAge: "+str(self.header["FlagAge"]) + " "*5
        messt+=line
        line="FlagMetals: "+str(self.header["FlagMetals"]) + " "*5
        messt+=line
        line="FlagCooling: "+str(self.header["FlagCooling"]) + " "*5
        messt+=line
        if len(messt)>check: check=len(messt)
        mess+=messt

        messt=""
        line="\nflag_entr_ics: "+str(self.header["flag_entr_ics"]) + " "*5
        mess+=line
        line="filename: "+str(self.header["filename"]) + " "*5
        mess+=line
        if len(messt)>check: check=len(messt)
        mess+=messt


        int="| Header |"
        mess2=int.center(check,"*")+"\n"

        mess2+=mess
        mess2+="\n"+"*"*check


        return mess2
    '''////////////////////////////////////////////////'''



    '''////////////////////////////////////////////////'''
    '''Method: read_header( string filename )
    this method load the header information from a snapshot or IC file with name filename and stores it
    in the variable header.
    '''
    def read_header(self,filename,end='<'):

        #Check format
        wformat={'=':'Native', '<':'little-endian', '>':'big-endian'}
        if end not in wformat:
            print('Warning: incorrect  format type, it will be set to default value (little-endian)')
            end='<'
        f = open(filename, "rb") ##Open file in read binary mode ("rb")
        print ("\nReading header in file %s......" % (filename))

        block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block check (C-type int)

        #Now read the binary file following the format=1 in Gadget.
        self.header['Npart'][0] = list( struct.unpack(end+"iiiiii", f.read(24)) )
        self.header['Massarr'][0] = list( struct.unpack(end+"dddddd", f.read(48)) )
        self.header['Time'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['Redshift'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['FlagSfr'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['FlagFeedback'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['Nall'] =  list(struct.unpack(end+"iiiiii", f.read(24)))
        self.header['FlagCooling'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['NumFiles'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['BoxSize'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['Omega0'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['OmegaLambda'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['HubbleParam'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['FlagAge'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['FlagMetals'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['NallHW'][0] =  list(struct.unpack(end+"iiiiii", f.read(24)))
        self.header['flag_entr_ics'] = struct.unpack(end+"i", f.read(4))[0]

        #Now read the last unused byte (the total dimension of the header is 256.)
        f.read(256 - 24 - 48 - 8 - 8 - 4 - 4 - 24 - 4 - 4 - 8 - 8 - 8 - 8 - 4 - 4 - 24 - 4)

        block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the end block check (C-type int)

        #Now check if the initial and final block check is equal:
        if block_check_start != block_check_end: utility.Continue_check("Warning: Control Block failed")

        f.close() #close the file

        self.header['Ntot'] = sum(self.header['Npart'][0]) #Store the total number of particle in the header, we sum in Npart instead Nall
                                                           #because is safer. Indeed, in some Ics file the Nall file is left to 0.

        self.header['filename'] = filename #Store the name of the read file in the header

        print ("header data loaded from file %s \n" % (filename))
    '''////////////////////////////////////////////////'''

    '''////////////////////////////////////////////////'''
    '''METHOD: read_header_multi( string filename )
    this method load the header information from a snapshot or IC distributed on more than one file
    with names filename.i and stores it in the variable header.
    Note that in the n files the header information are always the same except for the variable  Npart.
    '''
    def read_header_multi(self,filename,end='<'):

        #Check format
        wformat={'=':'Native', '<':'little-endian', '>':'big-endian'}
        if end not in wformat:
            print('Warning: incorrect write format type, it will be set to default value (little-endian)')
            end='<'

        buff = filename + ".0" #Name of the first file.
        f = open(buff, "rb")
        print ("\nReading header in file " + buff)


        block_check_start = struct.unpack("<i", f.read(4))[0] #read the initial block check (C-type int)

        #Now read the binary file following the format=1 in Gadget.
        self.header['Npart'][0] = list(struct.unpack(end+"iiiiii", f.read(24)))
        self.header['Massarr'][0] = list(struct.unpack(end+"dddddd", f.read(48)))
        self.header['Time'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['Redshift'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['FlagSfr'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['FlagFeedback'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['Nall'] = list(struct.unpack(end+"iiiiii", f.read(24)))
        self.header['FlagCooling'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['NumFiles'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['BoxSize'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['Omega0'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['OmegaLambda'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['HubbleParam'] = struct.unpack(end+"d", f.read(8))[0]
        self.header['FlagAge'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['FlagMetals'] = struct.unpack(end+"i", f.read(4))[0]
        self.header['NallHW'][0] = list(struct.unpack(end+"iiiiii", f.read(24)))
        self.header['flag_entr_ics'] = struct.unpack(end+"i", f.read(4))[0]

        #Now read the last unused byte (the total dimension of the header is 256.)
        f.read(256 - 24 - 48 - 8 - 8 - 4 - 4 - 24 - 4 - 4 - 8 - 8 - 8 - 8 - 4 - 4 - 24 - 4)

        block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the end block check (C-type int)

        #Now check if the initial and final block check is equal:
        if block_check_start != block_check_end: utility.Continue_check("Warning: Control Header Block failed at file" + "." + str(i))

        f.close()


        #Now add the information of Nall from the other files
        for i in range(1, self.header["NumFiles"]):

            buff = filename + "." + str(i) #name of the current reading file
            f = open(buff, "rb")

            print ("Reading header in file " + buff)

            block_check_start = struct.unpack(end+"i", f.read(4))[0] #read the initial block check (C-type int)

            self.header["Npart"].append(list(struct.unpack(end+"iiiiii", f.read(24)))) #update the Npart information
            f.read(256 - 24) #Skip all the other byte that are equal among the files

            block_check_end = struct.unpack(end+"i", f.read(4))[0] #read the final block check (C-type int)

            #Now check if the initial and final block check is equal:
            if block_check_start != block_check_end: utility.Continue_check("Warning: Control Header Block failed at file" + "." + str(i))

            f.close()

        for i in range(self.header["NumFiles"]):
            self.header['Ntot'] += sum(self.header['Npart'][i]) #Store the total number of particle in the header, we sum in Npart instead Nall
                                                                   #because is safer. Indeed, in some Ics file the Nall file is left to 0.


        self.header['filename'] = filename #Store the name of the read file in the header

        print ("header data loaded from file %s. \n" % (filename))
    '''////////////////////////////////////////////////'''


    '''////////////////////////////////////////////////'''
    '''METHOD: write_header
    '''
    def write(self,stream,fclose=True,end='<'):
        """
        Write the header data on a file pointed by stream.
        :param stream: The name of the pointer to the open file (ex. f=open(..) stream=f)
        :param fclose: If True (default) after the writing routine the file (pointed by stream) will be closed.
        :return cbyte: Total byte written in the file
        """

        cbyte=0

        print ("\nWriting header in file " + stream.name)

        #Check format
        wformat={'=':'Native', '<':'little-endian', '>':'big-endian'}
        if end not in wformat:
            print('Warning: incorrect write format type, it will be set to default value (little-endian)')
            end='<'

        block_check_start = struct.pack(end+"i", 256) #write the initial block check (C-type int)
        stream.write(block_check_start)
        cbyte+=4 #integer start block check


        for i in range(6):
            buf=struct.pack(end+"i",self.header['Npart'][0][i])
            stream.write(buf)

        for i in range(6):
            buf=struct.pack(end+"d",self.header['Massarr'][0][i])
            stream.write(buf)

        buf=struct.pack(end+"d",self.header['Time'])
        stream.write(buf)

        buf=struct.pack(end+"d",self.header['Redshift'])
        stream.write(buf)

        buf=struct.pack(end+"i",self.header['FlagSfr'])
        stream.write(buf)

        buf=struct.pack(end,"i",self.header['FlagFeedback'])
        stream.write(buf)


        for i in range(6):
            buf=struct.pack(end+"i",self.header['Nall'][0][i])
            stream.write(buf)

        buf=struct.pack(end+"i",self.header['FlagCooling'])
        stream.write(buf)

        buf=struct.pack(end+"i",self.header['NumFiles'])
        stream.write(buf)

        buf=struct.pack(end+"d",self.header['BoxSize'])
        stream.write(buf)

        buf=struct.pack(end+"d",self.header['Omega0'])
        stream.write(buf)

        buf=struct.pack(end+"d",self.header['OmegaLambda'])
        stream.write(buf)

        buf=struct.pack(end+"d",self.header['HubbleParam'])
        stream.write(buf)

        buf=struct.pack(end+"i",self.header['FlagAge'])
        stream.write(buf)

        buf=struct.pack(end+"i",self.header['FlagMetals'])
        stream.write(buf)

        for i in range(6):
            buf=struct.pack(end+"i",self.header['NallHW'][0][i])
            stream.write(buf)

        buf=struct.pack(end+"i",self.header['flag_entr_ics'])
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

        if fclose==True: stream.close()

        print ("header write in file %s." % (stream.name))
        time.sleep(0.5)

        return cbyte
    '''////////////////////////////////////////////////'''

    '''////////////////////////////////////////////////'''
    '''METHOD: set_standard
    '''
    def set_standard(self, Gas=[0,0], Halo=[0,0], Disk=[0,0], Bulge=[0,0], Stars=[0,0], Bndry=[0,0]):
        """
        Set the header in the standard way.
        :param Gas: [0] Number of gas particles & [1] Mass of gas particles.
        :param Halo: [0] Number of halo particles & [1] Mass of halo particles.
        :param Disk: [0] Number of disk particles & [1] Mass of disk particles.
        :param Bulge: [0] Number of bulge particles & [1] Mass of bulge particles.
        :param Stars: [0] Number of star particles & [1] Mass of star particles.
        :param Bndry: [0] Number of bndry particles & [1] Mass of bndry particles.
        """



        self.header['Npart'][0] = [Gas[0],Halo[0],Disk[0],Bulge[0],Stars[0],Bndry[0]]
        self.header['Massarr'][0] = [Gas[1],Halo[1],Disk[1],Bulge[1],Stars[1],Bndry[1]]
        self.header['Nall'][0] =  self.header['Npart'][0]




        self.header['NumFiles'] = 1
        self.header['Time'] = 0.0
        self.header['Redshift'] = 0.0
        self.header['FlagSfr'] = 0
        self.header['FlagFeedback'] = 0
        self.header['FlagCooling'] = 0
        self.header['BoxSize'] = 0.0
        self.header['Omega0'] = 0.0
        self.header['OmegaLambda'] = 0.0
        self.header['HubbleParam'] = 1.0
        self.header['FlagAge'] = 0
        self.header['FlagMetals'] = 0
        self.header['NallHW'][0] =  [0,0,0,0,0,0]
        self.header['flag_entr_ics'] = 0
    '''////////////////////////////////////////////////'''

class Particle:
    def __init__(self,id=0,type=1,pos=(0,0,0),vel=(0,0,0),mass=0):
        self.Pos = pos  # (Cartesian) Coordinates of the particles (list of tuple) (C-type float)
        self.Vel = vel  # (Cartesian) componente of the particles velocities (list of tuple) (C-type float)
        self.Id = id  # Particle identification number (C-type int)
        self.Mass = mass  # Particle mass (C-type float)
        self.Type = type # Particle Type (C-type int): 0-Gas, 1-Halo, 2-Disk, 3-Bulge, 4-Stars, 5-Bndry

        # The following blocks exist in the snapshot only if enabled in the makefile
        self.Pot = 0  # Gravitational potential of the particles (C-type float)
        self.Acce = 0  # Acceleration of particles (C-type float)
        self.Tstp = 0  # Time at the start of simulation  (C-type float)

        # The following blocks exist only per SPH particles (Tyep-0)
        self.U = 0  # SPH Particles internal energy. (C-type float)
        self.Rho = 0  # SPH Particle density (C-type float)
        # The following blocks exist in the snapshot only if enabled in the makefile for SPH particles
        self.Hsml = 0  # Sph smoothing length
        self.Endt = 0  # Rate of change of entropy (if enabled in the makefile)

        # The following variables are not imported from file and are used only in internal method of the class
        self.Radius = 0  # Distance of the particle from the axes-origin.
        self.Vel_tot =0  # Particle total velocity
        self.Energy = 0  # Particle total energy (in the sense of epsilon in BT)
        self.Pcord = "Pos not stored"
        self.Vcord = "Vel not stored"

        self.setRadius()
        self.setVel()
    '''%%%%%%%%%%%%%%%%%%%'''
    ''' METHOD DEFINITION '''
    '''%%%%%%%%%%%%%%%%%%%'''

    '''////////////////////////////////////////////////'''
    ''' Method:reset() '''
    '''
    The method simply reset the object to the initial null values
    '''

    def reset(self):
        self.Pos = (0,0,0)
        self.Vel = (0,0,0)
        self.Id = 0
        self.Mass = 0
        self.Type = 1
        self.Potential = 0
        self.Acce = 0
        self.Tstp = 0
        self.Radius = 0
        self.Vel_tot = 0
        self.Energy = 0

        self.U = 0
        self.Rho = 0
        self.Endt = 0
        self.Hsml = 0

    '''////////////////////////////////////////////////'''

    def print_general(self):
        print('\n')
        print("General".center(35, "*"))
        print('Id:', self.Id, '|', 'Type:', self.Type, '|', 'Mass:', self.Mass)

    def print_pos(self):
        print("Position".center(50, "*"))
        st = "Coo system: " + self.Pcord
        print(st.center(50, " "))
        print('X: ', self.Pos[0], 'Y: ', self.Pos[1], 'Z: ', self.Pos[2])
        try:
            radius = sqrt((self.Pos[0]) ** 2 + (self.Pos[1]) ** 2 + (self.Pos[2]) ** 2)
        except:
            radius = None
        print('Radius:', self.Radius, 'Radius_calc:', radius)

    def print_vel(self):
        print("Velocity".center(50, "*"))
        st = "Coo system: " + self.Pcord
        print(st.center(50, " "))
        print('Vx: ', self.Vel[0], 'Vy: ', self.Vel[1], 'Vz: ', self.Vel[2])
        try:
            vel = sqrt((self.Vel[0]) ** 2 + (self.Vel[1]) ** 2 + (self.Vel[2]) ** 2)
        except:
            vel = None
        print('Veltot:', self.Vel_tot, 'Veltot_calc:', vel)

    def __str__(self):
        mess = ""

        line = '\n' + "General".center(35, "*") + '\n' + 'Id: ' + str(self.Id) + '|' + 'Type: ' + str(
            self.Type) + '| ' + 'Mass:' + str(self.Mass) + '\n'
        mess += line

        st = "Coo system: " + self.Pcord
        line = "Position".center(50, '*') + '\n' + st.center(50, " ") + '\n' + 'X: ' + str(self.Pos[0]) + ' Y: ' + str(
            self.Pos[1]) + ' Z: ' + str(self.Pos[2]) + '\n'
        rad = sqrt((self.Pos[0]) ** 2 + (self.Pos[1]) ** 2 + (self.Pos[2]) ** 2)
        line += 'Radius: ' + str(self.Radius) + ' Radius_calc: ' + str(rad) + '\n'
        mess += line

        st = "Coo system: " + self.Vcord
        line = "Velocity".center(50, "*") + '\n' + st.center(50, " ") + '\n' + 'Vx: ' + str(
            self.Vel[0]) + ' Vy: ' + str(self.Vel[1]) + ' Vz: ' + str(self.Vel[2]) + '\n'
        vel = sqrt((self.Vel[0]) ** 2 + (self.Vel[1]) ** 2 + (self.Vel[2]) ** 2)
        line += 'Veltot: ' + str(self.Vel_tot) + ' Veltot_calc: ' + str(vel)
        mess += line

        return mess

    def setRadius(self):
        self.Radius = sqrt(self.Pos[0] * self.Pos[0] + self.Pos[1] * self.Pos[1] + self.Pos[2] * self.Pos[2])

    def setVel(self):
        self.Vel_tot = sqrt(self.Vel[0] * self.Vel[0] + self.Vel[1] * self.Vel[1] + self.Vel[2] * self.Vel[2])

class Particles:
    """
    NB. In the case p is prenset, the header will be not counted and a new header from the data will be
    created
    """

    def __init__(self,p=None,h=None,N=1):

        if ((p is None) and (h is None)):
            self.n=N
            self._initialize_vars(self.n) #Create from scratch
            self._make_header()


        elif p is None:
            if (isinstance(h,Header)): self.header=h.header
            else: self.header=h

            self.n=np.sum(self.header['Nall'])
            self._initialize_vars(self.n)
            #id
            self._makeid()
            #make type and mass
            self._maketypemass()


        else:
            if not (isinstance(p, Particle) or ((isinstance(p, (list, tuple, np.ndarray))) and (isinstance(p[0], Particle)))): raise ValueError('Incorrect particle format')
            if isinstance(p, Particle): p=[p]
            self.n = (len(p))
            self._initialize_vars(self.n) #Create from scratch
            self._fill_from_particle(p)
            self._make_header()

        self.par_dic = {'id': self.Id, 'type':self.Type, 'mass':self.Mass, 'rad':self.Radius, 'velt':self.Vel_tot ,'x': self.Pos[:,0], 'y': self.Pos[:,1], 'z': self.Pos[:,2], 'vx': self.Vel[:,0], 'vy': self.Vel[:,1], 'vz': self.Vel[:,2], 'pos': self.Pos, 'vel':self.Vel}




    def _initialize_vars(self,n):
            self.Pos = np.zeros(shape=[n, 3], dtype=float)
            self.Vel = np.zeros(shape=[n, 3], dtype=float)
            self.Id = np.arange(1,n+1, dtype=int)
            self.Mass = np.zeros(n, dtype=float)
            self.Type = np.full(n, 1, dtype=object)

            # The following blocks exist in the snapshot only if enabled in the makefile
            self.Pot = np.full(n, 0, dtype=float)  # Gravitational potential of the particles (C-type float)
            self.Acce = np.zeros(shape=[n, 3], dtype=float) # Acceleration of particles (C-type float)
            self.Tstp = np.full(n, 0, dtype=float)  # Time at the start of simulation  (C-type float)


            # The following blocks exist only per SPH particles (Tyep-0)
            self.U = np.full(n, 0, dtype=float)  # SPH Particles internal energy. (C-type float)
            self.Rho = np.full(n, 0, dtype=float)  # SPH Particle density (C-type float)
            # The following blocks exist in the snapshot only if enabled in the makefile for SPH particles
            self.Hsml = np.full(n, 0, dtype=object)  # Sph smoothing length
            self.Endt = np.full(n, 0, dtype=object)  # Rate of change of entropy (if enabled in the makefile)

            # The following variables are not imported from file and are used only in internal method of the class
            self.Radius = np.full(n, 0, dtype=float)  # Distance of the particle from the axes-origin.
            self.Vel_tot = np.full(n, 0, dtype=float)  # Particle total velocity
            self.Energy = np.full(n, 0, dtype=float)  # Particle total energy (in the sense of epsilon in BT)
            self.Pcord = np.full(n, 'Pos not stored', dtype=object)
            self.Vcord = np.full(n, 'Vel not stored', dtype=object)

            # History variables: variables that store the operation made on the object
            self.hpos = [(0, 0, 0)]  # List of  position in wich the partciles have been moved (The last triple, are actually the coordinate of the reference system for the particles)
            self.hvel = [(0, 0,0)]  # List of vel position in wich the partciles have been moved (The last triple, are actually the velocity of the reference system for the particles)
            self.heul = [(0, 0, 0)]  # List of euler angles  in wich the partciles have been rotated
            self.order_var=None

    def _make_header(self):
        h=Header()
        self.header=h.header

        count=[ np.sum(self.Type==i) for i in range(7)]

        self.header['Ntot']=self.n
        self.header['Npart']=[count]
        self.header['Nall']=count

    def _fill_from_particle(self,p):
        for i in range(self.n):
                self.Pos[i] = p[i].Pos  # (Cartesian) Coordinates of the particles (list of tuple) (C-type float)
                self.Vel[i] = p[i].Vel  # (Cartesian) componente of the particles velocities (list of tuple) (C-type float)
                self.Id[i] = p[i].Id  # Particle identification number (C-type int)
                self.Mass[i] = p[i].Mass  # Particle mass (C-type float)
                self.Type[i] = p[i].Type  # Particle Type (C-type int): 0-Gas, 1-Halo, 2-Disk, 3-Bulge, 4-Stars, 5-Bndry

                # The following blocks exist in the snapshot only if enabled in the makefile
                self.Pot[i] = p[i].Pot  # Gravitational potential of the particles (C-type float)
                self.Acce[i] = p[i].Acce  # Acceleration of particles (C-type float)
                self.Tstp[i] = p[i].Tstp  # Time at the start of simulation  (C-type float)

                # The following blocks exist only per SPH particles (Tyep-0)
                self.U[i] = p[i].U  # SPH Particles internal energy. (C-type float)
                self.Rho[i] = p[i].Rho  # SPH Particle density (C-type float)
                # The following blocks exist in the snapshot only if enabled in the makefile for SPH particles
                self.Hsml[i] = p[i].Hsml  # Sph smoothing length
                self.Endt[i] = p[i].Endt  # Rate of change of entropy (if enabled in the makefile)

                # The following variables are not imported from file and are used only in internal method of the class
                self.Radius[i] = p[i].Radius  # Distance of the particle from the axes-origin.
                self.Vel_tot[i] = p[i].Vel_tot  # Particle total velocity
                self.Energy[i] = p[i].Energy  # Particle total energy (in the sense of epsilon in BT)
                self.Pcord[i] = p[i].Pcord
                self.Vcord[i] = p[i].Vcord

    def _makeid(self):
        self.Id=np.arange(0,self.n)

    def _maketypemass(self):
        i=0
        maxincum=0
        for maxin in self.header['Nall']:
            self.Type[maxincum:maxin+maxincum]=i
            self.Mass[maxincum:maxin+maxincum]=self.header['Massarr'][0][i]
            i+=1
            maxincum+=maxin

    def _maketype(self):
        i=0
        maxincum=0
        for maxin in self.header['Nall']:
            self.Type[maxincum:maxin+maxincum]=i
            i+=1
            maxincum+=maxin

    def order(self,key='Id'):

        allowed_order_keys=('Id','Mass','Type','Pot','Acce','U','Rho','Radius','Vel_tot','Energy')
        if key not in allowed_order_keys: raise ValueError('key: %s. Not supported order key',key)

        if key=='Id':           sort_idx=np.argsort(self.Id)
        elif key=='Radius':     sort_idx= np.argsort(self.Radius)
        elif key=='Vel_tot':    sort_idx= np.argsort(self.Vel_tot)
        elif key=='Energy':     sort_idx= np.argsort(self.Energy)
        elif key=='Mass':       sort_idx= np.argsort(self.Mass)
        elif key=='Type':       sort_idx= np.argsort(self.Type)
        elif key=='Pot':        sort_idx= np.argsort(self.Pot)
        elif key=='Acce':       sort_idx= np.argsort(self.Acce)
        elif key=='U':          sort_idx= np.argsort(self.U)
        elif key=='Rho':        sort_idx= np.argsort(self.Rho)

        self.Pos=self.Pos[sort_idx]
        self.Vel=self.Vel[sort_idx]
        self.Id=self.Id[sort_idx]
        self.Mass=self.Mass[sort_idx]
        self.Type=self.Type[sort_idx]
        self.Pot= self.Pot[sort_idx]
        self.Acce=self.Acce[sort_idx]
        self.Tstp=self.Tstp[sort_idx]
        self.U=self.U[sort_idx]
        self.Rho=self.Rho[sort_idx]
        self.Hsml=self.Hsml[sort_idx]
        self.Endt=self.Endt[sort_idx]
        self.Radius=self.Radius[sort_idx]
        self.Vel_tot=self.Vel_tot[sort_idx]
        self.Pcord=self.Pcord[sort_idx]
        self.Vcord=self.Vcord[sort_idx]
        self.order_var=key

        return sort_idx

    def setrad(self):
        m=self.Pos*self.Pos #Matrix whit the pos coordinate**2
        self.Radius=np.sqrt(m.sum(axis=1)) #Sum the quadratic coordinate for each particle (axis=1, row) ans square it

    def setvelt(self):
        m=self.Vel*self.Vel #Matrix whit the vel componentes**2
        self.Vel_tot=np.sqrt(m.sum(axis=1)) #Sum the quadratic vel compontens for each particle (axis=1, row) ans square it

    def rtmove(self, mpos=(0,0,0), mvel=(0,0,0), eangle=(0,0,0)):
        '''
        This function applies a rototraslation for the particle system.
        The traslations are applied before of the rotation.
        :param mpos: position vector mx to change position coordinate (xnew=x-mpos)
        :param mvel: velocity angle mv to change position coordinate (vnew=v-mvel)
        :param eangle: Euler angle to rotate the system following the zyz convenction
        :return:
        '''

        mr=mt.sqrt(mpos[0]*mpos[0]+mpos[1]*mpos[1]+mpos[2]*mpos[2])
        mv=mt.sqrt(mvel[0]*mvel[0]+mvel[1]*mvel[1]+mvel[2]*mvel[2])
        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
        print('Class particles: Function rototranslation')

        if (mr==0): print('The frame of reference will be no moved')
        else:
            print('The frame of reference will be moved following the vector:')
            print('{0:<s} {1:>6.3f} {2:<s} {3:>6.3f} {4:<s} {5:>6.3f} {6:<s} {7:>6.3f}'.format('Xm:', mpos[0] ,'Ym:', mpos[1], 'Zm:',mpos[2], 'Module:', mr))
            self.Pos=self.Pos-mpos #Apply the traslation
            print('Done')

        if (mv==0): print('The linear motion of the frame of reference will be no changed')
        else:
            print('The linear motion of the frame of reference will be changed following the vector:')
            print('{0:<s} {1:>6.3f} {2:<s} {3:>6.3f} {4:<s} {5:>6.3f} {6:<s} {7:>6.3f}'.format('Vxm:', mvel[0] ,'Vym:', mvel[1], 'Vzm:',mvel[2], 'Module:', mv))
            self.Vel=self.Vel-mvel #Apply the traslation in velocity
            print('Done')

        if (eangle[0]==0 and eangle[1]==0 and eangle[2]==0): print('The frame of reference will be not rotated')
        else:
            print('The frame of reference will be rotated following the y-convention:')
            print('{0:<40s} {1:>6.3f} {2:>s}'.format('First rotation around the z-axis:',eangle[0],'deg'))
            print('{0:<40s} {1:>6.3f} {2:>s}'.format('Second rotation around the y\'-axis:',eangle[1],'Deg'))
            print('{0:<40s} {1:>6.3f} {2:>s}'.format('Third rotation around the z\'\'-axis:',eangle[2],'Deg'))

        #Define the the euler angle in radiant
        a=eangle[0]*2*pi/360
        b=eangle[1]*2*pi/360
        g=eangle[2]*2*pi/360

        #Define the matrix element
        ca=cos(a)
        sa=sin(a)
        cb=cos(b)
        sb=sin(b)
        cg=cos(g)
        sg=sin(g)

        #Define the rotation matrix following the y-convention
        R=np.array([[cg*cb*ca - sg*sa, cg*cb*sa + sg*ca  , -cg*sb ],[ -sg*cb*ca -cg*sa , -sg*cb*sa + cg*ca , sg*sb ],[ sb*ca ,  sb*sa ,  cb ]])

        '''
        NB: as defined above the rotation matrix works on column vector v:
        v'=R*v.
        But the position and velocity coordinate of the class particles have been defined
        as row vector w. In this case the transformation is:
        w'=w*R^T, where R^T is the transpose matrix of R
        '''

        #Rotate the system
        self.Pos=np.dot(self.Pos,R.T) #matrix algebric multiplication
        self.Vel=np.dot(self.Vel,R.T)

        print('Done')

        #Set new radius and total velocity
        self.setrad()
        self.setvelt()


        #Store this rototraslation in the hystorical variables
        self.hpos.append(mpos)
        self.hvel.append(mvel)
        self.heul.append(eangle)

        print('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

    def __getitem__(self, key):
        """ Get method

        :param key: int, slice, str or boolean array
        :return: if key is int return a Particle object, if slice or boolena array return a  list of Particle objects, if
        str return the quantities depending on the attribute par_dic
        """

        if isinstance(key, int):
            return self._getitem(key)
        elif isinstance(key, slice):
            return self._getslice(key)
        elif isinstance(key,list) or isinstance(key,np.ndarray) or isinstance(key,tuple):
            return self._getboolean(np.array(key))
        elif isinstance(key,str):
            return self._getdic(key)
        else:
            raise TypeError('Index must be int, not {}'.format(type(key).__name__))

    def _getitem(self,id):
        """Get a single element

        :param id: int
        :return:  a Particle element with the properties of the ith particle
        """

        pos=self.Pos[id]
        vel=self.Vel[id]
        mass=self.Mass[id]
        type=self.Type[id]

        part=Particle(id=id,type=type,pos=pos,vel=vel,mass=mass)

        return part

    def _getslice(self,idslice):
        """Return a slice of Particles

        :param idslice:  slice object
        :return: a list of Particle objects sliced from the original one
        """



        if idslice.start is None: start=0
        else: start=idslice.start

        if idslice.stop is None: stop=self.n
        else: stop=idslice.stop

        if idslice.step is None: step=1
        else: step=idslice.step

        plist=[]

        for i in range(start,stop,step):
            plist.append(Particle(id=self.Id[i], type=self.Type[i], pos=self.Pos[i], vel=self.Vel[i], mass=self.Mass[i]))


        return np.array(plist)

    def _getboolean(self,bollist):
        """Return a list of Particle objects where bollist is True

        :param bollist: a boolean array with the dimension of Particles
        :return:  a list of Particle objects
        """

        pos=self.Pos[bollist]
        vel=self.Vel[bollist]
        mass=self.Mass[bollist]
        type=self.Type[bollist]
        id=self.Id[bollist]

        leni=len(id)

        plist=np.zeros(leni,dtype=object)

        for i in range(leni):
            plist[i]=Particle(id=id[i],type=type[i],pos=pos[i],vel=vel[i],mass=mass[i])

        return plist

    def _getdic(self,strg):
        """Return the attribute following the dic in par_dic

        :param strg: strg
        :return:  an array with attribute defined in strg
        """

        if strg.lower() in self.par_dic:
            return self.par_dic[strg.lower()]
        else:
            raise ValueError('par %s not in par_dic'%strg)

    def extract(self,*args,mode='and'):
        """This method create a new Particles or Sky_Particles objects filtering the original one
        with condition in args. If there are no particles that fulfill all the conditions the function return None.
        The args need to be str object with the name of attribute, the condition (<,>,=,<=,>=) and the numerical value, e.g.:

        'Id>10', 'type=2', 'mass>30'

        All the condition

        :param args: list of str
        :return: a new Particles or Sky_Particles objects or None
        """

        if mode.lower()=='and': idx=np.ones(self.n,dtype=bool)
        elif mode.lower()=='or':  idx=np.zeros(self.n,dtype=bool)
        else: raise ValueError('Invalid mode, it hsould be and or or')

        for conditions in args:

            if isinstance(conditions,str)==False: raise ValueError('Condition need to be a str')

            if mode.lower()=='and': idx=self._check_and(conditions,idx)
            elif mode.lower()=='or': idx=self._check_or(conditions,idx)

        truecount=np.sum(idx)

        if truecount==0:

            print('Warning 0 Filtered object.. returning None')
            return None

        part_list=self._getboolean(idx)

        return Particles(p=part_list)


    def _check_and(self,cond,idx):
        """

        :param cond: str object
        :param idx: boolean numpy array
        :return:
        """

        cond, key, num = find_symbol(cond)
        check_list = self._getdic(key)
        num = float(num)

        if cond == '<=': idx *= check_list <= num
        elif cond == '<': idx *= check_list < num
        elif cond == '>=': idx *= check_list >= num
        elif cond == '>': idx *= check_list > num
        elif cond == '=': idx *= check_list == num
        else: raise ('Invalid condition')

        return idx

    def _check_or(self, cond, idx):
        """

        :param cond: str object
        :param idx: boolean numpy array
        :return:
        """

        cond, key, num = find_symbol(cond)
        check_list = self._getdic(key)
        num = float(num)

        if cond == '<=': idx += check_list <= num
        elif cond == '<': idx += check_list < num
        elif cond == '>=': idx += check_list >= num
        elif cond == '>': idx += check_list > num
        elif cond == '=': idx += check_list == num
        else: raise ('Invalid condition')

        return idx



if __name__=='main':
    p1=Particle(pos=(1,0,0))
    p2=Particle(id=1,pos=(3,0,0),type=4)
    p=np.array([p1,p2])
    pp=Particles(N=10)
    print(pp.header)
    pp=Particles(p=p)
    print(pp.header)
    h=Header()
    h.header['Nall']=[1,4,2,0,0,1]
    h.header['Massarr']=[[1,2,3,0,0,1]]
    pp=Particles(h=h,p=p)
    print(pp.Id)
    print(pp.Type)
    print(pp[1])
