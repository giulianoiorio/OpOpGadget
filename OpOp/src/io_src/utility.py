from __future__ import  division, print_function

import numpy as np

length_unit_cm={'mm':1e-1,'cm':1,'dm':1e1,'m':1e2,'hm':1e3,'km':1e4,'Mm':1e7,'Gm':1e10,'UA':1.496e+13, 'pc':3.085678e18, 'hpc':3.085678e20, 'kpc':3.085678e21, 'Mpc':3.085678e24, 'Gpc':3.085678e27 }
velocity_unit_cms={'mms':1e-1,'cms':1,'dms':1e1,'ms':1e2,'hms':1e3,'kms':1e4,'Mms':1e7,'Gms':1e10,'UAs':1.496e+13, 'pcs':3.085678e18, 'hpcs':3.085678e20, 'kpcs':3.085678e21, 'Mpcs':3.085678e24, 'Gpcs':3.085678e27 }
mass_unity_g={'g':1, 'msun':1.989e33, '1msun':1.989e34, '2msun':1.989e35, '3msun':1.989e36, '4msun':1.989e37, '5msun':1.989e38, '6msun':1.989e39, '7msun':1.989e40, '8msun':1.989e41, '9msun':1.989e42, '10msun':1.989e43    }

def Gadget_parfil(tfin,outfile='param.param',tini=0,lunit='kpc',vunit='kms',munit='msun'):


    vu=velocity_unit_cms[vunit]
    mu=mass_unity_g[munit]
    lu=length_unit_cm[lunit]
    tu=lu/vu #in sec
    sec_per_year=3.1536e7
    sec_per_Gyear=3.1536e16

    print(vu,mu,lu,tu/sec_per_Gyear)

    h=''

    #Sure
    h+='\nIntiCondFile'
    h+='\nOutputDir'


    #Code options
    h+='\n%Code options'
    h+='\n{:30} {:d}'.format('ICFormat',1)
    h+='\n{:30} {:d}'.format('SnapFormat', 1)
    h+='\n{:30} {:d}'.format('ComovingIntegrationOn', 0)
    h+='\n{:30} {:d}'.format('TypeOfTimestepCriterion',0)
    h+='\n{:30} {:d}'.format('OutputListOn', 0)
    h+='\n{:30} {:d}'.format('PeriodicBoundariesOn', 0)

    h+='\n\n%Code options'


    return h

if __name__ == "__main__":

    print(Gadget_parfil(2))