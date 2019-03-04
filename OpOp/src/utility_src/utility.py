from __future__ import  division, print_function


import os
import time
import numpy as np
import astropy.io.fits as ft

def radec_to_xieta(ra, dec, ra_c, dec_c):
    
    degtrad=np.pi/180
    sd=np.sin(dec*degtrad)
    cd=np.cos(dec*degtrad)
    sdc=np.sin(dec_c*degtrad)
    cdc=np.cos(dec_c*degtrad)
    
    
    denominator = sd*sdc + cd*cdc*np.cos( (ra-ra_c)*degtrad )
    xin= cd * np.sin( (ra-ra_c)*degtrad  ) / denominator
    xi=xin/degtrad
    
    etan= sd*cdc - cd*sdc*np.cos( (ra-ra_c)*degtrad )
    eta = etan/degtrad

    return xi*3600, eta*3600
    
def cartesian_to_spherical_vector(Ax, Ay, Az, theta, phi, degree=True):
    """
    Convert the component of a vector from cartesian to spherical coordinates
    :param Ax: x component of the vector
    :param Ay: y component of the vector
    :param Az: z component of the vector
    :param theta: zenithal  angle (Arccos(z/r)) of the point
    :param phi:  azimuthal angle (Arctan2(y/x))  of the point
    :param degree:  If True, theta and phi are in degrees, if False in radians.
    :return:
    """
    if degree: degree_to_rad=np.pi/180.
    else: degree_to_rad=1

    stheta=np.sin(theta*degree_to_rad)
    ctheta=np.cos(theta*degree_to_rad)
    sphi=np.sin(phi*degree_to_rad)
    cphi=np.cos(phi*degree_to_rad)
    
    Ar     =  Ax*stheta*cphi + Ay*stheta*sphi + Az*ctheta
    Atheta =  Ax*ctheta*cphi + Ay*ctheta*sphi - Az*stheta
    Aphi   = -Ax*sphi        + Ay*cphi
    
    return Ar, Atheta, Aphi



def list_check(object_to_be_checked):
    
    return isinstance(object_to_be_checked,(list, tuple, np.ndarray))



def Continue_check(warning):
    valid_choice = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}

    check=True
    while check==True:
        choice = input(warning + "\nContinue?" + "[Y/n]").lower()


        if choice in valid_choice or choice=='':
            check=False
            if choice=='' or valid_choice[choice]==True:
                    print ("Continue running.....")
                    pass
            elif valid_choice[choice]==False:
                    print ("Abort.....")
                    exit()
        else:
            print ("Invalid input: please respond yes or no (ye, y and n allowed)")

def file_check(filename):
    """
    :param filename:
    :return:
    """

    check_file=os.path.isfile(filename)

    if check_file:

        valid_choice = ("a","n","e")

        check=True
        while check==True:
            choice=input("The output file " + "\"" + str(filename) + "\" " +  "already exists. Waiting for instruction..\n(a for append, n for change name, e to exit, enter to overwrite)\n").lower()

            if choice in valid_choice or choice=='':
                if choice=="n": choice=str(input("New name:"))
                elif choice=="e":
                    print ("Abort.....")
                    exit()

                check=False



            else:
                print ("!!!!!Invalid input!!!!!!")
                time.sleep(0.3)

    else: choice=''

    return choice

def nparray_check(value_to_check):

    if isinstance(value_to_check,np.ndarray): return  value_to_check
    elif (isinstance(value_to_check,list) or isinstance(value_to_check,tuple)): return np.array(value_to_check)
    elif (isinstance(value_to_check,int) or isinstance(value_to_check,float)): return np.array([value_to_check,])
    else: raise TypeError('Non valid format, use an nparray, a list a tyle a int or a float')

def find_symbol(strg):

    for s in ('>=','<=','<','>','='):
        a = strg.split(s)
        if len(a)==2:
            return (s,a[0],a[1])
    raise ValueError('logic symbol error, it should be >, <, >=, <=, =')


def make_fits(dict, outname=None, header_key={}):
    '''
    Make a fits table from a numpy array
    args must be dictionary containing the type  and the columnf of the table, e.g.
    {'l':(col1,'D'),'b':(col2,'D')}
    '''

    col = []
    for field in dict:
        if len(dict[field]) == 2:
            format = dict[field][1]
            array = dict[field][0]
        else:
            format = 'D'
            array = dict[field]

        col.append(ft.Column(name=field, format=format, array=array))

    cols = ft.ColDefs(col)
    tab = ft.BinTableHDU.from_columns(cols)
    for key in header_key:
        item = header_key[key]
        if item is None:
            tab.header[key] = str(item)
        else:
            tab.header[key] = item

    if outname is not None: tab.writeto(outname, clobber=True)

    return tab


def stat_error(y, yerr, err_type='prop', n=1000, bootstrap_error=False, vsys_fix=None):
    """
    """

    N = len(y)
    err_type == err_type.lower()

    if err_type[0] == 'b' or err_type[0] == 'j':

        if err_type[0] == 'b':

            if bootstrap_error:
                yerr_boot = yerr
            else:
                yerr_boot = None
            data_sample = bootstrap_resampling(y, nsample=n, dataerr=yerr_boot)

        elif err_type[0] == 'j':

            data_sample = jackknife_resampling(y)

        if vsys_fix is None:

            mean_list = np.mean(data_sample, axis=1)
            std_list = np.std(data_sample, axis=1)

            mean = np.mean(mean_list)
            std = np.mean(std_list)

            mean_error = np.std(mean_list)
            std_error = np.std(std_list)

        else:

            mean = vsys_fix
            mean_error = np.nan

            std_list = std_fix(data_sample, mean, axis=1)
            std = np.mean(std_list)
            std_error = np.std(std_list)


    else:
        y = np.array(y)
        yerr = np.array(yerr)

        if vsys_fis is None:

            mean = np.mean(y)
            std = np.std(y)

            mean_error = 1. / N * np.sqrt((np.sum(yerr * yerr)))

        else:

            mean = vsys_fix
            mean_error = np.nan

            std = std_fix(y, mean)

    if err_type[0] == 'p':
        diff = y - mean
        std_err_a = np.sum(diff * diff * yerr * yerr)

        sumdiff = np.sum(diff)
        std_err_b = sumdiff * sumdiff * mean_error * mean_error

        std_error = (1 / (N * std)) * np.sqrt(std_err_a + std_err_b)

    elif err_type[0] == 'g':
        err_mean = np.mean(yerr)
        sigma_tot2 = std * std + err_mean * err_mean

        deltavar = np.sqrt(2 * sigma_tot2 * sigma_tot2) / N
        std_error = deltavar / (2 * std)

    return mean, mean_error, std, std_error


def bootstrap_resampling(data, fraction=1, dataerr=None, nsample=1):
    N = len(data)

    N_resample = int(N * fraction)

    if dataerr is not None:
        w = 1 / dataerr
        w = w / np.sum(w)
    else:
        w = None

    data_resample = np.random.choice(data, size=N_resample * nsample, replace=True, p=w)

    out = data_resample.reshape((nsample, N_resample))

    return out