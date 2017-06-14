__author__ = 'Giuliano'

import os
import time
import numpy as np

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

