from __future__ import  division, print_function

import numpy as np
from scipy.interpolate import  UnivariateSpline
from numpy.ctypeslib import ndpointer
import ctypes as ct
import  os


'''
def df_isotropic(dens,pot,egrid=None,normalize=False,use_c=False, check_negative=True):
    """
    Calcola la funzione di distribuzione per un sistema sferico isotropico con f(e) dove
    e è l'energia efficace normalizzata.
    :param dens: griglia con i valori di densità della componente, normalizzata a 1 nel valore al raggio piu basso
    :param pot: griglia con i valori del potenziale totale, normalizzata ad 1 nel raggio piu basso
    :param egrid: griglia di energia dove calcolare la df numerica. se lasciata ==none, la griglia in energia è uguale a quella del potenziale.
    :param plot: if True print the df and the logaritmic df
    :return: griglia usata, df numerica sulla griglia, df numerica funzione
    """

    #1 Definisci griglia pot-dens normalizzata con l'aggiunta del punto 0-0
    potnew=np.zeros(len(pot)+1,dtype=float,order='C')
    potnew[-1]=0
    potnew[:-1]=pot
    densnew=np.zeros(len(dens)+1,dtype=float,order='C')
    densnew[-1]=0
    densnew[:-1]=dens
    if normalize==True:
        potnew=potnew/np.max(potnew)
        densnew=densnew/np.max(densnew)


    if use_c==False:
        dphi=InterpolatedUnivariateSpline(potnew[::-1],densnew[::-1],k=3,ext=2)
        dphider=dphi.derivative() #derivata di tho rispetto a Psi



        #Semplicemente abbiamo aggiunto un punto a pot=0 dato che in questo punto accade sempre che dens(0)=0 (All infinito sia la densità che il potenziale per definizione valgono 0)

        #definisci griglia di energia se npn è data in input
        if egrid==None: egrid=np.copy(potnew)
        elif (np.max(egrid)>1) | (np.min(egrid)<0): raise ValueError('The energy grid allowed value in the range 0-1')

        inte=np.zeros_like(egrid) #Vettore per accogliere i risultati

        #Definisci funzione integranda:
        func= lambda potential,energy: dphider(potential)/np.sqrt(energy-potential)




        #Integra per ogni e con il metodo delle quadature di gauss
        i=0
        for e in egrid:
            if e==0: inte[i]=0
            else: inte[i]=integrate.fixed_quad(func,0,e,args=(e,))[0]
            i+=1


        #Interpoalte the results
        dfe=InterpolatedUnivariateSpline(egrid[::-1],inte[::-1],k=3,ext=2)



        dfeder=dfe.derivative() #derivative

        dfgrid=dfeder(egrid)
        dfgrid[-1]=0 #Impongo a mano che a e=0 la df=0

        #dfe_norm=InterpolatedUnivariateSpline(egrid[::-1],dfeder(egrid[::-1])/np.max(dfeder(egrid[::-1])),k=2,ext=2)
        #dfe_norm=InterpolatedUnivariateSpline(egrid[::-1],dfgrid[::-1]/np.max(dfgrid[::-1]),k=2,ext=2)
        dfe_norm=InterpolatedUnivariateSpline(egrid[::-1],dfgrid[::-1],k=2,ext=2)

        if check_negative:
            if (np.any((dfe_norm(egrid)<0) & (np.abs(dfe_norm(egrid))>1e-3))): raise ValueError('negative df')

        return (egrid,dfe_norm(egrid),dfe_norm)

    else:
        dfe_norm=np.zeros_like(potnew,dtype=float,order='C')
        lib=ct.cdll.LoadLibrary("/Users/Giuliano/PycharmProjects/OpOp/OpOp/df/df_c_ext/df_spherical.so")
        df_func=lib.df_spherical
        df_func.restype=None
        df_func.argtypes=[ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ct.c_int,ndpointer(ct.c_double, flags="C_CONTIGUOUS")]
        df_func(densnew,potnew,len(densnew),dfe_norm)
        dfe_norm=dfe_norm/np.max(dfe_norm)
        dfe_norm_func=InterpolatedUnivariateSpline(egrid[::-1],dfe_norm[::-1]/np.max(dfe_norm[::-1]),k=2,ext=2)

        return (potnew,dfe_norm,dfe_norm_func)
'''


def df_isotropic(dens,pot,**kwargs):
    """
    Calculate the distribution function. This is scaled by some constant C, but this is not
    important given that the df is always normalized to the value of the max energy

    :param dens: Array with density
    :param pot:  Array with potential
    :return:
            -pot: grid of energy
            -df: distribution function sampled in the grid pot
            -df_func: df function of variable e (interpolation of pot-df)
    """



    df=np.zeros(len(dens),dtype=float,order='C')

    if 'use_c' in kwargs: use_c=kwargs['use_c']
    else: use_c=True



    if use_c==True:
        pot=np.ascontiguousarray(pot,dtype=float)
        dens=np.ascontiguousarray(dens, dtype=float)



        #add to path to use relative path
        dll_name='df_c_ext/df_spherical.so'
        dllabspath = os.path.dirname(os.path.abspath(__file__)) + os.path.sep + dll_name
        lib = ct.CDLL(dllabspath)
        #add to path to use relativ path



        df_func=lib.df_spherical
        df_func.restype=None
        df_func.argtypes=[ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ndpointer(ct.c_double, flags="C_CONTIGUOUS"),ct.c_int,ndpointer(ct.c_double, flags="C_CONTIGUOUS")]

        df_func(dens,pot,len(dens),df)



    else:
        drpsi=np.zeros(len(dens),dtype=float,order='C')
        inte=np.zeros(len(dens),dtype=float,order='C')


        for i in range(len(dens)-1):
            drpsi[i]=(dens[i+1]-dens[i])/(pot[i+1]-pot[i])
        drpsi[-1]=drpsi[-2]



        for i in range(len(dens)-1):
            dqi=0
            for j in np.arange(i,len(dens)-1):
                dqs=np.sqrt(pot[i]-pot[j+1])
                inte[i]+=drpsi[j]*(dqs-dqi)
                dqi=dqs
        inte[-1]=inte[-2]



        for i in range(len(dens)-1):
            df[i]=(inte[i+1]-inte[i])/(pot[i+1]-pot[i])




    idx=np.isfinite(df)
    pot=pot[idx]
    df=df[idx]
    df_func=UnivariateSpline(pot[::-1],df[::-1],k=1,s=0,ext=1) #Set always df=0, outside the grid




    return pot,df,df_func

