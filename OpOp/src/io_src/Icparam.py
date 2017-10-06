from __future__ import  division, print_function
import numpy as np

def write_icparam_fvfps(id_new=0, tmax=1, max_dt=0.001, poscm=(0,0,0), velcm=(0,0,0), theta=0.5, epsc=0.01, isoft=1, noutput=1, iene=120, new=0, epotential='J95', outfolder='.' **kwargs):
    """
    see fvfps documentation
    :param id_new:
    :param tmax:
    :param max_dt:
    :param poscm:
    :param velcm:
    :param theta:
    :param epsc:
    :param isoft:
    :param noutput:
    :param iene:
    :param new:
    :param epotential:
    :param kwargs:  lscale, mscale, vscale, tscale
    :return:
    """
    #tmax Gy
    #max_dt Gy

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

    #if id_new<0: raise ValueError('id_new lower than 0 not allowed (write_icparam)')
    #elif id_new<10:  outname='out0'+str(id_new)
    #else: outname='out'+str(id_new)

    outname=outfolder+'./input.data'



    fo=open(outname,'w')

    print('theta,  epsc, isoft',file=fo)
    print('%4.3f,  %4.3f,  %i'%(theta,epsc,isoft),file=fo)
    print('itest,ibound,ivir,pout',file=fo)
    print('0,    0,     0,   0',file=fo)
    print('nout, iene', file=fo)
    print('%4i, %4i'%(noutput,iene), file=fo)
    print('imd, ctstep ', file=fo)
    print('100,  2', file=fo)
    print('tmax,  max_dt ', file=fo)
    print('%5.4f, %5.4f'%(tmax/tscale, max_dt/tscale), file=fo)
    print('new,  id_new ', file=fo)
    print('%i, %i'%(new,id_new), file=fo)

    #epot
    if epotential is None:
        print('iext,  mdisk, mspher, vhalo !No external potential', file=fo)
        print('0,  0,   0,   0', file=fo)
        print('adisk,  bdisk,  cspher, dhalo', file=fo)
        print('0,  0,   0,   0', file=fo)
    elif epotential=='J95':
        print('iext,  mdisk, mspher, vhalo !J95', file=fo)
        print('1,     10.,   3.4,    128.', file=fo)
        print('adisk,  bdisk,  cspher, dhalo', file=fo)
        print('6.5,	0.26,	0.7,	12.', file=fo)
    else:
        raise NotImplementedError('External potential %s not implemented'%epotential)

    #cm
    print('xcm, ycm, zcm, vxcm, vycm, vzcm',file=fo)
    print('%4.3f, %4.3f, %4.3f, %4.3f, %4.3f, %4.3f'%(poscm[0],poscm[1],poscm[2],velcm[0],velcm[1],velcm[2]),file=fo)

    fo.close()

