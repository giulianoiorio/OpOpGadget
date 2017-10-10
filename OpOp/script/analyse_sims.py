#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from OpOp.io import load_snap
from OpOp.analysis import Observe, Profile, Analysis
# import matplotlib as mpl
label_size = 20
# mpl.rcParams.update({'figure.autolayout':True})
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
mpl.rcParams['mathtext.default'] = 'regular'
import glob
import os
import sys
import importlib

class Param:

    def __init__(self, filename=None, default='long'):


        default_longrange={'proffile':None, 'folder':None, 'radmax':1.8, 'psun':(8.0,0.,0.),
                       'vsun':(-11.1, 12.24, 7.25), 'vrot':218, 'gpos':(4.913,-9.772,-85.387),
                       'gvel':(-36.774, 163.875, -96.074), 'skyposg':(287.534, -83.156), 'skypos':(15.0392, -33.7092),
                       'rh_obs': 0.283, 'vdisp_obs':8.69, 'vdisp_obs_tot':9.43, 'outdir': None, 'Nresample':100000,
                       'file_vdisp':'Scl_binned_profile_s3_rc.txt', 'dist':86, 'Vlos':111.05, 'pmot':(0.011, 0.143),
                       'mstar':4.6e6, }

        default_shortrange={'proffile':None, 'folder':None, 'radmax':1.2, 'psun':(8.0,0.,0.),
                       'vsun':(-11.1, 12.24, 7.25), 'vrot':218, 'gpos':(4.913,-9.772,-85.387),
                       'gvel':(-36.774, 163.875, -96.074), 'skyposg':(287.534, -83.156), 'skypos':(15.0392, -33.7092),
                       'rh_obs': 0.283, 'vdisp_obs':8.69, 'vdisp_obs_tot':9.40, 'outdir': None, 'Nresample':100000,
                       'file_vdisp':'Scl_binned_profile_s3_rc.txt', 'dist':86, 'Vlos':111.05, 'pmot':(0.011, 0.143),
                       'mstar':4.6e6, }

        default_exttrange={'proffile':None, 'folder':None, 'radmax':1.9, 'psun':(8.195,0.,0.),
                       'vsun':(-9.833, 11.945, 8.107), 'vrot':219.386, 'gpos':(4.881,-10.490,-91.658),
                       'gvel':(-40.843, 48.310, -81.362), 'skyposg':(287.534, -83.156), 'skypos':(15.0392, -33.7092),
                       'rh_obs': 0.303, 'vdisp_obs':8.69, 'vdisp_obs_tot':9.43, 'outdir': None, 'Nresample':100000,
                       'file_vdisp':'Scl_binned_profile_s3_rc.txt', 'dist':92.316, 'Vlos':110.714, 'pmot':(-0.05843, 0.39278),
                       'mstar':4.6e6, }

        default_exttranges={'proffile':None, 'folder':None, 'radmax':1.76, 'psun':(8.008,0.,0.),
                       'vsun':(-9.088, 10.286, 6.859), 'vrot':213.936, 'gpos':(4.983, -9.5487, -83.441),
                       'gvel':(-27.76599, 134.551590, -93.913965), 'skyposg':(287.534, -83.156), 'skypos':(15.0392, -33.7092),
                       'rh_obs': 0.276, 'vdisp_obs':8.69, 'vdisp_obs_tot':9.43, 'outdir': None, 'Nresample':100000,
                       'file_vdisp':'Scl_binned_profile_s3_rc.txt', 'dist':84.040, 'Vlos':110.914, 'pmot':(-0.0231058, 0.196979),
                       'mstar':4.6e6, }


        if default.lower()[0]=='l':
            self.default=default_longrange
        elif default.lower()[0]=='s':
            self.default=default_shortrange
        elif default.lower()[0]=='e':
            self.default=default_exttrange
        elif default.lower()[0]=='k':
            self.default=default_exttrange

        self.description={'proffile':'File with the velocity dispersion', 'folder':'?', 'radmax':'maximum radius to consider', 'psun':'Position of the Sun (X,Y,Z)',
                       'vsun':'Local velocty of the Sun (Vx, Vy, Vz)', 'vrot':'Velocity of LSR', 'gpos':'Galactic position of the object (Xg, Yg, Zg)',
                       'gvel':'Galactic velocity of the object (Xg, Yg, Zg)', 'skyposg': 'Position in sky coordinates (l [deg], b[deg])', 'skypos':'Position in equatorial coordinates (ra [deg], dec[deg])',
                       'rh_obs': 'Observed half light radius', 'vdisp_obs':'Observed velocity dispersion inside half-light radius', 'vdisp_obs_tot':'Observed velocity dispersion inside Rmax', 'outdir': 'Name of the output folder', 'Nresample':'DM particle to plot',
                       'file_vdisp':'File containing the velcoty dispersion profile', 'dist':'distance from the Sun ', 'Vlos':'Vlos wrt the Sun', 'pmot':'Proper motion (mul, mub) in mas/yr',
                       'mstar':'Stellar mass inside Rmax in solar masses', }

        self.used={}

        if filename is None:
            for key in self.default:
                #print(key)
                setattr(self, key, self.default[key])
                self.used[key]=self.default[key]
        else:


            directory, module_name = os.path.split(filename)
            module_name = os.path.splitext(module_name)[0]
            old_path = list(sys.path)
            if directory == '':
                sys.path.insert(0, os.getcwd())
            else:
                sys.path.insert(0, directory)

            try:
                fp = importlib.import_module(module_name)
            finally:
                sys.path[:] = old_path  # restore

            for key in self.default:
                try:
                    val=eval('fp.'+key)
                    setattr(self, key, val)
                    self.used[key] = val
                except AttributeError:
                    setattr(self, key, self.default[key])
                    self.used[key]=self.default[key]

    def save(self, filename='used_param'):

        fw=open(filename+'.py', 'w')

        for key in self.used:
            if self.used[key] is None:
                fw.write('%-20s = %-s #%s\n' % (key, str(self.used[key]), self.description[key]))
            elif key=='proffile' or key=='folder' or key=='outdir' or key=='file_vdisp':
                fw.write('%-20s = \'%-s\' #%s\n'%(key, str(self.used[key]), self.description[key]))
            else:
                fw.write('%-20s = %-s #%s\n' % (key, str(self.used[key]), self.description[key]))

        fw.close()


if len(sys.argv)>1:
    if sys.argv[1]=='-d' or sys.argv[1]=='-dl':
        par = Param()
        par.save('default_param')
        exit()
    elif sys.argv[1]=='-ds':
        par = Param(default='short')
        par.save('default_param')
        exit()
    elif sys.argv[1]=='-de':
        par = Param(default='ext')
        par.save('default_param')
        exit()
    elif sys.argv[1]=='-des':
        par = Param(default='k')
        par.save('default_param')
        exit()
    else:
        print(sys.argv)
        filename=sys.argv[1]
else:
    filename=None



#READ
par=Param(filename=filename)
proffile=par.proffile
folder=par.folder
radmax=par.radmax
psun=par.psun
vsun=par.vsun
vrot=par.vrot
gpos=par.gpos
gvel=par.gvel
skyposg=par.skyposg
skypos=par.skypos
rh_obs=par.rh_obs
vdisp_obs=par.vdisp_obs
vdisp_obs_tot=par.vdisp_obs_tot
#vdisp_obs=None
#vdisp_obs_tot=None

outdir=par.outdir
Nresample=par.Nresample
file_vdisp=par.file_vdisp
dist=par.dist
Vlos=par.Vlos
pmot=par.pmot
mstar=par.mstar



#START
if folder is None: folder = '.'
if outdir is None: outdir = './analysis'
if not os.path.exists(outdir):
    os.makedirs(outdir)

if file_vdisp is None:
    pass
elif not os.path.exists(file_vdisp):
    file_vdisp = None



olist = []
alist = []
proflist = []
simfiles = glob.glob(folder + '/*.bin')
simfiles.sort()
print('File read:', flush=True)
print(simfiles, flush=True)
Nfiles = len(simfiles)

figorbit, axarr = plt.subplots(3, 3, sharex=True, sharey=True)
figantd, axtd = plt.subplots(1, 3, tight_layout=True)
figanstd, axstd = plt.subplots(2, 3, tight_layout=True)
figobs, axobs = plt.subplots(2, 2, tight_layout=True, figsize=(10, 10))
figstat, axstat = plt.subplots(1, 1, tight_layout=True)
figev2D, axev2d = plt.subplots(2, 2, tight_layout=True, figsize=(10, 10))
figev3D, axev3d = plt.subplots(1, 2, tight_layout=True, figsize=(10, 5))
figvrot, axvrot = plt.subplots(1, 3, tight_layout=True, figsize=(15, 5))

colortd = ('blue', 'darkgreen', 'red')
check_td = 0

idx_plot_orbit = (int(Nfiles / 3.) - 1, int(2 * Nfiles / 3.) - 1, Nfiles - 1)
print(idx_plot_orbit)

log = ''
log += 'Setting:\n'
log += 'Xsun=%.3f Ysun=%.3f Zsun=%.3f\n' % psun
log += 'Vxsun=%.3f Vysun=%.3f Vzsun=%.3f Vrot=%.3f\n' % (vsun[0], vsun[1], vsun[2], vrot)
log += 'Distance=%.3f\n' % dist
if gpos is not None: log += 'Current position: Xg=%.3f Yg=%.3f Zg=%.3f\n' % gpos
if gvel is not None: log += 'Current velocity: Vxg=%.3f Vyg=%.3f Vzg=%.3f\n' % gvel
if skypos is not None: log += 'Current Sky position: ra=%.5f  dec=%.5f\n' % skypos
if skyposg is not None: log += 'Currrent Sky position2: l=%.5f b=%.5f\n' % skyposg
if Vlos is not None: log += 'Current Vlos: %.3f\n' % Vlos
if pmot is not None: log += 'Current pmot: mul=%.5f mub=%.5f\n' % pmot
if rh_obs is not None: log += 'Rh_obs:%.3f\n' % rh_obs
if vdisp_obs is not None: log += 'Vdisp_obs(R<Rh_obs):%.3f\n' % vdisp_obs
if vdisp_obs_tot is not None: log += 'Vdisp_obs(R<radmax):%.3f\n' % vdisp_obs_tot
if radmax is not None: log += 'Rmax=%.3f\n' % radmax
if mstar is not None: log += 'M*(R<Rmax)=%.3e\n' % mstar
log += '\n'

#save used param
par.save(filename=outdir+'/used_param')



fin_array = np.zeros((len(simfiles), 11))
fin_array_header = '0-T[Gyr] 1-Xg[kpc] 2-Yg[kpc] 3-Zg[kpc] 4-Vxg[kpc] 5-Vyg[kpc] 6-Vzg[kpc] 7-Rh[kpc] 8-Vdisp(R<Rh)[km/s] 9-Mass(R<Rmax) 10-Vdisp(R<Rmax)[km/s]'

fin_array_intrinsic = np.zeros((len(simfiles), 19))
fin_array_intrinsic_header = '0-T[Gyr] 1-rg[kpc] 2-Vg[km/s] 3-Rh_3d[kpc] 4-Rh_x[kpc] 5-Rh_y[kpc] 6-Rh_z[kpc] 7-Vdisp_x(R<Rhx)[kpc] 8-Vdisp_y(R<Rhy)[kpc] 9-Vdisp_z(R<Rhz)[kpc] 10-Vdisp_x(R<radmax)[kpc] 11-Vdisp_y(R<radmax)[kpc] 12-Vdisp_z(R<radmax)[kpc] 13-Mass(r<radmax) [Msun], 14-Mass_x (R<radmax) [Msun], 15-Mass_y (R<radmax) [Msun], 16-Mass_x (R<radmax) [Msun], 17-p, 18-q'

i = 0
iplot = 0
for file in simfiles:
    p_tmp = load_snap(file)
    time_tmp = p_tmp.header['Time']
    fin_array[i, 0] = time_tmp
    fin_array_intrinsic[i, 0] = time_tmp

    log += '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
    log += 'Time= %.3f\n' % time_tmp

    # Plot orbit
    if (i == idx_plot_orbit[0]) or (i == idx_plot_orbit[1]) or (i == idx_plot_orbit[2]):
        idh = p_tmp.Type == 1
        ids = ~idh
        Nhalos = np.sum(idh)
        Nstars = np.sum(ids)
        print('Plottng', flush=True)
        if Nhalos > 0:
            hp = p_tmp.Pos[idh]
            if (Nresample is not None) and (Nresample < len(hp)):
                idh = np.random.choice(len(hp), Nresample, replace=False)
            else:
                idh = np.array([True, ] * len(hp))
            axarr[iplot, 0].scatter(hp[idh, 0], hp[idh, 1], s=0.00005, c='black')
            axarr[iplot, 1].scatter(hp[idh, 0], hp[idh, 2], s=0.00005, c='black')
            axarr[iplot, 2].scatter(hp[idh, 1], hp[idh, 2], s=0.00005, c='black')

        if Nstars > 0:
            sp = p_tmp.Pos[ids]
            axarr[iplot, 0].scatter(sp[:, 0], sp[:, 1], s=0.00005, c='red', zorder=1000)
            axarr[iplot, 1].scatter(sp[:, 0], sp[:, 2], s=0.00005, c='red', zorder=1000)
            axarr[iplot, 2].scatter(sp[:, 1], sp[:, 2], s=0.00005, c='red', zorder=1000)
            axarr[iplot, 0].set_xlim(-150, 150)

        axarr[iplot, 0].set_ylim(-150, 150)
        axarr[iplot, 1].set_xlim(-150, 150)
        axarr[iplot, 1].set_ylim(-150, 150)
        axarr[iplot, 2].set_xlim(-150, 150)
        axarr[iplot, 2].set_ylim(-150, 150)
        axarr[iplot, 0].text(-150, 130, 'T=%.2f Gyr' % time_tmp, fontsize=20)
        if i == idx_plot_orbit[0]:
            axarr[0, 0].text(0, 150, 'XY', fontsize=25)
            axarr[0, 1].text(0, 150, 'XZ', fontsize=25)
            axarr[0, 2].text(0, 150, 'YZ', fontsize=25)

        if gpos is not None:
            axarr[iplot, 0].plot([gpos[0], gpos[0]], [-1000, 1000], color='cyan', zorder=2000, lw=0.5)
            axarr[iplot, 0].plot([-1000, 1000], [gpos[1], gpos[1]], color='cyan', zorder=2000, lw=0.5)
            axarr[iplot, 1].plot([gpos[0], gpos[0]], [-1000, 1000], color='cyan', zorder=2000, lw=0.5)
            axarr[iplot, 1].plot([-1000, 1000], [gpos[2], gpos[2]], color='cyan', zorder=2000, lw=0.5)
            axarr[iplot, 2].plot([gpos[1], gpos[1]], [-1000, 1000], color='cyan', zorder=2000, lw=0.5)
            axarr[iplot, 2].plot([-1000, 1000], [gpos[2], gpos[2]], color='cyan', zorder=2000, lw=0.5)

        if (i == idx_plot_orbit[2]):
            axarr[2, 0].set_xlabel('kpc', fontsize=20)
            axarr[2, 1].set_xlabel('kpc', fontsize=20)
            axarr[2, 2].set_xlabel('kpc', fontsize=20)
            axarr[0, 0].set_ylabel('kpc', fontsize=20)
            axarr[1, 0].set_ylabel('kpc', fontsize=20)
            axarr[2, 0].set_ylabel('kpc', fontsize=20)
            figorbit.set_size_inches(15, 15, forward=True)
            # plt.setp([a.get_xticklabels() for a in figorbit.axes[:-1]], visible=False)
            figorbit.savefig(outdir + '/orbit.png')
            del axarr
            del figorbit
        iplot += 1

        print('Done', flush=True)

    # Plot analysis
    print('Find COM', flush=True)
    a_tmp = Analysis(particles=p_tmp, safe=True, auto_centre=True, iter=True, single=False)
    com_tmp = a_tmp.pcom
    vcom_tmp = a_tmp.pvcom
    log += 'COM: X=%.3f Y=%.3f Z=%.3f\n' % tuple(com_tmp)
    log += 'VCOM: VX=%.3f VY=%.3f VZ=%.3f\n' % tuple(vcom_tmp)
    print('COM', com_tmp, flush=True)
    print('VCOM', vcom_tmp, flush=True)
    print('Done', flush=True)

    # Some calculations
    if radmax is None: radmax = 100

    # 3D
    rh_3d = a_tmp.qmass(q=50, rad_max=radmax, type=2)
    mass_3d = a_tmp.mass(rad=radmax, type=2)
    # Shape
    evalue = a_tmp.inertia_tensor(eig=True, maxrad=radmax, type=2)[1]
    ytox_ratio = evalue[1]
    ztox_ratio = evalue[0]

    # Proj-X
    rh_px = a_tmp.qmassup(q=50, pax='x', rad_max=radmax, type=2)
    rh_tmp_x = rh_px
    vdisp_rh_px = a_tmp.vdisp(rh_px, pax='x', type=2)
    vdisp_rh_px_tot = a_tmp.vdisp(radmax, pax='x', type=2)
    mass_px = a_tmp.massup(radmax, pax='x', type=2)

    # Proj-X
    rh_py = a_tmp.qmassup(q=50, pax='y', rad_max=radmax, type=2)
    rh_tmp_y = rh_py
    vdisp_rh_py = a_tmp.vdisp(rh_py, pax='y', type=2)
    vdisp_rh_py_tot = a_tmp.vdisp(radmax, pax='y', type=2)
    mass_py = a_tmp.massup(radmax, pax='y', type=2)

    # Proj-X
    rh_pz = a_tmp.qmassup(q=50, pax='z', rad_max=radmax, type=2)
    rh_tmp_z = rh_pz
    vdisp_rh_pz = a_tmp.vdisp(rh_pz, pax='z', type=2)
    vdisp_rh_pz_tot = a_tmp.vdisp(radmax, pax='z', type=2)
    mass_pz = a_tmp.massup(radmax, pax='z', type=2)

    rg_tmp = np.sum(np.sqrt(com_tmp * com_tmp))
    vg_tmp = np.sum(np.sqrt(vcom_tmp * vcom_tmp))

    fin_array_intrinsic[i, 1] = rg_tmp
    fin_array_intrinsic[i, 2] = vg_tmp
    fin_array_intrinsic[i, 3] = rh_3d
    fin_array_intrinsic[i, 4] = rh_px
    fin_array_intrinsic[i, 5] = rh_py
    fin_array_intrinsic[i, 6] = rh_pz
    fin_array_intrinsic[i, 7] = vdisp_rh_px
    fin_array_intrinsic[i, 8] = vdisp_rh_py
    fin_array_intrinsic[i, 9] = vdisp_rh_pz
    fin_array_intrinsic[i, 10] = vdisp_rh_px_tot
    fin_array_intrinsic[i, 11] = vdisp_rh_py_tot
    fin_array_intrinsic[i, 12] = vdisp_rh_pz_tot
    fin_array_intrinsic[i, 13] = mass_3d
    fin_array_intrinsic[i, 14] = mass_px
    fin_array_intrinsic[i, 15] = mass_py
    fin_array_intrinsic[i, 16] = mass_pz
    fin_array_intrinsic[i, 17] = ytox_ratio
    fin_array_intrinsic[i, 18] = ztox_ratio

    if (i == 0) or (i == idx_plot_orbit[1]) or (i == idx_plot_orbit[2]):
        # dens
        try:
            prof_tmp_h = Profile(particles=a_tmp.p, xmin=0.01, xmax=100, ngrid=100, kind='log', type=1)
            arr = prof_tmp_h.dens()[0]
            r = arr[:, 0]
            d = arr[:, 1]
            axtd[0].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
        except:
            pass

        prof_tmp_s = Profile(particles=a_tmp.p, xmin=0.0001, xmax=5, ngrid=100, kind='lin', type=2)
        arr = prof_tmp_s.dens()[0]
        r = arr[:, 0]
        d = arr[:, 1]
        axtd[1].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
        arr = prof_tmp_s.vdisp3d()[0]
        r = arr[:, 0]
        vd = arr[:, 1]
        axtd[2].plot(r, vd, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])

        arr = prof_tmp_s.supdens(pax='z')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        dminS0 = np.min(d)
        dmaxS0 = np.max(d)
        axstd[0, 0].plot([rh_tmp_z, rh_tmp_z], [np.min(d), np.max(d)], color=colortd[check_td])
        axstd[0, 0].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])

        arr = prof_tmp_s.supdens(pax='y')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        dminS1 = np.min(d)
        dmaxS1 = np.max(d)
        axstd[0, 1].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
        axstd[0, 1].plot([rh_tmp_y, rh_tmp_y], [np.min(d), np.max(d)], color=colortd[check_td])

        arr = prof_tmp_s.supdens(pax='x')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        dminS2 = np.min(d)
        dmaxS2 = np.max(d)
        axstd[0, 2].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
        axstd[0, 2].plot([rh_tmp_x, rh_tmp_x], [np.min(d), np.max(d)], color=colortd[check_td])  # vdisp

        arr = prof_tmp_s.vdisp2d(pax='z')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        axstd[1, 0].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
        arr = prof_tmp_s.vdisp2d(pax='y')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        axstd[1, 1].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
        arr = prof_tmp_s.vdisp2d(pax='x')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        axstd[1, 2].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])

        #Vsys
        arr = prof_tmp_s.vsys(pax='z')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        vsysmax=[np.max(d),]
        vsysmin=[np.min(d),]
        axvrot[0].plot(r, d, label='T=%.2f Gyr'%time_tmp, color=colortd[check_td])
        arr = prof_tmp_s.vsys(pax='y')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        vsysmax.append(np.max(d))
        vsysmin.append(np.min(d))
        axvrot[1].plot(r, d, label='T=%.2f Gyr'%time_tmp, color=colortd[check_td])
        arr = prof_tmp_s.vsys(pax='x')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        vsysmax.append(np.max(d))
        vsysmin.append(np.min(d))
        axvrot[2].plot(r, d, label='T=%.2f Gyr'%time_tmp, color=colortd[check_td])
        axvrot[0].axhline(0, color='black', ls='--')
        axvrot[1].axhline(0, color='black', ls='--')
        axvrot[2].axhline(0, color='black', ls='--')

        check_td += 1

        if check_td == 3:

            #Vsys
            axvrot[0].set_xlabel('$R=\\sqrt{X^2+Y^2} [kpc]$', fontsize=20)
            axvrot[1].set_xlabel('$R=\\sqrt{X^2+Z^2}  [kpc]$', fontsize=20)
            axvrot[2].set_xlabel('$R=\\sqrt{Y^2+Z^2}  [kpc]$', fontsize=20)
            axvrot[0].set_ylabel('$V_{sys} [km/s]$', fontsize=20)
            axvrot[0].set_title('Z-Projection', fontsize=20)
            axvrot[1].set_title('Y-Projection', fontsize=20)
            axvrot[2].set_title('X-Projection', fontsize=20)
            figvrot.savefig(outdir+'/Vsys.png')
            #axvrot_min=np.min(vsysmin)
            #axvrot_max=np.max(vsysmax)
            #axvrot[0].set_xlim(axvrot_min*1.5,ax)


            # 3D
            axtd[0].set_xlabel('r [kpc]', fontsize=20)
            axtd[1].set_xlabel('r [kpc]', fontsize=20)
            axtd[2].set_xlabel('r [kpc]', fontsize=20)
            axtd[0].set_ylabel('$\\rho_h \\ [M_\\odot / kpc^3]$', fontsize=20)
            axtd[1].set_ylabel('$\\rho_* \\ [M_\\odot / kpc^3]$', fontsize=20)
            axtd[2].set_ylabel('$\\sigma \\ [km/s]$', fontsize=20)
            axtd[0].set_xlim(0.1, 100)
            axtd[0].set_ylim(1e-2, 1e10)
            axtd[0].set_xscale('log')
            axtd[0].set_yscale('log')
            axtd[1].set_xlim(0.01, 10)
            axtd[1].set_xscale('log')
            axtd[1].set_yscale('log')
            axtd[2].set_xlim(0, 2)
            axtd[0].legend(loc='upper right')
            figantd.set_size_inches(15, 5, forward=True)
            figantd.savefig(outdir + '/3Danalysis.png')
            del axtd
            del figantd

            # 2D
            if rh_obs is not None:
                axstd[0, 0].plot([rh_obs, rh_obs], [dminS0, dmaxS0], '--', color='black', zorder=5000,
                                 label='$R^{obs}_h$')
                axstd[0, 1].plot([rh_obs, rh_obs], [dminS1, dmaxS1], '--', color='black', zorder=5000)
                axstd[0, 2].plot([rh_obs, rh_obs], [dminS2, dmaxS2], '--', color='black', zorder=5000)
            if vdisp_obs is not None:
                axstd[1, 0].plot([0, 2], [vdisp_obs, vdisp_obs], '--', color='black', zorder=5000,
                                 label='$\\sigma^{los}(R<R_h)$')
                axstd[1, 1].plot([0, 2], [vdisp_obs, vdisp_obs], '--', color='black', zorder=5000)
                axstd[1, 2].plot([0, 2], [vdisp_obs, vdisp_obs], '--', color='black', zorder=5000)
            if file_vdisp is not None:
                data = np.loadtxt(file_vdisp)
                x = dist * np.tan(data[:, 0] * (np.pi) / 180)
                ex = dist * np.tan(data[:, 1] * (np.pi) / 180)
                axstd[1, 0].errorbar(x, data[:, 4], data[:, 5], ex, fmt='o', c='black', label='Observed', zorder=8000)
                axstd[1, 1].errorbar(x, data[:, 4], data[:, 5], ex, fmt='o', c='black', label='Observed', zorder=8000)
                axstd[1, 2].errorbar(x, data[:, 4], data[:, 5], ex, fmt='o', c='black', label='Observed', zorder=8000)

            axstd[1, 0].set_xlabel('$R=\\sqrt{X^2+Y^2} [kpc]$', fontsize=20)
            axstd[1, 1].set_xlabel('$R=\\sqrt{X^2+Z^2}  [kpc]$', fontsize=20)
            axstd[1, 2].set_xlabel('$R=\\sqrt{Y^2+Z^2}  [kpc]$', fontsize=20)
            axstd[0, 0].set_xscale('log')
            axstd[0, 0].set_yscale('log')
            axstd[0, 1].set_xscale('log')
            axstd[0, 1].set_yscale('log')
            axstd[0, 2].set_xscale('log')
            axstd[0, 2].set_yscale('log')
            axstd[1, 0].set_xlim(0, 2)
            axstd[1, 1].set_xlim(0, 2)
            axstd[1, 2].set_xlim(0, 2)
            axstd[1, 0].set_ylim(0, 20)
            axstd[1, 1].set_ylim(0, 20)
            axstd[1, 2].set_ylim(0, 20)
            axstd[1, 0].set_ylabel('$\\sigma \\ [km/s]$', fontsize=20)
            axstd[0, 0].set_ylabel('$\\Sigma \\ [M_\\odot / kpc^2]$', fontsize=20)
            axstd[0, 0].legend(loc='upper right')
            axstd[1, 0].legend(loc='upper left')
            # axstd[0, 0].text(0.5, 1.0, 'Projection ax: X', fontsize=25, transform=axstd[0, 0].transAxes)
            # axstd[0, 1].text(0.5, 1.0, 'Projection ax: Y', fontsize=25, transform=axstd[0, 1].transAxes)
            # axstd[0, 2].text(0.5, 1.0, 'Projection ax: Z', fontsize=25, transform=axstd[0, 2].transAxes)
            figanstd.set_size_inches(15, 10, forward=True)
            figanstd.savefig(outdir + '/2Danalysis.png')

    # observe
    o_tmp = Observe(particles=p_tmp, type=2)
    s, c = o_tmp.observe(psun=psun, vsun=vsun, vrot=vrot, mq=50)
    a = Analysis(s, safe=True, auto_centre=False, iter=False, single=False)

    if radmax is None:
        rh_sim = a.qmassup(q=50, pax='obs')
        vdisp_sim = a.vdisp(rh_sim, pax='obs')
        vdisp_sim_tot = a.vdisp(100, pax='obs')
        mass_sim = a.massup(rad=100, pax='obs')
    else:
        rh_sim = a.qmassup(q=50, pax='obs', rad_max=radmax)
        vdisp_sim = a.vdisp(rh_sim, pax='obs')
        vdisp_sim_tot = a.vdisp(radmax, pax='obs')
        mass_sim = a.massup(radmax, pax='obs')

    fin_array[i, 1:4] = c.Pos[:]
    fin_array[i, 4:7] = c.Vel[:]
    fin_array[i, 7] = rh_sim
    fin_array[i, 8] = vdisp_sim
    fin_array[i, 9] = mass_sim
    fin_array[i, 10] = vdisp_sim_tot

    log += '--Observed--\n'
    log += 'Rh_sim=%.3f\n' % rh_sim
    log += 'Vdisp(R<Rh_sim)=%.3f\n' % vdisp_sim
    log += 'Vdisp(R<Rad_max)=%.3f\n' % vdisp_sim_tot
    log += 'Mass(R<Rmax)=%.3e\n' % mass_sim
    log += 'Stellar centre:\n'
    log += c.__str__()
    log += '\n'
    log += '\n \n'

    if i == idx_plot_orbit[2]:

        prof_obs = Profile(particles=s, xmin=0.00001, xmax=10, ngrid=100, kind='lin')

        axobs[0,0].scatter(s.xi[:] / 3600., s.eta[:] / 3600., s=0.005, c='red')
        axobs[0,0].set_xlabel('$\\xi [deg]$', fontsize=20)
        axobs[0,0].set_ylabel('$\\eta [deg]$', fontsize=20)
        axobs[0,0].set_xlim(-2, 2)
        axobs[0,0].set_ylim(-2, 2)

        arr = prof_obs.vdisp2d(pax='obs')[0]
        r = arr[:, 0]
        vd = arr[:, 1]
        axobs[1,0].plot(r, vd, label='Vlos', lw=3, color='red', zorder=2000)

        if file_vdisp is not None:
            try:
                data = np.loadtxt(file_vdisp)
                x = dist * np.tan(data[:, 0] * (np.pi) / 180)
                ex = dist * np.tan(data[:, 1] * (np.pi) / 180)
                axobs[1,0].errorbar(x, data[:, 4], data[:, 5], ex, fmt='o', c='black', label='Observed', zorder=1000)

                obins = np.zeros(len(data) + 1)
                obins[:-1] = data[:, 6]
                obins[-1] = data[-1, 7]
                Nperbin = data[:, 8]
                obins = dist * np.tan(obins * np.pi / 180.)

                color_disp = ('blue', 'orange', 'magenta')

                for j in range(3):
                    b = a.binned_dispersion(bins=obins, pax='obs', Nperbin=Nperbin, bins_kind='lin', velocity_err=None,
                                            err_distibution='uniform', nboot=10000)

                    axobs[1,0].errorbar(b[0], b[4], b[5], b[1], fmt='o', c=color_disp[j], mfc='white')
            except FileNotFoundError:
                print('File %s not found.. skipping' % file_vdisp)

        axobs[1,0].set_xlabel('$R [kpc]$', fontsize=20)
        axobs[1,0].set_ylabel('$\\sigma_{los} [km/s]$', fontsize=20)
        axobs[1,0].set_xlim(0, 2)
        axobs[1,0].set_ylim(0, 15)

        # densup
        arr = prof_obs.supdens(pax='obs')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        axobs[0,1].plot(r, d, lw=3, color='red')
        axobs[0,1].plot([rh_sim, rh_sim], [np.min(d), np.max(d)], color='red', label='$R^{sim}_h$')
        if rh_obs is not None:
            axobs[0,1].plot([rh_obs, rh_obs], [np.min(d), np.max(d)], '--', color='black', label='$R^{obs}_h$')
        axobs[0,1].set_xlabel('$R [kpc]$', fontsize=20)
        axobs[0,1].set_ylabel('$\\Sigma_{los} [M_\\odot/kpc^2]$', fontsize=20)
        # axobs[1].set_xlim(0.001,10)
        axobs[0,1].set_xscale('log')
        axobs[0,1].set_yscale('log')
        axobs[0,1].legend(loc='upper right')

        #Vsys
        arr = prof_obs.vsys(pax='obs')[0]
        r = arr[:, 0]
        d = arr[:, 1]
        axobs[1,1].plot(r, d, lw=3, color='red')
        if Vlos is not None:
            axobs[1,1].axhline(Vlos, color='black', ls='--', label='Observed')
        else:
            medVlos=np.median(d)
            axobs[1,1].axhline(medVlos, color='black', ls='--', label='Median')
        axobs[1,1].legend(loc='best')
        axobs[1,1].set_xlabel('$R [kpc]$', fontsize=20)
        axobs[1,1].set_ylabel('$V_{sys} [km/s]$', fontsize=20)
        axobs[1,1].set_xlim(0,2)

        #figobs.set_size_inches(15, 5, forward=True)
        figobs.savefig(outdir + '/Obs_analysis.png')
    # a_tmp=Analysis(p_tmp,safe=True, auto_centre=True, iter=True, single=False)
    # o_tmp=Observe(particles=p_tmp, type=2)
    # s_tmp,c_tmp=o_tmp.observe(psun=psun, vsun=vsun, vrot=vrot, mq=50)
    i += 1

T = fin_array[1:, 0]
rhp = fin_array[1:, 7]
if rh_obs is None:
    rhp_ini = rhp[0]
    label = '$R^{sim}_h$ wrt T=0'
else:
    rhp_ini = rh_obs
    label = '$R^{sim}_h$ wrt obs'
rel_diff = 100 * (rhp - rhp_ini) / rhp_ini
axstat.plot(T, rel_diff, '-o', c='red', label=label, zorder=3000)

vdp = fin_array[1:, 8]
if vdisp_obs is None:
    vdp_ini = vdp[0]
    label = '$\\sigma^{sim}_{los}(R<R_h)$ wrt T=0'
else:
    vdp_ini = vdisp_obs
    label = '$\\sigma^{sim}_{los}(R<R_h)$ wrt obs'
rel_diff = 100 * (vdp - vdp_ini) / vdp_ini
axstat.plot(T, rel_diff, '-s', c='blue', label=label, zorder=3000)

mp = fin_array[1:, 9]
mp_ini = mp[0]
if mstar is None:
    mp_ini = mp[0]
    label = '$M^{sim}_*$ wrt T=0'
else:
    mp_ini = mstar
    label = '$M^{sim}_*$ wrt obs'
rel_diff = 100 * (mp - mp_ini) / mp_ini
axstat.plot(T, rel_diff, '-^', c='darkgreen', label=label, zorder=3000)

vdp_tot = fin_array[1:, 10]
if vdisp_obs_tot is None:
    vdp_tot_ini = vdp_tot[0]
    label = '$\\sigma^{sim}_{los}(R<R_{max})$ wrt T=0'
else:
    vdp_tot_ini = vdisp_obs_tot
    label = '$\\sigma^{sim}_{los}(R<R_{max})$ wrt obs'
rel_diff = 100 * (vdp_tot - vdp_tot_ini) / vdp_ini
axstat.plot(T, rel_diff, '-s', c='cyan', label=label, zorder=3000)
axstat.legend(loc='best', fontsize=12)

axstat.axhline(0, color='black', ls='--')
axstat.set_xlabel('T [Gyr]', fontsize=20)
axstat.set_ylabel('Rel diff [%]', fontsize=20)

figstat.savefig(outdir + '/evolution_los.png')

# EV2D
# Vdisp
T2 = fin_array_intrinsic[:, 0]
vv_x = fin_array_intrinsic[:, 7]
vv_x0 = vv_x[0]
vv_x = 100 * (vv_x - vv_x0) / (vv_x0)
vv_y = fin_array_intrinsic[:, 8]
vv_y0 = vv_y[0]
vv_y = 100 * (vv_y - vv_y0) / (vv_y0)
vv_z = fin_array_intrinsic[:, 9]
vv_z0 = vv_z[0]
vv_z = 100 * (vv_z - vv_z0) / (vv_z0)
vv_los = fin_array[1:, 8]
vv_los0 = vv_los[0]
vv_los = 100 * (vv_los - vv_los0) / (vv_los0)

axev2d[0, 0].plot(T2, vv_x, '-o', label='X-projection', lw=3, markersize=10)
axev2d[0, 0].plot(T2, vv_y, '-^', label='Y-projection', lw=3, markersize=10)
axev2d[0, 0].plot(T2, vv_z, '-s', label='Z-projection', lw=3, markersize=10)
axev2d[0, 0].plot(T, vv_los, '-d', label='los-projection', lw=3, markersize=10)
axev2d[0, 0].axhline(0, color='black', ls='--')

# Vdisp_tot
T2 = fin_array_intrinsic[:, 0]
vv_x = fin_array_intrinsic[:, 10]
vv_x0 = vv_x[0]
vv_x = 100 * (vv_x - vv_x0) / (vv_x0)
vv_y = fin_array_intrinsic[:, 11]
vv_y0 = vv_y[0]
vv_y = 100 * (vv_y - vv_y0) / (vv_y0)
vv_z = fin_array_intrinsic[:, 12]
vv_z0 = vv_z[0]
vv_z = 100 * (vv_z - vv_z0) / (vv_z0)
vv_los = fin_array[1:, 10]
vv_los0 = vv_los[0]
vv_los = 100 * (vv_los - vv_los0) / (vv_los0)

axev2d[0, 1].plot(T2, vv_x, '-o', label='X-projection', lw=3, markersize=10)
axev2d[0, 1].plot(T2, vv_y, '-^', label='Y-projection', lw=3, markersize=10)
axev2d[0, 1].plot(T2, vv_z, '-s', label='Z-projection', lw=3, markersize=10)
axev2d[0, 1].plot(T, vv_los, '-d', label='los-projection', lw=3, markersize=10)
axev2d[0, 1].axhline(0, color='black', ls='--')

# Rh
T2 = fin_array_intrinsic[:, 0]
vv_x = fin_array_intrinsic[:, 4]
vv_x0 = vv_x[0]
vv_x = 100 * (vv_x - vv_x0) / (vv_x0)
vv_y = fin_array_intrinsic[:, 5]
vv_y0 = vv_y[0]
vv_y = 100 * (vv_y - vv_y0) / (vv_y0)
vv_z = fin_array_intrinsic[:, 6]
vv_z0 = vv_z[0]
vv_z = 100 * (vv_z - vv_z0) / (vv_z0)
vv_los = fin_array[1:, 7]
vv_los0 = vv_los[0]
vv_los = 100 * (vv_los - vv_los0) / (vv_los0)

axev2d[1, 0].plot(T2, vv_x, '-o', label='X-projection', lw=3, markersize=10)
axev2d[1, 0].plot(T2, vv_y, '-^', label='Y-projection', lw=3, markersize=10)
axev2d[1, 0].plot(T2, vv_z, '-s', label='Z-projection', lw=3, markersize=10)
axev2d[1, 0].plot(T, vv_los, '-d', label='los-projection', lw=3, markersize=10)
axev2d[1, 0].axhline(0, color='black', ls='--')

# Mass
T2 = fin_array_intrinsic[:, 0]
vv_x = fin_array_intrinsic[:, 14]
vv_x0 = vv_x[0]
vv_x = 100 * (vv_x - vv_x0) / (vv_x0)
vv_y = fin_array_intrinsic[:, 15]
vv_y0 = vv_y[0]
vv_y = 100 * (vv_y - vv_y0) / (vv_y0)
vv_z = fin_array_intrinsic[:, 16]
vv_z0 = vv_z[0]
vv_z = 100 * (vv_z - vv_z0) / (vv_z0)
vv_los = fin_array[1:, 9]
vv_los0 = vv_los[0]
vv_los = 100 * (vv_los - vv_los0) / (vv_los0)

axev2d[1, 1].plot(T2, vv_x, '-o', label='X-projection', lw=3, markersize=10)
axev2d[1, 1].plot(T2, vv_y, '-^', label='Y-projection', lw=3, markersize=10)
axev2d[1, 1].plot(T2, vv_z, '-s', label='Z-projection', lw=3, markersize=10)
axev2d[1, 1].plot(T, vv_los, '-d', label='los-projection', lw=3, markersize=10)
axev2d[1, 1].axhline(0, color='black', ls='--')

axev2d[0, 0].set_title('$\\sigma(R<R_h)$', fontsize=20)
axev2d[0, 1].set_title('$\\sigma(R<R_{max})$', fontsize=20)
axev2d[1, 0].set_title('$R_h$', fontsize=20)
axev2d[1, 1].set_title('$Mass(R<R_{max})$', fontsize=20)

axev2d[0, 0].legend(loc='best', fontsize=12)
axev2d[1, 0].set_xlabel('T [Gyr]', fontsize=20)
axev2d[1, 1].set_xlabel('T [Gyr]', fontsize=20)
axev2d[0, 0].set_ylabel('Rel diff [%]', fontsize=20)
axev2d[1, 0].set_ylabel('Rel diff [%]', fontsize=20)
axev2d[0, 0].set_ylim(-30, 30)
axev2d[0, 1].set_ylim(-30, 30)
axev2d[1, 0].set_ylim(-30, 30)
axev2d[1, 1].set_ylim(-30, 30)

# figev2D.set_size_inches(20, 15, forward=True)
figev2D.savefig(outdir + '/evolution_2D.png')

# EV3D
pf = fin_array_intrinsic[:, 17]
qf = fin_array_intrinsic[:, 18]
axev3d[0].plot(T2, pf, '-o', label='Y-to-X axial ratio', lw=3, markersize=10)
axev3d[0].plot(T2, qf, '-s', label='Z-to-X axial ratio', lw=3, markersize=10)

vv_x = fin_array_intrinsic[:, 3]
vv_x0 = vv_x[0]
vv_x = 100 * (vv_x - vv_x0) / (vv_x0)
vv_y = fin_array_intrinsic[:, 13]
vv_y0 = vv_y[0]
vv_y = 100 * (vv_y - vv_y0) / (vv_y0)

axev3d[1].plot(T2, vv_x, '-o', label='$r_h$', lw=3, markersize=10, zorder=3000)
axev3d[1].plot(T2, vv_y, '-s', label='$Mass(r<r_{max}}$', lw=3, markersize=10)

axev3d[0].legend(loc='best', fontsize=12)
axev3d[1].legend(loc='best', fontsize=12)
axev3d[1].set_xlabel('T [Gyr]', fontsize=20)
axev3d[0].set_xlabel('T [Gyr]', fontsize=20)
axev3d[0].set_ylabel('Axial ratio', fontsize=20)
axev3d[1].set_ylabel('Rel diff [%]', fontsize=20)

figev3D.savefig(outdir + '/evolution_3D.png')

np.savetxt(outdir + '/data_obs.txt', fin_array, fmt='%.3e', header=fin_array_header)
np.savetxt(outdir + '/data_intrinsic.txt', fin_array_intrinsic, fmt='%.3e', header=fin_array_intrinsic_header)
file = open(outdir + '/info.log', mode='w')
file.write(log)
file.close()
