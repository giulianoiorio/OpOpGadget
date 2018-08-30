#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from OpOp.io import load_snap
from OpOp.analysis import Observe, Profile, Analysis
from OpOp.utility import radec_to_xieta
# import matplotlib as mpl
label_size = 20
# mpl.rcParams.update({'figure.autolayout':True})
mpl.rcParams['xtick.labelsize'] = label_size
mpl.rcParams['ytick.labelsize'] = label_size
mpl.rcParams['mathtext.default'] = 'it'
import glob
import os
import sys
import importlib

class Param:

	def __init__(self, filename=None, default='h'):


		default_helmi={'proffile':None, 'folder':None, 'radmax':1.5, 'psun':(8.0,0.,0.),
					   'vsun':(-11.1, 12.24, 7.25), 'vrot':218, 'gpos':(4.913,-9.772,-85.387),
					   'gvel':(-17.15, 155.07, -95.78), 'skyposg':(287.534, -83.156), 'skypos':(15.0392, -33.7092),
					   'rh_obs': 0.283, 'vdisp_obs':8.4, 'vdisp_obs_tot':8.94, 'outdir': None, 'Nresample':100000,
					   'file_vdisp':None, 'dist':86, 'Vlos':111.05, 'pmot':(0.082, -0.131),
					   'mstar':4.6e6, 'comtype': 2, 'orbit_track': None, 'Sigma_scale':1 }

		default_longrange={'proffile':None, 'folder':None, 'radmax':1.8, 'psun':(8.0,0.,0.),
					   'vsun':(-11.1, 12.24, 7.25), 'vrot':218, 'gpos':(4.913,-9.772,-85.387),
					   'gvel':(-36.774, 163.875, -96.074), 'skyposg':(287.534, -83.156), 'skypos':(15.0392, -33.7092),
					   'rh_obs': 0.283, 'vdisp_obs':8.69, 'vdisp_obs_tot':9.43, 'outdir': None, 'Nresample':100000,
					   'file_vdisp':'Scl_binned_profile_s3_rc.txt', 'dist':86, 'Vlos':111.05, 'pmot':(0.011, 0.143),
					   'mstar':4.6e6, 'comtype': 2, 'orbit_track': None, 'Sigma_scale':1 }

		default_shortrange={'proffile':None, 'folder':None, 'radmax':1.2, 'psun':(8.0,0.,0.),
					   'vsun':(-11.1, 12.24, 7.25), 'vrot':218, 'gpos':(4.913,-9.772,-85.387),
					   'gvel':(-36.774, 163.875, -96.074), 'skyposg':(287.534, -83.156), 'skypos':(15.0392, -33.7092),
					   'rh_obs': 0.283, 'vdisp_obs':8.69, 'vdisp_obs_tot':9.40, 'outdir': None, 'Nresample':100000,
					   'file_vdisp':'Scl_binned_profile_s3_rc.txt', 'dist':86, 'Vlos':111.05, 'pmot':(0.011, 0.143),
					   'mstar':4.6e6, 'comtype': 2, 'orbit_track': None, 'Sigma_scale':1 }

		default_exttrange={'proffile':None, 'folder':None, 'radmax':1.9, 'psun':(8.195,0.,0.),
					   'vsun':(-9.833, 11.945, 8.107), 'vrot':219.386, 'gpos':(4.881,-10.490,-91.658),
					   'gvel':(-40.843, 48.310, -81.362), 'skyposg':(287.534, -83.156), 'skypos':(15.0392, -33.7092),
					   'rh_obs': 0.303, 'vdisp_obs':8.69, 'vdisp_obs_tot':9.43, 'outdir': None, 'Nresample':100000,
					   'file_vdisp':'Scl_binned_profile_s3_rc.txt', 'dist':92.316, 'Vlos':110.714, 'pmot':(-0.05843, 0.39278),
					   'mstar':4.6e6, 'comtype': 2, 'orbit_track': None, 'Sigma_scale':1 }

		default_exttranges={'proffile':None, 'folder':None, 'radmax':1.76, 'psun':(8.008,0.,0.),
					   'vsun':(-9.088, 10.286, 6.859), 'vrot':213.936, 'gpos':(4.983, -9.5487, -83.441),
					   'gvel':(-27.76599, 134.551590, -93.913965), 'skyposg':(287.534, -83.156), 'skypos':(15.0392, -33.7092),
					   'rh_obs': 0.276, 'vdisp_obs':8.69, 'vdisp_obs_tot':9.43, 'outdir': None, 'Nresample':100000,
					   'file_vdisp':'Scl_binned_profile_s3_rc.txt', 'dist':84.040, 'Vlos':110.914, 'pmot':(-0.0231058, 0.196979),
					   'mstar':4.6e6, 'comtype': 2, 'orbit_track': None, 'Sigma_scale':1}


		if default.lower()[0]=='l':
			self.default=default_longrange
		elif default.lower()[0]=='s':
			self.default=default_shortrange
		elif default.lower()[0]=='e':
			self.default=default_exttrange
		elif default.lower()[0]=='k':
			self.default=default_exttranges
		elif default.lower()[0]=='h':
			self.default=default_helmi

		self.description={'proffile':'File with the velocity dispersion', 'folder':'?', 'radmax':'maximum radius to consider in deg', 'psun':'Position of the Sun (X,Y,Z)',
					   'vsun':'Local velocty of the Sun (Vx, Vy, Vz)', 'vrot':'Velocity of LSR', 'gpos':'Galactic position of the object (Xg, Yg, Zg)',
					   'gvel':'Galactic velocity of the object (Xg, Yg, Zg)', 'skyposg': 'Position in sky coordinates (l [deg], b[deg])', 'skypos':'Position in equatorial coordinates (ra [deg], dec[deg])',
					   'rh_obs': 'Observed half light radius', 'vdisp_obs':'Observed velocity dispersion inside half-light radius', 'vdisp_obs_tot':'Observed velocity dispersion inside Rmax', 'outdir': 'Name of the output folder', 'Nresample':'DM particle to plot',
					   'file_vdisp':'File containing the velcoty dispersion profile', 'dist':'distance from the Sun ', 'Vlos':'Vlos wrt the Sun', 'pmot':'Proper motion (mul, mub) in mas/yr',
					   'mstar':'Stellar mass inside Rmax in solar masses', 'comtype': 'Use this type particles to calculate COM, if None use all particles', 'orbit_track':'File with orbits', 'Sigma_scale': 'Factor to rescale the observe Sigma profile'}

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
			elif key=='proffile' or key=='folder' or key=='outdir' or key=='file_vdisp' or key=='orbit_track':
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
orbit_track=par.orbit_track
Sigma_scale=par.Sigma_scale
#vdisp_obs=None
#vdisp_obs_tot=None

outdir=par.outdir
Nresample=par.Nresample
file_vdisp=par.file_vdisp
dist=par.dist
Vlos=par.Vlos
pmot=par.pmot
mstar=par.mstar

comtype=par.comtype

radmaxobs=radmax
radmax=dist*np.tan(radmaxobs*np.pi/180.)

#START
if folder is None: folder = '.'
if outdir is None: outdir = './analysis'
outdirdata=outdir+'/data'
if not os.path.exists(outdir):
	os.makedirs(outdir)
if not os.path.exists(outdirdata):
	os.makedirs(outdirdata)




if file_vdisp is not None:
	#ADJUST WITH SPACES
	file_vdisp=file_vdisp.rstrip()
	file_vdisp=file_vdisp.lstrip()
	#check whether file exists
	if  not os.path.exists(file_vdisp):
		print('file_vdisp %s not found, setting velocity dispersion profile to None'%file_vdisp)
		sys.stdout.flush()
		file_vdisp = None



if proffile is not None:
	#ADJUST WITH SPACES
	proffile=proffile.rstrip()
	proffile=proffile.lstrip()
	#check whether file exists
	if  not os.path.exists(proffile):
		print('proffile %s not found, setting surface density  profile to None'%file_vdisp)
		sys.stdout.flush()
		proffile = None
		
if orbit_track is not None:
	#ADJUST WITH SPACES
	orbit_track=orbit_track.rstrip()
	orbit_track=orbit_track.lstrip()
	#check whether file exists
	if  not os.path.exists(proffile):
		print('orbit_track %s not found, setting surface density  profile to None'%file_vdisp)
		sys.stdout.flush()
		orbit_track = None


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
figobs, axobs = plt.subplots(2, 2, tight_layout=True, figsize=(11.0, 10))
figstat, axstat = plt.subplots(1, 1, tight_layout=True)
figev2D, axev2d = plt.subplots(2, 2, tight_layout=True, figsize=(10, 10))
figev3D, axev3d = plt.subplots(1, 2, tight_layout=True, figsize=(10, 5))
figvrot, axvrot = plt.subplots(1, 3, tight_layout=True, figsize=(15, 5))

colortd = ('blue', 'darkgreen', 'red')
check_td = 0

idx_plot_orbit = (int(Nfiles / 3.) - 1, int(2 * Nfiles / 3.) - 1, Nfiles - 1)
#print(idx_plot_orbit)

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
if radmax is not None: log += 'Rmax=%.3f kpc\n' % radmax
if radmaxobs is not None: log += 'Rmax=%.3f deg \n' % radmaxobs
if mstar is not None: log += 'M*(R<Rmax)=%.3e\n' % mstar
if comtype is not None:  log += 'ComType=%i\n'%comtype
else: log += 'ComType=None'
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
		axarr[iplot, 0].text(-145, 125, '$t$=%.2f Gyr' % time_tmp, fontsize=20)
		#if i == idx_plot_orbit[0]:
		#	axarr[0, 0].text(0, 150, '$xy$', fontsize=25)
		#	axarr[0, 1].text(0, 150, '$xz$', fontsize=25)
		#	axarr[0, 2].text(0, 150, '$tz$', fontsize=25)

		if gpos is not None:
			axarr[iplot, 0].plot([gpos[0], gpos[0]], [-1000, 1000], color='darkgreen', zorder=2000, lw=0.5)
			axarr[iplot, 0].plot([-1000, 1000], [gpos[1], gpos[1]], color='darkgreen', zorder=2000, lw=0.5)
			axarr[iplot, 1].plot([gpos[0], gpos[0]], [-1000, 1000], color='darkgreen', zorder=2000, lw=0.5)
			axarr[iplot, 1].plot([-1000, 1000], [gpos[2], gpos[2]], color='darkgreen', zorder=2000, lw=0.5)
			axarr[iplot, 2].plot([gpos[1], gpos[1]], [-1000, 1000], color='darkgreen', zorder=2000, lw=0.5)
			axarr[iplot, 2].plot([-1000, 1000], [gpos[2], gpos[2]], color='darkgreen', zorder=2000, lw=0.5)
		
		'''  
		if orbit_track is not None: 
			dorbit=np.loadtxt(orbit_track)
			idxtorb=np.searchsorted(dorbit[:,0], time_tmp)
			if idxtorb!=0:
				dorbit=dorbit[:idxtorb]
				#print(dorbit.shape)
				axarr[iplot, 0].plot(dorbit[:, 1], dorbit[:, 2],lw=0.8)
				axarr[iplot, 1].plot(dorbit[:, 1], dorbit[:, 3],lw=0.8)
				axarr[iplot, 2].plot(dorbit[:, 2], dorbit[:, 3],lw=0.8)
			del dorbit
		'''
				
		if (i == idx_plot_orbit[2]):
			
			axarr[0, 0].set_ylabel('y [kpc]', fontsize=20)
			axarr[0, 0].set_xlabel('x [kpc]', fontsize=20)
			axarr[1, 0].set_ylabel('y [kpc]', fontsize=20)
			axarr[1, 0].set_xlabel('x [kpc]', fontsize=20)
			axarr[2, 0].set_ylabel('y [kpc]', fontsize=20)
			axarr[2, 0].set_xlabel('x [kpc]', fontsize=20)
			
			axarr[0, 1].set_ylabel('z [kpc]', fontsize=20)
			axarr[0, 1].set_xlabel('x [kpc]', fontsize=20)
			axarr[1, 1].set_ylabel('z [kpc]', fontsize=20)
			axarr[1, 1].set_xlabel('x [kpc]', fontsize=20)
			axarr[2, 1].set_ylabel('z [kpc]', fontsize=20)
			axarr[2, 1].set_xlabel('x [kpc]', fontsize=20)
			
			axarr[0, 2].set_ylabel('z [kpc]', fontsize=20)
			axarr[0, 2].set_xlabel('y [kpc]', fontsize=20)
			axarr[1, 2].set_ylabel('z [kpc]', fontsize=20)
			axarr[1, 2].set_xlabel('y [kpc]', fontsize=20)
			axarr[2, 2].set_ylabel('z [kpc]', fontsize=20)
			axarr[2, 2].set_xlabel('y [kpc]', fontsize=20)
			
			#axarr[2, 0].set_xlabel('[kpc]', fontsize=20)
			#axarr[2, 1].set_xlabel('[kpc]', fontsize=20)
			#axarr[2, 2].set_xlabel('[kpc]', fontsize=20)
			#axarr[0, 0].set_ylabel('x [kpc]', fontsize=20)
			#axarr[1, 0].set_ylabel('[kpc]', fontsize=20)
			#axarr[2, 0].set_ylabel('[kpc]', fontsize=20)
			figorbit.set_size_inches(15, 15, forward=True)
			# plt.setp([a.get_xticklabels() for a in figorbit.axes[:-1]], visible=False)
			figorbit.savefig(outdir + '/orbit.png')
			del axarr
			del figorbit
			


			
		iplot += 1

		print('Done', flush=True)

	# Plot analysis
	print('Find COM', flush=True)
	a_tmp = Analysis(particles=p_tmp, safe=True, auto_centre=True, iter=True, single=False, type_com=comtype)
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
		have_halo=False        

        
        
		# dens
		try:
			prof_tmp_h = Profile(particles=a_tmp.p, xmin=0.01, xmax=100, ngrid=100, kind='log', type=1)
			arr = prof_tmp_h.dens()[0]
			r = arr[:, 0]
			d = arr[:, 1]
			axtd[0].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
			np.savetxt(outdirdata+'/3Dprofile_DM_T%.2f.txt'%time_tmp,arr,fmt='%.3e',header='0-r [kpc] 1-dens [Msun kpc^-3]') 
            
			arr_m_h, f_m_h = prof_tmp_h.mass()
			np.savetxt(outdirdata+'/MassprofileDM_T%.2f.txt'%time_tmp,arr_m_h,fmt='%.3e',header='0-r [kpc] 1-Mhalo [Msun]')
			have_halo=True
		except:
			pass

		prof_tmp_s = Profile(particles=a_tmp.p, xmin=0.0001, xmax=5, ngrid=100, kind='lin', type=2)
		arr = prof_tmp_s.dens()[0]
		r = arr[:, 0]
		d = arr[:, 1]
		axtd[1].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
		np.savetxt(outdirdata+'/3Dprofile_stars_T%.2f.txt'%time_tmp,arr,fmt='%.3e',header='0-r [kpc] 1-dens [Msun kpc^-3]')     

		#mass
		arr_m_s = prof_tmp_s.mass()[0]
		arr_m   = np.zeros(shape=(len(arr_m_s),4))
		arr_m[:,0] = arr_m_s[:,0]   
		arr_m[:,1] = arr_m_s[:,1]
		if have_halo:  arr_m[:,2] = f_m_h(arr_m_s[:,0])   
		arr_m[:,3] = arr_m_s[:,1] +  arr_m[:,2]      
		np.savetxt(outdirdata+'/Massprofile_T%.2f.txt'%time_tmp,arr_m,fmt='%.3f %.3e %.3e %.3e',header='0-r [kpc] 1-Mstar [Msun] 2-Mhalo [Msun] 3-Mtot [Msun]') 
        
		prof_tmp_s2 = Profile(particles=a_tmp.p, xmin=0.0001, xmax=5, ngrid=20, kind='lin', type=2)
		arrx = prof_tmp_s2.vdisp3d(ax='x')[0]
		arry = prof_tmp_s2.vdisp3d(ax='y')[0]
		arrz = prof_tmp_s2.vdisp3d(ax='z')[0]
		r = arrx[:, 0]
		vdx = arrx[:, 1]
		vdy = arry[:, 1]
		vdz = arrz[:, 1]
		vd = np.sqrt(vdx*vdx + vdy*vdy + vdz*vdz)
		arr=np.zeros((len(arrx),5))
		arr[:,0]=r
		arr[:,1]=vdx
		arr[:,2]=vdy
		arr[:,3]=vdz
		arr[:,4]=vd
		axtd[2].plot(r, vd, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
		np.savetxt(outdirdata+'/3DVdisp_stars_T%.2f.txt'%time_tmp,arr,fmt='%.3f %.3f %.3f %.3f %.3f',header='0-r [kpc] 1-Vdispx 2-Vdispy 3-Vdispz 4-Vdisptot [km/s]') 
        
		arr = prof_tmp_s.supdens(pax='z')[0]
		arr_output= np.zeros(shape=(len(arr),4))
		r = arr[:, 0]
		d = arr[:, 1]
		arr_output[:,1]= d
		dminS0 = np.min(d)
		dmaxS0 = np.max(d)
		axstd[0, 0].plot([rh_tmp_z, rh_tmp_z], [np.min(d), np.max(d)], color=colortd[check_td])
		axstd[0, 0].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])

		arr = prof_tmp_s.supdens(pax='y')[0]
		r = arr[:, 0]
		d = arr[:, 1]
		arr_output[:,2]= d
		dminS1 = np.min(d)
		dmaxS1 = np.max(d)
		axstd[0, 1].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
		axstd[0, 1].plot([rh_tmp_y, rh_tmp_y], [np.min(d), np.max(d)], color=colortd[check_td])

		arr = prof_tmp_s.supdens(pax='x')[0]
		r = arr[:, 0]
		d = arr[:, 1]
		arr_output[:,3]= d
		dminS2 = np.min(d)
		dmaxS2 = np.max(d)
		axstd[0, 2].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
		axstd[0, 2].plot([rh_tmp_x, rh_tmp_x], [np.min(d), np.max(d)], color=colortd[check_td])  # vdisp
		arr_output[:,0]= r
		np.savetxt(outdirdata+'/Sigma_stars_T%.2f.txt'%time_tmp,arr_output,fmt='%.3e',header='0-r [kpc] 1-Sigma_x [Msun / kpc^2] 2-Sigma_y [Msun / kpc^2] 3-Sigma_z [Msun / kpc^2]')  

		arr = prof_tmp_s.vdisp2d(pax='z')[0]
		arr_output= np.zeros(shape=(len(arr),4))
		r = arr[:, 0]
		d = arr[:, 1]
		arr_output[:,1]= d
		axstd[1, 0].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
		arr = prof_tmp_s.vdisp2d(pax='y')[0]
		r = arr[:, 0]
		d = arr[:, 1]
		arr_output[:,2]= d
		axstd[1, 1].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
		arr = prof_tmp_s.vdisp2d(pax='x')[0]
		r = arr[:, 0]
		d = arr[:, 1]
		arr_output[:,3]= d
		axstd[1, 2].plot(r, d, lw=2, label='T=%.2f Gyr' % time_tmp, color=colortd[check_td])
		arr_output[:,0]= r
		np.savetxt(outdirdata+'/2DVdisp_stars_T%.2f.txt'%time_tmp,arr_output,fmt='%.3f',header='0-r [kpc] 1-Vdisp_x [km/s] 2-Vdisp_y [km/s] 3-Vdisp_z [km/s]')   
        
        
        
		#Vsys
		arr = prof_tmp_s.vsys(pax='z')[0]
		r = arr[:, 0]
		d = arr[:, 1]
		vsysmax=[np.max(d),]
		vsysmin=[np.min(d),]
		axvrot[0].plot(r, d, label='$t$=%.2f Gyr'%time_tmp, color=colortd[check_td])
		arr = prof_tmp_s.vsys(pax='y')[0]
		r = arr[:, 0]
		d = arr[:, 1]
		vsysmax.append(np.max(d))
		vsysmin.append(np.min(d))
		axvrot[1].plot(r, d, label='$t$=%.2f Gyr'%time_tmp, color=colortd[check_td])
		arr = prof_tmp_s.vsys(pax='x')[0]
		r = arr[:, 0]
		d = arr[:, 1]
		vsysmax.append(np.max(d))
		vsysmin.append(np.min(d))
		axvrot[2].plot(r, d, label='$t$=%.2f Gyr'%time_tmp, color=colortd[check_td])
		axvrot[0].axhline(0, color='black', ls='--')
		axvrot[1].axhline(0, color='black', ls='--')
		axvrot[2].axhline(0, color='black', ls='--')

		check_td += 1

		if check_td == 3:

			#Vsys
			axvrot[0].set_xlabel('$R=\\sqrt{x^2+y^2} \$ [kpc]', fontsize=20)
			axvrot[1].set_xlabel('$R=\\sqrt{x^2+z^2} \$ [kpc]', fontsize=20)
			axvrot[2].set_xlabel('$R=\\sqrt{y^2+z^2} \$ [kpc]', fontsize=20)
			axvrot[0].set_ylabel('$V_{sys}  \  \mathrm{[km \ s^{-1}]}$', fontsize=20)
			axvrot[0].set_title('$z$-Projection', fontsize=20)
			axvrot[1].set_title('$y$-Projection', fontsize=20)
			axvrot[2].set_title('$x$-Projection', fontsize=20)
			figvrot.savefig(outdir+'/Vsys.png')
			#axvrot_min=np.min(vsysmin)
			#axvrot_max=np.max(vsysmax)
			#axvrot[0].set_xlim(axvrot_min*1.5,ax)


			# 3D
			axtd[0].set_xlabel('$r$ [kpc]', fontsize=20)
			axtd[1].set_xlabel('$r$ [kpc]', fontsize=20)
			axtd[2].set_xlabel('$r$ [kpc]', fontsize=20)
			axtd[0].set_ylabel('$\\rho_h \\  \mathrm{[M_\\odot \ kpc^{-3}]}$', fontsize=20)
			axtd[1].set_ylabel('$\\rho_* \\  \mathrm{[M_\\odot \ kpc^{-3}]}$', fontsize=20)
			axtd[2].set_ylabel('$\\sigma_{3D} \\  \mathrm{[km \ s^{-1}]}$', fontsize=20)
			axtd[0].set_xlim(0.1, 100)
			axtd[0].set_ylim(1e-2, 1e10)
			axtd[0].set_xscale('log')
			axtd[0].set_yscale('log')
			axtd[1].set_xlim(0.01, 10)
			axtd[1].set_xscale('log')
			axtd[1].set_yscale('log')
			axtd[2].set_xlim(0, 2)
			axtd[0].legend(loc='upper right',ncol=3,fontsize=16)
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
				if proffile is not None:
					datapr = np.loadtxt(proffile)
					x = dist * np.tan((datapr[:, 0]/60.) * (np.pi) / 180)
					axstd[0, 0].errorbar(x, datapr[:, 1]*Sigma_scale, datapr[:, 2]*Sigma_scale, fmt='o', c='black', ms=4, capsize=4,  zorder=1000,alpha=0.5,label='DATA (rescaled)')
					axstd[0, 1].errorbar(x, datapr[:, 1]*Sigma_scale, datapr[:, 2]*Sigma_scale, fmt='o', c='black', ms=4, capsize=4,  zorder=1000,alpha=0.5,label='DATA (rescaled)')
					axstd[0, 2].errorbar(x, datapr[:, 1]*Sigma_scale, datapr[:, 2]*Sigma_scale, fmt='o', c='black', ms=4, capsize=4,  zorder=1000,alpha=0.5,label='DATA (rescaled)')

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

			axstd[1, 0].set_xlabel('$R=\\sqrt{X^2+Y^2} \$  [kpc]', fontsize=20)
			axstd[1, 1].set_xlabel('$R=\\sqrt{X^2+Z^2} \$  [kpc]', fontsize=20)
			axstd[1, 2].set_xlabel('$R=\\sqrt{Y^2+Z^2} \$  [kpc]', fontsize=20)
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
			axstd[1, 0].set_ylabel('$\\sigma \\  \mathrm{[km \ s^{-1}]}$', fontsize=20)
			axstd[0, 0].set_ylabel('$\\Sigma \\  \mathrm{[M_\\odot \ kpc^{-2}]}$', fontsize=20)
			axstd[0, 0].legend(loc='upper right')
			axstd[1, 0].legend(loc='upper left')
			# axstd[0, 0].text(0.5, 1.0, 'Projection ax: X', fontsize=25, transform=axstd[0, 0].transAxes)
			# axstd[0, 1].text(0.5, 1.0, 'Projection ax: Y', fontsize=25, transform=axstd[0, 1].transAxes)
			# axstd[0, 2].text(0.5, 1.0, 'Projection ax: Z', fontsize=25, transform=axstd[0, 2].transAxes)
			figanstd.set_size_inches(15, 10, forward=True)
			figanstd.savefig(outdir + '/2Danalysis.png')

	# observe
	o_tmp = Observe(particles=p_tmp, type=2)
	s, c = o_tmp.observe(psun=psun, vsun=vsun, vrot=vrot, com=(com_tmp, vcom_tmp))
	st, ct = o_tmp.observe(psun=psun, vsun=vsun, vrot=vrot, com=(gpos, gvel))
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

		#axobs[0,0].axis('equal')
		#axobs[0,1].axis('equal')
		#axobs[1,0].axis('equal')
		#axobs[1,1].axis('equal')
		xict,etact=radec_to_xieta(ct.mura, ct.mudec, c.mura, c.mudec)


		prof_obs_large = Profile(particles=s, xmin=0.01, xmax=10, ngrid=512, kind='lin')
		prof_obs_small = Profile(particles=s, xmin=0.01, xmax=10, ngrid=50, kind='lin')

	   # ifilter=(np.abs(s.xi/3600.)<3)&(np.abs(s.eta/3600.)<3)

		H, xedges, yedges = np.histogram2d(s.xi[:] / 3600., s.eta[:] / 3600., bins=(30,30),range=((-2,2),(-2.,2.)))

		pixel_area=np.abs(xedges[1]-xedges[0])*np.abs(yedges[1]-yedges[0])*60*60
		H=H/pixel_area
		levels=(0.01, 0.05, 0.1,1.0,5,10,50,100)
		#axobs[0,0].imshow(H.T,origin='lower',extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]),aspect='auto',cmap='YlOrBr')
		axobs[0,0].contour(H.T,origin='lower',extent=(xedges[0],xedges[-1],yedges[0],yedges[-1]),zorder=10000,levels=levels,colors='black')
		theta=np.linspace(0,2*np.pi,1000)
		RRR=radmaxobs
		xxx=RRR*np.cos(theta)
		yyy=RRR*np.sin(theta)
		axobs[0,0].plot(xxx,yyy,'--', color='blue',lw=2,zorder=30000,label='$R=%.1f^\circ$'%(RRR))
		axobs[0,0].scatter(s.xi[:] / 3600., s.eta[:] / 3600., s=0.005, c='red')
		axobs[0,0].scatter(xict/ 3600., etact  / 3600., s=100, c='blue', marker='X', zorder=1000000, label='Sculptor centre')
		
		axobs[0,0].quiver(xict / 3600., etact  / 3600., ct.mura, ct.mudec, angles='uv', scale=0.5,zorder=300010, linewidth=2.5,  label='PM Observed',linestyle='dashed',color='blue')
		axobs[0,0].quiver(0, 0, c.mura, c.mudec, angles='uv', scale=0.5, width=0.01, headwidth=6, headlength=5,zorder=100000,label='PM Simulation')
		axobs[0,0].scatter(1e6 / 3600.,0, c='red',label='Star particles')
		axobs[0,0].set_xlabel('$\\xi$  [deg]', fontsize=20)
		axobs[0,0].set_ylabel('$\\eta$  [deg]', fontsize=20)
		axobs[0,0].set_xlim(-2.0, 2.0)
		axobs[0,0].set_ylim(-2.0, 2.0)
		axobs[0,0].plot([1e6,1e6],[1e6,1e9],color='black',label='Iso-density')
		#axobs[0,0].legend(ncol=2)
		
		#Vdisp
		prof_obs=prof_obs_small
		arr = prof_obs.vdisp2d(pax='obs')[0]
		r = arr[:, 0]
		vd = arr[:, 1]
		axobs[1,0].plot(r, vd, lw=3, color='red', zorder=2000,label='Simulation')
		np.savetxt(outdirdata+'/Vdisplos_T%.2f.txt'%time_tmp,arr,fmt='%.3e %.3f')        

		if file_vdisp is not None:
			try:
				data = np.loadtxt(file_vdisp)
				x = dist * np.tan(data[:, 0] * (np.pi) / 180)
				ex = dist * np.tan(data[:, 1] * (np.pi) / 180)
				axobs[1,0].errorbar(x, data[:, 4], data[:, 5], ex, fmt='o', c='black', zorder=1000,label='Data (this work)')

				obins = np.zeros(len(data) + 1)
				obins[:-1] = data[:, 6]
				obins[-1] = data[-1, 7]
				Nperbin = data[:, 8]
				obins = dist * np.tan(obins * np.pi / 180.)

				color_disp = ('blue', 'orange', 'magenta')

				for j in range(3):
					b = a.binned_dispersion(bins=obins, pax='obs', Nperbin=Nperbin, bins_kind='lin', velocity_err=None,
											err_distibution='uniform', nboot=10000)

					axobs[1,0].errorbar(b[0], b[4], b[5], b[1], fmt='o', c=color_disp[j], mfc='white',label='Realisation %i'%j)
					outbinrel=np.zeros((len(b[0]),4))
					outbinrel[:,0]=b[0]
					outbinrel[:,1]=b[1]
					outbinrel[:,2]=b[4]
					outbinrel[:,3]=b[5]
					np.savetxt(outdirdata+'/vdisprel_%i.txt'%j,outbinrel,fmt='%.3e %.3e %.3f %.3f',header='0-R 1-eR 2-V 3-eV Vdisp_tot=%.3f'%b[-1])
			except FileNotFoundError:
				print('File %s not found.. skipping' % file_vdisp)

		axobs[1,0].set_xlabel('$R$   [kpc]', fontsize=20)
		axobs[1,0].set_ylabel('$\\sigma_\mathrm{los}  \  \mathrm{[km \ s^{-1}]}$', fontsize=20)
		axobs[1,0].set_xlim(0, 2)
		axobs[1,0].set_ylim(0, 15)

		# densup
		#Vdisp
		prof_obs=prof_obs_large

		arr,supdensfunc = prof_obs.supdens(pax='obs',ret=True,func=True,s=0)
		r = arr[:, 0]
		d = arr[:, 1]
		np.savetxt(outdirdata+'/Sigmalos.txt',arr,fmt='%.3e',header='0-R [kpc] 1-Sigma [Msun/kpc^2] Rh=%.3f'%rh_sim)
		axobs[0,1].plot(r, d, lw=3, color='red')
		axobs[0,1].plot([rh_sim, rh_sim], [np.min(d), np.max(d)], color='magenta', lw=2, label='$R^{sim}_h$')
		if rh_obs is not None:
			axobs[0,1].plot([rh_obs, rh_obs], [np.min(d), np.max(d)], '--',  lw=1.5, color='black', label='$R^{obs}_h$  (McConnachie12)')
		axobs[0,1].set_xlabel('$R$  [kpc]', fontsize=20)
		axobs[0,1].set_ylabel('$\\Sigma_* \  \mathrm{[M_\\odot \ kpc^{-2}]}$', fontsize=20)
		axobs[0,1].set_xlim(0.01, 2)
		#axobs[0,1].set_ylim(1e3, 4e7)
		# axobs[1].set_xlim(0.001,10)
		axobs[0,1].set_xscale('log')
		axobs[0,1].set_yscale('log')

		if proffile is not None:
			try:
				datapr = np.loadtxt(proffile)
				x = dist * np.tan((datapr[:, 0]/60.) * (np.pi) / 180)

				#idxmed=int(len(x)/2)
				#xmed=x[idxmed]
				#nnorm=supdensfunc(xmed)

				#dnorm=datapr[idxmed,1]
				#fnorm=nnorm/dnorm
				#datapr[:,1]=fnorm*datapr[:,1]
				#datapr[:,2]=fnorm*datapr[:,2]

				axobs[0,1].errorbar(x, datapr[:, 1]*Sigma_scale, datapr[:, 2]*Sigma_scale, fmt='o', c='black', ms=4, capsize=4,  zorder=1000,alpha=0.5,label='DATA (rescaled)')
			except FileNotFoundError:
				print('File %s not found.. skipping' % file_vdisp)
		proffile=None

		#Vsys
		prof_obs=prof_obs_small
		arr = prof_obs.vsys(pax='obs')[0]
		r = arr[:, 0]
		d = arr[:, 1]
		axobs[1,1].plot(r, d, lw=3, color='red', label='Simulation')
		if Vlos is not None:
			axobs[1,1].axhline(Vlos, color='black', ls='--', label='Observed (this work)')
			Vlose=Vlos
		else:
			medVlos=np.median(d)
			axobs[1,1].axhline(medVlos, color='black', ls='--', label='Median')
			Vlose=medVlos
		axobs[1,1].set_xlabel('$R$  [kpc]', fontsize=20)
		axobs[1,1].set_ylabel('$V_\mathrm{los} \   \mathrm{[km \ s^{-1}]}$', fontsize=20)
		axobs[1,1].set_xlim(0,2)
		axobs[1,1].set_ylim(Vlose-20,Vlose+20)

		#axobs[0,0].legend(loc='upper center',fontsize=14,ncol=2)
		axobs[0,1].legend(loc='best',fontsize=14)
		axobs[1,0].legend(loc='lower right',fontsize=14)
		axobs[1,1].legend(loc='best',fontsize=14)
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
	label = '$M^{sim}_*(<R_{*})$ wrt $t$=0'
else:
	mp_ini = mstar
	label = '$M^{sim}_*(<R_{*})$ wrt obs'
rel_diff = 100 * (mp - mp_ini) / mp_ini
axstat.plot(T, rel_diff, '-^', c='darkgreen', label=label, zorder=3000)

vdp_tot = fin_array[1:, 10]
if vdisp_obs_tot is None:
	vdp_tot_ini = vdp_tot[0]
	label = '$\\sigma^{sim}_{los}(<R_{*})$ wrt $t$=0'
else:
	vdp_tot_ini = vdisp_obs_tot
	label = '$\\sigma^{sim}_{los}(<R_{*})$ wrt obs'
rel_diff = 100 * (vdp_tot - vdp_tot_ini) / vdp_ini
axstat.plot(T, rel_diff, '-s', c='cyan', label=label, zorder=3000)
axstat.legend(loc='best', fontsize=12)

axstat.axhline(0, color='black', ls='--')
axstat.set_xlabel('$t$ [Gyr]', fontsize=20)
axstat.set_ylabel('Rel. diff [%]', fontsize=20)

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

axev2d[1, 0].plot(T2, vv_x, '-o', label='$x$-projection', lw=3, markersize=10)
axev2d[1, 0].plot(T2, vv_y, '-^', label='$y$-projection', lw=3, markersize=10)
axev2d[1, 0].plot(T2, vv_z, '-s', label='$z$-projection', lw=3, markersize=10)
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
axev2d[1, 1].set_title('$M(R<R_{max})$', fontsize=20)

axev2d[0, 0].legend(loc='best', fontsize=12)
axev2d[1, 0].set_xlabel('$t$ [Gyr]', fontsize=20)
axev2d[1, 1].set_xlabel('$t$ [Gyr]', fontsize=20)
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
axev3d[0].plot(T2, pf, '-o', label='$y$-to-$x$ axial ratio', lw=3, markersize=10)
axev3d[0].plot(T2, qf, '-s', label='$z$-to-$x$ axial ratio', lw=3, markersize=10)

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
axev3d[1].set_xlabel('$t$  [Gyr]', fontsize=20)
axev3d[0].set_xlabel('$t$  [Gyr]', fontsize=20)
axev3d[0].set_ylabel('Axial ratio', fontsize=20)
axev3d[1].set_ylabel('Rel diff [%]', fontsize=20)

figev3D.savefig(outdir + '/evolution_3D.png')

np.savetxt(outdir + '/data_obs.txt', fin_array, fmt='%.3e', header=fin_array_header)
np.savetxt(outdir + '/data_intrinsic.txt', fin_array_intrinsic, fmt='%.3e', header=fin_array_intrinsic_header)
file = open(outdir + '/info.log', mode='w')



#Correction for the next Run
#Stellar Mass
Mstarnew_slow= (mstar/mass_sim)**0.5
Mstarnew_fast= (mstar/mass_sim)

Mvirnew_slow= (vdisp_obs_tot/vdisp_sim_tot)**0.5
Mvirnew_fast= (vdisp_obs_tot/vdisp_sim_tot)

rcnew_slow= (vdisp_sim/vdisp_obs)**0.5
rcnew_fast= (vdisp_sim/vdisp_obs)

bnew_slow= (rh_obs/rh_sim)**0.5
bnew_fast= (rh_obs/rh_sim)


if orbit_track is not None and gpos is not None and gvel is not None:
	dorbit=np.loadtxt(orbit_track)
	difftot=np.abs(dorbit[:,1]-com_tmp[0])
	difftot+=np.abs(dorbit[:,2]-com_tmp[1])
	difftot+=np.abs(dorbit[:,3]-com_tmp[2])
	difftot+=np.abs(dorbit[:,5]-vcom_tmp[0])
	difftot+=np.abs(dorbit[:,6]-vcom_tmp[1])
	difftot+=np.abs(dorbit[:,7]-vcom_tmp[2])
	idxmineff=np.argmin(difftot)
	Teff=dorbit[idxmineff]
	
	com_tmpp=gpos
	vcom_tmpp=gvel
	difftot=np.abs(dorbit[:,1]-com_tmpp[0])
	difftot+=np.abs(dorbit[:,2]-com_tmpp[1])
	difftot+=np.abs(dorbit[:,3]-com_tmpp[2])
	difftot+=np.abs(dorbit[:,5]-vcom_tmpp[0])
	difftot+=np.abs(dorbit[:,6]-vcom_tmpp[1])
	difftot+=np.abs(dorbit[:,7]-vcom_tmpp[2])   
	idxminobj=np.argmin(difftot)	
	Tobj=dorbit[idxminobj]
	
	Toffset=Tobj[0]-Teff[0]
	
	del dorbit

	
if gpos is not None and gvel is not None:
	log+='\n\n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
	log+='FINAL ANALYSIS\n'
	log+='Observed Pos: X=%.2f Y=%.2f Z=%.2f\n'%gpos
	log+='Last Simulation Pos: X=%.2f Y=%.2f Z=%.2f\n'%tuple(com_tmp)
	log+='Observed Vel: VX=%.2f VY=%.2f VZ=%.2f\n'%gvel
	log+='Last Simulation Vel: VX=%.2f VY=%.2f VZ=%.2f\n'%tuple(vcom_tmp)
	log+='%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'

if mstar is not None and vdisp_obs is not None and vdisp_obs_tot is not None and rh_obs is not None:
	log+='\n\n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'
	log+='CORRECTION FACTORS\n'
	log+='f(Rh)= %.6f (slow) %.6f (fast)\n'%(bnew_slow, bnew_fast)
	log+='f(Mstar)= %.6f (slow) %.6f (fast)\n'%(Mstarnew_slow, Mstarnew_fast)
	log+='f(rc)= %.6f (slow) %.6f (fast)\n'%(rcnew_slow, rcnew_fast)
	log+='f(Mvir)= %.6f (slow) %.6f (fast)\n'%(Mvirnew_slow, Mvirnew_fast)
if orbit_track is not None:
	log+='TORB OFFSET\n'
	log+='Target time: %.4f (idx: %i)\n'%(Tobj[0], idxminobj)
	log+='Effective time: %.4f (idx: %i)\n'%(Teff[0], idxmineff)
	log+='Toffset: %.4f\n'%Toffset
	log+='\n'
	log+='Properties at Ttarget: Pos(%.3f, %.3f, %.3f) Vel(%.3f, %.3f, %.3f)\n'%(Tobj[1], Tobj[2], Tobj[3], Tobj[5], Tobj[6], Tobj[7])
	log+='Properties at Teffective: Pos(%.3f, %.3f, %.3f) Vel(%.3f, %.3f, %.3f)\n'%(Teff[1], Teff[2], Teff[3], Teff[5], Teff[6], Teff[7])
log+='%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n'

	


	
file.write(log)
file.close()
