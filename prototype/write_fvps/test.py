import numpy as np
from OpOp.particle import Particle, Particles
from OpOp.io import write_snap, load_snap, write_icparam_fvfps

#Prop Iniziali Sculptor dopo avet integrato all indietro per 8 Gyr
Xg=4.78628
Yg=127.32181
Zg=43.99764
Vx=-11.92509
Vy=-45.72985
Vz=128.48708

ps=Particle(mass=1,id=0,type=1,pos=(Xg,Yg,Zg),vel=(Vx,Vy,Vz))
pf=Particle(mass=1,id=1,type=1,pos=(0,0,0),vel=(0,0,0))
print(ps)
print(pf)

plist=(ps,pf)
p=Particles(p=plist)
p.header['Massarr']=[[0, 1, 0, 0, 0, 0]]
print(p.header)
print(p[0])
print(p[1])

write_snap(p,'prova.bin',kind='fvfps')

pl=load_snap('prova.bin')
print(pl.header)
print(pl[0])
print(pl[1])

write_icparam_fvfps()