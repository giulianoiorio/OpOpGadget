from OpOp.particle import Particle,Particles, Sky_Particle, Sky_Particles
from OpOp.analysis import Observe
import matplotlib.pyplot as plt
import numpy as np

"""
p=Sky_Particle(l=60,b=30,distance=1,mul=0,mub=0,vlos=10)
print(p)

p=Sky_Particle(l=60,b=30,distance=1,mul=0,mub=0,vlos=10,centre_loc=(60,0))
print(p)

p=Sky_Particle(l=60,b=30,distance=1,mul=0,mub=0,vlos=10,centre_loc=(60,30))
print(p)
"""

p1=Sky_Particle(l=60,b=0,distance=1,mul=0,mub=0,vlos=10,id=0)
p2=Sky_Particle(l=30,b=0,distance=1,mul=0,mub=0,vlos=10,id=1)
p3=Sky_Particle(l=60,b=30,distance=1,mul=0,mub=0,vlos=10,id=2)

p=Sky_Particles(p=[p1,p2,p3])

con='l>30'
p2=p.extract('l>30','b=0',mode='or')
print(p2['id'])
p2=p.extract('l>30','b=0',mode='and')
print(p2['id'])
cod=('l=30','b=30')
p2=p.extract(*cod,mode='or')
print(p2['id'])
