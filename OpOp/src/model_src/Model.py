import numpy as np


class Model:

    def __init__(self):
        self.use_c=use_c
        self._use_nparray=False

    def dens(self,R):
        return self._evaluatedens(R)

    def sdens(self,R,*args):
        return self._evaluatesdens(R,*args)

    def mass(self,R,*args):
        if self.use_c==True: mass=self._evaluatemassc
        elif self._use_nparray==True: mass=self._evaluatemass
        else: mass=np.vectorize(self._evaluatemass)
        return mass(R,*args)

    def pot(self,R,*args):
        if self.use_c==True: pot=self._evaluatepotc
        elif self._use_nparray==True: pot=self._evaluatepot
        else:pot=np.vectorize(self._evaluatepot)
        return pot(R,*args)

    def __add__(self, other):
        res=Multimodel(self,other)

        return res

class Multimodel(Model):

    def __init__(self,*args):

        self.mcomp=args
        self.use_c=False
        self._use_nparray=False

    def _evaluatedens(self,R):

        if isinstance(R,(float,int)):
            ret=0
            for comp in self.mcomp:
                ret+=comp.dens(R)
        else:
            ret=np.zeros(len(R))
            for comp in self.mcomp:
                ret+=comp.dens(R)

        return ret

    def _evaluatesdens(self,R):

        if isinstance(R,(float,int)):
            ret=0
            for comp in self.mcomp:
                ret+=comp.sdens(R)
        else:
            ret=np.zeros(len(R))
            for comp in self.mcomp:
                ret+=comp.sdens(R)

        return ret

    def _evaluatemass(self,R):

        if isinstance(R,(float,int)):
            ret=0
            for comp in self.mcomp:
                ret+=comp.mass(R)
        else:
            ret=np.zeros(len(R))
            for comp in self.mcomp:
                ret+=comp.mass(R)

        return ret

    def _evaluatepot(self,R):

        if isinstance(R,(float,int)):
            ret=0
            for comp in self.mcomp:
                ret+=comp.pot(R)
        else:
            ret=np.zeros(len(R))
            for comp in self.mcomp:
                ret+=comp.pot(R)

        return ret






