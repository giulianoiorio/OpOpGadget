import numpy as np

class Model:

    def __init__(self):
        self.use_c=use_c

    def dens(self,R):
        return self._evaluatedens(R)

    def sdens(self,R):
        return self._evaluatesdens(R)

    def mass(self,R):
        if self.use_c==True: mass=self._evaluatemassc
        elif self._use_nparray==True: mass=self._evaluatemass
        else: mass=np.vectorize(self._evaluatemass)
        return mass(R)

    def pot(self,R):
        if self.use_c==True: pot=self._evaluatepotc
        elif self._use_nparray==True: pot=self._evaluatepot
        else:pot=np.vectorize(self._evaluatepot)
        return pot(R)







