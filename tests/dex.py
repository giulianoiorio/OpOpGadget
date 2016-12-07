from functools import wraps
from scipy.integrate import quad

def cvar(rs):
    def cvartmp(func):
        #@wraps(func)
        def func_wrapper(r,*args):
            return func(r*rs,*args)
        return func_wrapper
    return cvartmp

@cvar(10)
def f(r):
    return r

print(quad(f,0,1))
print(f.__name__)
print(f.__doc__)