import numpy as np
from scipy.optimize import brentq as _brentq

_iter = 1000
_xtol = 2e-12
_rtol = 4 * np.finfo(float).eps

def brentq(f, x0, dx, args=(), xtol=_xtol, rtol=_rtol, max_loops=_iter, full_output=False, disp=True):
    if f(x0+dx) * f(x0-dx) > 0:
        for iter in range(max_loops):
            dx *= 2
            if f(x0+dx) * f(x0-dx) < 0:
                break
    return _brentq(f, x0-dx, x0+dx, args, xtol, rtol, max_loops, full_output, disp)

