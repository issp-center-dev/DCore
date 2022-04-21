import numpy

def isfloating(x):
    """ Return whether x is a floating number (real) """
    return isinstance(x, (numpy.floating, float))