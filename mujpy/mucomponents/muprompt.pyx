%%cython
# Muprompt now written in Cython
# The %%cython IPython does the following:
# * Call Cython to generate C code for a Python C extension.
# * Compile it into a Python C extension (a shared library)
# * Load it into the current namespace
# If you don't understand these things, don't worry, it basically means:
# * Get full-metal speed easily
cimport cython
from libc.math cimport exp, M_PI, sqrt, erf
#from scipy.special cimport erf
#from numpy cimport sqrt, pi, exp, erf 
@cython.binding(True) # IMPORTANT: this tells Cython to dump the function signature    

def muprompt(int x, double a, double x0, double dx, double ak1, double ak2): 
    # fit function for a PSI prompt, 
    # data  first row is bin number, second row optional is histo content, third is error on data
    # a gaussian peak coincident with the edge betwee two plateaus (a constant + an erf)
    # par contains peak height, peak center, peak std width, first plateau, second plateau
    cdef double f
    f = a/(sqrt(2.*M_PI)*dx) * exp(-.5*((x-x0)/dx)**2)+ak2/2.+ak1+ak2/2.*erf((x-x0)/sqrt(2.)/dx)
    return f 
