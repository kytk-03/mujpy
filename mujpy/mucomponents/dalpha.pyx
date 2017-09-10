# We want to speed things up with Cython
%load_ext Cython

%%cython
# damu now written in Cython
# The %%cython IPython does the following:
# * Call Cython to generate C code for a Python C extension.
# * Compile it into a Python C extension (a shared library)
# * Load it into the current namespace
# If you don't understand these things, don't worry, it basically means:
# * Get full-metal speed easily
cimport cython
import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.float
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.float_t DTYPE_t
@cython.binding(True) # IMPORTANT: this tells Cython to dump the function signature    

def dalpha_cython(np.array x, double dalpha): 
    # fit function for a PSI prompt, 
    # data  first row is bin number, second row optional is histo content, third is error on data
    # a gaussian peak coincident with the edge betwee two plateaus (a constant + an erf)
    # par contains peak height, peak center, peak std width, first plateau, second plateau
    return np.array np.zeros(x.len())

