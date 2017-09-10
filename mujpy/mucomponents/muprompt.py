from numpy import sqrt, pi, exp
from scipy.special import erf

def muprompt(x, a, x0, dx, ak1, ak2): 
    # fit function for a PSI prompt, sluggish python version
    # data  first row is bin number, second row optional is histo content, third is error on data
    # a gaussian peak coincident with the edge betwee two plateaus (a constant + an erf)
    # par contains peak height, peak center, peak std width, first plateau, second plateau
    f = a/(sqrt(2.*pi)*dx) * exp(-.5*((x-x0)/dx)**2)+ak2/2.+ak1+ak2/2.*erf((x-x0)/sqrt(2.)/dx)
    return f 
