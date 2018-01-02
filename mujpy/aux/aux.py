########################
# FFT AUTO PHASE METHODS
########################
def autops(data, fn, p0=0.0, p1=0.0):
    """
    Automated phase correction from NMRglue by https://github.com/jjhelmus
    These functions provide support for automatic phasing of NMR data. 
    ----------
    Automatic linear phase correction
    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    fn : str or function
        Algorithm to use for phase scoring. Built in functions can be
        specified by one of the following strings: "acme", "peak_minima"
    p0 : float
        Initial zero order phase in degrees.
    p1 : float
        Initial first order phase in degrees.
    Returns
    -------
    ndata : ndarray
        Phased NMR data.
    """
    import numpy as np
    import scipy.optimize

    if not callable(fn):
        fn = {
            'peak_minima': _ps_peak_minima_score,
            'acme': _ps_acme_score,
        }[fn]
    opt = [p0, p1]
    opt = scipy.optimize.fmin(fn, x0=opt, args=(data, ))
    return ps(data, p0=opt[0], p1=opt[1]), opt[0], opt[1]


def _ps_acme_score(ph, data):
    """
    Phase correction using ACME algorithm by Chen Li et al.
    Journal of Magnetic Resonance 158 (2002) 164-168
    Parameters
    ----------
    pd : tuple
        Current p0 and p1 values
    data : ndarray
        Array of NMR data.
    Returns
    -------
    score : float
        Value of the objective function (phase score)
    """
    import numpy as np

    stepsize = 1

    phc0, phc1 = ph

    s0 = ps(data, p0=phc0, p1=phc1)
    data = np.real(s0)

    # Calculation of first derivatives
    ds1 = np.abs((data[1:]-data[:-1]) / (stepsize*2))
    p1 = ds1 / np.sum(ds1)

    # Calculation of entropy
    p1[p1 == 0] = 1

    h1 = -p1 * np.log(p1)
    h1s = np.sum(h1)

    # Calculation of penalty
    pfun = 0.0
    as_ = data - np.abs(data)
    sumas = np.sum(as_)

    if sumas < 0:
        pfun = pfun + np.sum((as_/2) ** 2)

    p = 1000 * pfun

    return h1s + p


def _ps_peak_minima_score(ph, data):
    """
    Phase correction using simple minima-minimisation around highest peak
    This is a naive approach but is quick and often achieves reasonable
    results.  The optimisation is performed by finding the highest peak in the
    spectra (e.g. TMSP) and then attempting to reduce minima surrounding it.
    Parameters
    ----------
    pd : tuple
        Current p0 and p1 values
    data : ndarray
        Array of NMR data.
    Returns
    -------
    score : float
        Value of the objective function (phase score)
    """

    phc0, phc1 = ph

    s0 = ps(data, p0=phc0, p1=phc1)
    data = np.real(s0)

    i = np.argmax(data)
    mina = np.min(data[i-100:i])
    minb = np.min(data[i:i+100])

    return np.abs(mina - minb)

def ps(data, p0=0.0, p1=0.0, inv=False):
    """
    Linear phase correction
    Parameters
    ----------
    data : ndarray
        Array of NMR data.
    p0 : float
        Zero order phase in degrees.
    p1 : float
        First order phase in degrees.
    inv : bool, optional
        True for inverse phase correction
    Returns
    -------
    ndata : ndarray
        Phased NMR data.
    """
    import numpy as np

    p0 = p0 * np.pi / 180.  # convert to radians
    p1 = p1 * np.pi / 180.
    size = data.shape[-1]
    apod = np.exp(1.0j * (p0 + (p1 * np.arange(size) / size))).astype(data.dtype)
    if inv:
        apod = 1 / apod
    return apod * data

##############
# MUGUI AUX
##############

def derange(string):
    '''
    derange(string) 
    reads string 
    assuming 2, 3 or 4 csv values, integers or floats,
    or 2, 3 or 4 space separated values, integers or floats
    returns 2, 3 or 4 floats, or 
    two negative values, if the string does not 
    contain commas or spaces
    '''
#    print('In derange')
    try:  
        try:
            values = [float(x) for x in string.split(',')]
        except:
            values = [float(x) for x in string.split('')]
        if len(values)==5:
            if values[3]<values[1] or values[4]>values[3] or values[2]>values[1] or values[1]<values[0] or sum(n<0 for n in values)>0:
                return -1,-1
            return values[0],values[1],values[2],values[3],values[4] # start, stop, packe, last, packl
        if len(values)==4:
            return values[0],values[1],1,values[2],values[3] # start, stop,packe, last, packl
        elif len(values)==3:
            return values[0],values[1],values[2] # start, stop, pack
        elif len(values)==2:
            return values[0],values[1] # start, stop
    except:
        return -1,-1 

def derange_int(string):
    '''
    derange(string) 
    reads string 
    assuming 2, 3 or 4 csv values, integers or floats,
    or 2, 3 or 4 space separated values, integers or floats
    returns 2, 3 or 4 floats, or 
    two negative values, if the string does not 
    contain commas or spaces
    '''
#    print('In derange_int')
    try:  
        try:
            values = [int(x) for x in string.split(',')]
        except:
            values = [int(x) for x in string.split('')]
        if len(values)==5:
            if values[3]<values[1] or values[4]>values[3] or values[2]>values[1] or values[1]<values[0] or sum(n<0 for n in values)>0:
                return -1,-1
            return values[0],values[1],values[2],values[3],values[4] # start, stop, packe, last, packl
        if len(values)==4:
            if values[3]<values[1] or values[2]>values[1] or values[1]<values[0] or sum(n<0 for n in values)>0:
                return -1,-1
            return values[0],values[1],1,values[2],values[3] # start, stop, packe, last, packl
        elif len(values)==3:
            return values[0],values[1],values[2] # start, stop, pack
        elif len(values)==2:
            return values[0],values[1] # start, stop
    except:
        return -1,-1 

def derun(string):
    '''
    parses string, producing a list of runs; 
    expects comma separated items
    looks for 'l','l:m','l+n+m' 
    where l, m, n are integers
    rejects all other characters
    returns a list of lists of integer
    '''
    s = []
    try:
    # systematic str(int(b[])) to check that b[] ARE integers
        for b in string.split(','): # csv
            kcolon = b.find(':') # ':' and '+' are mutually exclusive
            kplus = b.find('+')
            if kcolon>0: # value might be a range
                for j in range(int(b[:kcolon]),int(b[kcolon+1:])+1):
                    s.append([str(j)]) # strings
            elif kplus>0:
                ss = []
                k0 = 0
                while kplus>0: # str(int(b[]))
                    ss.append(int(b[k0:kplus])) 
                    k0 = kplus+1
                    kplus = b.find('+',k0)
                ss.append(int(b[k0:]))
                s.append([str(q) for q in ss])
            else:# value should be an integer
                int(b) # produces an Error if b is not an integer
                s.append([b]) # added as a string
    except Exception as e:
        s = []
        errmsg = e
    if 'errmsg' not in locals():
        errmsg = None
    return s, errmsg

def findall(p, s):
    '''Yields all the positions of
    the pattern p in the string s.'''
    i = s.find(p)
    while i != -1:
        yield i
        i = s.find(p, i+1)

def find_nth(haystack, needle, n):
           start = haystack.rfind(needle)
           while start >= 0 and n > 1:
               start = haystack.rfind(needle, 1, start-1)
               n -= 1
           return start

def get_grouping(groupcsv):
    """
    name = 'forward' or 'backward'
    grouping(name) is an np.array wth detector indices
    group.value[k] for k=0,1 is a shorthand csv like '1:3,5' or '1,3,5' etc.
    index is present mugui.mainwindow.selected_index
    out is mugui._output_ for error messages
    returns
      grouping, group, index
         group and index are changed only in case of errors
    """
    import numpy as np

    # two shorthands: either a list, comma separated, such as 1,3,5,6 
    # or a pair of integers, separated by a colon, such as 1:3 = 1,2,3 
    # only one column is allowed, but 1, 3, 5 , 7:9 = 1, 3, 5, 7, 8, 9 
    # or 1:3,5,7 = 1,2,3,5,7  are also valid
    #       get the shorthand from the gui Text 
    groupcsv = groupcsv.replace('.',',') # can only be a mistake: '.' means ','
    try:
        if groupcsv.find(':')==-1:
            # colon not found, csv only
            grouping = np.array([int(s) for s in groupcsv.split(',')])
        else:
            # colon found                 
            if groupcsv.find(',')+groupcsv.find(':')==-2:
                grouping = np.array([int(groupcsv)])
            elif groupcsv.find(',')+1: # True if found, False if not found (not found yields -1)    
                firstcomma = groupcsv.index(',')
                lastcomma = groupcsv.rindex(',')
                if firstcomma < groupcsv.find(':'): # first read csv then range
                    partial = np.array([int(s) for s in groupcsv[:lastcomma].split(',')])
                    fst = int(groupcsv[lastcomma:grouping.find(':')])
                    lst = int(groupcsv[groupcsv.find(':')+1:])
                    grouping[name] = np.concatenate((partial,arange(fst,lst+1,dtype=int)))
                else: # first read range then csv
                    partial = np.array([int(s) for s in groupcsv[:lastcomma].split(',')])
                    fst = int(groupcsv[:groupcsv.find(':')])
                    lst = int(groupcsv[groupcsv.find(':')+1:firstcomma])
                    grouping = np.concatenate((np.arange(fst,lst+1,dtype=int),partial))
            else: # only range
                fst = int(groupcsv[:groupcsv.find(':')])
                lst = int(groupcsv[groupcsv.find(':')+1:])
                grouping = np.arange(fst,lst+1,dtype=int)
            grouping -=1 # python starts at 0
    except:
        grouping = np.array([-1])
        
    return grouping

def get_title(run):
    '''
    form standard psi title
    '''
    return '{} {} {} {}'.format(run.get_sample(),run.get_field(),run.get_orient(),run.get_temp())    

def muvalid(string,out,index):
    '''
    parse function CHECK WITH MUCOMPONENT, THAT USES A DIFFERENT SCHEME
    accepted functions are RHS of agebraic expressions of parameters p[i], i=0...ntot  
    '''
    import re

    pattern = re.compile(r"\p\[(\d+)\]") # find all patterns p[*] where * is digits
    test = pattern.sub(r"a",string) # substitute "a" to "p[*]" in s
    #           strindices = pattern.findall(string)
    #           indices = [int(strindices[k]) for k in range(len(strindices))] # in internal parameter list
    #           mindices = ... # produce the equivalent minuit indices  
    valid = True
    message = ''
    try: 
        safetry(test) # should select only safe use (although such a thing does not exist!)
    except Exception as e:
        with out:
            print('function: {}. Tested: {}. Wrong or not allowed syntax: {}'.format(string,test,e))
        index = 3
        valid = False
    return valid

def muvaluid(string):
    '''
    Run suite fits: muvaluid returns True/False
    checks the syntax for string function 
    corresponding to flag='l'. Meant for pars
    displaying large changes across the run suite,
    requiring different migrad start guesses.
    '''
    # string syntax: e.g. "0.2*3,2.*4,20."
    # means that for the first 3 runs value = 0.2,
    #            for the next 4 runs value = 2.0
    #            from the 8th run on value = 20.0
    try:
        value_times_list = string.split(',')
        last = value_times_list.pop()
        for value_times in value_times_list:
            value,times = value_times.split('*')
            dum, dum = float(value),int(times)
        dum = float(last)
        return True
    except:
        return False

def muvalue(lrun,string):
    '''
    Run suite fits: muvalue returns the value 
    for the nint-th parameter of the lrun-th run
    according to string (corresponding flag='l').
    Large parameter change across the run suite
    requires different migrad start guesses.
    '''
    # string syntax: e.g. "0.2*3,2.*4,20."
    # means that for the first 3 runs value = 0.2,
    #            for the next 4 runs value = 2.0
    #            from the 8th run on value = 20.0

    value = []
    for value_times in string.split(','):
        try:
            value,times = value_time.split('*') 
            for k in range(int(times)):
                value.append(float(value))
        except:
            for k in range(len(value),lrun):
                value.append(float(value_time))
    return value[lrun]

def muzeropad(runs,out):
    '''
    muzeropad(runs,out)
        runs is a string containing the run number
    utility of the suite tab, not a method
    future:
    1) determine how many leading zeros for padding
       read a file from data dir
       check number of digits before '.'
       count number of digits in run
       zero pad
    now:
    0) zeroth version pads a fixed number of zeros to 4 digits
    out in mugui is _output_, Textarea for printing messages
    '''
    zeros='0000'
    if len(runs)<len(zeros):
        return zeros[:len(zeros)-len(runs)]+runs
    elif len(runs)==len(zeros):
        return runs
    else:
        with out:
            print('Too long run number!')
        return []

def norun_msg(out):
    with out:
         print('No run loaded yet! Load one first (select suite tab).')


def plotile(x,xdim=0,offset=0):
    '''
    xt = plotile(x,xdim,xoffset=0)
        e.g. x.shape = (1,1000) y.shape = (4,1000)
        xt = plotile(x,4)
        yt = plotile(y,offset=0.1) 
    '''
    # x is an array(x.shape[0],x.shape[1])
    # xoffset is a step offset
    # xdim = x.shape[0] if xdim == 0 else xdim
    # each row is shifted by xoffset*n, where n is the index of the row  
    # 
    # 
    from copy import deepcopy
    from numpy import tile, arange
    xt = deepcopy(x)
    if xdim != 0: # x is a 1D array, must be tiled to xdim
        xt = tile(xt,(int(xdim),1))
    if offset != 0:
        xt += tile(offset*arange(xt.shape[0]),(x.shape[1],1)).transpose()
    return xt

def rebin(x,y,strstp,pack,e=None):
    '''
    x,y[,e] are 2D arrays to be rebinned
    pack is the rebinning factor, e.g:
        xb = array([x[0:pack].sum()/pack])
    strstp = [start,stop] is a list of indices
        rebinning is done on slices of x,y[,e]
        such as x[start:stop]

    use either 
    xb,yb = rebin(x,y,strstp,pack)
    or
    xb,yb,eyb = rebin(x,y,strstp,pack,ey) # the 5th is y error
    Works also with pack = 1

    * Now works for 1d array x and 2d ndarrays y, ey
    xb, yb, eb are 2D arrays, e.g. xb.shape = (1,25000), yb.shape = (nruns,25000) 
    '''
    from numpy import floor, sqrt,empty, array
    start,stop = strstp
    if pack==1: # no rebinning, just a slice of 2D arrays
        xx = x[:,start:stop] 
        yy = y[:,start:stop]
        # print('in rebin: shape xx {}, yy {}'.format(xx.shape,yy.shape)) 
        if e is None:
            return xx,yy
        else:
            ee = e[:,start:stop]
            return xx,yy,ee
    else:
        m = int(floor((stop-start)/pack)) # length of rebinned xb
        mn = m*pack # length of x slice 
        xx =x[:,start:start+mn] # slice of the first 2D array
        xx = xx.reshape(m,pack) # temporaty 2d array
        xb = array([xx.sum(1)/pack]) # rebinned first ndarray
        nruns = y.shape[0] # number of runs
        yb = empty((nruns,m))
        if e is not None:
            eb = empty((nruns,m))
        for k in range(nruns): # each row is a run
            yy = y[k][start:start+mn]  # slice row
            yy = yy.reshape(m,pack)  # temporaty 2d
            yb[k] = yy.sum(1)/pack # rebinned row
            if e is not None:
                ey = e[k][start:start+mn]   # slSice row
                ey = ey.reshape(m,pack)  # temporaty 2d
                eb[k] = sqrt((ey**2).sum(1))/pack  # rebinned row
        if e is not None:
            return xb,yb,eb
        else:
            return xb,yb

def safetry(string):
    from math import acos,asin,atan,atan2,ceil,cos,cosh,degrees,e,exp,floor,log,log10,pi,pow,radians,sin,sinh,sqrt,tan,tanh
    safe_list = ['a','acos', 'asin', 'atan', 'atan2', 'ceil', 'cos', 'cosh', 'degrees', 'e', 
                 'exp', 'floor', 'log', 'log10', 'pi', 'pow', 'radians', 'sin', 'sinh', 'sqrt', 'tan', 'tanh']
    # 	use the list to filter the local namespace
    a = 0.3
    safe_dict={}
    for k in safe_list:
        safe_dict[k]=locals().get(k)
    #    print(safe_dict[k])
    return eval(string,{"__builtins__":None},safe_dict)

def set_bar(n,b):
    '''
    service to animate histograms
    e.g. in the fit tab
    extracted from matplotlib animate 
    histogram example
    '''
    from numpy import array, zeros, ones
    import matplotlib.path as path

    # get the corners of the rectangles for the histogram
    left = array(b[:-1])
    right = array(b[1:])
    bottom = zeros(len(left))
    top = bottom + n
    nrects = len(left)

    # here comes the tricky part -- we have to set up the vertex and path
    # codes arrays using moveto, lineto and closepoly

    # for each rect: 1 for the MOVETO, 3 for the LINETO, 1 for the
    # CLOSEPOLY; the vert for the closepoly is ignored but we still need
    # it to keep the codes aligned with the vertices
    nverts = nrects*(1 + 3 + 1)
    verts = zeros((nverts, 2))
    codes = ones(nverts, int) * path.Path.LINETO
    codes[0::5] = path.Path.MOVETO
    codes[4::5] = path.Path.CLOSEPOLY
    verts[0::5, 0] = left
    verts[0::5, 1] = bottom
    verts[1::5, 0] = left
    verts[1::5, 1] = top
    verts[2::5, 0] = right
    verts[2::5, 1] = top
    verts[3::5, 0] = right
    verts[3::5, 1] = bottom
    xlim = [left[0], right[-1]]
    return verts, codes, bottom, xlim

def translate(nint,lmin,function):
    string = function[nint].value
    # search for integers between '[' and ']'
    start = [i+1 for i in  findall('[',string)]
    stop = [i for i in  findall(']',string)]
    nints = [string[i:j] for (i,j) in zip(start,stop)]
    nmins = [lmin[int(string[i:j])] for (i,j) in zip(start,stop)]
    for lstr,m in zip(nints,nmins):
        string = string.replace(lstr,str(m))
    return string

def value_error(value,error):
    '''
    value_error(v,e)
    returns a string of the format v(e) 
    '''
    from numpy import floor, log10
    exponent = int(floor(log10(error)))  
    most_significant = int(round(error/10**exponent))
    if most_significant>9:
        exponent += 1
        most_significant=1
    exponent = -exponent if exponent<0 else 0
    form = '"{:.'
    form += '{}'.format(exponent)
    form += 'f}({})".format(value,most_significant)'
    return eval(form)
