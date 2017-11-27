def muzeropad(runs,out):
    '''
    muzeropad(runs,out)
    utility of the suite tab, not a method
    future:
    1) look for '+' and separate runs to be added
    2) determine how many leading zeros for padding
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

def rebin(x,y,strstp,pack,e=None):
    '''
    use either 
    xb,yb = rebin(x,y,strstp,pack)
    or
    xb,yb,eyb = rebin(x,y,strstp,pack,ey) # the 5th is an error
    Works also with pack = 1
    Must allow for each run being a row in 2d ndarrays y, ey

    '''
    from numpy import floor, sqrt,empty, array
    start,stop = strstp
    ydim = y.ndim
    if pack==1:
        xx = array(x[start:stop])
        yy = array(y[start:stop]) if ydim==1 else array(y[:,start:stop])
        # print('in rebin: shape xx {}, yy {}'.format(xx.shape,yy.shape)) 
        if e is None: # no rebinning, just a slice
            return xx,yy
        else:
            ee = array(e[start:stop]) if ydim==1 else array(e[:,start:stop])
            return xx,yy,ee
    else:
        m = int(floor((stop-start)/pack)) # length of rebinned xb
        mn = m*pack # length of x slice 
        xx = x[start:start+mn] # slice of the first 1-dimensional ndarray
        xx = xx.reshape(m,pack) # temporaty 2d array
        xb = xx.sum(1)/pack # rebinned first ndarray
        if ydim == 1: # single 
            yy = y[start:start+mn]  # slice 1-d
            yy = yy.reshape(m,pack)  # temporaty 2d 
            yb = yy.sum(1)/pack   # rebinned 1-d
            if e is not None:
                ey = e[start:start+mn]  # slice 1-d
                ey = ey.reshape(m,pack)  # temporaty 2d
                eb = sqrt((ey**2).sum(1))/pack  # rebinned 1-d
        else: # suite
            yb = empty((ydim,m))
            if e is not None:
                eb = empty((ydim,m))
            for k in range(y.shape[0]): # each row is a run
                yy = y[k][start:start+mn]  # slice row
                yy = yy.reshape(m,pack)  # temporaty 2d
                yb[k] = yy.sum(1)/pack # rebinned row
                if e is not None:
                    ey = e[k][start:start+mn]   # slice row
                    ey = ey.reshape(m,pack)  # temporaty 2d
                    eb[k] = sqrt((ey**2).sum(1))/pack  # rebinned row
        if e is not None:
            return xb,yb,eb
        else:
            return xb,yb

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


