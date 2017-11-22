##########################
# REBIN
##########################
def rebin(x,y,pack,e=None):
    '''
    use either 
    xb,yb = rebin(x,y,pack)
    or
    xb,yb,eyb = rebin(x,y,pack,ey) # the 4rd is an error
    Works also with pack = 1
    '''
    from numpy import floor, sqrt
    if pack==1:
        if e is None:
            return x,y
        else:
            return x,y,e
    else:
        m = int(floor(len(x)/pack))
        mn = m*pack
        xx = x[:mn]
        xx = xx.reshape(m,pack)
        yy = y[:mn]
        yy = yy.reshape(m,pack)
        xb = xx.sum(1)/pack
        yb = yy.sum(1)/pack
        if e is not None:
            ey = e[:mn]
            ey = ey.reshape(m,pack)
            eb = sqrt((ey**2).sum(1))/pack
            return xb,yb,eb
        else:
            return xb,yb

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

