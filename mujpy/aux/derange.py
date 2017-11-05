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
    try:  
        try:
            values = [float(x) for x in string.split(',')]
        except:
            values = [float(x) for x in string.split('')]
        if len(values)==4:
            return values[0],values[1],values[2],values[3] # start, stop, last, pack
        elif len(values)==3:
            return values[0],values[1],values[2] # start, stop, pack
        elif len(values)==2:
            return values[0],values[1] # start, stop
    except:
        return -1,-1 

