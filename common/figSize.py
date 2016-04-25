def convertFigSize(keyword, pagetype):
    '''
    will return a tuple of (x-inches, y-inches) to give to the figsize argument based on keywords specifying the desired figure size.
    We will assume for now that the figures should all be square.
    keywords=['1','2','3', '4'] #this corresponds to the number of figures there will be across the page. i.e. 4= this figure will be a quarter page.
    pagetypes=['full', 'half']
    '''
    conversion=float(1)/25.4 #multiply this by the mm values to get the number of inches that we need
    num= float(keyword) #the number of plots there will be across the page
    
    if pagetype=='full':
        max_width= 183
        max_height=247 #these are Nature-specific values in mm, change if necessary
    else:
        max_width=89
        max_height=247
        
    #sizeD={} #calculate necessary width for 1x, 0.33x, 0.5x or 0.25x width figure
    margin=3 #change as necessary, this amount of space on left and right to enable placement of the a, b, c.. labels
    spacing=1 #this is the factor between the figures to expect (multiply marging times this factor)
        
    availw= max_width-2*margin
    deadw= (num-1)*(margin*spacing)
    availw= availw-deadw
    wperplot= availw/num
    
    w_in= wperplot*conversion
    
    return (w_in, w_in)

def getLinesize(keyword, pagetype):
    '''
    As above, but get recommended linesize for this size of figure
    For now, it will return as follows:
    (1, 'full')= 0.75
    (2, 'full')= 0.5
    (3, 'full')= 0.5
    (4, 'full')= 0.25
    (1, 'half')=0.5
    '''
    keyword= int(keyword)
    
    if pagetype=='full':
        if keyword==1:
            return llinesize
        elif (keyword==2) or (keyword==3):
            return mlinesize
        elif keyword==4:
            return slinesize

    elif pagetype=='half':
        if keyword==1:
            return mlinesize
        elif (keyword==2) or (keyword==3):
            return slinesize
        elif keyword==4:
            return slinesize
