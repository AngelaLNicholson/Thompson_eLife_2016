'''
140330, Boris Zinshteyn and Mary K. Thompson

Purpose: Converts a .csv file output by the BioComp gradient master software, and automatically detects peaks and integrates areas between minima
Input:
    <CSV file> - the input data
    <outputprefix> - the location and base filename of the output
    <skip> - # of mm at start of gradient to ignore
Output:
    a pdf of the polysome plot with integration values
'''

import sys, matplotlib, numpy
matplotlib.use('agg')
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
import os

'''
Colorblind safe colors from Bang Wong, Nature Methods 8. 441 (2011)
'''
black = (0,0,0)
orange = (230/255.0,159/255.0,0)
skyBlue = (86/255.0,180/255.0,233/255.0)
bluishGreen = (0,158/255.0,115/255.0)
yellow = (240/255.0,228/255.0,66/255.0)
blue = (0,114/255.0,178/255.0)
vermillion = (213/255.0,94/255.0,0)
reddishPurple = (204/255.0,121/255.0,167/255.0)
colors = [black, orange, skyBlue, bluishGreen, vermillion, blue, reddishPurple, yellow]

def derivative(x, y, windowHalfSize = 30):
    """
    input:
        x - an array of x values
        y - an array of y values f(x)
    output:
        fpX - returns an array of the approximate derivative at each position i.
            as approximated by a linear fit
    """
    fpX=[]
    for i in range(len(x)):
        if i<windowHalfSize:
            dx=x[0:i+windowHalfSize]
            dy=y[0:i+windowHalfSize]
        elif i>=(len(x)-windowHalfSize):
            dx=x[i-windowHalfSize:i]
            dy=y[i-windowHalfSize:i]
        else:
            dx=x[i-windowHalfSize:i+windowHalfSize]
            dy=y[i-windowHalfSize:i+windowHalfSize]
        slope, b=numpy.polyfit(dx, dy, 1)
        fpX.append(slope)
    return fpX

def findExtremes(x, y, firstDerivative, secondDerivative, d):
    """
    purpose:
        find all minima in y(x), by looking for a sign flip in the first derivative, then checking the sign of the 2nd derivative to see if it's a max or min.
    input:
        x - x values
        y - y(x) values
        firstDerivative - the first derivative of y(x)
        secondDerivative - the second derivative of y(x)
        skip - how many mm at beginning of data to ignore (to not get the top of the gradient)
        mode- 'force' or 'natural'. Force will force the assignment of a minima between 34 and 38. Use for cases where 80S is overlapping with disome so
        don't have baseline separation
    output:
        mins - indices of all minima
        maxs - indices of all maxima
    """
    
    skip= d['skip']
    mode= d['mode']
    mbounds= d['mbounds']
    
    mins=[]
    maxs=[]
    i=1
    lastExtreme=0
    minSeperation=2
    if mode=='force':
        while i < len(x)-1:
            if mbounds[0]<x[i]<mbounds[1]: #special case to help it find the min on the side of the large 80S peak. Adjust as necessary
                if x[i]>=skip and x[i]-lastExtreme>=minSeperation:#this requires that 1 value be positive, and one negative, thus there must be a zero in there somewhere
                    if secondDerivative[i]>0: #concave up, this is a minimum
                        mins.append(i)
                    elif secondDerivative[i]<0: #concave down, this is a maximum
                        maxs.append(i)
                    lastExtreme=x[i] #Used to impose minimum distance between extremes
                    i+=1 #this prevents 2 neighboring datapoints from being called as peaks when then change happened between them
                i+=1
            else:
                if x[i]>=skip and x[i]-lastExtreme>=minSeperation and firstDerivative[i-1]*firstDerivative[i+1] <=0:#this requires that 1 value be positive, and one negative, thus there must be a zero in there somewhere
                    if secondDerivative[i]>0: #concave up, this is a minimum
                        mins.append(i)
                    elif secondDerivative[i]<0: #concave down, this is a maximum
                        maxs.append(i)
                    lastExtreme=x[i] #Used to impose minimum distance between extremes
                    i+=1 #this prevents 2 neighboring datapoints from being called as peaks when then change happened between them
                i+=1
    elif mode=='natural':
        while i < len(x)-1:
            if x[i]>=skip and x[i]-lastExtreme>=minSeperation and firstDerivative[i-1]*firstDerivative[i+1] <=0:#this requires that 1 value be positive, and one negative, thus there must be a zero in there somewhere
                if secondDerivative[i]>0: #concave up, this is a minimum
                    mins.append(i)
                elif secondDerivative[i]<0: #concave down, this is a maximum
                    maxs.append(i)
                lastExtreme=x[i] #Used to impose minimum distance between extremes
                i+=1 #this prevents 2 neighboring datapoints from being called as peaks when then change happened between them
            i+=1
                
    return mins, maxs

def integrate(x, y, index1, index2, m, b):
    """
    Purpose:
        integrates from x[index1] to x[index2]
    output:
        volume - the area under the curve from x[index1] (inclusive) to x[index2](exclusive)
    
    finds the baseline for each i using the line with slope m and intercept b
    
    """
    volume=0
    for i in range(index1, index2-1):
        baseline=evalLine(m, b, index1)
        yVal=numpy.mean(y[i:i+1])-baseline
        dX=x[i+1]-x[i]
        dV=yVal*dX
        volume+=dV
    return volume

def plotAllCurves(dataTraces, firstPeaks, order, outPrefix):
    fig = plt.figure(figsize = (8,8))
    plot = fig.add_subplot(111)
    i = 0
    referencePeak = firstPeaks[order[0]]
    for dataset in order:
        offset = referencePeak-firstPeaks[dataset]
        plot.plot(numpy.array(dataTraces[dataset][0])+offset, dataTraces[dataset][1], label = dataset, color = colors[i], lw =2)
        i+=1
    plot.legend(frameon = False)
    plt.xlabel('distance(mm)')
    plt.ylabel('absorbance (OD254)')
    plt.savefig(outPrefix+'combined_plots.pdf', transparent='True', format='pdf')
    plt.clf()

def printAllRatios(ratios, datasetOrder, ratioOrder, outPrefix):
    f = open(outPrefix+'_ratios.txt', 'w')
    header = '\t'.join(['dataset']+ratioOrder) +'\n'
    f.write(header)
    for dataset in datasetOrder:
        line = '\t'.join([dataset]+[str(ratios[dataset][ratio]) for ratio in ratioOrder])+'\n'
        f.write(line)
    f.close()

def parseParams(paramfile, indir):
    
    '''Parse a parameter file containing the input files and the parameters for each run'''
    
    csvfiles=[]
    params={}
    
    f= open(paramfile, 'r')
    for line in f:
        fields= line.strip('\n').split()
        filename= fields[0]
        fullname=os.path.join(indir, filename)
        csvfiles.append(fullname)
        params[fullname]={}
        outname= fields[1]
        mode_fields= fields[2].split('=')
        mode= mode_fields[0]
        if len(mode_fields)>1:
            mbounds= mode_fields[1][1:-1].split(',')
            mbounds= (float(mbounds[0]), float(mbounds[1]))
        else:
            mbounds=''
        skip= fields[3]
        params[fullname]['outname']= outname
        params[fullname]['mode']= mode
        params[fullname]['mbounds']= mbounds
        params[fullname]['skip']= float(skip)
        
    f.close()
    return csvfiles, params

def findBaseline(params, mins, absorbances, distances):
    '''Given the absorbances, find a baseline that goes through the minimum in the first half (typically between 40S and 60S)
    and the minimum in the second half (typically at end of the trace)'''
    
    #get the index corresponding to the start (skip)
    target= params[csvFile]['skip']
    startI=''
    for i in range(0, len(distances)):
        if distances[i]>target:
            startI=i
            break
        
    min1= min(absorbances[startI:len(absorbances)/2]) #find min in first half of the list
    x1=absorbances.index(min1)
    min2= absorbances[mins[-1]] #find last min
    x2= absorbances.index(min2)
    if distances[x2]<65: #probably didn't find a min at the end
        min2= absorbances[-1]
        x2= len(absorbances)-1
    #fit the line:
    m=(min2-min1)/(x2-x1) #this will be the slope (Abs/index)
    #find intercept using the first point
    b= min1-m*x1
    return m, b, x2

def evalLine(m, b, x):
    '''Return the y-value at a given x value and for a line with intercept b and slope m'''
    y= m*x+b
    return y
    
outputPrefix, indir, paramfile = sys.argv[1:4]
#skip = int(skip)
if not os.path.exists(outputPrefix):
    os.makedirs(outputPrefix)
    
csvFiles, params= parseParams(paramfile, indir)
datasetOrder =[]
dataTraces = {} #Maps csvFileName to tuple of X and Y data
firstPeaks = {} #Sample name to x pos of first peak
ratios = {} #Sample name to essential ratios
for csvFile in csvFiles:
    print 'csvfile', csvFile
    f=open(csvFile)
    #First skip all the headers to get to the data, take down experiment info meanwhile
    expInfo={}
    lines=f.readlines()
    i=0
    while not lines[i].startswith('Distance(mm)'):
        if not lines[i].split(',')[0].strip()=='':
            ll=lines[i].split(',')[0].strip().split(':')
            if len(ll)>1:
                expInfo[ll[0]]=ll[1]
        i+=1
        
    dataset= os.path.basename(expInfo['Gradient Profiler Text Data File'])
    datasetOrder.append(dataset)
    
    #Get OD254 data
    distances=[]
    absorbances=[]
    fracNums=[]
    fracVolumes=[]
    i+=1
    while i < len(lines):
        ll=lines[i].strip().strip(',').split(',')
        if float(ll[0])>0:#the datafile starts with a bunch of entries at distance zero, this presents problems for the derivative, and so is removed
            distances.append(float(ll[0]))
            absorbances.append(float(ll[1]))
            if len(ll)>=4:
                fracNums.append(float(ll[2]))
                fracVolumes.append(float(ll[3]))
            else:
                fracNums.append(None)
                fracVolumes.append(None)
        i+=1
    f.close()
    dataTraces[dataset] = (distances, absorbances)
   
    
    firstDir=derivative(distances, absorbances)
    secondDir=derivative(distances, firstDir)
    mins, maxs=findExtremes(distances, absorbances, firstDir, secondDir, params[csvFile])
    
    #turn baseline into a linear function now and then will just solve for the corresponding value during integration:
    m,b, finalI= findBaseline(params, mins, absorbances, distances)
    
    vol40s = integrate(distances, absorbances, mins[0], mins[1], m, b)
    vol60s = integrate(distances, absorbances, mins[1], mins[2], m, b)
    vol80s = integrate(distances, absorbances, mins[2], mins[3], m, b)
    volPoly = integrate(distances, absorbances, mins[3], finalI, m, b)

    #change this here to try integrating for 32 mm from the disome min so that all of them will be more equalized and won't have the issue of very weird things at the bottom
    ratios[dataset] = {'40s':vol40s, '60s':vol60s, '80s':vol80s, 'poly':volPoly, 'poly_40':volPoly/vol40s,
                       'pol_60':volPoly/vol60s, 'pol_80':volPoly/vol80s}
    ratioOrder = ['40s','60s','80s','poly','poly_40','pol_60','pol_80']
    
    minis= [distances[i] for i in mins]
    
    peak40s = maxs[0]
    firstPeaks[dataset] = distances[peak40s]
    peak60s = maxs[1]
    peak80s = maxs[2]
    #peakDisome = maxs[3]
    plt.plot(distances, absorbances, color='Black', label='f(x)')
    #plt.plot(distances, firstDir, color='red', label='f\'(x)')
    #plt.plot(distances, secondDir, color='green', label='f\'\'(x)')
    for i in mins:
        plt.axvline(distances[i], color='blue', linestyle='--')
    for i in maxs:
        plt.axvline(distances[i], color='pink', linestyle='--')
        
    #plot the baseline according to the slope of the line:
    #This plots a y=x line
    x_points=[]
    y_points=[]
    for i in range(0, len(distances)):
        x_points.append(distances[i])
        y1= evalLine(m, b, i)
        y_points.append(y1)
        
    plt.plot(x_points, y_points, linestyle='-', color='r')
    
    plt.axvline(params[csvFile]['skip'], color='black', linestyle='--')
    plt.annotate('40s=%f' % (vol40s), xy=(distances[peak40s], absorbances[peak40s]+0.1))
    plt.annotate('60s=%f' % (vol60s), xy=(distances[peak60s], absorbances[peak60s]+0.1))
    plt.annotate('80s=%f' % (vol80s), xy=(distances[peak80s], absorbances[peak80s]+0.1))
    #plt.annotate('polysome=%f' % (volPoly), xy=(distances[peakDisome], absorbances[peakDisome]+0.1))
    ratio=volPoly/vol80s
    #plt.annotate('polysome/80s=%f' % (ratio), xy=(distances[peakDisome], absorbances[peakDisome]+0.05))
    plt.title(expInfo['Gradient Profiler Text Data File'].split('.')[0]+' '+expInfo['Date'])
    plt.xlabel('distance(mm)')
    plt.ylabel('absorbance (OD254)')
    #lg=plt.legend()
    #lg.get_frame().set_linewidth(0)
    plt.savefig(outputPrefix+dataset+'.pdf', transparent='True', format='pdf')
    plt.clf()
plotAllCurves(dataTraces, firstPeaks, datasetOrder, outputPrefix)
printAllRatios(ratios,datasetOrder, ratioOrder, outputPrefix)