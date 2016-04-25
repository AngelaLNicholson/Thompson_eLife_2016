#140519 Mary K. Thompson
#purpose: To add plot polysome gradients and add P/M values computed elsewhere.

import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
import argparse
import math
import numpy as np

fontsize= 16

#fontsizes permitted
lfsize=7
mfsize=6
sfsize=5

#linesize permitted
slinesize=0.25
mlinesize=0.5
llinesize=0.75

def convertFigSize(keyword, pagetype):
    '''
    will return a tuple of (x-inches, y-inches) to give to the figsize argument based on keywords specifying the desired figure size.
    We will assume for now that the figures should all be square.
    keywords=['1','2','3', '4'] #this corresponds to the number of figures there will be across the page. i.e. 4= this figure will be a quarter page.
    pagetypes=['full', 'half']
    '''
    conversion=float(1)/25.4 #multiply this by the mm values to get the number of inches that we need
    num= int(keyword) #the number of plots there will be across the page
    
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

def read(file, mark='ON'):

	distance= []
	Abs= []
	Mark= []

	begin= 0

	f= open(file, 'r')
	for line in f:
		if line.startswith('Distance(mm)'):
			#we're at the start of the values
			begin= 1
		elif begin== 1:
			#add to arrays
			fields= line.strip().split(',')
			if fields== ['']:
				continue 
			else:
				distance.append(float(fields[0]))
				Abs.append(float(fields[1]))
				#add the fraction mark if present
				try:
					a= int(fields[2])
					Mark.append((fields[0],fields[2])) #this tells you where the mark occurred
				except:
					continue
		else:
			continue
	f.close()
	return distance, Abs, Mark

def readLabels(labelfile):
	
	labels=[]
	f= open(labelfile, 'r')
	for line in f:
		name= line.strip()
		labels.append(name)
	f.close()
	return labels

def addDead(x, deadvol):
	'''
	adjust data based on the deadvol since it by default doesn't add the deadvol recordings. All of these will be marked with dist=0
	'''
	
	#in order to replace the 0s with a distance, we will assume that the samples are being taken at even rate and assign the next number
	#samples to the samples in the dead volume
	
	num=0
	for i in range(0, len(x)):
		if x[i]==0.0:
			num=i+1#this will set num to the index of the first non-zero point collected
	
	lendead= num
	#reset the first xvalues to the next x values
	x[0:num]= x[num:num+lendead]
	
	#add the deadvol to the rest of the x values
	for i in range(num, len(x)):
		x[i]=x[i]+deadvol
	
	return x

def zeroAbs(Abs, distance):
	'''
	return abs with minimum set to zero, only the minimum between 2 and 70 mm
	'''
	
	#find indices of 2 and 70 mm
	mini=0
	maxi=0
	for i in range(0, len(distance)):
		if distance[i]>2:
			mini= i
			break
	for i in range(0, len(distance)):
		if distance[i]>70:
			maxi=i
			break
	
	minVal=min(Abs[mini:maxi])
	
	#subtract the minVal from all the numbers
	newAbs=[]
	for j in Abs:
		newAbs.append(j-minVal)
	return newAbs

def clip(x, y, lbound, rbound):
	
	#find indices of 2 and 70 mm
	mini=0
	maxi=0
	for i in range(0, len(x)):
		if x[i]>lbound:
			mini= i
			break
	for i in range(0, len(x)):
		if x[i]>rbound:
			maxi=i
			break
		
	#therefore range(mini, maxi) will give you all the points between the min and the max (non-inclusive)
	
	return x[mini:maxi], y[mini:maxi]
	
def plot(allData, outname, args):
	#plot all the traces:
	#if you don't specifiy xdim and ydim, make it square.
	num_plots=len(allData)
	dead_vol= 1.69 #distance due to mm. note that this is actually the distance (in mm) that it doesn't add to the distance b/c it considers it in the dead volume
	
	if args['xdim']=='':
				
		dim= int(math.ceil(math.sqrt(num_plots)))
		
		x_num= dim
		y_num= dim
	else:
		x_num= args['xdim']
		y_num= args['ydim']
				
	leftmost=[]
	for l in range(0, num_plots):
		if l%x_num==0:	
			leftmost.append(l+1)
	
	fig=plt.figure(figsize=convertFigSize(1, 'full'))

	for p in range(0, len(allData)):
		ax = fig.add_subplot(y_num,x_num,p+1, frameon=False)
		if args['include_dead']==True:
			allData[p]['distance']= addDead(allData[p]['distance'], dead_vol)
		if args['zero']==True:
			allData[p]['Abs']= zeroAbs(allData[p]['Abs'], allData[p]['distance'])
			
		ax.plot(allData[p]['distance'],allData[p]['Abs'], 'k', linewidth=args['lw'])
		start, end= ax.get_ylim()
		
		newID= allData[p]['Name']
		plt.text(0.5, 1.15, newID, horizontalalignment='center',verticalalignment='center', transform=ax.transAxes, fontsize=7)
		
		if args['mark']==True:
			for i in allData[p]['Mark']:
	
				fracx= float(i[0])#-dead_vol #the dead_vol is already accounted for by the way we parse the file
				if int(i[1])<1:
					continue
				ax.axvline(fracx)
				plt.text(fracx, args['maxA']*0.9, str(i[1]), horizontalalignment='right', fontsize=5, color= 'r')
		
		#try this from the perspective of only labeling the leftmost ones rather than unlabeling the others:

		if args['noAxes']!=True:
			
			if args['sparseAx']==True: #only label the leftmost plots
				if (p+1) in leftmost:
					plt.ylabel('A254', fontsize=6)
					plt.xlabel('Distance (mm)', fontsize= 6)
					ax.yaxis.set_ticks(np.arange(start, end, 0.1))
					ax.xaxis.set_ticks_position('bottom')
					ax.yaxis.set_ticks_position('left')
					ax.xaxis.set_tick_params(labelsize=6)
					ax.yaxis.set_tick_params(labelsize=6) #changing these from 8 to 7 to see if it makes a difference when opening in Illustrator
				else:
					plt.gca().axison=False
			else:
				plt.ylabel('A254', fontsize=6)
				plt.xlabel('Distance (mm)', fontsize= 6)
				ax.yaxis.set_ticks(np.arange(start, end, 0.1))
				ax.xaxis.set_ticks_position('bottom')
				ax.yaxis.set_ticks_position('left')
				ax.xaxis.set_tick_params(labelsize=6)
				ax.yaxis.set_tick_params(labelsize=6)
						
		else:
			plt.gca().axison = False #this is still drawing the left axis. Why???
			
		plt.ylim(ymax= args['maxA'])
		plt.xlim(xmin=args['clip_left'], xmax= args['clip_right']) #try setting xmin to 2 to get rid of that spike
		#actually though I think it would look better if I deal with it at the data collection step rather than the xlimit step
	
	if args['hpad']!=None:
		plt.subplots_adjust(hspace=args['hpad']) #I had it on 0.4 before
		
	plt.savefig('%s.png'% outname, transparent=True)
	plt.savefig('%s.pdf'% outname, transparent=True)	

def main(argList):
	
	parser= argparse.ArgumentParser(description= 'parse command line args')
	parser.add_argument('outfilename', help= 'outfileprefix')
	parser.add_argument('samplenamefile', help= 'plateTemplate showing plate setup')
	parser.add_argument('csvfiles', nargs='+', help='.csv files for each plate, must be in order')
	parser.add_argument('-maxA', type=float, default=2.0, help='max absorbance (in A254 units) for the y-axis')
	parser.add_argument('--noAxes', action= 'store_true', help='add option if you only want to display the names, no axes, axes labels or tick marks')
	  #also add options allowing manual override of default plotting positions
	parser.add_argument('-xdim', type=int, help='num of spaces to plot in the x-dimension')
	parser.add_argument('-ydim', type=int, help='num of spaces to plot in the y-dimension')
	parser.add_argument('-lw', type=float, default=0.75, help= 'linewidth for plotting polysomes')
	parser.add_argument('--sparseAx', action='store_true', help='turn off all axes except for leftmost ones')
	parser.add_argument('-hpad', type=float, help='the amount of vertical spacing to add between subplots, if needed')
	parser.add_argument('-mark', action='store_true', help='add this if you want to mark where the fractions were collected')
	parser.add_argument('-include_dead', action='store_true', help='add this to plot the points collected during the dead volume and readjust other fractions accordingly')
	parser.add_argument('-clip_left', type=float, default=0, help='the leftmost coordinate to plot in the graph, 2 mm is good usually')
	parser.add_argument('-clip_right', type=float, default=72, help='the rightmost coordinate to plot, 70 mm is good')
	parser.add_argument('--zero', action='store_true', help='set zero on the y-axis equal to the minimum value from the data, theyll be more square this way')
	#Note that it marks the output file where the UV monitor is when the fraction switches and not the estimated spot based on
	#the length of the tubing (i.e. 329 mm)
    
	ar=parser.parse_args(args=argList)
	args= vars(ar)
	
	outname=args['outfilename']
	labelfile=args['samplenamefile']
	csvFiles=args['csvfiles']
	
	#for now, keep csvFiles as a list you input, later, give it a directory and then let it plot all the data in a specific directory
	#SampleNames will be a separate, line-delimited file that contains the names of the samples, as will appear on the graph in same order as files are input.
	labels= readLabels(labelfile)
	
	allData=[] #store data from csvFiles here in order as dictionaries, i.e. allData=[{sampelID:ID, Abs:.., distance:..., Mark...},{},{}
	
	for k in range(0, len(csvFiles)):
		before_sampleID= csvFiles[k].split('/')[-1].split('.')[0:-1]
		sampleID= ''.join(before_sampleID)
		pathandName= csvFiles[k].split('.')[0]
		distance, Abs, Mark= read(csvFiles[k])
		d={}
		d['ID']= sampleID
		d['path']= pathandName
		d['distance']= distance
		d['Abs']= Abs
		d['Mark']= Mark
		d['Name']= labels[k]
		allData.append(d)
		
	plot(allData, outname, args)

main(sys.argv[1:])

