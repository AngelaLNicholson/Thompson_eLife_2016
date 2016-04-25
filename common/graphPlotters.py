#140317 Mary K. Thompson and Pavan Vaidyanathan
#purpose: Creates common types of plots using the analysesRunner objects

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as plticker
plt.rcParams['pdf.fonttype'] = 42
import numpy as np
import scipy.stats as stats
import argparse
import operator

from scipy.stats import gaussian_kde
from numpy.random import normal
from numpy import arange
from matplotlib.patches import Polygon
import sys
import cPickle as pickle
from analysesRunner import *
from pylab import *
import figSize
import pandas as pd
from scipy.stats import scoreatpercentile

black = (0,0,0)
orange = (230/255.0,159/255.0,0)
skyBlue = (86/255.0,180/255.0,233/255.0)
bluishGreen = (0,158/255.0,115/255.0)
yellow = (240/255.0,228/255.0,66/255.0)
blue = (0,114/255.0,178/255.0)
vermillion = (213/255.0,94/255.0,0)
reddishPurple = (204/255.0,121/255.0,167/255.0)
colors = [black, blue, reddishPurple, vermillion, orange, skyBlue, yellow, blue, yellow]
box_colors=[blue, reddishPurple, skyBlue, orange]
colorConvert={'black':black, 'orange':orange, 'skyBlue':skyBlue, 'bluishGreen':bluishGreen, 'yellow':yellow, 'blue':blue, 'vermillion':vermillion, 'reddishPurple':reddishPurple}

#default args, unless overriden by the command line
args={'numcats': None, 'barkey': None, 'axisfont': None, 'flipAxes': False, 'pickledObjects': None, 'colors': None, 'xmin': None, 'figsizey': None, 'ymin': None, 'draw_diagonal': False, 'legend_loc': 'topleft', 'legendhandles': None, 'figsizex': None, 'ymax': None, 'title': '', 'percent': False, 'legendfont': None, 'binsize': 0.1, 'draw_zeros': False, 'ylabel': '', 'test': False, 'figsize': 'full', 'ax_lw': 0.5, 'linehistfill': False, 'beta': 4, 'barplotfile': None, 'fig_scalerx': 1, 'fig_scalery': 1, 'legends': None, 'tickfont': None, 'drawbar': False, 'draw_slope': False, 'smooth': False, 'scattersize': 0.2, 'otype': None, 'xlabel': '', 'GOcatfile': None, 'xmax': None, 'panels': 1, 'barformat': None}
LegendLocs={'topright':1, 'topleft':2, 'bottomleft':3, 'bottomright':4}

class GraphPlotter(object):
	''' Base class for all graph plotting classes '''

	def __init__(self):
		super(GraphPlotter, self).__init__()
		
		pagetype=args['figsize']
		npans= args['panels']
		thisFigSize= figSize.convertFigSize(npans, pagetype)
		
		#override these if figsizex and figsizey are specified directly
		if args['figsizex']!= None:
			x= args['figsizex']/25.4 #convert mm->in
			y= args['figsizey']/25.4
			thisFigSize=(x,y)
			
		x=thisFigSize[0]*args['fig_scalerx']
		y=thisFigSize[1]*args['fig_scalery']
		
		thisFigSize=(x,y)
		
		self.title = ""
		self.xlabel = ""
		self.ylabel = ""
		self.fig = plt.figure(figsize=thisFigSize) #w, h (x,y)
		self.axes={}
		
	def set_basic_properties(self, t, xl, yl):
		''' sets the basic properties that all graphs should contain '''
		self.title = t
		self.xlabel = xl
		self.ylabel = yl
	
	def draw_subplot(self, obj, pos, fsize=7):
		'''
		implement the specific drawing function for each child class in this method.
		all subclasses should extend this method. The _multi version will draw a subplot in a given area of the canvas.
		i= the number of the plot that you are on. num_plots= the total number of plots that there will be. pos=[x,y,n] where
		x and y are the number of x and y grid spots that you want and n is the iteration that we're on. This is equivalent to draw_plot if pos=(1,1,1)
		If you are drawing a plot with the figure coordinates instead of subplotting (e.g. with the HistScatter), do not call this
		'''
		
		if pos!= None: #set pos= None if we don't want to use the subplot system but rather figure coordinates
			x,y,num= pos[0], pos[1], pos[2]
			
			self.currentax=num
			
			self.axes[num]= self.fig.add_subplot(x, y, num)
	
		else:
			self.currentax=1 #if just making one plot then set that ax to num=1
	
	def label_axes(self, thisAxis, fsize=7):
		'''
		this is an alternative to draw_subplot() that will label the specified axis
		'''
		thisAxis.set_title(self.title, fontsize=fsize)
		thisAxis.set_xlabel(self.xlabel, labelpad=2, fontsize=fsize)
		thisAxis.set_ylabel(self.ylabel, labelpad=2, fontsize=fsize)
		
	def set_limits(self, thisAxis, xmin='', xmax='', ymin='', ymax=''):
		'''Set the limits on a plot. We will assume that if you provide a min you must also provide a max and vice versa'''
		
		
		if xmin!='':
			thisAxis.set_xlim(xmin=xmin)
		if xmax!='':
			thisAxis.set_xlim(xmax=xmax)
			
		if ymin!='':
			thisAxis.set_ylim(ymin=ymin)
		
		if ymax!='':
			thisAxis.set_ylim(ymax=ymax)
	
	def set_tick_params(self, thisAxis, axcode='LB', lsize=7, multi=False):
		'''modify the properties of the axes
			use a code to specifiy which axes should be present. i.e. LB= left and bottom'''
		
		T,B,L,R= getaxsides(axcode)
		
		#use these parameters for the smaller tick marks on the multiScatter objects
		if multi==True:
			thisAxis.tick_params(axis='y', length=0.8, labelsize=lsize, pad=2) #format tick length and distance to labels (for some reason it only formats the x axis if you leave it open)
			thisAxis.tick_params(axis='x', length=0.8, labelsize=lsize, pad=2)

		else:
			thisAxis.tick_params(axis='y', length=2, labelsize=lsize, pad=2) #format tick length and distance to labels (for some reason it only formats the x axis if you leave it open)
			thisAxis.tick_params(axis='x', length=2, labelsize=lsize, pad=2)

		thisAxis.tick_params(top=T, bottom=B, right=R, left=L) #remove unwanted ticks
		

	def set_big_ax_params(self, labels, axcode='LB', lsize=7, fsize=7):
		'''set the parameters for the outside axes, for example the axes that go around all the subplots.
			fsize is the fontsize for the labels
			'''
		num_datasets= len(labels)
	
		bigAxes= self.fig.add_subplot(1,1,1, frameon=False)
		bigAxes.set_xticks(np.arange(0.0, 1.0, 1.0/num_datasets)+0.5/num_datasets)
		bigAxes.set_yticks(np.arange(0.0, 1.0, 1.0/num_datasets)+0.5/num_datasets)
		
		bigAxes.set_xticklabels([dataset for dataset in labels], fontsize=fsize)
		bigAxes.set_yticklabels([dataset for dataset in labels][::-1], rotation=90, fontsize=fsize)
		
		T,B,L,R= getaxsides(axcode)
		
		bigAxes.tick_params(top=T, bottom=B, right=R, left=L, length=1, labelsize=fsize, pad=20)
		
	def fix_spacing(self):
		self.fig.tight_layout()
		
	def save_plots(self, out_prefix, pdf=True):
		''' Saves the plot to files '''
		if pdf:
			plt.savefig(out_prefix + '.pdf', dpi=600, facecolor='w', edgecolor='w', format='pdf', transparent=True)
		else:
			plt.savefig(out_prefix + '.png', dpi=600, facecolor='w', edgecolor='w', format='png', transparent=True)
			plt.savefig(out_prefix + '.svg', dpi=600, facecolor='w', edgecolor='w', format='svg', transparent=True)
		
		plt.clf()
		#if you want to prevent clipping, add to savefig(bbox_inches='tight', however, this changes the size of the figure)

	def draw_extra_text(self, thisAxis, text_list, xy_list, transform='axis', fsize=7, ha='center'):
		''' Adds extra text information to the plot at the given x,y coordinates Note: xy_list is a list of x,y tuples '''
		if len(text_list) <> len(xy_list):
			print "GraphPlotter: draw_extra_text: the text and xy coordinate lists are not of the same length"
		else:
			if transform=='axis':
				transform=thisAxis.transAxes
			else:
				transform=thisAxis.transData
				
			for i in range(len(text_list)):
				thisAxis.text(xy_list[i][0], xy_list[i][1], text_list[i], transform=transform, fontsize=fsize, horizontalalignment=ha)
	
	def draw_legends(self, thisAxis, legenddata, legends, location, scatterpoints=1, markerscale=1, fsize=5):
		'''
		Draw and format legends
		'''
		if args['legendhandles']!=None:
			matplotlib.rcParams['legend.handlelength']=args['legendhandles'] #this seems like the appropriate scaling for the dash style I've set for the TE histograms, 1.2
		
		if legends==None:
			return
		leg=thisAxis.legend(legenddata, legends, loc=location, scatterpoints=scatterpoints, markerscale=markerscale, fontsize=fsize)

		leg.get_frame().set_linewidth(0.25)
		leg.draw_frame(False)

		
	def clear_extra_axes(self, thisAxis, axcode='LB'):
		'''
		Get rid of the top and right axes and tick marks
		'''
		
		T,B,L,R= getaxsides(axcode)
		
		if T=='off':
			thisAxis.spines['top'].set_visible(False)
		if B=='off':
			thisAxis.spines['bottom'].set_visible(False)
		if L=='off':
			thisAxis.spines['left'].set_visible(False)
		if R=='off':
			thisAxis.spines['right'].set_visible(False)
					

def getBoxParams(num_datasets, num_legends, w, find_width=True):
	'''
	calculate where to draw the boxes. Assume that we want a full width between groups and a half width between boxes in the same group.
	'''
	
	if find_width==True:
		w=float(2)/(3*num_datasets+1)
	
	groupspace= w*num_datasets+ (num_datasets-1)*w/float(2)
	boxspace= w+w/float(2)
	
	lpoints=[]
	for i in range(1, num_legends+1):
		lpoints.append(i-groupspace/2)
		
	positions=[]
	for i in lpoints:
		for j in range(0, num_datasets):
			mid_box=i+boxspace*j+w/float(2)
			positions.append(mid_box)
	
	return positions, w

def getViolinParams(num_datasets, num_legends, w, find_width=True):
	'''
	calculate where to draw the boxes. Assume that we want a full width between groups and a half width between boxes in the same group.
	Change this to be two widths between groups since the violins are fat!
	'''
	
	if find_width==True:
		w=float(2)/(3*num_datasets+3)
	
	groupspace= w*num_datasets+ (num_datasets-1)*w/float(2)
	boxspace= w+w/float(2)
	
	lpoints=[]
	for i in range(1, num_legends+1):
		lpoints.append(i-groupspace/2) #this is the leftmost point of the group
		
	positions=[]
	for i in lpoints:
		for j in range(0, num_datasets):
			mid_box=i+boxspace*j+w/float(2) #this is the midpoint for each box
			positions.append(mid_box)
	
	return positions, w

class BoxWhiskerPlotter(GraphPlotter):
#http://matplotlib.org/examples/pylab_examples/boxplot_demo2.html ##useful example
	def __init__(self):
		super(BoxWhiskerPlotter, self).__init__()
	
	def draw_subplot(self, histObjects, pos, **kwargs):
		''' Draws a box and whisker plot of the data (error_bars parameter is ignored).
		->box_width is calculated from the number of datasets provided
		->Note here that there are two types of legends. The ones plotted to the right of the graph
		will be the datasets (i.e. WT, rack1, etc.). Set these as self.legends and will be formatted with the formatIt() function
		The other labels are categories (i.e. mito, ribo, etc.) and I'll call those xcats. Will be dealt with in this function'''
		
		super(BoxWhiskerPlotter, self).draw_subplot(histObjects, pos, **kwargs)
		
		#check if we're getting a list of objects or just one object:
		if type(histObjects)!=list:
			histObjects=[histObjects] #nest this so that it can be treated the same
			
		num_datasets= len(histObjects[0].allHists)
		num_spaces=num_datasets+1
		
		currentax= self.currentax
		ax1= self.axes[currentax]
		
		dataseries=[]
		toPlot=[] #keep track of list of lists here to plot
		for i in range(0, len(histObjects[0].legends)): #for each category, i.e. bg, mito, ribo
			for j in range(0, len(histObjects)): 
				toPlot.append(histObjects[j].allVals[i]) #add that mutant's values corresponding to this category
		
		box_width=0.15
		positions, box_width= getBoxParams(len(histObjects), len(histObjects[0].legends), box_width)
	
		bp = ax1.boxplot(toPlot, positions=positions, widths=box_width, notch=0, sym='')
		
		plt.setp(bp['boxes'], color='black') #ax doesn't have .setp
		plt.setp(bp['whiskers'], color='black', linestyle='-')
		#plt.setp(bp['fliers'], color='black', marker='o')
		plt.setp(bp['fliers'], color='black', linestyle='-')
		plt.setp(bp['medians'], color='black')
		
		ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
		
		#DEAL WITH THE LEGENDS ISSUE. NOTE THAT HERE WE WANT TO MAKE THE LEGENDS THE MUTANTS AND USE THE LEGENDS TO LABEL THE XCATS
		###Get the gene counts for each group and make it to the newlegends
		xcats=histObjects[0].legends
		genecounts=[len(p) for p in histObjects[0].allVals]		
		newxcats= modifyLegends(xcats, genecounts) #get other data to add for legend annotation
		self.legends= histObjects[0].mutants
		
		#CREATE THE BOX COORDINATES
		numBoxes = len(toPlot)
		medians = range(numBoxes)
		for i in range(numBoxes):
			box = bp['boxes'][i]
			boxX = []
			boxY = []
			for j in range(5):
				boxX.append(box.get_xdata()[j])
				boxY.append(box.get_ydata()[j])
			boxCoords = zip(boxX,boxY)
			k=i%len(histObjects)
			boxPolygon = Polygon(boxCoords, facecolor=box_colors[k])
			if i<len(histObjects):
				dataseries.append(boxPolygon)
			ax1.add_patch(boxPolygon)
			
			med = bp['medians'][i]
			medianX = []
			medianY = []
			for j in range(2):
				medianX.append(med.get_xdata()[j])
				medianY.append(med.get_ydata()[j])
				ax1.plot(medianX, medianY, 'k', zorder=10)
				medians[i] = medianY[0]
		
		#LABEL THE CATEGORIES
		ticklocations1= range(1, len(xcats)+1)
		ticks1=ax1.set_xticks(ticklocations1)
		ax1.set_xticklabels(newxcats, rotation=45, fontsize=8)
		ax1.set_axisbelow(True)

		ax1.axhline(y=0, linestyle='-', color='k', zorder=1) #tried zorder to get the line to plot below, doesn't seem to work:(

		self.dataseries= dataseries


class ScatterPlotter(GraphPlotter):
	
	def __init__(self):
		super(ScatterPlotter, self).__init__()
	
	def draw_subplot(self, scatterObject, pos, **kwargs):
		''' Draws a scatter plot of the data (legends and error_bars parameter are ignored)'''

		super(ScatterPlotter, self).draw_subplot(scatterObject, pos, **kwargs)
		
		currentax= self.currentax
		ax1= self.axes[currentax]
		
		'''
		#uncomment this if you want lines at x=0 and y=0. Could add this as an option later
		ax1.axhline(y=0, linestyle='-', color='gray')
		ax1.axvline(x=0, linestyle='-', color='gray')
		'''
		
		if not scatterObject.otype=='multiScattero':
			self.legends=['All'] #don't want to do this if you're not looking at genesets
			if scatterObject.legends!= None:
				self.legends.extend(scatterObject.legends)
			self.legends= modifyLegends(self.legends, rvalues=scatterObject.pear_rvalues)
		
		allPlotData= scatterObject.allPlotData
				
		#keep track of objects to assign legends to
		dataseries=[]
		
		#plot the first data in black
		x1= scatterObject.allPlotData[0][0]
		y1= scatterObject.allPlotData[0][1]
		
		if scatterObject.flipAxes==True:
			x1= scatterObject.allPlotData[0][1]
			y1= scatterObject.allPlotData[0][0]

		
		d= ax1.scatter(x1, y1, s=args['scattersize'], color='k') #what does it mean again to put a comma after this?
		dataseries.append(d)
		
		#get the slope of the best fit line
		if args['draw_slope']==True:
			fit= polyfit(x1, y1, 1)
			slope= fit[0]
			fit_fn= poly1d(fit)
			ax1.plot(x1, fit_fn(x1), color='k')
			self.draw_extra_text(ax1, ['slope=%1.3f' % slope], [[0.8, 0.8]])		
		
		#if you want to plot an x=y line
		if args['draw_diagonal']==True:
			view=ax1.xaxis.get_view_interval()
			now= view[0]
			maxx= view[1]
			x_points=[]
			while now< maxx:
				x_points.append(now)
				now+= 0.1
			y_points=[i for i in x_points]
			ax1.plot(x_points, y_points, linestyle='-', color='r', linewidth=0.75)
		
		for i in range(1, len(allPlotData)):
			
			x=allPlotData[i][0]
			y=allPlotData[i][1]
			
			if scatterObject.flipAxes==True:
				x= allPlotData[i][1]
				y= allPlotData[i][0]
			
			d = ax1.scatter(x, y, s=args['scattersize'], color=colors[i])
			dataseries.append(d)
		
		
		#ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
		#ax1.xaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
		if args['draw_zeros']==True:
			ax1.axhline(y=0, linestyle='--', color='k')
			ax1.axvline(x=0, linestyle='--', color='k')
		self.dataseries= dataseries

class ScatterHistPlotter(GraphPlotter):
	'''
        I'm making this for the Hist_Scatter. Doesn't need to be subplottable.
        '''
        
	def __init__(self):
		super(ScatterHistPlotter, self).__init__()
	
	def draw_subplot(self, scatterObject, pos, **kwargs):
		''' Draws a scatter plot of the data (legends and error_bars parameter are ignored)'''
		
		super(ScatterHistPlotter, self).draw_subplot(scatterObject, pos, **kwargs)

		spotsize=0.2
		self.legends=['All'] #don't want to do this if you're not looking at genesets	
		self.legends.extend(scatterObject.legends)
		self.legends= modifyLegends(self.legends, rvalues=scatterObject.pear_rvalues)
		allPlotData= scatterObject.allPlotData
		
		#keep track of objects to assign legends to
		dataseries=[]
		
		#plot the first data in black
		x1= scatterObject.allPlotData[0][0]
		y1= scatterObject.allPlotData[0][1]
		
		if scatterObject.flipAxes==True:
			x1= scatterObject.allPlotData[0][1]
			y1= scatterObject.allPlotData[0][0]

		#manually set up where the axes are:
		left, width = 0.1, 0.65
                bottom, height = 0.1, 0.65
                bottom_h = left_h = left+width+0.0
                
                rect_scatter = [left, bottom, width, height]
                rect_histx = [left, bottom_h, width, 0.2]
                rect_histy = [left_h, bottom, 0.2, height]
                
                axScatter = plt.axes(rect_scatter)
		self.axes[self.currentax]= axScatter
		
                axHistx = plt.axes(rect_histx)
                axHisty = plt.axes(rect_histy)
		
		self.axes[2]= axHistx
		self.axes[3]= axHisty #assign x to 2 and y to 3 so that can manipulate these later
		
		#this clears away all the extra ticks and axes for the histogram portion
		self.clear_extra_axes(axHistx, axcode='B')
		self.set_tick_params(axHistx, axcode='', lsize=0)
		self.clear_extra_axes(axHisty, axcode='L')
		self.set_tick_params(axHisty, axcode='', lsize=0)
                
		#plot the background distribution
		d= axScatter.scatter(x1, y1, s=spotsize, color='k') #changed s to smaller number (6 to 4) for plotting on less of the page
		dataseries.append(d)
		
                #now determine nice limits by hand:
		binwidth=0.1
                xymax = np.max( [np.max(np.fabs(x1)), np.max(np.fabs(y1))] )
                lim = ( int(xymax/binwidth) + 1) * binwidth
                
                axScatter.set_xlim( (-lim, lim) )
                axScatter.set_ylim( (-lim, lim) )
                
                bins = np.arange(-lim, lim + binwidth, binwidth)
                axHistx.hist(x1, bins=bins, histtype="stepfilled", alpha=0.2, color='k', normed=True)
                axHistx.hist(x1, bins=bins, histtype="step", color='k', alpha=0.5, normed=True)

                axHisty.hist(y1, bins=bins, orientation='horizontal', histtype="stepfilled", alpha=0.2, color='k', normed=True)
                axHisty.hist(y1, bins=bins, orientation='horizontal', histtype="step", color='k', alpha=0.5, normed=True)
		
		#This plots a y=x line
		view=axScatter.xaxis.get_view_interval()
		now= view[0]
		maxx= view[1]
		maxx=3 ##manually assign
		x_points=[]
		while now< maxx:
			x_points.append(now)
			now+= 0.1
		y_points=[i for i in x_points]
		axScatter.plot(x_points, y_points, linestyle='-', color='r', linewidth=0.75)
		
		#get the slope of the best fit line
		if args['draw_slope']==True:
			fit= polyfit(x1, y1, 1)
			slope= fit[0]
			fit_fn= poly1d(fit)
			axScatter.plot(x1, fit_fn(x1), color='k')
			self.draw_extra_text(axScatter, ['slope=%1.3f' % slope], [[0.8, 0.8]])
			
		for i in range(1, len(allPlotData)):
			
			x=allPlotData[i][0]
			y=allPlotData[i][1]
			
			if scatterObject.flipAxes==True:
				x= allPlotData[i][1]
				y= allPlotData[i][0]
			
			d = axScatter.scatter(x, y, s=spotsize, color=colors[i])
			dataseries.append(d)
                        
                        axHistx.hist(x, bins=bins, histtype="stepfilled", alpha=0.5, color=colors[i], normed=True)
                        axHistx.hist(x, bins=bins, histtype="step", color=colors[i], normed=True)

                        axHisty.hist(y, bins=bins, orientation='horizontal', histtype="stepfilled", alpha=0.5, color=colors[i], normed=True)
                        axHisty.hist(y, bins=bins, orientation='horizontal', histtype="step", color=colors[i], normed=True)
                
		#ax1.yaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
		#ax1.xaxis.grid(True, linestyle='-', which='major', color='lightgrey', alpha=0.5)
		axScatter.axhline(y=0, linestyle='--', color='k', linewidth=0.75)
		axScatter.axvline(x=0, linestyle='--', color='k', linewidth=0.75)
		self.dataseries=dataseries
			    
class CDFPlotter(GraphPlotter):
	'''
	Input supported:
	One hist object, not a list of histObjects-> would need to link this to grid subplotting. Maybe in the future.
	Add the number of genes and the pvalues on the graph now
	'''
	def __init__(self):
		super(CDFPlotter, self).__init__()
	
	def draw_subplot(self, histObject, pos, **kwargs):
		super(CDFPlotter, self).draw_subplot(histObject, pos, **kwargs)
		
		currentax= self.currentax
		ax1=self.axes[currentax]
		
		arrays= histObject.allCdfs
		edges= histObject.allEdges
		legends= histObject.legends
		genecounts= histObject.genecounts
		pvals= histObject.pvals
		
		#modify the legends with the number of genes
		newlegends= modifyLegends(legends, genecounts, pvals)
		newlegends=[newlegends[k] for k in range(len(newlegends)) if k%2==0]
					
		currentax= self.currentax
		ax1=self.axes[currentax]
		
		dataseries=[]
		
		bg, =ax1.plot(edges[0][1:], arrays[0], label=legends[0], linewidth=.75, color=colors[0], linestyle=colorstyles[0])
		dataseries.append(bg)
		
		for i in range(1, len(arrays)):
			s,=ax1.plot(edges[i][1:], arrays[i], label=legends[i], linewidth=.75, color=colors[i], linestyle=colorstyles[i])
			dataseries.append(s)
			
		self.legends= newlegends
		self.dataseries=dataseries		

class LineHistPlotter(GraphPlotter):
	'''
	Make a histogram by drawing a line. Can choose bin width and smoothing parameters
	'''
	def __init__(self):
		super(LineHistPlotter, self).__init__()
		
	def draw_subplot(self, histObject, pos, **kwargs):
		super(LineHistPlotter, self).draw_subplot(histObject, pos, **kwargs)
		
		#http://glowingpython.blogspot.com/2012/02/convolution-with-numpy.html
		#comments: smoothing seems to be dependent both on the bin size you choose and beta. Higher beta = more smooth
		#in this way smoothing definitely disguises the shape of the data a little bit depending on what you choose--> maybe a little artificial, be careful
		currentax= self.currentax
		ax1= self.axes[currentax]
		
		#check if we're getting a list of objects or just one object:
		if type(histObject)!=list:
			histObject=[histObject] #nest this so that it can be treated the same
		
		allVals= histObject[0].allVals
		legends= histObject[0].legends
		if args['legends']!=None:
			legends= args['legends']
		
		dataseries=[]
	
		#now determine nice limits by hand:
		binwidth=args['binsize'] #the smoothing effect is fairly
		beta= args['beta']
                xymax = np.max( [np.max(np.fabs(allVals[0])), np.max(np.fabs(allVals[0]))] )
                lim = ( int(xymax/binwidth) + 1) * binwidth
		bins = np.arange(-lim, lim + binwidth, binwidth)
		
		#TRY IMPLEMENTING A DIFFERENT NORMALIZATION SCHEME TO MAKE IT LOOK LIKE A FREQUENCY HISTOGRAM
		for i in range(0, len(histObject)):
			ny, bin_edges = np.histogram(np.array(histObject[i].allVals), bins=bins, density=False)
			#use the middle of each xbin for plotting
			xs= [np.mean(bin_edges[m:m+2]) for m in range(0, len(bin_edges)-1)]
			
			totaly= float(sum(ny))
			ny_frac=[]
			for p in ny:
			    frac= p/totaly
			    ny_frac.append(frac)
		   #uncomment if you want the first line to be grey instead of black
			if i==0:
				fillalphalevel=0.3
				linealphalevel=0.5
			else:
				fillalphalevel=0.5
				linealphalevel=1
			
			
			#KAISER WINDOW SMOOTHING
			if args['smooth']==True:
				yy=self.KaiserSmooth(ny_frac,beta)
			else:
				yy= ny_frac
			
			s,=ax1.plot(xs, yy, alpha=linealphalevel, color=colors[i], linestyle=colorstyles[i], linewidth=0.5)
			if colorstyles[i]=='--':
				s.set_dashes([1])
			if args['linehistfill']==True:
				ax1.fill_between(xs, yy, alpha=fillalphalevel, color=colors[i])
				testrect= plt.Rectangle((bin_edges[0],0), 0,0, facecolor=colors[i],alpha=fillalphalevel)
				dataseries.append(testrect)
			else:
				dataseries.append(s)
				
			    
		self.legends= legends
		self.dataseries= dataseries
		
	def KaiserSmooth(self, x, beta):
		""" kaiser window smoothing
		from http://glowingpython.blogspot.co.uk/2012/02/convolution-with-numpy.html"""
		window_len=11
		# extending the data at beginning and at the end
		# to apply the window at the borders
		s = np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
		w = np.kaiser(window_len,beta)
		y = np.convolve(w/w.sum(),s,mode='valid')
		return y[5:len(y)-5]

class BarHistPlotter(GraphPlotter):
	'''
	Try to make a prettier version than HistPlotter.
	For first pass, normalize after the hist and plot bars to make a BarHist looking plot.-> could also incorporate
	this into the barHist function as an option--> perhaps that would be more cohesive. The only difference right now is how the limits are determined
	Then add smoothing to see how it looks
	'''
	def __init__(self):
		super(BarHistPlotter, self).__init__()
		
	def draw_subplot(self, histObject, pos, **kwargs):
		super(BarHistPlotter, self).draw_subplot(histObject, pos, **kwargs)
		
		currentax= self.currentax
		ax1= self.axes[currentax]
		
		hists= histObject.allHists
		edges= histObject.allEdges
		legends= histObject.legends
		genecounts= histObject.genecounts
		pvalues= histObject.pvals #note that pvalues start one back since the first one is a comparison of geneset1 to bg
		allVals= histObject.allVals
		
		dataseries=[]
		#now determine nice limits by hand:
		binwidth=args['binsize']
                xymax = np.max( [np.max(np.fabs(allVals[0])), np.max(np.fabs(allVals[0]))] )
                lim = ( int(xymax/binwidth) + 1) * binwidth
		bins = np.arange(-lim, lim + binwidth, binwidth)
		
		#TRY IMPLEMENTING A DIFFERENT NORMALIZATION SCHEME TO MAKE IT LOOK LIKE A FREQUENCY HISTOGRAM
		for i in range(0, len(allVals)):
			ny, bin_edges = np.histogram(np.array(allVals[i]), bins=bins, density=False)
			totaly= float(sum(ny))
			ny_frac=[]
			for p in ny:
			    frac= p/totaly
			    ny_frac.append(frac)
			
			if i==0:
				fillalphalevel=0.2
				linealphalevel=0.5
			else:
				fillalphalevel=0.5
				linealphalevel=1
		
			if args['drawbar']==True: #draw the actual bars
				testrect=ax1.bar(bin_edges[0:-1], ny_frac, width=binwidth, facecolor=colors[i], alpha=fillalphalevel)
				ax1.bar(bin_edges[0:-1], ny_frac, width=binwidth, facecolor='none')

			else: #just outline the top of the bars
				newbins=[bin_edges[q] for q in range(1, len(bin_edges)) for j in range(2)]
				#try duplicating the one to the right first so that left side and right side have the same y-value
				doubledbins= [bin_edges[0]]
				doubledbins.extend(newbins)
				doubledbins=doubledbins[0:-1]
				newfrac=[ny_frac[q] for q in range(len(ny_frac)) for j in range(2)]
				
				ax1.plot(doubledbins, newfrac, alpha=linealphalevel, color=colors[i])
	
				ax1.fill_between(doubledbins, newfrac, alpha=fillalphalevel, color=colors[i])
				testrect= plt.Rectangle((doubledbins[0],0), 0,0, facecolor=colors[i],alpha=fillalphalevel)
	
			dataseries.append(testrect)
			    
		self.legends= legends
		self.dataseries= dataseries
	
class StackedHistPlotter(GraphPlotter):
	'''
	Plot stacked histogram, where the subsets make up part of the full histogram
	'''
	def __init__(self):
		super(StackedHistPlotter, self).__init__()
		
	def draw_subplot(self, histObject, pos, **kwargs):
		super(StackedHistPlotter, self).draw_subplot(histObject, pos, **kwargs)
		
		currentax= self.currentax
		ax1= self.axes[currentax]
		
		hists= histObject.allHists
		edges= histObject.allEdges
		legends= histObject.legends
		genecounts= histObject.genecounts
		pvalues= histObject.pvals #note that pvalues start one back since the first one is a comparison of geneset1 to bg
		allVals= histObject.allVals
		
		dataseries=[]
		#now determine nice limits by hand:
		binwidth=args['binsize']
                xymax = np.max( [np.max(np.fabs(allVals[0])), np.max(np.fabs(allVals[0]))] )
                lim = ( int(xymax/binwidth) + 1) * binwidth
		bins = np.arange(-lim, lim + binwidth, binwidth)
		
		totaly=0
		#TRY IMPLEMENTING A DIFFERENT NORMALIZATION SCHEME TO MAKE IT LOOK LIKE A FREQUENCY HISTOGRAM
		for i in range(0, len(allVals)):
			ny, bin_edges = np.histogram(np.array(allVals[i]), bins=bins, density=False)
			if totaly==0:
				totaly= float(sum(ny)) #by doing this we ensure that it's only calculated 1x/round
			ny_frac=[]
			for p in ny:
			    frac= p/totaly
			    ny_frac.append(frac)
			
			if i==0:
				fillalphalevel=0.2
				linealphalevel=0.5
			else:
				fillalphalevel=0.5
				linealphalevel=1
		
			testrect=ax1.bar(bin_edges[0:-1], ny_frac, width=binwidth, facecolor=colors[i], alpha=fillalphalevel)
			ax1.bar(bin_edges[0:-1], ny_frac, width=binwidth, facecolor='none')
			dataseries.append(testrect)
			    
		self.legends= legends
		self.dataseries= dataseries
	
class BarLineHistPlotter(GraphPlotter):
	'''
	This version will plot the first series as a bar or stephist and the second series as a line on top of that.
	'''
	def __init__(self):
		super(BarLineHistPlotter, self).__init__()
		
	def draw_subplot(self, histObjects, pos, **kwargs):
		super(BarLineHistPlotter, self).draw_subplot(histObjects, pos, **kwargs)
		
		barlinecolors=[black, blue]
		drawline=False
		currentax= self.currentax
		ax1= self.axes[currentax]
		
		#check if we're getting a list of objects or just one object:
		if type(histObjects)!=list:
			histObjects=[histObjects] #nest this so that it can be treated the same
		
		allVals= histObjects[0].allVals
		legends= histObjects[0].legends
		
		dataseries=[]
		#now determine nice limits by hand:
		binwidth=args['binsize']
                xymax = np.max( [np.max(np.fabs(allVals[0])), np.max(np.fabs(allVals[0]))] )
                lim = ( int(xymax/binwidth) + 1) * binwidth
		bins = np.arange(-lim, lim + binwidth, binwidth)
		
		#TRY IMPLEMENTING A DIFFERENT NORMALIZATION SCHEME TO MAKE IT LOOK LIKE A FREQUENCY HISTOGRAM
		for i in range(0, 1):
			allVals=histObjects[i].allVals
			ny, bin_edges = np.histogram(np.array(allVals[0]), bins=bins, density=False)
			totaly= float(sum(ny))
			ny_frac=[]
			for p in ny:
			    frac= p/totaly
			    ny_frac.append(frac)
			
			if i==0:
				fillalphalevel=0.2
				linealphalevel=0.5
			else:
				fillalphalevel=0.5
				linealphalevel=1
		
			if args['drawbar']==True: #draw the actual bars
				testrect=ax1.bar(bin_edges[0:-1], ny_frac, width=binwidth, facecolor=barlinecolors[i], alpha=fillalphalevel)
				ax1.bar(bin_edges[0:-1], ny_frac, width=binwidth, facecolor='none', alpha=linealphalevel, edgecolor= barlinecolors[i])

			else: #just outline the top of the bars
				newbins=[bin_edges[q] for q in range(1, len(bin_edges)) for j in range(2)]
				#try duplicating the one to the right first so that left side and right side have the same y-value
				doubledbins= [bin_edges[0]]
				doubledbins.extend(newbins)
				doubledbins=doubledbins[0:-1]
				newfrac=[ny_frac[q] for q in range(len(ny_frac)) for j in range(2)]
				
				ax1.plot(doubledbins, newfrac, alpha=linealphalevel, color=barlinecolors[i])
	
				ax1.fill_between(doubledbins, newfrac, alpha=fillalphalevel, color=barlinecolors[i])
				testrect= plt.Rectangle((doubledbins[0],0), 0,0, facecolor=barlinecolors[i],alpha=fillalphalevel)
			
		for i in range(1, 2):
			if drawline==True: #this will just produce a line without the bar shape
				allVals= histObjects[i].allVals
				ny, bin_edges = np.histogram(np.array(allVals[0]), bins=bins, density=False)
				#use the middle of each xbin for plotting
				xs= [np.mean(bin_edges[m:m+2]) for m in range(0, len(bin_edges)-1)]
				
				totaly= float(sum(ny))
				ny_frac=[]
				for p in ny:
				    frac= p/totaly
				    ny_frac.append(frac)
				    
				p,=ax1.plot(xs, ny_frac, alpha=linealphalevel, color=barlinecolors[i])
				#ax1.fill_between(xs, yy, alpha=fillalphalevel, color=colors[i])
				dataseries.append(p)
			else: #this line will be the outline of the bars
				allVals=histObjects[i].allVals
				ny, bin_edges = np.histogram(np.array(allVals[0]), bins=bins, density=False)
				totaly= float(sum(ny))
				ny_frac=[]
				for p in ny:
				    frac= p/totaly
				    ny_frac.append(frac)
				
				
				newbins=[bin_edges[q] for q in range(1, len(bin_edges)) for j in range(2)]
				#try duplicating the one to the right first so that left side and right side have the same y-value
				doubledbins= [bin_edges[0]]
				doubledbins.extend(newbins)
				doubledbins=doubledbins[0:-1]
				newfrac=[ny_frac[q] for q in range(len(ny_frac)) for j in range(2)]
				
				p, =ax1.plot(doubledbins, newfrac, alpha=1, color=barlinecolors[i])

				dataseries.append(p)    
		
			    
		self.legends= ['replicate error']
		self.dataseries= dataseries
		
class BarPlotter(GraphPlotter):
	
	def __init__(self):
		super(BarPlotter, self).__init__()
	
	def draw_plot(self, data, legends, error_bars, categories, **kwargs):
		
		#here data will be a list of series values:
		#i.e. data=[[series1 data],[series2 data]]
		#legends will be a list of names corresponding to each series:
		#i.e. legends=['series1name', 'series2name', etc.]
		#error bars= same dimensions as the data, error to plot for each point
	
		ax= self.ax1
		graphcolors=[blue, vermillion, orange, reddishPurple]

		n= len(data[0]) #this is the number of plots that we're putting on the same axis
		m= len(data)
		width=.5
		bar_width= width/m
		
		ind= np.arange(n)
				
		#figure out how wide to make the other bars based on the number of datasets:
		num_bars= n
		
		half_width= width/2
	       
		mid_points=[]
		for i in ind:
			mid_points.append(ind[i]+0.5)
		    
		left_points=[] #this is where to put the leftmost point of the first bars for the bar graph
		for i in mid_points:
			left_points.append(i-half_width)
		    
		series_lpts={} #the leftmost bar line for each data series
		for i in range(0, n): #for each series
			pts=[]
			for j in range(0, len(left_points)):
				p=left_points[j]+bar_width*i
				pts.append(p)
			series_lpts[i]= pts
		
		allRects=[]
		for s in range(0, len(data)):
			print 'data', data[s]
			print 'series_lpts', series_lpts[s]
			print 'err', error_bars[s]
			
			rects=ax.bar(series_lpts[s], data[s], bar_width, color=graphcolors[s], yerr=error_bars[s], ecolor='k')
			allRects.append(rects)
			    
		ax.set_xticks(mid_points) # removing this line doesn't make the ticks go away and it puts them at default positions (i.e. starting at zero)
		ax.set_xticklabels(categories, rotation=45, fontsize=7) ##maybe nice to remove x-ticks later? They are not doing anything???
		
		#a hacky way to get rid of the ugly ticks:
		#which color is the in middle?
		if n%2==0: #even number, so tick will overlap with bar boundary anyway
			tickcolor= black
		else:
			tickcolor= graphcolors[int(math.floor(n/float(2)))]
		    
		axis= ax.xaxis
		ticklines= axis.get_ticklines()
		for line in ticklines:
			matplotlib.artist.setp(line, color= tickcolor)
		    
		
		#add legends for the series
		legenddata=[i[0] for i in allRects]
		ax.legend(legenddata, legends, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0., fontsize=7)

def convertToPercent(data, logbase):
	'''Given a list of lists, convert to percentage points based on the log base'''
	transformeddata=[]
	for d in data:
		newvs=[]
		for i in d:
			v= (logbase**i)*100-100
			newvs.append(v)
			
		transformeddata.append(newvs)
	return transformeddata
	
class ViolinPlotter(GraphPlotter):
	''' Adapted from Flavio Coelho: http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html'''
	def __init__(self):
		super(ViolinPlotter, self).__init__()
	
	def draw_subplot(self, histObjects, pos, **kwargs):
		
		super(ViolinPlotter, self).draw_subplot(histObjects, pos, **kwargs)
		plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.25)
		
		#check if we're getting a list of objects or just one object:
		if type(histObjects)!=list:
			histObjects=[histObjects] #nest this so that it can be treated the same
		
		ax1=self.axes[self.currentax]
		
		dataseries=[]
		data=[] #keep track of list of lists here to plot
		for i in range(0, len(histObjects[0].legends)): #for each category, i.e. bg, mito, ribo
			for j in range(0, len(histObjects)): 
				data.append(histObjects[j].allVals[i]) #add that mutant's values corresponding to this category
		
		if args['percent']==True: #convert the y-axis scale to percentage points
			data= convertToPercent(data, 2)
				
		box_width=0.15
		pos, w= getViolinParams(len(histObjects), len(histObjects[0].legends), box_width)

		dist = max(pos)-min(pos)
		w=w*.6
		#w = min(0.15*max(dist,1.0),0.5)
		#w = min(0.15*max(dist,1.0),0.5/(len(histObjects)+1)) #determined this empirically
		#w = min(0.04*max(dist,1.0),0.5/(len(histObjects)+1)) #determined this empirically

		zippeddata= zip(data, pos)
		for i in range(len(data)):
			d,p= zippeddata[i]
			d = [float(dm) for dm in d]
			k = gaussian_kde(d) #calculates the kernel density
			
			#the tukey method for gettting the whiskers
			q1= scoreatpercentile(data[i], 25)
			q3= scoreatpercentile(data[i], 75)
			iqr= q3-q1
			
			whiskrange= 1.5*iqr
			mw= q1- whiskrange
			Mw= q3+ whiskrange
			data[i].sort()
			
			#find lower whisker
			tmw=''
			tMw=''
			for j in range(0, len(data[i])):
				if data[i][j]>= mw:
					tmw= data[i][j]
					break
			for l in range(j+1, len(data[i])):
				if data[i][l]> Mw:
					tMw= data[i][l-1]
					break
			if tmw=='': #didn't get reassigned, use min
				tmw= data[i][0]
			if tMw=='':
				tMw= data[i][-1]
				
			m= tmw
			M= tMw
			
			x = arange(m,M,(M-m)/100.) # support for violin
			v = k.evaluate(x) #violin profile (density curve)
			median= np.median(d)
			colori= i%len(histObjects)
			v = v/v.max()*w #scaling the violin to the available space
			###below fills both sides separately
			plt.fill_betweenx(x,p,v+p,facecolor=box_colors[colori],alpha=1, linewidth=0.75) ##this is what makes the violin plots
			plt.fill_betweenx(x,p,-v+p,facecolor=box_colors[colori],alpha=1, linewidth=0.75)
			
			if i<len(histObjects):
				testrect= plt.Rectangle((p,0), 0,0, facecolor=box_colors[colori],alpha=1, linewidth=0.75)
				dataseries.append(testrect)
			ax1.plot([p-w, p+w],[median, median], 'k', linewidth=0.75) #x chould be beginning coordinate, end coordinate, y should be the median
	
		
		#LABEL THE CATEGORIES
		xcats= histObjects[0].legends
		genecounts=[len(p) for p in histObjects[0].allVals]		
		newxcats= modifyLegends(xcats, genecounts) #get other data to add for legend annotation
		
		ticklocations1= range(1, len(xcats)+1)
		ticks1=ax1.set_xticks(ticklocations1)
		ax1.set_xticklabels(newxcats, rotation=45, fontsize=8)
		ax1.set_axisbelow(True)
		
		ax1.axhline(y=0, linestyle='-', color='k', zorder=1, linewidth=0.75) #tried zorder to get the line to plot below, doesn't seem to work:(
		
		self.dataseries=dataseries
		self.legends= histObjects[0].mutants

class ViolinPlotterSig(GraphPlotter):
	''' Adapted from Flavio Coelho: http://pyinsci.blogspot.com/2009/09/violin-plot-with-matplotlib.html
	This version I'm adapting to perform the Mann Whitney U test on the data and indicate the pvalues on the plot'''
	
	def __init__(self):
		super(ViolinPlotterSig, self).__init__()
	
	def draw_subplot(self, histObjects, pos, **kwargs):
		
		super(ViolinPlotterSig, self).draw_subplot(histObjects, pos, **kwargs)
		plt.subplots_adjust(left=0.1, right=0.95, top=0.9, bottom=0.25)
		
		#check if we're getting a list of objects or just one object:
		if type(histObjects)!=list:
			histObjects=[histObjects] #nest this so that it can be treated the same
		
		ax1=self.axes[self.currentax]
		
		dataseries=[]
		data=[] #keep track of list of lists here to plot
		for i in range(0, len(histObjects[0].legends)): #for each category, i.e. bg, mito, ribo
			for j in range(0, len(histObjects)): 
				data.append(histObjects[j].allVals[i]) #add that mutant's values corresponding to this category
		
		if args['percent']==True: #convert the y-axis scale to percentage points
			data= convertToPercent(data, 2)
				
		box_width=0.15
		pos, w= getViolinParams(len(histObjects), len(histObjects[0].legends), box_width)

		dist = max(pos)-min(pos)
		w=w*.6
		#w = min(0.15*max(dist,1.0),0.5)
		#w = min(0.15*max(dist,1.0),0.5/(len(histObjects)+1)) #determined this empirically
		#w = min(0.04*max(dist,1.0),0.5/(len(histObjects)+1)) #determined this empirically

		graphMaxes=[]
		
		zippeddata= zip(data, pos)
		ps=[] #keep track of the x position of the plot to use later for plotting the pvalues
		for i in range(len(data)):
			d,p= zippeddata[i]
			ps.append(p)
			d = [float(dm) for dm in d]
			k = gaussian_kde(d) #calculates the kernel density
	
			#the tukey method for gettting the whiskers
			q1= scoreatpercentile(data[i], 25)
			q3= scoreatpercentile(data[i], 75)
			iqr= q3-q1
			
			whiskrange= 1.5*iqr
			mw= q1- whiskrange
			Mw= q3+ whiskrange 
			data[i].sort()
			
			#find the whiskers within the range
			tmw=''
			tMw=''
			for j in range(0, len(data[i])):
				if data[i][j]>= mw:
					tmw= data[i][j]
					break
			for l in range(j+1, len(data[i])):
				if data[i][l]> Mw:
					tMw= data[i][l-1]
					break
				
			if tmw=='': #didn't get reassigned, use min
				tmw= data[i][0]
			if tMw=='':
				tMw= data[i][-1]
				
			m= tmw
			M= tMw
			#append M to the maxes so that I can use this for plotting the asterices
			graphMaxes.append(M)
			
			x = arange(m,M,(M-m)/100.) # support for violin
			v = k.evaluate(x) #violin profile (density curve)
			
			median= np.median(d)
			colori= i%len(histObjects)
			v = v/v.max()*w #scaling the violin to the available space
			
			if args['colors']!=None:
				box_colors= [colorConvert[q] for q in args['colors']]
			
			#edgecolor= box_colors[colori] #if you want to make the line disappear
			plt.fill_betweenx(x,v+p,-v+p,facecolor=box_colors[colori],alpha=1, linewidth=0.25)

			if i<len(histObjects):
				testrect= plt.Rectangle((p,0), 0,0, facecolor=box_colors[colori],alpha=1, linewidth=0.25)
				dataseries.append(testrect)
			ax1.plot([p-w, p+w],[median, median], 'k', linewidth=0.25) #x chould be beginning coordinate, end coordinate, y should be the median
			#bp = ax1.boxplot(data, positions=pos, widths=w, notch=0, sym='') #this plots a box plot on top, looks kind of ugly.

		#LABEL THE CATEGORIES
		xcats= histObjects[0].legends
		genecounts=[len(p) for p in histObjects[0].allVals]		
		newxcats= modifyLegends(xcats, genecounts) #get other data to add for legend annotation
		
		ticklocations1= range(1, len(xcats)+1)
		ticks1=ax1.set_xticks(ticklocations1)
		newxcats=[]
		for i in range(0, len(xcats)):
			newlabel= '%s\nn=%s' % (xcats[i], genecounts[i])
			newxcats.append(newlabel)
			
		ax1.set_xticklabels(newxcats, fontsize=6)
		ax1.set_axisbelow(True)

		ax1.axhline(y=0, linestyle='-', color='k', zorder=1, linewidth=0.5, alpha=0.3) #tried zorder to get the line to plot below, doesn't seem to work:(

		self.dataseries=dataseries
		self.legends= histObjects[0].mutants
		
		#collect all pvalues into list before drawing them onto plot:
		pvalstoPlot=[]
			
		#histObjects[0]
		for i in range(0, len(histObjects[0].legends)-1): #iterate through the number of categories, start at 0 b/c there are n-1 pvalues
			print 'legend num', i
			for j in range(0, len(histObjects)): #add the jth pvalue for that category
				print 'series num', j
				pval= histObjects[j].pvals[i]
				pvalstoPlot.append(pval)
			
		print 'pvalstoPlot', pvalstoPlot #note that these are not going to be in the right order now with respect to the plot coordinates, need to add iteratively
		
		#exclude the bg plot positions, no pvalues associated with these
		skip= len(histObjects)#the number of series is = to the number of background positions that there are, skip these
		for k in range(0, len(ps)-skip): #can change this to skip the first ones (all genes) that won't have an associated pvalues
			thispos= ps[k+skip]/3 -0.5/3#-(float(1)/3)*(1/3.5) #this makes it so that center is at the left
			thisPval= pvalstoPlot[k]
				
			#add stars instead of pvalues on top of the violin plot
			yrange=args['ymax']-args['ymin']
			aspad= yrange*0.02
		
			if thisPval<1e-18:
				ax1.text(ps[k+skip], graphMaxes[k+skip]+aspad, '***', fontsize=5, ha='center')
			elif thisPval<1e-9:
				ax1.text(ps[k+skip], graphMaxes[k+skip]+aspad, '**', fontsize=5, ha='center')
			elif thisPval<1e-3:
				ax1.text(ps[k+skip], graphMaxes[k+skip]+aspad, '*', fontsize=5, ha='center')

def parseBarFile(infile, formatting='GO'):
	'''parse the input file for the barplot'''
	if formatting=='GO':
		
		d={}
		f=open(infile, 'r')
		header= f.readline().strip('\n').split('\t')
		d['header']= header[1:] #add this to keep track of the order that we want them to appear in
		d['xcats']=[]
		for p in header:
			d[p]={}
		for line in f:
			fields= line.strip('\n').split('\t')
			d['xcats'].append(fields[0])
			for i in range(1, len(fields)):
				try:
					d[header[i]][fields[0]]= -float(fields[i]) #make this a negative numbers so that it will graph the pvalues as positive
				except:
					d[header[i]][fields[0]]= ''

		f.close()
		
	else: #make this format for a series with error bars in the next column
		d={}
		f=open(infile, 'r')
		header= f.readline().strip('\n').split('\t')
		d['header']= header[1:] #add this to keep track of the order that we want them to appear in
		d['xcats']=[]
		for p in header:
			d[p]={}
		for line in f:
			fields= line.strip('\n').split('\t')
			d['xcats'].append(fields[0])
			for i in range(1, len(fields)):
				try:
					d[header[i]][fields[0]]= float(fields[i])
				except:
					d[header[i]][fields[0]]= ''

		f.close()
		
	return d

def getBarValues(d, key, numcats, GOcatsfile=None):
	'''Return the values to use for making the barplot.
	Can use to sort by value for a specific mutant and extract the corresponding values for the other mutants
	key should be the name of the key in the dictionary (i.e. the dataseries) that you want to sort'''
	
	legends= d['header']
	
	if GOcatsfile==None:	
		thisD= d[key]
		sortedthisD= sorted(thisD.iteritems(), key=operator.itemgetter(1))
		xcats=[]
		for i in range(numcats): #change to however many categories that you want to plot
			xcats.append(sortedthisD[i][0])
		
	else:
		xcats=[]
		f=open(GOcatsfile, 'r')
		for line in f:
			cat= line.strip('\n')
			xcats.append(cat)
		f.close()
		
	toPlot=[]
	for i in legends:
		thisSeries=[]
		for j in xcats:
			val= d[i][j]
			thisSeries.append(val)
		toPlot.append(thisSeries)
		
		
	
	return legends, xcats, toPlot
	
		
class BarPlotterSimple(GraphPlotter):
	'''This one plots a bar graph given an input file with the values already calculated'''
	def __init__(self):
		super(BarPlotterSimple, self).__init__()
	def draw_subplot(self, infile, pos, **kwargs):
		''' Draws a scatter plot of the data (legends and error_bars parameter are ignored)'''
		super(BarPlotterSimple, self).draw_subplot(infile, pos, **kwargs)
		
	
		graphcolors= colors[1:]
		
		#PARSE THE INPUT FILE:
		bardict= parseBarFile(infile, formatting=args['barformat'])
				
		legends=bardict['header']
		if 'std dev' in legends:
			end=1 #this is just one series, stop here
		else:
			end= len(legends)

		xcats= bardict['xcats']
		
		toPlot=[] #make list of lists, where each list are the values from one mutant
		for i in range(0, end):
			l=[]
			for j in range(0, len(xcats)):
				l.append(bardict[legends[i]][xcats[j]])
		
			toPlot.append(l)
			
		errortoPlot=[]
		if end==1:
			e=[]
			for j in range(0, len(xcats)):
				e.append(bardict[legends[end]][xcats[j]])
		
			errortoPlot.append(e)
		
		#JUST DEAL WITH SOME OF THE FORMATTING BELOW:
		ax= self.axes[self.currentax]
		
		n= len(legends) #this is the number of bars that we're putting in each group
		
		width=.5
		bar_width= width/n
		
		ind= np.arange(len(xcats))
				
		#figure out how wide to make the other bars based on the number of datasets:
		num_bars= len(toPlot)
		
		half_width= width/2
	       
		mid_points=[]
		for i in ind:
			mid_points.append(ind[i]+0.5)
		    
		left_points=[] #this is where to put the bars
		for i in mid_points:
			left_points.append(i-half_width)
		    
		series_lpts={} #the leftmost bar line for each data series
		for i in range(0, n): #for each series
			pts=[]
			for j in range(0, len(left_points)):
				p=left_points[j]+bar_width*i
				pts.append(p)
		    
			series_lpts[i]= pts
		
		allRects=[]
		for s in range(0, len(toPlot)):
			rects=ax.bar(series_lpts[s], toPlot[s], bar_width, color=colors[s], linewidth=0.25, yerr=errortoPlot[s], ecolor='k')
			allRects.append(rects)
		
		ax.set_xticks(mid_points) # removing this line doesn't make the ticks go away and it puts them at default positions (i.e. starting at zero)
		ax.set_xticklabels(xcats, rotation=45, fontsize=5)
		ax.axhline(y=0, color='k')
		legenddata=[i[0] for i in allRects]
		
		self.dataseries=legenddata
		self.legends=legends
		
	    
def modifyLegends(legends, genecounts=[], pvalues=[], rvalues=[]):
		
	newlegends=[]
	
	for i in range(0, len(legends)):
		basename= legends[i]
		if pvalues!=[]:
			if i==0:
				newname= '%s (%s)' %(basename, genecounts[i])
			else:
				newname='%s (%s) p=%1.2e' %(basename, genecounts[i], pvalues[i-1])
		
			newlegends.append(newname)
		if rvalues!=[]:
			newname= '%s (r=%1.2f)' % (basename, rvalues[i])
			newlegends.append(newname)
			
		else: #no pvalues, just add the gene count data
			newname= '%s (%s)' %(basename, genecounts[i])
			newlegends.append(newname)
		
	return newlegends
			
def getaxsides(axcode):
	
	if 'L' in axcode:
		L='on'
	else:
		L='off'
	
	if 'R' in axcode:
		R='on'
	else:
		R='off'
	
	if 'T' in axcode:
		T='on'
	else:
		T='off'
	
	if 'B' in axcode:
		B= 'on'
	else:
		B='off'
	
	
	return T,B,L,R

def unPickle(pickledObjects):
	'''this returns either an unpickled object or a list of unpickled object
	One tiny issue that this creates is that now otype is not listed for all things. We will assume that otype
	needs to be the same for a list of objects.
	What if we are giving it a pickled list? Then there will still be one object but it will actually be a list'''
	
	if len(pickledObjects)==1:
		pickledObject= pickledObjects[0]
		f= open(pickledObject, 'r')
		d= pickle.load(f)
		f.close()
	
	else:
		d=[]
		for i in pickledObjects:
			f= open(i, 'r')
			o= pickle.load(f)
			f.close()
			d.append(o)
	return d

def getOrder(d):
	#return a list of paired indices to be plotted together
	num= len(d.legends)
		 #change here to deal with multiple pairs of libraries:
	plotPositions=[] #make this into a list of tuples of the form ([x,y], n) where x are the x values, y values and n is the spot in the grid
	    
	range_y= range(0, num)
	range_x= range(0, num)
	    
	counter=1
	for i in range_y:
		for j in range_x:
			pair= [i,j]
			if i==j:
				counter+= len(range_x)-j
				break
			else:
				plotPositions.append((pair, counter))
				counter+=1
	return plotPositions


def formatIt(d, plot, method):
	'''
	Here define the conventional formatting for the different types of plots that I'm creating
	'''
	
	if type(d)==list: #this is a list of objects
		keyobj= d[0]
	else:
		keyobj=d
	
	if method== 'MultiScatter':
		return
	
	thisAxis= plot.axes[plot.currentax]
	
	#change all the axes lines to lw 0.5 and the ticks to 0.25 for consistency
	for axis in ['top', 'bottom', 'left', 'right']:
		thisAxis.spines[axis].set_linewidth(0.5)
		
	thisAxis.tick_params(width=0.25) #don't understand why I can't get the tick lines to be shorter. The width is changing OK?

	legendloc= LegendLocs[args['legend_loc']]
	fs= args['axisfont']
	ls= args['tickfont']
	legf= args['legendfont']
	
	if method=='Box':
		plot.label_axes(thisAxis, fsize=fs)
		plot.draw_legends(thisAxis, plot.dataseries, plot.legends, legendloc, fsize=legf)
		plot.clear_extra_axes(thisAxis, axcode='L')
		plot.set_tick_params(thisAxis, axcode='L', lsize=ls)
		plot.set_limits(thisAxis, ymin=args['ymin'], ymax=args['ymax'])

	elif method=='ScatterHist': #for this one, plot.currentax is actually the axis, not the reference to the axis. too confusing?
		plot.label_axes(thisAxis, fsize=fs)
		plot.set_limits(thisAxis, ymin=args['ymin'], ymax=args['ymax'], xmin=args['xmin'], xmax=args['xmax'])
		xaxis= plot.axes[2]
		yaxis= plot.axes[3]
		plot.set_limits(xaxis, xmin=args['xmin'], xmax=args['xmax'])
		plot.set_limits(yaxis, ymin=args['ymin'], ymax=args['ymax'])
		plot.draw_legends(xaxis, plot.dataseries, plot.legends, legendloc, scatterpoints=1, markerscale=5, fsize=legf)
		plot.set_tick_params(thisAxis, axcode='LBTR', lsize=ls)
	
	elif method=='Scatter':
		plot.label_axes(thisAxis, fsize=fs)
		plot.clear_extra_axes(thisAxis, axcode='LB')
		plot.set_tick_params(thisAxis, axcode='LB', lsize=ls)
		if args['autolimits']==False:
			plot.set_limits(thisAxis, xmin=args['xmin'], xmax=args['xmax'], ymin=args['ymin'], ymax=args['ymax'])

		if args['hidelegenddots']==False:
			plot.draw_legends(thisAxis, plot.dataseries, plot.legends, legendloc, scatterpoints=1, markerscale=3, fsize=legf)
		else:
			plot.draw_legends(thisAxis, plot.dataseries, plot.legends, legendloc, scatterpoints=1, markerscale=0, fsize=legf)
		
	elif method=='SimpleBar':
		plot.label_axes(thisAxis, fsize=fs)
		plot.set_limits(thisAxis, ymin=args['ymin'], ymax=args['ymax'])
		plot.clear_extra_axes(thisAxis, axcode='LB')
		plot.draw_legends(thisAxis, plot.dataseries, plot.legends, legendloc, fsize=legf)
		plot.set_tick_params(thisAxis, axcode='LB', lsize=ls)
		loc= plticker.MultipleLocator(base=10)
		thisAxis.yaxis.set_major_locator(loc)
	
	elif (method=='Violin') or (method=='ViolinSig'):
		plot.label_axes(thisAxis, fsize=fs)
		plot.set_limits(thisAxis, ymin=args['ymin'], ymax=args['ymax'])
		plot.clear_extra_axes(thisAxis, axcode='LB')
		loc= plticker.MultipleLocator(base=15)
		thisAxis.yaxis.set_major_locator(loc)
		plot.set_tick_params(thisAxis, axcode='LB', lsize=ls)
		plot.draw_legends(thisAxis, plot.dataseries, plot.legends, legendloc, fsize=legf)
		
	elif method=='CDF':
		plot.label_axes(thisAxis, fsize=fs)
		plot.set_limits(thisAxis, xmin=args['xmin'], xmax=args['xmax'])
		plot.clear_extra_axes(thisAxis, axcode='LB')
		plot.set_tick_params(thisAxis, axcode='LB', lsize=ls)
		plot.draw_legends(thisAxis, plot.dataseries, plot.legends, legendloc, scatterpoints=1, markerscale=2, fsize=legf)
		
	elif (method=='BarHist') or (method=='LineHist') or (method=='BarLineHist') or (method=='StackedHist'):
		plot.label_axes(thisAxis, fsize=fs)
		plot.set_limits(thisAxis, xmin=args['xmin'], xmax=args['xmax'])
		plot.clear_extra_axes(thisAxis, axcode='LB')
		plot.set_tick_params(thisAxis, axcode='LB', lsize=ls)
		plot.draw_legends(thisAxis, plot.dataseries, plot.legends, legendloc, scatterpoints=1, markerscale=2, fsize=legf)
	
	elif method=='BarHist':
		plot.label_axes(thisAxis, fsize=fs)
		plot.set_limits(thisAxis, xmin=args['xmin'], xmax=args['xmax'])
		plot.clear_extra_axes(thisAxis, axcode='LB')
		plot.set_tick_params(thisAxis, axcode='LB', lsize=ls)
		plot.draw_legends(thisAxis, plot.dataseries, plot.legends, legendloc, scatterpoints=1, markerscale=2, fsize=legf)
		
	elif method== 'MultiScatter':
		pass
	
	else:
		plot.draw_legends(thisAxis, plot.dataseries, plot.legends, legendloc, fsize=legf)
		plot.label_axes(thisAxis, fsize=fs)


def graphIt(d, method, outname, args):
	'''
	Make the plot. Different plots need slightly different call patterns.
	'''
	
	if type(d)==list: #this is a list of objects
		
		keyobj= d[0]
		#test if the legends of the key object are the same as the other, if not, merge
		firstleg= keyobj.legends
		if firstleg!=None:	
			newlegends= firstleg
			for i in range(1, len(d)):
				if d[i].legends!= firstleg:
					newlegends.extend(d[i].legends)
			d[0].legends=newlegends
		
		#test if the mutants of the key object are the same as the other, if not, merge
		firstleg= keyobj.mutants
		if firstleg!=None:
			newlegends= firstleg
			for i in range(1, len(d)):
				if d[i].legends!= firstleg:
					newlegends.extend(d[i].mutants)
			d[0].mutants=newlegends
	else:
		keyobj=d
	
	function= Methods[method]
	plot= function()
	xstring=' '.join(args['xlabel'])
	ystring=' '.join(args['ylabel'])
	titlestring=' '.join(args['title'])
	
	plot.set_basic_properties(titlestring, xstring, ystring) #this should be called first, that way draw_plot_multi will draw it, call this separately for methods below
	
	if method in AbsAxesTypes: #these ones are made without using the subplot function
		plot.draw_subplot(d, None)
	
	elif method=='SimpleBar':
		return plot
	
	elif keyobj.otype=='multiScattero': #I'm having to distinguish otype and method here since I'm using the same method for some different objects
		pairs= d.order
		num= len(d.legends)
		plot.set_big_ax_params(d.legends, axcode='')
		loc= plticker.MultipleLocator(base=args['xmax'])

		for i in range(0, len(pairs)):
			#make a new transient object to mimic a single plot
			pos= pairs[i][1]
			rval= d.pear_rvalues[i]
			rvaltext='%s=%1.2f' % ('r', rval)
			newo= SingleScatterObject(d, pairs[i])
			newo.otype='multiScattero'
			plot.draw_subplot(newo, (num,num,pos))
			thisAxis= plot.axes[plot.currentax]
			plot.clear_extra_axes(thisAxis, axcode='')
			plot.draw_extra_text(thisAxis, [rvaltext], [(0.5,1.1)], fsize=args['axisfont'])
			plot.set_tick_params(thisAxis, axcode='LB', lsize=args['tickfont'], multi=True)
			thisAxis.tick_params(width=0.25)
			plot.set_limits(thisAxis, ymin=args['ymin'], ymax=args['ymax'], xmin=args['xmin'], xmax=args['xmax'])
			thisAxis.xaxis.set_major_locator(loc)
			thisAxis.yaxis.set_major_locator(loc)
	
		plot.set_big_ax_params(d.legends, axcode='', lsize= args['tickfont'], fsize= args['axisfont'])
		plot.fix_spacing()

	else:
		plot.draw_subplot(d, (1,1,1))
		
	formatIt(d, plot, method)
	plot.save_plots(outname)
  
Methods={'CDF':CDFPlotter, 'BarHist':BarHistPlotter, 'Scatter':ScatterPlotter, 'Box':BoxWhiskerPlotter, 'ScatterHist':ScatterHistPlotter, 'SimpleBar':BarPlotterSimple, 'Violin':ViolinPlotter, 'LineHist':LineHistPlotter, 'BarLineHist':BarLineHistPlotter, 'StackedHist':StackedHistPlotter, 'ViolinSig':ViolinPlotterSig, 'MultiScatter':ScatterPlotter}
AbsAxesTypes= set(['ScatterHist'])
MultiAxesTypes= set(['MultiScatter'])
SubplotAxesTypes= set(['CDF', 'Hist', 'BarHist', 'Scatter', 'Box'])
LegendLocs={'topright':1, 'topleft':2, 'bottomleft':3, 'bottomright':4}

def main(argList):
	
	global args
	
	parser= argparse.ArgumentParser(description= 'parse command line args')
	parser.add_argument('method', help='the method to use, see Methods dictionary for acceptable ones')
	parser.add_argument('outname', help= 'the output prefix')
	parser.add_argument('-pickledObjects',nargs='+', help='pickled Objects to use. Right now it supports more than one')
	parser.add_argument('--test', action='store_true', help='dont expect pickledObject input, it will generate and plot some test data')
	parser.add_argument('-barplotfile', help='a tab-delimited file with values to use for producing a barplot')
	parser.add_argument('-GOcatfile', help='a tab-delimited file containing GO categories or other barplot key values to plot')
	parser.add_argument('-barkey', help='a key to use for sorting the barplot values, such as M1X')
	parser.add_argument('-numcats', help='the number of categories to plot for the bar graph')
	parser.add_argument('-otype', help='object type to make for the test, otherwise the object will already have this attribute')
	parser.add_argument('-title', default='', help='graph title')
	parser.add_argument('-xlabel', nargs='+', default='', help='x-axis label')
	parser.add_argument('-ylabel', nargs='+', default='', help='y-axis label')
	parser.add_argument('-flipAxes', action='store_true', help='option to add if you want to swap the given x and y-axis')
	parser.add_argument('-xmin', type=float, default=-2, help='provide min for x-axis')
	parser.add_argument('-xmax', type=float, default=2, help='provide max for x-axis')
	parser.add_argument('-ymin', type=float, default=-2, help='provide min for y-axis')
	parser.add_argument('-ymax', type=float, default=2, help='provide max for y-axis')
	parser.add_argument('--autolimits', action='store_true', help='set this to override the limits for the axes, currently implmented only for scatterplot')
	parser.add_argument('--hidelegenddots', action='store_true', help='set this to remove the dots in the legends')
	parser.add_argument('-binsize', type=float, default=0.1, help='provide bin width for a histogram')
	parser.add_argument('-beta', type=int, default=4, help='parameter for the smoothing function')
	parser.add_argument('--smooth', action='store_true', help='add if you want the histogram to be smoothed')
	parser.add_argument('--drawbar', action='store_true', help='add if you want the bar edges to be drawn on the histogram')
	parser.add_argument('-figsize', default='full', help='enter figure size (full or half page by Nature Standards)')
	parser.add_argument('-panels', type=float, default=1, help='enter the number of subpanels that will fit in the page the figure will go on')
	parser.add_argument('--percent', action='store_true', help='add this if you want the y-axis to be converted from log2 to percentage change. Only implemented for bar and violin plots')
	parser.add_argument('--draw_diagonal', action='store_true', help='add this to draw a y=x line on a scatter plot')
	parser.add_argument('--draw_zeros', action= 'store_true', help='add an x=0 and y=0 dashed line on a scatter plot')
	parser.add_argument('--draw_slope', action='store_true', help='draw a best fit line to the scatter data')
	parser.add_argument('-tickfont', type=int, help='set the fontsize of the tickmarks')
	parser.add_argument('-axisfont', type=int, help='set the size of the x and ylabel fonts')
	parser.add_argument('-legend_loc', default='topleft', help='keyword specifying a built-in legend location, one of these (topright, topleft, bottomleft, bottomright)')
	parser.add_argument('-ax_lw', type=float, default=0.5, help='linewidth for the main axes')
	parser.add_argument('-colors', nargs='+', help='a list of colors with names matching the ones defined at top. This will override the built-in colors')
	parser.add_argument('-legendfont', type=int, help='fontsize for the legend')
	parser.add_argument('-scattersize', type=float, help='input the scattersize')
	parser.add_argument('-linehistfill', action='store_true', help='add this if you want it to fill the linehist')
	parser.add_argument('-legends', nargs='+', help='allow this to overwrite the previous legends of the object')
	parser.add_argument('-legendhandles', type=float, help='modifies the the rcParams[legend.handlelength]')
	parser.add_argument('-fig_scalerx', default=1, type=float, help='adds an extra scaling factor to the figure size. This is useful if you want the actal sizes between different figures to be similar')
	parser.add_argument('-fig_scalery', default=1, type=float, help='adds an extra scaling factor to the figure size. This is useful if you want the actal sizes between different figures to be similar')
	parser.add_argument('-figsizex', type=int, help='specify the width of figure in mm directly')
	parser.add_argument('-figsizey', type=int, help='specify the height of figure in mm directly')
	parser.add_argument('-barformat', help='write GO if the values should be made negative')
	
	ar=parser.parse_args(args=argList)
	args= vars(ar)
	
	method= args['method']
	outname= args['outname']
	
	if args['colors']!=None:
		colors= [colorConvert[i.split('_')[0]] for i in args['colors']]
		colorstyles=[]
		for i in args['colors']:
			if len(i.split('_'))>1:
				if i.split('_')[1]=='dashed':
					colorstyles.append('--')
			else:
				colorstyles.append('-')
				
		global colors
		global colorstyles
		
	if args['pickledObjects']!=None:
		pickledObjects=args['pickledObjects']
		d= unPickle(pickledObjects)
	
	if args['barplotfile']!=None: #just plot from the barplot file directly
		d= args['barplotfile']
		plot= graphIt(d, method, outname, args)
		plot.draw_subplot(args['barplotfile'], (1,1,1))
		formatIt(d, plot, 'SimpleBar')
		plot.save_plots(outname)
		sys.exit()
		
	if args['test']==True:
		d= generateTestData(args['otype'])
	
	if args['flipAxes']==True:
		d.flipAxes=True
	
	graphIt(d, method, outname, args)
	
if __name__ == '__main__':	
	main(sys.argv[1:])