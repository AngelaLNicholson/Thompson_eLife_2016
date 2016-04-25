#Author: Boris Zinshteyn
#purpose: Takes a pickled dictionary of ORFs, as well as windows around the start and end of genes, and compiles metagene data in those windows, and plots 2 historams, for start and stop

#version history:
#created 5/25/2010

import sys, pickle, numpy, math
from scipy import *
from scipy import fftpack
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt

################################################
def compile_metagene(pickledorfs, startLeft, startRight, stopLeft, stopRight):
################################################

    f= open(pickledorfs,'r')
    orfs= pickle.load(f) #load the dictionary of {orfName, {'hits', {Location, #hits}}}
    startIndices=range(startLeft, startRight+1)
    allIndices=range(-12, 3000)
    startValues=[0]*len(startIndices)
    stopIndices=range(stopLeft, stopRight+1)
    stopValues=[0]*len(stopIndices)
    allValues=[0]*3000
    for i in range(len(startIndices)):
        index=startIndices[i]
        for orf in orfs.keys():
            if index in orfs[orf]['hits'].keys():
                startValues[i]+=orfs[orf]['hits'][index]
    for k in range(len(allValues)):
        index=allIndices[k]
        for orf in orfs.keys():
            if index in orfs[orf]['hits'].keys():
                allValues[k]+=orfs[orf]['hits'][index]
    for j in range(len(stopIndices)):
        index=stopIndices[j]
        for orf in orfs.keys():
            #dictionary is indexed from front, but we have index from end
            #print orfs[orf]['length']
            adjustedIndex=index+orfs[orf]['length']
            if adjustedIndex in orfs[orf]['hits'].keys():
                #print "index:"+str(index)
                #print "length:"+str(orfs[orf]['length'])
                #print "adjusted Index:"+str(adjustedIndex)
                stopValues[j]+=orfs[orf]['hits'][adjustedIndex]
    return 	startIndices, startValues, stopIndices, stopValues, allValues

def plot_histogram(xLabels, heights, outputPrefix):
    width=1.0 #bar width
    
    plot=plt.bar(xLabels, heights, color='blue', width=width, linewidth=0, align='edge')
    
    
    xLabels=numpy.arange(xLabels[0], xLabels[-1], 3)
    plt.ylabel('Total # of Read 5\' ends')
    plt.xticks(xLabels+width/2., xLabels)
    title=outputPrefix.split('/')[-1]
    plt.title(title)
    outName=outputPrefix+'.pdf'
    plt.savefig(outName, transparent=True, format='pdf')
    plt.clf()
    
def printFFT(vect, outFile):
    f=open(outFile, 'w')
    Y=fft(vect)
    n=len(Y)
    power = abs(Y[1:(n/2)])**2
    nyquist=1./2
    freq=array(range(n/2))/(n/2.0)*nyquist
    period=1./freq
    plt.semilogy(period[1:len(period)], power)
    plt.xlabel('Period [nt]')
    plt.xlim(0, 10)
    plt.xticks(arange(10), arange(10))
    plt.ylabel('|FFT|**2')
    plt.savefig(outFile+'.pdf', transparent=True, format='pdf')
    pickle.dump(Y, f)
    f.close()
    
def main():

    pickledorfs, startLeft, startRight, stopLeft, stopRight, outputPrefix= sys.argv[1:]

    startIndices, startValues, stopIndices, stopValues, allValues = compile_metagene(pickledorfs, int(startLeft), int(startRight), int(stopLeft), int(stopRight))
    
    plot_histogram(startIndices, startValues, outputPrefix+'_start')
    plot_histogram(stopIndices, stopValues, outputPrefix+'_stop')
    
    printFFT(allValues, outputPrefix+'.fft')
    

main()
                                                                                                               
                                                                                                                
