#131004 Mary K. Thompson
#script to build and pickle a dictionary using the .obo annotation (go-basic.obo) file provided by SGD as well as the gene->annotation file (gene_association.sgd.gz)
#this version can also build the tree for humans using the FuncAssociate2.0 file from the Roth lab

import sys
import gzip
import numpy
import cPickle as pickle
import datetime
import argparse

convertNamespace={'cellular_component':'C','molecular_function':'F','biological_process':'P'}
namespaces= convertNamespace.values()

def getGOgenes(stanzas, assocfile, file_format):
    '''
    go thru association file and get mapping genes.
    '''
    namespaces= stanzas.keys()
    
    if file_format=='SGD':    
        f=gzip.open(assocfile, 'r')
        for line in f:
            if line.startswith('!'):
                continue
            if line.startswith('SGD'):
                fields= line.strip().split('\t')
                GO_ID=fields[4].split(':')[-1]
                namespace= fields[8]
                gene_name= fields[10].split('|')[0]
                stanzas[namespace][GO_ID]['genes'].add(gene_name)
               
            else:
                print line #I don't think there's anything that doesn't start with these, but good to check
        
        f.close()
    
    elif file_format=='ROTH':
        f= open(assocfile, 'r')
        for line in f:
            if line.startswith('#'):
                continue
            else:
                fields= line.strip('\n').split('\t')
                GO_ID= fields[0].split(':')[-1]
                GO_description= fields[1]
                genes= fields[2].split(' ')
                #these don't have the namespace listed, so search for it:
                for namespace in namespaces:
                    if GO_ID in stanzas[namespace]:
                        stanzas[namespace][GO_ID]['genes']= set(genes)
                
    return stanzas

def initializeID(ID, namespace, stanzas):
    '''
    initialize a new ID
    '''
    
    stanzas[namespace][ID]={}
    stanzas[namespace][ID]['children']=set() #children will be added as they are found
    stanzas[namespace][ID]['genes']= set()
    return stanzas

def updateStanzas(ID, namespace, parents, name, toDelete, stanzas):
    '''
    update info given data in each new stanza pertaining to parent/child relationships
    '''
    
    if toDelete==1:
        return stanzas
    else:
        if ID not in stanzas[namespace]:
            stanzas=initializeID(ID, namespace, stanzas)
                    
        if 'name' not in stanzas[namespace][ID]: #have to add this here b/c we will initialize some IDs without knowing the names
            stanzas[namespace][ID]['name']= name
            
        stanzas[namespace][ID]['parents']= parents #parents only get added at this final step whereas children can get added constantly
            
        for parent in parents:
            if parent not in stanzas[namespace]:
                stanzas=initializeID(parent, namespace, stanzas)
            stanzas[namespace][parent]['children'].add(ID)
                    
    return stanzas

def parseObo(obofile):
    '''
    Main obo file parsing function. Note that I am treating 'part_of' and 'is a' relationships equivalently here-- both as parents of the category
    '''
    
    stanzas={}
    for k in convertNamespace.values():
        stanzas[k]={}
        
    thisID='' #use so that whenever we hit a new ID we clear the 'cache'
    parents=set() #we'll be able to keep track of each ID's parents
    toDelete=0 #if =1, then ID is obsolete, will delete it
    
    f= open(obofile, 'r')
    num_processed=0
    for line in f:
        
        if (line.startswith('id')) and (thisID!=''): #we hit another one and should update last one
            namespace=convertNamespace[cat]
            stanzas=updateStanzas(thisID, namespace, parents, name, toDelete, stanzas)
            
            GO_ID= line.strip().split(':')[-1] #this should be the number only
            thisID=GO_ID
            toDelete=0
            parents=set()
            
        elif line.startswith('namespace'): #this is the process function component part
            cat=line.strip().split(': ')[-1]
            
        elif line.startswith('name'):
            name= line.strip().split(': ')[-1]
            
        elif line.startswith('relationship: part_of'):
            parent=line.strip().split(' ! ')[0].split(':')[-1]
            parents.add(parent)
            
        elif line.startswith('is_a'):
            #check if they exist parents exist yet, if they do, then add children, if not, intialize then add children
            parent=line.strip().split(' ! ')[0].split(':')[-1]
            parents.add(parent)
            
        elif line.startswith('is_obsolete'): #delete this entry from the dictionary
            toDelete=1
        
        elif line.startswith('id'):
            thisID=line.strip().split(':')[-1] #this should only happen on the first initialization
        
        elif line.startswith('[Typedef]'): #we hit the end of the ids
            break

    f.close()
    return stanzas

def performCalcs(stanzas):
    '''
    find out how many GO categories have more than n genes and less than x genes
    '''
    
    for namespace in stanzas:
        print 'namespace', namespace
        l=[]
        #get average number of genes per category
        #get number with <10 genes and >300 genes
        for ID in stanzas[namespace]:
            num_genes= len(stanzas[namespace][ID]['genes'])
            l.append(num_genes)
        print 'total num genes', len(l)
        print 'average num genes', numpy.mean(l)
        print 'median num genes', numpy.median(l)
        
        i=0
        p=0
        j=0
        
        for item in l:
            if item <2:
                j+=1
            elif item <5:
                i+=1
            elif item >300:
                p+=1
        print 'less than 2', j
        print 'less than 5', i
        print 'greater than 300', p
    
def filterStanzas(stanzas, min_genes):
    '''
    right now implemented to delete GO_IDs with 1 or zero assigned genes. 
    '''
    filteredStanzas={}
    min_genes= int(min_genes)
    
    for namespace in stanzas:
        filteredStanzas[namespace]={}
        for ID in stanzas[namespace]:
            if 'genes' in stanzas[namespace][ID]:
                num_genes= len(stanzas[namespace][ID]['genes'])
                if num_genes>=min_genes:
                    filteredStanzas[namespace][ID]= stanzas[namespace][ID]
    
    return filteredStanzas

def addChildren(stanzas):
    '''
    starting with an ID, add the children and grandchildren genes, etc. to its gene set
    '''
    i=0
    for namespace in stanzas:
        for ID in stanzas[namespace]:
            i+=1
            if i%1000==0:
                print 'i', i
            children= stanzas[namespace][ID]['children']
            
            #first get all the IDs
            allchildren=children
            tolook= children #initialize with original children set
            while tolook != set(): #if set is not empty, continue iteractions:
                found=set() #keep track of new children found
                for child in tolook:
                    newchildren= stanzas[namespace][child]['children']
                    allchildren= allchildren.union(newchildren)
                    found= found.union(newchildren) #add new children to ones found this round
                    
                tolook= found #if we don't find new children then it will stop
                
            #now make a new set of genes from all the children:
            allgenes=stanzas[namespace][ID]['genes'] #initialize with genes mapping specifically to that ID
            
            for child in allchildren:
                genes= stanzas[namespace][child]['genes']
                allgenes= allgenes.union(genes)
            
            #print 'all genes', allgenes
            stanzas[namespace][ID]['genes']= allgenes
            
    return stanzas

def main(argList):
   
    parser= argparse.ArgumentParser(description= 'parse command line args')
    parser.add_argument('obofile', help='the go .obo file')
    parser.add_argument('association_file', help='the file containing GO category-> gene id')
    parser.add_argument('outfile', help= 'the output prefix for the GO tree dictionary')
    parser.add_argument('-format', default='SGD', help='the file format of the association file, SGD or ROTH')
        
    ar= parser.parse_args(args=argList)
    args= vars(ar)
    
    obofile= args['obofile']
    assocfile= args['association_file']
    outfile= args['outfile']
    fileformat= args['format']
    
    #obofile, assocfile, outfile
    print 'start time ', datetime.datetime.now().time()
    
    #eligibleIDs= getGOgenes(assocfile)
    #CANT DO THIS. IF WE REMOVE ONES WITHOUT ANNOTATED GENES THEN WE CAN'T TRACE THE GO TREE (LIKE A BREAK IN THE CHAIN)
    
    print 'parsing Obo'
    stanzas= parseObo(obofile)
    
    print 'getting GO genes'
    #get genes annotated to each GO category:
    stanzas= getGOgenes(stanzas, assocfile, fileformat)
    
    print 'merging all genes'
    #add children genes to genesets. This is necessary because otherwise when you query a specific GO category, you only get the direct children and no grandchildren, etc.
    stanzas= addChildren(stanzas)
    
    print 'end graph time ', datetime.datetime.now().time()

    print 'graph complete'
    print 'pickling...'
    f=open(outfile+'.p', 'w')
    pickle.dump(stanzas, f)
    f.close()
    
    print 'end time', datetime.datetime.now().time()
    

if __name__ == '__main__':
    main(sys.argv[1:])