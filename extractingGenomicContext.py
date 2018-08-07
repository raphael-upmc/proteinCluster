#! /usr/bin/env python


import os,sys,re
from collections import defaultdict


def kNearestNeighbors(orfName,geneList,k) :
    genomicContextList = list()
    sortedList = sorted(geneList,key=lambda x:x[1])
    for i in range( len(sortedList) ) :
        if sortedList[i][3] != orfName :
            continue
        else:
            if i-k < 0 :
                start = 0
            else :
                start = i-k
            
            if i+k+1 > len(sortedList) :
                end = len(sortedList)
            else :
                end = i+k+1

            for j in range(start,end) :
                centroid = sortedList[i][3]
                if j == i :
                    genomicContextList.append(sortedList[j])
#                    print('centroid\t'+str(sortedList[i]) )                
                else :
#                    print('neighbor\t'+str(sortedList[j]) )
                    neighbor = sortedList[j][2]
                    genomicContextList.append(sortedList[j])                    

        return genomicContextList



def genome2orfOrder(feature_filename) :
    orf2genome = dict()
    orf2scaffold = dict()    
    genome2scaffold2orfList = dict()
    file = open(feature_filename,'r')
    for line in file :
        line = line.rstrip()
        orfName,genome,scaffold,start,end,strand = line.split("\t")
        orfName = orfName.rstrip()
        genome = genome.rstrip()
        orf2genome[orfName] = genome
        orf2scaffold[orfName] = scaffold
        if genome not in genome2scaffold2orfList :
            genome2scaffold2orfList[genome] = dict()

        if scaffold not in genome2scaffold2orfList[genome] :
            genome2scaffold2orfList[genome][scaffold] = list()

        genome2scaffold2orfList[genome][scaffold].append([ int(start) , int(end) , strand , orfName ])
    file.close()
    
    return genome2scaffold2orfList,orf2genome,orf2scaffold




if __name__ == "__main__":
    orfList_filename = sys.argv[1]
    feature_filename = sys.argv[2]
    output_filename = sys.argv[3]

    if len(sys.argv) == 4 :
        k=20
    else:
        k = int(sys.argv[4])

    print(k)
    orfSet = set()
    file = open(orfList_filename,'r')
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        orfSet.add(liste[0])
    file.close()

    genome2scaffold2orfList,orf2genome,orf2scaffold = genome2orfOrder(feature_filename)

    
    
    orf2genomicContext = dict()
    for orf in orfSet :
#        print(orf)
        genome = orf2genome[orf]
        scaffold = orf2scaffold[orf]
        orf2genomicContext[orf] = kNearestNeighbors(orf,genome2scaffold2orfList[genome][scaffold],k)        
        

    output = open(output_filename,'w')
    output.write('centroid'+'\t'+'genome'+'\t'+'scaffold'+'\t'+'start'+'\t'+'end'+'\t'+'strand'+'\t'+'orfname'+'\n')
    for orf in orfSet :
#        print(orf)
        genome = orf2genome[orf]
        scaffold = orf2scaffold[orf]
        for elt in orf2genomicContext[orf] :
            output.write(orf+'\t'+genome+'\t'+scaffold+'\t'+str(elt[0])+'\t'+str(elt[1])+'\t'+elt[2]+'\t'+elt[3]+'\n')

    output.close()
        
