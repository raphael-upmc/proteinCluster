#! /usr/bin/python

""" this script detects families operons """

import scipy.stats as st
import os,sys,re
from collections import defaultdict
import random
import numpy as np
import multiprocessing as mp


def kNearestNeighbors(geneList,k,pair2count) :
    sortedList = sorted(geneList,key=lambda x:x[1])
    for i in range( len(sortedList) ) :
        if i-k < 0 :
            start = 0
        else :
            start = i-k


        if i+k+1 > len(sortedList) :
            end = len(sortedList)
        else :
            end = i+k+1

        for j in range(start,end) :
            centroid = sortedList[i][2]
            if centroid == None :
                continue
            
            if j == i :
                continue
            else :
                neighbor = sortedList[j][2]
                if neighbor == None :
                    continue
                
                pair = '\t'.join( sorted([centroid,neighbor]) )
                pair2count[ pair ] += 0.5


def randomizingfamiliesOrder(scaffold2strand2geneList) :
    ''' how to rendomize families??? '''

    newListRandomized = list()
    for scaffold,strand2geneList in scaffold2strand2geneList.items() :
        for strand,geneList in strand2geneList.items() :
            newListRandomized.extend(geneList)
            
    cpt = 0
    randomizedFamiliesOrder = dict()
    random.shuffle(newListRandomized)
    for scaffold,strand2geneList in scaffold2strand2geneList.items() :
        if scaffold not in randomizedFamiliesOrder :
            randomizedFamiliesOrder[ scaffold ] = defaultdict(list)

        for strand,geneList in strand2geneList.items() :
            for gene in geneList :
                fakeGene = gene
                fakeGene[2] = newListRandomized[cpt][2]
                randomizedFamiliesOrder[scaffold][strand].append(fakeGene)
                cpt += 1
                
    return randomizedFamiliesOrder


def orf2familyFunction() :
    orf2family = dict()
    orf2bin_filename = '/data7/proteinfams/3.6k.PF/3600genomes.4pub.all.prot_bin_subFam_fam_clan.tab'
    file = open(orf2bin_filename,'r')
    for line in file :
        line = line.rstrip()
        orfName,genome,subfamily,family,clan = line.split("\t")
        orf2family[orfName] = family
    file.close()
    return orf2family

def genome2familyOrder(familySet,genomeSet) :
    orf2family = orf2familyFunction()
    genome2scaffold2strand2geneList = dict()
    feature_filename = "3600genomes.4pub.all.genomeFeature"
    file = open(feature_filename,'r')
    for line in file :
        line = line.rstrip()
        orfName,genome,scaffold,start,end,strand = line.split("\t")

        if len(genomeSet) != 0 and genome not in genomeSet :
            continue
        
        if orfName not in orf2family :
            family = None
        else :
            if len(familySet) != 0 and orf2family[orfName] not in familySet :                
                family = None
            else :
                family = orf2family[orfName]

                    
        if genome not in genome2scaffold2strand2geneList :
            genome2scaffold2strand2geneList[genome] = dict()

        if scaffold not in genome2scaffold2strand2geneList[genome] :
            genome2scaffold2strand2geneList[genome][scaffold] = defaultdict(list)

#        print(start,end,family,sep='\t')
        genome2scaffold2strand2geneList[genome][scaffold][strand].append([ int(start) , int(end) , family ])
    file.close()
    return genome2scaffold2strand2geneList

def normalizing(genome2scaffold2strand2geneList) :
    pair2maxWeight = defaultdict(float)
    """ for each pair of families, the max number of co-occurrences expected  """
    for genome,scaffold2strand2geneList in genome2scaffold2strand2geneList.items() :
        for scaffold,strand2geneList in scaffold2strand2geneList.items() :
            for strand,geneList in strand2geneList.items() :
                for i1 in range(len(geneList)) :
                    for i2 in range(len(geneList)) :
                        if i1 == i2 :
                            continue
                        else :
                            gene1 = geneList[i1]
                            gene2 = geneList[i2]
                            if gene1[2] == None or gene2[2] == None :
                                continue
                            else :
                                pair = '\t'.join( sorted([gene1[2],gene2[2]]) )
                                pair2maxWeight[pair] += 0.5
    return pair2maxWeight


def simulation(output_filename,genome2scaffold2strand2geneList) :
    pair2weight = defaultdict(float)
    genome2scaffold2strand2randomizedGeneList = dict()
    for genome,scaffold2strand2geneList in genome2scaffold2strand2geneList.items() :
        genome2scaffold2strand2randomizedGeneList[genome] = randomizingfamiliesOrder(scaffold2strand2geneList)
        for scaffold,strand2geneList in genome2scaffold2strand2randomizedGeneList[genome].items() :
            for strand,geneList in strand2geneList.items() :
                kNearestNeighbors(geneList,k,pair2weight)

    pair2MaxWeightExpected = normalizing(genome2scaffold2strand2randomizedGeneList)
    output = open(output_filename,'w')
    for pair,maxWeightExpected in pair2MaxWeightExpected.items() :
        weight = pair2weight[ pair ]
        output.write(pair+'\t'+str(weight/maxWeightExpected)+'\n')
    output.close()


if __name__ == "__main__":

    genomeSet = set()    
    familySet = set()
    file = open('/home/meheurap/proteinCluster/coreCPR/coreFamiliesAnnotation.tsv','r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        family = liste[1]
        module = liste[0]
        if module != 'specificNonCprBacteriaCore' :
            familySet.add(family)
    file.close()
    print('number of families: '+str(len(familySet))+'\n')

    print('reading genomic architectures...')
    genome2scaffold2strand2geneList = genome2familyOrder(familySet,genomeSet)
    print('done\n')

    ######################
    # observation module #
    ######################    
    print('observation...')
    k=2
    observation_filename = 'output.txt'
    
    pair2weight = defaultdict(float)

    for genome,scaffold2strand2geneList in genome2scaffold2strand2geneList.items() :
        for scaffold,strand2geneList in scaffold2strand2geneList.items() :
            for strand,geneList in strand2geneList.items() :
                kNearestNeighbors(geneList,k,pair2weight)
    print(len(pair2weight))

    pair2MaxWeightExpected = normalizing(genome2scaffold2strand2geneList)
    print(len(pair2MaxWeightExpected))

    output = open('output.txt','w')
    for pair,weight in pair2weight.items() :
        maxWeightExpected = pair2MaxWeightExpected[ pair ]
        output.write(pair+'\t'+str(weight/maxWeightExpected)+'\n')
    output.close()
    print('done\n')

    
    #####################
    # simulation module #
    #####################    
    print("simulation...")
    N = 10
    results = list()
    pool = mp.Pool(processes=10,maxtasksperchild=1) # start 20 worker processes and 1 maxtasksperchild in order to release memory
    for i in range(N) :        
        print(i,flush=False)
        output_filename = str(i)+'.txt'
        results.append( pool.apply_async( simulation, args= (output_filename,genome2scaffold2strand2geneList) ))
    pool.close() # Prevents any more tasks from being submitted to the pool
    pool.join() # Wait for the worker processes to exit

    
    #########################
    # analyzing the results #
    #########################
    print("analizing the results...")
    
    pair2distribution = defaultdict(list)
    for i in range(N) :
        print(i,flush=False)
        simulation_filename = str(i)+'.txt'
        file = open(simulation_filename,'r')
        for line in file :
            line = line.rstrip()
            fam1,fam2,weight = line.split('\t')
            pair = fam1+'\t'+fam2
            pair2distribution[pair].append(float(weight))
        file.close()


    output_filename = 'graph.txt'
    output = open(output_filename,'w')
    file = open(observation_filename,'r')
    for line in file :
        line = line.rstrip()
        fam1,fam2,weight = line.split('\t')
        pair = fam1+'\t'+fam2
        mean = np.mean(pair2distribution[pair])
        std = np.std(pair2distribution[pair])
        zscore  = ( ( float(weight) - float(mean) ) / float(std) )
        pvalue = st.norm.cdf(zscore)
        if pvalue > 0.95 :
            output.write(pair+'\t'+str(weight)+'\t'+str(mean)+'\t'+str(std)+'\t'+str(zscore)+'\t'+str(pvalue)+'\n')
    file.close()
