#! /usr/bin/python3.6

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import pickle
from annotation import *


if __name__ == "__main__":


    print("loading pickle...")
    t1 = datetime.now()   
    dataset = pickle.load( open( "/home/meheurap/script/annotation.p", "rb" ) ) # Load the dictionary back from the pickle file.
    t2 = datetime.now()
    t3 = t2 - t1 
    print("done in "+str(t3.seconds)+" sec.")


    cpt = 0
    orf2accessions = defaultdict(list)
    hmmsearch_filename = sys.argv[1]
    print(hmmsearch_filename)
    file = open(hmmsearch_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.search(r'^#',line) :
            continue

        cpt += 1
        
        liste = line.split()
        orf = liste[0]
        length = liste[2]
        hmm = liste[3]
        ko = liste[3].split(".")[0]
        hmmLength = liste[5]
        evalue = liste[6]
        bitscore = liste[7]
        cEvalue = liste[11] # conditional Evalue
        iEvalue = liste[12] # independant Evalue

        hmmS = liste[15]
        hmmE = liste[16]

        aliS = liste[17]
        aliE = liste[18]

        envS = liste[19]
        envE = liste[20]

        orfCover = float( int(envE) - int(envS) + 1 ) / float(length)
        hmmCover = float( int(hmmE) - int(hmmS) + 1 ) / float(hmmLength)
        if orfCover > 0.4 and hmmCover > 0.4 :
            orf2accessions[orf].append( ( hmm , float(evalue) , float(bitscore) , orfCover , hmmCover , length, hmmLength ) )
    file.close()


    orf_filename = hmmsearch_filename+'.orfList'
    output = open(orf_filename,'w')
    for orf,liste in orf2accessions.items() :
        output.write(orf+'\n')
    output.close()

    basename = hmmsearch_filename
    context_filename = basename+'.genomicContext'
    feature_filename = '/home/meheurap/proteinCluster/genomicContext/3600genomes.4pub.all.genomeFeature'
    cmd = 'extractingGenomicContext.py '+orf_filename+' '+feature_filename+' '+context_filename+' -k 5'
    print(cmd)
    os.system(cmd)


    
    annotation_filename = basename+'.genomicContext.annotation'
    output = open(annotation_filename,'w')

    centroidSet = set()
    file = open(context_filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        centroid = liste[0]
        genome = liste[1]
        orf = liste[6]
        if orf in orf2accessions :
            tag = 'target'
        else:
            tag = '-'
            
        if orf not in dataset.orfName2orfObjet :
            output.write(genome+'\t'+'Na'+'\t'+orf+'\t'+coord+'\t'+'Na'+'\t'+'Na'+'\t'+'Na'+'\t'+tag+'\n')
            continue


        coord = liste[3]+'..'+liste[4]+' ('+liste[5]+')'    
        orfObject = dataset.orfName2orfObjet[orf]
        genome = orfObject.bin
        taxonomy = dataset.bin2taxonomy[ genome ]
        family = orfObject.family

        if len(orfObject.pfam) == 0 :
            pfam = 'Na'
        else :
            architecture = list()
            for domainInfoList in sorted( orfObject.pfam,key=itemgetter(0,1) ) :
                architecture.append(domainInfoList[2].name+" ("+domainInfoList[2].accession+")")
            pfam = ' + '.join(architecture)
        
        if centroid not in centroidSet :
            centroidSet.add(centroid)
            if len(centroidSet) > 0 :
                output.write('\n')
        
        if orfObject.kegg == 'Na' :
            output.write(genome+'\t'+taxonomy+'\t'+orf+'\t'+coord+'\t'+family+'\t'+'Na'+'\t'+pfam+'\t'+tag+'\n')
        else:
            output.write(genome+'\t'+taxonomy+'\t'+orf+'\t'+coord+'\t'+family+'\t'+orfObject.kegg[0].description+' ('+orfObject.kegg[0].accession+')'+'\t'+pfam+'\t'+tag+'\n')
    file.close()
    output.close()
