#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import numpy as np
import argparse

''' this script parses the domtblout kegg hmm file and keep the best hit  '''

def bestHit(accessionList,hmm2nc,coverage) :
    maxi = maxi_hit = 0
    range2score = dict()
    for i in range( len(accessionList) ) :
        hit = accessionList[i]
        queryCover = hit[3]
        hmmCover = hit[4]
        bitscore = hit[2]
        hmm = hit[0]
        nc = hmm2nc[ hmm ]

        if coverage :
            score = float(hmmCover) * float(queryCover) * float(bitscore)
        else:
            score = float(bitscore)
            
        range2score[ i ] = score
        if score > maxi :
            maxi = score
            maxi_hit = i
            besthit = hit
    return range2score,maxi_hit,besthit



if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='this script parses the domtblout kegg hmm file and keep the best hit')
    parser.add_argument('hmmsearch_filename', help='the path of the HMMSEARCH_FILENAME')
    parser.add_argument('output_filename', help='the path of the OUTPUT_FILENAME')
    parser.add_argument('--coverage',action='store_true',default=False,help='taking into account the coverage')
    
    args = parser.parse_args()

    if os.path.exists(args.hmmsearch_filename) :
        hmmsearch_filename = os.path.abspath(args.hmmsearch_filename)
    else:
        sys.exit(args.hmmsearch_filename+' does not exist, exit')

    
    kegg_filename = '/groups/banfield/projects/multienv/proteinfams/duduKeggHmm/20150610/ggTables/ko_name_def.tsv'
    ko2pathway_filename = '/groups/banfield/projects/multienv/proteinfams/duduKeggHmm/20150610/ko_pathway_modules.flat.filtered.tab'
    ko2desc_filename = '/groups/banfield/projects/multienv/proteinfams/duduKeggHmm/20150610/ko_name_def.tab'

    hmm2nc = dict()
    hmm2length = dict()
    hmm2sd = dict()
    info_filename = '/groups/banfield/projects/multienv/proteinfams/duduKeggHmm/20150610/ko.gCluster.hmm.thresholds.len_mean_sd.tab'
    file = open(info_filename,'r')
    for line in file :
        line = line.rstrip()
        hmm,NC,length_mean,sd_mean = line.split('\t')
        hmm2nc[ hmm ] = NC
        hmm2length[ hmm ] = length_mean
        hmm2sd[ hmm ] = sd_mean
    file.close()

    
    print('reading '+hmmsearch_filename)
    orf2accessions = defaultdict(list)
    cpt = 0
    file = open(hmmsearch_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.search(r'^#',line) :
            continue

        if cpt % 1000000 == 0 :
            print(cpt)

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
    
        orf2accessions[orf].append( ( hmm , float(evalue) , float(bitscore) , orfCover , hmmCover , length, hmmLength ) )
    file.close()
    print('done\n')


    print('writting ouput....')
    print(str(len(orf2accessions))+' orfs with kegg annotation')
    cpt = 0    
    output = open(args.output_filename,'w')
    output.write('orf'+'\t'+'KEGG'+'\t'+'bitscore'+'\t'+'NC'+'\t'+'evalue'+'\t'+'raw_bitscore'+'\t'+'orfCover'+'\t'+'hmmCover'+'\n')
    for orf,accessionList in orf2accessions.items() :
        if cpt % 100000 == 0 :
            print(cpt)
        cpt += 1
        range2score,maxi,hit = bestHit(accessionList,hmm2nc,args.coverage)
        output.write(orf+'\t'+hit[0]+'\t'+str(range2score[maxi])+'\t'+str(hmm2nc[ hit[0] ])+'\t'+str(hit[1])+'\t'+str(hit[2])+'\t'+str(hit[3])+'\t'+str(hit[4])+'\n')
    output.close()
    print('done')
    sys.exit()

