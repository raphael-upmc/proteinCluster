#! /usr/bin/env python


import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import argparse

def bestHit(accessionList) :
    maxi = maxi_hit = 0
    range2score = dict()
    for i in range( len(accessionList) ) :
        hit = accessionList[i]
        queryCover = hit[3]
        hmmCover = hit[4]
        bitscore = hit[2]
        hmm = hit[0]
        score = float(bitscore)
            
        range2score[ i ] = score
        if score > maxi :
            maxi = score
            maxi_hit = i
            besthit = hit
    return range2score,maxi_hit,besthit


def readingHMM(domtblout_filename,bitscoreThreshold) :
    orf2accessions = defaultdict(list)
    file = open(domtblout_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.search(r'^#',line) :
            continue

        liste = line.split()
        orf = liste[0]
        length = liste[2]
        hmm = liste[3]
        ko = liste[3].split(".")[0]
        hmmLength = liste[5]
        evalue = liste[6]
        bitscore = liste[7]

        if float(bitscore) < float(bitscoreThreshold) :
            continue
        
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
    return orf2accessions
    

def runningHMM(domtblout_filename,hmm_filename,fasta_filename) :
    orf2accessions = defaultdict(list)
    orf2besthit = dict()
    cmd = 'hmmsearch -E 1e-3 --cpu 6 --domtblout '+domtblout_filename+' '+hmm_filename+' '+fasta_filename+' >/dev/null 2>/dev/null'
    print(cmd)
    status = os.system(cmd)
    return cmd,status





if __name__ == "__main__":    
    parser = argparse.ArgumentParser(description='this script runs hmmsearch, parses the resulting domtblout output file and extracts the fasta sequences')
    parser.add_argument('hmm_filename', help='the path of the HMM_FILENAME')
    parser.add_argument('--bitscoreThreshold',type=float,default=0,help='bitscore cutoff')
    
    args = parser.parse_args()

    if os.path.exists(args.hmm_filename) :
        hmm_filename = os.path.abspath(args.hmm_filename)
    else:
        sys.exit(args.hmm_filename+' does not exist, exit')

    if args.bitscoreThreshold < 0 :
        sys.exit('bitscoreCutoff has to be a positive value')
    else :
        bitscoreThreshold = float(args.bitscoreThreshold)
        

    ggkbase_filename = '/data7/proteinfams/ggkbase_201907/raph-all-proteins_20190722.txt.cleaned'
    domtblout_filename = '/data7/proteinfams/ggkbase_201907_domtblout/'+os.path.basename(hmm_filename)+'_ggkbase.domtblout'
    fasta_domtblout_filename = '/data7/proteinfams/ggkbase_201907_domtblout/'+os.path.basename(hmm_filename)+'_ggkbase.domtblout.faa'


    print('ggkbase fasta filename: '+ggkbase_filename)
    print('hmm filename: '+hmm_filename)
    print('hmmsearch output filename: '+domtblout_filename)
    print('fasta output filename: '+fasta_domtblout_filename)
    print('bitscore cutoff: '+str(bitscoreThreshold))
    print()

    cmd,status = runningHMM(domtblout_filename,hmm_filename,ggkbase_filename)
    print('command line: '+cmd)
    print('status: '+str(status))

    if status != 0 :
        sys.exit('\nERROR EXIT\n')

    orf2accessions = readingHMM(domtblout_filename,bitscoreThreshold)
    print()


    print('extracting fasta sequences...')
    seqList = list()
    for record in SeqIO.parse(ggkbase_filename,'fasta') :
        if record.id in orf2accessions :
            seqList.append(record)
    SeqIO.write(seqList,fasta_domtblout_filename,'fasta')
    print('done')
