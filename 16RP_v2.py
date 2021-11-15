#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor,wait

def HMMsearch() :
    hmm_database = ''
    directory = ''
    for directory in directoryList :
        for root, dirs, files in os.walk(directory):
            for filename in files :
                basename = filename.replace('.hmm','')
                hmm_filename = root+'/'+filename
                domtblout_filename = '/data7/proteinfams/genomicContext/domtblout/'+basename+'.domtblout'
                if not os.path.exists(domtblout_filename) :
                    cmd = 'hmmsearch -E 1e-10 --cpu 1 --domtblout '+domtblout_filename+' '+hmm_filename+' '+fasta_filename+' >/dev/null 2>/dev/null'
                    status = os.system(cmd)
                    print(str(status)+'\t'+cmd)

def readingHMM(domtblout_filename) :
    orf2hmm = defaultdict(list)
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
        orf2hmm[ orf ].append( hmm , float(evalue) , float(bitscore) , orfCover , hmmCover , length, hmmLength )
    file.close()
    return orf2hmm


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='extracting the 16 ribosomal proteins')
    parser.add_argument('protein_filename', help='the path of the FASTA_PROTEIN_FILE')
    parser.add_argument('orf2bin_filename',help='the path of the ORF2BIN_FILE, this file is a tab-separated file, first col is the orf, second col is the genome. First line is skipped')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILE, results will be stored in this file')
    parser.add_argument('output_summary_filename',help='the path of the OUTPUT_SUMMARY_FILE, results will be stored in this file')
    parser.add_argument('--cpu',type=int,default=6,help='number of CPUs (default: 6)')
    parser.add_argument('--rp14',action='store_true',default=False,help='only consider 14RPs for Archaea (default: False)')
    parser.add_argument('--low-memory',action='store_true',default=False,help='low-memory (longer) (default: False)')
    
    args = parser.parse_args()
    
    if os.path.exists(args.protein_filename) :
        protein_filename = os.path.abspath(args.protein_filename)
    else:
        sys.exit(args.protein_filename+' does not exist, exit')

    if os.path.exists(args.orf2bin_filename) :
        orf2bin_filename = os.path.abspath(args.orf2bin_filename)
    else:
        sys.exit(args.orf2bin_filename+' does not exist, exit')

    cpu = args.cpu
    
    output_summary_filename = os.path.abspath(args.output_summary_filename)
    
    output_aln_filename = os.path.abspath(args.output_filename)
    cwd = '/'.join(output_aln_filename.split('/')[:-1])
    
    print('protein_filename: '+protein_filename)
    print('orf2bin_filename: '+orf2bin_filename)
    print('output_aln_filename: '+output_aln_filename)
    print('output_summary_filename: '+output_summary_filename)
    print('number of CPUs: '+str(cpu))
    print('current working directory: '+cwd)
    print('14RP: '+str(args.rp14))
    print('low-memory: '+str(args.low_memory))
    
    folder = cwd+"/16RP_results"
    if os.path.exists(folder) :
        sys.exit(folder+" already exists, remove it first")

    os.mkdir(folder)
