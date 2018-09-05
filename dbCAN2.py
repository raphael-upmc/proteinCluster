#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from annotation import *
from Bio import SeqIO
import argparse

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='annotating protein sequences with CAZY using dbCAN2')
    parser.add_argument('protein_filename', help='the path of the FASTA_PROTEIN_FILE')
    parser.add_argument('orf2bin_filename',help='the path of the ORF2BIN_FILE, this file is a tab-separated file, first col is the orf, second col is the genome. First line is skipped')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILE, results will be stored in this file')
    parser.add_argument('--cpu',type=int,default=6,help='number of CPUs (default: 6)')

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

    orf2bin = dict()
    file = open(orf2bin_filename,'r')
    for line in file :
        line = line.rstrip()
        orf,genome = line.split('\t')
        orf2bin[ orf ] = genome
    file.close()

    ## ! when no results in hotpep, no file is created....
    genome2seqList = defaultdict(list)
    for record in SeqIO.parse(protein_filename,'fasta') :
        genome = orf2bin [ record.id ]
        genome2seqList[ genome ].append( record )
    print(len(genome2seqList))

    index2genome = dict()
    cpt = 0
    for genome, seqList in genome2seqList.items() :
        cpt += 1
        output_filename = cwd+'/'+'tmp'+'.faa'
        SeqIO.write(seqList,output_filename,'fasta')
        print(str(cpt)+'\t'+genome)
        index2genome[ str(cpt) ] = genome
        cmd = "/data7/proteinfams/dbCAN2/Tools/run_dbcan.py "+output_filename+" protein --dia_eval 0.001 --dia_cpu "+str(cpu)+" --hmm_eval 0.001 --hmm_cpu "+str(cpu)+" --hotpep_cpu "+str(cpu)+" --out_dir "+cwd+"/Results --db_dir /data7/proteinfams/dbCAN2/Tools/example/db --out_pre "+str(cpt)+'.'
        print(cmd)
        status = os.system(cmd)
        print('done with status: '+str(status))

        if cpt == 5 :
            break
    
    sys.exit()
