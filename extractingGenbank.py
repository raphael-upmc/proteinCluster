#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
from operator import itemgetter
import argparse


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='extracting the proteins, the contigs and/or the features of a genbank file')
    parser.add_argument('genbank_filename', help='the path of the GBK_FILE')
    parser.add_argument('--proteins',,action='store_true',default=False,help='extracting the proteins')
    parser.add_argument('--contigs',,action='store_true',default=False,help='extracting the contigs')
    parser.add_argument('--features',,action='store_true',default=False,help='extracting the features')

    args = parser.parse_args()


    if os.path.exists(args.genbank_filename) :
        gbk_filename = os.path.abspath(args.genbank_filename)
    else:
        sys.exit(args.genbank_filename+' does not exist, exit')

    protein2sequence = dict()
    seqList = list()
    for record in SeqIO.parse(gbk_filename,'genbank') :
        print(record.name)
        seqList.append(record)
    
        if record.name not in scaffold2nb2orf:
            scaffold2nb2orf[record.name] = dict()
    
        for seq_record in record.features: 
            if seq_record.type == 'CDS':
                if 'protein_id' in seq_record.qualifiers :
                    protein2sequence[ '>'+seq_record.qualifiers['protein_id'][0]+'|'+seq_record.qualifiers['product'][0] ] = seq_record.qualifiers['translation'][0]
                    #scaffold2nb2orf[record.name][seq_record.location.start] = seq_record.qualifiers['protein_id'][0]
                if seq_record.type == 'ORIGIN':
                    print(seq_record.type)


    if contigs :
        SeqIO.write(seqList,gbk_filename+'.contigs.fna','fasta')

    if proteins :
        fasta_filename = gbk_filename+'.proteins.faa'
        output = open(fasta_filename,'w')
        for defline,sequence in protein2sequence.items() :
            output.write(defline+'\n')
            output.write(sequence+'\n')
        output.close()

    
    sys.exit()

