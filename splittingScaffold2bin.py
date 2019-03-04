#! /usr/bin/env python


import os,sys,re
from collections import defaultdict
from Bio import SeqIO


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='from a scaffold2bin file from ggkbase and a fasta file with scaffolds, the script creates one scaffolds file per bin')
    parser.add_argument('scaffold2bin_filename', help='the path of the scaffold2bin_FILENAME in GGKBASE format')
    parser.add_argument('fasta_filename', help='the path of the FASTA_FILENAME that contains the nt scaffolds sequences')
    parser.add_argument('output_directory',help='the path of the OUTPUT_DIRECTORY where the bin files will be stored')
    parser.add_argument('--force',action='store_true',default=False,help='force the files creation if the directory already exists (default: False)')

    args = parser.parse_args()
    
    # checking arguments
    if os.path.exists(args.scaffold2bin_filename) :
        scaffold2bin_filename = os.path.abspath(args.scaffold2bin_filename)
    else:
        sys.exit(args.scaffold2bin_filename+' does not exist, exit')


    if os.path.exists(args.fasta_filename) :
        fasta_filename = os.path.abspath(args.fasta_filename)
    else:
        sys.exit(args.fasta_filename+' does not exist, exit')

        
    if os.path.exists(args.output_directory) :
        directory = os.path.abspath(args.output_directory)
        if not args.force :
            sys.exit(directory+' already exists, remove it first or use option --force to force creating files into')
    else:
        directory = os.path.abspath(args.output_directory)
        os.mkdir(directory)

    genome2taxonomy = dict()
    scaffold2bin = dict()
    file = open(scaffold2bin_filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        scaffold,genome,lineage = line.split('\t')
        scaffold2bin[ scaffold ] = genome
        genome2taxonomy[ genome ] = lineage
    file.close()

    genome2seqList = defaultdict(list)
    for record in SeqIO.parse(fasta_filename,'fasta') :
        genome = scaffold2bin[record.id].replace(' ','_')
        genome2seqList[ genome ].append(record)

    for genome,seqList in genome2seqList.items() :
        fasta_filename = directory+'/'+genome+'.fna'
        SeqIO.write(seqList,fasta_filename,'fasta')
        print(genome+'\t'+genome2taxonomy[ genome ]+'\t'+fasta_filename)
    
