#! /usr/bin/env python


import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import argparse



if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='splitting a fasta file into k files')
    parser.add_argument('fasta_filename', help='the path of the FASTA_FILE')
    parser.add_argument('nb',type=int,help='the number of sub fasta files you want to create')

    args = parser.parse_args()
    
    if os.path.exists(args.fasta_filename) :
        fasta_filename = os.path.abspath(args.fasta_filename)
    else:
        sys.exit(args.fasta_filename+' does not exist, exit')

    nb = args.nb
    if nb < 1 :
        sys.exit('error! '+str(nb)+' must be a positive')

    if os.path.exists('fasta') :
        sys.exit('fasta already exists, remove it first')
    else :
        os.mkdir('fasta')

    print('Fasta filename: '+fasta_filename)
    print('Number of bins: '+str(nb))

    cpt = 0
    for seq_record in SeqIO.parse(fasta_filename , "fasta") : # get the number of sequences
        cpt += 1
    print('\nNumber of sequences in '+fasta_filename+': '+str(cpt)+'\n')

    ###################
    # splitting fasta #
    ###################

    binSeq = ( cpt / nb ) + 1
    cpt = 0
    binNb = 1
    output_fasta_filename = "fasta"+"/"+str(binNb)+".faa"
    output = open(output_fasta_filename,"w")
    for seq_record in SeqIO.parse(fasta_filename , "fasta") :
        if cpt < binSeq :
            cpt += 1
            SeqIO.write(seq_record,output,"fasta")
        else :
            output.close()
            print('Bin: '+str(binNb)+'\t'+str(cpt)+' sequences')
            binNb += 1
            cpt = 1
            output_fasta_filename = "fasta"+"/"+str(binNb)+".faa"
            output = open(output_fasta_filename,"w")
            SeqIO.write(seq_record,output,"fasta")
    output.close()
    print('Bin: '+str(binNb)+'\t'+str(cpt)+' sequences')
