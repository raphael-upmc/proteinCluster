#! /usr/bin/env python

from Bio import SeqIO
import os,sys,re
from collections import defaultdict
import argparse

def checkingMSA(msa_filename,coverage) :
    cpt = total = 0
    seqList = list()
    for seq_record in SeqIO.parse(msa_filename,'fasta') :
        total += 1
        l = len(seq_record)
        nb = seq_record.seq.count('-')
        if ( float(nb) / float(l) ) >= coverage :
            print(seq_record.id+'\t'+str(float(nb) / float(l)) )
            cpt += 1
        else :
            seqList.append(seq_record)
            continue
    SeqIO.write(seqList,msa_filename+'.removingGappedSeq','fasta')


if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='from a fasta file, perform a phylogenetic tree with FastTree')
    parser.add_argument('fasta_filename', help='the path of the FASTA_FILENAME')
    parser.add_argument('--directory',type=str,default='./',help='DIRECTORY where the results will be stored (default: current directory)')
    parser.add_argument('--coverage',type=float,default=0.5,help='the maximum ratio of gaps a sequence should have (default: 0.5)')
    parser.add_argument('--cpu',type=int,default=1,help='number of CPUs (default: 1)')

    args = parser.parse_args()
    
    if os.path.exists(args.fasta_filename) :
        fasta_filename = os.path.abspath(args.fasta_filename)
    else:
        sys.exit(args.fasta_filename+' does not exist, exit')

    if os.path.exists(args.directory) :
        directory = os.path.abspath(args.directory)
    else:
        directory = os.path.abspath('./')
        
    print(fasta_filename)        
    print(directory)
    print(args.cpu)

    #########################
    ######## MAFFT ##########
    #########################

    mafft_filename = directory+'/'+os.path.basename(fasta_filename)+'.mafft'
    cmd = 'mafft --auto --thread '+str(args.cpu)+' '+fasta_filename+' > '+mafft_filename+' 2>/dev/null'
    print(cmd)
    os.system(cmd)

    ##########################
    ######## TRIMAL ##########
    ##########################


    trimal_filename = mafft_filename+'.trimal'
    html_filename = trimal_filename+'.html'
    cmd = 'trimal -fasta -gappyout -in '+mafft_filename+' -out '+trimal_filename+' >/dev/null' #+' -htmlout '+html_filename
    print(cmd)
    os.system(cmd)


    print('checkingMsa')
    checkingMSA(trimal_filename,args.coverage)
    trimal_filename = trimal_filename+'.removingGappedSeq'

    ############
    # FastTree #
    ############

    fasttree_filename = trimal_filename+'.fasttree'
    cmd = 'FastTree -out '+fasttree_filename+' '+trimal_filename+' 2>/dev/null'
    print(cmd)
    os.system(cmd)

    sys.exit()
