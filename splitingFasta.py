#! /usr/bin/env python


import os,sys,re
from collections import defaultdict
from Bio import SeqIO

fasta_filename = sys.argv[1]
nb = int(sys.argv[2])


if os.path.exists('fasta') :
    sys.exit('psort_tmp already exists, remove it first')
else :
    os.mkdir('fasta')

    

cpt = 0
for seq_record in SeqIO.parse(fasta_filename , "fasta") : # get the number of sequences
    cpt += 1


###################
# splitting fasta #
###################

binSeq = ( cpt / nb ) + 1
cpt = 0
binNb = 1
output_fasta_filename = "fasta"+"/"+str(binNb)+".fasta"
output = open(output_fasta_filename,"w")
for seq_record in SeqIO.parse(fasta_filename , "fasta") :
    cpt += 1
    if cpt < binSeq : 
        SeqIO.write(seq_record,output,"fasta")
    else :
        output.close()
        binNb += 1
        cpt = 1
        output_fasta_filename = "fasta"+"/"+str(binNb)+".fasta"
        output = open(output_fasta_filename,"w")
        SeqIO.write(seq_record,output,"fasta")
output.close()
