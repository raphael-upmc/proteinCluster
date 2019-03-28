#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO

familyList_filename = sys.argv[1]
fasta_filename = sys.argv[2]
orf2family_filename = sys.argv[3]

print('reading the list of families you want to extract...')
family_set = set()
file = open(familyList_filename,"r")
for line in file :
    family = line.rstrip()
    family_set.add(family)
file.close()


print('reading the orf2family filename...')
orfName2family = dict()
file = open(orf2family_filename,"r")
for line in file :
    line = line.rstrip()
    orf,family = line.split("\t")
    if family in family_set :
        orfName2family[ orf ] = family
file.close()

print('reading the fasta filename...')
family2orfSequence = defaultdict(list)
for seq_record in SeqIO.parse(fasta_filename, "fasta"):
    if seq_record.id in orfName2family :
        family = orfName2family[seq_record.id]
        family2orfSequence[ family ].append(seq_record)
    else :
        continue


print('writting the outputs...')

for family,my_records in family2orfSequence.items() :
    output_fasta_sequence = family+".faa"
    SeqIO.write(my_records, output_fasta_sequence, "fasta")
    
