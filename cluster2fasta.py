#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO

directory = "/data7/proteinfams/3.6k.PF"

family_set = set()
family_selected_filename = sys.argv[1]
file = open(family_selected_filename,"r")
for line in file :
    family = line.rstrip()
    family_set.add(family)
file.close()




output_directory = sys.argv[2]
# creation iBBiGbicluster folder
if os.path.exists(output_directory) :
    sys.exit(output_directory+" already exists! remove it first")
else :
    os.mkdir(output_directory)


orfName2family = dict()
filename = directory+"/"+"3600genomes.4pub.all.prot_bin_subFam_fam_clan.tab"
file = open(filename,"r")
for line in file :
    line = line.rstrip()
    orf,genome,subfamily,family,clan = line.split("\t")
    if family in family_set :
        orfName2family[ orf ] = family
file.close()


family2orfSequence = defaultdict(list)
fasta_filename = directory+"/"+"3600genomes.4pub.faa"
for seq_record in SeqIO.parse(fasta_filename, "fasta"):
    if seq_record.id in orfName2family :
        family = orfName2family[seq_record.id]
        family2orfSequence[ family ].append(seq_record)
    else :
        continue




for family,my_records in family2orfSequence.items() :
    output_fasta_sequence = output_directory+"/"+family+".fa"
    SeqIO.write(my_records, output_fasta_sequence, "fasta")
    
