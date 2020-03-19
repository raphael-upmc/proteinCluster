#! /usr/bin/env python


import os,sys,re
from collections import defaultdict

url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt'

summary_filename = 'assembly_summary_genbank.txt'
print(summary_filename)
file = open(summary_filename,'r')
for line in file :
    line = line.rstrip()
    if re.match('#',line) :
        continue
    liste = line.split('\t')
    accession = liste[0]
    if accession not in asmSet :
        continue
        
    output.write(line+'\n')
    taxId = liste[5]
    speciesTaxId = liste[6]
    genome = liste[7]
    genome_rep = liste[13]
    ftp_path = liste[19]
    basename = ftp_path.split('/')[-1]
    protein_filename = ftp_path+'/'+basename+'_protein.faa.gz'
    accession2filename[ accession ].append(protein_filename)
    accession2line[ accession ] = accession+'\t'+genome
file.close()
print(len(accession2filename))
output.close()


accessionError = set()
for accession,liste in accession2filename.items() :
    for filename in liste :
        cmd = 'wget -P '+output_directory+' '+filename
        basename = os.path.basename(filename)
        if not os.path.exists(output_directory+'/'+basename) :
            print(accession)
            print(cmd)
            status = os.system(cmd)
            if status != 0 :
                print(status)
                accessionError.add(accession)
            else :
                continue

for accession in accessionError :
    print(accession2line[ accession ])
