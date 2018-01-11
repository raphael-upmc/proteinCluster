#! /usr/bin/python

import os,sys,re
from collections import defaultdict

familySet = set(sys.argv[1:])
print(familySet)

bin2taxonomy = dict()
filename = "/home/meheurap/proteinCluster/taxonomy/bin2nearestTaxaGroup.txt"
file = open(filename,"r")
for line in file :
    line = line.rstrip()
    genome,nearestTaxa,lineage,lineageNorm = line.split("\t")
    taxonomy = lineageNorm.split(',')[-3]
    bin2taxonomy[ genome ] = lineageNorm
file.close()



family2bin = defaultdict(set)
orf2family_filename = "/data7/proteinfams/3.6k.PF/3600genomes.4pub.all.prot_bin_subFam_fam_clan.tab"
file = open(orf2family_filename,"r")
for line in file :
    line = line.rstrip()            
    orfName,genome,subfamily,family,clan = line.split("\t")
    if family not in familySet :
        continue
    family2bin[ family ].add(genome)
file.close()


for family,binList in family2bin.items() :
    taxa2count = defaultdict(int)
    for genome in binList :
        lineage = bin2taxonomy[ genome ]
        try :
            taxonomy = lineage.split(',')[-4]
        except :
            taxonomy = lineage.split(',')[-3]
        taxa2count[taxonomy] += 1
    print(family+'\t'+str(sorted(taxa2count.items(),key=lambda x:x[1],reverse=True)))
