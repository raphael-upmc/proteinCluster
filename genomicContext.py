#! /usr/bin/python

import os,sys,re
from collections import defaultdict

family_set = set()
file = open('coreFamiliesAnnotation.tsv','r')
header = next(file)
for line in file :
    line = line.rstrip()
    liste = line.split('\t')
    family = liste[1]
    module = liste[0]
    if module != 'specificNonCprBacteriaCore' :
        family_set.add(family)
file.close()
print(family_set)
print(len(family_set))

scaffold2orfs = defaultdict(set)
genome2scaffold2families = dict()
orf2bin_filename = '/data7/proteinfams/3.6k.PF/3600genomes.4pub.all.prot_bin_subFam_fam_clan.tab'
file = open(orf2bin_filename,'r')
for line in file :
    line = line.rstrip()
    orfName,genome,subfamily,family,clan = line.split("\t")
    scaffold = '_'.join(orfName.split('_')[:-1])
    
    if family in family_set :
        scaffold2orfs[scaffold].add(orfName)
        if genome in genome2scaffold2families :
            genome2scaffold2families[ genome ][ scaffold ].add(family)
        else :
            genome2scaffold2families[ genome ] = defaultdict(set)
            genome2scaffold2families[ genome ][ scaffold ].add(family)
file.close()


pair2weight = defaultdict(int)
for genome,scaffold2families in genome2scaffold2families.items() :
    print(genome)
    for scaffold,familyList in sorted(scaffold2families.items(),key=lambda x:len(x[1]),reverse=True) :
        for family1 in familyList :
            for family2 in familyList :
                if family1 == family2 :
                    continue
                key = '_'.join(sorted([family1,family2]))
                pair2weight[ key ] += 1

output = open('scaffoldCooccurrence.txt','w')
output.write('\t'+'\t'.join(list(family_set))+'\n')
for family1 in family_set :
    output.write(family1)
    for family2 in family_set :
        if family1 == family2 :
            output.write('\t'+'NA')
        else :
            key = '_'.join(sorted([family1,family2]))
            if key in pair2weight :
                output.write('\t'+str(pair2weight[key]/2))
            else :
                output.write('\t'+'0')
    output.write('\n')
output.close()
