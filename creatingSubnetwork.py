#! /usr/bin/env python

import os,sys,re
from collections import defaultdict



familySet = set(['fam000603','fam000096','fam000839','fam000152','fam000058'])

subfam2fam = dict()
orf2subfam2family_filename = '/data7/proteinfams/3.6k.PF/mmseqsProteinClustering/20190129/orf2subfamily2family.tsv'
file = open(orf2subfam2family_filename,'r')
for line in file :
    line = line.rstrip()
    orfname,subfam,fam = line.split('\t')
    if fam in familySet :
        subfam2fam[subfam] = fam
file.close()


output = open('output.dash.cytoscape','w')
network_filename = '/data7/proteinfams/3.6k.PF/mmseqsProteinClustering/20190129/hhblits/hhr.network'
file = open(network_filename,'r')
for line in file :
    line = line.rstrip()
    subfam1,subfam2,weight = line.split('\t')
    fam1 = subfam2fam[subfam1]
    fam2 = subfam2fam[subfam2]
    if subfam1 in subfam2fam and subfam2 in subfam2fam :
        output.write(subfam1+ '('+fam1+')'+'\t'+subfam2+ '('+fam2+')'+'\t'+weight+'\n')
file.close()
output.close()
