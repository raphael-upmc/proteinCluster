#! /usr/bin/python

import os,sys,re
from collections import defaultdict
from operator import itemgetter


class Pfam :
    """ a pfam clan """

    def __init__(self,accession,name,clanAccession,clanName,description):
        self.accession = accession
        self.name = name
        self.clanAccession = clanAccession
        self.clanName = clanName
        self.description = description


def isOverlapping(s1,e1,s2,e2,cover_threshold) :
    """
    check if two intervalle overlap over a cover_threshold
    """
    if s2 <= e1 :
        cover1 = ( min([e1,e2]) - max([s1,s2]) + 1.0 ) / ( float(e1) - float(s1 + 1.0) )
        cover2 = ( min([e1,e2]) - max([s1,s2]) + 1.0 ) / ( float(e2) - float(s2 + 1.0) )
        if cover1 > cover_threshold or cover2 > cover_threshold :
            return True
        else :
            return False
    else :
        return False
            
    
pfamAccession2pfamObject = dict()
filename = "/groups/banfield/projects/multienv/proteinfams/NCBI_balanced_dataset/ncbiGenomeDbComprehensive20171212/PFAM/Pfam-A.clans.tsv"
file = open(filename,"r")
for line in file :
    line = line.rstrip()
    accession,clanAccession,clanName,name,description = line.split("\t")
    pfamAccession2pfamObject[accession] = Pfam(accession,name,clanAccession,clanName,description)
file.close()    
        


orfName2proteinArch = defaultdict(list)
dama_filename = sys.argv[1]
file = open(dama_filename,"r")
for line in file :
    line = line.rstrip()
    liste = line.split("\t")
    evalue = liste[0]
    start = liste[1]
    end = liste[2]
    orfName = liste[3]
    pfamAccession = liste[4]
    tag = liste[5]
    orfName2proteinArch[ orfName ].append( (int(start),int(end),pfamAccession,evalue,tag) )
file.close()

output_filename = sys.argv[2]
output = open(output_filename,"w")

count = 0
for orfName,domainList in orfName2proteinArch.items() :
    if len(domainList) == 1 :
        count += 1
        domain = domainList[0]
        output.write(orfName+"\t"+domain[2]+"\t"+str(domain[0])+"\t"+str(domain[1])+"\t"+str(domain[3])+"\n")
    else : # multidomains protein
#        print(orfName)
        coords = sorted( domainList , key=itemgetter(0,1) ) # sorting coordinates

        cpt = 0
        node2cc = dict()
        for node in coords :
            node2cc[cpt]=cpt
            cpt += 1
            
        i = 0
        while i < len(coords) - 1 :
            j = i + 1
            while j < len(coords) :
                coord1 = coords[i]
                start1 = float(coord1[0])
                end1 = float(coord1[1])

                coord2 = coords[j]
                start2 = float(coord2[0])
                end2 = float(coord2[1])

                if isOverlapping(start1,end1,start2,end2,0.70) :
                    node2cc[j] = node2cc[i]
                    j += 1
                else :
                    j += 1
            i += 1

        cpt2domain = defaultdict(list)
        cpt = 0
#        print("\t"+"with overlaps")
        for domain in coords :
            cpt2domain[node2cc[cpt]].append(domain)
            clan = pfamAccession2pfamObject[domain[2]].clanAccession
            clanName = pfamAccession2pfamObject[domain[2]].clanName
#            print("\t\t"+str(domain)+"\t"+str(node2cc[cpt])+"\t"+clan+"\t"+clanName)
            cpt += 1
#        print("\t"+"without overlaps")
        for cpt,domainList in sorted(cpt2domain.items()) :
            domain = sorted(domainList,key=lambda x:float(x[3]))[0]
            clan = pfamAccession2pfamObject[domain[2]].clanAccession
            clanName = pfamAccession2pfamObject[domain[2]].clanName
            output.write(orfName+"\t"+domain[2]+"\t"+str(domain[0])+"\t"+str(domain[1])+"\t"+str(domain[3])+"\n")
#            print("\t\t"+str(domain)+"\t"+str(node2cc[cpt])+"\t"+clan+"\t"+clanName)

print(count)
