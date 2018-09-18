#! /usr/bin/python

""" a script to convert a hmmsearch domtblout into a domainsHitFile recquired for DAMA software, the output file is a 9 columns tab-separated file (evalue,startPro,endPro,namePro,pfamAcc,startPfam,endPfam,bitscore,namePfam """

import os,sys,re
from collections import defaultdict


def convert(pfam_filename,output_filename) :
    output = open(output_filename,"w")
    cpt = 0
    orf2architecture = defaultdict(list)
    file = open(pfam_filename,"r")
    for line in file :
        line = line.rstrip()
        if re.match("#",line) :
            continue
        cpt += 1
        liste = line.split()
        orfName = liste[0]
        length = liste[2]
        pfamName = liste[3]
        pfamAccession = liste[4].split(".")[0]
        pfamLength = liste[5]
        evalue = liste[6]
        bitscore = liste[7]
        cEvalue = liste[11] # conditional Evalue
        iEvalue = liste[12] # independant Evalue

        hmmS = liste[15]
        hmmE = liste[16]

        aliS = liste[17]
        aliE = liste[18]

        envS = liste[19]
        envE = liste[20]
        newList = [cEvalue,aliS,aliE,orfName,pfamAccession,hmmS,hmmE,bitscore,pfamName]
        output.write("\t".join(newList)+"\n")
    file.close()
    output.close()
    
pfam_filename = sys.argv[1]
output_filename = sys.argv[2]
convert(pfam_filename,output_filename)
