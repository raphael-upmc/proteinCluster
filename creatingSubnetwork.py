#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import json
import argparse


def getFamilies(family_filename) :
    familySet = set()
    file = open(family_filename,'r')
    for line in file :
        family = line.rstrip()
        familySet.add(family)
    file.close()
    return familySet



def readingOrf2subfam2fam(orf2subfam2family_filename,familySet) :
    subfam2fam = dict()
    file = open(orf2subfam2family_filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        orfname,subfam,fam = line.split('\t')
        if fam in familySet :
            subfam2fam[subfam] = fam
    file.close()
    return subfam2fam


def writtingOutput(network_filename,subfamily2nb,subfam2fam,output_filename) :
    output = open(output_filename,'w')

    for subfam in subfam2fam :
        output.write('#'+subfam+'\t'+subfamily2nb[subfam]+'\n')

    file = open(network_filename,'r')
    for line in file :
        line = line.rstrip()
        subfam1,subfam2,weight = line.split('\t')
        if subfam1 in subfam2fam and subfam2 in subfam2fam :
            fam1 = subfam2fam[subfam1]
            fam2 = subfam2fam[subfam2]
            output.write(subfam1+ ' ('+fam1+')'+'\t'+subfam2+ ' ('+fam2+')'+'\t'+weight+'\n')
    file.close()
    output.close()



if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='from a list of subfamilies, exctract the subnetwork of hmm-hmm comparison')
    parser.add_argument('directory', help='the path to the  protein clustering directory (that contains the config.json file)')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILENAME where the subnetwork will be store')
    parser.add_argument('familyList_filename',help='the path of the FAMILYLIST_FILENAME')

    args = parser.parse_args()


    if os.path.exists(args.familyList_filename) :
        familyList_filename = os.path.abspath(args.familyList_filename)
    else:
        sys.exit(args.familyList_filename+' does not exist, exit')


    network_filename = args.directory+'/'+'hhblits/hhr.network'
    if os.path.exists(network_filename) :
        network_filename = os.path.abspath(network_filename)
    else:
        sys.exit(network_filename+' does not exist, exit')

    config_filename = args.directory+'/'+'config.json'
    if os.path.exists(config_filename) :
        config_filename = os.path.abspath(config_filename)
    else:
        sys.exit(config_filename+' does not exist, exit')



    orf2subfam2family_filename = args.directory+'/'+'orf2subfamily2family.tsv'
    if os.path.exists(orf2subfam2family_filename) :
        orf2subfam2family_filename = os.path.abspath(orf2subfam2family_filename)
    else:
        sys.exit(orf2subfam2family_filename+' does not exist, exit')

        
    with open(config_filename) as f:
        data = json.load(f)
    subfamily2nb = data['clusters']

    print('reading '+familyList_filename+'...')
    familySet = getFamilies(familyList_filename)

    print('reading '+orf2subfam2family_filename+'...')
    subfam2fam = readingOrf2subfam2fam(orf2subfam2family_filename,familySet)


    print('reading '+network_filename+'...'+' and writting output '+args.output_filename+'...')
    writtingOutput(network_filename,subfamily2nb,subfam2fam,args.output_filename)
    print('done')
    sys.exit()
