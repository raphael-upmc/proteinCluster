#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import gzip
import argparse


def functionTaxId2taxonName() :
    taxId2taxName = dict()
    taxName2taxId = dict()
    filename = "/home/meheurap/NCBI/taxdump/names.dmp"
    file = open(filename,"r")
    for line in file :
        line = line.rstrip()
        liste = line.split("\t|\t")
        taxId = liste[0].rstrip()
        taxName = liste[1].rstrip()
        taxName2taxId[taxName] = taxId
        if liste[3] == "scientific name\t|" :
            taxId2taxName[taxId] = taxName
    file.close()
    return taxId2taxName,taxName2taxId

def functionTaxId2taxon() :
    taxId2rank = dict()
    taxId2parent = dict()
    filename = "/home/meheurap/NCBI/taxdump/nodes.dmp"
    file = open(filename,"r")
    for line in file :
        line = line.rstrip()
        liste = line.split("\t|\t")
        taxId = liste[0]
        parentTaxId = liste[1]
        rank = liste[2]
        taxId2rank[taxId] = rank
        taxId2parent[taxId] = parentTaxId
    file.close()
    return taxId2rank,taxId2parent

def gettingRank(taxId,rank,taxId2taxRank,taxId2parent) :
    taxRank = taxId2taxRank[taxId]
    while taxRank != rank and taxId != '1' :
        taxRank = taxId2taxRank[ taxId2parent[taxId] ]
        taxId = taxId2parent[taxId]
        
    return taxId

def gettingFullLineage(taxId,taxId2taxName,taxId2parent) :
    lineage = list()
    try :
        while taxId != '1' :        
            lineage.append(taxId2taxName[ taxId ])
            taxId = taxId2parent[taxId]

        
        return ','.join(list(reversed(lineage)))
    except:
        print('error with'+taxId)
        return 'Na'
    
def gi2taxIdFunction() :
    gi2taxId = dict()
    gi2taxId_filename = '/home/meheurap/NCBI/taxdump/gi_taxid_prot.dmp'
    cpt = 0
    print('reading '+gi2taxId_filename+'...')
    file = open(gi2taxId_filename,'r')
    for line in file :
        line = line.rstrip()
        if cpt % 100000000 == 0 :
            print(cpt)
        cpt += 1
        gi,taxId = line.split()
        gi2taxId[gi] = taxId
    file.close()
    print('done ('+str(len(gi2taxId))+' GIs)')
    return gi2taxId

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='from a list of gi accessions, retrieve the taxonomy of the organism')
    parser.add_argument('taxId_filename', help='the path of the TAXID_FILENAME, one ncbi taxonomy accession per line')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILENAME where the results will be store')

    args = parser.parse_args()
    
    if os.path.exists(args.taxId_filename) :
        taxId_filename = os.path.abspath(args.taxId_filename)
    else:
        sys.exit(args.taxId_filename+' does not exist, exit')

    
    print('loading taxDump database...')
    taxId2taxName,taxName2taxId = functionTaxId2taxonName()
    taxId2taxRank,taxId2parent = functionTaxId2taxon()
    print('done')


    print('writting the results in '+args.output_filename)
    output = open(args.output_filename,'w')
    output.write('taxId'+'\t'+'ncbi_lineage'+'\n')
    file = open(taxId_filename,'r')
    for line in file :
        taxId = line.rstrip()
        lineage = gettingFullLineage(taxId,taxId2taxName,taxId2parent)
        output.write(taxId+'\t'+lineage+'\n')
    file.close()
    output.close()
    print('done')
    
    sys.exit()
    # gi = '1084748677'
    
    # taxId = gi2taxId[gi]
    # lineage = gettingFullLineage(taxId,taxId2taxName,taxId2parent)
    # print(gi+'\t'+taxId+'\t'+lineage)
    # genusId = gettingRank(taxId,'genus',taxId2taxRank,taxId2parent)
    # genusName = taxId2taxName[ genusId ]
