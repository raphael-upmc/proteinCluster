#! /usr/bin/env python


import os,sys,re
from collections import defaultdict
import urllib.request
from datetime import date
import argparse

def functionTaxId2taxonName() :
    taxId2taxName = dict()
    taxName2taxId = dict()
    filename = "/env/cns/proj/agc/home/rmeheust/scripts/usefulScripts/taxdump/names.dmp"
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
    filename = "/env/cns/proj/agc/home/rmeheust/scripts/usefulScripts/taxdump/nodes.dmp"
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

def gettingFullLineage(taxId,taxId2taxName,taxId2parent) :
    lineage = list()
    try :
        while taxId != '1' :        
            lineage.append(taxId2taxName[ taxId ])
            taxId = taxId2parent[taxId]        
        return ','.join(list(reversed(lineage)))
    except:

#        print('error with '+taxId)
        return 'Na'

def downloadingAssemblySummary(assembly_summary_filename) :
    url = 'ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/assembly_summary_genbank.txt'
    print(assembly_summary_filename)
    if not os.path.exists(assembly_summary_filename) :
        print('downloading assembly_summary_genbank.txt from ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS')
        liste = urllib.request.urlretrieve(url, assembly_summary_filename)
        #print(liste)

def readingAssemblySummary(assembly_summary_filename,taxId2taxName,taxId2parent,taxaSet,accessionSet,accession2gtdb) :
    accessionDetected = set()
    taxaDetected = set()

    accession2filename = dict()
    accession2info = dict()
    accession2ncbiTaxonomy = dict()
    file = open(assembly_summary_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.match('#',line) :
            continue
        liste = line.split('\t')
        accession = liste[0]
        taxId = liste[5]
        speciesTaxId = liste[6]
        genome = liste[7]
        genome_rep = liste[13]
        seq_rel_date = liste[14]
        submitter = liste[16]
        ftp_path = liste[19]
        basename = ftp_path.split('/')[-1]
        lineage_ncbi = gettingFullLineage(taxId,taxId2taxName,taxId2parent)
        if accession in accession2gtdb :
            lineage_gtdb = accession2gtdb[ accession ]
        else:
            lineage_gtdb = 'Na'
        genome_filename = ftp_path+'/'+basename+'_genomic.fna.gz'

        if accession in accessionSet or liste[17] in accessionSet :
            accession2filename[accession] = genome_filename
            accession2ncbiTaxonomy[accession] = lineage_ncbi
            accession2info[accession] = seq_rel_date+'\t'+submitter
            accessionDetected.add(accession)
            accessionDetected.add(liste[17])
        else :
            for taxon in taxaSet :
                if re.search(';'+taxon+';',lineage_gtdb) or re.search(','+taxon+',',lineage_ncbi):
                    #print(accession+'\t'+lineage_gtdb+'\t'+lineage_ncbi)
                    accession2filename[accession] = genome_filename
                    accession2ncbiTaxonomy[accession] = lineage_ncbi
                    accession2info[accession] = seq_rel_date+'\t'+submitter
                    taxaDetected.add(accession)
                else:
                    continue
    file.close()


    if len(accessionSet) > 0 :
        if len(accessionSet - accessionDetected) > 0 :
            print(str(len(accessionSet - accessionDetected))+' accessions were not detected: ')
            for accession in accessionSet - accessionDetected :
                print('\t'+accession)


    if len(taxaSet) > 0 :
        print(str(len(taxaDetected))+' accessions were detected using taxaSet: ')
        

    return accession2filename,accession2ncbiTaxonomy,accession2info

def readingListTaxa(filename) :
    taxaSet = set()
    file = open(filename,'r')
    for line in file :
        taxon = line.rstrip()
        if re.match(r'#',taxon) :
            continue
        taxaSet.add(taxon)
    file.close()
    return taxaSet

def readingAccessionList(filename) :
    accessionSet = set()
    file = open(filename,'r')
    for line in file :
        accession = line.rstrip()
        if re.match(r'#',accession) :
            continue
        accessionSet.add(accession)
    file.close()
    return accessionSet


def readingGtdb() :
    accession2gtdb = dict()
    filenameList = ['/env/export/cns_n02_agc/agc/bank/GTDB/release202/ar122_taxonomy_r202.tsv','/env/export/cns_n02_agc/agc/bank/GTDB/release202/bac120_taxonomy_r202.tsv']
    for filename in filenameList :
        file = open(filename,'r')
        for line in file :
            line = line.rstrip()
            liste = line.split('\t')
            accession = liste[0].split('_')[-1]
            #print(accession)
            lineage_gtdb = liste[1]
            accession2gtdb[ 'GCA_'+accession ] = lineage_gtdb
            accession2gtdb[ 'GCF_'+accession ] = lineage_gtdb
        file.close()
    return accession2gtdb




if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='downloading genome assemblies from NCBI')
    parser.add_argument('directory', help='the path of the DIRECTORY')
    parser.add_argument('output_summary_filename',help='the path of the OUTPUT_SUMMARY_FILE, results will be stored in this file')
    parser.add_argument('--accessionList',help='list of accessions')    
    parser.add_argument('--taxaList',help='list of taxa')    

    args = parser.parse_args()
    
    if os.path.exists(args.directory) :
        directory = os.path.abspath(args.directory)
    else:
        sys.exit(args.directory+' does not exist, exit')

    output_summary_filename = os.path.abspath(args.output_summary_filename)


    if args.taxaList == None and args.accessionList == None :
        sys.exit('taxaList_filename and/or accessionList_filename needed, exit')

    print('directory: '+directory)
    print('output_summary_filename: '+output_summary_filename)

    if args.taxaList != None :
        if os.path.exists(args.taxaList) :
            taxaList_filename = os.path.abspath(args.taxaList)
            print('taxaList_filename: '+taxaList_filename)
        else:
            sys.exit(args.taxaList+' does not exist, exit')

    if args.accessionList != None :
        if os.path.exists(args.accessionList) :
            accessionList_filename = os.path.abspath(args.accessionList)
            print('accessionList_filename: '+accessionList_filename)
        else:
            sys.exit(args.accessionList+' does not exist, exit')



    print()
    taxaSet = set()
    if args.taxaList != None :
        taxaSet = readingListTaxa(taxaList_filename)
        print(str(len(taxaSet))+' taxa in '+taxaList_filename)

    accessionSet = set()
    if args.accessionList != None :
        accessionSet = readingAccessionList(accessionList_filename)
        print(str(len(accessionSet))+' accessions in '+accessionList_filename)    
    print()

    today = date.today()
    assembly_summary_filename = 'assembly_summary_genbank_'+str(today)+'.txt'
    downloadingAssemblySummary(assembly_summary_filename)

    print('loading taxDump database...')
    taxId2taxName,taxName2taxId = functionTaxId2taxonName()
    taxId2taxRank,taxId2parent = functionTaxId2taxon()
    print('done')

    accession2gtdb = readingGtdb()

    print()
    accession2filename,accession2ncbiTaxonomy,accession2info = readingAssemblySummary(assembly_summary_filename,taxId2taxName,taxId2parent,taxaSet,accessionSet,accession2gtdb)
    print('number of assemblies found: '+str(len(accession2filename)))
    print()

    output = open(output_summary_filename,'w')
    output.write('assembly\tncbiTaxonomy\tgtdbTaxonomy\tftpUrl\n')

    cpt = 0
    accessionError = set()
    for accession,genome_filename in accession2filename.items() :
        cpt += 1
        output_genome_filename = directory+'/'+accession+'.fna.gz'
        if accession in accession2ncbiTaxonomy :
            ncbiTaxonomy = accession2ncbiTaxonomy[accession]
        else:
            ncbiTaxonomy = 'Na'

        if accession in accession2gtdb :
            gtdbTaxonomy = accession2gtdb[accession]
        else:
            gtdbTaxonomy = 'Na'
        
        output.write(accession+'\t'+ncbiTaxonomy+'\t'+gtdbTaxonomy+'\t'+accession2info[accession]+'\t'+genome_filename+'\n')

        if os.path.exists(output_genome_filename) :
            continue
        if os.path.exists(output_genome_filename.replace('.gz','')) :
            continue

        liste = urllib.request.urlretrieve(genome_filename, output_genome_filename)
        print( str(cpt)+'\t'+str(liste[0])+'\t'+str(liste) )

        if not os.path.exists(output_genome_filename) :
            accessionError.add(accession)

    output.close()


    print()
    print(str(len(accessionError))+' assemblies were not downloaded')
    for accession in accessionError :
        print(accession+'\t'+accession2filename[ accession ])
