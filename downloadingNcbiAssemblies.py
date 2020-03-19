#! /usr/bin/env python


import os,sys,re
from collections import defaultdict
import urllib.request
from datetime import date

def functionTaxId2taxonName() :
    taxId2taxName = dict()
    taxName2taxId = dict()
    filename = "/groups/banfield/users/meheurap/NCBI/taxdump/names.dmp"
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
    filename = "/groups/banfield/users/meheurap/NCBI/taxdump/nodes.dmp"
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
        print(liste)

def readingAssemblySummary(assembly_summary_filename,taxId2taxName,taxId2parent,taxaSet) :
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
        ftp_path = liste[19]
        basename = ftp_path.split('/')[-1]
        lineage = gettingFullLineage(taxId,taxId2taxName,taxId2parent)
        genome_filename = ftp_path+'/'+basename+'_genomic.fna.gz'
        for taxon in taxaSet :
            if re.search(','+taxon+',',lineage) :
                print(accession+'\t'+lineage)
            else:
                continue
    file.close()
    

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

def readingGtdb() :
    bacteria_filename = 'bac120_metadata.tsv'
    file = open(bacteria_filename,'r')
    #header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        print(liste)
        accession = liste[54]
        lineage_16S = liste[37]
        lineage_gtdb = liste[16]
        lineage_ncbi = liste[78]
        print(accession+'\t'+lineage_ncbi+'\t'+lineage_gtdb)
    file.close()

readingGtdb()
sys.exit()

filename = 'taxaList'
taxaSet = readingListTaxa(filename)

today = date.today()
assembly_summary_filename = 'assembly_summary_genbank_'+str(today)+'.txt'
downloadingAssemblySummary(assembly_summary_filename)

print('loading taxDump database...')
taxId2taxName,taxName2taxId = functionTaxId2taxonName()
taxId2taxRank,taxId2parent = functionTaxId2taxon()
print('done')



liste = readingAssemblySummary(assembly_summary_filename,taxId2taxName,taxId2parent,taxaSet)


sys.exit()

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
