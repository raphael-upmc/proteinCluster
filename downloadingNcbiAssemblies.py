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
        #print(liste)

def readingAssemblySummary(assembly_summary_filename,taxId2taxName,taxId2parent,taxaSet,accession2gtdb) :
    accession2filename = dict()
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
        ftp_path = liste[19]
        basename = ftp_path.split('/')[-1]
        lineage_ncbi = gettingFullLineage(taxId,taxId2taxName,taxId2parent)
        if accession in accession2gtdb :
            lineage_gtdb = accession2gtdb[ accession ]
        else:
            lineage_gtdb = 'Na'
        genome_filename = ftp_path+'/'+basename+'_genomic.fna.gz'
        for taxon in taxaSet :
            if re.search(';'+taxon+';',lineage_gtdb) or re.search(','+taxon+',',lineage_ncbi):
                #print(accession+'\t'+lineage_gtdb+'\t'+lineage_ncbi)
                accession2filename[accession] = genome_filename
                accession2ncbiTaxonomy[accession] = lineage_ncbi
            else:
                continue
    file.close()
    return accession2filename,accession2ncbiTaxonomy

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
    accession2gtdb = dict()
    bacteria_filename = 'bac120_metadata.tsv'
    file = open(bacteria_filename,'r')
    #header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        # print(liste)
        accession = liste[54]
        lineage_16S = liste[37]
        lineage_gtdb = liste[16]
        lineage_ncbi = liste[78]
        #print(accession+'\t'+lineage_ncbi+'\t'+lineage_gtdb)
        accession2gtdb[ accession ] = lineage_gtdb
    file.close()
    return accession2gtdb






filename = 'taxaList'
taxaSet = readingListTaxa(filename)

today = date.today()
assembly_summary_filename = 'assembly_summary_genbank_'+str(today)+'.txt'
downloadingAssemblySummary(assembly_summary_filename)

print('loading taxDump database...')
taxId2taxName,taxName2taxId = functionTaxId2taxonName()
taxId2taxRank,taxId2parent = functionTaxId2taxon()
print('done')

accession2gtdb = readingGtdb()

accession2filename,accession2ncbiTaxonomy = readingAssemblySummary(assembly_summary_filename,taxId2taxName,taxId2parent,taxaSet,accession2gtdb)
print('number of assemblies: '+str(len(accession2filename)))


output_filename = 'accession2taxonomy.txt'
output = open(output_filename,'w')
output.write('assembly\tncbiTaxonomy\tgtdbTaxonomy\tftpUrl\n')
directory = '/groups/banfield/users/meheurap/Paula/genomes'
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
        
    output.write(accession+'\t'+ncbiTaxonomy+'\t'+gtdbTaxonomy+'\t'+genome_filename+'\n')
    if os.path.exists(output_genome_filename) :
        continue
    liste = urllib.request.urlretrieve(genome_filename, output_genome_filename)
    print( str(cpt)+'\t'+str(liste[0])+'\t'+str(liste) )

    if not os.path.exists(output_genome_filename) :
        accessionError(accession)

output.close()



print(str(len(accessionError))+' assemblies were not downloaded')
for accession in accessionError :
    print(accession+'\t'+accession2filename[ accession ])
