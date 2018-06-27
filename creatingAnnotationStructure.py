#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import json
from Bio import SeqIO

class Genome :
    """ a genome class """

    def __init__(self,name) :
        self.name = name
        self.lineage = 'Na'
        self.protein = 'Na'
        self.size = 'Na'
        self.completeness = 'Na'
        self.contamination = 'Na'
    
class Pfam :
    """ a pfam class """

    def __init__(self,accession,name,clanAccession,clanName,description):
        self.accession = accession
        self.name = name
        self.clanAccession = clanAccession
        self.clanName = clanName
        self.description = description

class Kegg :
    """ a kegg class """

    def __init__(self,accession,name,description,bigCategory,category,pathway):
        self.accession = accession
        self.name = name
        self.description = description
        self.bigCategory = bigCategory
        self.category = category
        self.pathway = pathway

class Orf:
    """ classe definissant une orf """

    def __init__(self,name) :
        self.name = name
        self.genome = "Na"
        self.lineage = "Na"
        self.seq = None
        self.family = "Na"
        self.subfamily = "Na"
        self.clan = "Na"
        self.signalp = None
        self.psort = "Na"
        self.pfam = list()
        self.tmhmm = -1
        self.kegg = "Na"
        self.psort = "Na"
        self.cazy = "Na"
        
    def __len__(self):
        return len(self.seq)

    def __str__(self):
        string = self.name+"\n"+"length: "+str(len(self))+"\n"+"family: "+self.family+"\n"+"signalP: "+str(self.signalp)+"\n"+"ggkbase: "+self.ggkbase+"\n"+"Psortb: "+self.psort+"\n"+"Pfam: "+str(self.pfam)+"\n"+"TMHMM: "+str(self.tmhmm)
        return string


class DatasetAnnotation:
    """ classe annotation qui contient toutes les annotations necessaires! """

    def __init__(self,name) :
        """ notre constructeur """
        self.name = name
        self.genomeList = dict()
        self.orfList = dict()
        
    def addingTaxonomy(self,filename) :
        file = open(filename,"r")
        for line in file :
            line = line.rstrip()
            genome,lineage = line.split("\t")
            self.genomeList[ genome ].lineage = lineage
        file.close()

    def addingFastaFile(self,filename) :
        for seq_record in SeqIO.parse(filename, "fasta"):
            self.orfList[seq_record.id] = Orf(seq_record.id)
            self.orfList[seq_record.id].seq = seq_record

    def addingGenome(self,filename) :
        file = open(filename,'r')
        for line in file :
            line = line.rstrip()
            orf,genome = line.split('\t')
            self.orfList[orf].genome = genome
            self.genomeList[genome] = Genome(genome)
        file.close()
    
    def addingSignalP(self,filename) :
        file = open(filename,'r')
        for line in file :
            line = line.rstrip()

            if re.match("#",line) :
                continue

            orfName,result = line.split('\t')

            if self.orfList[orfName].signalp != None :
                sys.exit(orfName+" already have signalP prediction!")
            if result == "Y" :
                self.orfList[orfName].signalp = True
            else :
                self.orfList[orfName].signalp = False
        file.close()

                
    def addingTMHMM(self,filename) :
        file = open(filename,"r")
        header = next(file)
        for line in file :
            line = line.rstrip()            
            liste = line.split()
            orfName = liste[0]
            result = int(liste[4].replace("PredHel=",""))
            if not self.orfName2orfObjet[orfName].tmhmm == -1 :
                sys.exit(orfName+" already have TMHMM prediction!")
            else :
                self.orfName2orfObjet[orfName].tmhmm = result 
        file.close()

    
    def addingPFAM(self,filename) :

        print("\treading Pfam info...")
        pfamAccession2pfamObject = dict()
        info_filename = "/data7/proteinfams/PFAM/Pfam-A.clans.tsv"
        file = open(info_filename,"r")
        for line in file :
            line = line.rstrip()
            accession,clanAccession,clanName,name,description = line.split("\t")
            pfamAccession2pfamObject[accession] = Pfam(accession,name,clanAccession,clanName,description)
        file.close()

        print("\treading Pfam filename...")
        cpt = 0
        orf2architecture = defaultdict(list)
        file = open(filename,"r")
        for line in file :
            line = line.rstrip()
            liste = line.split()
            orfName = liste[0]
            pfamAccession = liste[1]
            start = liste[2]
            end = liste[3]
            cEvalue = liste[4] # conditional Evalue
            self.orfList[orfName].pfam.append( ( int(start) , int(end) , pfamAccession2pfamObject[pfamAccession] , cEvalue ) )
        file.close()
        



if __name__ == "__main__":


    config_filename = '/home/meheurap/script/creatingAnnotationStructure.json'
    name = 'Elusimicrobia'
    
    if not os.path.exists(config_filename) :
        sys.exit(config_filename+' does not exist, exit')
    else:
        config_filename = os.path.abspath(config_filename)
        
    with open(config_filename) as f:
        data = json.load(f)


    dataset = DatasetAnnotation(name)
    print(dataset)
    print(dataset.name)


    # reading config_filename     
    if 'fasta_filename' in data and os.path.exists(data['fasta_filename']) :
        print("reading fasta file...")
        dataset.addingFastaFile(data['fasta_filename'])
    else :
        sys.exit('error no fasta_filename in '+config_filename)

    # if sequence2genome_filename available 
    if 'sequence2genome_filename' in data and os.path.exists( data['sequence2genome_filename'] ) :
        print("adding Genome...")
        dataset.addingGenome(data['sequence2genome_filename'])

    # if genome2taxonomy_filename available 
    if 'genome2taxonomy_filename' in data and os.path.exists( data['genome2taxonomy_filename'] ) :
        print("adding Taxonomy...")
        dataset.addingTaxonomy(data['genome2taxonomy_filename'])

    # if signalP_filename available
    if 'signalp_filename' in data and os.path.exists( data['signalp_filename'] ) :
        print("adding signalP...")
        dataset.addingSignalP(data['signalp_filename'])

    # if pfam_filename available
    if 'pfam_filename' in data and os.path.exists( data['pfam_filename'] ) :
        print("adding Pfam...")
        dataset.addingPFAM(data['pfam_filename'])
    
    
    # print(len(dataset.orfName2orfObjet))
    # print("saving dataset object in pickle...")
    # pickle.dump( dataset, open( "annotation.p","wb" ) ) # Save a dictionary into a pickle file.
    
