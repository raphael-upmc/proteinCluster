#! /usr/bin/env python

import os,sys,re
from collections import defaultdict


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

    def __init__(self,name,config_filename) :
        """ notre constructeur """
        self.name = name
        self.genomeList = dict()
        self.orfList = dict()
        
    def addingTaxonomy(self) :
        filename = "/home/meheurap/proteinCluster/taxonomy/bin2taxonomy.txt"
        file = open(filename,"r")
        for line in file :
            line = line.rstrip()
            genome,nearestTaxa,lineageNorm = line.split("\t")
            self.bin2taxonomy[ genome ] = lineageNorm
        file.close()

    def addingFastaFile(self) :
        for seq_record in SeqIO.parse(self.fasta_filename, "fasta"):
            self.orfName2orfObjet[seq_record.id] = Orf(seq_record.id)
            self.orfName2orfObjet[seq_record.id].seq = seq_record
        
    def addingSignalP(self,signalP_filename) :
        file = open(signalP_filename,'r')
        for line in file :
            line = line.rstrip()

            if re.match("#",line) :
                continue

            orfName,result = line.split('\t')

            if self.orfName2orfObjet[orfName].signalp != None :
                sys.exit(orfName+" already have signalP prediction!")
            if result == "Y" :
                self.orfName2orfObjet[orfName].signalp = True
            else :
                self.orfName2orfObjet[orfName].signalp = False
        file.close()

                
    def addingTMHMM(self) :
        tmhmm_filename = "/data7/proteinfams/TMHMM/3600genomes.4pub.faa.tmhmm"
        file = open(tmhmm_filename,"r")
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

    def addingPFAM(self,pfamAccession2pfamObject) :

        print("reading Pfam info...")
        pfamAccession2pfamObject = dict()
        filename = "/data7/proteinfams/PFAM/Pfam-A.clans.tsv"
        file = open(filename,"r")
        for line in file :
            line = line.rstrip()
            accession,clanAccession,clanName,name,description = line.split("\t")
            pfamAccession2pfamObject[accession] = Pfam(accession,name,clanAccession,clanName,description)
        file.close()


        cpt = 0
        orf2architecture = defaultdict(list)
        pfam_filename = "/data7/proteinfams/PFAM/3600genomes.4pub.faa.domtblout.domainsHit.dama.removingOverlapping"
        file = open(pfam_filename,"r")
        for line in file :
            line = line.rstrip()
            liste = line.split()
            orfName = liste[0]
            pfamAccession = liste[1]
            start = liste[2]
            end = liste[3]
            cEvalue = liste[4] # conditional Evalue
            self.orfName2orfObjet[orfName].pfam.append( ( int(start) , int(end) , pfamAccession2pfamObject[pfamAccession] , cEvalue ) )
        file.close()
        



if __name__ == "__main__":

    pfam_filename = '/data7/proteinfams/Elusimicrobia/annotation/PFAM/elusimicrobia.faa.domtblout.domainHit.dama.removingOverlapping'

    signalp_filename = '/data7/proteinfams/Elusimicrobia/annotation/SignalP/elusimicrobia.faa.signalP'

    fasta_filename = '/data7/proteinfams/Elusimicrobia/proteinClustering/elusimicrobia.faa'
    sequence2genone_filename = '/data7/proteinfams/Elusimicrobia/proteinClustering/elusimicrobia.txt'
    genome2taxonomy_filename = ''




    dataset = DatasetAnnotation('Elusimicrobia')
    print(dataset)
    print(dataset.name)

    # print("adding Taxonomy...")
    # dataset.addingTaxonomy()
    
    print("reading fasta file...")
    dataset.addingFastaFile()

    print("adding signalP...")
    dataset.addingSignalpAnnotation()

    print("adding Pfam...")
    dataset.addingPfamAnnotation(pfamAccession2pfamObject)
    
    
    # print(len(dataset.orfName2orfObjet))
    # print("saving dataset object in pickle...")
    # pickle.dump( dataset, open( "annotation.p","wb" ) ) # Save a dictionary into a pickle file.
    
