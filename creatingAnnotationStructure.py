#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import json
from Bio import SeqIO
import pickle
import numpy as np
from datetime import date, time, datetime
from operator import itemgetter

class Family :
    def __init__(self,name) :
        self.name = name
        self.orfList = dict()

    def __len__(self) :
        """ return the number of ORFs in the protein cluster """
        return len(self.orfList)

    def oneLineAnnotation(self) :
        """ sumarizing ProteinCluster annotation with one line """

        annot2count = {"sequenceLength" : [] , "signalP" : 0 , "TMHMM" : [] , "pfam" : defaultdict(int) , "keggBigCategory" : defaultdict(int) , "keggCategory" : defaultdict(int) , "keggPathway" : defaultdict(int) , "keggAccession" : defaultdict(int) , "keggDescription" : defaultdict(int) , "psort" : defaultdict(int), "cazy" : defaultdict(int) }

        unknown_set = set( [ "hypothetical","Uncharacterized","seg","unannotated" ] )
        for orfName,orfObject in self.orfList.items() :            
            # sequence length
            annot2count["sequenceLength"].append(len(orfObject))

            # psort annotation
            annot2count["psort"][orfObject.psort] += 1

            # Cazy annotation
            annot2count["cazy"][orfObject.cazy] += 1

            # signalP
            if orfObject.signalp :
                annot2count["signalP"] += 1
                
            # TMHMM
            annot2count["TMHMM"].append(orfObject.tmhmm)

            # Pfam
            if len(orfObject.pfam) == 0 :
                annot2count["pfam"]["No domains detected"] += 1
            else :
                architecture = list()
                for domainInfoList in sorted( orfObject.pfam,key=itemgetter(0,1) ) :
                    architecture.append(domainInfoList[2].name+" ("+domainInfoList[2].accession+")")
                annot2count["pfam"][" + ".join(architecture)] += 1

            # Kegg
            if orfObject.kegg == 'Na' :
                annot2count["keggBigCategory"]["No KEGG annotation"] += 1
                annot2count["keggCategory"]["No KEGG annotation"] += 1
                annot2count["keggPathway"]["No KEGG annotation"] += 1
                annot2count["keggAccession"]["No KEGG annotation"] += 1
                annot2count["keggDescription"]["No KEGG annotation"] += 1
            else :
                annot2count["keggBigCategory"][ orfObject.kegg[0].bigCategory ] += 1
                annot2count["keggCategory"][ orfObject.kegg[0].category ] += 1
                annot2count["keggPathway"][ orfObject.kegg[0].pathway ] += 1
                annot2count["keggAccession"][ orfObject.kegg[0].accession ] += 1
                annot2count["keggDescription"][ orfObject.kegg[0].description ] += 1

        # sum-up the results
        psort = sorted(annot2count["psort"].items(),key=itemgetter(1),reverse=True)[0]
        psort = psort[0]+" ("+str( format( float(psort[1])/float(len(self)) , '.2f' ) )+")"

        cazy = sorted(annot2count["cazy"].items(),key=itemgetter(1),reverse=True)[0]
        cazy = cazy[0]+" ("+str( format( float(cazy[1])/float(len(self)) , '.2f' ) )+")"

        pfam = sorted(annot2count["pfam"].items(),key=itemgetter(1),reverse=True)[0]
        pfam = pfam[0]+" ("+str( format( float(pfam[1])/float(len(self)) , '.2f' ) )+")"

        keggBigCategory = sorted(annot2count["keggBigCategory"].items(),key=itemgetter(1),reverse=True)[0]
        keggBigCategory = keggBigCategory[0]+" ("+str( format( float(keggBigCategory[1])/float(len(self)) , '.2f' ) )+")"

        keggCategory = sorted(annot2count["keggCategory"].items(),key=itemgetter(1),reverse=True)[0]
        keggCategory = keggCategory[0]+" ("+str( format( float(keggCategory[1])/float(len(self)) , '.2f' ) )+")"

        keggPathway = sorted(annot2count["keggPathway"].items(),key=itemgetter(1),reverse=True)[0]
        keggPathway = keggPathway[0]+" ("+str( format( float(keggPathway[1])/float(len(self)) , '.2f' ) )+")"

        keggAccession = sorted(annot2count["keggAccession"].items(),key=itemgetter(1),reverse=True)[0]
        keggAccession = keggAccession[0]+" ("+str( format( float(keggAccession[1])/float(len(self)) , '.2f' ) )+")"

        keggDescription = sorted(annot2count["keggDescription"].items(),key=itemgetter(1),reverse=True)[0]
        keggDescription = keggDescription[0]+" ("+str( format( float(keggDescription[1])/float(len(self)) , '.2f' ) )+")"

        result = self.name+"\t"+str(int(len(self)))+"\t"+str(int(np.median(annot2count["sequenceLength"])))+"\t"+str(format( float(annot2count["signalP"])/float(len(self)) , '.2f'))+"\t"+str(int(np.median(annot2count["TMHMM"])))+"\t"+psort+"\t"+pfam+"\t"+keggAccession+"\t"+keggDescription+'\t'+keggPathway+"\t"+keggCategory+"\t"+keggBigCategory+"\t"+cazy
        return result
    

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
        self.familyList = dict()
        self.subfamilyList = dict()
        
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

    def addingFamily(self,filename) :
        file = open(filename,'r')
        header = next(file)
        for line in file :
            line = line.rstrip()
            orf,family = line.split('\t')
            self.orfList[orf].family = family
            if family not in self.familyList :
                self.familyList[ family ] = Family(family)
                self.familyList[ family ].orfList[orf] = self.orfList[orf]
            else:
                self.familyList[ family ].orfList[orf] = self.orfList[orf]
        file.close()

    def addingSubfamily(self,filename) :
        file = open(filename,'r')
        header = next(file)
        for line in file :
            line = line.rstrip()
            orf,subfamily = line.split('\t')
            self.orfList[orf].subfamily = subfamily
            if subfamily not in self.subfamilyList :
                self.subfamilyList[ subfamily ] = Family(subfamily)
                self.subfamilyList[ subfamily ].orfList[orf] = self.orfList[orf]
            else:
                self.subfamilyList[ subfamily ].orfList[orf] = self.orfList[orf]

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
            if not self.orfList[orfName].tmhmm == -1 :
                sys.exit(orfName+" already have TMHMM prediction!")
            else :
                self.orfList[orfName].tmhmm = result 
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


    # reading fasta_filename     
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
        
    # if orf2subfamily_filename available 
    if 'orf2subfamily_filename' in data and os.path.exists( data['orf2subfamily_filename'] ) :
        print("adding subfamilies...")
        dataset.addingSubfamily(data['orf2subfamily_filename'])

    # if orf2family_filename available 
    if 'orf2family_filename' in data and os.path.exists( data['orf2family_filename'] ) :
        print("adding families...")
        dataset.addingFamily(data['orf2family_filename'])

    # if signalP_filename available
    if 'signalp_filename' in data and os.path.exists( data['signalp_filename'] ) :
        print("adding signalP...")
        dataset.addingSignalP(data['signalp_filename'])

    # if pfam_filename available
    if 'pfam_filename' in data and os.path.exists( data['pfam_filename'] ) :
        print("adding Pfam...")
        dataset.addingPFAM(data['pfam_filename'])

    # if tmhmm_filename available
    if 'tmhmm_filename' in data and os.path.exists( data['tmhmm_filename'] ) :
        print("adding TMHMM...")
        dataset.addingTMHMM(data['tmhmm_filename'])

    
    # print(len(dataset.orfName2orfObjet))
    # print("saving dataset object in pickle...")
    # pickle.dump( dataset, open( "annotation.p","wb" ) ) # Save a dictionary into a pickle file.
    


    output_filename = "family2annotation_test.txt"
    output = open(output_filename,"w")
    output.write("\t".join( ["family","nb","seqLength","signalP","TMHMM","psort","Pfam",'keggAccession','keggDescription','keggPathway','keggCategory','keggBigCategory','Cazy'])+"\n")
    for family,familyObject in dataset.familyList.items() :
        output.write(familyObject.oneLineAnnotation()+"\n")
    output.close()
    
