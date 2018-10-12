#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import json
from Bio import SeqIO
import pickle
import numpy as np
from datetime import date, time, datetime
from operator import itemgetter
import argparse

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
        self.protein = 0
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
        next(file)
        for line in file :
            line = line.rstrip()
            orf,genome = line.split('\t')
            if orf not in self.orfList :
                continue
            
            self.orfList[orf].genome = genome
            if genome not in self.genomeList :
                self.genomeList[genome] = Genome(genome)
            self.genomeList[genome].protein += 1
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

    def addingCAZY(self,filename) :
        file = open(filename,"r")
        header = next(file)
        for line in file :
            line = line.rstrip()            
            orfName,cazy = line.split('\t')
            self.orfList[orfName].cazy = cazy
        file.close()


    def addingKEGG(self,filename):
        ko2kegg = dict()        
        json_filename = '/data7/proteinfams/3.6k.PF/annotation/keggHMM/ko00000.json'

        with open(json_filename) as f:
            data = json.load(f)

        for kegg in data['children'] :
            keggBigCategory = kegg['name']
            for kegg1 in kegg['children'] :
                keggCategory = kegg1['name']
                for kegg2 in kegg1['children'] :
                    keggPathway = kegg2['name']
                    if 'children' in kegg2 :            
                        for kegg3 in kegg2['children'] :
                            ko = kegg3['name'].split()[0]
                            name = kegg3['name'].split('  ')[1].split(';')[0]
                            description = kegg3['name'].split('  ')[1].split(';')[1].strip()
                            if ko not in ko2kegg :
                                ko2kegg[ko] = Kegg(ko,name,description,keggBigCategory,keggCategory,keggPathway)
                            else:
                                if ko2kegg[ko].accession != ko :
                                    ko2kegg[ko].accession = 'Multiple'

                                if ko2kegg[ko].name != name :
                                    ko2kegg[ko].name = 'Multiple'

                                if ko2kegg[ko].description != description :
                                    ko2kegg[ko].description = 'Multiple'

                                if ko2kegg[ko].bigCategory != keggBigCategory :
                                    ko2kegg[ko].bigCategory = 'Multiple'
                            
                                if ko2kegg[ko].category != keggCategory :
                                    ko2kegg[ko].category = 'Multiple'

                                if ko2kegg[ko].pathway != keggPathway :
                                    ko2kegg[ko].pathway = 'Multiple'

                        
                    else:
                        continue


        koSet = set()
        file = open(filename,'r')
        header = next(file)
        for line in file :
            line = line.rstrip()
            liste = line.split('\t')
            orfName = liste[0]
            KO = liste[1].split('.')[0]
            evalue = liste[4]

            hmmCover = float(liste[7])
            orfCover = float(liste[6])

            if orfCover > 0.7 :
                if hmmCover > 0.7 :
                    reliability = 'complete'
                else:
                    reliability = 'truncatedOrf??'
            else:
                if hmmCover > 0.7 :
                    reliability = 'fusion??'
                else:
                    reliability = 'not_reliable'
                
                    
            try :
                self.orfList[orfName].kegg = ( ko2kegg[KO] , evalue , reliability )
            except :
                koSet.add(KO)
        file.close()
        print(koSet)

        
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
        

    def addingPsort(self,filename) :
        file = open(filename,"r")
        header = next(file)
        for line in file :
            line = line.rstrip()            
            liste = line.split('\t')
            orfName = liste[0].split()[0]
            localisation = liste[-5]
            finalScore = float(liste[-3])
            cytoplasmicScore = float(liste[-9])
            if not self.orfList[orfName].psort == "Na" :
                print(orfName+" already have psort prediction! "+self.orfList[orfName].psort+' '+localisation)

            if self.orfList[orfName].psort != "Na" and self.orfList[orfName].psort != localisation :
                print(orfName+" already have psort prediction and it's different! "+self.orfList[orfName].psort+' '+localisation)
                    
#                    sys.exit(orfName+" already have psort prediction!")
            else :
                if localisation == 'Unknown' and finalScore > 4 and cytoplasmicScore < 4 :
                    localisation = 'Exported'
                    self.orfList[orfName].psort = localisation
                else :
                    self.orfList[orfName].psort = localisation
        file.close()

        


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run all-vs-all hhblits on the subfamilies')
    parser.add_argument('config_filename', help='the path of the CONFIG_FILENAME that contains the paths of the annotation files')
    parser.add_argument('pickle_filename', help='the pickle filename where the results will be store')
    parser.add_argument('--name',default='Unknow',help='the name of the dataset (default: Unknown)')
    parser.add_argument('--family2annotation_filename',help='the output filename where the family2annotation will be store (default: ./FASTA_FILENAME_proteinClutering)')
    
    args = parser.parse_args()

    
    if not os.path.exists(args.config_filename) :
        sys.exit(args.config_filename+' does not exist, exit')
    else:
        config_filename = os.path.abspath(args.config_filename)

    if os.path.exists(args.pickle_filename) :
        sys.exit(args.pickle_filename+' already exists, remove it first')
    else:
        pickle_filename = os.path.abspath(args.pickle_filename)


    if args.family2annotation_filename != None :
        if os.path.exists(args.family2annotation_filename) :
            sys.exit(args.family2annotation_filename+' already exists, remove it first, exit')
        else:
            family2annotation_filename = os.path.abspath(args.family2annotation_filename)


        
    with open(config_filename) as f:
        data = json.load(f)
    dataset = DatasetAnnotation(args.name)

    



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

    # if cazy_filename available
    if 'cazy_filename' in data and os.path.exists( data['cazy_filename'] ) :
        print("adding CAZY...")
        dataset.addingCAZY(data['cazy_filename'])

    # if kegg_filename available
    if 'kegg_filename' in data and os.path.exists( data['kegg_filename'] ) :
        print("adding KEGG...")
        dataset.addingKEGG(data['kegg_filename'])


    # if psort_filename available
    if 'psort_filename' in data and os.path.exists( data['psort_filename'] ) :
        print("adding Psort...")
        dataset.addingPsort(data['psort_filename'])
        
    

    print("saving dataset object in pickle file"+pickle_filename+"...")
    pickle.dump( dataset, open(pickle_filename,"wb" ) ) # Save a dictionary into a pickle file.
    

    if args.family2annotation_filename != None :
        print('creating '+family2annotation_filename+'...')
        output = open(family2annotation_filename,"w")
        output.write("\t".join( ["family","nb","seqLength","signalP","TMHMM","psort","Pfam",'keggAccession','keggDescription','keggPathway','keggCategory','keggBigCategory','Cazy'])+"\n")
        for family,familyObject in dataset.familyList.items() :
            output.write(familyObject.oneLineAnnotation()+"\n")
        output.close()
