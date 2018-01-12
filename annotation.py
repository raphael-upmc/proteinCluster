#! /usr/bin/python
""" annotation class to store annotation of each protein cluster """

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import pickle
import numpy as np
from datetime import date, time, datetime
from operator import itemgetter



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
        self.bin = "Na"
        self.lineage = "Na"
        self.seq = None
        self.family = "Na"
        self.subfamily = "Na"
        self.clan = "Na"
        self.signalp = None
        self.ggkbase = "Na"
        self.psort = "Na"
        self.pfam = list()
        self.tmhmm = -1
        self.kegg = "Na"
        self.psort = "Na"
        
    def __len__(self):
        return len(self.seq)

    def __str__(self):
        string = self.name+"\n"+"length: "+str(len(self))+"\n"+"family: "+self.family+"\n"+"signalP: "+str(self.signalp)+"\n"+"ggkbase: "+self.ggkbase+"\n"+"Psortb: "+self.psort+"\n"+"Pfam: "+str(self.pfam)+"\n"+"TMHMM: "+str(self.tmhmm)
        return string
    
class ProteinCluster :
    """ classe definissant un cluster de proteine """

    def __init__(self,name):
        self.name=name
        self.orfDict = dict()

    def __len__(self) :
        """ return the number of ORFs in the protein cluster """
        return len(self.orfDict)

    def __str__(self) :
        """ sumarizing ProteinCluster annotation """
        annot2count = {"sequenceLength" : [] , "ggkbase" : defaultdict(int) , "signalP" : 0 , "TMHMM" : [] , "pfam" : defaultdict(int) }
        for orfName,orfObject in self.orfDict.items() :            
            # sequence length
            annot2count["sequenceLength"].append(len(orfObject))

            # ggkbase annotation
            annot2count["ggkbase"][orfObject.ggkbase] += 1

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

        result = "-----------------------------------------------"+"\n"+ \
        self.name+"\t"+str(len(self))+"\n"+ \
        "\t"+"sequenceLength :"+str(np.median(annot2count["sequenceLength"]))+"\n"+ \
        "\t"+"signalP: "+str(float(annot2count["signalP"])/float(len(self)))+"\n"+ \
        "\t"+"TMHMM :"+str(np.median(annot2count["TMHMM"]))+"\n"+ \
        "\t"+"top5 ggkbase annotations:"+"\n"
        cpt = 0
        for annot,count in sorted(annot2count["ggkbase"].items(),key=lambda x:x[1],reverse=True) :
            cpt += 1
            result += "\t"+"\t"+annot+"\t"+str(count)+"\n"
            if cpt == 5 :
                break

        result += "\t"+"top5 Pfam architectures:"+"\n"
        cpt = 0
        for annot,count in sorted(annot2count["pfam"].items(),key=lambda x:x[1],reverse=True) :
            cpt += 1
            result += "\t"+"\t"+annot+"\t"+str(count)+"\n"
            if cpt == 5 :
                break        
        result += "-----------------------------------------------"
        return result

    def oneLineAnnotation(self) :
        """ sumarizing ProteinCluster annotation with one line """
        annot2count = {"sequenceLength" : [] , "ggkbase" : defaultdict(int) , "signalP" : 0 , "TMHMM" : [] , "pfam" : defaultdict(int) , "keggBigCategory" : defaultdict(int) , "keggCategory" : defaultdict(int) , "keggPathway" : defaultdict(int) , "keggAccession" : defaultdict(int) , "keggDescription" : defaultdict(int) , "psort" : defaultdict(int), "CPR" : set() , "non-CPR-Bacteria" : set() , "non-DPANN-Archaea" : set() , "DPANN" : set() }
        unknown_set = set( [ "hypothetical","Uncharacterized","seg","unannotated" ] )
        for orfName,orfObject in self.orfDict.items() :            
            # sequence length
            annot2count["sequenceLength"].append(len(orfObject))

            # psort annotation
            annot2count["psort"][orfObject.psort] += 1

            
            # ggkbase annotation
            if orfObject.ggkbase in unknown_set :
                annot2count["ggkbase"]["unknown"] += 1
            else :
                annot2count["ggkbase"][orfObject.ggkbase] += 1

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


            # taxonomy
            if orfObject.lineage == 'Na' :
                print(orfObject.bin+'\t'+orfObject.lineage)
            else :
                lineageNorm = orfObject.lineage
                taxonomy = lineageNorm.split(',')[-3]

            if taxonomy not in annot2count :
                print(taxonomy+'\t'+lineageNorm)
            else :
                annot2count[taxonomy].add(orfObject.bin)
            
        ggkbase = sorted(annot2count["ggkbase"].items(),key=itemgetter(1),reverse=True)[0]
        ggkbase = ggkbase[0]+" ("+str( format( float(ggkbase[1])/float(len(self)) , '.2f' ) )+")"

        psort = sorted(annot2count["psort"].items(),key=itemgetter(1),reverse=True)[0]
        psort = psort[0]+" ("+str( format( float(psort[1])/float(len(self)) , '.2f' ) )+")"

        
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

        cpr = str(len(annot2count['CPR']))
        nonCprBacteria = str(len(annot2count['non-CPR-Bacteria']))
        dpann = str(len(annot2count['DPANN']))
        nonDpannArchaea = str(len(annot2count['non-DPANN-Archaea']))
        
        result = self.name+"\t"+str(int(len(self)))+"\t"+str(int(np.median(annot2count["sequenceLength"])))+"\t"+str(format( float(annot2count["signalP"])/float(len(self)) , '.2f'))+"\t"+str(int(np.median(annot2count["TMHMM"])))+"\t"+psort+"\t"+ggkbase+"\t"+pfam+"\t"+keggAccession+"\t"+keggDescription+'\t'+keggPathway+"\t"+keggCategory+"\t"+keggBigCategory+"\t"+cpr+"\t"+nonCprBacteria+"\t"+dpann+"\t"+nonDpannArchaea
        return result
    
    def writtingFasta(self,output_filename):
        output = open(output_filename,"w")       
        for orfName,orfObject in self.orfDict.items() :
            SeqIO.write(orfObject.seq,output,"fasta")
        output.close()
            
class DatasetAnnotation:
    """ classe annotation qui contient toutes les annotations necessaires! """

    def __init__(self) :
        """ notre constructeur """
        self.name = "All"
        self.fasta_filename = "/data7/proteinfams/3.6k.PF/3600genomes.4pub.faa"
        self.orf2family_filename = "/data7/proteinfams/3.6k.PF/3600genomes.4pub.all.prot_bin_subFam_fam_clan.tab"
        self.orfName2orfObjet = dict()
        self.clusterName2clusterObject = dict()
        self.bin2taxonomy = dict()

    def addingTaxonomy(self) :
        filename = "/home/meheurap/proteinCluster/taxonomy/bin2nearestTaxaGroup.txt"
        file = open(filename,"r")
        for line in file :
            line = line.rstrip()
            genome,nearestTaxa,lineage,lineageNorm = line.split("\t")
            taxonomy = lineageNorm.split(',')[-3]
            self.bin2taxonomy[ genome ] = lineageNorm
        file.close()

        
    def fillingOrf(self) :
        for seq_record in SeqIO.parse(self.fasta_filename, "fasta"):
            self.orfName2orfObjet[seq_record.id] = Orf(seq_record.id)
            self.orfName2orfObjet[seq_record.id].seq = seq_record
            liste = seq_record.description.split()
            del(liste[0])
            try :
                for i in range(len(liste)) :
                    if re.search("=",liste[i]) :
                        break
                    if i != 0 :
                        annot = " ".join(liste[:i])
                    else :
                        annot = "unannotated"
            except :
                annot = "unannotated"
            self.orfName2orfObjet[seq_record.id].ggkbase = annot
                
            m = re.search('bin=\"([\w\s\.\-]*)\"\s',seq_record.description)
            try :
                genome = m.group(1).rstrip()
                self.orfName2orfObjet[seq_record.id].bin = genome.rstrip()
                self.orfName2orfObjet[seq_record.id].lineage = self.bin2taxonomy[ genome.rstrip() ]
            except :
                print(self.orfName2orfObjet[seq_record.id].bin+'\t'+self.orfName2orfObjet[seq_record.id].lineage)
                continue

                

            
        file = open(self.orf2family_filename,"r")
        for line in file :
            line = line.rstrip()            
            orfName,genome,subfamily,family,clan = line.split("\t")
            self.orfName2orfObjet[orfName].subfamily = subfamily.rstrip()
            self.orfName2orfObjet[orfName].family = family.rstrip()
            self.orfName2orfObjet[orfName].clan = clan.rstrip()

            if self.orfName2orfObjet[orfName].bin == 'Na' :
                self.orfName2orfObjet[orfName].bin = genome.rstrip()
                
            if family not in self.clusterName2clusterObject :
                self.clusterName2clusterObject[family] = ProteinCluster(family)
                self.clusterName2clusterObject[family].orfDict[orfName] = self.orfName2orfObjet[orfName]
            else :
                self.clusterName2clusterObject[family].orfDict[orfName] = self.orfName2orfObjet[orfName]
        file.close()



        
    def addingSignalpAnnotation(self) :
        for root, dirs, files in os.walk("/data7/proteinfams/SignalP/output"):
            for filename in files :
                file = open(root+"/"+filename,"r")
                for line in file :
                    line = line.rstrip()
                    if re.match("#",line) :
                        continue
            
                    liste = line.split()
                    orfName = liste[0]
                    result = liste[9]
                    if self.orfName2orfObjet[orfName].signalp != None :
                        sys.exit(orfName+" already have signalP prediction!")
                    if result == "Y" :
                        self.orfName2orfObjet[orfName].signalp = True
                    else :
                        self.orfName2orfObjet[orfName].signalp = False
                file.close()

                
    def addingTmhmmAnnotation(self) :
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

    def addingPsortAnnotation(self) :
        psort_filename = "/home/meheurap/proteinCluster/coreCPR/cprOnly.psort"
        file = open(psort_filename,"r")
        header = next(file)
        for line in file :
            line = line.rstrip()            
            liste = line.split()
            orfName = liste[0].split()[0]
            localisation = liste[-5]
            if not self.orfName2orfObjet[orfName].psort == "Na" :
                sys.exit(orfName+" already have psort prediction!")
            else :
                self.orfName2orfObjet[orfName].psort = localisation
        file.close()

        
    def addingPfamAnnotation(self,pfamAccession2pfamObject) :
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
        
    def addingKeggAnnotation(self,keggAccession2keggObject) :
        koSet = set()
        kegg_filename = "/data7/proteinfams/3.6k.PF/annotation/keggHMM/3600genomes.keggHMM.report.tab"
        file = open(kegg_filename,'r')
        header = next(file)
        for line in file :
            line = line.rstrip()
            liste = line.split()
            orfName = liste[0]
            KO = liste[2]
            if KO == '-' :
                continue
            name = liste[3]
            description = liste[4]
            reliability = liste[6]
            evalue = liste[7]
            try :
                self.orfName2orfObjet[orfName].kegg = ( keggAccession2keggObject[KO] , evalue , reliability )
            except :
                koSet.add(KO)
        file.close()
        print(koSet)

if __name__ == "__main__":

    print("reading Pfam info...")
    pfamAccession2pfamObject = dict()
    filename = "/data7/proteinfams/PFAM/Pfam-A.clans.tsv"
    file = open(filename,"r")
    for line in file :
        line = line.rstrip()
        accession,clanAccession,clanName,name,description = line.split("\t")
        pfamAccession2pfamObject[accession] = Pfam(accession,name,clanAccession,clanName,description)
    file.close()

    print("reading Kegg info...")
    keggAccession2keggObject = dict()
    filename = "/data7/proteinfams/3.6k.PF/annotation/keggHMM/kegg.info"
    file = open(filename,"r")
    for line in file :
        line = line.rstrip()
        accession,description,bigCategory,category,pathway = line.split("\t")
        liste = description.split(';')
        
        name = liste[0].rstrip()
        description = ''.join(liste[1:]).rstrip()
        
        keggAccession2keggObject[accession] = Kegg(accession,name,description,bigCategory,category,pathway)
    file.close()    

    
    dataset = DatasetAnnotation()
    print(dataset)
    print(dataset.name)

    print("adding Taxonomy...")
    dataset.addingTaxonomy()
    
    print("fillingOrf...")
    dataset.fillingOrf()

    print("adding Psort...")
    dataset.addingPsortAnnotation()
    
    print("adding TMHMM...")
    dataset.addingTmhmmAnnotation()

    print("adding signalP...")
    dataset.addingSignalpAnnotation()

    print("adding Pfam...")
    dataset.addingPfamAnnotation(pfamAccession2pfamObject)

    print("adding Kegg...")
    dataset.addingKeggAnnotation(keggAccession2keggObject)
    
    print(len(dataset.orfName2orfObjet))
    print("saving dataset object in pickle...")
    pickle.dump( dataset, open( "annotation.p","wb" ) ) # Save a dictionary into a pickle file.

    output_filename = "family2annotation.txt"
    output = open(output_filename,"w")
    output.write("\t".join( ["family","nb","seqLength","signalP","TMHMM","psort","Ggkbase","Pfam",'keggAccession','keggDescription','keggPathway','keggCategory','keggBigCategory','CPR','non-CPR-Bacteria','DPANN','non-DPANN-Archaea'])+"\n")
    for family,familyObject in dataset.clusterName2clusterObject.items() :
        output.write(familyObject.oneLineAnnotation()+"\n")
    output.close()
