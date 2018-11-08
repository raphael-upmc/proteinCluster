#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
import multiprocessing as mp
import traceback
from copy import deepcopy

class Taxonomy:
    """ classe definissant la taxonomy d'un genome (herite de la classe genome) """

    def __init__(self,name,color):
        self.name=name
        self.color=color

class ProteinCluster:
    """ classe definissant un cluster de proteine """
    def __init__(self,name):
        self.name=name
        self.orfList = list()
        self.annotation = list()
        self.genomeList = set()

    
class Genome:        
    """ classe definissant un genome """
    def __init__(self,name):
        self.name=name
        self.orfList = list()
        self.taxonomy = Taxonomy

        
class Community:
    """ classe definissant une communaute detectee par la methode de Louvain. Cette classe contient obligatoirement :
    - un nom
    - la liste des familles
    """
    
    def __init__(self,name) :
        self.name = name
        self.familyList = set()
        self.genomeList = set()

    def creatingMatrixImage(self,matrixObject,image_filename,genomeOrder_filename,familyOrder_filename) :
        matrix_filename = "/tmp/"+self.name+".matrix"
        output = open(matrix_filename,"w")
        output.write("\t"+"\t".join(matrixObject.familyList)+"\n")
        for genome in matrixObject.genomeList :
            output.write(genome)
            for family in matrixObject.familyList :
                if family in matrixObject.matrix[ genome ] :
                    if family in self.familyList :
                        output.write("\t"+"2")
                    else :
                        output.write("\t"+"1")
                else :
                    output.write("\t"+"0")
            output.write("\n")
        output.close()
        cmd = "/home/meheurap/script/lowMemoryHeatmap.r "+matrix_filename+" "+image_filename+" "+genomeOrder_filename+" "+familyOrder_filename #+" 2>/dev/null"
        os.system(cmd)
        os.remove(matrix_filename)

    def creatingSubmatrixImage(self,matrixObject,submatrix_filename,image_filename) :      
        output = open(submatrix_filename,"w")
        
        output.write("\t"+"\t".join(list(self.familyList))+"\n")
        for genome in self.genomeList :
            output.write(genome)
            for family in self.familyList :
                if family in matrixObject.matrix[ genome ] :
                    output.write("\t"+"1")
                else :
                    output.write("\t"+"0")
            output.write("\n")
        output.close()

        cmd = "/home/meheurap/script/heatmap.R "+submatrix_filename+" "+image_filename #+" 2>/dev/null"
        os.system(cmd)
        print(cmd)

        
class  Matrix:
    """ classe definissant une matrice presence/absence genome x family """

    def __init__(self):
        self.genomeList = list()
        self.familyList = list()
        self.matrix = defaultdict(dict)
        self.family2genomes = defaultdict(set)
        self.genome2families = defaultdict(set)
        
    def remplir(self,matrix_filename) :
        """ methode qui remplit l'objet Matrix """
        file = open(matrix_filename,"r")
        header = next(file).rstrip().split("\t")
        del(header[0])
        self.familyList = header
        for line in file :
            line = line.rstrip()
            liste = line.split("\t")
            genome = liste[0]
            self.genomeList.append(genome)
            for i in range(1,len(liste)) :
                if liste[i] == "1" :
                    family = self.familyList[i-1]
                    self.matrix[genome][family] = 1
                    self.family2genomes[ family ].add(genome)
                    self.genome2families[ genome ].add(family)
        file.close()
        
    # def effacer(self):
    #     """ method qui efface l'objet matrix """
    #     self.genomeList = del(self.genomeList[:])
    #     self.familyList = del(self.familyList[:])
    #     self.matrix = self.matrix.clear()
        

if __name__ == "__main__":

    min_size = 10 # opt arg 4
    cpu = 6
    genomeOrder_filename = sys.argv[3] # opt arg 3
    familyOrder_filename = sys.argv[4] # opt arg 3
    
    # creating folders
    print("creating folders...")
    directory = sys.argv[5] # opt arg 3
    tmp_dir = directory
    if os.path.exists(tmp_dir) :
        sys.exit(tmp_dir+" already exists, remove it")
    else :
        os.mkdir(tmp_dir)
        os.mkdir(tmp_dir+"/"+"heatmapAll")
        os.mkdir(tmp_dir+"/"+"heatmapModule")
        os.mkdir(tmp_dir+"/"+"heatmapModule"+"/"+"figure")
        os.mkdir(tmp_dir+"/"+"heatmapModule"+"/"+"file")


    matrix_filename = sys.argv[1] # arg1
    genomeFamilyMatrix = Matrix()
    genomeFamilyMatrix.remplir(matrix_filename)

        
    # family2module

    fakeCommunity = Community("fakeModule")
    print("reading family2module file...")
    community2objetCommunity = dict()
    family2cluster_filename = sys.argv[2] # arg2
    file = open(family2cluster_filename,"r")
    for line in file :
        line = line.rstrip()
        family,community = line.split("\t")
        family = family.replace("X","")
        fakeCommunity.familyList.add(family)
        fakeCommunity.genomeList |= genomeFamilyMatrix.family2genomes[family]
        if community not in community2objetCommunity :
            community2objetCommunity[ community ] = Community(community)
            community2objetCommunity[ community ].familyList.add(family)
            community2objetCommunity[ community ].genomeList |= genomeFamilyMatrix.family2genomes[family]
        else :
            community2objetCommunity[ community ].familyList.add(family)
            community2objetCommunity[ community ].genomeList |= genomeFamilyMatrix.family2genomes[family]
    file.close()
    print(str(len(fakeCommunity.familyList))+" families are present in a module")

    
    results = list()
    pool = mp.Pool(processes=cpu,maxtasksperchild=1) # start 20 worker processes and 1 maxtasksperchild in order to release memory

    cpt = 1
    for community,communityObject in sorted(community2objetCommunity.items(),key=lambda x:len(x[1].familyList),reverse=True) :
        if len(communityObject.familyList) < min_size :
            continue
        
        print(str(cpt)+"\t"+community+" "+str(len(communityObject.genomeList))+"x"+str(len(communityObject.familyList)))

        submatrix_image_filename = tmp_dir+"/"+"heatmapModule"+"/"+"figure"+"/"+communityObject.name+"_heatmap.pdf"
        submatrix_filename = tmp_dir+"/"+"heatmapModule"+"/"+"file"+"/"+communityObject.name+".matrix"
        ARGS = (genomeFamilyMatrix,submatrix_filename,submatrix_image_filename)
        results.append( pool.apply_async( communityObject.creatingSubmatrixImage, args=deepcopy(ARGS) ))

        matrix_image_filename = tmp_dir+"/"+"heatmapAll"+"/"+communityObject.name+"_image.tiff"
        ARGS = (genomeFamilyMatrix,matrix_image_filename,genomeOrder_filename,familyOrder_filename)
        results.append( pool.apply_async( communityObject.creatingMatrixImage, args=deepcopy(ARGS) ))

        
        cpt += 1
        
    pool.close() # Prevents any more tasks from being submitted to the pool
    print("pool has been closed")
    pool.join() # Wait for the worker processes to exit
    print("pool has been joined")
    output = list()
    for p in results :
        print(p.ready())
        print(p.successful())
        print(p.wait())
        output.append( p.get() )
    print(output)
