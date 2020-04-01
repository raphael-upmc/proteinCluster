#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from ete3 import Tree
import argparse
import random

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='creating a colorstrip annotation file to drag and drop in Itol')
    parser.add_argument('tree_filename', help='the path of the TREE_FILE')
    parser.add_argument('annotation_filename',help='the path of the ANNOTATION_FILE, this file is a tab-separated file, first col is the genome, second col is the annotation. First line is skipped')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILE, results will be stored in this file')
    parser.add_argument('--color',help='the path of the COLOR2ANNOTATION_FILE, this file is a tab-separated file, first col is the annotation, second col is the color. First line is skipped')
    parser.add_argument('--name',help='name of the dataset in itol')
    parser.add_argument('--nb',type=int,default=1,help='minimum number of annotations to be shown on itol (default: 1)')

    args = parser.parse_args()

    if os.path.exists(args.tree_filename) :
        tree_filename = os.path.abspath(args.tree_filename)
    else:
        sys.exit(args.tree_filename+' does not exist, exit')

    if os.path.exists(args.annotation_filename) :
        annotation_filename = os.path.abspath(args.annotation_filename)
    else:
        sys.exit(args.annotation_filename+' does not exist, exit')


    if args.color != None :
        if os.path.exists(args.color) :
            color2annotation_filename = os.path.abspath(args.color)
        else:
            sys.exit(args.color+' does not exist, exit')

    if args.name != None :
        name = args.name
    else:
        name = 'Unknown'


    randomColorList = ['#FFFFFF', '#AEFF7D', '#F7A541', '#000000','#e6194B', '#B7B180', '#3F9600', '#717c75', '#000084', '#14e557', '#006B01', '#7D2726', '#C9FFBA', '#F36579', '#B79EFF', '#00f9ff', '#4D0057', '#D15800', '#B4FF00', '#F00000', '#6D6D6D', '#BBD181', '#ABBFE0', '#B8BF4F', '#ff0000', '#EA15D6', '#6081FF', '#89AFFF', '#BFBF4F', '#5EFFA1', '#F6CAF2', '#FFFF00', '#2444FF', '#7529FF', '#A8FFDA', '#ff2828', '#DBA9A8', '#4FFF00', '#FFFF92', '#E3EC5F', '#FEFF81', '#4DB6FF', '#ff00e7', '#FFFFFF', '#FDE4A5', '#FFFF0F','#800000','#9A6324','#808000','#469990','#000075','#f58231','#ffe119','#bfef45','#3cb44b','#42d4f4','#4363d8','#911eb4','#f032e6','#fabebe','#ffd8b1','#fffac8','#aaffc3','#e6beff']
               
    # read annotation
    genome2annotation = dict()
    file = open(annotation_filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        genome,annot = line.split('\t')
        genome2annotation[ genome ] = annot
    file.close()


    annotation2count = defaultdict(int)
    t = Tree(tree_filename,1)
    for leaf in t:
        if leaf.name not in genome2annotation :
            continue
        annotation = genome2annotation[ leaf.name ]
        annotation2count[ annotation ] += 1

    annotationList = list( set(annotation2count.keys()) )
    for annotation in annotationList :
        count = annotation2count[annotation]
        if count <= args.nb :
            print(annotation)
            del( annotation2count[annotation] )
        else:
            continue

        
    # annotation2color
    annotation2color = dict()
    if args.color != None : # color file provided
        file = open(color2annotation_filename,'r')
        header = next(file)
        for line in file :
            line = line.rstrip()
            annot,color = line.split('\t')
            annotation2color[annot] = color
        file.close()
    else:
        annotationList = list( set(annotation2count.keys()) )
        k = len(annotationList)
        selectedColorList = random.sample(randomColorList, k)
        for i in range(k) :
            annotation = annotationList[i]
            color = selectedColorList[i]
            annotation2color[annotation] = color

    # Load a tree structure from a newick file.
    annotationMissing = set()
    genomeMissing = set()
    annotation2colorFinal = dict()

    t = Tree(tree_filename,1)
    otuSet = set()
    for leaf in t:

        if leaf.name not in genome2annotation :
            print(leaf.name)
            genomeMissing.add(leaf.name)
            continue
        
        annotation = genome2annotation[ leaf.name ]
        if annotation in annotation2color :
            otuSet.add(leaf.name)
            color = annotation2color[annotation]
            annotation2colorFinal[annotation] = color
        else:
            print(annotation)
            annotationMissing.add(annotation)


        
    output = open(args.output_filename,'w')
    output.write('DATASET_COLORSTRIP'+'\n')
    output.write('SEPARATOR TAB'+'\n')
    output.write('DATASET_LABEL\t'+name+'\n')
    output.write('COLOR\t#ff0000'+'\n')

    output.write('\n\n')
    output.write('LEGEND_TITLE\tDataset_legend'+'\n')

    output.write('LEGEND_SHAPES\t'+'\t'.join(['1'] * len(annotation2colorFinal))+'\n')
    output.write('LEGEND_LABELS\t'+'\t'.join(annotation2colorFinal.keys())+'\n')
    output.write('LEGEND_COLORS\t'+'\t'.join(annotation2colorFinal.values())+'\n')
    
    output.write('\n\n')
    output.write('DATA'+'\n')
    
    for otu in otuSet :
        if otu not in genome2annotation :
            continue
        
        annotation = genome2annotation[ otu ]
        
        if annotation in annotation2color :
            color = annotation2color[annotation]
            output.write(otu+'\t'+color+'\n')
        else:
            continue
    output.close()

    print(annotationMissing)
    print('\n\n')
    print(genomeMissing)
    
