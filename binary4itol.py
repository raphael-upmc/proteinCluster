#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from ete3 import Tree
import argparse
import random

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='creating a binary annotation file to drag and drop in Itol')
    parser.add_argument('tree_filename', help='the path of the TREE_FILE')
    parser.add_argument('matrix_filename',help='the path of the MATRIX_FILE, this file is a tab-separated file, cols are functions/traits, lines are genomes')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILE, results will be stored in this file')
    parser.add_argument('--name',help='name of the dataset in itol')
    parser.add_argument('--color',help='the path of the COLOR2ANNOTATION_FILE, this file is a tab-separated file, first col is the annotation, second col is the color. First line is skipped')
    args = parser.parse_args()

    if os.path.exists(args.tree_filename) :
        tree_filename = os.path.abspath(args.tree_filename)
    else:
        sys.exit(args.tree_filename+' does not exist, exit')

    if os.path.exists(args.matrix_filename) :
        matrix_filename = os.path.abspath(args.matrix_filename)
    else:
        sys.exit(args.matrix_filename+' does not exist, exit')


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
    annotationSize = 0
    genome2annotation = dict()
    file = open(matrix_filename,'r')
    annotationList = next(file).rstrip().split('\t')
    del(annotationList[0])
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        genome = liste[0]
        annot = '\t'.join( liste[1:] )
        genome2annotation[ genome ] = annot
        annotationSize = len(liste[1:])
    file.close()
    print('# of annotations: '+str(annotationSize))

    # annotation2color
    annotation2color = dict()
    color = randomColorList[random.randint(0,len(randomColorList)-1)]
    if args.color != None : # color file provided
        file = open(color2annotation_filename,'r')
        header = next(file)
        for line in file :
            line = line.rstrip()
            annot,color = line.split('\t')
            annotation2color[annot] = color
        file.close()
    else:
        for annotation in annotationList :
            annotation2color[ annotation ] = color
            
    # Load a tree structure from a newick file.
    genomeMissing = set()
    t = Tree(tree_filename)
    otuSet = set()
    for leaf in t:
        
        if leaf.name not in genome2annotation :
            genomeMissing.add(leaf.name)
            continue
        else:
            otuSet.add(leaf.name)
        
    output = open(args.output_filename,'w')
    output.write('DATASET_BINARY'+'\n')
    output.write('SEPARATOR TAB'+'\n')
    output.write('DATASET_LABEL\t'+name+'\n')
    output.write('COLOR\t'+color+'\n')

    output.write('\n\n')
    output.write('SHOW_LABELS\t1'+'\n')

    output.write('FIELD_SHAPES\t'+'\t'.join(['1'] * annotationSize)+'\n')
    output.write('FIELD_LABELS\t'+'\t'.join(annotationList)+'\n')
    output.write('FIELD_COLORS')
    for annot in annotationList :
        if annot in annotation2color :
            output.write('\t'+annotation2color[annot])
        else:
            output.write('\t'+'#000000')
    output.write('\n')
    
    output.write('\n\n')
    output.write('DATA'+'\n')
    
    for otu in otuSet :

        if otu not in genome2annotation :
            continue
        
        annotation = genome2annotation[ otu ]
        output.write(otu+'\t'+annotation.replace('0','-1')+'\n')

    output.close()


    print('\n\n')
    print(genomeMissing)
    
