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
    parser.add_argument('--name',help='name of the dataset in itol')

    args = parser.parse_args()

    if os.path.exists(args.tree_filename) :
        tree_filename = os.path.abspath(args.tree_filename)
    else:
        sys.exit(args.tree_filename+' does not exist, exit')

    if os.path.exists(args.annotation_filename) :
        annotation_filename = os.path.abspath(args.annotation_filename)
    else:
        sys.exit(args.annotation_filename+' does not exist, exit')


    if args.name != None :
        name = args.name
    else:
        name = 'Unknown'


    randomColorList = ['#FFFFFF', '#AEFF7D', '#F7A541', '#000000', '#B7B180', '#3F9600', '#717c75', '#000084', '#14e557', '#006B01', '#7D2726', '#C9FFBA', '#F36579', '#B79EFF', '#00f9ff', '#4D0057', '#D15800', '#B4FF00', '#F00000', '#6D6D6D', '#BBD181', '#ABBFE0', '#B8BF4F', '#ff0000', '#EA15D6', '#6081FF', '#89AFFF', '#BFBF4F', '#5EFFA1', '#F6CAF2', '#FFFF00', '#2444FF', '#7529FF', '#A8FFDA', '#ff2828', '#DBA9A8', '#4FFF00', '#FFFF92', '#E3EC5F', '#FEFF81', '#4DB6FF', '#ff00e7', '#FFFFFF', '#FDE4A5', '#FFFF0F']
               
    # read annotation
    genome2annotation = dict()
    file = open(annotation_filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        genome,size = line.split('\t')
        genome2annotation[ genome ] = float(size)
    file.close()


    # Load a tree structure from a newick file.
    annotationMissing = set()
    genomeMissing = set()
    annotation2colorFinal = dict()

    annotation2count = defaultdict(int)
    t = Tree(tree_filename)
    otuSet = set()
    for leaf in t:

        if leaf.name not in genome2annotation :
            genomeMissing.add(leaf.name)
            continue
        else:
            otuSet.add(leaf.name)
        
    output = open(args.output_filename,'w')
    output.write('DATASET_SIMPLEBAR'+'\n')
    output.write('SEPARATOR TAB'+'\n')
    output.write('DATASET_LABEL\t'+name+'\n')
    output.write('COLOR\t#ff0000'+'\n')

    # output.write('\n\n')
    # output.write('LEGEND_TITLE\tDataset_legend'+'\n')

    # output.write('LEGEND_SHAPES\t'+'\t'.join(['1'] * len(annotation2colorFinal))+'\n')
    # output.write('LEGEND_LABELS\t'+'\t'.join(annotation2colorFinal.keys())+'\n')
    # output.write('LEGEND_COLORS\t'+'\t'.join(annotation2colorFinal.values())+'\n')
    
    output.write('\n\n')
    output.write('DATA'+'\n')
    
    for otu in otuSet :
        size = str( genome2annotation[ otu ] )
        output.write(otu+'\t'+size+'\n')

    output.close()

    print(annotationMissing)
    print('\n\n')
    print(genomeMissing)
    
