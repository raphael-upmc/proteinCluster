#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio.SeqFeature import SeqFeature, FeatureLocation


orf2color = dict()

scaffold2orf2coord = dict()
scaffold2min = dict()
scaffold2max = dict()

#feature_filename = '/data7/proteinfams/nitrogenase/myxobacteria/genomicContext/myxococcales_nitrogenase.annotation.txt'
feature_filename = '/data7/proteinfams/Elusimicrobia/metabolism/nitrogenase.genomicContext.annotation'
feature_filename = '/data7/proteinfams/BD1-5/genomicContext/BD1-5_nitrogenase.annotation.txt'
feature_filename = '/home/meheurap/YcaO/allison_YcaO/genomicDiagram.txt'
file = open(feature_filename,'r')
for line in file :
    line = line.rstrip()
    if line == '' :
        continue
    #genome,scaffold,coord,orf,length,kegg_acc,kegg,pfam = line.split('\t')
    #genome,scaffold,orf,coord,length,keggAcc,kegg,pfam = line.split('\t')
    scaffold,orf,coord,length,keggAcc,kegg,pfam = line.split('\t')g

    if re.search(r'trna',orf):
        print(orf)
        continue
    
    #scaffold = '_'.join( orf.split('_')[:-1] )
    print(coord)
    start,endStrand = coord.split('..')
    end,strand = endStrand.split()
    print(strand)

    
    # if not ( lineage == '-' or lineage == 'Elusimicrobia' ) :
    #     continue

    # if strand == '(u)' :
    #     strand = '1'
    # else:
    #     strand = '-1'


    if scaffold not in scaffold2min :
        if int(start) < int(end) :
            scaffold2min[scaffold] = int(start)
            scaffold2max[scaffold] = int(end)
        else:
            scaffold2min[scaffold] = int(end)
            scaffold2max[scaffold] = int(start)
            
    if int(start) < scaffold2min[scaffold] :
        scaffold2min[scaffold] = int(start)

    if int(end) < scaffold2min[scaffold] :
        scaffold2min[scaffold] = int(end)


        
    if int(end) > scaffold2max[scaffold] :
        scaffold2max[scaffold] = int(end)

    if int(start) > scaffold2max[scaffold] :
        scaffold2max[scaffold] = int(start)


        
#    if kegg == ' nitrogenase molybdenum-iron protein alpha chain' :
    if pfam == 'YcaO (PF02624)':
        orf2color[orf] = 'green'

#    if kegg == ' nitrogenase molybdenum-iron protein beta chain' :
    if pfam == 'Nitroreductase (PF00881)' :
        orf2color[orf] = 'red'

    if re.search(r'PF04055',pfam) :
        if pfam == 'Radical_SAM (PF04055)' :
            orf2color[orf] = 'blue'
        else:
            orf2color[orf] = 'purple'
            
    if pfam == 'Fer4_NifH (PF00142)' :
        orf2color[orf] = 'gold'

    if pfam == 'Fer4_NifH (PF00142) + Oxidored_nitro (PF00148)' :
        orf2color[orf] = 'darkorange'

    if pfam == 'Na' :
        orf2color[orf] = 'silver'

        
    #scaffold2features[scaffold].append( SeqFeature( FeatureLocation( int( int(start) ) , int( int(end) ) ) , strand=int(strand), type = 'CDS' ) )
    #print(scaffold+'\t'+start+'\t'+end+'\t'+strand+'\t'+kegg)

    if scaffold not in scaffold2orf2coord :
        scaffold2orf2coord[ scaffold ] = dict()

    scaffold2orf2coord[ scaffold ][orf] = (int(start),int(end),int(strand))
    print(scaffold2orf2coord[ scaffold ][orf])
file.close()



from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO


for scaffold,orf2coord in scaffold2orf2coord.items() :
    print(scaffold)
    scaffold_filename = scaffold+'_genomicDiagram.pdf'
    gd_diagram = GenomeDiagram.Diagram("a new type of nitrogenase")
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()

    
    # print('\t'+str(scaffold2max[scaffold]))
    # print('\t'+str(scaffold2min[scaffold]))

    for orf,coord in orf2coord.items() :
        start = coord[0]
        end = coord[1]
        strand = coord[2]
        feature = SeqFeature( FeatureLocation( int( int(start) ) , int( int(end) ) ) , strand=int(strand), type = 'CDS' )
        if orf in orf2color :
            color = orf2color[orf]
        else:
            color = 'white'
        gd_feature_set.add_feature(feature, color=color, name="toto", label=False,sigil="BIGARROW",arrowhead_length=0.25,arrowshaft_height=0.5)

    intervalle = int ( ( scaffold2max[scaffold] - scaffold2min[scaffold] + 1 ) * 0.025 )
    size = int( scaffold2max[scaffold] + intervalle ) - int( scaffold2min[scaffold] - intervalle )
    print('\t'+scaffold+'\t'+str(size) )

    if float(size / 1000) > 1 :
        l = float(size / 1000) * 1.5
    else:
        print(scaffold)
        l = float(size / 100) * 1.5
        
    print( '\t'+str(l) )
    if int( scaffold2min[scaffold] - intervalle ) < 1 :
        start = 1
    else:
        start = int( scaffold2min[scaffold] - intervalle )
        
    gd_diagram.draw(format="linear", orientation="landscape", pagesize=(3*cm,l*cm),fragments=1, start=start , end=int( scaffold2max[scaffold] + intervalle ) )
    gd_diagram.write(scaffold_filename, "PDF")


