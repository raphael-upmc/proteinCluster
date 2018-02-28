#! /home/meheurap/.pyenv/shims/python

import os,sys,re
from collections import defaultdict

# make network object and load file
from clustergrammer import Network

if __name__ == "__main__":

    # genome annotation
    bin2taxonomy = dict()
    filename = sys.argv[3]
    file = open(filename,"r")
    next(file)
    for line in file :
        line = line.rstrip()
        genome,ggkbaseTaxonomy,normalizedTaxonomy = line.split("\t")
        taxonomy = normalizedTaxonomy.split(',')[1]
        phylum = normalizedTaxonomy.split(',')[0]
        bin2taxonomy[ genome ] = 'Taxonomy: '+taxonomy+'\t'+'Phylum: '+phylum
    file.close()


    # family enrichment
    family2enrichment = dict()
    family2enrichmentOddRatio = dict()    
    filename = '/home/meheurap/proteinCluster/coreCPR/enrichedFam/enrichFam.txt'
    file = open(filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        family = liste[0]
        enrichment = liste[1]
        family2enrichment[ family ] = enrichment        
        if enrichment == 'enriched' :
            oddratio = float(liste[4])
            if oddratio < 2 :
                enrichment = 'equally_distributed'
        elif enrichment == 'depleted' :
            oddratio = float(liste[7])
            if oddratio < 2 :
                enrichment = 'equally_distributed'

        family2enrichmentOddRatio[ family ] = enrichment        
    file.close()

    # family annotation      
    family2desc = dict()
    family2annot = defaultdict(list)
    filename = '/home/meheurap/script/family2annotation.txt'
    file = open(filename,'r')
    header = next(file)
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')

        family = liste[0]
        keggDesc = liste[9].split(' (')[0]
        family2desc[family] = keggDesc
        keggAnnotBig = liste[12].split(' (')[0]
        keggAnnot = liste[11].split(' (')[0]
        pathway = liste[10].split(' (')[0]
        psort = liste[5].split(' (')[0]
#        family2annot[family].append('module: '+module)
        family2annot[family].append('KEGG_big: '+keggAnnotBig)
        family2annot[family].append('psort: '+psort)
        family2annot[family].append('KEGG: '+keggAnnot)
        family2annot[family].append('pathway: '+pathway)

        if family not in family2enrichment :
            family2annot[family].append( 'enrichment: '+'Unknown' )
        else :
            family2annot[family].append( 'enrichment: '+family2enrichment[ family ] )

        if family not in family2enrichmentOddRatio :
            family2annot[family].append( 'enrichOddRatio: '+'Unknown' )
        else :
            family2annot[family].append( 'enrichOddRation: '+family2enrichmentOddRatio[ family ] )
            
            
    file.close()


    print('writting matrix...')
    matrix_filename = sys.argv[1]

    output_filename = os.path.basename(matrix_filename)+'.tmp'
    output = open(output_filename,'w')

    file = open(matrix_filename,'rt')

    familyList = next(file).rstrip().split('\t')
    del(familyList[0])
    header = '\t\t'
    for family in familyList :
        header += '\t'+'family: '+family+'|'+family2desc[family]
    output.write(header+'\n')
        
    header = '\t\t'
    for family in familyList :
        header += '\t'+family2annot[family][0]
    output.write(header+'\n')

    header = '\t\t'
    for family in familyList :
        header += '\t'+family2annot[family][1]
    output.write(header+'\n')

    header = '\t\t'
    for family in familyList :
        header += '\t'+family2annot[family][2]
    output.write(header+'\n')

    header = '\t\t'
    for family in familyList :
        header += '\t'+family2annot[family][3]
    output.write(header+'\n')

    header = '\t\t'
    for family in familyList :
        header += '\t'+family2annot[family][4]
    output.write(header+'\n')

    header = '\t\t'
    for family in familyList :
        header += '\t'+family2annot[family][5]
    output.write(header+'\n')
    
    # header = '\t\t'
    # for family in familyList :
    #     header += '\t'+family2annot[family][3]
    #     output.write(header+'\n')


    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        genome = liste[0]
        taxonomy = bin2taxonomy[ genome ]
        output.write('genome: '+genome+'\t'+taxonomy+'\t'+'\t'.join(liste[1:])+'\n')
    file.close()
    output.close()
    print('done')
    print(output_filename)


    print('loading file...')
    net = Network()
    # load matrix file
    net.load_file(output_filename)
    print('done')

    
    # cluster using default parameters
    print('clustering the matrix...')
    net.cluster(dist_type='jaccard',linkage_type='complete')
    print('done')
    
    # save visualization JSON to file for use by front end
    print('saving results in json file...')
    json_filename = output_filename+'.json'
    net.write_json_to_file('viz', json_filename)
    print('done')


    # creating the html page
    print('creating the html page...')
    network_data = ''
    file = open(json_filename,'rt')
    for line in file :
        network_data += line
    file.close()
    print(len(network_data))
    
    html_output_filename = sys.argv[2]
    print(html_output_filename)
    output = open(html_output_filename,'w')

    html_filename = '/home/meheurap/script/clustergrammer_model.html'
    file = open(html_filename,'rt')
    for line in file :
        line = line.rstrip()
        if re.search('network_data = JSON',line) :
            output.write('network_data = '+network_data+'\n')
        else :
            output.write(line+'\n')
    file.close()
    output.close()
    print('done')
    

    sys.exit()
