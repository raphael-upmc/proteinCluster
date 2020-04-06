#! /usr/bin/env python

import os,sys,re
from collections import defaultdict

# make network object and load file
from clustergrammer import Network

if __name__ == "__main__":

    matrix_filename = sys.argv[1]
    html_output_filename = sys.argv[2]
    
    print('loading file...')
    net = Network()
    # load matrix file
    net.load_file(matrix_filename)
    print('done')

    
    # cluster using default parameters
    print('clustering the matrix...')
    net.cluster(dist_type='jaccard',linkage_type='complete')
#    net.cluster(run_clustering=False) 
    print('done')
    
    # save visualization JSON to file for use by front end
    print('saving results in json file...')
    json_filename = matrix_filename+'.json'
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

    load_viz_new_filename = '/home/meheurap/scripts/proteinCluster/load_viz_new.js'
    load_viz_new = ''
    file = open(load_viz_new_filename,'rt')
    for line in file :
        load_viz_new += line
    file.close()


    clustergrammer_js_filename = '/home/meheurap/scripts/proteinCluster/clustergrammer.js'
    cg_new = ''
    file = open(clustergrammer_js_filename,'rt')
    for line in file :
        cg_new += line
    file.close()

    
    print(html_output_filename)
    output = open(html_output_filename,'w')

    html_filename = '/home/meheurap/scripts/proteinCluster/clustergrammer_model.html'
    file = open(html_filename,'rt')
    for line in file :
        line = line.rstrip()
        if re.search('network_data = JSON',line) :
            output.write('\n\n\n'+cg_new+'\n\n\n')
            output.write('\n\n\n'+load_viz_new+'\n\n\n')
            output.write('network_data = '+network_data+'\n')
        else :
            output.write(line+'\n')
    file.close()
    output.close()
    print('done')


    
    # # creating the html page
    # print('creating the html page...')
    # network_data = ''
    # file = open(json_filename,'rt')
    # for line in file :
    #     network_data += line
    # file.close()
    # print(len(network_data))

    # load_viz_new_filename = '/home/meheurap/script/load_viz_new.js'
    # load_viz_new = ''
    # file = open(load_viz_new_filename,'rt')
    # for line in file :
    #     load_viz_new += line
    # file.close()

    

    # print(html_output_filename)
    # output = open(html_output_filename,'w')

    # html_filename = '/home/meheurap/script/clustergrammer_model.html'
    # file = open(html_filename,'rt')
    # for line in file :
    #     line = line.rstrip()
    #     if re.search('network_data = JSON',line) :
    #         output.write('\n\n\n'+load_viz_new+'\n\n\n')
    #         output.write('network_data = '+network_data+'\n')
    #     else :
    #         output.write(line+'\n')
    # file.close()
    # output.close()
    # print('done')
    

    sys.exit()
