#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor,wait

def HMMsearch() :
    hmm_database = ''
    directory = ''
    for directory in directoryList :
        for root, dirs, files in os.walk(directory):
            for filename in files :
                basename = filename.replace('.hmm','')
                hmm_filename = root+'/'+filename
                domtblout_filename = '/data7/proteinfams/genomicContext/domtblout/'+basename+'.domtblout'
                if not os.path.exists(domtblout_filename) :
                    cmd = 'hmmsearch -E 1e-10 --cpu 1 --domtblout '+domtblout_filename+' '+hmm_filename+' '+fasta_filename+' >/dev/null 2>/dev/null'
                    status = os.system(cmd)
                    print(str(status)+'\t'+cmd)

def readingHMM(domtblout_filename) :
    orf2hmm = defaultdict(list)
    file = open(domtblout_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.search(r'^#',line) :
            continue

        liste = line.split()
        orf = liste[0]
        length = liste[2]
        hmm = liste[3]
        ko = liste[3].split(".")[0]
        hmmLength = liste[5]
        evalue = liste[6]
        bitscore = liste[7]
        cEvalue = liste[11] # conditional Evalue
        iEvalue = liste[12] # independant Evalue

        hmmS = liste[15]
        hmmE = liste[16]

        aliS = liste[17]
        aliE = liste[18]

        envS = liste[19]
        envE = liste[20]

        orfCover = float( int(envE) - int(envS) + 1 ) / float(length)
        hmmCover = float( int(hmmE) - int(hmmS) + 1 ) / float(hmmLength)
        orf2hmm[ orf ].append( hmm , float(evalue) , float(bitscore) , orfCover , hmmCover , length, hmmLength )
    file.close()
    return orf2hmm


def running16RP(genome,cpt,cwd,seqList):
    fasta_filename = cwd+'/'+genome+'.faa'
    result_filename = cwd+'/'+genome+'.tsv'
    SeqIO.write(seqList,fasta_filename,'fasta')
    cmd = "/home/meheurap/.pyenv/shims/rp16.py -f "+fasta_filename+" -d /shared/db/ribosomal_proteins/hug/r04_11_2016 -t 1 "+"1>"+result_filename+" 2>>"+cwd+'/rp16.error'
    print(cmd)
    status = os.system(cmd)
    print(str(cpt)+'\t'+genome+'\t'+str(status))

    positionList = list()
    if status == 0 :
        file = open(result_filename,'r')
        next(file)
        for line in file :
            positionList.append(genome+'\t'+line)
        file.close()
            
    for usearch_filename in [genome+'-usearch_prot-rpL14_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpS8_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpL16_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpL4_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpS3_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpL5_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpS19_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpL24_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpL22_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpL18_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpL15_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpL3_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpS10_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpS17_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpL6_JGI_MDM.filtered.b6',genome+'-usearch_prot-rpL2_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL14_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpS8_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL16_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL4_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpS3_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL5_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpS19_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL24_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL22_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL18_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL15_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL3_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpS10_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpS17_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL6_JGI_MDM.filtered.b6',genome+'-usearch_nucl-rpL2_JGI_MDM.filtered.b6'] :

        if os.path.exists(usearch_filename) :
            os.remove(usearch_filename)

    if os.path.exists(fasta_filename) :
        os.remove(fasta_filename)
        
    if os.path.exists(result_filename) :
        os.remove(result_filename)
        
    return(status,genome,positionList)


def running16RP_low_memory(genome,cpt,cwd,fasta_filename):
    result_filename = cwd+'/'+str(cpt)+'.tsv'
    cmd = "/home/meheurap/.pyenv/shims/rp16.py -f "+fasta_filename+" -d /shared/db/ribosomal_proteins/hug/r04_11_2016 -t 1 "+"1>"+result_filename+" 2>/dev/null"
    status = os.system(cmd)
    #print(str(cpt)+'\t'+genome+'\t'+str(status))

    positionList = list()
    if status == 0 :
        file = open(result_filename,'r')
        next(file)
        for line in file :
            positionList.append(genome+'\t'+line)
        file.close()
            
    for usearch_filename in [str(cpt)+'-usearch_prot-rpL14_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpS8_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpL16_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpL4_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpS3_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpL5_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpS19_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpL24_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpL22_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpL18_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpL15_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpL3_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpS10_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpS17_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpL6_JGI_MDM.filtered.b6',str(cpt)+'-usearch_prot-rpL2_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL14_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpS8_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL16_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL4_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpS3_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL5_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpS19_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL24_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL22_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL18_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL15_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL3_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpS10_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpS17_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL6_JGI_MDM.filtered.b6',str(cpt)+'-usearch_nucl-rpL2_JGI_MDM.filtered.b6'] :

        if os.path.exists(usearch_filename) :
            os.remove(usearch_filename)

    if os.path.exists(fasta_filename) :
        os.remove(fasta_filename)
        
    if os.path.exists(result_filename) :
        os.remove(result_filename)
        
    return(status,genome,positionList)


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='extracting the 16 ribosomal proteins')
    parser.add_argument('protein_filename', help='the path of the FASTA_PROTEIN_FILE')
    parser.add_argument('orf2bin_filename',help='the path of the ORF2BIN_FILE, this file is a tab-separated file, first col is the orf, second col is the genome. First line is skipped')
    parser.add_argument('output_filename',help='the path of the OUTPUT_FILE, results will be stored in this file')
    parser.add_argument('output_summary_filename',help='the path of the OUTPUT_SUMMARY_FILE, results will be stored in this file')
    parser.add_argument('--cpu',type=int,default=6,help='number of CPUs (default: 6)')
    parser.add_argument('--rp14',action='store_true',default=False,help='only consider 14RPs for Archaea (default: False)')
    parser.add_argument('--low-memory',action='store_true',default=False,help='low-memory (longer) (default: False)')
    
    args = parser.parse_args()
    
    if os.path.exists(args.protein_filename) :
        protein_filename = os.path.abspath(args.protein_filename)
    else:
        sys.exit(args.protein_filename+' does not exist, exit')

    if os.path.exists(args.orf2bin_filename) :
        orf2bin_filename = os.path.abspath(args.orf2bin_filename)
    else:
        sys.exit(args.orf2bin_filename+' does not exist, exit')

    cpu = args.cpu
    
    output_summary_filename = os.path.abspath(args.output_summary_filename)
    
    output_aln_filename = os.path.abspath(args.output_filename)
    cwd = '/'.join(output_aln_filename.split('/')[:-1])
    
    print('protein_filename: '+protein_filename)
    print('orf2bin_filename: '+orf2bin_filename)
    print('output_aln_filename: '+output_aln_filename)
    print('output_summary_filename: '+output_summary_filename)
    print('number of CPUs: '+str(cpu))
    print('current working directory: '+cwd)
    print('14RP: '+str(args.rp14))
    print('low-memory: '+str(args.low_memory))
    
    folder = cwd+"/16RP_results"
    if os.path.exists(folder) :
        sys.exit(folder+" already exists, remove it first")

    os.mkdir(folder)


    
    if args.low_memory :
        genomeSet = set()
        file = open(orf2bin_filename,'r')
        header = next(file)
        for line in file :
            line = line.rstrip()
            orf,genome = line.split('\t')
            genomeSet.add(genome)
        file.close()
        print(len(genomeSet))
        
        bin2genomeSet = defaultdict(set)
        nb = 4500
        print(nb)
        cpt_bin = 0
        cpt_genome = 0
        for genome in genomeSet :
            bin2genomeSet[ cpt_bin ].add(genome)
            if cpt_genome > nb :
                print(str(cpt_bin)+'\t'+str(len(bin2genomeSet[ cpt_bin ])))
                cpt_genome = 0
                cpt_bin += 1
            else:
                cpt_genome += 1
        
        print(str(cpt_bin)+'\t'+str(len(bin2genomeSet[ cpt_bin ])))
        print('number of bins: '+str(len(bin2genomeSet)))
        print('\n\n')

        print('starting 16rp...')
        bin2results = defaultdict(list)
        for cpt_bin in sorted(bin2genomeSet) :
            pool = ProcessPoolExecutor(cpu) # start 20 worker processes and 1 maxtasksperchild in order to release memory    
            cpt = 0

            orf2bin = dict()
            file = open(orf2bin_filename,'r')
            header = next(file)
            for line in file :
                line = line.rstrip()
                orf,genome = line.split('\t')
                if genome in bin2genomeSet[ cpt_bin ] :
                    orf2bin[ orf ] = genome
                else:
                    continue
            file.close()

            
            print(str(cpt_bin)+'\t'+str(len(bin2genomeSet[cpt_bin])))

            genome2seqList = defaultdict(list)
            for record in SeqIO.parse(protein_filename,'fasta') :                
                if record.id not in orf2bin :
                    continue
                else:
                    genome = orf2bin[ record.id ]
                    if genome in bin2genomeSet[cpt_bin] :
                        genome2seqList[genome].append(record)
                    else:
                        sys.exit('error')

            
            for genome in bin2genomeSet[cpt_bin] :
                cpt += 1
                print('\t'+str(cpt)+'\t'+genome)
                
                genome_filename = folder+'/'+str(cpt)+'.faa'
                SeqIO.write(genome2seqList[genome],genome_filename,'fasta')
                del[ genome2seqList[genome][:] ]

                future = pool.submit( running16RP_low_memory,genome,cpt,cwd,genome_filename )
                bin2results[cpt_bin].append(future)
            wait(bin2results[cpt_bin])
            
            print('\t'+str(len(bin2results[cpt_bin])))
            orf2bin.clear()
            genome2seqList.clear()            
            pool.shutdown()
            
        rp16_table_filename = cwd+'/'+'rp16_table.tsv'
        output = open(rp16_table_filename,'w')
        output.write('genome\tscaffold\tL15\tL18\tL6\tS8\tL5\tL24\tL14\tS17\tL16\tS3\tL22\tS19\tL2\tL4\tL3\tS10'+'\n')

        print('\n\n')
        orfSet = set()
        error = 0
        for cpt_bin,results in bin2results.items() :
            for elt in results :
                status,genome,positionResult = elt.result()
                if status != 0 :
                    print('\t'+genome+' '+'==>'+' '+'Error' )
                    error += 1
                else:
                    for line in positionResult :
                        line = line.rstrip()
                        output.write(line+'\n')
                        liste = line.split('\t')
                        del(liste[0])
                        for orf in liste :
                            if orf == '-' :
                                continue
                            else:
                                orfSet.add(orf)
        output.close()
        print(len(orfSet))

        
        if error != 0 :
            print('\n'+str(error)+' genomes failed to run 16RP.py\n')
        print('done')            
        
        orf2bin = dict()
        file = open(orf2bin_filename,'r')
        header = next(file)
        for line in file :
            line = line.rstrip()
            orf,genome = line.split('\t')
            if orf in orfSet :
                orf2bin[ orf ] = genome
        file.close()
        print(len(orf2bin))
        genome2seqList = defaultdict(list)
        orf2seq = dict()
        for record in SeqIO.parse(protein_filename,'fasta') :                
            if record.id in orfSet :
                genome = orf2bin[ record.id ]
                genome2seqList[ genome ].append( record )
                orf2seq[record.id] = record
        print(len(orf2seq))

    else:
        genomeSet = set()
        orf2bin = dict()
        file = open(orf2bin_filename,'r')
        header = next(file)
        for line in file :
            line = line.rstrip()
            orf,genome = line.split('\t')
            orf2bin[ orf ] = genome
            genomeSet.add(genome)
        file.close()
        print(len(genomeSet))

        genome2seqList = defaultdict(list)
        orf2seq = dict()
        for record in SeqIO.parse(protein_filename,'fasta') :
            genome = orf2bin[ record.id ]
            genome2seqList[ genome ].append( record )
            orf2seq[record.id] = record
        print(str(len(genome2seqList))+' genomes')


        #################################################
        # running the script of ctb to detect the 16RPs #
        #################################################

        error = 0
        results = list()
        pool = ProcessPoolExecutor(cpu) # start 20 worker processes and 1 maxtasksperchild in order to release memory    
        cpt = 0
        for genome,seqList in genome2seqList.items() :
            cpt += 1
            future = pool.submit( running16RP,genome,cpt,cwd,seqList )
            results.append(future)
            # if cpt == 100 :
            #     break
        wait(results)


        rp16_table_filename = cwd+'/'+'rp16_table.tsv'
        output = open(rp16_table_filename,'w')
        output.write('genome\tscaffold\tL15\tL18\tL6\tS8\tL5\tL24\tL14\tS17\tL16\tS3\tL22\tS19\tL2\tL4\tL3\tS10'+'\n')

        print('\n\n')
        for elt in results :
            status,genome,positionResult = elt.result()
            if status != 0 :
                print('\t'+genome+' '+'==>'+' '+'Error' )
                error += 1
            else:
                for line in positionResult :
                    output.write(line)
        output.close()
        pool.shutdown()
    
        if error != 0 :
            print('\n'+str(error)+' genomes failed to run 16RP.py\n')
        print('done')

    
    ###########################
    # reading the 16RP output #
    ###########################
    
    genome2scaffold2rp = dict()
    file = open(rp16_table_filename,'r')
    header = next(file).rstrip().split('\t')
    for line in file :
        line = line.rstrip()
        liste = line.split('\t')
        genome = liste[0]
        scaffold = liste[1]

        if genome not in genome2scaffold2rp :
            genome2scaffold2rp[ genome ] = defaultdict(dict)

        if scaffold in genome2scaffold2rp[ genome ] :
            sys.exit('error! should be impossible...')
        
        genome2scaffold2rp[ genome ][scaffold] = dict()

        for i in range(2,len(liste)) :
            rp = header[i]

            if args.rp14 :
                if rp == 'L16' or rp == 'S10' :
                    continue
       
            if liste[i] == '-' :
                continue
            else :
                orf = liste[i]
                genome2scaffold2rp[ genome ][scaffold][ rp ] = orf
    file.close()
    print(len(genome2scaffold2rp))


    ###############################
    # selecting the best scaffold #
    ###############################
    
    genome2summary = dict()
    genome2scaffold = dict()
    for genome,scaffold2rp in genome2scaffold2rp.items() :
        nb = 0
        contaminationSet = set()
        contamination = 'no'
        rpSet = set()
        for scaffold in scaffold2rp :
            nb += len(scaffold2rp[scaffold])
            for rp in scaffold2rp[scaffold] :
                if rp not in rpSet :
                    rpSet.add(rp)
                else :
                    contaminationSet.add(rp)
                    contamination = 'yes' # the same RP is present in several copy

        scaffold_nb = len(scaffold2rp)
        scaffold = sorted(scaffold2rp.items(),key=lambda x:len(x[1]), reverse = True)[0][0] #selecting the scaffold with the most of 16RP along
        best = len( scaffold2rp[scaffold] )
        genome2scaffold[genome] = scaffold
        if len(contaminationSet) == 0 :
            result = '-'
        else:
            result = ','.join(list(contaminationSet))
        genome2summary[ genome ] = genome+'\t'+str(scaffold_nb)+'\t'+str(nb)+'\t'+str(best)+'\t'+contamination+'\t'+result


    ##########################################################
    # extracting the fasta sequences and performing the MSAs #
    ##########################################################
    
    rp2seq = defaultdict(list)
    for genome,scaffold2rp in genome2scaffold2rp.items() :
        if genome not in genome2scaffold :
            continue
        scaffold = genome2scaffold[genome]
        
        for rp,seq in scaffold2rp[scaffold].items() :
            rp2seq[rp].append( SeqRecord(seq=orf2seq[seq].seq,id=genome,description="") )


    print('performing MSA...')
    for rp,seqList in rp2seq.items() :
        print(rp+'\t'+str( len(seqList) ) )

        output_filename = folder+'/'+rp+'.fa'
        SeqIO.write(seqList,output_filename,'fasta')
        mafft_filename = output_filename.replace('.fa','.mafft')
        cmd = '/shared/software/bin/mafft --auto --thread '+str(cpu)+' '+output_filename+' > '+mafft_filename+' 2>/dev/null'
        print(cmd)
        os.system(cmd)
        
        trimal_filename = mafft_filename.replace('.mafft','.trimal')
        cmd = '/groups/banfield/users/meheurap/programs/trimal-trimAl/source/trimal -fasta -gappyout -in '+mafft_filename+' -out '+trimal_filename
        print(cmd)
        os.system(cmd)
    
    print('done')

    print('creating final output...')
    genome2aln = defaultdict(str)
    for (path, dirs, files) in os.walk(folder):
        for filename in files :
            if not re.search(r'.trimal$',filename) :
                continue
            print(filename)
            rpGenomeSet = set()
            lengthSet = set()
            for seq_record in SeqIO.parse(path+'/'+filename, "fasta"):
                l = len(seq_record)
                lengthSet.add(l)
                genome2aln[ seq_record.description ] += str(seq_record.seq)
                rpGenomeSet.add(seq_record.description)
            if len(lengthSet) != 1 :
                sys.exit('error : '+str(lengthSet) )

            fakeSeq = ''
            for i in range(l) :
                fakeSeq += '-'

            for genome in genome2scaffold2rp :
                if genome not in rpGenomeSet :
                    genome2aln[ genome ] += fakeSeq
                else :
                    continue

                
    ###############################
    # Concatenating the 16RP MSAs #
    ###############################

    print('genome2aln: '+str(len(genome2aln)))
    output1 = open(output_summary_filename,'w')
    output1.write('genome'+'\t'+'nb_of_scaffolds'+'\t'+'nb_of_RPs'+'\t'+'nb_of_RPs_on_the_best_scaffold'+'\t'+'are_RPs_duplicated'+'\t'+'list_of_RPs_duplicated'+'\t'+'size'+'\n')

    lengthSet = set()
    output = open(output_aln_filename,'w')
    for genome,aln in genome2aln.items() :

        l = float(len(aln.replace('-',''))) / float(len(aln))
        genome2summary[ genome ] += '\t'+str(l)
        output1.write(genome2summary[ genome ]+'\n')

        output.write('>'+genome+'\n')
        output.write(aln+'\n')
    output.close()
    output1.close()
    print('done')




