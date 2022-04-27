#! /usr/bin/env python

import os,sys,re
from collections import defaultdict
from Bio import SeqIO
import argparse
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ProcessPoolExecutor,wait

def kNearestNeighbors(markerSet,geneList,k) :
    last = None
    cluster2markers = dict()
    cluster = 0
    sortedList = sorted(geneList.items(),key=lambda x:x[1][0])
    for i in range( len(sortedList) ) :
        if sortedList[i][0] not in markerSet :
            continue
        else:
            if cluster not in cluster2markers :
                last = i
                cluster2markers[cluster] = [ sortedList[i][0] ]
            else:
                delta = i - last - 1
                if delta <= k : 
                    cluster2markers[cluster].append(sortedList[i][0])
                    last = i
                else:
                    cluster += 1
                    cluster2markers[cluster] = [ sortedList[i][0] ]
                    last = i
    return cluster2markers


def rp2pfam(rp_nb) :
    if rp_nb :
        rp2pfam = {'RPL14' : 'PF00238' , 'RPL15' : 'PF00828' , 'RPL18' : 'PF00861' , 'RPL22' : 'PF00237' , 'RPL24' : 'PF17136' , 'RPL2' : 'PF03947' , 'RPL3' : 'PF00297' , 'RPL4' : 'PF00573' , 'RPL5' : 'PF00673' , 'RPL6': 'PF00347' , 'RPS17' : 'PF00366' , 'RPS19' : 'PF00203' , 'RPS3' : 'PF00189' , 'RPS8' : 'PF00410' }
    else:
        rp2pfam = {'RPL14' : 'PF00238' , 'RPL15' : 'PF00828' , 'RPL16' : 'PF00252' , 'RPL18' : 'PF00861' , 'RPL22' : 'PF00237' , 'RPL24' : 'PF17136' , 'RPL2' : 'PF03947' , 'RPL3' : 'PF00297' , 'RPL4' : 'PF00573' , 'RPL5' : 'PF00673' , 'RPL6': 'PF00347' , 'RPS10' : 'PF00338' , 'RPS17' : 'PF00366' , 'RPS19' : 'PF00203' , 'RPS3' : 'PF00189' , 'RPS8' : 'PF00410' }

    pfam2rp = dict()
    for rp,pfam in rp2pfam.items() :
        pfam2rp[pfam] = rp

    return rp2pfam, pfam2rp

def buildingHmmDb(pfam2rp,hmm_filename) :
    output = open(hmm_filename,'w')

    pfamList = list()
    name2accession = dict()
    pfam2desc = dict()
    pfamList_filename = '/env/cns/db/Pfam/Pfam_latest/Pfam-A.hmm'
    file = open(pfamList_filename,'r')
    for line in file :
        line = line.rstrip()
        pfamList.append(line)

        liste = line.split()
        if re.match(r'NAME',line) :
            name = ' '.join(liste[1:])
            accession = ''
            desc = ''

        if re.match(r'ACC',line) :
            accession = ' '.join(liste[1:])
            name2accession[ name ] = accession            
            
        if re.match(r'DESC',line) :
            desc = ' '.join(liste[1:])
            pfam2desc[accession] = desc+' ('+accession+')'

        if line == '//' :
            if accession.split('.')[0] in pfam2rp :
                print(pfam2desc[accession])
                output.write('\n'.join(pfamList)+'\n')
                pfamList = []
            else:
                pfamList = []
                continue
                print(pfam2desc[accession])
    file.close()
    output.close()


def runningHMM(domtblout_filename,hmm_filename,fasta_filename,cpu) :
    cmd = 'hmmsearch -E 1e-3 --cpu '+str(cpu)+' --domtblout '+domtblout_filename+' '+hmm_filename+' '+fasta_filename+' >/dev/null 2>/dev/null' 
    status = os.system(cmd)
    return cmd,status

def readingHMM(domtblout_filename) :
    orf2hmm = dict()
    file = open(domtblout_filename,'r')
    for line in file :
        line = line.rstrip()
        if re.search(r'^#',line) :
            continue

        liste = line.split()
        orf = liste[0]
        length = liste[2]
        hmm = liste[3]
        pfam = liste[4].split(".")[0]
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
        orf2hmm[ orf ] = [ pfam , float(evalue) , float(bitscore) , orfCover , hmmCover , length, hmmLength ]
    file.close()
    return orf2hmm


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='extracting the 16 ribosomal proteins')
    parser.add_argument('protein_filename', help='the path of the FASTA_PROTEIN_FILE')
    parser.add_argument('feature_filename',help='the path of the FEATURE_FILE, this file is a tab-separated file.')
    parser.add_argument('--cpu',type=int,default=1,help='number of CPUs (default: 1)')
    parser.add_argument('--rp14',action='store_true',default=False,help='only consider 14RPs for Archaea (default: False)')
    
    args = parser.parse_args()
    
    if os.path.exists(args.protein_filename) :
        protein_filename = os.path.abspath(args.protein_filename)
    else:
        sys.exit(args.protein_filename+' does not exist, exit')

    if os.path.exists(args.feature_filename) :
        feature_filename = os.path.abspath(args.feature_filename)
    else:
        sys.exit(args.feature_filename+' does not exist, exit')


    cpu = args.cpu
    cwd = os.getcwd()
    
    print('protein_filename: '+protein_filename)
    print('feature_filename: '+feature_filename)
    print('number of CPUs: '+str(cpu))
    print('current working directory: '+cwd)
    print('14RP: '+str(args.rp14))
    
    folder = cwd+"/16RP_results"
    if os.path.exists(folder) :
        sys.exit(folder+" already exists, remove it first")

    os.mkdir(folder)


    output_aln_filename = folder+'/'+'16RP.aln'
    output_summary_filename = folder+'/'+'16RP.summary'
    matrix_filename = folder+'/'+'16RP.matrix'
    missing_genomes_filename = folder+'/'+'16RP.missingGenomes'

    rp2pfam, pfam2rp = rp2pfam(args.rp14)
    hmm_filename = folder+'/'+'16RP.hmm'
    if not os.path.exists(hmm_filename) :
        buildingHmmDb(pfam2rp,hmm_filename)
    
    domtblout_filename = folder+'/'+'16RP.domtblout'
    if not os.path.exists(domtblout_filename) :
        cmd,status = runningHMM(domtblout_filename,hmm_filename,protein_filename,cpu)
        print(cmd)
        print(status)
    
 
    orf2hmm = readingHMM(domtblout_filename)
    
    # get fasta sequences of ORFs matching a genetic markers
    orf2seq = dict()
    for record in SeqIO.parse(protein_filename,'fasta') :
        if record.id in orf2hmm :
            orf2seq[record.id] = record
        else:
            continue

    # get the genome and scaffold names of genetic markers
    totalGenomeSet = set()
    scaffoldSet = set()
    file = open(feature_filename,'r')
    for line in file :
        line = line.rstrip()
        orf,genome,scaffold,start,end,strand = line.split('\t')
        totalGenomeSet.add(genome)
        if orf in orf2hmm :
            scaffoldSet.add(genome+'\t'+scaffold)        
    file.close()


    # get the genomic context of scaffolds encoding a genetic markers
    genome2scaffold2orf = dict()
    orf2coordinate = dict()
    file = open(feature_filename,'r')
    for line in file :
        line = line.rstrip()
        orf,genome,scaffold,start,end,strand = line.split('\t')
        if genome+'\t'+scaffold in scaffoldSet :
            orf2coordinate[ orf ] = [ int(start),int(end),strand ]
            if genome not in genome2scaffold2orf :
                genome2scaffold2orf[ genome ] = defaultdict(list)
            
            if scaffold not in genome2scaffold2orf[ genome ] :
                genome2scaffold2orf[ genome ][scaffold] = dict()
            genome2scaffold2orf[ genome ][scaffold][ orf ] = [ int(start) , int(end) , strand ]
    file.close()

    # partioning the markers into gene clusters using the parameter k
    genome2scaffold2cluster2orf = dict()
    for genome,scaffold2orf in genome2scaffold2orf.items() :
        print(genome)
        if genome not in genome2scaffold2cluster2orf :
            genome2scaffold2cluster2orf[genome] = dict()
        for scaffold,orf2coordinate in scaffold2orf.items() :
            print('\t'+scaffold)
            cluster2markers = kNearestNeighbors(orf2hmm,orf2coordinate,2)
            genome2scaffold2cluster2orf[ genome ][scaffold] = cluster2markers


    # writing the matrix output
    output = open(matrix_filename,'w')
    output.write('genome'+'\t'+'scaffold'+'\t'+'cluster')
    for rp,pfam in sorted(rp2pfam.items()) :
        output.write('\t'+rp+' ('+pfam+')')
    output.write('\n')

    for genome, scaffold2cluster2orf in genome2scaffold2cluster2orf.items() :
        print(genome)
        for scaffold,cluster2markers in scaffold2cluster2orf.items() :
            for cluster,liste in cluster2markers.items() :
                output.write(genome+'\t'+scaffold+'\t'+str(cluster))
                rp2orfs = defaultdict(set)
                for orf in liste :
                    pfam = orf2hmm[orf][0]
                    rp = pfam2rp[pfam]
                    rp2orfs[rp].add(orf)

                print('\t\t'+str(cluster))
                for rp,pfam in sorted(rp2pfam.items()) :
                    if rp in rp2orfs :
                        output.write('\t'+','.join(list(rp2orfs[rp])))
                    else:
                        output.write('\t'+'-')
                output.write('\n')
    output.close()


    ###############################
    # selecting the best scaffold #
    ###############################
    
    genome2summary = dict()
    genome2scaffold = dict()

    for genome,scaffold2cluster2orf in genome2scaffold2cluster2orf.items() :
        cluster_nb = 0
        scaffold_nb = len(scaffold2cluster2orf)
        contaminationSet = set()
        contamination = 'no'
        rpSet = set()
        nb = 0
        best = 0        
        for scaffold in scaffold2cluster2orf :
            cluster_nb += len(scaffold2cluster2orf[scaffold])
            for cluster in scaffold2cluster2orf[scaffold] :
                cluster_contamination = 'no'
                clusterRpSet = set()

                nb += len(scaffold2cluster2orf[scaffold][cluster])


                for orf in scaffold2cluster2orf[scaffold][cluster] : 
                    pfam = orf2hmm[orf][0]
                    rp = pfam2rp[pfam]

                    if rp not in rpSet :
                        rpSet.add(rp)
                    else :
                        contaminationSet.add(rp)
                        contamination = 'yes' # the same RP is present in several copy

                    if rp not in clusterRpSet :
                        clusterRpSet.add(rp)
                    else :
                        cluster_contamination = 'yes' # the same RP is present in several copy in the same cluster


                if len(scaffold2cluster2orf[scaffold][cluster]) > best and cluster_contamination == 'no':
                    best = len(scaffold2cluster2orf[scaffold][cluster])
                    scaffoldBest = [scaffold,cluster]
                    
        genome2scaffold[genome] = scaffoldBest 

        if len(contaminationSet) == 0 :
            result = '-'
        else:
            result = ','.join(list(contaminationSet))

        genome2summary[ genome ] = genome+'\t'+str(scaffold_nb)+'\t'+str(cluster_nb)+'\t'+str(nb)+'\t'+str(best)+'\t'+contamination+'\t'+result


    ##########################################################
    # extracting the fasta sequences and performing the MSAs #
    ##########################################################
    
    rp2seq = defaultdict(list)
    for genome in genome2scaffold :
        scaffold,cluster = genome2scaffold[genome]
        for orf in genome2scaffold2cluster2orf[genome][scaffold][cluster] : 
            pfam = orf2hmm[orf][0]
            rp = pfam2rp[pfam]        
            rp2seq[rp].append( SeqRecord(seq=orf2seq[orf].seq,id=genome,description="") )


    print('performing MSA...')
    for rp,seqList in rp2seq.items() :
        print(rp+'\t'+str( len(seqList) ) )

        output_filename = folder+'/'+rp+'.fa'
        SeqIO.write(seqList,output_filename,'fasta')
        mafft_filename = output_filename.replace('.fa','.mafft')
        cmd = 'mafft --auto --thread '+str(cpu)+' '+output_filename+' > '+mafft_filename+' 2>/dev/null'
        print(cmd)
        os.system(cmd)
        
        trimal_filename = mafft_filename.replace('.mafft','.trimal')
        cmd = 'trimal -fasta -gappyout -in '+mafft_filename+' -out '+trimal_filename
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

            for genome in genome2scaffold :
                if genome not in rpGenomeSet :
                    genome2aln[ genome ] += fakeSeq
                else :
                    continue

                
    ###############################
    # Concatenating the 16RP MSAs #
    ###############################

    print('genome2aln: '+str(len(genome2aln)))
    output1 = open(output_summary_filename,'w')
    output1.write('genome'+'\t'+'nb_of_scaffolds'+'\t'+'nb_of_clusters'+'\t'+'nb_of_RPs'+'\t'+'nb_of_RPs_on_the_best_scaffold'+'\t'+'are_RPs_duplicated'+'\t'+'list_of_RPs_duplicated'+'\t'+'size'+'\n')

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

    output = open(missing_genomes_filename,'w')
    for genome in totalGenomeSet :
        if genome not in genome2aln :
            output.write(genome+'\n')
        else:
            continue
    output.close()

