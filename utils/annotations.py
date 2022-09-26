# coding: utf-8

''' a list of useful functions for genome annotation '''


import os
import re
from collections import defaultdict
from operator import itemgetter
from Bio import SeqIO

PFAM_HMM_FILENAME = '/env/cns/proj/projet_CSX/scratch/bank/pfam/release_35.0/Pfam-A.hmm'
PFAM_INFO_FILENAME = '/env/cns/proj/projet_CSX/scratch/bank/pfam/release_35.0/Pfam-A.clans.tsv'
KEGG_HMM_FILENAME = '/env/cns/proj/projet_CSX/scratch/bank/kegg/release_100.0/KOFAM.hmm'
KEGG_INFO_FILENAME = '/env/cns/proj/projet_CSX/scratch/bank/kegg/release_100.0/ko_list'


def running_homemade_hmm(protein_filename,hmm_filename,domtblout_filename,cpu,genome) :
    cmd1 = 'hmmsearch -E 1e-3 --cpu ' \
           +str(cpu)+' --domtblout ' \
           +domtblout_filename+' ' \
           +hmm_filename+' '+protein_filename \
           +' >/dev/null 2>/dev/null'
    status1 = os.system(cmd1)
    besthit_filename = domtblout_filename+'.besthit'
    cmd2 = '/env/cns/proj/agc/home/rmeheust/scripts/' \
           'usefulScripts/cath-resolve-hits.ubuntu14.04' \
           '--input-format hmmer_domtblout '+ \
           domtblout_filename+ \
           ' >'+besthit_filename
    status2 = os.system(cmd2)
    return status1,cmd1,status2,cmd2,genome


def running_tmhmm(protein_filename,tmhmm_filename,genome):
    cmd = '/env/cns/proj/agc/home/rmeheust/programs/tmhmm-2.0c/bin/tmhmm' \
          '--short '+protein_filename+'  > '+tmhmm_filename
    status = os.system(cmd)
    return status,cmd,genome

def running_pfam(protein_filename,pfam_hmm_filename,domtblout_filename,cpu,genome):
    cmd1 = 'hmmsearch -E 1e-3 --cpu '+str(cpu) \
           +' --domtblout '+domtblout_filename \
           +' '+pfam_hmm_filename+' '+protein_filename \
           +' >/dev/null 2>/dev/null'
    status1 = os.system(cmd1)
    besthit_filename = domtblout_filename+'.besthit'
    cmd2 = '/env/cns/proj/agc/home/rmeheust/scripts/usefulScripts/' \
           'cath-resolve-hits.ubuntu14.04' \
           ' --input-format hmmer_domtblout ' \
           +domtblout_filename+' >'+besthit_filename
    status2 = os.system(cmd2)
    os.remove(domtblout_filename)
    return status1,cmd1,status2,cmd2,genome

def running_kegg(protein_filename,kegg_hmm_filename,domtblout_filename,genome):
    cmd1 = 'hmmsearch -E 1e-3 --cpu '+str(1) \
           +' --domtblout '+domtblout_filename \
           +' '+kegg_hmm_filename+' '+protein_filename \
           +' >/dev/null 2>/dev/null'
    status1 = os.system(cmd1)
    besthit_filename = domtblout_filename+'.besthit'
    cmd2 = '/env/cns/proj/agc/home/rmeheust/scripts/usefulScripts/' \
           'cath-resolve-hits.ubuntu14.04' \
           ' --input-format hmmer_domtblout ' \
           +domtblout_filename+' >'+besthit_filename
    status2 = os.system(cmd2)
    os.remove(domtblout_filename)
    return status1,cmd1,status2,cmd2,genome

def running_signalp6(protein_filename,signalp_filename,signalp_directory,genome) :
    log_filename = '/dev/null'
    cmd = 'source /env/cns/proj/agc/scratch/conda/miniconda.profile' \
          '&& conda activate signalp-6.0 && signalp6 -m fast -fmt txt' \
          ' -ff '+protein_filename+' -wp 1 -od '+signalp_directory \
          +' > '+log_filename+' 2>>'+log_filename
    status = os.system(cmd)
    if status == 0 :
        with open(signalp_filename,'w',encoding="utf8") as output :
            with open(signalp_directory+'/'+'prediction_results.txt','r',encoding="utf8") as file:
                for line in file :
                    line = line.rstrip()
                    if re.match("#",line) :
                        continue
                    liste = line.split('\t')
                    defline = liste[0].split()[0]
                    result = liste[1]
                    output.write(defline+"\t"+result+"\n")

    return status,cmd,genome

def running_signalp5(protein_filename,signalp_filename,log_filename,genome) :
    cmd = 'cd /env/cns/proj/agc/home/rmeheust/programs/signalp-5.0b/bin ;' \
          ' ./signalp -format short -org gram+ -fasta ' \
          +protein_filename+' -prefix '+signalp_filename \
          +' > '+log_filename+' 2>>'+log_filename
    print(cmd)
    status = os.system(cmd)
    return status, cmd, genome

def domtblout_is_ok(domtblout_filename) :

    with open(domtblout_filename,'r',encoding="utf8") as file :
        line = ''
        for line in file :
            line = line.rstrip()
    return line == '# [ok]'



def getting_tmhmm(genome_set,orf_set,tmhmm_directory) :
    print('reading tmhmm....')
    orf_to_tmhmm = {}
    for genome in genome_set :
        tmhmm_filename = tmhmm_directory+'/'+genome+'_protein.faa'+'.tmhmm'
        with open(tmhmm_filename,'r',encoding="utf8") as file :
            for line in file :
                line = line.rstrip()
                liste = line.split('\t')
                orf = liste[0]
                if orf not in orf_set :
                    continue
                tmhmm = int(liste[4].replace('PredHel=',''))
                orf_to_tmhmm[orf] = tmhmm
    return orf_to_tmhmm

def getting_signalp(genome_set,orf_set,signalp_directory) :
    orf_to_signalp = {}
    for genome in genome_set :
        signalp_filename = signalp_directory+'/'+genome+'_protein.faa.signalp_summary.signalp5'
        with open(signalp_filename,'r',encoding="utf8") as file :
            for line in file :
                line = line.rstrip()
                if re.match('#',line) :
                    continue
                liste = line.split('\t')
                orf = liste[0]
                if orf not in orf_set :
                    continue
                signalp = liste[1]
                orf_to_signalp[orf] = signalp
    return orf_to_signalp

def getting_kegg(genome_set,orf_set,kegg_directory,kegg_info_filename) :
    print('reading kegg....')
    accession_to_kegg_description = {}
    with open(kegg_info_filename,"r",encoding="utf8") as file :
        for line in file :
            line = line.rstrip()
            liste = line.split("\t")
            accession = liste[0]
            description = liste[-1]
            accession_to_kegg_description[accession] = description

    orf_to_kegg = defaultdict(list)
    for genome in genome_set :
        kegg_filename = kegg_directory+'/'+genome+'_protein.faa.kegg.domtblout.besthit'
        with open(kegg_filename,'r',encoding="utf8") as file :
            next(file)
            for line in file :
                line = line.rstrip()
                if re.match('#',line) :
                    continue
                liste = line.split()
                orf = liste[0]
                if orf not in orf_set :
                    continue
                accession = liste[1]
                if accession not in accession_to_kegg_description :
                    description =  accession
                else:
                    description = accession_to_kegg_description[accession]
                start,end = liste[3].split('-')
                c_evalue = liste[5] # conditional Evalue
                orf_to_kegg[orf].append( (  \
                                            int(start) , \
                                            int(end) , \
                                            accession , \
                                            description , \
                                            c_evalue \
                                        ) )

    orf_to_kegg_architecture = {}
    for orf,liste in orf_to_kegg.items() :
        architecture = []
        for domain_info_list in sorted( liste ,key=itemgetter(0,1) ) :
            architecture.append(domain_info_list[3]+" ("+domain_info_list[2]+")")
        orf_to_kegg_architecture[orf] = " + ".join(architecture)

    return orf_to_kegg_architecture

def getting_taxonomy(gtdb_taxonomy_filename,genome_set = None) :
    genome_to_taxonomy = {}
    with open(gtdb_taxonomy_filename,'r',encoding="utf8") as file :
        for line in file :
            line = line.rstrip()
            genome,lineage = line.split('\t')
            if genome_set  is None :
                genome_to_taxonomy[ genome ] = lineage
            else:
                if genome in genome_set :
                    genome_to_taxonomy[ genome ] = lineage
    return genome_to_taxonomy

def getting_pfam(pfam_filename,pfam_info_filename) :
    print('reading pfam....')
    name_to_pfam = {}
    pfam_to_name = {}
    with open(pfam_info_filename,"r",encoding="utf8") as file :
        for line in file :
            line = line.rstrip()
            accession,clan_accession,clan_name,name,description = line.split("\t")
            pfam_to_name[accession] = name
            name_to_pfam[name] = accession

    orf_to_pfam_architecture = {}
    orf_to_pfam = defaultdict(list)
    with open(pfam_filename,'r',encoding="utf8") as file :
        for line in file :
            line = line.rstrip()
            if re.match(r'#',line) :
                continue
            liste = line.split()
            orf_name = liste[0]
            pfam_name = liste[1]
            pfam_accession = name_to_pfam[pfam_name]
            start,end = liste[3].split('-')
            c_evalue = liste[5] # conditional Evalue
            orf_to_pfam[orf_name].append( ( \
                                            int(start) , \
                                            int(end) , \
                                            pfam_accession , \
                                            pfam_to_name[pfam_accession] , \
                                            c_evalue \
                                        ) )

    for orf,domain_list in orf_to_pfam.items() :
        architecture = []
        for domain_info_list in sorted( domain_list ,key=itemgetter(0,1) ) :
            architecture.append(domain_info_list[3]+" ("+domain_info_list[2]+")")
        orf_to_pfam_architecture[orf] = " + ".join(architecture)

    return orf_to_pfam_architecture

def getting_homemade_hmm(besthit_filename) :
    orf_to_pfam_architecture = {}
    orf_to_pfam = defaultdict(list)
    with open(besthit_filename,'r',encoding="utf8") as file :
        for line in file :
            line = line.rstrip()
            if re.match(r'#',line) :
                continue
            liste = line.split()
            orf_name = liste[0]
            pfam_name = liste[1]
            start,end = liste[3].split('-')
            c_evalue = liste[5] # conditional Evalue
            orf_to_pfam[orf_name].append( ( int(start) , int(end) , pfam_name , c_evalue ) )

    for orf,liste in orf_to_pfam.items() :
        architecture = []
        for domain_info_list in sorted( liste ,key=itemgetter(0,1) ) :
            architecture.append(domain_info_list[2])
        orf_to_pfam_architecture[orf] = " + ".join(architecture)

    return orf_to_pfam_architecture

def prodigal_to_feature_file(genome_to_path,feature_filename) :
    with open(feature_filename,'w',encoding="utf8") as output :
        for genome in genome_to_path :
            protein_filename = genome_to_path[genome]
            for record in SeqIO.parse(protein_filename,'fasta') :
                liste = record.description.split(' # ')
                orf = liste[0].split('>')[-1]
                scaffold = '_'.join( orf.split('_')[:-1] )
                start = liste[1]
                end = liste[2]
                strand = liste[3]
                output.write(orf+'\t'+genome+'\t'+scaffold+'\t'+start+'\t'+end+'\t'+strand+'\n')



if __name__ == "__main__":
    print('Hello world')
