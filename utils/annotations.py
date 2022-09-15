# coding: utf-8

from Bio import SeqIO
import os,sys


pfam_hmm_filename = '/env/cns/proj/projet_CSX/scratch/bank/pfam/release_35.0/Pfam-A.hmm'
pfam_info_filename = '/env/cns/proj/projet_CSX/scratch/bank/pfam/release_35.0/Pfam-A.clans.tsv'
kegg_hmm_filename = '/env/cns/proj/projet_CSX/scratch/bank/kegg/release_100.0/KOFAM.hmm'
kegg_info_filename = '/env/cns/proj/projet_CSX/scratch/bank/kegg/release_100.0/ko_list'

def runningTmhmm(protein_filename,tmhmm_filename,genome):
    cmd = '/env/cns/proj/agc/home/rmeheust/programs/tmhmm-2.0c/bin/tmhmm --short '+protein_filename+'  > '+tmhmm_filename
    status = os.system(cmd)
    return status,cmd,genome

def runningPfam(protein_filename,pfam_hmm_filename,domtblout_filename,genome):
    cmd1 = 'hmmsearch -E 1e-3 --cpu '+str(1)+' --domtblout '+domtblout_filename+' '+pfam_hmm_filename+' '+protein_filename+' >/dev/null 2>/dev/null' 
    status1 = os.system(cmd1)
    
    besthit_filename = domtblout_filename+'.besthit'
    cmd2 = '/env/cns/proj/agc/home/rmeheust/scripts/usefulScripts/cath-resolve-hits.ubuntu14.04 --input-format hmmer_domtblout '+domtblout_filename+' >'+besthit_filename
    status2 = os.system(cmd2)
    os.remove(domtblout_filename)
    return status1,cmd1,status2,cmd2,genome

def runningKegg(protein_filename,kegg_hmm_filename,domtblout_filename,genome):
    cmd1 = 'hmmsearch -E 1e-3 --cpu '+str(1)+' --domtblout '+domtblout_filename+' '+kegg_hmm_filename+' '+protein_filename+' >/dev/null 2>/dev/null' 
    status1 = os.system(cmd1)
    
    besthit_filename = domtblout_filename+'.besthit'
    cmd2 = '/env/cns/proj/agc/home/rmeheust/scripts/usefulScripts/cath-resolve-hits.ubuntu14.04 --input-format hmmer_domtblout '+domtblout_filename+' >'+besthit_filename
    status2 = os.system(cmd2)
    os.remove(domtblout_filename)
    return status1,cmd1,status2,cmd2,genome

def runningSignalp6(protein_filename,signalp_filename,signalp_directory,genome) :
    log_filename = '/dev/null'
    cmd = 'source /env/cns/proj/agc/scratch/conda/miniconda.profile && conda activate signalp-6.0 && signalp6 -m fast -fmt txt  -ff '+protein_filename+' -wp 1 -od '+signalp_directory+' > '+log_filename+' 2>>'+log_filename
    status = os.system(cmd)
    if status == 0 :
        output=open(signalp_filename,'w')
        file = open(signalp_directory+'/'+'prediction_results.txt','r')
        for line in file :
            line = line.rstrip()
            if re.match("#",line) :
                continue
            else :
                liste = line.split('\t')
                defline = liste[0].split()[0]
                result = liste[1]
                output.write(defline+"\t"+result+"\n")
        file.close()
        output.close()
    return status,cmd,genome

def runningSignalP5(protein_filename,signalp_filename,log_filename,genome) :
    cmd = 'cd /env/cns/proj/agc/home/rmeheust/programs/signalp-5.0b/bin ; ./signalp -format short -org gram+ -fasta '+protein_filename+' -prefix '+signalp_filename+' > '+log_filename+' 2>>'+log_filename
    print(cmd)
    status = os.system(cmd)
    return status, cmd, genome

def domtblout_is_ok(domtblout_filename) :
    file = open(domtblout_filename)
    line = ''
    for line in file :
        line = line.rstrip()
    file.close()
    if line == '# [ok]' :
        return True
    else:
        return False


def gettingTMHMM(genomeSet,orfSet,tmhmm_directory) :
    print('reading tmhmm....')
    orf2tmhmm = dict()
    for genome in genomeSet :
        tmhmm_filename = tmhmm_directory+'/'+genome+'_protein.faa'+'.tmhmm'
        file = open(tmhmm_filename,'r')
        for line in file :
            line = line.rstrip()
            liste = line.split('\t')
            orf = liste[0]
            if orf not in orfSet :
                continue
            tmhmm = int(liste[4].replace('PredHel=',''))
            orf2tmhmm[orf] = tmhmm
        file.close()
    return orf2tmhmm

def gettingSignalP(genomeSet,orfSet,signalp_directory) :
    orf2signalp = dict()
    for genome in genomeSet :
        signalp_filename = signalp_directory+'/'+genome+'_protein.faa.signalp_summary.signalp5'
        file = open(signalp_filename,'r')
        for line in file :
            line = line.rstrip()
            if re.match('#',line) :
                continue
            liste = line.split('\t')
            orf = liste[0]
            if orf not in orfSet :
                continue
            signalp = liste[1]
            orf2signalp[orf] = signalp
        file.close()
    return orf2signalp

def gettingKegg(genomeSet,orfSet,kegg_directory,kegg_info_filename) :
    print('reading kegg....')
    accession2keggDescription = dict()
    file = open(kegg_info_filename,"r")
    for line in file :
        line = line.rstrip()
        liste = line.split("\t")
        accession = liste[0]
        description = liste[-1]
        accession2keggDescription[accession] = description
    file.close()

    orf2keggArchitecture = dict()
    orf2kegg = defaultdict(list)
    for genome in genomeSet :
        kegg_filename = kegg_directory+'/'+genome+'_protein.faa.kegg.domtblout.besthit'
        file = open(kegg_filename,'r')
        next(file)
        for line in file :
            line = line.rstrip()
            if re.match('#',line) :
                continue
            liste = line.split()
            orf = liste[0]
            if orf not in orfSet :
                continue
            accession = liste[1]
            if accession not in accession2keggDescription :
                description =  accession
            else:
                description = accession2keggDescription[accession]
            start,end = liste[3].split('-')
            cEvalue = liste[5] # conditional Evalue
            orf2kegg[orf].append( ( int(start) , int(end) , accession , description , cEvalue ) )
        file.close()

    orf2keggArchitecture = dict()
    for orf,liste in orf2kegg.items() :
        architecture = list()
        for domainInfoList in sorted( liste ,key=itemgetter(0,1) ) :
            architecture.append(domainInfoList[3]+" ("+domainInfoList[2]+")")
        orf2keggArchitecture[orf] = " + ".join(architecture)
    return orf2keggArchitecture


def gettingTaxonomy(genomeSet,gtdb_taxonomy_filename) :
    genome2taxonomy = dict()
    file = open(gtdb_taxonomy_filename,'r')
    for line in file :
        line = line.rstrip()
        genome,lineage = line.split('\t')
        if genome in genomeSet :
            genome2taxonomy[ genome ] = lineage
    file.close()
    return genome2taxonomy

def gettingPfam(genomeSet,orfSet,pfam_directory,pfam_info_filename) :
    print('reading pfam....')
    name2pfam = dict()
    pfam2name = dict()
    file = open(pfam_info_filename,"r")
    for line in file :
        line = line.rstrip()
        accession,clanAccession,clanName,name,description = line.split("\t")
        pfam2name[accession] = name
        name2pfam[name] = accession
    file.close()

    orf2pfamArchitecture = dict()
    orf2pfam = defaultdict(list)
    for genome in genomeSet :
        pfam_filename = pfam_directory+'/'+genome+'_protein.faa.pfam.domtblout.besthit'
        file = open(pfam_filename,'r')
        for line in file :
            line = line.rstrip()
            if re.match(r'#',line) :
                continue
            liste = line.split()
            orfName = liste[0]
            if orfName not in orfSet :
                continue
            pfamName = liste[1]
            pfamAccession = name2pfam[pfamName]
            start,end = liste[3].split('-')
            cEvalue = liste[5] # conditional Evalue
            orf2pfam[orfName].append( ( int(start) , int(end) , pfamAccession , pfam2name[pfamAccession], cEvalue ) )
        file.close()

    for orf,liste in orf2pfam.items() :
        architecture = list()
        for domainInfoList in sorted( liste ,key=itemgetter(0,1) ) :
            architecture.append(domainInfoList[3]+" ("+domainInfoList[2]+")")
        orf2pfamArchitecture[orf] = " + ".join(architecture)
    return orf2pfamArchitecture


def prodigal2featureFile(prodigal_filename,feature_filename,genome2path) :
    output = open(feature_filename,'w')
    for genome in genome2path :
        protein_filename = genome2path[genome]
        for record in SeqIO.parse(protein_filename,'fasta') :
            liste = record.description.split(' # ')
            orf = liste[0].split('>')[-1]        
            scaffold = '_'.join( orf.split('_')[:-1] )
            start = liste[1]
            end = liste[2]
            strand = liste[3]
            output.write(orf+'\t'+genome+'\t'+scaffold+'\t'+start+'\t'+end+'\t'+strand+'\n')
    output.close()



if __name__ == "__main__":
    print(annotations.__name__)
