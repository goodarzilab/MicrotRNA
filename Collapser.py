import numpy as np
import distance,os,time,itertools,subprocess,operator,statistics,random, pysam
from Bio import pairwise2
from sklearn.decomposition import PCA
from multiprocessing import Pool
import pandas as pd
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from numpy import array, linspace
import kmeans1d

def call_calc_hamming_distance(names_list,n_threads):
        iterlist = names_list
        caller = Pool(n_threads)
        caller.map(calc_hamming_distance,iterlist)

def Collapse(listforcollapse):
    returnlist =[]
    trna = {}
    score = 0
    current_collapse = ""
    for seq in listforcollapse:
        trna[seq.split("#")[1]] = seq  # count : seq
    cntlist = [int(x) for x in sorted(trna.keys())]
    clusters, centroids = kmeans1d.cluster(cntlist, 2)
    comparetor = max(centroids)
    cseqs_keys = [cnt for cnt in cntlist if cnt/comparetor > 2/3]
    cseqs = [trna[str(k)] for k in cseqs_keys]           ##important: if we need the to collapsed reads count we need add stuff here  but since we don't use them in downstream analyses I ignored it
    return cseqs

def calc_hamming_distance(clustername):
    tRNA_repertoire = []
    tRNAseq = open("MSA_fastas/"+clustername[0])
    tRNAseqAligned = open("MSA_fastas/"+clustername[1])
    seqdict , aligneddict = {},{}
    for i,line in enumerate(tRNAseq.readlines()):
        if i%2 == 0:
            tRNAname = line[1:].strip()
        else:
            seqdict[tRNAname] = line.strip()
    for i,line in enumerate(tRNAseqAligned.readlines()):
        if i%3 == 0:
            tRNAname = line[1:].strip()
        elif i%3 == 1:
            firstline = line.strip()
        elif i%3 == 2:
            fullseq = firstline+line.strip()
            aligneddict[tRNAname] = fullseq
    keys = list(seqdict.keys())
    had_collapsed = [False for key in keys]
    for i in range(len(keys)):
        if had_collapsed[i]:
            continue
        had_collapsed[i] = True
        list_for_collapse = []
        list_for_collapse.append(seqdict[keys[i]]+"#"+keys[i].split(":")[1])
        for j in range(i+1, len(keys)):
            if had_collapsed[j]:
                continue
            try:
                dis = distance.hamming(aligneddict[keys[i]],aligneddict[keys[j]])
                if dis < 5 :
                    list_for_collapse.append(seqdict[keys[j]]+"#"+keys[j].split(":")[1])
                    had_collapsed[j] = True
            except:
                continue
        collapsed_Seq = Collapse(list_for_collapse)
        for tRNA in collapsed_Seq:
            tRNA_repertoire.append(tRNA)
    filename = "tRNAs_after_collapse/tRNA-" + str(random.randint(1, 10000000))+".rfa"
    fl= open(filename, "w")
    for t in tRNA_repertoire:
        st = t+"\n"
        fl.write(st)
    fl.close()


def P_Alignment(cluster):
    fname = "MSA_fastas/tRNA-" + str(random.randint(1, 1000000))+"_seq_.fa"
    with open(fname, "w") as file:
        for i,l in enumerate(cluster):
            entry = ">tRNA-"+str(i)+"-count:"+l.split("#")[1] + "\n"+ l.split('#')[0] +"\n"
            file.write(entry)
        file.close()
        faligned = fname+"_aligned.fa"
        cline = "muscle -in %s -out %s " %(fname,faligned)
        subprocess.call(cline,shell=True)
def call_alignment(wholeclusters,n_threads=4):
    iterlist = wholeclusters
    caller = Pool(n_threads)
    caller.map(P_Alignment, iterlist)
def RunCollapser(log,clustersToOpen,tRNAlistFileToOpen,n_threads):
    n_threads = n_threads
    clusters = np.load(clustersToOpen, allow_pickle=True)
    wholecluster = []
    with open(tRNAlistFileToOpen, "r+") as file:
        reads = file.readlines()
    trna = [x.strip() for x in reads]
    for cluster in clusters:
        tcluster = []
        for ix in cluster:
            rd =trna[ix]
            tcluster.append(rd)
        wholecluster.append(tcluster)
    log.write("#### Running Multiple sequencing alignment with MUSCLE#####\n")
    log.write("\n\n")
    print("#### Running Multiple sequencing alignment with MUSCLE#####\n")
    call_alignment(wholecluster, n_threads)
    filenames = list(set(name.split("_")[0] for name in set(os.listdir('MSA_fastas/'))))
    cluster_names = [(name+"_seq_.fa", name+"_seq_.fa_aligned.fa") for name in filenames]
    print("#### Collabsing the tRNAs in same cluster based on their Hamming distance #####\n")
    call_calc_hamming_distance(cluster_names,n_threads)
