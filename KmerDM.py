import os, sys, time, subprocess, pysam
import numpy as np
from multiprocessing import Pool
import pandas as pd
from scipy.spatial import distance_matrix
from sklearn.metrics import pairwise_distances
from scipy.sparse import dok_matrix
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.neighbors import kneighbors_graph
import leidenalg as la
from igraph import *
class MakeKmerDM():
    def __init__(self,fqfile,log):
        self.file = fqfile
        self.log = log
        self.seqlist = []
        self.klist = []
        self.reads = {}
        self.matrix = None
        self.distace_matrix = None
        self.clustered= {}
    def readfile(self):
        print(self.file)
        log = self.log
        start = time.time()
        log.write("####making the read-kmer dictionary####\n")
        with pysam.FastxFile(self.file) as fin:
            for entry in fin:
                seqcount = entry.name.split("#")[1]
                nseq = entry.sequence+"#"+seqcount
                self.reads[nseq] = {}
                for p in range(1+len(entry.sequence)-6):               ##slide over the sequence to get k-mers,  p as position
                    kmer = entry.sequence[p:p+6]
                    if kmer not in self.klist :
                        self.klist.append(kmer)
                    if kmer not in self.reads[nseq].keys():
                        self.reads[nseq][kmer] = 1
                    elif kmer in self.reads[nseq].keys():
                        self.reads[nseq][kmer] += 1            ##making this dict >>> {sequence#n : { kmer : kmercount}

        log.write("####read-kmer dictionary is done. now making the matrix####\n")
        log.write("\n\n")
    def make_distance_matrix(self):
        log = self.log
        coded_tRNAs, coded_kmers = {},{}
        tfile = open(self.file.split(".")[0]+"_tRNAlist.txt", "w")
        for i,tRNA in enumerate(self.reads.keys()):
            coded_tRNAs[tRNA] = i
            tfile.write(tRNA)
            tfile.write("\n")
        for j,kmer in enumerate(self.klist):
            coded_kmers[kmer] = j
        self.matrix = dok_matrix((len(coded_kmers),len(coded_tRNAs)), dtype=np.int)
        for read,kmers in self.reads.items():
            for kmer in kmers.keys():
                self.matrix[coded_kmers[kmer],coded_tRNAs[read]] = self.reads[read][kmer]
        self.distance_matrix = pairwise_distances(self.matrix.T, metric='euclidean' ,n_jobs= -1)
        np.save(self.file.split(".")[0]+"_DistanceMatrix", self.distance_matrix)
        log.write("####read-kmer space distance matrix made####\n")
        log.write("\n\n")
class Leiden():
    def __init__(self,log,matrix,n_threads):
        self.log = log
        self.sampname = matrix
        self.NGraph = np.load(matrix)
        self.n_threads = n_threads

    def kneighbors_graph(self):
        log = self.log
        self.log = log
        log.write("####Leiden clustering is running on Kmer space of sample - this might tale a while####\n")
        print("Leiden clustering is running on Kmer space of {} - this might take a while".format(self.sampname))
        grf = Graph()
        self.NGraph = kneighbors_graph(self.NGraph ,n_neighbors= 14,  mode = 'distance', n_jobs= self.n_threads )
        self.NGraph = self.NGraph.toarray().tolist()
        self.NGraph = grf.Adjacency(self.NGraph)
        self.NGraph = la.find_partition(self.NGraph,  la.ModularityVertexPartition)
        self.NGraph = np.array(self.NGraph)
        name_to_save = self.sampname.split("_")[0] + "_Final_leiden_clusters"
        np.save(name_to_save, self.NGraph)
        log.write("\n\n")
        log.write("Leiden clustering is done!")
        return self.NGraph
