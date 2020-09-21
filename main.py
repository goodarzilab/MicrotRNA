import os,sys,re,argparse
import pandas as pd
import shutil
from ProcessFastqs import Process
from KmerDM import MakeKmerDM, Leiden
from Collapser import *
from pathlib import Path

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument("path", help=" directory path to the raw data and and metadatafile ", type = str)
  parser.add_argument("-p", "--process_fastqs", help="genral processing of raw files", dest='process', action='store_true')
  parser.add_argument("-mk", "--make_kmer_distance_mat", help="making the kmer distance matrix for the samples", dest='kmer', action='store_true')
  parser.add_argument("-lc", "--leiden_clustering", help="making the kmer distance matrix for the samples", dest='clusterit', action='store_true')
  parser.add_argument("-cls", "--Collabse_clusters", help="Collabse the the tRNAs in clusters", dest='collabse', action='store_true')
  parser.add_argument("-th", "--Number_of_Threads", help="number of threads for MSA and collabse part (Default is 4)", dest='threads', type= int)
  parser.set_defaults(threads=4)
  args = parser.parse_args()

  data_path=args.path
  os.chdir(data_path)
  meta = pd.read_csv("metadata", sep=",", header=0, index_col=0)
  log = open("pipeline.txt", "w+")
  if args.process :
      pro = Process(meta,log=log)
      pro.UTM()
      pro.Collapser()
      pro.Filter_nontRNAs()
  if args.kmer: ##change the output names for this part
      for sample in meta.index:
          print("making the kmer space of sample {}".format(sample))
          orgsample = sample
          processedsamp = orgsample + ".post.c.filtered.fastq.gz"
          try:
              run = MakeKmerDM(processedsamp,log)
              run.readfile()
              run.make_distance_matrix()
              del(run)
          except:
              continue
  if args.clusterit:
      for sample in meta.index:
          print("loading sample {} \n".format(sample))
          sample_to_cluster = sample + "_DistanceMatrix.npy"
          ld = Leiden(log,sample_to_cluster,n_threads=args.threads)
          ld.kneighbors_graph()
          del(ld)

  if args.collabse:
      for sample in meta.index:
          clustersToOpen = sample + "_Final_leiden_clusters.npy"
          tRNAlistFileToOpen = sample + "_tRNAlist.txt"
          Path(sample).mkdir(parents=True, exist_ok=True)
          shutil.copy(os.getcwd()+"/"+clustersToOpen, os.getcwd()+"/"+sample)
          shutil.copy(os.getcwd()+"/"+tRNAlistFileToOpen, os.getcwd()+"/"+sample)
          os.chdir(sample)
          Path("MSA_fastas").mkdir(parents=True, exist_ok=True)
          Path("tRNAs_after_collapse").mkdir(parents=True, exist_ok=True)
          RunCollapser(log,clustersToOpen,tRNAlistFileToOpen,n_threads=args.threads)
          os.chdir("..")
  log.close()
