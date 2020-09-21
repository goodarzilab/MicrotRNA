# MicrotRNA
Using k-mer deconvolution for rapid clustering of sequences into highly similar group and collapsing each group to tRNA families.


### Run MicrotRNA.
#### Create the metadata file
The metadata file resides in the working directory and lists the required information for each sample. 
For example:

| Sample.name  | R1  | R2  | 
|---|---|---|
Ecoli | 05_E-coli_S5_L001_R1_001.fastq.gz | 05_E-coli_S5_L001_R2_001.fastq.gz
Efaecalis | 07_E-faecalis_S7_L001_R1_001.fastq.gz | 07_E-faecalis_S7_L001_R2_001.fastq.gz
Etarda | 06_E-tarda_S6_L001_R1_001.fastq.gz | 06_E-tarda_S6_L001_R2_001.fastq.gz


#### Run the analysis
PWD=path to the working directory(required) 
```bash
python main.py $PWD -p -mk -lc -cls -th 20
```
#### Options
Run `python main.py` for usage.
The following are the options:
1. `-p` or `--process_fastqs` : preprocessing of raw fastq files : UMI extraction, adapter removal, mergeing and quality filtering
2. `-mk` or `--make_kmer_distance_mat`: making the kmer distance matrix for the samples
3. `-lc` or `--leiden_clustering` : clustering reads based on the kmer space of samples
4. `-cls` or `--Collabse_clusters`: Collabsing each read clusters to further simialr groups
5. `-th` or `--Number_of_Threads` : Number of cpu threads to use for clustering and collabsing ( default : 4)

