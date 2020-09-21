import subprocess, re, os, tempfile, pysam, glob
import pandas as pd

class Process():
    def __init__(self, metadata, log, isPaired=True, useRead="R2", hasUMI=True, runMode=True):
        self.metadata = metadata
        self.log = log
        self.paired = isPaired
        self.read = useRead
        self.umi = hasUMI
        self.rm = runMode

    def Filter_nontRNAs(self):
        meta = self.metadata
        log = self.log
        log.write("####Filtering out the sequence errors and non tRNA reads#####\n")
        for file in glob.iglob('**/*.post.c.fastq.gz', recursive=True):
            fq = file
            fout = re.sub(".fastq.gz",".filtered.fastq.gz",fq)
            with pysam.FastxFile(fq) as fin, open(fout, mode='w') as fout:
                for entry in fin:
                    if int(entry.name.split("#")[1]) > 2 and entry.sequence[-3:] =="CCA":
                        eout = str(entry)+"\n"
                        fout.write(eout)
        log.write("\n\n")
        meta.loc[file, 'fastq'] = fout
    def Collapser(self):
        meta = self.metadata
        log = self.log
        log.write("####Collapsing identical reads#####\n")
        files = os.listdir()
        fqfiles =  sorted(filter(lambda x: x[-14:] == ".post.fastq.gz", files))
        for sample in fqfiles:
          print(sample)
          fq =sample
          tmpdir = tempfile.mkdtemp()
          cmd = "gunzip -c %s | awk \'{ if(NR%%4==1) {print $1} } \' > %s/id" % (fq, tmpdir)
          print(cmd)
          log.write("%s\n" % (cmd))
          if self.rm: subprocess.call(cmd,shell=True)
          cmd = "gunzip -c %s | awk \'{ if(NR%%4==2) {print $1} } \' > %s/seq" % (fq, tmpdir)
          print(cmd)
          log.write("%s\n" % (cmd))
          if self.rm: subprocess.call(cmd,shell=True)
          cmd = "gunzip -c %s | awk \'{ if(NR%%4==0) {print $1} } \' > %s/qual" % (fq, tmpdir)
          print(cmd)
          log.write("%s\n" % (cmd))
          if self.rm: subprocess.call(cmd,shell=True)
          out = re.sub(".fastq.gz", ".c.fastq.gz", fq)
          meta.loc[sample, 'fastq'] = out
          cmd = "paste %s/id %s/qual %s/seq | sort -k 3 | uniq -f 2 -c | awk '{print $2 \"#\" $1 \"\\n\" $4 \"\\n+\\n\" $3 }' | gzip -c > %s"  %(tmpdir,tmpdir,tmpdir,out)
          print(cmd)
          log.write("%s\n" % (cmd))
          if self.rm: subprocess.call(cmd,shell=True)
          log.write("\n\n")

    def UTM(self): ##umi-extract,trimming,and,merging

        if (not self.umi): return True

        else:

            meta = self.metadata
            log = self.log
            log.write("####Extracting UMIs#####\n")

            if (self.paired):
                for sample in meta.index:
                    print(sample)
                    r1 = meta.loc[sample, 'R1']
                    r2 = meta.loc[sample, 'R2']
                    o1 = re.sub(".fastq.gz", ".u.fastq.gz", r1)
                    o2 = re.sub(".fastq.gz", ".u.fastq.gz", r2)
                    meta.loc[sample, 'R1'] = o1
                    meta.loc[sample, 'R2'] = o2
                    cmd = 'umi_tools extract --stdin={} --read2-in={} --bc-pattern=NNNN --bc-pattern2=NNNN -L log.out --stdout={} --read2-out={}'.format(r1, r2, o1, o2)
                    print(cmd)
                    log.write("%s\n" % (cmd))
                    if self.rm: subprocess.call(cmd,shell=True)
            else:
                for sample in meta.index:
                    r1 = meta.loc[sample, 'R1']
                    o1 = re.sub(".fastq.gz", ".c.fastq.gz", r1)
                    meta.loc[sample, 'R1'] = o1
                    cmd = 'umi_tools extract --stdin={} --bc-pattern=NNNN -L log.out --stdout={}'.format(r1, o1)
                    print(cmd)
                    log.write("%s\n" % (cmd))
                    if self.rm: subprocess.call(cmd,shell=True)
            log.write("\n\n")

        log.write("####Removing adaptors#####\n")

        if (self.paired and self.read=="both"):
            for sample in meta.index:
                r1 = meta.loc[sample, 'R1']
                r2 = meta.loc[sample, 'R2']
                o1 = re.sub(".fastq.gz", ".trim.fastq.gz", r1)
                o2 = re.sub(".fastq.gz", ".trim.fastq.gz", r2)
                cmd = 'cutadapt -e 0.12 -m 20 -q 15 -a ^ACTGGATACTGGN...GTATCCAGT -A ^ACTGGATAC...NCCAGTATCCAGT -o {} -p {} {} {}'.format(o1, o2, r1, r2)
                meta.loc[sample, 'R1'] = o1
                meta.loc[sample, 'R2'] = o2
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        elif (self.paired and self.read=="R2"):
            for sample in meta.index:
                r2 = meta.loc[sample, 'R2']
                print(r2)
                o2 = re.sub(".fastq.gz", ".trim.fastq.gz", r2)
                cmd = 'cutadapt -j 8 --trimmed-only  -e 0.12 -m 55 -q 15 -a ^ACTGGATAC...NCCAGTATCCAGT -o {} {}'.format(o2,r2)
                meta.loc[sample, 'R2'] = o2
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        else:
            for sample in meta.index:
                r1 = meta.loc[sample, 'R1']
                o1 = re.sub(".fastq.gz", ".trim.fastq.gz", r1)
                cmd = 'cutadapt -j 8 --trimmed-only  -e 0.12 -m 55 -q 15 -a ^ACTGGATACTGGN...GTATCCAGT -o {} {}'.format(o1,r1)
                meta.loc[sample, 'R1'] = o1
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        log.write("\n\n")

        log.write("####Merging reads if required#####\n")

        if (self.paired and self.read=="both"):
            for sample in meta.index:
                r1 = meta.loc[sample, 'R1']
                r2 = meta.loc[sample, 'R2']
                cmd = 'pear -n 30 -f {} -r {} -o {} '.format(r1, r2, sample)
                cmd = cmd + ' ; rm *.discarded.fastq; rm *.unassembled.*; '
                cmd = cmd + ' ; mv {}.assembled.fastq {}.post.fastq; gzip {}.post.fastq'.format(sample, sample, sample)
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        elif (self.paired and self.read=="R2"):
            for sample in meta.index:
                r2 = meta.loc[sample, 'R2']
                print(r2)
                cmd = 'zcat {} | fastx_reverse_complement -z -i - -o {}.post.fastq.gz'.format(r2,sample)
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        else:
            for sample in meta.index:
                r1 = meta.loc[sample, 'R1']
                cmd = 'mv {} {}.post.fastq.gz'.format(r1, sample)
                print(cmd)
                log.write("%s\n" % (cmd))
                if self.rm: subprocess.call(cmd,shell=True)
        log.write("\n\n")
