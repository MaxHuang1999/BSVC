#######################################################
#Data 2024-1-28
#Author Max Huang
#E-mail 1849732267@qq.com
#脚本用于处理二维SHAPE-MaP数据
#######################################################
#

import os,argparse
import sys
import subprocess

def rnaseq(unmod1,unmod2,mod1,mod2,title,refseq):
    os.system("bowtie2-build %s %s"%(refseq,title))
    os.system("rf-map -p 50 -wt 20 -bi %s -ctn -cmn 0 -cqo -cq5 20 -b2 -ow -o rf_map ./unmod/%s.fastq ./unmod/%s.fastq ./mod/%s.fastq ./mod/%s.fastq"%(title,unmod1,unmod2,mod1,mod2))
    os.system("bedtools genomecov -ibam ./rf_map/%s.bam  -dz -split > %s.coverage.txt"%(unmod1,unmod1))
    os.system("bedtools genomecov -ibam ./rf_map/%s.bam  -dz -split > %s.coverage.txt"%(unmod2,unmod2))
    os.system("bedtools genomecov -ibam ./rf_map/%s.bam  -dz -split > %s.coverage.txt"%(mod1,mod1))
    os.system("bedtools genomecov -ibam ./rf_map/%s.bam  -dz -split > %s.coverage.txt"%(mod2,mod2))
    os.system("samtools view -c ./rf_map/%s.bam" % (unmod1))
    os.system("samtools view -c ./rf_map/%s.bam" % (unmod2))
    os.system("samtools view -c ./rf_map/%s.bam" % (mod1))
    os.system("samtools view -c ./rf_map/%s.bam" % (mod2))
    os.system("rf-count -m -rd -p 50 -o rf_count -ow -f %s rf_map/*" % (refseq))
    os.system("mkdir rf_rctools")
    os.system("rf-rctools view -t ./rf_count/%s.rc > ./rf_rctools/%s.txt" % (unmod1,unmod1))
    os.system("rf-rctools view -t ./rf_count/%s.rc > ./rf_rctools/%s.txt" % (unmod2,unmod2))
    os.system("rf-rctools view -t ./rf_count/%s.rc > ./rf_rctools/%s.txt" % (mod1,mod1))
    os.system("rf-rctools view -t ./rf_count/%s.rc > ./rf_rctools/%s.txt" % (mod2,mod2))
    os.system("rf-norm -u ./rf_count/%s.rc ./rf_count/%s.rc -t ./rf_count/%s.rc ./rf_count/%s.rc -sm 3 -nm 3 -rb AC -o rf_norm -ow -p 50"%(unmod1,unmod2,mod1,mod2))
    os.system("rf-fold ./rf_norm/%s.xml -ct -w -fw 3000 -fo 300 -wt 200 -pw 1000 -po 250 -dp -sh -nlp -md 600 -o rf_fold -g -p 50 -ow"%(title))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="SHAPE-MaP data analysis pipeline. First-writtern by Max.")
    parser.add_argument("-u1", "--unmod1", type=str, required=True, help="Please type in the file of unmod1 id.")
    parser.add_argument("-u2", "--unmod2", type=str, required=True,help="Please type in the file of unmod2 id .")
    parser.add_argument("-m1", "--mod1", type=str, required=True,help="Please type in the file of mod1 id .")
    parser.add_argument("-m2", "--mod2", type=str, required=True,help="Please type in the file of mod2 id.")
    parser.add_argument("-t", "--title", type=str, required=True,help="Please type in the name of cov id;for example:NC_045512.2.")
    parser.add_argument("-r", "--refseq", type=str, required=True,help="Please type in the file of refseq.")
    Args = parser.parse_args()
    rnaseq(Args.unmod1,Args.unmod2,Args.mod1,Args.mod2,Args.title,Args.refseq)