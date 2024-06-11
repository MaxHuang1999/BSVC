######################################################
#Author Max Huang
#E-mail 1849732267@qq.com
#2023.10.30
#SHAPE-MaP上游数据分析
######################################################
import os,argparse
import sys
import subprocess
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_ind
import numpy as np
#mapping_and_pairs

def pipeline1(refseq,DMSOinput,NAIinput,dotfile):
    refseqid=str(refseq).split('.')[0]
    DMSOid=str(DMSOinput).split('.')[0]
    NAIid=str(NAIinput).split('.')[0]
    print("===reads before filter===")
    os.system("expr $(cat %s | wc -l) / 4"%(DMSOinput))
    os.system("expr $(cat %s | wc -l) / 4"%(NAIinput))
    os.system("bowtie2-build %s %s > indexlog.txt"%(refseq,refseqid))
    os.system("rf-map -p 50 -wt 20 -bi %s -ctn -cmn 0 -cqo -cq5 20 -b2 -ow -o rf_map_refstructure %s %s > maplog.txt"%(refseqid,DMSOinput,NAIinput))
    print("===reads after filter===")
    os.system("samtools view -c rf_map_refstructure/%s.bam"%(DMSOid))
    os.system("samtools view -c rf_map_refstructure/%s.bam"%(NAIid))
    os.system("rf-count -m -rd -p 50 -o rf_count_refstructure -ow -f %s rf_map_refstructure/* > countlog.txt"%(refseq))
    os.system("mkdir rf_rctools_refstructure")
    os.system("rf-rctools view -t rf_count_refstructure/%s.rc > rf_rctools_refstructure/%s.txt"%(DMSOid,DMSOid))
    os.system("rf-rctools view -t rf_count_refstructure/%s.rc > rf_rctools_refstructure/%s.txt" % (NAIid,NAIid))
    os.system("rf-norm -u rf_count_refstructure/%s.rc -t rf_count_refstructure/%s.rc -sm 3 -nm 3 -rb AC -o rf_norm_refstructure -ow > normlog.txt"%(DMSOid,NAIid))
    os.system("rf-jackknife -r %s -p 50 -rp '-md 600 -nlp' -x -ow -o rf_jackknife_refstructure rf_norm_refstructure/* > jackknifelog.txt"%(dotfile))


# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="vRICseq data analysis pipeline-RNAseq. First-writtern by Max.")
#     parser.add_argument("-r", "--refseq", type=str, required=True, help="Please type in the file of refseq.")
#     parser.add_argument("-f", "--fdata", type=str, required=True,
#                         help="Please type in the file of first sequencing data")
#     parser.add_argument("-s", "--sdata", type=str, required=True,
#                         help="Please type in the file of second sequencing data")
#     parser.add_argument("-t", "--title", type=str, required=True,
#                         help="Please type in the name of outputfile;for example:SARSCOV2-1")
#     Args = parser.parse_args()
#     shapelikepipeline(os.path.abspath(Args.refseq), os.path.abspath(Args.fdata), os.path.abspath(Args.sdata),
#                  os.path.abspath(Args.title))







