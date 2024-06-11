#######################################################
#Data 2024-1-23
#Author Max Huang
#E-mail 1849732267@qq.com
#脚本用于处理一维RNA-seq数据
#######################################################
#提前准备好fastq文件、用于比对的病毒基因组序列、病毒和宿主合并的序列、病毒基因组注释文件、
#合并冠状病毒和宿主基因组序列
#(cat SARS-CoV-2.fa ; zcat genome.fa.gz) | fasta_formatter -w 72 | bgzip -@ 20 -c /dev/stdin > SARS2_host.genome.fa.gz

import os,argparse
import sys
import subprocess

def rnaseq(refseq,fdata,sdata,title):
    os.system("fastp -i %s -I %s -o %s.read1.clean.fq -O %s.read2.clean.fq -z 4 -q 20 -u 30 -f 10" % (fdata, sdata, title, title))
    os.system("mkdir -p %s_host.genome.star"%(title))
    os.system("STAR --runThreadN 50 --runMode genomeGenerate --genomeDir %s_host.genome.star --genomeFastaFiles %s" % (title,refseq))
    os.system("STAR --runThreadN 50 --genomeDir %s_host.genome.star --readFilesIn %s.read1.clean.fq %s.read2.clean.fq --outFilterType BySJout --outFilterMultimapNmax 20 \
              --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outSJfilterOverhangMin 12 12 12 12 --outSJfilterCountUniqueMin 1 1 1 1 --outSJfilterCountTotalMin 1 1 1 1 --outSJfilterDistToOtherSJmin 0 0 0 0 \
              --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --scoreGapNoncan -4 --scoreGapATAC -4 --chimOutType Junctions WithinBAM HardClip --chimScoreJunctionNonGTAG 0 \
              --alignSJstitchMismatchNmax -1 -1 -1 -1 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --chimSegmentMin 20 --outFileNamePrefix %s_toGenome --outSAMmultNmax 32 \
              --outSAMtype BAM SortedByCoordinate  --outSAMunmapped Within KeepPairs"%(title,title, title,title))
    os.system("ln -s %s_toGenomeAligned.sortedByCoord.out.bam %s.genome.bam"%(title,title))
    os.system("samtools index %s.genome.bam"%(title))
    os.system("samtools view -@ 50 -b -o %s.virus.bam %s.genome.bam KT336560.1"%(title,title))
    os.system("samtools view -c -f 1 -F 12 %s.virus.bam"%(title))
    command1 ='''samtools view %s.virus.bam | awk 'BEGIN {{ OFS="\t"; }} {{ if ($6 ~ /N/) print $4, $6; }}' | bgzip -c -@ 4 /dev/stdin > %s.jump-cigars.txt.gz'''%(title,title)
    subprocess.call(command1, shell=True)
    command2 = '''zcat %s.jump-cigars.txt.gz | python convert-cigars.py | sort -k1,2n | uniq -c | awk '{{ OFS="\t"; }} {{ print $2, $3, $1; }}' | bgzip -c /dev/stdin > %s.jumps.txt.gz''' % (title, title)
    subprocess.call(command2, shell=True)
    os.system("gunzip %s.jumps.txt.gz" % (title))
    os.system("bedtools genomecov -ibam %s.virus.bam  -dz -split > %s.virus.coverage.txt"%(title,title))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNAseq data analysis pipeline. First-writtern by Max.")
    parser.add_argument("-r", "--refseq",type=str,required=True,help="Please type in the file of refseq.")
    parser.add_argument("-f","--fdata",type=str,required=True,help="Please type in the file of first sequencing data")
    parser.add_argument("-s", "--sdata",type=str,required=True,help="Please type in the file of second sequencing data")
    parser.add_argument("-t", "--title", type=str, required=True,help="Please type in the name of cov id;for example:NC_045512.2")
    Args = parser.parse_args()
    rnaseq(os.path.abspath(Args.refseq),os.path.abspath(Args.fdata),os.path.abspath(Args.sdata),Args.title)