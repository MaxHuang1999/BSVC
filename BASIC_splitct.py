#######################################################
#Data 2024-1-24
#Author Max Huang
#E-mail 1849732267@qq.com
#该脚本用于从全基因组ct文件中提取出目的区域
#######################################################
import pandas as pd
import numpy as np
import os,argparse
import sys
import subprocess
import shutil

def select_strc(ct_f, s, e, types = 'rna_strc'):
    strc_selet = ct_f[s-1:e].copy()
    mask = (strc_selet['5'] < s) | (strc_selet['5'] >= e)
    strc_selet['5'][mask] = 0
    if types == 'rna_strc':
        strc_selet['1'] = strc_selet['1'] - s +1
        strc_selet['3'] = strc_selet['3'] - s +1
        strc_selet['4'] = strc_selet['4'] - s +1
        strc_selet['5'][strc_selet['5']>0] = strc_selet['5'][strc_selet['5']>0] - s +1
        strc_selet['6'] = strc_selet['6'] - s +1
    return strc_selet

def splitct(inputfile,startpos,endpos,outputdir):
    startpos=int(startpos)
    endpos=int(endpos)
    title=str(inputfile)[:-3]
    strc_ric = pd.read_csv(inputfile,names=['1', '2', '3', '4', '5', '6'], sep='\t', skiprows=1)
    ric_1 = select_strc(strc_ric, startpos, endpos)
    filename=str(title)+'-'+ str(startpos) + '-' + str(endpos) + '.ct'
    ric_1.to_csv(filename, sep='\t', index=False, header=False)
    currentfile_directory = filename
    shutil.move(currentfile_directory, outputdir)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="RNAseq data analysis pipeline. First-writtern by Max.")
    parser.add_argument("-i", "--inputfile", type=str, required=True, help="Please type in the ct file.")
    parser.add_argument("-s", "--startpos", type=str, required=True, help="Please type in the start position.")
    parser.add_argument("-e", "--endpos", type=str, required=True, help="Please type in the end position.")
    parser.add_argument("-o", "--outputdir", type=str, required=True, help="Please type in the outputdir.")
    Args = parser.parse_args()
    splitct(os.path.abspath(Args.inputfile),Args.startpos,Args.endpos,os.path.abspath(Args.outputdir))