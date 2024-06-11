# -*- coding: utf-8 -*-
# """
# Created on 2023.2.26
#Covariation analysis
# @author: Max

import os
os.system("export PYTHONPATH=/store/max/IPyRSSA:$PYTHONPATH")
import General,Covariation
os.system("mafft --auto ./data/SARSCOV2seq.fasta > ./data/SARSCOV2seq.aln.fasta")
data="/store/max/RNA_Covariation-main/data"
workdir="/store/max/RNA_Covariation-main/data"## 这里设置你自己的数据目录和工作目录
ORF10_seq,ORF10_dot = General.load_dot(os.path.join(data, "SARSCOV2_ORF10.dot"))['ORF10']# 读入数据
input_sto = os.path.join(workdir, 'input.sto')
Covariation.dot2sto({'input':[ORF10_seq,ORF10_dot]}, "input", input_sto, refSeq=None, GS_DE=None, mode='w')# 把序列和二级结构存储为stockholm文件
#### cmbuild
out_CM = os.path.join(workdir, 'input.cm')
Covariation.cmbuild(input_sto, out_CM, verbose=False, showCMD=True)
#### cmcalibrate
h = Covariation.cmcalibrate(out_CM, use_LSF=False, LSF_parameters={'cpu': 5})
#### 4.3 cmsearch
seqDB = os.path.join(data, 'SARSCOV2seq.fasta')
out_txt = os.path.join(workdir, 'cmsearch_out.txt')
out_sto = os.path.join(workdir, 'cmsearch_out.sto')
h = Covariation.cmsearch(out_CM, seqDB, out_txt, out_sto, cpu=5, toponly=True, verbose=True, showCMD=True, use_LSF=False, LSF_parameters={'cpu': 5})

# #### cmalign
# filtered_fasta = os.path.join(workdir, 'filtered.fasta')
# out_sto = os.path.join(workdir, 'cmalign_out.txt')
# h = Covariation.cmalign(input_cm, filtered_fasta, out_sto, cpu=5, verbose=True, showCMD=True, use_LSF=False, LSF_parameters={'cpu': 5})
# h.wait()
