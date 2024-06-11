# -*- coding: utf-8 -*-
# """
# Created on 2023.2.10
# 将ct结构文件转换为bedpe格式，方便在热图中显示。可以同时将shapemap和vRICseq数据展示在HiC热图上。
# @author: Max
import numpy as np
file_dir = 'D:/360MoveData/Users/dgwei/Desktop/文章/Scientific_data-Sunju/result/'
# file_dir = '/store/max/'
# strc_ric = np.loadtxt( file_dir + 'PEDV.phase4_final.whole-2.ct', usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<U2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)
# strc_shape = np.loadtxt( file_dir + 'MK584552.1.incell.ct', usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<U2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)
def ct2bedpe(dirs , ct_f, out_f):
    rna_strc = np.loadtxt(dirs + ct_f, usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<U2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)
    outs = np.zeros(rna_strc.shape[0], dtype=[('c1','<U4'),('s1','<i4'),('e1','<i4'),('c2','<U4'),('s2','<i4'),('e2','<i4')])
    outs['c1'] = 'HS1'
    outs['c2'] = 'HS1'
    outs['s1'] = rna_strc['1']
    outs['s2'] = rna_strc['5']
    outs['e1'] = rna_strc['1'] + 1
    outs['e2'] = rna_strc['5'] + 1
    np.savetxt(dirs + out_f, outs[outs['s2']>0], fmt = '%s', delimiter='\t')
ct2bedpe(file_dir,"virus.phase4_final.whole.ct","virus.phase4_final.whole.bedpe")
