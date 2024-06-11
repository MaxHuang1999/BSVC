# -*- coding: utf-8 -*-
# """
# Created on 2023.4.9
#Shpe-map data analysis
#Local median SHAPE reactivity and Shannon Entropy were calculated in 55nt sliding windows.
# The global median SHAPE reactivity or Shannon Entropy were subtracted from calculated values to aid in data visualization.
# Regions with local SHAPE and Shannon Entropy signals
# 1) below the global median 2) for stretches longer than 40 nucleotides 3) that appear in both replicate datasets were considered well-folded.
# Disruptions, or regions where local SHAPE or Shannon Entropy rose above the global median,
# are not considered to disqualify well-folded regions if they extended for less than 40 nucleotides
# @author: Max
# """

#从xml文件中提取shape-react值
# import xml.etree.ElementTree as ET
# count=0
# with open("D:/360MoveData/Users/dgwei/Desktop/BASIC/two-dimensional data analysis/PEDV/incell/rf_norm_rep1/MK584552.1.txt","w")as F:
#     tree = ET.parse('D:/360MoveData/Users/dgwei/Desktop/BASIC/two-dimensional data analysis/PEDV/incell/rf_norm_rep1/MK584552.1.xml')
#     root = tree.getroot()
#     for line in root:
#         a=root[0][0].text
#         b=root[0][1].text
#         for i in b.split(","):
#             ii=i.strip()
#             count+=1
#             print(ii)
#             if ii=="NaN":
#                 ii=0
#             F.write(str(count)+"\t"+str(ii)+"\n")

#利用R代码将wig文件中缺少的shannon值补齐
# setwd("D:/360MoveData/Users/dgwei/Desktop/工作/2023.4/")
# df.shannon=read.csv2("./rf_fold1/shannon/MK584552.1.shannon.csv",header=F,sep=",",stringsAsFactors = F)
# df.shannon.new=data.frame(seq(from=1,to=28044,by=1))
# colnames(df.shannon.new)="V1"
# df.shannon.newest=merge(x=df.shannon.new,y=df.shannon,by="V1",all=T)
# df.shannon.newest[is.na(df.shannon.newest)]=0
# df.shannon.newest$V2=as.numeric(df.shannon.newest$V2)

#设置窗口在PEDV基因组上滑动
# import os
# os.chdir("D:/360MoveData/Users/dgwei/Desktop/工作/2023.4/shapemapdata/AJ1102_incell/")
# win_size=25
# step_size=1
# with open("PEDVwin25","w")as f:
#     with open("PEDV.chrom.sizes","r")as F:
#         for line in F:
#             line=line.rstrip()
#             ele=line.split("\t")
#             stt=1
#             while(stt+win_size <= int(ele[1])):#得到所有窗口
#                 end=stt+win_size
#                 print(ele[0]+"\t"+str(stt)+"\t"+str(end))
#                 f.write(ele[0]+"\t"+str(stt)+"\t"+str(end)+'\n')
#                 stt=stt+step_size

# #处理shape-react
# import numpy as np
# val=[];spl_val=[]
# with open("D:/360MoveData/Users/dgwei/Desktop/SARS-CoV-2_SHAPE_MaP_structure-master/SHAPE-MaP_data/SARS-CoV-2.winReact","w")as ff:
#     with open("D:/360MoveData/Users/dgwei/Desktop/SARS-CoV-2_SHAPE_MaP_structure-master/SHAPE-MaP_data/SARS-CoV-2-MC.csv","r")as F:
#         for line in F:
#             val.append(float(line.split(',')[1].rstrip()))
#         val_array=np.array(val)
#         val_median=np.median(val_array)
#         print(val_median)
#         with open("D:/360MoveData/Users/dgwei/Desktop/工作/2023.4/shapemapdata/SARS-CoV-2-Danny/SARSCoV2win50", "r") as f:
#             for l in f:
#                 spl_val.append(str(np.median(val_array[int(l.split('\t')[1])-1:int(l.split('\t')[2])])-val_median))
#             ff.write(str('\n'.join(spl_val)))

#处理病毒粒子PEDV-shapemap数据
# import os
# os.chdir("D:/360MoveData/Users/dgwei/Desktop/工作/2023.6/6.28/")
# win_size=25
# step_size=1
# with open("PEDVwin25.txt","w")as f:
#     with open("PEDV.chrom.sizes","r")as F:
#         for line in F:
#             line=line.rstrip()
#             ele=line.split("\t")
#             stt=1
#             while(stt+win_size <= int(ele[1])):#得到所有窗口
#                 end=stt+win_size
#                 print(ele[0]+"\t"+str(stt)+"\t"+str(end))
#                 f.write(ele[0]+"\t"+str(stt)+"\t"+str(end)+'\n')
#                 stt=stt+step_size

#处理shape-react
# import numpy as np
# val=[];spl_val=[]
# with open("D:/360MoveData/Users/dgwei/Desktop/工作/2023.6/6.28/PEDV.invitro.win25React.txt","w")as ff:
#     with open("D:/360MoveData/Users/dgwei/Desktop/工作/2023.6/6.28/AJ1102-in-vitro-reactivity.txt","r")as F:
#         for line in F:
#             val.append(float(line.split('\t')[1].rstrip()))
#         val_array=np.array(val)
#         val_median=np.median(val_array)
#
#         with open("D:/360MoveData/Users/dgwei/Desktop/工作/2023.6/6.28/PEDVwin25.txt", "r") as f:
#             for l in f:
#                 spl_val.append(str(np.median(val_array[int(l.split('\t')[1])-1:int(l.split('\t')[2])])-val_median))
#             ff.write(str('\n'.join(spl_val)))



#处理shape-shannon
# import numpy as np
# val=[];spl_val=[]
# with open("D:/360MoveData/Users/dgwei/Desktop/工作/2023.6/6.28/AJ1102-in-vitro-win25shannon.txt","w")as ff:
#     with open("D:/360MoveData/Users/dgwei/Desktop/工作/2023.6/6.28/AJ1102-in-vitro-shannon.txt","r")as F:
#         for line in F:
#             val.append(float(line.rstrip()))
#         val_array=np.array(val)
#         val_median=np.median(val_array)
#         print(val_median)
#         with open("D:/360MoveData/Users/dgwei/Desktop/工作/2023.6/6.28/PEDVwin25.txt", "r") as f:
#             for l in f:
#                 spl_val.append(str(np.median(val_array[int(l.split('\t')[1])-1:int(l.split('\t')[2])])))
#             ff.write(str('\n'.join(spl_val)))


#找出winReact和winShannon中都小于globa median的位点
# pos=0;a=[];s=[]
# with open("./rep2/MK584552.1.wellfold",'w') as ff:
#     with open("./rep2/MK584552.1.winReact",'r')as F,open("./rep2/MK584552.1.winShannon",'r') as f:
#         for line1 in F:
#             line2=f.readline()
#             pos+=1
#             if float(line1) < 0.232 and float(line2) < 0.0086:
#                 a.append(pos)
#         for i in sorted(set(a)):
#             if len(s) == 0 or s[-1] + 1 == i:
#                 s.append(i)  # 入栈
#             else:
#                 if len(s) >= 2:
#                     print(s)
#                     ff.write(str(s[0]) + '\t' + str(s[-1]) + "\n")
#                 s = []  # 清空
#                 s.append(i)  # 入栈
#         # 最后一轮，需判断下
#         if len(s) >= 2:
#             print(s)

#找出rep1和rep2中共同的wellfold区域
import os
import numpy as np
os.chdir("D:/360MoveData/Users/dgwei/Desktop/工作/2023.4/shapemapdata/AJ1102_incell/")
pos=0;b=[];c=[];a=[];s=[];val1=[];val2=[];val3=[];val4=[]
with open("D:/360MoveData/Users/dgwei/Desktop/工作/PEDV全基因组功能性结构的探究/2023.4.25-ldh/MK584552.1.50.lowrect.highsh.txt",'w')as FF:
    with open("./rep1/MK584552.1.mod.winReact",'r')as f1,open("./rep1/MK584552.1.mod.winShannon",'r') as f2,\
            open("./rep2/MK584552.1.mod.winReact",'r') as f3,open("./rep2/MK584552.1.mod.winShannon",'r') as f4:
        for l1 in f1:
            l2=f2.readline();l3=f3.readline();l4=f4.readline()
            val1.append(float(l1.split('\t')[1].rstrip()));val2.append(float(l2.split('\t')[1].rstrip()))
            val3.append(float(l3.split('\t')[1].rstrip()));val4.append(float(l4.split('\t')[1].rstrip()))
        val_array1 = np.array(val1);val_array2 = np.array(val2)
        val_array3 = np.array(val3);val_array4 = np.array(val4)
        val_median1 = np.percentile(val_array1,75);val_median2 = np.percentile(val_array2,75)
        val_median3 = np.percentile(val_array3,75);val_median4 = np.percentile(val_array4,75)

        for ll1, ll2, ll3 ,ll4 in zip(val1, val2, val3, val4):
            pos += 1
            if float(ll1) < val_median1 and float(ll2) > val_median2:
                b.append(pos+25)
            if float(ll3) < val_median3 and float(ll4) > val_median4:
                c.append(pos+25)
        a=list (set(b).intersection (set(c)))

        for i in sorted(set(a)):
            if len(s) == 0 or s[-1] + 1 == i:
                s.append(i)  # 入栈
            else:
                if len(s) >= 2:
                    print(s)
                    FF.write(str(s[0]) + '\t' + str(s[-1]) + "\n")
                s = []  # 清空
                s.append(i)  # 入栈
#
# with open("./rep1/MK584552.1.mod.winReact",'r')as f1,open("./rep1/MK584552.1.mod.winShannon",'r') as f2,\
#         open("./rep2/MK584552.1.mod.winReact",'r') as f3,open("./rep2/MK584552.1.mod.winShannon",'r') as f4:
#     for l1 in f1:
#         l2=f2.readline();l3=f3.readline();l4=f4.readline()
#         val1.append(float(l1.split(',')[0].rstrip()));val2.append(float(l2.split(',')[0].rstrip())+0.0083)
#         val3.append(float(l3.split(',')[0].rstrip()));val4.append(float(l4.split(',')[0].rstrip())+0.0086)
#     val_array1 = np.array(val1);val_array2 = np.array(val2)
#     val_array3 = np.array(val3);val_array4 = np.array(val4)
#     val_median1 = np.median(val_array1);val_median2 = np.median(val_array2)
#     val_median3 = np.median(val_array3);val_median4 = np.median(val_array4)
#     print(val_median1,val_median2,val_median3,val_median4)
#     val_mean1 = np.mean(val_array1);val_mean2 = np.mean(val_array2)
#     val_mean3 = np.mean(val_array3);val_mean4 = np.mean(val_array4)
#     print(val_mean1, val_mean2, val_mean3, val_mean4)


#SARSCOV2-Danny
# import numpy as np
# pos=0;b=[];c=[];a=[];s=[];val1=[];val2=[]
# with open("SARSCoV2.wellfold.txt",'w')as FF:
#     with open("SARS-CoV-2.winReact",'r')as f1,open("SARS-CoV-2.winShannon",'r') as f2:
#         for l1 in f1:
#             l2=f2.readline()
#             val1.append(float(l1.split(',')[0].rstrip()))
#             val2.append(float(l2.split(',')[0].rstrip()))
#         val_array1 = np.array(val1)
#         val_array2 = np.array(val2)
#         val_median1 = np.median(val_array1)
#         val_median2 = np.median(val_array2)
#         print(val_median1,val_median2)
#         for ll1,ll2 in zip(val1,val2):
#             pos+=1
#             if float(ll1) < val_median1 and float(ll2) < val_median2:
#                 b.append(pos+25)
#         print(b)
#         for i in sorted(set(b)):
#             if len(s) == 0 or s[-1] + 1 == i:
#                 s.append(i)  # 入栈
#             else:
#                 if len(s) >= 2:
#                     print(s)
#                     FF.write(str(s[0]) + '\t' + str(s[-1]) + "\n")
#                 s = []  # 清空
#                 s.append(i)  # 入栈


###统计basepair和wellfold占比
#AJ1102 shapemap incell
# a=27702
# b=28044
# basepair=0
# import os
# os.chdir("D:/360MoveData/Users/dgwei/Desktop/工作/PEDV全基因组功能性结构的探究/")
# with open("AJ1102-rep1-structure.txt","r")as F:
#     for line in F:
#         seq=line[a-1:b]
#         print(seq)
#         for l in seq:
#             if l == "(" or l == ")" :
#                basepair+=1
#         print(basepair)

###统计 wellfold nulce 在不同基因中的占比
##AJ1102 shapemap incell
# import os
# os.chdir("D:/360MoveData/Users/dgwei/Desktop/工作/PEDV全基因组功能性结构的探究/")
# gene=[27702,28044]
# list=[]
# with open("wellfold.40.txt","r")as F:
#     for line in F:
#         start=line.split("\t")[0]
#         end=line.split("\t")[1].rstrip()
#         for i in range(int(start),int(end)):
#             list.append(i)
#     pos=[i for i in list if i>=gene[0] and i<=gene[1]]
#     print(len(pos))

# import os
# import numpy as np
# start=27702
# end=28044
# os.chdir("D:/360MoveData/Users/dgwei/Desktop/工作/2023.4/shapemapdata/AJ1102_incell/rep1/")
# pos=0;b=[];c=[];a=[];s=[];val1=[];val2=[]
# with open("MK584552.1.winReact",'r')as f1,open("MK584552.1.winShannon",'r') as f2:
#     for l1 in f1:
#         l2=f2.readline()
#         val1.append(float(l1))
#         val2.append(float(l2))
#     val_array1 = np.array(val1)
#     val_array2 = np.array(val2)
#     val_median1 = np.median(val_array1)
#     val_median2 = np.median(val_array2)
#     print(val_median1,val_median2)
#     for ll1,ll2 in zip(val1,val2):
#         pos+=1
#         if float(ll1) < val_median1 and float(ll2) < val_median2:
#             b.append(pos+25)
#     print(b)
# c=[i for i in b if i>=start and i<=end]
# print(len(c))


#SARSCOV2和PEDV全基因组碱基配对和wellfold距离比较
#获取文件夹中所有的文件名
# import os
# path = '/storx/zwen/shapemap_max/AJ1102-incell-SHAPE/SARACOV2/Selected_structures/'
# os.chdir(path)
# filenames = os.listdir()
#
# for i in range(len(filenames)):
#     filenames[i] = filenames[i].split('.')[0].split('_')[1]
#
# filenames = sorted(filenames) # sort the filenames in dictionary order
# print(filenames)
#SARS-COV-2-Danny
# with open('D:/360MoveData/Users/dgwei/Desktop/工作/PEDV全基因组功能性结构的探究/SARSCOV2-Danny-wellfold.txt','w')as f:
#     with open('D:/360MoveData/Users/dgwei/Desktop/工作/PEDV全基因组功能性结构的探究/SARSCOV2-Danny-wellfold.fasta','r')as F:
#         for l in F:
#             if l.startswith('>'):
#                 line=l.split('_')[1]
#                 s=line.split('-')[0];e=line.split('-')[1].rstrip()
#                 f.write(str(s)+"\t"+str(e)+"\n")

#从PEDV shapemap得到的ct结构中批量提取wellfold区域
# import numpy as np
# import os
# os.chdir("D:/360MoveData/Users/dgwei/Desktop/工作/PEDV全基因组功能性结构的探究/")
# strc_ric = np.loadtxt('D:/360MoveData/Users/dgwei/Desktop/工作/2023.4/shapemapdata/AJ1102_incell/rep1/rf_fold1_fw3000/structures/MK584552.1.ct', usecols=[0,1,2,3,4,5] , dtype=[('1','<i4'),('2','<S2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)
#
# def select_strc(ct_f, s, e, types = 'varna'):
#     strc_selet = ct_f[s-1:e].copy()
#     mask = (strc_selet['5'] < s) | (strc_selet['5'] >= e)
#     strc_selet['5'][mask] = 0
#     if types == 'rna_strc':
#         strc_selet['1'] = strc_selet['1'] - s
#         strc_selet['3'] = strc_selet['3'] - s
#         strc_selet['4'] = strc_selet['4'] - s
#         strc_selet['5'][strc_selet['5']>0] = strc_selet['5'][strc_selet['5']>0] - s
#         strc_selet['6'] = strc_selet['6'] - s
#     return strc_selet
# with open('D:/360MoveData/Users/dgwei/Desktop/工作/PEDV全基因组功能性结构的探究/wellfold.40.txt','r')as F:
#     for l in F:
#         s=l.split('\t')[0];e=l.split('\t')[1].rstrip()
#         ric_1=select_strc(strc_ric, int(s),int(e))
#         np.savetxt('4.' + str(s) + '-' + str(e) + '.ct', ric_1, fmt='%s')

# list1=[];list2=[]
# with open("D:/360MoveData/Users/dgwei/Desktop/工作/PEDV全基因组功能性结构的探究/wellfold.40.txt",'r')as f1, \
#         open("D:/360MoveData/Users/dgwei/Desktop/工作/PEDV全基因组功能性结构的探究/SARSCOV2-Danny-wellfold.txt",'r')as f2:
#     for l1 in f1:
#         l2=f2.readline()
#         ll1=l1.split('\t')[0];ll2=l1.split('\t')[1].rstrip();ll3=l2.split('\t')[0];ll4=l2.split('\t')[1].rstrip()
#         list1.append(ll1+':'+ll2)
#         list2.append(ll3+':'+ll4)
#     print(','.join(list1))
#     print(','.join(list2))

# #######################################################
# #Data 2023-4-27
# #Author Max Huang
# #E-mail 1849732267@qq.com
# #Secondary structure coloring by normalized reactivity
# #######################################################
#读入db文件并获取结构
# with open("D:/360MoveData/Users/dgwei/Desktop/工作/2023.4/shapemapdata/AJ1102_incell/rep1/rf_fold1_fw3000/structures/MK584552.1.sd.dot","w")as f:
#     with open("D:/360MoveData/Users/dgwei/Desktop/工作/2023.4/shapemapdata/AJ1102_incell/rep1/rf_fold1_fw3000/structures/MK584552.1.dot","r")as F:
#         for line in F:
#             seq=["Single-Stranded" if l=="." else "Double-Stranded" for l in line ]
#             print(seq)
#         f.write("\n".join(seq))

# -*- coding: utf-8 -*-
# """
# Created on 2023.5.6
#计算不同类型碱基在全基因组中的占比
# import os
# ll1=[];ll2=[];ll3=[];ll4=[]
# os.chdir("D:/360MoveData/Users/dgwei/Desktop/工作/2023.4/shapemapdata/AJ1102_incell/")
# with open("MK584552.1.wellfold.txt","r")as f1,open("MK584552.1.lowrect.highsh.txt","r")as f2,\
#     open("MK584552.1.highrect.lowsh.txt","r")as f3,open("MK584552.1.highrect.highsh.txt","r")as f4:
#     for l1 in f1:
#         l2=f2.readline()
#         l3=f3.readline()
#         l4=f4.readline()
#         ll1+=[i for i in range(l1.split("\t")[0],l1.split("\t")[1])]
#         ll2+=[i for i in range(l2.split("\t")[0],l2.split("\t")[1])]
#         ll3+=[i for i in range(l3.split("\t")[0],l3.split("\t")[1])]
#         ll4+=[i for i in range(l4.split("\t")[0],l4.split("\t")[1])]
#     print(len(ll2)/28044)


# #######################################################
# -*- coding: utf-8 -*-
# #Data 2023-6-13
# #Author Max Huang
# #E-mail 1849732267@qq.com
# #协变分析从well-folded区域中找出潜在的保守功能元件
# #######################################################
#1）协变分析输入格式
# import os
# os.chdir("D:/360MoveData/Users/dgwei/Desktop/SHAPE-PEDV/")
# dict={};startpot=[];endpot=[];count=0
# with open("MK584552.fasta","r")as f1,open("wellfold.40.txt","r")as f2,open("MK584552.1.dot","r")as f3:
#     for l in f1:
#         if l.startswith(">"):
#             id=l.split(" ")[0].rstrip()
#             dict[id]=""
#         else:
#             dict[id]+=l.rstrip()
#     for ll in f2:
#         startpot.append(ll.split("\t")[0].rstrip())
#         endpot.append(ll.split("\t")[1].rstrip())
# # structure1=l3[int(startpot[count])-1:int(endpot[count])]
#     for lll in f3:
#         dot=lll
#
# print(dict[id][int(startpot[0])-1:int(endpot[0])])
# for i in range(len(startpot)):
#     with open("element%s.dot"%(i),"w")as F:
#         F.write(">element"+str(i)+"\n")
#         F.write(dict[id][int(startpot[i])-1:int(endpot[i])]+"\n")
#         F.write(dot[int(startpot[i])-1:int(endpot[i])])

#2)批量生成well-folded区域的结构（ct）文件
# import os
# os.chdir("D:/360MoveData/Users/dgwei/Desktop/6.17/")
# startpot=[];endpot=[]
# with open("PEDVstr.6.17.txt","r")as f1:
#     for ll in f1:
#         startpot.append(ll.split("\t")[0].rstrip())
#         endpot.append(ll.split("\t")[1].rstrip())
# import numpy as np
# strc_ric = np.loadtxt('MK584552.1.ct', usecols=[0,1,2,3,4,5], dtype=[('1','<i4'),('2','<S2'),('3','<i4'),('4','<i4'),('5','<i4'), ('6','<i4')], skiprows=1)
# def select_strc(ct_f, s, e, types = 'varna'):
#     strc_selet = ct_f[s-1:e].copy()
#     mask = (strc_selet['5'] < s) | (strc_selet['5'] >= e)
#     strc_selet['5'][mask] = 0
#     if types == 'rna_strc':
#         strc_selet['1'] = strc_selet['1'] - s
#         strc_selet['3'] = strc_selet['3'] - s
#         strc_selet['4'] = strc_selet['4'] - s
#         strc_selet['5'][strc_selet['5']>0] = strc_selet['5'][strc_selet['5']>0] - s
#         strc_selet['6'] = strc_selet['6'] - s
#     return strc_selet
#
# for i in range(len(startpot)):
#     s=int(startpot[i])
#     e=int(endpot[i])
#     ric_1 = select_strc(strc_ric, s, e)
#     np.savetxt('4.'+str(s)+'-'+str(e)+'.ct', ric_1, fmt = '%s')
#     with open('5.'+str(s)+'-'+str(e)+'.ct', 'w') as f:
#         with open('4.'+str(s)+'-'+str(e)+'.ct','r')as F:
#             for l in F:
#                 a=l.split(' ')[1];b=list(a);c=b[2]
#                 if a.startswith('b'):
#                     d=l.replace(a,c)
#                     f.write(d)

#然后借助VARNA将切片得到的结构改为正确的结构
#3)批量使用协变分析从wellfolded区域中找出潜在的功能元件
# import General,Covariation,os
# for i in range(74):#74 elements
#     DATA_DIR='/storx/max/workspace/covariant_analysis/RNA_Covariation-main/data'
#     WORK_DIR='/storx/max/workspace/covariant_analysis/RNA_Covariation-main/output'
#     seq,dot=General.load_dot(os.path.join(DATA_DIR,"SHAPE-PEDV/element%s.dot"%(i)))['element%s'%(i)]
#     input_sto=os.path.join(WORK_DIR,"input%s.sto"%(i))
#     Covariation.dot2sto({"input":[seq,dot]},"input",input_sto, refSeq=None, GS_DE=None, mode='w')
#     out_CM = os.path.join(WORK_DIR, 'input%s.cm'%(i))
#     Covariation.cmbuild(input_sto, out_CM, verbose=False, showCMD=True)
#     h = Covariation.cmcalibrate(out_CM, use_LSF=False)
#     seqDB = os.path.join(DATA_DIR, 'CoVrefseq.fasta')
#     out_txt = os.path.join(WORK_DIR, 'cmsearch_out%s.txt'%(i))
#     out_sto = os.path.join(WORK_DIR, 'cmsearch_out%s.sto'%(i))
#     h = Covariation.cmsearch(out_CM, seqDB, out_txt, out_sto, cpu=5, toponly=True, verbose=True, showCMD=True, use_LSF=False)
#     os.system("mkdir result%s"%(i))
#     os.system("R-scape -s --cacofold --outdir result%s cmsearch_out%s.sto"%(i,i))