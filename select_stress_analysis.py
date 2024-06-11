# -*- coding: utf-8 -*-
# """
# Created on 2023.6.21
# 用于选择压力分析的脚本
# @author: Max

#1）提取基因组文件的id号下载注释文件
# dict={}
# with open("D:/360MoveData/Users/dgwei/Desktop/AlphaCoV.fasta","r")as F:
#     for line in F:
#         if line.startswith(">"):
#             name=line.split(" ")[0].rstrip()
#             dict[name]=""
#         else:
#             dict[name]+=line.rstrip()
# with open("D:/360MoveData/Users/dgwei/Desktop/AlphaCoV.id.txt","w")as f:
#     for key,value in dict.items():
#         print(key)
#         f.write(key[1:]+"\n")
#
#2）从不同病毒序列中提取出不同的基因
import pandas as pd
import numpy as np
data=pd.read_excel("D:/360MoveData/Users/dgwei/Desktop/6.24/AlphaCoV.genepos.xlsx","genepos",header=0,skiprows=0)
posorf1a=data.loc[data.gene=="ORF1a",];posorf1b=data.loc[data.gene=="ORF1b",];s=data.loc[data.gene=="S"];e=data.loc[data.gene=="E"];
m=data.loc[data.gene=="M"];n=data.loc[data.gene=="N"];posorf1ab=data.loc[data.gene=="ORF1ab",]
dict={}
with open("D:/360MoveData/Users/dgwei/Desktop/6.24/AlphaCoV.fasta","r")as F:
    for line in F:
        if line.startswith(">"):
            name=line.split(" ")[0].rstrip()
            dict[name]=""
        else:
            dict[name]+=line.rstrip()
with open("D:/360MoveData/Users/dgwei/Desktop/6.24/posorf1a.alphaCoV.fasta","w")as f:
    for key,value in dict.items():
        start=posorf1a[posorf1a.id==key].iloc[:,2]
        end=posorf1a[posorf1a.id==key].iloc[:,-1]
        seq=value[int(start)-1:int(end)]
        f.write(key+"\n")
        f.write(seq+"\n")
# with open("D:/360MoveData/Users/dgwei/Desktop/6.24/posorf1b.alphaCoV.fasta", "w") as f:
#     for key, value in dict.items():
#         start = posorf1b[posorf1b.id == key].iloc[:, 2]
#         end = posorf1b[posorf1b.id == key].iloc[:, -1]
#         seq = value[int(start) - 1:int(end)-3]
#         f.write(key + "\n")
#         f.write(seq+ "\n")
# with open("D:/360MoveData/Users/dgwei/Desktop/6.24/s.alphaCoV.fasta", "w") as f:
#     for key, value in dict.items():
#         start = s[s.id == key].iloc[:, 2]
#         end = s[s.id == key].iloc[:, -1]
#         seq = value[int(start) - 1:int(end)-3]
#         f.write(key + "\n")
#         f.write(seq+ "\n")
# with open("D:/360MoveData/Users/dgwei/Desktop/6.24/e.alphaCoV.fasta", "w") as f:
#     for key, value in dict.items():
#         start = e[e.id == key].iloc[:, 2]
#         end = e[e.id == key].iloc[:, -1]
#         seq = value[int(start) - 1:int(end)-3]
#         f.write(key + "\n")
#         f.write(seq+ "\n")
# with open("D:/360MoveData/Users/dgwei/Desktop/6.24/m.alphaCoV.fasta", "w") as f:
#     for key, value in dict.items():
#         start = m[m.id == key].iloc[:, 2]
#         end = m[m.id == key].iloc[:, -1]
#         seq = value[int(start) - 1:int(end)-3]
#         f.write(key + "\n")
#         f.write(seq+ "\n")
# with open("D:/360MoveData/Users/dgwei/Desktop/6.24/n.alphaCoV.fasta", "w") as f:
#     for key, value in dict.items():
#         start = n[n.id == key].iloc[:, 2]
#         end = n[n.id == key].iloc[:, -1]
#         seq = value[int(start) - 1:int(end)-3]
#         f.write(key + "\n")
#         f.write(seq+ "\n")
# with open("D:/360MoveData/Users/dgwei/Desktop/6.24/posorf1ab.alphaCoV.fasta", "w") as f:
#     for key, value in dict.items():
#         start = posorf1ab[posorf1ab.id == key].iloc[:, 2]
#         end = posorf1ab[posorf1ab.id == key].iloc[:, -1]
#         seq = value[int(start) - 1:int(end)-3]
#         f.write(key + "\n")
#         f.write(seq+ "\n")



