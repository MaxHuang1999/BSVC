#######################################################
#Data 2023-10-15
#Author Max Huang
#E-mail 1849732267@qq.com
#RNA three-dimensional structure
#######################################################
# import numpy as np
# import os
#
# os.chdir("D:/360MoveData/Users/dgwei/Desktop/")
# examplematrix=np.load("counts.npy")
# print(examplematrix)
# # print(examplematrix.shape)
#
# from mayavi import mlab
#
# # Using numpy array may be better
# x = [1, 2, 3, 4, 5]
# y = [1, 2, 3, 4, 5]
# z = [1, 2, 3, 4, 5]
#
# s = [1.5, 4.7, -0.112, -3]
# def f(x, y, z):
# 	return x + y - z
#
# mlab.points3d(x, y, z)
# mlab.points3d(x, y, z, s)
# mlab.points3d(x, y, z, f)
# mlab.show()

####################################################################
#将matrix改为矩阵
import os
import numpy as np

os.chdir("D:/360MoveData/Users/dgwei/Desktop/")
data=[]
with open("./工作/Vis_of_Covdata/SARS-cov2-all-rep-rna.1nt.none.matrix","r")as file:
    for line in file:
        line=line.rstrip("\n")
        values=line.split(" ")
        start=int(float(values[0]))
        end=int(float(values[1]))
        value=int(float(values[2]))
        data.append([start,end,value])

max_start=max(row[0] for row in data)
max_end=max(row[1] for row in data)
max_size=(len(data),max_end-max_start+1)
print(max_size,max_end,max_end)
matrix_size=(max_end+1,max_end+1)
matrix=np.zeros(matrix_size)

for i,row in enumerate(data):
    start,end,value=row
    matrix[start,end]=value

np.save('SARSCoV2.1nt.npy',matrix)

