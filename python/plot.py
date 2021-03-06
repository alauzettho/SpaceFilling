########################################################################################
###
###	Created by Thomas Alauzet on Feb 17, 2019.
###	Copyright 2019. All rights reserved.
###
########################################################################################


import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from   mpl_toolkits.mplot3d import Axes3D


myfile		= open('../output.txt', 'r')
lineList	= myfile.readlines()
M			= len(lineList)
myfile.close()


F 			= 0
ndim 		= 0
npoint 		= 0

for i in range(M) :
	if (lineList[i][0:3] == "F :") :
		F = lineList[i][4:15]
	if (lineList[i][0:3] == "Par") :
		ndim	= M - i - 1
		npoint	= len(lineList[i + 1].split("\t")) - 1


colors		= np.random.rand(npoint)

plt.figure(figsize = (10, 15), num = 'Space Filling : ' + str(ndim) + " | " + str(npoint) + " | " + str(F))

if (ndim == 3) :
	ax		= plt.axes(projection='3d')
else :
	ax		= plt.subplot(1, 1, 1)

area		= 50
pause		= 4
prec		= 1
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)


def plot(a, b) :
	X		= [float(i) for i in lineList[a].strip().split('\t')]
	Y		= [float(i) for i in lineList[b].strip().split('\t')]
	plt.scatter(X, Y, s = area, c = colors, alpha = prec)
	plt.xlabel("Dim " + str(ndim + a - M + 1))
	plt.ylabel("Dim " + str(ndim + b - M + 1))


def plop(a, b, c) :
	X		= [float(i) for i in lineList[a].strip().split('\t')]
	Y		= [float(i) for i in lineList[b].strip().split('\t')]
	Z		= [float(i) for i in lineList[c].strip().split('\t')]
	ax.scatter(X, Y, Z, s = area, c = colors, alpha = prec)
	plt.xlabel("Dim " + str(ndim + a - M + 1))
	plt.ylabel("Dim " + str(ndim + b - M + 1))


if (ndim == 2) :
	plot(M - 2, M - 1)

if (ndim == 3) :
	plop(M - 3, M - 2, M - 1)

if (ndim >= 4) :
	plt.subplot(221)
	plot(M - 1, M - 2)
	plt.subplot(222)
	plot(M - 1, M - 3)
	plt.subplot(223)
	plot(M - 1, M - 4)
	plt.subplot(224)
	plot(M - 3, M - 2)

plt.show()