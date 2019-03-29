########################################################################################
###
###	Created by Thomas Alauzet on March 07, 2019.
###	Copyright 2019. All rights reserved.
###
########################################################################################


import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from   mpl_toolkits.mplot3d import Axes3D

file		= sys.argv[1]
myfile		= open(str(file), 'r')
lineList	= myfile.readlines()
ndim		= len(lineList)
npoint		= len((lineList[0].split('\t'))) - 1
myfile.close()



plt.figure(figsize = (10, 15), num = 'Space Filling : ' + str(ndim) + " | " + str(npoint))

if (ndim == 3) :
	ax		= plt.axes(projection='3d')
else :
	ax		= plt.subplot(1, 1, 1)

area	= 50
pause	= 4
prec	= 1
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)


def plot(a, b) :
	X		= [float(i) for i in lineList[a].strip().split('\t')]
	Y		= [float(i) for i in lineList[b].strip().split('\t')]
	plt.scatter(X, Y, s = area, alpha = prec)
	plt.xlabel("Dim " + str(a - 1))
	plt.ylabel("Dim " + str(b - 1))

def plop(a, b, c) :
	X		= [float(i) for i in lineList[a].strip().split('\t')]
	Y		= [float(i) for i in lineList[b].strip().split('\t')]
	Z		= [float(i) for i in lineList[c].strip().split('\t')]
	ax.scatter(X, Y, Z, s = area, alpha = prec)
	plt.xlabel("Dim " + str(a - 1))
	plt.ylabel("Dim " + str(b - 1))

if (ndim == 2) :
	plot(ndim - 2, ndim - 1)

if (ndim == 3) :
	plop(ndim - 3, ndim - 2, ndim - 1)

if (ndim >= 4) :
	plt.subplot(221)
	plot(ndim - 1, ndim - 2)
	plt.subplot(222)
	plot(ndim - 1, ndim - 3)
	plt.subplot(223)
	plot(ndim - 1, ndim - 4)
	plt.subplot(224)
	plot(ndim - 3, ndim - 2)


plt.show()