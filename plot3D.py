########################################################################################
###
###	Created by Thomas Alauzet on January 29, 2019.
###	Copyright 2019. All rights reserved.
###
###	For 3D models only
###
########################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D

myfile		= open('../input.txt', 'r')
lineList	= myfile.readlines()
N 			= int(lineList[1])
colors		= np.random.rand(N)
myfile.close()
myfile		= open('../output.txt', 'r')
lineList	= myfile.readlines()
M			= len(lineList)
myfile.close()


plt.figure(figsize = (10, 15))
ax			= plt.axes(projection='3d')
area		= 50
pause		= 4
prec		= 1
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)


def plot(a, b, c, f) :
	X			= [float(i) for i in lineList[a].strip().split('\t')]
	Y			= [float(i) for i in lineList[b].strip().split('\t')]
	Z			= [float(i) for i in lineList[c].strip().split('\t')]
	red_patch	= mpatches.Patch(color = 'blue', label = 'F = %s'%(f))
	ax.legend(handles = [red_patch])
	ax.scatter(X, Y, Z, s = area, c = colors, alpha = prec)
	

plot(M - 3, M - 2, M - 1, lineList[M - 6][4:23])

plt.show()