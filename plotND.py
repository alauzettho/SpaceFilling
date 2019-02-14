########################################################################################
###
###	Created by Thomas Alauzet on January 29, 2019.
###	Copyright 2019. All rights reserved.
###
###	For ND models only
###
########################################################################################

import numpy as np
import matplotlib.pyplot as plt

myfile		= open('../input.txt', 'r')
lineList	= myfile.readlines()
ndim 		= int(lineList[0])
npoint 		= int(lineList[1])
colors		= np.random.rand(npoint)
myfile.close()
myfile		= open('../output.txt', 'r')
lineList	= myfile.readlines()
M			= len(lineList)
myfile.close()


plt.figure(figsize = (10, 15), num = 'Projection')
ax			= plt.subplot(1, 1, 1)
area		= 50
pause		= 4
prec		= 1
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)


def plot(a, b) :
	X = [float(i) for i in lineList[a].strip().split('\t')]
	Y = [float(i) for i in lineList[b].strip().split('\t')]
	plt.scatter(X, Y, s = area, c = colors, alpha = prec)
	plt.xlabel("Dim " + str(ndim + a - M + 1))
	plt.ylabel("Dim " + str(ndim + b - M + 1))

	
plt.subplot(221)
plot(M - 1, M - 2)
plt.subplot(222)
plot(M - 1, M - 3)
plt.subplot(223)
plot(M - 1, M - 4)
plt.subplot(224)
plot(M - 3, M - 2)
plt.show()