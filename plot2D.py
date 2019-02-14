########################################################################################
###
###	Created by Thomas Alauzet on November 27, 2018.
###	Copyright 2018. All rights reserved.
###
###	For 2D models only
###
########################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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
ax			= plt.subplot(1, 1, 1)
area		= 50
pause		= 4
prec		= 1
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)


def plot(a, b, f) :
	X			= [float(i) for i in lineList[a].strip().split('\t')]
	Y			= [float(i) for i in lineList[b].strip().split('\t')]
	red_patch	= mpatches.Patch(color = 'blue', label = 'F = %s'%(f))
	plt.legend(handles = [red_patch])
	# prec		= float(b) / float(M)
	plt.scatter(X, Y, s = area, c = colors, alpha = prec)


# for i in range(0, M) :
# 	if (lineList[i][0] == "F") :
# 		#plt.clf()
# 		plot(i + 3, i + 4, lineList[i][4:23])
# 		plt.pause(pause)

plot(12, 13, lineList[9][4:23])

plt.pause(pause)
plt.clf()

plot(M - 2, M - 1, lineList[M - 5][4:23])

plt.show()