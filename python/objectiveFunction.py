from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib import cm
import numpy as np


fig = plt.figure("OBJECTIVE FUNCTION")
ax  = fig.gca(projection = '3d')


mq = 0.5
x1 = [0, 0]
x2 = [1, 0]
x3 = [0, 1]
x4 = [1, 1]


def monteCarlo(param) :
    m_npoint, m_ndim  = np.shape(param)
    coeff1 = - 72 * pow(m_npoint, 2.0 / (m_ndim + 4)) / m_ndim

    fvalue = 0.0
    for i in range(m_npoint) :

        coeff2 = 0.0
        for j in range(i, m_npoint) :
            norm = np.linalg.norm(param[i] - param[j])
            coeff2 += np.exp(coeff1 * norm * norm)

        fvalue += pow(coeff2, mq - 1)

    if (mq < 1) :
        fvalue = -fvalue

    return(fvalue)


def nearestNeighbor(param) :
    m_npoint, m_ndim  = np.shape(param)
    coeff1 = m_ndim * (1 - mq)

    fvalue = 0.0
    for i in range(m_npoint) :

        dist = 1000000
        for j in range(m_npoint) :
            if (i != j) :
                norm = np.linalg.norm(param[i] - param[j])
                if (norm < dist) :
                    dist = norm

        fvalue += pow(dist, coeff1)

    if (mq < 1) :
        fvalue = -fvalue

    return(fvalue)




X = Y = np.arange(0, 1, 0.01)
X, Y  = np.meshgrid(X, Y)
a, b  = np.shape(X)
Z     = np.zeros([a, b])

for i in range(a) :
    for j in range(b) :
        x5      = [X[i, j], Y[i, j]]
        Z[i, j] = monteCarlo(np.array([x1, x2, x3, [0.5, Y[i, j]], [X[i, j], 1]]))
        # Z[i, j] = monteCarlo(np.array([x1, x2, x3, x4, x5]))



surf = ax.plot_surface(X, Y, Z, cmap = cm.coolwarm, linewidth = 0.7)
ax.set(xlabel = r'$X_1^5$', ylabel = r'$X_2^4$')
fig.colorbar(surf, shrink = 0.5, aspect = 5)
ax.set_yticklabels([])
ax.set_zticklabels([])
ax.set_xticklabels([])

plt.savefig('monteCarlo_0,5.svg')