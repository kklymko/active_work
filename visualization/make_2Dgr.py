#!/usr/bin/python

from numpy.random import uniform, seed
from matplotlib.mlab import griddata
import matplotlib.pyplot as plt
import numpy as np
# make up data.
#npts = int(raw_input('enter # of random points to plot:'))
plt.figure(figsize=(6,5))
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.tick_params(labelsize=20, pad = 10, length=8, width=2,labeltop='off',labelright='off')

#X_dat,Y_dat,Z_dat = np.loadtxt('corr_dat_dil', usecols=(0,1,3),unpack=True)
X_dat,Y_dat,Z_dat = np.loadtxt('coefficient/CLONEScorr_dat_dil_2000_200', usecols=(0,1,3),unpack=True)
#XX,YY = np.loadtxt('Data_test',usecols=(1,0),unpack=True)
# Convert from pandas dataframes to numpy arrays
# Convert from pandas dataframes to numpy arrays
x, y, z, = np.array([]), np.array([]), np.array([])
for i in range(len(X_dat)):
        x = np.append(x,X_dat[i])
        y = np.append(y,Y_dat[i])
        z = np.append(z,400*Z_dat[i])
# define grid.
xi = np.linspace(x.min()-.1, x.max()+.001, 500)
yi = np.linspace(y.min()-1, y.max()+.001, 500)
# grid the data.
zi = griddata(x, y, z, xi, yi, interp='linear')
# contour the gridded data, plotting dots at the nonuniform data points.
CS = plt.pcolor(xi, yi, zi, cmap='jet')
#plt.plot(XX,YY,linewidth=2, color = 'w',linestyle='dashed')
#CS = plt.contourf(xi, yi, zi, 15,
                  #vmax=abs(zi).max(), vmin=-abs(zi).max())
cb = plt.colorbar(ticks=[0.0, 1.0, 2.0, 3.0, 4.0])  # draw colorbar
# plot data points.
#plt.xlim(x.min()+.001, x.max()+.001)
plt.ylim(y.min(), y.max())
cb.ax.tick_params(labelsize=20)

plt.xlim(0.8,5)
plt.xticks(np.arange(0,5,1))
#plt.yticks(np.arange(0,25.01,5))
#plt.xticks(np.arange(min(x), max(x)+.001, 0.2))
#plt.title('griddata test (%d points)' % npts)
plt.ylabel(r'$\phi$', fontsize=25, rotation='vertical')
plt.xlabel(r'$r/\sigma$', fontsize = 25, rotation='horizontal')

plt.tight_layout()
plt.savefig('clones.pdf')

plt.show()
