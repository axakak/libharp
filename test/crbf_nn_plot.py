#!/usr/bin/env python3.4
import yaml
import argparse

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, art3d
from numpy import *

parser = argparse.ArgumentParser(description="To visualize trained crbf neural network")
parser.add_argument('-z','--zdata', action='store_true', help='Plot z spacial data on the z-axis')
parser.add_argument('crbfFilename', metavar='[filename]')
args = parser.parse_args()

#setup 3D subplot
x,y = plt.figaspect(.65)*1.5
fig = plt.figure(figsize=(x,y), tight_layout=True)
ax = fig.add_subplot(111, projection='3d')

# set axes properties
ax.set_title('Spatio-temporal Layer')
ax.view_init(elev=30, azim=-70)

ax.set_xlabel('x')
ax.set_xlim(0, 1)
ax.set_xticks([0, 0.5, 1])

ax.set_ylabel('y')
ax.set_ylim(0, 1)
ax.set_yticks([0, 0.5, 1])

if args.zdata:
    zidx = 2
    ax.set_zlabel('z')
    ax.set_zlim3d(0, 1)
    ax.set_zticks([0, 0.5, 1])
else:
    zidx = 3
    ax.set_zlabel('time')
    ax.set_zlim(0, 2*np.pi)
    ax.set_zticks([0, np.pi, 2*np.pi])
    ax.set_zticklabels(['0', '$\pi$','2$\pi$'])

print('Plotting CRBF spatio-temporal neurons')

#open file into stream
yamlSream = open(args.crbfFilename, 'r')

#parse yaml stream
yamlDoc = yaml.load(yamlSream)

#load st-neuron weights into numpy array, skip data key (index 0) by [1:]
stnw = np.array((yamlDoc['spatio-temporal-neuron-weights'][1:]))

#load st-neuron edges (neuron indexed) into numpy array, skip data key by [1:]
stne = np.array((yamlDoc['spatio-temporal-neuron-edges'][1:]))

#convert edges from neuron index pairs to coordinate pairs
edges = np.array([vstack((stnw[eidx[0]], stnw[eidx[1]])).tolist() for eidx in stne])

#build 3D line collection from
lc = art3d.Line3DCollection([[tuple(e[0,[0,1,zidx]]), tuple(e[1,[0,1,zidx]])] for e in edges])

#set edge properties
lc.set_linewidth(0.5)
lc.set_color('silver')

#add neuron edges to axes
lines = ax.add_collection(lc)

cWeights = yamlDoc['class-layer']['class-neurons'][0]['weights'][1:]


if (len(stnw) != len(cWeights)):
    print('\x1B[31merror\x1B[0m: neuron count mismatch ({}:{})'.format(len(stnw),len(cWeights)))
    exit()

cSizes = [0] * len(cWeights)
cColor = [0] * len(cWeights)

for cw in cWeights:
    cSizes[cw[0]] = (cw[1]*20)**3+20
    cColor[cw[0]] = cw[2]

nScatter = ax.scatter3D(stnw[:,0],stnw[:,1],stnw[:,zidx], s=cSizes,
                        linewidth=0.1, edgecolor='gray',
                        c=cColor, cmap=plt.cm.bwr, vmin=-np.pi, vmax=np.pi,
                        marker='o', depthshade=False)

fig.colorbar(nScatter, ax=ax, shrink=0.7)

#legend = ax.legend([tdScatter, (lines, nScatter)] ,['Training Data','Spatio-temporal Neurons'], fontsize='medium', loc='lower right')
#legend.get_frame().set_edgecolor('darkgray')

plt.show()
