#!/usr/bin/env python3.4
import yaml
import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from mpl_toolkits.mplot3d import axes3d, art3d
from numpy import *

###############################################################################
# Argument Parser Setup
###############################################################################
parser = argparse.ArgumentParser(description="Plot libharp output")

parser.add_argument('-z','--zdata',
                    action='store_true',
                    help='Plot z spacial data on the z-axis')

parser.add_argument('crbfFilename',
                    metavar='[filename]')

args = parser.parse_args()

###############################################################################
# YAML Document Setup
###############################################################################
#open file into stream
yamlStream = open(args.crbfFilename, 'r')

docCount = 0

#count number of documents
for yamlDoc in yaml.load_all(yamlStream):
    docCount += 1

print('Document count: {}'.format(docCount))

#yaml.load* consumes the stream, must reopen
yamlStream = open(args.crbfFilename, 'r')

#parse yaml stream
if docCount > 1:
    yamlDocs = yaml.load_all(yamlStream)
else:
    yamlDoc = yaml.load(yamlStream)

#TODO: impliment doc type detection


###############################################################################
# Matplotlib Global Figure Setup
###############################################################################
#setup 3D subplot
x,y = plt.figaspect(.65)*1.5
fig = plt.figure(figsize=(x,y), tight_layout=True)
ax = fig.add_subplot(111, projection='3d')

#set axes properties
ax.view_init(elev=30, azim=-70)

#set x-axis properties
ax.set_xlabel('x')
ax.set_xlim(0, 1)
ax.set_xticks([0, 0.5, 1])

#set y-axis properties
ax.set_ylabel('y')
ax.set_ylim(0, 1)
ax.set_yticks([0, 0.5, 1])

#set z-axis properties
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

title = ''

###############################################################################
# Plotters
###############################################################################
if 'events' in yamlDoc:
    print('Plotting trace data')

    title = 'Trace Data (Normalized)'

    #load training trace date yaml file
    td = np.array((yamlDoc['events'][1:]))

    if yamlDoc['coordinate-space'] != "normalized":
        print('Normalizing trace data')
        tdMax = td.max(0)
        tdMin = td.min(0)
        tdScale = np.array([1,1,1,2*np.pi]) / (tdMax - tdMin)
        td = (td - tdMin) * tdScale

    tdScatter = ax.scatter(td[::25,0],td[::25,1],td[::25,zidx], s=5,
                           marker='o', c='grey', edgecolor='paleturquoise',
                           depthshade=False, zorder=0, alpha=0.8)


if 'spatio-temporal-neuron-weights' in yamlDoc:
    print('Plotting CRBF spatio-temporal neurons')

    title = 'Spatio-temporal Neurons'

    #load st-neuron weights into numpy array, skip data key (index 0) by [1:]
    stnw = np.array((yamlDoc['spatio-temporal-neuron-weights'][1:]))

    print('Coloring coding weights for class {}'.format(
          yamlDoc['class-layer']['class-neurons'][0]['class-group']))

    cWeights = yamlDoc['class-layer']['class-neurons'][0]['weights'][1:]

    cSizes = [20] * len(stnw)
    cColor = [0] * len(stnw)

    for cw in cWeights:
        cSizes[cw[0]] = (cw[1]*20)**3+20
        cColor[cw[0]] = cw[2]

    nScatter = ax.scatter3D(stnw[:,0],stnw[:,1],stnw[:,zidx], s=cSizes,
                            linewidth=0.1, edgecolor='gray',
                            c=cColor, cmap=plt.cm.bwr, vmin=-np.pi, vmax=np.pi,
                            marker='o', depthshade=False)

    fig.colorbar(nScatter, ax=ax, shrink=0.7)

if 'spatio-temporal-neuron-edges' in yamlDoc:
    print('Plotting CRBF spatio-temporal neuron edges')

    title += ' with Edges'

    #load st-neuron edges (neuron indexed) into numpy array, skip header [1:]
    stne = np.array((yamlDoc['spatio-temporal-neuron-edges'][1:]))

    #convert edges from neuron index pairs to coordinate pairs
    edges = np.array([vstack((stnw[eidx[0]],
                              stnw[eidx[1]])).tolist() for eidx in stne])

    #build 3D line collection from
    lc = art3d.Line3DCollection([[tuple(e[0,[0,1,zidx]]),
                                  tuple(e[1,[0,1,zidx]])] for e in edges])

    #set edge properties
    lc.set_linewidth(0.5)
    lc.set_color('silver')

    #add neuron edges to axes
    lines = ax.add_collection(lc)




# legend = ax.legend([tdScatter, (lines, nScatter)],
#                    ['Training Data','Spatio-temporal Neurons'],
#                    fontsize='medium', loc='lower right')

#legend.get_frame().set_edgecolor('darkgray')

ax.set_title(title)

plt.show()
