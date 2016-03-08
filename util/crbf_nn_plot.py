#!/usr/bin/env python3
#NOTE: the "harp.py plot" command replaces this script
import yaml

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from mpl_toolkits.mplot3d import axes3d, art3d
from numpy import *



###############################################################################
# Plotters
###############################################################################
# if 'events' in yamlDoc:

def plotTrace():
    print('Plotting trace data')

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



# if 'spatio-temporal-neuron-weights' in yamlDoc:
#     print('Plotting CRBF spatio-temporal neurons')
#
#     title = 'Spatio-temporal Neurons'
#
#     #load st-neuron weights into numpy array, skip data key (index 0) by [1:]
#     stnw = np.array((yamlDoc['spatio-temporal-neuron-weights'][1:]))
#     nScatter = ax.scatter3D(stnw[:,0],stnw[:,1],stnw[:,zidx],
#                             s=30, c='crimson', edgecolor='coral',
#                             marker='o', depthshade=False)
#
#     # print('Coloring coding weights for class: {}'.format(
#     #       yamlDoc['class-layer']['class-neurons'][0]['class-group']))
#     #
#     # cWeights = yamlDoc['class-layer']['class-neurons'][0]['weights'][1:]
#     #
#     # cSizes = [20] * len(stnw)
#     # cColor = [0] * len(stnw)
#     #
#     # for cw in cWeights:
#     #     cSizes[cw[0]] = (cw[1]*20)**3+20
#     #     cColor[cw[0]] = cw[2]
#     #
#     # nScatter = ax.scatter3D(stnw[:,0],stnw[:,1],stnw[:,zidx], s=cSizes,
#     #                         linewidth=0.1, edgecolor='gray',
#     #                         c=cColor, cmap=plt.cm.bwr, vmin=-np.pi, vmax=np.pi,
#     #                         marker='o', depthshade=False)
#     #
#     # fig.colorbar(nScatter, ax=ax, shrink=0.7)
#
# if 'spatio-temporal-neuron-edges' in yamlDoc:
#     print('Plotting CRBF spatio-temporal neuron edges')
#
#     title += ' with Edges'
#
#     #load st-neuron edges (neuron indexed) into numpy array, skip header [1:]
#     stne = np.array((yamlDoc['spatio-temporal-neuron-edges'][1:]))
#
#     #convert edges from neuron index pairs to coordinate pairs
#     edges = np.array([np.vstack((stnw[eidx[0]],
#                                  stnw[eidx[1]])).tolist() for eidx in stne])
#
#     #build 3D line collection from
#     lc = art3d.Line3DCollection([[tuple(e[0,[0,1,zidx]]),
#                                   tuple(e[1,[0,1,zidx]])] for e in edges])
#
#     #set edge properties
#     lc.set_linewidth(0.5)
#     lc.set_color('coral')
#
#     #add neuron edges to axes
#     lines = ax.add_collection3d(lc)
#
#
#
#
# # legend = ax.legend([tdScatter, (lines, nScatter)],
# #                    ['Training Data','Spatio-temporal Neurons'],
# #                    fontsize='medium', loc='lower right')
#
# #legend.get_frame().set_edgecolor('darkgray')
#
# ax.set_title(title)

def initAxis():
    ###############################################################################
    # Matplotlib Global Figure Setup
    ###############################################################################
    #setup 3D subplot
    x,y = plt.figaspect(.75)*1.5
    fig = plt.figure(figsize=(x,y), tight_layout=True)
    ax = fig.add_subplot(111, projection='3d')

    #set axes properties
    # ax.view_init(elev=40, azim=-115)#stc-3
    ax.view_init(elev=36, azim=-70)#stc-2

    #set x-axis properties
    ax.set_xlabel('x')
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])

    #set y-axis properties
    ax.set_ylabel('y')
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])

    #set z-axis properties
    if False:#args.zdata:
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

    return [ax, fig]

def td_arrayfromyaml(trace):
    skip = 200
    yamlDoc = yaml.load(open(trace))

    td = np.array((yamlDoc['events'][1::skip]))

    if yamlDoc['coordinate-space'] != 'normalized':
        tdMax = td.max(0)
        tdMin = td.min(0)
        tdScale = np.array([1,1,1,2*np.pi]) / (tdMax - tdMin)
        td = (td[:,:] - tdMin) * tdScale

    return td

if __name__ == "__main__":
    import argparse
    ###############################################################################
    # Argument Parser Setup
    ###############################################################################
    parser = argparse.ArgumentParser(description="Plot libharp output")
    #
    # parser.add_argument('-z','--zdata',
    #                     action='store_true',
    #                     help='Plot z spacial data on the z-axis')

    commandParser = argparse.ArgumentParser(description="harp plot utility")

    commandParser.add_argument('command',
                               action='store',
                               choices=['trace', 'net'],
                               help='Select plot command')

    commandParser.add_argument('args', nargs=argparse.REMAINDER)

    commandArgs = commandParser.parse_args()

    if commandArgs.command == 'trace':
        traceParser = argparse.ArgumentParser("plot trace files and lists")

        traceParser.add_argument('files',
                                nargs=2,
                                type=open,
                                help='one or more files to be plotted')

        traceArgs = traceParser.parse_args(commandArgs.args)

        [axis, figure] = initAxis()
        axis.set_title('Trace Data Set (normalized)')

        tracePlotsNon = []
        tracePlotsDom = []

        # for f in traceArgs.files:
        #     #if yaml file, plot trace
        #     #else if list, plot all traces in the list
        #     if '.txt' in f.name:
        for trace in traceArgs.files[0].read().splitlines():
            print(trace)
            td = td_arrayfromyaml(trace)

            tracePlotsNon.append(axis.scatter(td[:,0],td[:,1],td[:,3],
                                s=5, c='grey', edgecolor='blue',
                                   marker='o', depthshade=False, zorder=0,
                                   alpha=0.8))

        for trace in traceArgs.files[1].read().splitlines():
            print(trace)
            td = td_arrayfromyaml(trace)

            tracePlotsDom.append(axis.scatter(td[:,0],td[:,1],td[:,3],
                                s=5, c='grey', edgecolor='red',
                                   marker='o', depthshade=False, zorder=0,
                                   alpha=0.8))

        legend = axis.legend([tracePlotsDom[0], tracePlotsNon[0]],
                           ['dominant hand','non-dominant hand'],
                           fontsize='medium', loc='lower right')

        legend.get_frame().set_edgecolor('darkgray')

        figure.savefig('stc_3_data_set.png')

            #else print error


    # ###############################################################################
    # # YAML Document Setup
    # ###############################################################################
    # #open file into stream
    # yamlStream = open(args.crbfFilename)
    #
    # docCount = 0
    #
    # #count number of documents
    # for yamlDoc in yaml.load_all(yamlStream):
    #     docCount += 1
    #
    # print('Document count: {}'.format(docCount))
    #
    # #yaml.load* consumes the stream, must reopen
    # yamlStream = open(args.crbfFilename)
    #
    # #parse yaml stream
    # if docCount > 1:
    #     yamlDocs = yaml.load_all(yamlStream)
    # else:
    #     yamlDoc = yaml.load(yamlStream)
    #
    # #TODO: impliment doc type detection




    plt.show()
