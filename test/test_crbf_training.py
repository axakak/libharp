#!/usr/bin/env python3.4
import os
import signal
import sys
import argparse
import subprocess as sub
import yaml

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axes3d, art3d
from numpy import *

parser = argparse.ArgumentParser(description="To test crbf training")
parser.add_argument('-t','--train', action='store', default='', metavar='<file>', help='Train neural network using trace list <file>')
parser.add_argument('-p','--plot', action='store_true', help='Plot data from previous training')
parser.add_argument('-z','--zdata', action='store_true', help='Plot z spacial data on the z-axis')
parser.add_argument('--gif', action='store_true')
parser.add_argument('--png', action='store_true')
args = parser.parse_args()

if not (args.train or args.plot):
    print('no command given, please use -h for a list of options')
    quit()

origDir = os.getcwd()
rootDir = os.path.dirname(os.path.dirname(origDir))

if args.train:
    crbfTrain = os.path.join(rootDir,'bin/crbfTrainer')

    # check that build is up to date
    makeCMD = "make -q --directory='{}'".format(rootDir);
    print(makeCMD)
    if sub.call(makeCMD, shell=True):
        print('\x1B[33mwarning:\x1B[0m build not up to date')

    # run c-rbf training
    crbfTrainCMD = "'{}' {}".format(crbfTrain, args.train)
    print(crbfTrainCMD)
    sub.call(crbfTrainCMD,shell=True,stdout=sys.stdout)


if args.plot:
    # display results
    print('\x1B[34m==> \x1B[0m Plotting results')

    #setup 3D subplot
    x,y = plt.figaspect(.75)*1.5
    fig = plt.figure(figsize=(x,y), tight_layout=True)
    ax = fig.add_subplot(111, projection='3d')

    # set axes properties
    ax.set_title('Spatio-temporal Layer Training')
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

    print('Plotting normalized training data')
    #load training trace date yaml file
    for yamlDoc in yaml.load_all(open('normalizedTrainingData.yaml', 'r')):
        ntd = np.array((yamlDoc['events'][1:]))
        tdScatter = ax.scatter(ntd[::60,0],ntd[::60,1],ntd[::60,zidx], s=5, c='grey', edgecolor='paleturquoise', marker='o', depthshade=False, zorder=0, alpha=0.8)

    ims = []

    print('Plotting neural network and links')
    for yamlDoc in yaml.load_all(open('spatioTemporalTrain.yaml', 'r')):
        stnw = np.array((yamlDoc['spatio-temporal-neuron-weights'][1:]))
        stnc = np.array((yamlDoc['spatio-temporal-neuron-edges'][1:]))
        links = np.array([vstack((stnw[cidx[0]], stnw[cidx[1]])).tolist() for cidx in stnc])
        lc = art3d.Line3DCollection([[tuple(con[0,[0,1,zidx]]), tuple(con[1,[0,1,zidx]])] for con in links], lw=0.5 )
        lc.set_color('coral')
        lines = ax.add_collection(lc)
        nScatter = ax.scatter(stnw[:,0],stnw[:,1],stnw[:,zidx], s=30, c='crimson', edgecolor='coral', marker='o', depthshade=False)
        ims.append([lines,nScatter])


    legend = ax.legend([tdScatter, (lines, nScatter)] ,['Training Data','Spatio-temporal Neurons'], fontsize='medium', loc='lower right')
    legend.get_frame().set_edgecolor('darkgray')

    ani = animation.ArtistAnimation(fig, ims, interval=1500, repeat_delay=3000)

    if args.gif:
        ani.save('neural_training.gif', writer='imagemagick');

    if args.png:
        ani.save('neural_training.png', writer='imagemagick');


    plt.show()
