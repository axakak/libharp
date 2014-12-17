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
parser.add_argument('-t','--train', action='store_true', help='Train neural network without plotting results')
parser.add_argument('-p','--plot', action='store_true', help='Plot data from previous training')
args = parser.parse_args()

rootDir = os.path.dirname(os.getcwd())
os.chdir(os.path.join(rootDir, 'test/stc_1_sample'))

if not args.plot:
    crbfTrain = os.path.join(rootDir,'bin/crbfTrainer')
    traceDataList = 'stc_1_list.txt'

    # check that build is up to date
    makeCMD = "make -q -C '{}'".format(rootDir);
    print(makeCMD)
    if sub.call(makeCMD, shell=True):
        print('\x1B[33mwarning:\x1B[0m build not up to date')

    # run c-rbf training
    crbfTrainCMD = "'{}' {}".format(crbfTrain, traceDataList)
    print(crbfTrainCMD)
    sub.call(crbfTrainCMD,shell=True,stdout=sys.stdout)


if not args.train:
    # display results
    print('\x1B[34m==> \x1B[0m Plotting results')

    #setup 3D subplot
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # set axes properties
    ax.set_title('Spatio-temporal Layer Training')

    ax.set_xlabel('x')
    ax.set_xlim3d(0, 1)

    ax.set_ylabel('y')
    ax.set_ylim3d(0, 1)

    zidx = 2

    if zidx == 3:
        ax.set_zlabel('time')
        ax.set_zlim3d(0, 2*np.pi)
    elif zidx == 2:
        ax.set_zlabel('z')
        ax.set_zlim3d(0, 1)


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


    ax.legend([tdScatter, (lines, nScatter)] ,['Training Data','Spatio-temporal Neurons'])
    
    ani = animation.ArtistAnimation(fig, ims, interval=1500, repeat_delay=3000)

    #ani.save('neural_training.mp4')

    plt.show()
