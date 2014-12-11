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
    #os.chdir(rootDir)
    if sub.call(["make","-q","-C {}".format(rootDir)]):
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

    # Setting the axes properties
    ax.set_xlabel('x')
    ax.set_xlim3d(0, 1)

    ax.set_ylabel('y')
    ax.set_ylim3d(0, 1)

    ax.set_zlabel('time')
    ax.set_zlim3d(0, 2*np.pi)

    print('Plotting normalized training data')
    #load training trace date yaml file
    for yamlDoc in yaml.load_all(open('normalizedTrainingData.yaml', 'r')):
        ntd = np.array((yamlDoc['events'][1:]))
        ax.scatter(ntd[::100,0],ntd[::100,1],ntd[::100,3], s=5, c='grey', edgecolor='paleturquoise', marker='o', depthshade=False, zorder=0, label='Training Data', alpha=0.8)

    ims = []

    print('Plotting neural network and links')
    for yamlDoc in yaml.load_all(open('spatioTemporalTrain.yaml', 'r')):
        stnw = np.array((yamlDoc['spatio-temporal-neuron-weights'][1:]))
        stnc = np.array((yamlDoc['spatio-temporal-neuron-edges'][1:]))
        links = np.array([vstack((stnw[cidx[0]], stnw[cidx[1]])).tolist() for cidx in stnc]) 
        lc = art3d.Line3DCollection([[tuple(con[0,[0,1,3]]), tuple(con[1,[0,1,3]])] for con in links], lw=0.5 )
        lc.set_color('coral')
        lines = ax.add_collection(lc)
        scatter = ax.scatter(stnw[:,0],stnw[:,1],stnw[:,3], s=25, c='crimson', edgecolor='coral', marker='o', depthshade=False)
        ims.append([lines,scatter])

    ani = animation.ArtistAnimation(fig, ims, interval=1000, repeat_delay=2000)

    #ani.save('neural_training.mp4')

    plt.show()
