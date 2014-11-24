#!/usr/bin/env python3.4
import os
import sys
import subprocess as sub
import yaml

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import axis3d

rootDir = os.path.dirname(os.getcwd())

crbfTrain = os.path.join(rootDir,'bin/crbfTrainer')
traceDataList = 'stc_1_list.txt'

# check that build is up to date
os.chdir(rootDir)
if sub.call(["make","-q"]):
    print('\x1B[33mwarning:\x1B[0m build not up to date')


# run c-rbf training
os.chdir(os.path.join(rootDir, 'test/stc_1_sample'))
crbfTrainCMD = "'{}' {}".format(crbfTrain, traceDataList)
print(crbfTrainCMD)
sub.call(crbfTrainCMD,shell=True,stdout=sys.stdout)


# display results
print('\x1B[34m==> \x1B[0m Plotting results')

fig2 = plt.figure()
ax2 = fig2.add_subplot(111, projection='3d')

ax2.set_xlabel('x')
ax2.set_xlim3d(0, 1)
ax2.set_ylabel('y')
ax2.set_ylim3d(0, 1)
ax2.set_zlabel('time')
ax2.set_zlim3d(0, 2*np.pi)

#load training trace date yaml file
for yamlDoc in yaml.load_all(open('normalizedTrainingData.yaml', 'r')):
    ntd = np.array((yamlDoc['events'][1:]))
    ax2.scatter(ntd[::100,0],ntd[::100,1],ntd[::100,3], c='g', edgecolor='c', marker='.', depthshade=False)

ims = []
for yamlDoc in yaml.load_all(open('spatioTemporalTrain.yaml', 'r')):
    stn = np.array((yamlDoc['spatio-temporal-neurons'][1:]))
    scatter = ax2.scatter(stn[:,0],stn[:,1],stn[:,3], c='r', edgecolor='m', marker='o', depthshade=False)
    ims.append([scatter])


ani = animation.ArtistAnimation(fig2, ims, interval=1000, repeat_delay=2000)

ani.save('neural_training.mp4')


plt.show()
