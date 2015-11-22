#!/usr/bin/env python3


def find_bin(name):
    import os

    cwdDir = os.getcwd()

    # get libharp position in path
    i = cwdDir.find('libharp')

    if i:
        libDir = os.path.join(cwdDir[:i],'libharp')
    else:
        print('Error: utility only works when run from inside project' )

    binPath = os.path.join(libDir,'build','util', name)

    if os.path.isfile(binPath):
        return binPath
    else:
        print('Error: unable to locate bin file: {}'.format(name))
        quit()

###############################################################################
# Train Neural Network
###############################################################################
def train(in_file, out_file='crbfNeuralNet.yaml'):
    import subprocess as sub
    from datetime import datetime

    harpTrain = find_bin('harptrain')

    out = ["'{}' {} {}\n".format(harpTrain, in_file, out_file)]

    # run c-rbf training
    harpTrainArgs = [harpTrain, in_file, out_file]
    print(harpTrainArgs)
    p = sub.Popen(harpTrainArgs, stdout=sub.PIPE, stderr=sub.STDOUT, bufsize=1)

    for line in iter(p.stdout.readline, b''):
        so = line.decode()
        out.append(so)
        if '==>' in so:
            so = so.replace('\n','') + datetime.today().strftime(" [%H:%M:%S]\n")

        print(so, end='')

    return out

###############################################################################
# Evaluate Neural Network
###############################################################################
def eval(neural_net, trace):
    import subprocess as sub

    harpEval = find_bin('harpevaluate')

    out = ["'{}' {} {}\n".format(harpEval, neural_net, trace)]

    # run c-rbf evaluation
    harpEvalArgs = [harpEval, neural_net, trace]
    print(harpEvalArgs)
    p = sub.Popen(harpEvalArgs, stdout=sub.PIPE, stderr=sub.STDOUT, bufsize=1)

    for line in iter(p.stdout.readline, b''):
        so = line.decode()
        out.append(so)
        print(so, end='')

    return out


###############################################################################
# Plot Training output
###############################################################################
def plot(train_list, st_neuron_train, zdata, gif, png, show):
    import yaml
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from mpl_toolkits.mplot3d import axes3d, art3d

    # display results
    print('\x1B[34m==> \x1B[0m Plotting ')

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

    if zdata:
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

    print('Plotting normalized training data in: {}'.format(train_list.name))
    #load training trace date yaml file
    trainList = train_list.read().splitlines()
    skip = 250
    for trace in trainList:
        print(trace)
        yamlDoc = yaml.load(open(trace))
        td = np.array((yamlDoc['events'][1:]))
        if yamlDoc['coordinate-space'] != 'normalized':
            tdMax = td.max(0)
            tdMin = td.min(0)
            tdScale = np.array([1,1,1,2*np.pi]) / (tdMax - tdMin)
            td = (td[::skip,:] - tdMin) * tdScale
        tdScatter = ax.scatter(td[:,0],td[:,1],td[:,zidx],
                               s=5, c='grey', edgecolor='paleturquoise',
                               marker='o', depthshade=False, zorder=0,
                               alpha=0.8)

    ims = []

    print('Plotting neural network and links in : {}'.format(st_neuron_train.name))
    for yamlDoc in yaml.load_all(st_neuron_train):
        stnw = np.array((yamlDoc['spatio-temporal-neuron-weights'][1:]))
        stnc = np.array((yamlDoc['spatio-temporal-neuron-edges'][1:]))
        links = np.array([np.vstack((stnw[cidx[0]], stnw[cidx[1]])).tolist()
                                 for cidx in stnc])

        lc = art3d.Line3DCollection([[tuple(con[0,[0,1,zidx]]),
                                      tuple(con[1,[0,1,zidx]])] for con in links], lw=0.5)

        lc.set_color('coral')
        lines = ax.add_collection(lc)
        nScatter = ax.scatter(stnw[:,0],stnw[:,1],stnw[:,zidx],
                              s=30, c='crimson', edgecolor='coral', marker='o',
                              depthshade=False)

        ims.append([lines,nScatter])


    legend = ax.legend([tdScatter, (lines, nScatter)],
                       ['Training Data','Spatio-temporal Neurons'],
                       fontsize='medium', loc='lower right')

    legend.get_frame().set_edgecolor('darkgray')

    ani = animation.ArtistAnimation(fig, ims, interval=1500, repeat_delay=3000)

    if gif:
        ani.save('neural_training.gif', writer='imagemagick');

    if png:
        ani.save('neural_training.png', writer='imagemagick');

    if show:
        print('Showing plot')
        plt.show()


if __name__ == "__main__":
    import argparse
    import time

    commandParser = argparse.ArgumentParser(description="harp command line utility")

    commandParser.add_argument('command',
                               action='store',
                               choices=['train', 'eval', 'plot'],
                               help='Select operating mode')

    commandParser.add_argument('args', nargs=argparse.REMAINDER)

    commandArgs = commandParser.parse_args()

    if commandArgs.command == 'train':
        trainParser = argparse.ArgumentParser("Train neural network")

        trainParser.add_argument('trace_list_file',
                                 help='list of trace files to train neural network')

        trainParser.add_argument('out_file',
                                 nargs='?',
                                 action='store',
                                 default='crbfNeuralNet.yaml',
                                 help='Export trained nueral network to <ofile>')

        trainArgs = trainParser.parse_args(commandArgs.args)

        sTime = time.time()
        train(trainArgs.trace_list_file, trainArgs.out_file)
        eTime = time.time()
        dTime = eTime-sTime
        print('Training completed in {:n}m{:.0f}s'.format(dTime//60, dTime%60))


    elif commandArgs.command == 'eval':
        evalParser = argparse.ArgumentParser("Evaluate trace")

        evalParser.add_argument('neural_net_file',
                                help='trained neural network')

        evalParser.add_argument('trace_file',
                                help='trace to be evaluated')

        evalArgs = evalParser.parse_args(commandArgs.args)

        sTime = time.time()
        eval(evalArgs.neural_net_file, evalArgs.trace_file)
        eTime = time.time()
        dTime = eTime-sTime
        print('Evaluation completed in {:n}m{:.0f}s'.format(dTime//60, dTime%60))


    elif commandArgs.command == 'plot':
        plotParser = argparse.ArgumentParser("visualize the input and output of harp")

        plotParser.add_argument('train_list',
                                type=open)

        plotParser.add_argument('st_neuron_train',
                                type=open)

        plotParser.add_argument('-z','--zdata',
                                action='store_true',
                                help='Plot z spacial data on the z-axis')

        plotParser.add_argument('--gif',
                                action='store_true')

        plotParser.add_argument('--png',
                                action='store_true')

        plotParser.add_argument('--show',
                                action='store_true')

        plotArgs = plotParser.parse_args(commandArgs.args)

        plot(plotArgs.train_list, plotArgs.st_neuron_train,
             plotArgs.zdata, plotArgs.gif, plotArgs.png, plotArgs.show)
