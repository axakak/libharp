#!/usr/bin/env python3


def find_bin(name):
    import os

    cwdDir = os.getcwd()

    # get libharp position in path
    i = cwdDir.find('libharp')

    if i:
        libDir = os.path.join(cwdDir[:i], 'libharp')
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
# Plot HARP data
###############################################################################
def td_arrayfromyaml(trace):
    import yaml
    import numpy as np

    skip = 200
    yamlDoc = yaml.load(open(trace))

    td = np.array((yamlDoc['events'][1::skip]))

    if yamlDoc['coordinate-space'] != 'normalized':
        tdMax = td.max(0)
        tdMin = td.min(0)
        tdScale = np.array([1,1,1,2*np.pi]) / (tdMax - tdMin)
        td = (td[:,:] - tdMin) * tdScale

    return td


def plot(harp_files, zdata, gif, png, show):
    import yaml
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    from mpl_toolkits.mplot3d import axes3d, art3d

    plt.rc('font', family = 'serif', serif = 'CMU Serif')

    trace_files = []
    trace_file_lists = []
    neural_net_files = []
    train_file = []

    print('\x1B[34m==> \x1B[0m Parsing Files')

    #TODO: if yaml files use an explicit file type key, parsing would be much simpler
    for fname in harp_files:
        print(fname)
        f = open(fname)
        if f.readline().startswith('%YAML'):# if YAML file
            docCount = f.read().count('---')#count number of documents
            f.seek(0) #return to start of stream
            if docCount == 1:
                yamlDoc = f.read(1000)
                if 'events' in yamlDoc:
                    trace_files.append(fname)
                elif 'spatio-temporal-neuron-weights' in yamlDoc:# crbf neural net
                    neural_net_files.append(fname)
            else: #multi doc yaml file, assume st training
                train_file.append(fname)
        else: # assume trace list text file
            trace_file_lists.append(fname)# add to trace data file list

    # display results
    print('\x1B[34m==> \x1B[0m Plotting')

    # vtk plot #################################################################

    import vtk


    if len(trace_files) == 1:

        td = td_arrayfromyaml(trace_files[0])

        view = vtk.vtkContextView()
        view.GetRenderer().SetBackground(1.0, 1.0, 1.0)
        view.GetRenderWindow().SetSize(400, 300)

        chart = vtk.vtkChartXYZ()
        view.GetScene().AddItem(chart)
        chart.SetShowLegend(True)
#TODO: construct vtkPolyData()

        points = chart.AddPlot(vtk.vtkChart.POINTS)
        points.SetInputData(table, 0, 1)
        points.SetColor(0, 0, 0, 255)
        points.SetWidth(1.0)
        points.SetMarkerStyle(vtk.vtkPlotPoints.CROSS)

        points = chart.AddPlot(vtk.vtkChart.POINTS)
        points.SetInputData(table, 0, 2)
        points.SetColor(0, 0, 0, 255)
        points.SetWidth(1.0)
        points.SetMarkerStyle(vtk.vtkPlotPoints.PLUS)

        points = chart.AddPlot(vtk.vtkChart.POINTS)
        points.SetInputData(table, 0, 3)
        points.SetColor(0, 0, 255, 255)
        points.SetWidth(1.0)
        points.SetMarkerStyle(vtk.vtkPlotPoints.CIRCLE)

        view.GetRenderWindow().SetMultiSamples(0)
        view.GetInteractor().Initialize()
        view.GetInteractor().Start()


    #setup 3D subplot
    x,y = plt.figaspect(.75)*1.5
    fig = plt.figure(figsize=(x,y), tight_layout=True)
    ax = fig.add_subplot(111, projection='3d')

    # set axes properties
    # ax.view_init(elev=40, azim=-115)#stc-3
    ax.view_init(elev=36, azim=-70)#stc-2

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

    leg_art = []
    leg_lable = []

    if len(trace_files) > 0:
        print('normalized trace data')
        #load training trace date yaml file
        for trace in trace_files:
            print(trace)
            td = td_arrayfromyaml(trace)

            tdScatter = ax.scatter(td[:,0],td[:,1],td[:,zidx],
                                   s=10, c='grey', edgecolor='black',
                                   marker='o', depthshade=False, zorder=0,
                                   alpha=0.8)

        leg_art.append(tdScatter)
        leg_lable.append('Trace Data')

    # for trace in f.read().splitlines():
    if len(trace_file_lists) > 0:
        print('normalized trace data groups')
        #load training trace date yaml file
        for trace_list in trace_file_lists:
            print(trace_list)

            if 'dom' in trace_list:
                eColor = 'red'
                label = 'dominant hand trace data'
            elif 'non' in trace_list:
                eColor = 'blue'
                label = 'non-dominant hand trace data'
            else:
                eColor = 'green'
                label = trace_list + ' trace data'

            for trace in open(trace_list).read().splitlines():
                print(' ',trace)
                td = td_arrayfromyaml(trace)

                tdScatter = ax.scatter(td[:,0],td[:,1],td[:,zidx],
                                       s=10, c='grey', edgecolor=eColor,
                                       marker='o', depthshade=False, zorder=0,
                                       alpha=0.8)

            leg_art.append(tdScatter)
            leg_lable.append(label)

    frames = []

    if len(train_file) > 0:
        print('Plotting neural network and links in...')
        print(train_file[0])
        for yamlDoc in yaml.load_all(open(train_file[0])):
            stnw = np.array((yamlDoc['spatio-temporal-neuron-weights'][1:]))
            stnc = np.array((yamlDoc['spatio-temporal-neuron-edges'][1:]))
            links = np.array([np.vstack((stnw[cidx[0]], stnw[cidx[1]])).tolist()
                                     for cidx in stnc])

            lc = art3d.Line3DCollection([[tuple(con[0,[0,1,zidx]]),
                                          tuple(con[1,[0,1,zidx]])] for con in links], lw=0.5)

            lc.set_color('coral')
            ntLines = ax.add_collection(lc)
            ntScatter = ax.scatter(stnw[:,0],stnw[:,1],stnw[:,zidx],
                                  s=30, c='crimson', edgecolor='coral', marker='o',
                                  depthshade=False)

            frames.append([ntLines,ntScatter])

        leg_art.append((ntLines,ntScatter))
        leg_lable.append('Spatio-temporal Neurons')

    elif len(neural_net_files) > 0:
        print('Plotting neural network and links...')
        for neural_net in neural_net_files:
            print(neural_net)
            yamlDoc = yaml.load(open(neural_net))
            stnw = np.array((yamlDoc['spatio-temporal-neuron-weights'][1:]))
            stnc = np.array((yamlDoc['spatio-temporal-neuron-edges'][1:]))
            links = np.array([np.vstack((stnw[cidx[0]], stnw[cidx[1]])).tolist()
                                     for cidx in stnc])

            lc = art3d.Line3DCollection([[tuple(con[0,[0,1,zidx]]),
                                          tuple(con[1,[0,1,zidx]])] for con in links], lw=0.5)

            lc.set_color('coral')

            nLines = ax.add_collection(lc)
            nScatter = ax.scatter(stnw[:,0],stnw[:,1],stnw[:,zidx],
                                  s=30, c='crimson', edgecolor='coral', marker='o',
                                  depthshade=False)
        leg_art.append((nLines,nScatter))
        leg_lable.append('Spatio-temporal Neurons')


    legend = ax.legend(leg_art,leg_lable,
                       fontsize='large', loc='lower right')

    legend.get_frame().set_edgecolor('darkgray')

    # ax.set_title('Spatio-temporal Layer Training')

    if len(frames) > 0:
        ani = animation.ArtistAnimation(fig, frames, interval=1500, repeat_delay=3000)

    if gif:
        if len(frames) > 0:
            ani.save('harp_plot.gif', writer='imagemagick')
        else:
            fig.savefig('harp_plot.gif')

    if png:
        if len(frames) > 0:
            ani.save('harp_plot.png', writer='imagemagick')
        else:
            fig.savefig('harp_plot.png')

    if show or not (gif or png):
        print('Showing plot')
        plt.show()


###############################################################################
# Data manipulation
###############################################################################
def data(harp_file_list, output):
    import yaml
    import numpy as np
    import csv
    import os

    #make csv directory, if it exists ignore error
    if not os.path.isdir(output):
        os.makedirs(output)

    for f in open(harp_file_list).read().splitlines():
        print(f)

        yamlDoc = yaml.load(open(f))

        (csvdir, csvfile) = os.path.split(f)

        csvdir = os.path.join(output, csvdir)
        csvfile = os.path.splitext(csvfile)[0] + '.' + output

        if not os.path.isdir(csvdir):
            os.makedirs(csvdir)

        csvpath = os.path.join(csvdir,csvfile)

        with open(csvpath, 'w') as csvDoc:
            wr = csv.writer(csvDoc)

            for row in yamlDoc['events']:
                wr.writerow(row)


###############################################################################
# Command parser
###############################################################################
if __name__ == "__main__":
    import argparse
    import time

    commandParser = argparse.ArgumentParser(description="harp command line utility")

    commandParser.add_argument('command',
                               action='store',
                               choices=['train', 'eval', 'plot','data'],
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

        plotParser.add_argument('-z','--zdata',
                                action='store_true',
                                help='Plot z spacial data on the z-axis')

        plotParser.add_argument('--gif',
                                action='store_true')

        plotParser.add_argument('--png',
                                action='store_true')

        plotParser.add_argument('--show',
                                action='store_true')

        plotParser.add_argument('harp_files',
                                nargs='+',
                                help='One or more harp files to be plotted')

        plotArgs = plotParser.parse_args(commandArgs.args)

        plot(plotArgs.harp_files, plotArgs.zdata,
             plotArgs.gif, plotArgs.png, plotArgs.show)


    elif commandArgs.command == 'data':
        dataParser = argparse.ArgumentParser("manipulate harp data")

        dataParser.add_argument('-o','--output',
                                action='store',
                                default='csv',
                                choices=['csv'],
                                help='Select output format')

        dataParser.add_argument('harp_file_list',
                                # nargs='+',
                                help='One or more harp files')

        dataArgs = dataParser.parse_args(commandArgs.args)

        data(dataArgs.harp_file_list, dataArgs.output)
