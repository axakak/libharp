#!/usr/bin/env python3
import harp
import os
import argparse
import sys
import time

from pathlib import Path
from datetime import datetime
from random import shuffle

# Argument parsing
parser = argparse.ArgumentParser(description="Run random permutations crbf test ")

# number of traces per group to use in evaluation
parser.add_argument('-n',
                    action='store',
                    metavar='<count>',
                    type=int,
                    default=5,
                    help='use <count> traces from each group to evaluate training')

# one list file for each group
parser.add_argument('group_list_files',
                    nargs='+',
                    type=open,
                    help='two or more files listing groups of traces')

args = parser.parse_args()


# check for a minimum of two groups
if len(args.group_list_files) < 2:
    print('error: group_list_files requiers two or more list files')
    quit()

#TODO: create better output directory name
timeStr = datetime.today().strftime("%Y-%m-%d_%H%M")
outDirStr = 'harp_random_{}_groups_{}'.format(len(args.group_list_files),timeStr)

# create directory for current test results
outDir = Path(outDirStr)
outDir.mkdir(exist_ok=True)
os.chdir(outDirStr)

evalListStr = 'eval_list.txt'
trainListStr = 'train_list.txt'

evalListFile = open(evalListStr, mode='w')
trainListFile = open(trainListStr, mode='w')

# create train and eval lists for processing
for group in args.group_list_files:
    # create list of trace files
    traceFiles = group.read().splitlines(keepends=True)
    traceFiles = [os.path.join('..',f) for f in traceFiles]

    # shuffle list of traces
    shuffle(traceFiles)

    # split group into evaluate and train subgroups
    evalListFile.writelines(traceFiles[0:args.n])
    trainListFile.writelines(traceFiles[args.n:])

# close list files
trainListFile.close()
evalListFile.close()

# train neural net using train_list
sTime = time.time()
harp.train(trainListStr)
eTime = time.time()
dTime = eTime-sTime
print('Training completed in {:n}m{:.0f}s'.format(dTime//60, dTime%60))


# evaluate neural net using eval_list
sTime = time.time()
harp.eval('crbfNeuralNet.yaml', evalListStr)
eTime = time.time()
dTime = eTime-sTime
print('Evaluations completed in {:n}m{:.0f}s'.format(dTime//60, dTime%60))
