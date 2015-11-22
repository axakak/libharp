#!/usr/bin/env python3
import harp
import os
import argparse
import sys
import time
import re

from pathlib import Path
from datetime import datetime
from random import shuffle

# Argument parsing
parser = argparse.ArgumentParser(description="Run random permutations crbf test ")

# one list file for each class
parser.add_argument('class_list_files',
                    nargs='+',
                    type=open,
                    help='two or more files with line separated trace file paths')

args = parser.parse_args()

# check for a minimum of two classes
if len(args.class_list_files) < 2:
    print('error: class_list_files requiers two or more list files')
    quit()

traceFiles = [f.read().splitlines(keepends=True) for f in args.class_list_files]
traceFiles = [[os.path.join('..', '..', f) for f in c] for c in traceFiles]

#TODO: Remove after testing
traceFiles = [traceFiles[0][:2], traceFiles[1][:2]]

# flatten list of file lists
traceFiles = [tFile for cList in traceFiles for tFile in cList]

print(traceFiles)

eval_csv = open('eval_out_all.csv', mode='w', buffering=1)

for i in range(len(traceFiles)):

    print()

    outDirStr = 'harp_{}_class_{}_of_{}'.format(len(args.class_list_files),i+1,len(traceFiles))

    # create directory for current test run
    outDir = Path(outDirStr)
    outDir.mkdir(exist_ok=True)
    os.chdir(outDirStr)

    # create log file
    log = open('out.log', mode='w')

    print(os.getcwd())
    log.write(os.getcwd()+'\n')

    # create eval and train list files
    evalListStr = 'eval_list.txt'
    evalListFile = open(evalListStr, mode='w')

    trainListStr = 'train_list.txt'
    trainListFile = open(trainListStr, mode='w')

    # create train and eval lists for processing
    evalListFile.writelines(traceFiles[i])
    trainListFile.writelines(traceFiles[:i])
    trainListFile.writelines(traceFiles[i+1:])#HACK: might cause bounds error

    # close list files
    trainListFile.close()
    evalListFile.close()

    #HACK: harptrain has a sporadic seg fault, try 10 times then stop
    for k in range(10):
        # train neural net using train_list
        sTime = time.time()
        output = ''.join(harp.train(trainListStr))
        eTime = time.time()
        dTime = eTime-sTime
        if 'Exporting neural network' in output:
            print('Training completed in {:n}m{:.0f}s'.format(dTime//60, dTime%60))
            break
        else:
            print('\nError: harptrain crashed, will retry')
    else:
        print('\n\n****** Error: retry max reached **********\n\n')


    # clean up and log terminal output
    output = re.sub(r'\x1b\[[0-9]{1,2}m', '', output)
    output = re.sub(r'.*\x1b\[1K\x1b\[40D', '', output)
    log.write(output)

    # for line in output:
    #     if 'Detecting classes:' in line:
    #         classList = re.findall(r'[0-9]+', line)

    # evaluate neural net using eval_list
    sTime = time.time()
    output = ''.join(harp.eval('crbfNeuralNet.yaml', evalListStr))
    eTime = time.time()
    dTime = eTime-sTime
    print('Evaluations completed in {:n}m{:.0f}s'.format(dTime//60, dTime%60))

    # clean up terminal output
    output = re.sub(r'\x1b\[[0-9]{1,2}m', '', output)
    log.write(output)

    # write csv header on initial loop
    if i is 0:
        entry = re.findall(r'Known class: [0-9]+.*\n', output, re.DOTALL)[-1]
        entry = re.sub(r'Known class: [0-9]+', r'', entry)
        entry = re.sub(r'\n([0-9]+): [01]\.[0-9]+', r',acc. error (class \1)', entry)
        entry = 'neural net,trace,true class'+ entry
        eval_csv.write(entry)

    # for each trace eval form csv entry and write to file
    for m in re.findall(r'==> Evaluating trace.*\n', output, re.DOTALL):
        entry = re.sub(r'==> Evaluating trace.*/(.*\.yaml).*\n', r'\1', m)
        entry = re.sub(r'Known class: ([0-9]+)', r',\1', entry)
        entry = re.sub(r'\n[0-9]+: ([01]\.[0-9]+)', r',\1', entry)
        eval_csv.write(outDirStr+','+entry)

    os.chdir('..')
