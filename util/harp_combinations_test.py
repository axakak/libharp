#!/usr/bin/env python3
import harp
import os
import argparse
import time
import re

from pathlib import Path
from datetime import datetime
from multiprocessing import Pool

def harpTrainEval(i):

    outDirStr = 'harp_{}_class_{:0=3}_of_{:0=3}'.format(len(args.class_list_files),i+1,len(traceFiles))

    # create directory for current test run
    Path(outDirStr).mkdir(exist_ok=True)
    os.chdir(outDirStr)

    # create log file
    logStr = 'out.log'
    log = open(logStr, mode='w')

    log.write(os.getcwd()+'\n')

    # create eval and train list files
    evalListStr = 'eval_list.txt'
    evalListFile = open(evalListStr, mode='w')

    trainListStr = 'train_list.txt'
    trainListFile = open(trainListStr, mode='w')

    # create train and eval lists for processing
    evalListFile.writelines(traceFiles[i])
    trainListFile.writelines(traceFiles[:i])
    trainListFile.writelines(traceFiles[i+1:])

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

    # remove ASCII escape codes in returned stdout
    output = re.sub(r'\x1b\[[0-9]{1,2}m', '', output)
    output = re.sub(r'.*\x1b\[1K\x1b\[40D', '', output)
    log.write(output)

    # evaluate neural net using eval_list
    sTime = time.time()
    output = ''.join(harp.eval('crbfNeuralNet.yaml', evalListStr))
    eTime = time.time()
    dTime = eTime-sTime
    print('Evaluations completed in {:n}m{:.0f}s'.format(dTime//60, dTime%60))

    # remove ASCII escape codes in returned stdout
    output = re.sub(r'\x1b\[[0-9]{1,2}m', '', output)
    log.write(output)

    entries = ''

    # write csv header on initial loop
    if i is 0:
        entry = re.findall(r'Known class: [0-9]+.*\n', output, re.DOTALL)[-1]
        entry = re.sub(r'Known class: [0-9]+', r'', entry)
        entry = re.sub(r'\n([0-9]+): [01]\.[0-9]+', r',acc. error (class \1)', entry)
        entry = 'neural net,trace,true class'+ entry
        entries += entry

    # for each trace eval form csv entry and write to file
    for m in re.findall(r'==> Evaluating trace.*\n', output, re.DOTALL):
        entry = re.sub(r'==> Evaluating trace.*/(.*\.yaml).*\n', r'\1', m)
        entry = re.sub(r'Known class: ([0-9]+)', r',\1', entry)
        entry = re.sub(r'\n[0-9]+: ([01]\.[0-9]+)', r',\1', entry)
        entries +=  outDirStr + ',' + entry

    os.chdir('..')

    return entries


if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Run random permutations crbf test ")

    # one list file for each class
    parser.add_argument('class_list_files',
                        nargs='+',
                        type=open,
                        help='one or more files with line separated trace file paths')

    parser.add_argument('-p','--processes',
                        choices=range(1, os.cpu_count()+1),
                        default=os.cpu_count()/2,
                        type=int,
                        help='number of worker processes to use. default: (cpu count)/2')

    args = parser.parse_args()

    traceFiles = [f.read().splitlines(keepends=True) for f in args.class_list_files]
    traceFiles = [[os.path.join('..', '..', f) for f in c] for c in traceFiles]

    # flatten list of file lists
    traceFiles = [tFile for cList in traceFiles for tFile in cList]

    eval_csv = open('eval_out_all.csv', mode='w', buffering=1)

    # create pool of processes and write returned csv entries
    with Pool(processes=args.processes) as p:
        imap_it = p.imap(harpTrainEval, range(len(traceFiles)))

        for x in imap_it:
            eval_csv.write(x)
