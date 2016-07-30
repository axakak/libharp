#!/usr/bin/env python3
import harp
import os
import argparse
import time
import re

from pathlib import Path
from multiprocessing import Pool

def getOutDirStr(i):

    return 'harp_{}_class_{:0=3}_of_{:0=3}'.format(len(args.class_list_files),i+1,traceCount)

def harpTrainEval(i):

    outDirStr = getOutDirStr(i)

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

    _harpTrain(trainListStr, log)
    _harpEval(i, evalListStr, outDirStr, log)


def _harpTrain(trainListStr, log):
    #HACK: harptrain has a sporadic seg fault, try 10 times then stop
    for k in range(10):
        # train neural net using train_list
        sTime = time.time()
        output = ''.join(harp.harp_train(trainListStr))
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



def harpEval(outDirStr):
    # change directory for current index
    os.chdir(outDirStr)

    # create log file
    logStr = 'out.log'
    log = open(logStr, mode='a')

    entries = _harpEval('eval_list.txt', outDirStr, log)

    os.chdir('..')

    return entries

def _harpEval(evalListStr, outDirStr, log):
    # evaluate neural net using eval_list
    sTime = time.time()
    output = ''.join(harp.harp_eval('crbfNeuralNet.yaml', evalListStr))
    eTime = time.time()
    dTime = eTime-sTime
    print('Evaluations completed in {:n}m{:.0f}s'.format(dTime//60, dTime%60))

    # remove ASCII escape codes in returned stdout
    output = re.sub(r'\x1b\[[0-9]{1,2}m', '', output)

    log.write(output)

    entries = ''

    # write csv header on initial loop
    if '01_' in outDirStr:
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

    return entries



if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Run random permutations crbf test ")

    # one list file for each class

    parser.add_argument('-p','--processes',
                        choices=range(1, os.cpu_count()+1),
                        default=os.cpu_count()/2,
                        type=int,
                        help='number of worker processes to use. default: (cpu count)/2')

    parser.add_argument('-e','--eval',
                        action='store_true',
                        help='evaluate pre-trained neural nets (note: no lass_list_files needed)')

    parser.add_argument('class_list_files',
                        nargs='*',
                        type=open,
                        help='one or more files with line separated trace file paths')


    args = parser.parse_args()

    if args.eval:
        subdirs = [str(x) for x in Path('.').iterdir() if x.is_dir()]
    else:
        traceFiles = [f.read().splitlines(keepends=True) for f in args.class_list_files]
        traceFiles = [[os.path.join('..', '..', f) for f in c] for c in traceFiles]

        # flatten list of file lists into single list
        traceFiles = [tFile for cList in traceFiles for tFile in cList]

        traceCount = len(traceFiles)

    eval_csv = open('eval_out_all.csv', mode='w', buffering=1)

    # create pool of processes and write returned csv entries
    with Pool(processes=args.processes) as p:
        if args.eval:
            imap_it = p.imap(harpEval, subdirs)
        else:
            imap_it = p.imap(harpTrainEval, range(traceCount))

        for x in imap_it:
            eval_csv.write(x)
