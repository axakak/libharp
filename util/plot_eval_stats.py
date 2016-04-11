#!/usr/bin/env python3
import argparse

import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from numpy import *
from numpy.lib.recfunctions import append_fields

###############################################################################
# Argument Parser Setup
###############################################################################
parser = argparse.ArgumentParser(description="Plot harp evaluation data")

parser.add_argument('file')

args = parser.parse_args()

header = open(args.file).readline().split(sep=',')
h_names = []
for name in header:
    if name == 'neural net':
        h_names.append('neural_net')
    elif name == 'trace':
        h_names.append('trace')
    elif name == 'true class':
        h_names.append('true_class')
    elif '1' in name:
        h_names.append('err_dom')
    elif '0' in name:
        h_names.append('err_non')

evalStats = np.recfromcsv(args.file,skip_header=1,names=h_names)
evalStats = append_fields(evalStats,'err_diff', (evalStats['err_non']-evalStats['err_dom']))

#create arrayS of dominant and nondominant error differences
errDiff_dom = vstack(evalStats[evalStats['true_class']%2 == 1]['err_diff'])
errDiff_non = vstack(evalStats[evalStats['true_class']%2 == 0]['err_diff'])

# print(hstack([errDiff_dom,errDiff_dom>0]))

trace = evalStats['trace'][0].decode()[:6]

plt.rc('font', family = 'serif', serif = 'CMU Serif')

###############################################################################
# Confusion matrix calculations
###############################################################################
N = len(evalStats)
TP = len(errDiff_dom[errDiff_dom>0]) # true positive (A)
FN = len(errDiff_dom[errDiff_dom<0]) # false negative (C)
FP = len(errDiff_non[errDiff_non>0]) # false positive (B)
TN = len(errDiff_non[errDiff_non<0]) # true negative (D)

PP = TP+FP # predicted positives
PN = FN+TN # predicted negatives
RP = TP+FN # real positives
RN = FP+TN # real negatives

bias = PP/N # label bias
prev = RP/N # population prevalence

prevG2 = prev*(1-prev)
biasG2 = bias*(1-bias)

TPR = TP/(TP+FN) # true positive rate, sensitivity, recall
TNR = TN/(FP+TN) # true negative rate, specificity, inverse recall

PPV = TP/(TP+FP) # positive predictive value, precision
NPV = TN/(FN+TN) # negative predictive value, inverse precision

informedness = TPR+TNR-1 # DeltaP'
markedness =  PPV+NPV-1 # DeltaP

correlation = (TP*TN - FP*FN) / sqrt((TP+FN)*(FP+TN)*(TP+FP)*(FN+TN))

accuracy = (TP+TN) / N

print('\t+R\t-R')
print('+P\t{}\t{}\t{}'.format(TP,FP,PP))
print('-P\t{}\t{}\t{}'.format(FN,TN,PN))
print('')
print('\t{}\t{}\t{}'.format(RP,RN,N))
print('')

print('Informedness: {:.4f}'.format(informedness))
print('Markedness: {:.4f}'.format(markedness))
print('Correlation: {:.4f}'.format(correlation))
print('')

print('Bias: {:.4f}'.format(bias))
print('Prevalence: {:.4f}'.format(prev))

print('BiasG2: {:.4f}'.format(biasG2))
print('PrevG2: {:.4f}'.format(prevG2))
print('')

print('Recall: {:.4f}'.format(TPR))
print('Precision: {:.4f}'.format(PPV))
print('IRecall: {:.4f}'.format(TNR))
print('IPrecision: {:.4f}'.format(NPV))
print('')

print('Accuracy: {:.4f}'.format(accuracy))
print('')

obs = np.array([[TP, FP], [FN, TN]])

oddsratio, pvalue = stats.fisher_exact(obs)
print('Fisher exact pvalue: {:.8f}'.format(pvalue))
print('Fisher exact odds: {:.8f}'.format(oddsratio))
print('')

###############################################################################
# Box Plot
###############################################################################
x,y = plt.figaspect(.9)
plt.figure(figsize=(x,y))
ax = plt.gca()

ax.set_title('Class Error Difference Distribution')

#set x-axis properties
ax.set_xlabel('Trace Input (True Class)')

#set y-axis properties
ax.set_ylabel('Class Error Difference $(c_{non}-c_{dom})$')

ax.axhline(y=0, color='gray')

labels = [('Dominant'),('Nondominant')]
ax.boxplot([errDiff_dom,errDiff_non],labels=labels)

ymax = max([abs(y) for y in ax.get_ylim()])

yoffset = (ax.get_ybound()[1] - ax.get_ybound()[0])/24

ax.set_ylim(-ymax, ymax)
ax.text(1.5, ymax-yoffset, 'More Likely\nDominant', color='gray', ha='center', va='top')
ax.text(1.5, -ymax+yoffset, 'More Likely\nNondominant', color='gray', ha='center')


###############################################################################
# Scatter Plot
###############################################################################
plt.figure(figsize=(x,y))
ax = plt.gca()

ax.set_title('Class Errors')
ax.set_xlabel('Dominant Class Error')

ax.set_ylabel('Nondominant Class Error')

domXY = (evalStats[evalStats['true_class']%2 == 1]['err_dom'],
         evalStats[evalStats['true_class']%2 == 1]['err_non'])

nonXY = (evalStats[evalStats['true_class']%2 == 0]['err_dom'],
         evalStats[evalStats['true_class']%2 == 0]['err_non'])

ax.scatter(domXY[0],domXY[1],label='dominant',c='r',alpha=0.75,
                linewidth=0.1,edgecolors='gray')

ax.scatter(nonXY[0],nonXY[1],label='nondominant',c='b',alpha=0.75,
                linewidth=0.1,edgecolors='gray')

axmax = max([ax.get_ybound()[1],ax.get_xbound()[1]])
axmin = min([ax.get_ybound()[0],ax.get_xbound()[0]])

ax.set_ylim(axmin, axmax)
ax.set_xlim(axmin, axmax)
ax.plot(np.linspace(axmin,axmax),np.linspace(axmin,axmax), color='gray')

yoffset = (ax.get_ybound()[1] - ax.get_ybound()[0])/24

ax.text(axmin+yoffset, axmax-yoffset, r'More Likely Dominant', color='gray', va='top')
ax.text(axmax-yoffset, axmin+yoffset, r'More Likely Nondominant', color='gray', ha='right')

ax.legend(loc='upper right', fontsize='medium',scatterpoints=1,title='True Class',framealpha=.75)


###############################################################################
#Bar Chart
###############################################################################
plt.figure(figsize=(x,y))
ax = plt.gca()

ax.set_title('Predicted Class Count')
ax.set_ylabel('Neural Network Predicted Class Count')

ax.set_xlabel('Trace Input (True Class)')
ax.set_xlim(0.5, 2.5)
ax.set_xticks([1,2])
ax.set_xticklabels(['Dominant', 'Nondominant'])

ax.axhline(y=0, color='gray')


barCount = [TP,-FN,FP,-TN]

ax.bar((1,1,2,2), barCount, width=0.2, align='center',
             color=('r','b','r','b'), edgecolor='gray')

ymax = max([abs(y) for y in ax.get_ybound()])
yoffset = (ax.get_ybound()[1] - ax.get_ybound()[0])/24

ax.set_ylim(-ymax,ymax)
ax.set_yticklabels([int(abs(y)) for y in ax.get_yticks()])
ax.text(1.5, ymax-yoffset, 'Classified\nDominant', color='gray', ha='center', va='top')
ax.text(1.5, -ymax+yoffset, 'Classified\nNondominant', color='gray', ha='center')


for fig in plt.get_fignums():
	plt.figure(num=fig)
	plt.savefig('{}{}.pdf'.format(trace,plt.gca().get_title().replace(' ','_').lower()))

plt.show()
