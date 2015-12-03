#!/usr/bin/env python3
import argparse

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


evalStats = np.recfromcsv(args.file,skip_header=1,
                        names=('neural_net','trace','true_class','err_non','err_dom'))

evalStats = append_fields(evalStats,'err_diff', (evalStats['err_non']-evalStats['err_dom']))

#create 2 col array of dominant and nondominant error differences
errDiff = hstack((vstack(evalStats[evalStats['true_class']%2 == 1]['err_diff']),
           vstack(evalStats[evalStats['true_class']%2 == 0]['err_diff'])))

trace = evalStats['trace'][0].decode()[:6]

plt.rc('font', family = 'serif', serif = 'CMU Serif')

###############################################################################
# Confusion matrix calculations
###############################################################################
TP = len(errDiff[errDiff[:,0]>0,0]) #true positive
FN = len(errDiff[errDiff[:,0]<0,0]) #false negative
FP = len(errDiff[errDiff[:,1]>0,1]) #false positive
TN = len(errDiff[errDiff[:,1]<0,1]) #true negative

TPR = TP/(TP+FN)
TNR = TN/(FP+TN)

PPV = TP/(TP+FP)
NPV = TN/(FN+TN)


print('True positive: {}'.format(TP))
print('False negative: {}'.format(FN))

print('False positive: {}'.format(FP))
print('True negative: {}'.format(TN))

print('True positive rate (sensitivity): {}'.format(TPR))
print('True negative rate (specificity): {}'.format(TNR))

print('Positive predictive value: {:.2f}'.format(PPV))
print('Negative predictive value: {:.2f}'.format(NPV))


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
ax.boxplot(errDiff,labels=labels)

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
