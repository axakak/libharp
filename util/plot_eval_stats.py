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


# evalStats = mlab.csv2rec(args.file, skiprows=1)
evalStats = np.recfromcsv(args.file,skip_header=1,
                            names=('trace','class','err_non','err_dom'))

evalStats['err_non'] = evalStats['err_non']*100
evalStats['err_dom'] = evalStats['err_dom']*100

evalStats = append_fields(evalStats,'err_diff', (evalStats['err_non']-evalStats['err_dom']))

#create 2 col array of dominant and nondominant error differences
errDiff = hstack((vstack(evalStats[evalStats['class']%2 == 1]['err_diff']),
           vstack(evalStats[evalStats['class']%2 == 0]['err_diff'])))

# x,y = plt.figaspect(.65)*1.5
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(16,5),tight_layout=True)

# plt.rc('text', usetex=True)
# plt.rc('font', family='Computer Modern Roman')

###############################################################################
# Box Plot
###############################################################################
axes[0].set_title('Classification Certainty Distribution')

#set x-axis properties
axes[0].set_xlabel('Trace Input (known class)')

#set y-axis properties
axes[0].set_ylabel('Nerual Network Class Error Difference (%)')

axes[0].axhline(y=0, color='k')

labels = [('Dominant'),('Nondominant')]
axes[0].boxplot(errDiff,labels=labels)

ymax = max([abs(y) for y in axes[0].get_ylim()])

yoffset = (axes[0].get_ybound()[1] - axes[0].get_ybound()[0])/24

axes[0].set_ylim(-ymax, ymax)
axes[0].text(1.5, ymax-yoffset, 'More Likely\nDominant', color='gray', ha='center', va='top')
axes[0].text(1.5, -ymax+yoffset, 'More Likely\nNondominant', color='gray', ha='center')


###############################################################################
# Scatter Plot
###############################################################################
axes[1].set_title('Class Errors')
axes[1].set_xlabel('Dominant Class Error (%)')

axes[1].set_ylabel('Nondominant Class Error (%)')

domXY = (evalStats[evalStats['class']%2 == 1]['err_dom'],
         evalStats[evalStats['class']%2 == 1]['err_non'])

nonXY = (evalStats[evalStats['class']%2 == 0]['err_dom'],
         evalStats[evalStats['class']%2 == 0]['err_non'])

axes[1].scatter(domXY[0],domXY[1],label='dominant',c='r',alpha=0.75,
                linewidth=0.1,edgecolors='gray')

axes[1].scatter(nonXY[0],nonXY[1],label='nondominant',c='b',alpha=0.75,
                linewidth=0.1,edgecolors='gray')

axmax = max([axes[1].get_ybound()[1],axes[1].get_xbound()[1]])
axmin = min([axes[1].get_ybound()[0],axes[1].get_xbound()[0]])
axes[1].set_ylim(axmin, axmax)
axes[1].set_xlim(axmin, axmax)
axes[1].plot(np.linspace(axmin,axmax),np.linspace(axmin,axmax), color='k')

yoffset = (axes[1].get_ybound()[1] - axes[2].get_ybound()[0])/24

axes[1].text(axmin+yoffset, axmax-yoffset, r'More Likely Dominant', color='gray', va='top')
axes[1].text(axmax-yoffset, axmin+yoffset, r'More Likely Nondominant', color='gray', ha='right')

axes[1].legend(loc='upper right', fontsize='medium',scatterpoints=1,title='known class',framealpha=.5)


###############################################################################
#Bar Chart
###############################################################################
axes[2].set_title('Classification Count')
axes[2].set_ylabel('Neural Network Classification Count')

axes[2].set_xlabel('Trace Input (known class)')
axes[2].set_xlim(0.5, 2.5)
axes[2].set_xticks([1,2])
axes[2].set_xticklabels(['Dominant', 'Nondominant'])

axes[2].axhline(y=0, color='k')

barCount = [len(errDiff[errDiff[:,0]>0,0]),
            len(errDiff[errDiff[:,0]<0,0])*-1,
            len(errDiff[errDiff[:,1]>0,1]),
            len(errDiff[errDiff[:,1]<0,1])*-1]

axes[2].bar((1,1,2,2), barCount, width=0.2, align='center',
             color=('r','b','r','b'), edgecolor='gray')

ymax = max([abs(y) for y in axes[2].get_ybound()])
yoffset = (axes[2].get_ybound()[1] - axes[2].get_ybound()[0])/24

axes[2].set_ylim(-ymax,ymax)
axes[2].set_yticklabels([int(abs(y)) for y in axes[2].get_yticks()])
axes[2].text(1.5, ymax-yoffset, 'Classified\nDominant', color='gray', ha='center', va='top')
axes[2].text(1.5, -ymax+yoffset, 'Classified\nNondominant', color='gray', ha='center')

plt.show()
