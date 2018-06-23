#!/usr/bin/python
import matplotlib
matplotlib.use('agg')
from pylab import *
import sys
import os
pwd = os.path.abspath('.')
files = sys.argv[1:]
plotme = True
if '--noplot' in files:
    files.remove('--noplot') 
    plotme = False

for fn in files:
    xs, ys=[],[]
    best = {}
    with open(fn) as file:
        lines = file.readlines()
    for line in lines:
        xs.append(int(line.split()[1]))
        ys.append(float(line.split()[3]))
        if xs[-1] not in best:
            best[xs[-1]] = []
        best[xs[-1]].append(ys[-1])
    for i in best:
        best[i].sort()
    lowest = [best[x][0] for x in best]
    lowest.sort()
    lowest.reverse()
    if plotme and not os.path.isfile(fn+'.png'):
        Figure()
        plot(xs,ys,'.')
        ylim(0,5)
        savefig(fn+'.png')
        close()
    print fn, lowest[0], lowest[1], lowest[2], lowest[3], lowest[4], lowest[5], 'sum6', sum(lowest[0:6]), 'sum', sum(lowest)
    os.chdir(pwd)
