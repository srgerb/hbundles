import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

name1=str(sys.argv[1])
name2=str(sys.argv[2])
#print (name)

data1 = np.genfromtxt("%s"%(name1), skip_header=1, usecols=("rms","score"), names=True)
data2 = np.genfromtxt("%s"%(name2), skip_header=1, usecols=("rms","score"), names=True)

fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(data1['rms'], data1['score'], linestyle='None', markersize=1, marker='x', markeredgewidth=2, color='r', alpha=1)
ax.plot(data2['rms'], data2['score'], linestyle='None', markersize=1, marker='x', markeredgewidth=2, color='g', alpha=1)

ax.set_xlim(0,15)
#ax.set_ylim(min(data['y']),statistics.median(data['y'])+statistics.median(data['y'])*.5)
ax.set_ylim(min(data2['score'] - 10),np.median(data1['score']))
ax.set_xlabel('RMSD')
ax.set_ylabel('Rosetta Energy')
ax.set_title("%s"%(name1))

plt.savefig(filename="plots/%s.pdf"%(name1), format='pdf', transparent=True)

#plt.show()
