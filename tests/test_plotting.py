from  __future__ import division
import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
import matplotlib
import random
from sklearn.neighbors import KernelDensity

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + '/data/'

path = mydir + 'Mean_R2_Random_Sample.txt'
IN1 = pd.read_csv(path, sep='\t')
x = IN1['HMP_zipf'].values
#test = np.histogram(x, bins=50,range = ([0,1]), density=True)
weights = np.ones_like(x) / x.shape[-1]
fig = plt.figure()
ax = fig.add_subplot(111)
n, bins, patches = ax.hist(x, 50, normed=1, weights = weights, histtype='step', cumulative=False)
kde = stats.gaussian_kde(x)
xx = np.linspace(0, 1, len(x))
ax.plot(xx, kde(x))

plt.savefig('temp.png')

#x = np.random.randn(1000)
#print x
#n = 10
# Count the occurrences in the sample.
#b = np.bincount(x, minlength=n+1)

# p is the array of probabilities.
#p = b / float(b.sum())
#plt.bar(np.arange(len(b)) - 0.5, p, width=1, facecolor='white')
#plt.xlim(-0.5, n + 0.5)
#plt.xlabel("Number of heads (k)")
#plt.ylabel("P(k)")


#fig = plt.figure()
#ax = fig.add_subplot(111)
#n, bins, rectangles = ax.hist(x, 50, normed=False)

#fig.canvas.draw()
#plt.savefig('temp.png')
