from  __future__ import division
import generate_figs_zipf as gf
import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'

#datasets = [ 'HMP','EMPclosed','MGRAST']
datasets = ['HMP']
#methods = ['geom', 'mete','zipf']
methods = ['geom',]
#params = ['N','S', 'N/S']
colors = ['red', 'green', 'blue']


#import_NSR2_data(data_dir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')



def generate_kde_to_file(datasets, methods):
    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            IN = gf.import_NSR2_data(mydir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')
            r2s = ((IN["R2"]))
            r2_kde = gf.CV_KDE(r2s)
            r2_kde_table = pd.DataFrame({'Grid':r2_kde[0], 'PDF':r2_kde[1]})
            OUT_name = method + '_' + dataset + '_' + str(round(r2_kde[2], 3)) + '_KDEs.txt'
            OUT_dir = mydir + 'KDEs/'
            r2_kde_table.to_csv(os.path.join(OUT_dir, OUT_name),  sep='\t')


def N_r2_KDE(datasets, methods):
    bins = np.linspace(-1, 1, 100)
    #fig = plt.figure()
    #fig, ax = plt.subplots()
    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            path = mydir + 'KDEs/' + method+'_'+dataset+'_1.259_KDEs.txt'
            IN1 = pd.read_csv(path, sep='\t')
            #IN1 = IN1[(IN1.Grid >= -1) & (IN1.Grid <= 1)]
            IN2 = gf.import_NSR2_data(mydir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')
            r2s = ((IN2["R2"]))
            weights = np.ones_like(r2s)/float(len(r2s))
            #plt.hist(myarray, weights=weights)
            print r2s
            #plt.hist(r2s, bins, alpha=0.5, label = method, color = colors[j])
            plt.plot(IN1.Grid, IN1.PDF, linewidth=2, alpha=0.5, label='bw=%.2f' % 1.259, color =colors[j])
            plt.hist(r2s, 200, weights=weights, color = colors[j], histtype='stepfilled',  alpha=0.3)
            #ax.hist(r2s, 50,color = colors[j],  alpha=0.3)
            #plt.xlim(-1, 1)


            plt.savefig('test.png')
            #r2_kde = gf.CV_KDE(r2s)

            #N_total = ((IN["N"]))


generate_kde_to_file(datasets, methods)
N_r2_KDE(datasets, methods)
