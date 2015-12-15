from  __future__ import division
import generate_figs_zipf as gf
import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
import matplotlib

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'

datasets = [ 'HMP','EMPclosed','MGRAST']
#datasets = ['HMP']
methods = ['geom', 'mete','zipf']
#methods = ['geom',]
#params = ['N','S', 'N/S']
colors = ['red', 'green', 'blue']


#import_NSR2_data(data_dir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')



def generate_kde_to_file(datasets, methods):
    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            IN = gf.import_NSR2_data(mydir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')
            r2s = ((IN["R2"]))
            #r2s = r2s[(r2s >= -1) & (r2s <= 1)]
            r2_kde = gf.CV_KDE(r2s)
            r2_kde_table = pd.DataFrame({'Grid':r2_kde[0], 'PDF':r2_kde[1]})
            OUT_name = method + '_' + dataset + '_KDEs.txt'
            OUT_dir = mydir + 'KDEs/'
            r2_kde_table.to_csv(os.path.join(OUT_dir, OUT_name),  sep='\t')


def to_percent(y, position):
    # Ignore the passed in position. This has the effect of scaling the default
    # tick locations.
    s = str(100 * y)

    # The percent symbol needs escaping in latex
    if matplotlib.rcParams['text.usetex'] is True:
        return s + r'$\%$'
    else:
        return s + '%'


def N_r2_KDE(datasets, methods):
    bins = np.linspace(-1, 1, 100)
    #fig = plt.figure()
    for i, dataset in enumerate(datasets):
        fig, ax = plt.subplots()
        max_x = 0
        min_x = 0
        for j, method in enumerate(methods):
            path = mydir + 'KDEs/' + method+'_'+dataset+'_KDEs.txt'
            IN1 = pd.read_csv(path, sep='\t')
            #IN1 = IN1[(IN1.Grid >= -1) & (IN1.Grid <= 1)]
            IN2 = gf.import_NSR2_data(mydir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')
            r2s = ((IN2["R2"]))
            #r2s = r2s[(r2s >= -1) & (r2s <= 1)]
            #plt.hist(r2s, bins, alpha=0.5, label = method, color = colors[j])
            plt.plot(IN1.Grid, IN1.PDF, linewidth=2, alpha=0.5, color =colors[j])
            #plt.hist(r2s, 30, stacked = True, color = colors[j],   alpha=0.3, normed= True)
            #ax.hist(r2s, 50,color = colors[j],  alpha=0.3)
            #n, bins, patches = plt.hist(r2s, 50, normed=1, facecolor='green', alpha=0.75)
            if np.max(IN1.PDF) > max_x:
                max_x = np.max(IN1.PDF)
            if np.min(IN1.PDF) < min_x:
                min_x = np.min(IN1.PDF)
        print np.min(IN1.PDF)
        figure_name = '../figures/' + str(dataset) + '_KDE.png'
        plt.grid(True)
        plt.xlabel(r'$r^{2}_{m}$', fontsize = 14)
        plt.ylabel(r'$probability$', fontsize = 14)
        plt.xlim(-7, 1)
        plt.savefig(figure_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        plt.close()


#generate_kde_to_file(datasets, methods)
N_r2_KDE(datasets, methods)
