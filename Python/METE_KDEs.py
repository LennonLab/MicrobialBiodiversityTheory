from  __future__ import division
import generate_figs_zipf as gf
import os
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import FuncFormatter
import matplotlib
import random
from sklearn.neighbors import KernelDensity

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'

#datasets = [ 'HMP','EMPclosed','MGRAST']
datasets = ['HMP', 'EMPclosed']
#datasets = [ 'EMPclosed']

#methods = ['geom', 'mete','zipf']
methods = ['zipf',]
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
            if j == 0:
                ax.plot(IN1.Grid, IN1.PDF, linewidth=2, alpha=0.5, color =colors[j], label='Broken-stick')
            elif j == 1:
                ax.plot(IN1.Grid, IN1.PDF, linewidth=2, alpha=0.5, color =colors[j], label='METE')
            elif j == 2:
                ax.plot(IN1.Grid, IN1.PDF, linewidth=2, alpha=0.5, color =colors[j], label='Zipf')

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
        ax.legend(loc='upper left', shadow=False)
        plt.xlabel(r'$r^{2}_{m}$', fontsize = 14)
        plt.ylabel(r'$probability$', fontsize = 14)
        if i == 0:
            plt.title('HMP', fontsize = 18)
        if i == 1:
            plt.title('EMP', fontsize = 18)
        if i == 2:
            plt.title('MG-RAST', fontsize = 18)
        plt.xlim(-7, 1)
        plt.savefig(figure_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        plt.close()

def random_lines(datasets, methods, sample_size, iterations):
    df = pd.DataFrame( range(iterations))
    for i, dataset in enumerate(datasets):
        if (dataset != 'HMP') or (dataset != 'EMPclose'):
            pass
        for j, method in enumerate(methods):
            #mean_array = np.empty(iterations, dtype=float)
            mean_list = []
            name = str(dataset) + '_' + str(method)
            print name
            for k in range(iterations):
                IN = gf.import_NSR2_data(mydir + 'NSR2/' + method + '_' + dataset + '_NSR2.txt')
                r2s = ((IN["R2"]))
                sample_mean = np.mean( np.random.choice(r2s, sample_size) )
                mean_list.append( sample_mean )
                #mean_array = np.insert(mean_array, k, sample_mean)
                #np.append(mean_array,sample_mean )
            #data = pd.DataFrame({name: mean_list})
            df[name] = pd.Series(mean_list, index=df.index)
            #df = df.append(data, ignore_index = True)
    OUT_name = 'Mean_R2_Random_Sample.txt'
    OUT_dir = mydir
    df.to_csv(os.path.join(OUT_dir, OUT_name),  sep='\t')
    print df.shape



def r2_random_sample_kde(datasets):
    path = mydir + 'Mean_R2_Random_Sample.txt'
    IN1 = pd.read_csv(path, sep='\t')
    for x in datasets:
        figure_name = str(x) + '.png'
        fig = plt.figure()
        ax = fig.add_subplot(111)
        if x == 'EMPclosed':

            EMPclosed_zipf = IN1['EMPclosed_zipf'].values
            EMPclosed_mete = IN1['EMPclosed_mete'].values
            EMPclosed_geom = IN1['EMPclosed_geom'].values

            #EMPclosed_zipf_kde = gf.CV_KDE(EMPclosed_zipf)
            #EMPclosed_mete_kde = gf.CV_KDE(EMPclosed_mete)
            #EMPclosed_geom_kde = gf.CV_KDE(EMPclosed_geom)

            #ax.plot(EMPclosed_zipf_kde[0], EMPclosed_zipf_kde[1], linewidth=2, alpha=0.5, color = 'blue', label='Zipf')
            #ax.plot(EMPclosed_mete_kde[0], EMPclosed_mete_kde[1], linewidth=2, alpha=0.5, color = 'green', label='METE')
            #ax.plot(EMPclosed_geom_kde[0], EMPclosed_geom_kde[1], linewidth=2, alpha=0.5, color = 'red', label='Broken-stick')
            
            #kde = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(X.ravel()[:, None])
            #log_dens = kde.score_samples(X_plot)  # evaluate the density model on the data.
            #plt.plot(np.exp(log_dens))

            #ax.hist(EMPclosed_zipf, 15, linewidth=2, normed=True, color='blue',alpha=0.5)
            #ax.hist(EMPclosed_mete, 15, linewidth=2, normed=True, color='green',alpha=0.5)
            #ax.hist(EMPclosed_geom, 15, linewidth=2, normed=True, color='red',alpha=0.5)

            plt.savefig(figure_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

        #elif x == 'HMP':
        #    HMP_zipf = IN1['HMP_zipf'].values
        #    HMP_mete = IN1['HMP_mete'].values
        #    HMP_geom = IN1['HMP_geom'].values
        #    ax.hist(HMP_zipf, 15, linewidth=2, color='blue',alpha=0.5)
        #    ax.hist(HMP_mete, 15, linewidth=2, color='green',alpha=0.5)
        #    ax.hist(HMP_geom, 15, linewidth=2, color='red',alpha=0.5)

        #    plt.savefig(figure_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)





#generate_kde_to_file(datasets, methods)
#N_r2_KDE(datasets, methods)
sample_size = 1000
iterations = 10000
#random_lines(datasets, methods, sample_size, iterations)

r2_random_sample_kde(datasets)
