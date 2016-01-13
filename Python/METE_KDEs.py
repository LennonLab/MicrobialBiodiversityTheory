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
import signal
import mete
import macroecotools
import macroeco_distributions as md
from scipy import stats
import random
import collections
from scipy.stats import gaussian_kde
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity


mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'


#import_NSR2_data(data_dir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')

def CV_KDE(x):
    # remove +/- inf
    x = x[np.logical_not(np.isnan(x))]
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.logspace(0.1, 5.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(x[:, None])
    x_grid = np.linspace(np.min(x), np.max(x), len(x))
    kde = grid.best_estimator_
    print "bandwidth is " + str(kde.bandwidth)
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple



def generate_kde_to_file(datasets, methods):
    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            IN = gf.import_NSR2_data(mydir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')
            r2s = ((IN["R2"]))
            #r2s = r2s[(r2s >= -1) & (r2s <= 1)]
            r2_kde = CV_KDE(r2s)
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


def r2_KDE(datasets, methods):
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
        figure_name = '../figures/KDEs/' + str(dataset) + '_KDE.png'
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


def get_GeomSeries(N,S,zeros):

    rank = range(1,S+1)
    cdf = [(S-i+0.5)/S for i in rank]
    SNratio = S/N
    if zeros == False:
        abd = md.trunc_geom.ppf(np.array(cdf), SNratio, N)
    return abd


class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException



def sample_lines(datasets, SAD_number, iterations):
    percents = [0.500000, 0.250000, 0.125000, 0.062500, 0.031250, 0.015625]
    SAD_number = int(SAD_number)
    iterations = int(iterations)
    methods = ['geom', 'mete','zipf']
    # testing percents
    #percent = 0.25
    #df = pd.DataFrame( range(iterations))
    for i, dataset in enumerate(datasets):
        signal.signal(signal.SIGALRM, timeout_handler)
        if dataset == 'MGRAST':
            # fix subset l8r
            IN = mydir  + dataset + '-Data' + '/MGRAST/MGRAST-SADs.txt'
            nsr2_data_zipf = gf.import_NSR2_data(mydir + 'NSR2/' + 'zipf_MGRAST_NSR2.txt')
            nsr2_data_mete_geom = gf.import_NSR2_data(mydir + 'NSR2/' + 'mete_MGRAST_NSR2.txt')
        elif dataset == '95' or dataset == '97' or dataset == '99':
            IN = mydir  + dataset + '-Data/' + str(dataset) + '/MGRAST-' + str(dataset) + '-SADs.txt'
            nsr2_data_zipf = gf.import_NSR2_data(mydir + 'NSR2/' +'zipf_MGRAST'+dataset+'_NSR2.txt')
            nsr2_data_mete_geom = gf.import_NSR2_data(mydir + 'NSR2/' + 'mete_MGRAST'+dataset+'_NSR2.txt')
        elif dataset == 'HMP':
            IN = mydir  + dataset + '-Data' + '/' + dataset +'-SADs_NAP.txt'
            nsr2_data_zipf = gf.import_NSR2_data(mydir + 'NSR2/' +  'zipf_'+dataset+'_NSR2.txt')
            nsr2_data_mete_geom = gf.import_NSR2_data(mydir + 'NSR2/' + 'mete_'+dataset+'_NSR2.txt')
        else:
            IN = mydir  + dataset + '-Data' + '/' + dataset +'-SADs.txt'
            nsr2_data_zipf = gf.import_NSR2_data(mydir + 'NSR2/' +  'zipf_'+dataset+'_NSR2.txt')
            nsr2_data_mete_geom = gf.import_NSR2_data(mydir + 'NSR2/' + 'mete_'+dataset+'_NSR2.txt')

        #argsort()[::-1][:n]
        # Get lines with x largest values of N.
        nsr2_data_zipf_N_site = np.column_stack((nsr2_data_zipf["site"], nsr2_data_zipf["N"]))
        nsr2_data_mete_geom_N_site = np.column_stack((nsr2_data_mete_geom["site"], nsr2_data_mete_geom["N"]))
        # Sort these arrays
        nsr2_data_zipf_top100 = nsr2_data_zipf_N_site[nsr2_data_zipf_N_site[:,1].argsort()[::-1]][:SAD_number,]
        nsr2_data_mete_geom_top100 = nsr2_data_mete_geom_N_site[nsr2_data_mete_geom_N_site[:,1].argsort()[::-1]][:SAD_number,]
        # Get the SAD numbers
        zipf_numbers = nsr2_data_zipf_top100[:,0]
        mete_geom_numbers = nsr2_data_mete_geom_top100[:,0]
        zipf_numbers = zipf_numbers.astype(int)
        mete_geom_numbers = mete_geom_numbers.astype(int)

        OUT1 = open(mydir + 'SubSampled-Data' + '/' + dataset + 'geom_SubSampled_Data.txt', 'w+')
        OUT2 = open(mydir + 'SubSampled-Data' + '/' + dataset + 'mete_SubSampled_Data.txt', 'w+')
        OUT3 = open(mydir + 'SubSampled-Data' + '/' + dataset + 'zipf_SubSampled_Data.txt', 'w+')
        num_lines = sum(1 for line in open(IN))
        #done = False
        for j,line in enumerate(open(IN)):
            if (j not in zipf_numbers) and (j not in mete_geom_numbers):
                continue
            #elif j in zipf_numbers:
            #    print "yay"
            else:
                print j
                if dataset == "HMP":
                    line = line.strip().split(',')
                    line = [x.strip(' ') for x in line]
                    line = [x.strip('[]') for x in line]
                    site_name = line[0]
                    #line.pop(0)
                    line.pop(0)
                else:
                    line = eval(line)
            obs = map(int, line)
            # Calculate relative abundance of each OTU
            # Use that as weights
            N_0 = sum(obs)
            S_0 = len(obs)
            N_max = max(obs)
            if S_0 < 10 or N_0 <= S_0:
                continue
            line_ra = map(lambda x: x/N_0, obs)
            # Calculate relative abundance of each OTU
            # Use that as weights
            gm_lines = SAD_number
            zipf_lines = SAD_number
            #percent_model_dict = {}

            geom_means = [N_0, S_0, N_max]
            mete_means = [N_0, S_0, N_max]
            zipf_means = [N_0, S_0, N_max]
            print dataset, N_0, S_0, ' countdown: ', gm_lines, ' ', zipf_lines
            # separate this. get percents for Zipf and mete/geom
            # then go on with the sampling
            
            for percent in percents:
                print percent
                sample_size = round(percent * N_0)
                if sample_size <= 10:
                    continue
                sample_k = np.random.multinomial(sample_size, line_ra, size = None)
                sample_k = sample_k[sample_k != 0]
                sample_k[::-1].sort()
                N_k = sum(sample_k)
                S_k = sample_k.size
                if S_k < 10 or N_k <= S_k:
                    continue
                N_max_k = max(sample_k)

                mg_iter = iterations
                bad_zipf_count = 20
                zipf_iter = iterations + bad_zipf_count

                N_max_list_mg = []
                N_0_list_mg = []
                S_0_list_mg = []
                r2_list_BS = []
                r2_list_METE = []

                N_max_list_zipf = []
                N_0_list_zipf = []
                S_0_list_zipf = []
                r2_list_zipf = []
                gamma_list = []
                len_N_0_list_zipf = len(N_0_list_zipf)
                if j in mete_geom_numbers:
                    while (mg_iter > 0):
                        logSeries = mete.get_mete_rad(S_k, N_k)
                        pred_mete = logSeries[0]
                        r2_mete = macroecotools.obs_pred_rsquare(np.log10(sample_k), np.log10(pred_mete))
                        pred_BS = get_GeomSeries(N_k, S_k, False) # False mean no zeros allowed
                        r2_BS = macroecotools.obs_pred_rsquare(np.log10(sample_k), np.log10(pred_BS))
                        #r2_list = [r2_zipf, r2_mete, r2_BS]
                        r2_list = [r2_mete, r2_BS]

                        if any( (r2 == -float('inf') ) or (r2 == float('inf') ) or (r2 == float('Nan') ) for r2 in r2_list):
                            continue
                        N_max_list_mg.append(N_max_k)
                        N_0_list_mg.append(N_k)
                        S_0_list_mg.append(S_k)
                        r2_list_BS.append(r2_BS)
                        r2_list_METE.append(r2_mete)
                        mg_iter -= 1

                if j in zipf_numbers:
                    while (len_N_0_list_zipf <= iterations) and (zipf_iter > 0):
                        # This try/except loop ensures that
                        #   you'll catch TimeoutException when it's sent.
                        signal.alarm(2)
                        try:
                            # Whatever your function that might hang
                            #Zipf_solve_line = gf.zipf_solver(sample_k)
                            # use S
                            zipf_class = gf.zipf(sample_k)
                            pred_tuple = zipf_class.from_cdf()
                            Zipf_solve_line = zipf_class.zipf_solver(sample_k)
                            rv = stats.zipf(Zipf_solve_line)

                            pred_zipf = pred_tuple[0]
                            gamma = pred_tuple[1]
                            r2_zipf = macroecotools.obs_pred_rsquare(np.log10(sample_k), np.log10(pred_zipf))
                            if (r2_zipf == -float('inf') ) or (r2_zipf == float('inf') ) or (r2_zipf == float('Nan') ):
                                continue
                            else:
                                r2_list_zipf.append(r2_zipf)
                                gamma_list.append(gamma)
                                N_max_list_zipf.append(N_max_k)
                                N_0_list_zipf.append(N_k)
                                S_0_list_zipf.append(N_k)
                                zipf_iter -= 1
                        except TimeoutException:
                            #continue # continue the while loop if function takes more than x seconds
                            pass # rest of script happens as normal
                        else:
                            # Reset the alarm
                            signal.alarm(0)
                if j in mete_geom_numbers:
                    N_0_mg_mean = np.mean(N_0_list_mg)
                    geom_means.append(N_0_mg_mean)
                    mete_means.append(N_0_mg_mean)

                    S_0_mean = np.mean(S_0_list_mg)
                    geom_means.append(S_0_mean)
                    mete_means.append(S_0_mean)


                    N_max_mg_mean = np.mean(N_max_list_mg)
                    geom_means.append(N_max_mg_mean)
                    mete_means.append(N_max_mg_mean)

                    r2_BS_mg_mean = np.mean(r2_list_BS)
                    geom_means.append(r2_BS_mg_mean)
                    r2_METE_mg_mean = np.mean(r2_list_METE)
                    mete_means.append(r2_METE_mg_mean)

                if j in zipf_numbers:
                    N_0_zipf_mean = np.mean(N_0_list_zipf)
                    zipf_means.append(N_0_zipf_mean)
                    S_0_zipf_mean = np.mean(S_0_list_zipf)
                    zipf_means.append(S_0_zipf_mean)
                    N_max_zipf_mean = np.mean(N_max_list_zipf)
                    zipf_means.append(N_max_zipf_mean)
                    r2_zipf_mean = np.mean(r2_list_zipf)
                    zipf_means.append(r2_zipf_mean)
                    gamma_mean = np.mean(gamma_list)
                    zipf_means.append(gamma_mean)

                    print percent
                    print N_0_zipf_mean, r2_zipf_mean

            '''Now we check if the lists are the right length
            there are 6 iterations for the percentage
            mete/ geom, append four items each iteration.
            4*6 = 24, add three original = 27
            likewise, for zipf, (5*6) + 3 = 33 '''
            #print geom_means
            #print mete_means
            #print zipf_means
            print zipf_means
            print len(zipf_means)
            if len(geom_means) == 27:
                geom_means_str = ' '.join(map(str, geom_means))
                #OUT1.write(','.join(map(repr, geom_means_str[i]))
                print>> OUT1, j, geom_means_str
            if len(mete_means) == 27:
                mete_means_str = ' '.join(map(str, mete_means))
                print>> OUT2, j, mete_means_str
            if len(zipf_means) == 33:
                zipf_means_str = ' '.join(map(str, zipf_means))
                print>> OUT3, j, zipf_means_str
            print dataset, percent





####################################
##### random_lines is Deprecated####
####################################
def random_lines(datasets, methods, sample_size, iterations):
    '''This piece of code is deprecated.'''
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


#########################
#### r2_kde is Deprecated
#########################

def r2_kde(datasets):
    path = mydir + 'Mean_R2_Random_Sample.txt'
    IN1 = pd.read_csv(path, sep='\t')
    for x in datasets:
        figure_name = str(x) + '.png'
        fig = plt.figure()
        ax = fig.add_subplot(111)
        x_zipf = str(x) + '_zipf'
        x_mete = str(x) + '_mete'
        x_geom = str(x) + '_geom'
        zipf = IN1[x_zipf].values
        mete = IN1[x_mete].values
        geom = IN1[x_geom].values
        #n, bins, rectangles = ax.hist(zipf, 50, linewidth=2, normed=1, color='blue',alpha=0.5)
        #fig.canvas.draw()
        #weights = np.ones_like(zipf)/len(zipf)
        #print weights
        #ax.hist(zipf, 50, linewidth=2, normed=1, color='blue',alpha=0.5, density=True)
        #ax.hist(mete, 15, linewidth=2, normed=True, color='green',alpha=0.5)
        #ax.hist(geom, 15, linewidth=2, normed=True, color='red',alpha=0.5)
        #plt.savefig(figure_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        #zipf_kde = CV_KDE(zipf)
        zipf_mete = CV_KDE(mete)
        #zipf_geom = CV_KDE(geom)

        #ax.plot(zipf_kde[0], zipf_kde[1], linewidth=2, alpha=0.5, color = 'blue', label='Zipf')
        ax.plot(zipf_mete[0], zipf_mete[1], linewidth=2, alpha=0.5, color = 'green', label='METE')
        #ax.plot(zipf_geom[0], zipf_geom[1], linewidth=2, alpha=0.5, color = 'red', label='Broken-stick')

        plt.savefig(figure_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

            #EMPclosed_zipf = IN1['EMPclosed_zipf'].values
            #EMPclosed_mete = IN1['EMPclosed_mete'].values
            #EMPclosed_geom = IN1['EMPclosed_geom'].values
            #print EMPclosed_geom
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


        #elif x == 'HMP':
        #    HMP_zipf = IN1['HMP_zipf'].values
        #    HMP_mete = IN1['HMP_mete'].values
        #    HMP_geom = IN1['HMP_geom'].values
        #    ax.hist(HMP_zipf, 15, linewidth=2, color='blue',alpha=0.5)
        #    ax.hist(HMP_mete, 15, linewidth=2, color='green',alpha=0.5)
        #    ax.hist(HMP_geom, 15, linewidth=2, color='red',alpha=0.5)

        #    plt.savefig(figure_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)

#datasets = [ 'HMP','EMPclosed','MGRAST']
datasets = ['HMP']
#datasets = [ 'EMPclosed']


#methods = ['geom', 'mete','zipf']
methods = ['geom']
#methods = ['zipf',]
#params = ['N','S', 'N/S']
colors = ['red', 'green', 'blue']

#generate_kde_to_file(datasets, methods)
#N_r2_KDE(datasets, methods)
sample_size = 1000
#iterations = 10000
#random_lines(datasets, methods, sample_size, iterations)
#r2_KDE(datasets, methods)
sample_lines(datasets, 10, 1)
