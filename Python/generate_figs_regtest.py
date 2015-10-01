#!/usr/bin/env sage -python
from __future__ import division
import sys
import os
import pickle
#sys.path.append(mydir)
mydir = os.path.expanduser("~/github/MicroMETE/data/")
import matplotlib.gridspec as gridspec

import scipy as sp
import  matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import macroecotools
import mete
import macroeco_distributions as md
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
from scipy import stats


"""This code was written using MIT liscenced code from the following Weecology
repos: METE (https://github.com/weecology/METE) and macroecotools
(https://github.com/weecology/macroecotools).
We in no way assume ownership their code"""

def get_SADs_mgrast(path, threshold):
    path_list = []
    path  = path + 'MGRAST-data/' + threshold + '/'
    for subdir, dirs, files in os.walk(path):
        for file in files:
            file_path = os.path.join(subdir, file)
            if file_path.endswith("-data.txt"):
                path_list.append(file_path)
    for x in path_list:
        SADdict = {}
        with open(x) as f:
            for d in f:
                print d
                #print len(d)
                if d.strip():
                    d = d.split()
                    if 'BOVINE' in x:
                        site = d[0]
                        #species = d[1] # Dataset name plus species identifier
                        abundance = float(d[-1])

                    else:
                        site = d[0]
                        #year = d[1]
                        #if closedref == True:
                        #    for i in d:
                        #        if 'unclassified' in i:
                        #            #print 'unclassified'
                        #            continue
                        #        elif 'unidentified' in i:
                        #            #print 'unidentified'
                        #            continue

                        abundance = float(d[-1])


                    if abundance > 0:
                        if site in SADdict:
                            SADdict[site].append(abundance)
                        else:
                            SADdict[site] = [abundance]
    SADs = SADdict.values()
    filteredSADs = []
    for sad in SADs:
        if len(sad) >= 10:
            filteredSADs.append(sad)
    OUT =  open(path + 'MGRAST-' + threshold + '-SADs.txt', 'w')
    #with open(OUT,'wb') as f:
    #    pickle.dump(f,OUT)
    #print >> OUT, filteredSADs
    #return filteredSADs
    SAD_nested_list = list(SADdict.values())
    for SAD in SAD_nested_list:
        print>> OUT, SAD


def get_SADs(path, name, closedref=True):

    SADdict = {}
    DATA = path + name + '-data.txt'

    with open(DATA) as f:

        for d in f:
            print d
            #print len(d)

            if d.strip():
                d = d.split()

                if name == 'GENTRY':
                    site = d[0]
                    #species = d[1] # Dataset name plus species identifier
                    abundance = float(d[-1])

                else:
                    site = d[0]
                    #year = d[1]

                    if closedref == True:
                        for i in d:
                            if 'unclassified' in i:
                                #print 'unclassified'
                                continue
                            elif 'unidentified' in i:
                                #print 'unidentified'
                                continue

                    abundance = float(d[-1])


                if abundance > 0:
                    if site in SADdict:
                        SADdict[site].append(abundance)
                    else:
                        SADdict[site] = [abundance]

    SADs = SADdict.values()
    filteredSADs = []
    for sad in SADs:
        if len(sad) >= 10:
            filteredSADs.append(sad)


    return filteredSADs




def EMP_SADs(path, name, mgrast):

    minS = 10

    IN = path + '/' + name + '-SSADdata.txt'
    n = sum(1 for line in open(IN))

    SiteDict = {}

    with open(IN) as f:

        for d in f:

            n -= 1

            if d.strip():

                d = d.split()
                #species = d[0]
                sample = d[1]
                abundance = float(d[2])

                if abundance > 0:
                    if sample not in SiteDict:

                        SiteDict[sample] = [abundance]

                    else:
                        SiteDict[sample].append(abundance)


    SADs = SiteDict.values()
    filteredSADs = []
    for sad in SADs:
        if len(sad) >= minS:
            filteredSADs.append(sad)

    return filteredSADs







def get_GeomSeries(N,S,zeros):

    rank = range(1,S+1)
    cdf = [(S-i+0.5)/S for i in rank]
    SNratio = S/N
    if zeros == False:
        abd = md.trunc_geom.ppf(np.array(cdf), SNratio, N)
    return abd

def generate_obs_pred_data(datasets, methods, size):

    for method in methods:
        for dataset in datasets:
            #OUT1 = open(mydir + "ObsPred/" + method +'_'+dataset+'_obs_pred.txt','w+')
            #OUT2 = open(mydir + "NSR2/" + method +'_'+dataset+'_NSR2.txt','w+')
            #OUT1 = open(mydir + "ObsPred/" + method +'_'+dataset+'_obs_pred_subset.txt','w+')
            #OUT2 = open(mydir + "NSR2/" + method +'_'+dataset+'_NSR2_subset.txt','w+')

            if dataset == "HMP":
                IN = mydir  + dataset + '-Data' + '/' + dataset +'-SADs.txt'
                num_lines = sum(1 for line in open(IN))
                OUT1 = open(mydir + "ObsPred/" + method +'_'+dataset+'_obs_pred.txt','w+')
                OUT2 = open(mydir + "NSR2/" + method +'_'+dataset+'_NSR2.txt','w+')
            elif dataset == 'EMPclosed' or dataset == 'EMPpen':
                IN = mydir  + dataset + '-Data' + '/' + dataset +'-SADs.txt'
                num_lines = sum(1 for line in open(IN))
                random_sites = np.random.randint(num_lines,size=size)
                num_lines = size
                OUT1 = open(mydir + "ObsPred/" + method +'_'+dataset+'_obs_pred_subset.txt','w+')
                OUT2 = open(mydir + "NSR2/" + method +'_'+dataset+'_NSR2_subset.txt','w+')
                num_lines = sum(1 for line in open(IN))
            else:
                IN = mydir + 'MGRAST-Data/' + dataset +  '/' + 'MGRAST-' + dataset + '-SADs.txt'
                num_lines = sum(1 for line in open(IN))
                OUT1 = open(mydir + "ObsPred/" + method +'_'+ 'MGRAST' + dataset+'_obs_pred.txt','w+')
                OUT2 = open(mydir + "NSR2/" + method +'_'+ 'MGRAST' + dataset+'_NSR2.txt','w+')

            for j,line in enumerate(open(IN)):
                if dataset == "HMP":
                    line = line.split()
                elif size == 0:
                    line = eval(line)
                else:
                    line = eval(line)
                    if j not in random_sites:
                        continue
                #line.strip("[]")
                #line.split()
                obs = map(int, line)

                N = sum(obs)
                S = len(obs)

                if S < 10 or N <= S:
                    num_lines += 1
                    continue

                obs.sort()
                obs.reverse()
                print method, dataset, N, S, 'countdown:', num_lines,

                if method == 'geom': # Predicted geometric series
                    pred = get_GeomSeries(N, S, False) # False mean no zeros allowed

                elif method == 'mete': # Predicted log-series
                    logSeries = mete.get_mete_rad(S, N)
                    pred = logSeries[0]

                r2 = macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))
                print " r2:", r2
                if r2 == -float('inf') or r2 == float('inf') or r2 == float('Nan'):
                    print r2 + " is Nan or inf, removing..."
                    continue
                print>> OUT2, j, N, S, r2
                # write to file, by cite, observed and expected ranked abundances
                for i, sp in enumerate(pred):
                    print>> OUT1, j, obs[i], pred[i]


                num_lines -= 1

            OUT1.close()

        print dataset



def import_obs_pred_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8", names = ['site','obs','pred'], delimiter = " ")
    #test = data[0:10000]
    #return test
    return data


def hist_mete_r2(sites, obs, pred):  # TAKEN FROM Macroecotools or the mete_sads.py script used for White et al. (2012)
    """Generate a kernel density estimate of the r^2 values for obs-pred plots"""
    r2s = []
    for site in np.unique(sites):
        if int(site) >= 100:
            break
        else:
            obs_site = obs[sites==site]
            pred_site = pred[sites==site]
            r2 = macroecotools.obs_pred_rsquare(obs_site, pred_site)
            r2s.append(r2)

    hist_r2 = np.histogram(r2s, range=(0, 1))
    xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
    xvals = xvals[0:len(xvals)-1]
    yvals = hist_r2[0]
    plt.plot(xvals, yvals, 'k-', linewidth=2)
    plt.axis([0, 1, 0, 1.1 * max(yvals)])

def obs_pred_r2_multi(methods, datasets, data_dir= mydir): # TAKEN FROM THE mete_sads.py script
    print 'generating 1:1 line R-square values for dataset(s)'

    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + "_" + dataset + '_obs_pred.txt')
            #obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + "_" + dataset + '_obs_pred.txt')
            obs = ((obs_pred_data["obs"]))
            pred = ((obs_pred_data["pred"]))
            print method, dataset,' ',macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))

def import_NSR2_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8", names = ['site','N','S', 'R2'], delimiter = " ")
    #test = data[0:5000]
    #return test
    return data

def plot_obs_pred_sad(methods, datasets, data_dir= mydir, radius=2): # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    # Used for Figure 3 Locey and White (2013)        ########################################################################################
    """Multiple obs-predicted plotter"""
    fig = plt.figure()


    count = 0


    #ax = fig.add_subplot(111)
    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            #if method == 'mete' and dataset == 'EMP': continue
            #obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method+'_'+dataset+'_obs_pred_test.txt')
            if str(dataset) == 'EMPclosed' or str(dataset) == 'EMPopen':
                obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method+'_'+dataset+'_obs_pred_subset.txt')
            elif str(dataset) == 'HMP':
                obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method+'_'+dataset+'_obs_pred.txt')
            else:
                obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_' + 'MGRAST' + dataset +'_obs_pred.txt')
            print method, dataset
            site = ((obs_pred_data["site"]))
            obs = ((obs_pred_data["obs"]))
            pred = ((obs_pred_data["pred"]))
            axis_min = 0.5 * min(obs)
            axis_max = 2 * max(obs)
            ax = fig.add_subplot(3, 3, count+1)
            ax.set(adjustable='box-forced', aspect='equal')
            if j == 0:
                if i == 0:
                    ax.set_ylabel("Human Microbiome Project" + "\n" + "(HMP)", rotation=90, size=8)
                elif i == 1:
                    ax.set_ylabel("Earth Microbiome Project" + "\n" + "(EMP)", rotation=90, size=8)
                elif i == 2:
                    ax.set_ylabel("MG-RAST", rotation=90, size=8)
                #elif i == 3:
                #    ax.set_ylabel('MG-RAST', rotation=90, size=12)
            if i == 0 and j == 0:
                ax.set_title("Broken-stick")
            elif i == 0 and j == 1:
                ax.set_title("METE")

            macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=1,
                            plot_obj=plt.subplot(3,3,count+1))

            plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
            plt.xlim(axis_min, axis_max)
            plt.ylim(axis_min, axis_max)

            plt.tick_params(axis='both', which='major', labelsize=8)
            plt.subplots_adjust(wspace=0.5, hspace=0.3)

            #plt.text(xs[0][1],xs[0][0],dataset+'\n'+rs[0],fontsize=8)
            #xs.pop(0)
            #rs.pop(0)
            # Create inset for histogram of site level r^2 values
            axins = inset_axes(ax, width="30%", height="30%", loc=4)
            if str(dataset) == 'EMPclosed' or str(dataset) == 'EMPopen':
                INh2 = import_NSR2_data(data_dir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')
                r2s = ((INh2["R2"]))
                hist_r2 = np.histogram(r2s, range=(0, 1))
                xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
                xvals = xvals[0:len(xvals)-1]
                yvals = hist_r2[0]
                plt.plot(xvals, yvals, 'k-', linewidth=2)
                plt.axis([0, 1, 0, 1.1 * max(yvals)])
            else:
                hist_mete_r2(site, np.log10(obs), np.log10(pred))
            plt.setp(axins, xticks=[], yticks=[])
            count += 1
        count += 1
    #ax.set_xlabel(-8,-80,'Rank-abundance at the centre of the feasible set',fontsize=10)
    #ax.set_ylabel(-8.5,500,'Observed rank-abundance',rotation='90',fontsize=10)
    #ax.set_ylabel('Rank-abundance at the centre of the feasible set',rotation='90',fontsize=10)
    #fig.tight_layout()
    #fig.subplots_adjust(hspace=.5)
    fig.subplots_adjust(wspace=0.0001, left=0.1)


    fig.text(0.05, 0.5, 'Observed rank-abundance', ha='center', va='center', rotation='vertical')
    #fig.text(0.35, 0.04, 'Predicted rank-abundance', ha='center', va='center')
    fig.text(0.37, 0.04, 'Predicted rank-abundance', ha='center', va='center')
    #ax.set_xlabel('Observed rank-abundance',fontsize=10)
    plt.savefig('obs_pred_plots.png', dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()



# Make a function to generate the histogram.
def NSR2_regression(methods, datasets, data_dir= mydir):
    fig = plt.figure()
    count  = 0
    test_count = 0
    #st = fig.suptitle("Broken-stick", fontsize="x-large")
    #fig.text(0.02, 0.5, r'$r^{2}$', ha='center', va='center', rotation='vertical', size = 'x-large')
    #fig.text(0.04, 0.5, 'common Y', va='center', rotation='vertical')
    for i, dataset in enumerate(datasets):
        for k, param in enumerate(params):
            for j, method in enumerate(methods):
                nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')

                #nsr2_data[[~np.isnan(nsr2_data).any(axis=1)]]
                #nsr2_data[~np.isinf(nsr2_data).any(axis=1)]
                #nsr2_data[~np.isnan(nsr2_data).any(1)]
                list = ['nan', 'NAN', '-inf', 'inf']
                #for x in nsr2_data:
                #    print type(x)
                #    value = str(x[3])

                #if np.isinf(x[3]) == True:
                #    print "infinity"
                #mask = np.all(np.isinf(nsr2_data), axis=1)


                y = ((nsr2_data["R2"]))
                mean = np.mean(y)
                std_error = sp.stats.sem(y)
                print method, param, dataset
                #print "mean = " + str(mean)
                #print "standard error = " + str(std_error)

                #print method, dataset, param
                ax = fig.add_subplot(3, 2, count+1)
                if param == "N" or param == "S":
                    x = np.log10(((nsr2_data[param])))
                else:
                    N_count = ((nsr2_data["N"]))
                    S_count = ((nsr2_data["S"]))
                    print dataset, method
                    print "mean N is " + str(np.mean(N_count))
                    print "mean S is " + str(np.mean(S_count))
                    x = np.divide(N_count, S_count)
                    x = np.log10(x)

                #elif str(n_or_s).capitalize() == 'S'
                #x = ((nsr2_data["N"]))
                macroecotools.plot_color_by_pt_dens(x, y, 0.1, loglog=0,
                                plot_obj=plt.subplot(3, 2, count+1))
                slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
                #if param == 'N/S':
                #    plt.xlim(np.amin(x), 1000)
                #else:
                plt.xlim(np.amin(x), np.amax(x))
                plt.ylim(-1,1)

                predict_y = intercept + slope * x
                pred_error = y - predict_y
                degrees_of_freedom = len(x) - 2
                residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
                plt.plot(x, predict_y, 'k-')
                plt.axhline(linewidth=2, color='darkgrey',ls='--')
                #plt.hline(0, xmin, xmax, color="0.3", ls='--')
                plt.subplots_adjust(wspace=0.2, hspace=0.3)
                # Plotting
                plt.xlabel(param)
                if k == 1 and j ==0:
                    plt.ylabel(r'$r^{2}_{m}$', fontsize = 'xx-large')
                #if k == 1 and i ==0:
                #    plt.ylabel(r'$r^{2}_{m}$', fontsize = 'xx-large')
                #if k == 0 and i ==0 :
                #    plt.title('HMP', fontsize = 'large')
                #elif k == 0 and i ==1:
                #    plt.title('EMP Closed', fontsize = 'large')
                #elif k == 0 and i ==2:
                #    plt.title('EMP Open', fontsize = 'large')
                #elif k == 0 and i ==3:
                #    plt.title('MG-RAST, 97%', fontsize = 'large')
                if j == 0 and k == 0:
                    plt.title('Broken-stick', fontsize = 'large')
                elif j == 1 and k == 0:
                    plt.title('METE', fontsize = 'large')

                #plt.ylabel(r'$r^{2}$',fontsize=16)
                #r_2 = "r2 =" + str(round(r_value,2))
                #p_s = "p =" + str(round(p_value,2))
                #plt.text(0, 1, r'$p$'+ ' = '+str(round(p_value,2)), fontsize=12)
                #plt.text(0, 1, r'$r_{2}$'+ ' = '+str(round(r_value,2)), fontsize=12)
                #ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
                leg = plt.legend(loc=1,prop={'size':10})
                ax.tick_params(axis='x', labelsize=6)
                ax.tick_params(axis='y', labelsize=6)
                #tick.label.set_fontsize(14)
                #leg.draw_frame(False)
                #plt.legend(loc='upper left')
                print r_value, p_value

                count += 1
    #ax.set_ylabel('common ylabel')
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

    #fig.text(0.02, 0.5, r'$r^{2}$', ha='center', va='center', rotation='vertical', size = 'large')
    #st.set_y(0.95)
    fig.subplots_adjust(top=0.85)
    #ax.set_ylabel('common ylabel')
    #fig.text(-8,-80,'Rank-abundance at the centre of the feasible set',fontsize=10)
    #plt.suptitle(-8.5,500,r'$r^{2}$',rotation='90',fontsize=10)
    fig_name = 'NSR2_GeomMete' + str(dataset) + '.png'
    plt.savefig(fig_name)
    #plt.xscale()
    plt.close()

methods = ['geom', 'mete']
#methods = ['geom']
datasets = ['HMP', 'EMPclosed', '97']
#datasets = ['EMPclosed', 'EMPopen']
#datasets = ['95', '97','99']
#datasets = ['HMP']

params = ['N','S', 'N/S']
#params = ['N/S']
#get_SADs()
#generate_obs_pred_data(datasets, methods, 0)
#empclosed = ['EMPclosed']
#generate_obs_pred_data(empclosed, geommethods, 500)
plot_obs_pred_sad(methods, datasets)
#NSR2_regression(methods, datasets, data_dir= mydir)

#get_SADs_mgrast(mydir, '99')
