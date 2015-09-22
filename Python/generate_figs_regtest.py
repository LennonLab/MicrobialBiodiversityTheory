#!/usr/bin/env sage -python
from __future__ import division
import sys
import os

#mydir = os.path.expanduser("~/Desktop/Repos/HMP_EMP/")
#sys.path.append(mydir)
mydir = os.path.expanduser("~/github/MicroMETE/data/")

import  matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid.inset_locator import inset_axes

import macroecotools
import mete
import macroeco_distributions as md
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
from scipy import stats
def get_GeomSeries(N,S,zeros):

    rank = range(1,S+1)
    cdf = [(S-i+0.5)/S for i in rank]
    SNratio = S/N
    if zeros == False:
        abd = md.trunc_geom.ppf(np.array(cdf), SNratio, N)

    return abd

def generate_obs_pred_data(datasets, methods):

    for method in methods:
        for dataset in datasets:

            #OUT1 = open(mydir + "ObsPred/" + method +'_'+dataset+'_obs_pred.txt','w+')
            #OUT2 = open(mydir + "NSR2/" + method +'_'+dataset+'_NSR2.txt','w+')
            OUT1 = open(mydir + "ObsPred/" + method +'_'+dataset+'_obs_pred_subset.txt','w+')
            OUT2 = open(mydir + "NSR2/" + method +'_'+dataset+'_NSR2_subset.txt','w+')
            IN = mydir  + dataset + '-Data' + '/' + dataset +'-SADs.txt'
            num_lines = sum(1 for line in open(IN))
            # randomly pick 5000 linres
            random_sites = np.random.randint(300,size=num_lines)
            num_lines = 300
            for j,line in enumerate(open(IN)):
                if dataset == "HMP" :
                    line = line.split()
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



def plot_obs_pred_sad(methods, datasets, data_dir= mydir, radius=2): # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    # Used for Figure 3 Locey and White (2013)        ########################################################################################
    """Multiple obs-predicted plotter"""
    fig = plt.figure()

    #xs = [[60,1], [100,1], [20,1], [60,1], [40,1], [200,1], [800,1.5], [200,1.5]]
    #rs = ['0.93','0.77','0.84','0.81','0.78','0.83','0.58','0.76']
    count = 0
    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):

            #if method == 'mete' and dataset == 'EMP': continue

            print method, dataset
            if method == 'EMPclosed' or method == 'EMPopen':
                obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method+'_'+dataset+'_obs_pred_subset.txt')
            else:
                obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method+'_'+dataset+'_obs_pred.txt')
            print method, dataset
            site = ((obs_pred_data["site"]))
            obs = ((obs_pred_data["obs"]))
            pred = ((obs_pred_data["pred"]))
            #if method == 'EMPclosed' or method == 'EMPopen':
            #    unique = np.unique(((obs_pred_data["site"])))
            #    site_num = unique.size
            #    random_sites = np.random.randint(5000,size=site_num)
            #    obs_pred_data = obs_pred_data[random_sites, :]
            #    site = ((obs_pred_data["site"]))
            #    obs = ((obs_pred_data["obs"]))
            #    pred = ((obs_pred_data["pred"]))
            #else:
            #    site = ((obs_pred_data["site"]))
            #    obs = ((obs_pred_data["obs"]))
            #    pred = ((obs_pred_data["pred"]))
            axis_min = 0.5 * min(obs)
            axis_max = 2 * max(obs)
            ax = fig.add_subplot(3, 2, count+1)

            macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=1,
                            plot_obj=plt.subplot(3,2,count+1))

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
            hist_mete_r2(site, np.log10(obs), np.log10(pred))
            plt.setp(axins, xticks=[], yticks=[])
            count += 1
    plt.text(-8,-80,'Rank-abundance at the centre of the feasible set',fontsize=10)
    plt.text(-8.5,500,'Observed rank-abundance',rotation='90',fontsize=10)
    plt.savefig('obs_pred_plots.png', dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()


def import_NSR2_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8", names = ['site','N','S', 'R2'], delimiter = " ")
    #test = data[0:5000]
    #return test
    return data

# Make a function to generate the histogram.
def NSR2_regression(methods, datasets, data_dir= mydir):
    fig = plt.figure()
    for i, dataset in enumerate(datasets):
        count  = 0
        test_count = 0
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
                mask = np.all(np.isinf(nsr2_data), axis=1)


                y = ((nsr2_data["R2"]))
                print method, dataset, param
                ax = fig.add_subplot(2, 3, count+1)
                if param == "N" or param == "S":
                    x = np.log10(((nsr2_data[param])))
                else:
                    N_count = ((nsr2_data["N"]))
                    S_count = ((nsr2_data["S"]))
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
                plt.subplots_adjust(wspace=0.5, hspace=0.3)
                # Plotting
                #plt.text(-8,-80,'Rank-abundance at the centre of the feasible set',fontsize=10)
                #plt.text(-8.5,500,'Observed rank-abundance',rotation='90',fontsize=10)
                plt.xlabel(param)
                plt.ylabel(r'$r^{2}$',fontsize=16)
                r_2 = "r2 =" + str(round(r_value,2))
                p_s = "p =" + str(round(p_value,2))
                #plt.text(0, 1, r'$p$'+ ' = '+str(round(p_value,2)), fontsize=12)
                #plt.text(0, 1, r'$r_{2}$'+ ' = '+str(round(r_value,2)), fontsize=12)
                #ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14,
    #verticalalignment='top', bbox=props)
                leg = plt.legend(loc=1,prop={'size':10})
                #leg.draw_frame(False)
                #plt.legend(loc='upper left')
                print r_value, p_value
                count += 1
        plt.tight_layout()
        fig_name = 'NSR2_GeomMete' + str(dataset) + '.png'
        plt.savefig(fig_name)
        #plt.xscale()
        plt.close()

#methods = ['geom', 'mete']
methods = ['geom']
#datasets = ['HMP', 'EMPclosed', 'EMPopen']
datasets = ['EMPclosed']
params = ['N','S', 'N/S']
generate_obs_pred_data(datasets, methods)
plot_obs_pred_sad(methods, datasets)
#NSR2_regression(methods, datasets, data_dir= mydir)
