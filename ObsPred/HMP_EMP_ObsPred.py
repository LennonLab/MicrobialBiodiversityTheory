from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

mydir = os.path.expanduser("~/GitHub")
#sys.path.append(mydir)
sys.path.append(mydir + "/macroecotools/")
import macroecotools
import predRADs
import mete


def generate_obs_pred_data(datasets, methods):

    for method in methods:
        for dataset in datasets:

            gN = 0
            #OUT = open(mydir+'/data/'+method+'_'+dataset+'_obs_pred.txt','w+')
            IN = mydir+'/data/'+dataset+'_SADs.txt'
            num_lines = sum(1 for line in open(IN))

            for line in open(IN):

                line = line.split()
                obs = map(int, line)
                obs = list([x for x in obs if x > 1])

                N = sum(obs)
                gN += N
                print N
                S = len(obs)

                if S < 10:
                    continue

                obs.sort()
                obs.reverse()
                print method, dataset, N, S, 'countdown:', num_lines,

                if method == 'geom': # Predicted geometric series
                    pred = predRADs.get_GeomSeries(N, S, False) # False mean no zeros allowed

                elif method == 'mete': # Predicted log-series
                    logSeries = mete.get_mete_rad(S, N)
                    pred = logSeries[0]

                r2 = macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))
                print " r2:", r2

                # write to file, by cite, observed and expected ranked abundances
                #for i, sp in enumerate(pred):
                #    print>> OUT, obs[i], pred[i]

                num_lines -= 1

            print 'N(HMP): ',gN
            #OUT.close()

        print dataset




def transform_obs_data(datasets):

    for dataset in datasets:

        gN = 0
        OUT = open(mydir+'/data/'+dataset+'-data.txt','w+')
        IN = mydir+'/data/'+dataset+'_SADs.txt'
        num_lines = sum(1 for line in open(IN))
        ct = 1

        for line in open(IN):
            line = line.split()
            obs = map(int, line)
            obs = list([x for x in obs if x > 1])
            N = sum(obs)
            S = len(obs)

            if S < 10: continue

            gN += N

            obs.sort()
            obs.reverse()
            for i, sp in enumerate(obs):
                print>>OUT, ct, obs[i]

            print gN, 'countdown:', num_lines - ct
            ct += 1

    OUT.close()
    print 'finished'

    return



def import_obs_pred_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    data = np.genfromtxt(input_filename, dtype = "f8,f8", names = ['obs','pred'], delimiter = " ")
    return data



def hist_mete_r2(sites, obs, pred):  # TAKEN FROM Macroecotools or the mete_sads.py script used for White et al. (2012)
    """Generate a kernel density estimate of the r^2 values for obs-pred plots"""
    r2s = []
    for site in np.unique(sites):
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



def obs_pred_r2_multi(methods, datasets, data_dir='/home/kenlocey/data1/'): # TAKEN FROM THE mete_sads.py script
    print 'generating 1:1 line R-square values for dataset(s)'

    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            obs_pred_data = import_obs_pred_data(data_dir + dataset + '/' + dataset + '_obs_pred.txt')
            obs = ((obs_pred_data["obs"]))
            pred = ((obs_pred_data["pred"]))
            print method, dataset,' ', macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))



def plot_obs_pred_sad(methods, datasets, data_dir='/home/kenlocey/data1/', radius=2): # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    # Used for Figure 3 Locey and White (2013)        ########################################################################################

    """Multiple obs-predicted plotter"""
    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):

            #if method == 'mete' and dataset == 'EMP': continue

            obs_pred_data = import_obs_pred_data(mydir+'/data/truedata/'+method+'_'+dataset+'_obs_pred.txt')
            #site = ((obs_pred_data["site"]))
            obs = ((obs_pred_data["obs"]))
            pred = ((obs_pred_data["pred"]))

            axis_min = 0.5 * min(obs)
            axis_max = 2 * max(obs)

            macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=1,
                            plot_obj=plt.subplot(2, 2, j+1))

            plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
            plt.xlim(axis_min, axis_max)
            plt.ylim(axis_min, axis_max)

            plt.tick_params(axis='both', which='major', labelsize=8)
            plt.subplots_adjust(wspace=0.5, hspace=0.3)

            r2 = macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))
            print method, dataset, r2

            #Create inset for histogram of site level r^2 values
            #axins = inset_axes(ax, width="30%", height="30%", loc=4)
            #hist_mete_r2(site, np.log10(obs), np.log10(pred))
            #plt.setp(axins, xticks=[], yticks=[])

            if method == 'mete': plt.title("Log-series")
            else: plt.title("Geometric series")
            plt.text(1, 2000,  r'$R^2$' + '='+ str(round(r2,3)))
            plt.ylabel('Observed abundance',rotation='90',fontsize=12)
            plt.xlabel('Predicted abundance',fontsize=12)
    plt.savefig(mydir+'/obs_pred_plots.png', dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)



methods = ['mete']
datasets = ['HMP']

#generate_obs_pred_data(datasets, methods)
#plot_obs_pred_sad(methods, datasets)
transform_obs_data(['HMP'])
