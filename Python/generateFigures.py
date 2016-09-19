from __future__ import division
import imp, os, signal, datetime, random
import importData
import macroecotools
import mete
import scipy as sp
from scipy import stats, optimize
import pandas as pd
import numpy as np
import  matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import modelsAndMetrics as mo
import itertools
import matplotlib.ticker as mticker

mydir = os.path.expanduser("~/github/MicroMETE/")

def fig1(figname = 'Fig1', data_dir= mydir, saveAs = 'eps'):
    SAD = [10000, 8000, 6000, 5000, 1000, 200, 100,  20, 18, 16, 14, 12, 10, 4,5,
        4, 4, 3, 3, 2, 2, 2, 2, 2,2, 1, 1, 1, 1, 1,1,1,1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    SAD.sort()
    SAD.reverse()
    x = range(1, len(SAD) +1)
    N = sum(SAD)
    S = len(SAD)

    geom = np.log10(mo.get_GeomSeries(N, S, False))

    logSeries = np.log10(mete.get_mete_rad(S, N)[0])


    lognorm_pred = mo.lognorm(SAD, 'pln')
    lognorm_SAD = np.log10(lognorm_pred.get_rad_from_obs())
    zipf_class = mo.zipf(SAD, 'fmin')
    pred_tuple = zipf_class.from_cdf()
    zipf_SAD = np.log10(pred_tuple[0])
    gamma = pred_tuple[1]

    SAD = np.log10(SAD)
    fig = plt.figure()
    plt.plot()

    max_y = max(max(SAD),  max(zipf_SAD))

    plt.plot(x, SAD,color = '#A9A9A9', linestyle = '-', linewidth=2, label="Observed")
    plt.plot(x, geom,color = '#00008B', linestyle = '-', linewidth=2, label="Broken-stick")
    plt.plot(x, lognorm_SAD, color = '#0000CD',linestyle = '--', linewidth=2, label="Lognormal")
    plt.plot(x, logSeries, color = '#FF4500',linestyle = '-.', linewidth=2, label="Log-series")
    plt.plot(x, zipf_SAD, color = 'red',linestyle = '-',linewidth=2,  label="Zipf")

    plt.tight_layout()
    plt.xlabel('Rank Abundance', fontsize = 22)
    plt.ylabel('Abundance, ' +r'$log_{10}$', fontsize = 22)
    output = "dorm_fix_prob.png"
    plt.legend(loc='upper right')

    #plt.yscale('log')
    #plt.yscale('log')
    plt.xlim(1, len(SAD))
    plt.ylim(-0.25 , max_y)

    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.legend(frameon=False, fontsize= 18)

    fig_name = str(mydir + 'figures/' + figname + '_RGB.' + saveAs)
    plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600, \
        format = saveAs)
    plt.close()

def fig2(n=352899, figname = 'Fig2', data_dir=mydir, \
    stratify = True, radius=2, remove = 0, zipfType = 'mle', RGF = False, \
    lognormType = 'pln', saveAs = 'eps'):
    # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    # Used for Figure 3 Locey and White (2013)
    """Multiple obs-predicted plotter"""
    fig = plt.figure()
    count = 0
    plot_dim = 2
    methods = ['geom', 'lognorm', 'mete', 'zipf']
    fig.subplots_adjust(bottom= 0.30)
    for i, method in enumerate(methods):
        if method == 'zipf':
            obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/'+ method + '_'+  zipfType+'_obs_pred_stratify.txt')
            INh2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified/' + method + '_mle' + '_NSR2_stratify.txt')
        elif method == 'lognorm':
            obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/'+ method + '_'+  lognormType+'_obs_pred_stratify.txt')
            INh2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified/' + method + '_'+  lognormType + '_NSR2_stratify.txt')
        else:
            obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/'+ method +'_obs_pred_stratify.txt')
            INh2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified/' + method + '_NSR2_stratify.txt')
        obs = np.asarray(list(((obs_pred_data["obs"]))))
        pred = np.asarray(list(((obs_pred_data["pred"]))))
        site = np.asarray(list(((obs_pred_data["site"]))))

        obs2 = []
        pred2 = []
        site2 = []

        obs_all = np.asarray(obs)
        pred_all = np.asarray(pred)
        site_all = np.asarray(site)

        if n == 'all' or len(obs) <= n:
            obs2 = list(obs)
            pred2 = list(pred)
            site2 = list(site)

        else:
            if len(obs) > n:
                inds = np.random.choice(range(len(site)), size=n, replace=False)
                for ind in inds:
                    obs2.append(obs[ind])
                    pred2.append(pred[ind])
                    site2.append(site[ind])

        obs = np.asarray(obs2)
        pred = np.asarray(pred2)
        site =  np.asarray(site2)

        if method == 'zipf':
            axis_min = 0
            axis_max = 2  * max(pred)
        else:
            axis_min = 0
            axis_max = 2 * max(obs)
        ax = fig.add_subplot(plot_dim, plot_dim, count+1)
        if method == 'zipf':
            NSR2_BS = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/'+ method  + '_mle_NSR2_stratify.txt')
        elif method == 'lognorm':
            NSR2_BS = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/'+ method  + '_pln_NSR2_stratify.txt')
        else:
            NSR2_BS = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/'+ method  +'_NSR2_stratify.txt')

        if method == 'geom':
            ax.set_title("Broken-stick")
        elif method == 'lognorm':
            ax.set_title("Lognormal")
        elif method == 'mete':
            ax.set_title("Log-series")
        elif method == 'zipf':
            ax.set_title("Zipf")
        print len(pred), len(obs)
        macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=1,
                        plot_obj=plt.subplot(plot_dim,plot_dim,count+1))

        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
        if method == 'zipf':
            plt.xlim(0, axis_max)
            plt.ylim(0, axis_max)
        else:
            plt.xlim(0, axis_max)
            plt.ylim(0, axis_max)
        r2s = ((INh2["R2"]))
        r2s = r2s.astype(float)
        # insert r2 of all data
        r2_all = np.mean(((NSR2_BS["R2"])))
        print method + ' mean = ' + str(r2_all)
        print method + ' std dev = ' +str(np.std(r2_all))
        r2text = r"${}^{{2}}_{{m}} = {:.{p}f} $".format('r',r2_all , p=2)
        if method == 'geom':
            plt.text(0.25, 0.90, r2text,  fontsize=14,
                horizontalalignment='center',
                verticalalignment='center',transform = ax.transAxes)
        else:
            plt.text(0.22, 0.90, r2text,  fontsize=14,
                horizontalalignment='center',
                verticalalignment='center',transform = ax.transAxes)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.subplots_adjust(wspace=0.0000000001, hspace=0.5)

        axins = inset_axes(ax, width="30%", height="30%", loc=4)

        hist_r2 = np.histogram(r2s, range=(0, 1))
        xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
        xvals = xvals[0:len(xvals)-1]
        yvals = hist_r2[0]
        plt.plot(xvals, yvals, 'k-', linewidth=2)
        plt.axis([0, 1, 0, 1.1 * max(yvals)])

        ax.set(adjustable='box-forced', aspect='equal')
        plt.setp(axins, xticks=[], yticks=[])


        count += 1
    plt.tight_layout(pad=1.5, w_pad=0.8, h_pad=0.8)
    fig.text(0.50, 0.03, 'Predicted abundance', ha='center', va='center', fontsize=16)
    fig.text(0.08, 0.5, 'Observed abundance', ha='center', va='center', rotation='vertical', fontsize=16)
    fig_name = str(mydir + 'figures/' + figname + '_RGB.' + saveAs)
    plt.savefig(fig_name, dpi=600, format = saveAs)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()


def fig3(figname = 'Fig3', \
    zipfType = 'mle', lognormType = 'pln', Stratified = True, data_dir= mydir, \
    saveAs = 'eps'):
    methods = ['geom', 'lognorm', 'mete', 'zipf']
    fig = plt.figure()
    count  = 0
    params = ['N']
    removeSADs = []
    for i, param in enumerate(params):
        for j, method in enumerate(methods):
            if method == 'zipf':
                obs_pred_data = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified/'+ method + '_'+ zipfType + '_NSR2_stratify.txt')
            elif method == 'lognorm':
                obs_pred_data = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified/'+ method + '_'+ lognormType +'_NSR2_stratify.txt')
            else:
                obs_pred_data = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified/'+ method +'_NSR2_stratify.txt')
            site = np.asarray(list(((obs_pred_data["site"]))))
            y = np.asarray(list(((obs_pred_data["R2"]))))
            x = np.log10(np.asarray(list(((obs_pred_data[param])))))
            print "nmax" + str(np.mean(np.asarray(list(((obs_pred_data["NmaxObs"]))))))

            mean_x = np.mean(x)
            mean_y = np.mean(y)
            std_error = sp.stats.sem(y)
            print method, param
            print "mean modified r2 = " + str(mean_y)
            print "modified r2 standard error = " + str(std_error)
            print "mean " + param  + " is " + str(np.mean(np.asarray(list(((obs_pred_data[param]))))))
            ax = fig.add_subplot(2, 2, count+1)
            macroecotools.plot_color_by_pt_dens(x, y, 0.1, loglog=0,
                            plot_obj=plt.subplot(2, 2, count+1))
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            print "slope is " + str(slope)
            print "r2-value is " + str(r_value **2)
            print "p-value is " + str(p_value)

            print "NmaxPred ", method
            NmaxPred = np.log10(np.asarray(list(((obs_pred_data["NmaxPred"])))))
            slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(x,NmaxPred)

            diff1 = slope1  - 1
            p_diff1 = (np.absolute(diff1) /  ( 0.5 *(slope1 + 1) )) * 100
            print "percent difference " + str(p_diff1)

            print "evennessPred ",method
            evennessPred = np.log10(np.asarray(list(((obs_pred_data["evennessPred"])))))
            slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(x,evennessPred)
            diff2 = slope2  - (-0.31)
            p_diff2 = (np.absolute(diff2)  /  ( 0.5 *(slope1 + 1) )) * 100


            print "percent difference " + str(p_diff2)

            print "skewnessPred ",method
            skewnessPred = np.log10(np.asarray(list(((obs_pred_data["skewnessPred"])))))
            slope3, intercept3, r_value3, p_value3, std_err3 = stats.linregress(x,skewnessPred)
            diff3 = slope3  - 0.13
            p_diff3 = (np.absolute(diff3)  /  ( 0.5 *(slope1 + 1) )) * 100

            print "percent difference " + str(p_diff3)
            plt.xlim(np.amin(x), np.amax(x))

            plt.ylim(-1.5,1.5)
            predict_y = intercept + slope * x
            pred_error = y - predict_y
            degrees_of_freedom = len(x) - 2
            residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
            plt.plot(x, predict_y, 'k-')
            plt.axhline(linewidth=2, color='darkgrey',ls='--')
            plt.tick_params(axis='both', which='major', labelsize=10)

            if i == 0 and j == 0:
                ax.set_title("Broken-stick", fontsize = 15)
                ax.set_ylabel(r'$r^{2}_{m}$', fontsize = 22)
            if i == 0 and j == 1:
                ax.set_title("Lognormal", fontsize = 15)
            if i == 0 and j == 2:
                ax.set_title("Log-series", fontsize = 15)
                ax.set_xlabel('Abundance, ' +r'$log_{10}$', fontsize = 17)
                ax.set_ylabel(r'$r^{2}_{m}$', fontsize = 22)
            if i == 0 and j == 3:
                ax.set_title("Zipf", fontsize = 15)
                ax.set_xlabel('Abundance, ' +r'$log_{10}$', fontsize = 17)
            ax.tick_params(axis='x', labelsize=12)
            ax.tick_params(axis='y', labelsize=12)

            count += 1

    plt.tight_layout(pad=0.8, w_pad=0.8, h_pad=0.8)
    fig_name = str(mydir + 'figures/' + figname  + '_RGB.' + saveAs)
    plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600, \
        format = saveAs)
    plt.close()



def fig4(figname = 'Fig4', data_dir=mydir, radius=2, saveAs = 'eps'):
    fig = plt.figure()
    fig.subplots_adjust(bottom= 0.15)
    plot_dim = 1
    count = 0

    IN_Obs_Pred = importData.import_NSR2_data(mydir + \
        'data/NSR2/Stratified/lognorm_pln_NSR2_stratify.txt')
    N = np.asarray(list(((IN_Obs_Pred["N"]))))
    S = np.asarray(list(((IN_Obs_Pred["S"]))))
    NmaxObs = np.asarray(list(((IN_Obs_Pred["NmaxObs"]))))
    models = ['geom', 'lognorm', 'mete', 'zipf']
    modelSlopes = [0.647520323289, 0.942904468437, 0.769214774397, 0.954497727096]
    modelInterepts = [0.116508916992, 0.292527611072, 0.19240314275, 0.189954627996]
    for g, model in enumerate(models):
        NmaxPred = []
        SPred = []
        for i in range(len(N)):
            NmaxPred_i = mo.predictS(N[i], NmaxObs[i], \
                predictNmax=True).getNmax(b = modelInterepts[g], slope = modelSlopes[g])
            SPred_i = mo.predictS(N[i], NmaxObs[i], predictNmax=True).getS()
            NmaxPred.append(NmaxPred_i)
            SPred.append(SPred_i)
        NmaxPred = np.asarray(NmaxPred)
        SPred = np.asarray(SPred)
        axis_min = 0
        axis_max = 2 * max(NmaxObs)
        ax = fig.add_subplot(2, 2, count+1)
        if model == 'zipf':
            OUT2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/'+ model  + '_mle_NSR2_stratify.txt')
        elif model == 'lognorm':
            OUT2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/'+ model  + '_pln_NSR2_stratify.txt')
        else:
            OUT2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/'+ model  +'_NSR2_stratify.txt')

        NmaxObs_BS = np.asarray(list(((OUT2["NmaxObs"]))))
        NmaxPred_BS = np.asarray(list(((OUT2["NmaxPred"]))))

        if model == 'geom':
            ax.set_title("Broken-stick")
        elif model == 'lognorm':
            ax.set_title("Lognormal")
        elif model == 'mete':
            ax.set_title("Log-series")
        elif model == 'zipf':
            ax.set_title("Zipf")
        macroecotools.plot_color_by_pt_dens(NmaxPred, NmaxObs, radius, loglog=1,
                        plot_obj=plt.subplot(2,2,count+1))
        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
        plt.xlim(axis_min, axis_max)
        plt.ylim(0, axis_max)
        r2_all = macroecotools.obs_pred_rsquare(np.log10(NmaxObs), np.log10(NmaxPred))
        r2text = r"${}^{{2}}_{{m}} = {:.{p}f} $".format('r',r2_all , p=2)
        plt.text(0.72, 0.12, r2text,  fontsize=13,
            horizontalalignment='center',
            verticalalignment='center',transform = ax.transAxes)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.subplots_adjust(wspace=0.00001, hspace=0.3)
        ax.set(adjustable='box-forced', aspect='equal')

        count += 1
    fig.text(0.50, 0.055 , 'Predicted, ' +r'$log_{10}(N_{max})$', ha='center', va='center', fontsize = 19)
    fig.text(0.09, 0.5, 'Observed, ' +r'$log_{10}(N_{max})$', ha='center', va='center', rotation='vertical',\
        fontsize = 19)
    fig_name = str(mydir + 'figures/' + figname + '_RGB.' + saveAs)
    plt.savefig(fig_name, dpi=600, format = saveAs)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()

def figS1(n=35289, figname = 'FigS1', data_dir=mydir, radius=2, zipfType = 'mle', \
    saveAs = 'eps', lognormType = 'pln'):
    methods = ['geom', 'lognorm', 'mete', 'zipf']
    datasets = ['95', '97', '99']
    fig = plt.figure()
    count = 0
    rows = len(datasets)
    columns = len(methods)
    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            print count
            if method == 'zipf':
                obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/'+ method + '_'+  zipfType+'_obs_pred_stratify.txt')
                INh2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified/' + method + '_mle' + '_NSR2_stratify.txt')
            elif method == 'lognorm':
                obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/'+ method + '_'+  lognormType+'_obs_pred_stratify.txt')
                INh2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified/' + method + '_'+  lognormType + '_NSR2_stratify.txt')
            else:
                obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/'+ method +'_obs_pred_stratify.txt')
                INh2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified/' + method + '_NSR2_stratify.txt')
            obs = np.asarray(list(((obs_pred_data["obs"]))))
            pred = np.asarray(list(((obs_pred_data["pred"]))))
            site = np.asarray(list(((obs_pred_data["site"]))))


            obs2 = []
            pred2 = []
            site2 = []

            obs_all = np.asarray(obs)
            pred_all = np.asarray(pred)
            site_all = np.asarray(site)

            if n == 'all' or len(obs) <= n:
                obs2 = list(obs)
                pred2 = list(pred)
                site2 = list(site)

            else:
                if len(obs) > n:
                    inds = np.random.choice(range(len(site)), size=n, replace=False)
                    for ind in inds:
                        obs2.append(obs[ind])
                        pred2.append(pred[ind])
                        site2.append(site[ind])

            obs = np.asarray(obs2)
            pred = np.asarray(pred2)
            site =  np.asarray(site2)

            print "number of points " + str(len(obs))
            if method == 'zipf':
                axis_min = 0
                axis_max = 2  * max(pred)
            else:
                axis_min = 0
                axis_max = 2 * max(obs)
            ax = fig.add_subplot(rows, columns, count+1)

            if i == 0 and j == 0:
                ax.set_title("Broken-stick")
            elif i == 0 and j == 1:
                ax.set_title("Lognormal")
            elif i == 0 and j == 2:
                ax.set_title("Log-series")
            elif i == 0 and j == 3:
                ax.set_title("Zipf")

            if j == 0:
                if dataset == '95':
                    ax.set_ylabel("MG-RAST 95%", rotation=90, size=12)
                elif dataset == '97':
                    ax.set_ylabel("MG-RAST 97%", rotation=90, size=12)
                elif dataset == '99':
                    ax.set_ylabel("MG-RAST 99%", rotation=90, size=12)

            macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=1,
                            plot_obj=plt.subplot(rows,columns,count+1))

            plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
            if method == 'zipf':
                plt.xlim(0, axis_max)
                plt.ylim(0, axis_max)
            else:
                plt.xlim(0, axis_max)
                plt.ylim(0, axis_max)
            r2s = ((INh2["R2"]))
            r2s = r2s.astype(float)
            mean_r2s = np.mean(r2s)
            std_r2s = np.std(r2s)
            print method, dataset
            print "Mean r2 " + str(mean_r2s)
            print "Standard dev. " + str(std_r2s)
            if method == 'zipf':
                getR2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/SequenceSimilarity/'+ method  + '_mle_' +dataset +'_NSR2_stratify.txt')
            elif method == 'lognorm':
                getR2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/SequenceSimilarity/'+ method  + '_' + lognormType + '_' +dataset  +'_NSR2_stratify.txt')
            else:
                getR2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/SequenceSimilarity/'+ method + '_' +dataset +'_NSR2_stratify.txt')
            r2_mean = np.mean(((getR2["R2"])))
            if method == 'geom':
                r2text = r"${}^{{2}}_{{m}} = {:.{p}f} $".format('r',r2_mean , p=3)
            else:
                r2text = r"${}^{{2}}_{{m}} = {:.{p}f} $".format('r',r2_mean , p=2)
            if method == 'geom':
                plt.text(0.28, 0.90, r2text,  fontsize=10,
                    horizontalalignment='center',
                    verticalalignment='center',transform = ax.transAxes)
            else:
                plt.text(0.25, 0.90, r2text,  fontsize=10,
                    horizontalalignment='center',
                    verticalalignment='center',transform = ax.transAxes)
            plt.tick_params(axis='both', which='major', labelsize=8)
            plt.subplots_adjust(wspace=0.0000000001, hspace=0.5)

            axins = inset_axes(ax, width="30%", height="30%", loc=4)

            hist_r2 = np.histogram(r2s, range=(0, 1))
            xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
            xvals = xvals[0:len(xvals)-1]
            yvals = hist_r2[0]
            plt.plot(xvals, yvals, 'k-', linewidth=2)
            plt.axis([0, 1, 0, 1.1 * max(yvals)])

            ax.set(adjustable='box-forced', aspect='equal')
            plt.setp(axins, xticks=[], yticks=[])

            count += 1

    plt.tight_layout(pad=1.5, w_pad=0.8, h_pad=0.8)
    fig.subplots_adjust(left=0.1)
    fig.text(0.50, 0.02, 'Predicted abundance', ha='center', va='center', fontsize=14)
    fig.text(0.03, 0.5, 'Observed abundance', ha='center', va='center', rotation='vertical', fontsize=14)
    fig_name = str(mydir + 'figures/' + figname + '_RGB.' + saveAs)
    plt.savefig(fig_name, dpi=600, format = saveAs)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()

def figS2(n=35289, figname = 'FigS2', data_dir=mydir, \
    stratify = True, radius=2, remove = 1, zipfType = 'mle', RGF = False, \
    saveAs = 'eps', lognormType = 'pln'):
    # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    # Used for Figure 3 Locey and White (2013)
    """Multiple obs-predicted plotter"""
    fig = plt.figure()
    count = 0
    plot_dim = 2
    fig.subplots_adjust(bottom= 0.30)
    methods = ['geom', 'lognorm', 'mete', 'zipf']
    for i, method in enumerate(methods):
        if method == 'zipf':
            obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Remove_1s/Stratified/'+ method + '_'+  zipfType+'_obs_pred_1_stratify.txt')
            INh2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_1s/Stratified/' + method + '_mle' + '_NSR2_1_stratify.txt')
        elif method == 'lognorm':
            obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Remove_1s/Stratified/'+ method + '_'+  lognormType+'_obs_pred_1_stratify.txt')
            INh2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_1s/Stratified/' + method + '_'+  lognormType + '_NSR2_1_stratify.txt')
        else:
            obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Remove_1s/Stratified/'+ method +'_obs_pred_1_stratify.txt')
            INh2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_1s/Stratified/' + method + '_NSR2_1_stratify.txt')
        obs = np.asarray(list(((obs_pred_data["obs"]))))
        pred = np.asarray(list(((obs_pred_data["pred"]))))
        site = np.asarray(list(((obs_pred_data["site"]))))

        obs2 = []
        pred2 = []
        site2 = []

        obs_all = np.asarray(obs)
        pred_all = np.asarray(pred)
        site_all = np.asarray(site)

        if n == 'all' or len(obs) <= n:
            obs2 = list(obs)
            pred2 = list(pred)
            site2 = list(site)

        else:
            if len(obs) > n:
                inds = np.random.choice(range(len(site)), size=n, replace=False)
                for ind in inds:
                    obs2.append(obs[ind])
                    pred2.append(pred[ind])
                    site2.append(site[ind])

        obs = np.asarray(obs2)
        pred = np.asarray(pred2)
        site =  np.asarray(site2)

        if method == 'zipf':
            axis_min = 0
            axis_max = 2  * max(pred)
        else:
            axis_min = 0
            axis_max = 2 * max(obs)
        ax = fig.add_subplot(plot_dim, plot_dim, count+1)

        if method == 'zipf':
            NSR2_BS = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/Remove_1s/'+ method  + '_mle_NSR2_1_stratify.txt')
        elif method == 'lognorm':
            NSR2_BS = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/Remove_1s/'+ method  + '_pln_NSR2_1_stratify.txt')
        else:
            NSR2_BS = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/Remove_1s/'+ method  +'_NSR2_1_stratify.txt')

        if method == 'geom':
            ax.set_title("Broken-stick")
        elif method == 'lognorm':
            ax.set_title("Lognormal")
        elif method == 'mete':
            ax.set_title("Log-series")
        elif method == 'zipf':
            ax.set_title("Zipf")
        macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=1,
                        plot_obj=plt.subplot(plot_dim,plot_dim,count+1))

        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
        if method == 'zipf':
            plt.xlim(0, axis_max)
            plt.ylim(0, axis_max)
        else:
            plt.xlim(0, axis_max)
            plt.ylim(0, axis_max)
        r2s = ((INh2["R2"]))
        r2s = r2s.astype(float)
        # insert r2 of all data
        r2_all = np.mean(((NSR2_BS["R2"])))
        print method
        r2text = r"${}^{{2}}_{{m}} = {:.{p}f} $".format('r',r2_all , p=2)
        if method == 'geom':
            plt.text(0.25, 0.90, r2text,  fontsize=14,
                horizontalalignment='center',
                verticalalignment='center',transform = ax.transAxes)
        else:
            plt.text(0.22, 0.90, r2text,  fontsize=14,
                horizontalalignment='center',
                verticalalignment='center',transform = ax.transAxes)
        plt.tick_params(axis='both', which='major', labelsize=10)
        plt.subplots_adjust(wspace=0.0000000001, hspace=0.5)

        axins = inset_axes(ax, width="30%", height="30%", loc=4)

        hist_r2 = np.histogram(r2s, range=(0, 1))
        xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
        xvals = xvals[0:len(xvals)-1]
        yvals = hist_r2[0]
        plt.plot(xvals, yvals, 'k-', linewidth=2)
        plt.axis([0, 1, 0, 1.1 * max(yvals)])

        ax.set(adjustable='box-forced', aspect='equal')
        plt.setp(axins, xticks=[], yticks=[])
        count += 1
    plt.tight_layout(pad=1.5, w_pad=0.8, h_pad=0.8)
    fig.text(0.50, 0.02, 'Predicted abundance', ha='center', va='center', fontsize=16)
    fig.text(0.08, 0.5, 'Observed abundance', ha='center', va='center', rotation='vertical', fontsize=16)
    fig_name = str(mydir + 'figures/' + figname + '_RGB.' + saveAs)
    plt.savefig(fig_name, dpi=600, format = saveAs)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()


def statOutput( data_dir=mydir, lognormType = 'pln', remove =0, seqSim = False):
    methods = ['geom', 'lognorm', 'mete', 'zipf']
    if seqSim != False:
        print seqSim
    for method in methods:
        if remove == 0 and seqSim == False:
            if method == 'zipf':
                nsr2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/' + method + '_mle' + '_NSR2_stratify.txt')
            elif method == 'lognorm':
                nsr2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/' + method + '_'+  lognormType + '_NSR2_stratify.txt')
            else:
                nsr2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/' + method + '_NSR2_stratify.txt')
        elif remove != 0 and seqSim == False:
            if method == 'zipf':
                nsr2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/Remove_1s/' + method + '_mle' + '_NSR2_1_stratify.txt')
            elif method == 'lognorm':
                nsr2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/Remove_1s/' + method + '_'+  lognormType + '_NSR2_1_stratify.txt')
            else:
                nsr2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/Remove_1s/' + method + '_NSR2_1_stratify.txt')
        else:
            if method == 'zipf':
                nsr2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/SequenceSimilarity/'+ method  + '_mle_' +seqSim +'_NSR2_stratify.txt')
            elif method == 'lognorm':
                nsr2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/SequenceSimilarity/'+ method  + '_' + lognormType + '_' +seqSim  +'_NSR2_stratify.txt')
            else:
                nsr2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified_Test/SequenceSimilarity/'+ method + '_' +seqSim +'_NSR2_stratify.txt')

        print method
        print "r2 mean = " + str(np.mean(((nsr2["R2"]))))
        print "r2 std = " +  str(np.mean(((nsr2["R2std"]))))
        if seqSim == False:
            domSlope = np.mean(((nsr2["NmaxPredSlope"])))
            evenSlope = np.mean(((nsr2["evennessPredSlope"])))
            skewSlope = np.mean(((nsr2["skewnessPredSlope"])))
            print "Dominance slope " + str(domSlope)
            diffDom = domSlope  - 1
            p_diffDom = ((np.absolute(diffDom) /  ( 0.5 *(domSlope + 1) ))) * 100
            diffDom_norm = diffDom / 1
            P_errorDom = np.absolute(diffDom_norm) * 100
            print "percent difference " + str(p_diffDom)
            print "percent error " + str(P_errorDom)

            print "Evenness slope " + str( evenSlope )
            diffEven = evenSlope  - (-0.31)
            p_diffEven = ((np.absolute(diffEven)  /  ( 0.5 *(np.absolute(evenSlope + -0.31)) ))) * 100
            diffEven_norm = diffEven / -0.31
            P_errorEven = np.absolute(diffEven_norm) * 100
            print "percent difference " + str(p_diffEven)
            print "percent error " + str(P_errorEven)

            print "Skewness slope " + str(skewSlope )
            diffSkew = skewSlope  - 0.13
            p_diffSkew = ((np.absolute(diffSkew)  /  ( 0.5 *(skewSlope + 0.13) ))) * 100
            diffSkew_norm = diffSkew / 0.13
            P_errorSkew = np.absolute(diffSkew_norm) * 100
            print "percent difference " + str(p_diffSkew)
            print "percent error " + str(P_errorSkew)


def figS3(data_dir=mydir, saveAs = 'eps'):
    models = ['geom', 'lognorm', 'mete', 'zipf']
    #models = ['geom']
    fig = plt.figure()
    count = 0
    for count_model, model in enumerate(models):
        IN = pd.read_csv(data_dir + 'data/Subsample_N/' + model +'_SubSampled_Data.txt', sep = ' ',
        names = ['Percent', 'N', 'S', 'Nmax_obs', 'Nmax_pred', 'r2'], index_col =None)

        P_0_5 = IN.loc[IN['Percent'] == 0.5][['r2']].values.tolist()
        P_0_25 = IN.loc[IN['Percent'] == 0.250000][['r2']].values.tolist()
        P_0_125 = IN.loc[IN['Percent'] == 0.125000][['r2']].values.tolist()
        P_0_0625 = IN.loc[IN['Percent'] == 0.062500][['r2']].values.tolist()
        P_0_03125 = IN.loc[IN['Percent'] == 0.031250][['r2']].values.tolist()
        P_0_015625 = IN.loc[IN['Percent'] == 0.015625][['r2']].values.tolist()
        N_0_5 = list(itertools.repeat(50, len(P_0_5)))
        N_0_25 = list(itertools.repeat(25, len(P_0_25)))
        N_0_125 = list(itertools.repeat(12.5, len(P_0_125)))
        N_0_0625 = list(itertools.repeat(6.25, len(P_0_0625)))
        N_0_03125 = list(itertools.repeat(3.125, len(P_0_03125)))
        N_0_015625 = list(itertools.repeat(1.5625, len(P_0_015625)))

        x = np.asarray(list(itertools.chain(N_0_015625, N_0_03125, N_0_0625, N_0_125, N_0_25, N_0_5)))
        y = np.asarray(list(itertools.chain(P_0_015625, P_0_03125, P_0_0625, P_0_125, P_0_25, P_0_5)))
        x = np.asarray([float(i) for i in x])
        x = np.log10(x)
        y = np.asarray([float(i) for i in y])
        ax = fig.add_subplot(2, 2, count + 1)
        #plt.scatter(x, y, alpha=0.05)
        plt.hexbin(x, y, mincnt=1, gridsize = 20, bins='log', cmap=plt.cm.jet)


        slope, intercept, r_value, p_value, std_err = \
            stats.linregress(x,y)

        print "slope is " + str(slope)
        print "r-value is " + str(r_value)
        print "p-value is " + str(p_value)
        plt.xlim(np.amin(x), np.amax(x))
        plt.ylim(np.amin(y), np.amax(y))
        predict_y = intercept + slope * x
        print model
        pred_error = y - predict_y
        plt.plot(x, predict_y, 'k-')

        degrees_of_freedom = len(x) - 2
        residual_std_error = np.sqrt(np.sum(pred_error**2) / \
            degrees_of_freedom)
        plt.xticks(fontsize = 10) # work on current fig
        plt.yticks(fontsize = 10)

        if model == 'geom':
            ax.set_title("Broken-stick")
        elif model == 'lognorm':
            ax.set_title("Lognormal")
        elif model == 'mete':
            ax.set_title("Log-series")
        elif model == 'zipf':
            ax.set_title("Zipf")

        count += 1
    name = str(data_dir + '/figures' + '/FigS3_RGB.' + saveAs)
    plt.tight_layout(pad=1.5, w_pad=0.8, h_pad=0.8)

    fig.text(0.50, 0.001, 'Percent of ' + r'$N,\; log_{10}$', ha='center', va='center', fontsize=16)
    fig.text(0.02, 0.5, r'$r_m^2$', ha='center', va='center', rotation='vertical', fontsize=16)
    plt.savefig(name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600, format = saveAs)

def fig1_Presentation(figname = 'Fig1', data_dir= mydir, saveAs = 'png'):
    SAD = [10000, 8000, 6000, 5000, 1000, 200, 100,  20, 18, 16, 14, 12, 10, 4,5,
        4, 4, 3, 3, 2, 2, 2, 2, 2,2, 1, 1, 1, 1, 1,1,1,1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    SAD.sort()
    SAD.reverse()
    x = range(1, len(SAD) +1)
    N = sum(SAD)
    S = len(SAD)

    geom = np.log10(mo.get_GeomSeries(N, S, False))

    logSeries = np.log10(mete.get_mete_rad(S, N)[0])


    lognorm_pred = mo.lognorm(SAD, 'pln')
    lognorm_SAD = np.log10(lognorm_pred.get_rad_from_obs())
    zipf_class = mo.zipf(SAD, 'fmin')
    pred_tuple = zipf_class.from_cdf()
    zipf_SAD = np.log10(pred_tuple[0])
    gamma = pred_tuple[1]

    SAD = np.log10(SAD)
    fig = plt.figure()
    plt.plot()

    max_y = max(max(SAD),  max(zipf_SAD))

    plt.plot(x, SAD,color = '#A9A9A9', linestyle = '-', linewidth=4, label="Observed")
    plt.plot(x, geom,color = 'blue', linestyle = '-', linewidth=4, label="Broken-stick")
    plt.plot(x, lognorm_SAD, color = 'black',linestyle = '-', linewidth=4, label="Lognormal")
    plt.plot(x, logSeries, color = '#228B22',linestyle = '-', linewidth=4, label="Log-series")
    plt.plot(x, zipf_SAD, color = 'red',linestyle = '-',linewidth=4,  label="Zipf")

    plt.tight_layout()
    plt.xlabel('Rank Abundance', fontsize = 22)
    plt.ylabel('Abundance, ' +r'$log_{10}$', fontsize = 22)
    output = "dorm_fix_prob.png"
    plt.legend(loc='upper right')

    #plt.yscale('log')
    #plt.yscale('log')
    plt.xlim(1, len(SAD))
    plt.ylim(-0.25 , max_y)

    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.legend(frameon=False, fontsize= 18)

    fig_name = str(mydir + 'figures/' + figname + '_RGB_Presentation.' + saveAs)
    plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600, \
        format = saveAs)
    plt.close()



#352899
#fig4()
#fig2(352899, figname = 'Fig2', data_dir=mydir, \
#    stratify = True, radius=2, remove = 0, zipfType = 'mle', RGF = False, lognormType = 'pln')
#fig3()
#fig1()
#figS2(352899)
#test(seqSim = '95')
#test(seqSim = '97')
#test(seqSim = '99')
#test(remove = 0)
#figS3()
fig1_Presentation()
