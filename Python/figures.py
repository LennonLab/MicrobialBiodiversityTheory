
import imp, os, signal, datetime, random
import importData
import macroecotools
import mete
import scipy as sp
from scipy import stats, optimize

import numpy as np
import  matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import models as mo
#import generate_figs_zipf as gfz


mydir = os.path.expanduser("~/github/MicroMETE/")
importPredictS = imp.load_source('sim_lognormal', mydir + 'lognormal/predictS.py')

def Supp(figname = 'Supp', data_dir=mydir, radius=2):
    # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    # Used for Figure 3 Locey and White (2013)
    """Multiple obs-predicted plotter"""
    fig = plt.figure()
    count = 0

    plot_dim = 2
    IN_Obs_Pred = importData.import_obs_pred_data(mydir + \
        'data/ObsPred/Stratified/lognorm_75_25_obs_pred_stratify_test.txt')
    site = np.asarray(list(((IN_Obs_Pred["site"]))))
    obs = np.asarray(list(((IN_Obs_Pred["obs"]))))
    pred7525 = np.asarray(list(((IN_Obs_Pred["pred7525"]))))
    predPln = np.asarray(list(((IN_Obs_Pred["predPln"]))))
    toIterate = [pred7525, predPln]
    for x in range(2):
        axis_min = 0
        axis_max = 2 * max(obs)
        #print plot_dim
        ax = fig.add_subplot(plot_dim, plot_dim, count+1)
        if x == 0:
            ax.set_title(r"$\mathbf{75:25\, Simulation}$")
        else:
            ax.set_title(r"$\mathbf{Lognormal\, MLE}$")

        macroecotools.plot_color_by_pt_dens(toIterate[x], obs, radius, loglog=1,
                        plot_obj=plt.subplot(plot_dim,plot_dim,count+1))
        #
        #plt.text(0.1, 0.9,'matplotlib', ha='center', va='center', transform=ax.transAxes)


        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')

        plt.xlim(0, axis_max)
        plt.ylim(0, axis_max)
        #r2s = ((INh2["R2"]))
        #r2s = r2s.astype(float)
        # insert r2 of all data
        r2_all = macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(toIterate[x]))
        r2text = r"${}^{{2}}_{{m}} = {:.{p}f} $".format('r',r2_all , p=2)

        plt.text(0.18, 0.93, r2text,  fontsize=10,
            horizontalalignment='center',
            verticalalignment='center',transform = ax.transAxes)
        plt.tick_params(axis='both', which='major', labelsize=7)
        plt.subplots_adjust(wspace=0.5, hspace=0.3)

        axins = inset_axes(ax, width="30%", height="30%", loc=4)

        #hist_r2 = np.histogram(r2s, range=(0, 1))
        #xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
        #xvals = xvals[0:len(xvals)-1]
        #yvals = hist_r2[0]
        #plt.plot(xvals, yvals, 'k-', linewidth=2)
        #plt.axis([0, 1, 0, 1.1 * max(yvals)])
        ax.set(adjustable='box-forced', aspect='equal')
        #plt.setp(axins, xticks=[], yticks=[])

        count += 1
    fig.text(0.50, 0.04, r'$Predicted \; rank-abundance$', ha='center', va='center')
    fig.text(0.05, 0.5, r'$Observed \; rank-abundance$', ha='center', va='center', rotation='vertical')
    fig_name = str(mydir + 'figures/' + figname + '.png')
    plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()

def figSuppp(figname = 'SuppFig3', data_dir=mydir, radius=2):
    fig = plt.figure()
    plot_dim = 2
    count = 0

    IN_Obs_Pred = importData.import_NSR2_data(mydir + \
        'data/NSR2/Stratified/lognorm_pln_NSR2_stratify.txt')
    N = np.asarray(list(((IN_Obs_Pred["N"]))))
    S = np.asarray(list(((IN_Obs_Pred["S"]))))
    NmaxObs = np.asarray(list(((IN_Obs_Pred["NmaxObs"]))))
    NmaxPred = []
    SPred = []
    for i in range(len(N)):
        NmaxPred_i = importPredictS.predictS(N[i], NmaxObs[i], predictNmax=True).getNmax()
        SPred_i = importPredictS.predictS(N[i], NmaxObs[i], predictNmax=True).getS()
        NmaxPred.append(NmaxPred_i)
        SPred.append(SPred_i)
    NmaxPred = np.asarray(NmaxPred)
    SPred = np.asarray(SPred)
    toIteratePred = [NmaxPred, SPred]
    toIterateObs = [NmaxObs, S]
    for x in range(2):
        axis_min = 0
        axis_max = 2 * max(toIteratePred[x])
        #print plot_dim
        ax = fig.add_subplot(plot_dim-1, plot_dim, count+1)
        if x == 0:
            ax.set_title(r"$\mathbf{N_{max}}$")
        else:
            ax.set_title(r"$\mathbf{S}$")

        macroecotools.plot_color_by_pt_dens(toIteratePred[x], toIterateObs[x], radius, loglog=1,
                        plot_obj=plt.subplot(plot_dim-1,plot_dim,count+1))
        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
        plt.xlim(axis_min, axis_max)
        plt.ylim(0, axis_max)
        r2_all = macroecotools.obs_pred_rsquare(np.log10(toIterateObs[x]), np.log10(toIteratePred[x]))
        r2text = r"${}^{{2}}_{{m}} = {:.{p}f} $".format('r',r2_all , p=2)
        plt.text(0.18, 0.93, r2text,  fontsize=10,
            horizontalalignment='center',
            verticalalignment='center',transform = ax.transAxes)
        plt.tick_params(axis='both', which='major', labelsize=7)
        plt.subplots_adjust(wspace=0.5, hspace=0.3)

        #axins = inset_axes(ax, width="30%", height="30%", loc=4)

        ax.set(adjustable='box-forced', aspect='equal')
        #plt.setp(axins, xticks=[], yticks=[])

        count += 1
    fig.text(0.50, 0.04, r'$Predicted$', ha='center', va='center')
    fig.text(0.05, 0.5, r'$Observed$', ha='center', va='center', rotation='vertical')

    fig_name = str(mydir + 'figures/' + figname + '.png')
    plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()

def fig1(figname = 'Fig1', data_dir= mydir):
    SAD = [10000, 8000, 6000, 5000, 1000, 200, 100,  20, 18, 16, 14, 12, 10, 4,5,
        4, 4, 3, 3, 2, 2, 2, 2, 2,2, 1, 1, 1, 1, 1,1,1,1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    SAD.sort()
    SAD.reverse()
    x = range(1, len(SAD) +1)
    N = sum(SAD)
    S = len(SAD)
    geom = mo.get_GeomSeries(N, S, False)

    logSeries = mete.get_mete_rad(S, N)[0]


    lognorm_pred = mo.lognorm(SAD, 'pln')
    lognorm_SAD = lognorm_pred.get_rad_from_obs()
    zipf_class = mo.zipf(SAD, 'fmin')
    pred_tuple = zipf_class.from_cdf()
    zipf_SAD = pred_tuple[0]
    gamma = pred_tuple[1]

    fig = plt.figure()
    plt.plot()

    max_y = max(max(SAD),  max(zipf_SAD))

    plt.plot(x, SAD,color = '#A9A9A9', linestyle = '-', linewidth=2, label="Observed")
    plt.plot(x, geom,color = '#00008B', linestyle = '-', linewidth=2, label="Broken-stick")
    plt.plot(x, lognorm_SAD, color = '#0000CD',linestyle = '--', linewidth=2, label="Lognormal")
    plt.plot(x, logSeries, color = '#FF4500',linestyle = '-.', linewidth=2, label="Log-series")
    plt.plot(x, zipf_SAD, color = '#8B0000',linestyle = ':',linewidth=2,  label="Zipf")

    plt.tight_layout()
    #plt.xlabel(r'$Rank \; Abundance$', fontsize = 18)
    plt.xlabel('Rank Abundance', fontsize = 22)
    plt.ylabel('Abundance, ' +r'$log_{10}$', fontsize = 22)
    #plt.ylabel(r'$Abundance, \, log_{10}$', fontsize = 18)
    output = "dorm_fix_prob.png"
    plt.legend(loc='upper right')
    plt.yscale('log')
    plt.xlim(1, len(SAD))
    plt.ylim(1, max_y)
    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.legend(frameon=False, fontsize= 18)

    #fig.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig_name = str(mydir + 'figures/' + figname + '.png')
    plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    #plt.xscale()
    plt.close()

def fig2(n, figname = 'Fig2', data_dir=mydir, \
    stratify = True, radius=2, remove = 0, zipfType = 'mle', RGF = False, lognormType = 'pln'):
    # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    # Used for Figure 3 Locey and White (2013)
    """Multiple obs-predicted plotter"""
    fig = plt.figure()
    count = 0
    plot_dim = 2
    methods = ['geom', 'lognorm', 'mete', 'zipf']
    for i, method in enumerate(methods):
        if method == 'zipf':
            obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/'+ method + '_'+  zipfType+'_obs_pred_stratify.txt')
            INh2 = importData.import_NSR2_data(data_dir + 'data/NSR2/Stratified/' + method + '_mle' + '_NSR2_stratify.txt')
        #if method == 'rgf':
        #    obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/Stratified/'+ 'zipf' + '_'+  method+'_obs_pred_stratify.txt')
        #    INh2 = import_NSR2_data(data_dir + 'NSR2/Stratified/' +'zipf' + '_'+  method + '_NSR2_stratify.txt')
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
        #print plot_dim
        ax = fig.add_subplot(plot_dim, plot_dim, count+1)

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
        r2_all = macroecotools.obs_pred_rsquare(np.log10(obs_all), np.log10(pred_all))
        #r2text = r"${}^{{2}}_{{m}} = {:.{p}f} $".format('r',r2_all , p=2)
        r2text = r"${:.{p}f} $".format(r2_all , p=2)
        if method == 'geom':
            plt.text(0.18, 0.90, r2text,  fontsize=14,
                horizontalalignment='center',
                verticalalignment='center',transform = ax.transAxes)
        else:
            plt.text(0.15, 0.90, r2text,  fontsize=14,
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
    #fig.tight_layout()

    #plt.tight_layout(pad=0.4, w_pad=0.8, h_pad=0.5)
    plt.tight_layout(pad=1.5, w_pad=0.8, h_pad=0.8)
    #plt.subplots_adjust(wspace=0.2, hspace=0.1)
    fig.text(0.50, 0.02, 'Predicted rank-abundance', ha='center', va='center', fontsize=14)
    fig.text(0.08, 0.5, 'Observed rank-abundance', ha='center', va='center', rotation='vertical', fontsize=14)
    fig_name = str(mydir + 'figures/' + figname + '.png')
    plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()


def fig3(figname = 'Fig3', \
    zipfType = 'mle', lognormType = 'pln', Stratified = True, data_dir= mydir):
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
            #if i == 0:
            x = np.log10(np.asarray(list(((obs_pred_data[param])))))
            print "nmax" + str(np.mean(np.asarray(list(((obs_pred_data["NmaxObs"]))))))

            print len(x),  len(y)
            mean_x = np.mean(x)
            mean_y = np.mean(y)
            std_error = sp.stats.sem(y)
            print method, param
            print "mean modified r2 = " + str(mean_y)
            print "modified r2 standard error = " + str(std_error)
            print "mean " + param  + " is " + str(np.mean(np.asarray(list(((obs_pred_data[param]))))))
            ax = fig.add_subplot(2, 2, count+1)
            print len(x), len(y)
            macroecotools.plot_color_by_pt_dens(x, y, 0.1, loglog=0,
                            plot_obj=plt.subplot(2, 2, count+1))
            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            print "r-value is " + str(r_value)
            print "p-value is " + str(p_value)

            plt.xlim(np.amin(x), np.amax(x))
            if method == 'lognorm' or method == 'zipf':
                plt.ylim(0,2)
            else:
                plt.ylim(-1,1)

            predict_y = intercept + slope * x
            pred_error = y - predict_y
            degrees_of_freedom = len(x) - 2
            residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
            plt.plot(x, predict_y, 'k-')
            plt.axhline(linewidth=2, color='darkgrey',ls='--')
            plt.tick_params(axis='both', which='major', labelsize=10)

            #plt.subplots_adjust(wspace=0.5, hspace=0.5 )
            if i == 0 and j == 0:
                ax.set_title("Broken-stick", fontsize = 15)
                ax.set_ylabel(r'$r^{2}_{m}$', fontsize = 19)
            if i == 0 and j == 1:
                ax.set_title("Lognormal", fontsize = 15)
            if i == 0 and j == 2:
                ax.set_title("Log-series", fontsize = 15)
                ax.set_xlabel('Abundance, ' +r'$log_{10}$', fontsize = 17)
                ax.set_ylabel(r'$r^{2}_{m}$', fontsize = 19)
            if i == 0 and j == 3:
                ax.set_title("Zipf", fontsize = 15)
                ax.set_xlabel('Abundance, ' +r'$log_{10}$', fontsize = 17)
            ax.tick_params(axis='x', labelsize=12)
            ax.tick_params(axis='y', labelsize=12)

            count += 1

    #fig.subplots_adjust(wspace = 0.00001, hspace = 0.1, top=0.70)
    plt.tight_layout(pad=0.8, w_pad=0.8, h_pad=0.8)
    fig_name = str(mydir + 'figures/' + figname  + '.png')
    plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    #plt.xscale()
    plt.close()



def fig4(figname = 'Fig4', data_dir=mydir, radius=2):
    fig = plt.figure()
    plot_dim = 1
    count = 0

    IN_Obs_Pred = importData.import_NSR2_data(mydir + \
        'data/NSR2/Stratified/lognorm_pln_NSR2_stratify.txt')
    N = np.asarray(list(((IN_Obs_Pred["N"]))))
    S = np.asarray(list(((IN_Obs_Pred["S"]))))
    NmaxObs = np.asarray(list(((IN_Obs_Pred["NmaxObs"]))))
    # order
    models = ['geom', 'lognorm', 'mete', 'zipf']
    modelSlopes = [0.647520323289, 0.942904468437, 0.769214774397, 0.954497727096]
    modelInterepts = [0.116508916992, 0.292527611072, 0.19240314275, 0.189954627996]
    for g, model in enumerate(models):
        NmaxPred = []
        SPred = []
        for i in range(len(N)):
            NmaxPred_i = importPredictS.predictS(N[i], NmaxObs[i], \
                predictNmax=True).getNmax(b = modelInterepts[g], slope = modelSlopes[g])
            SPred_i = importPredictS.predictS(N[i], NmaxObs[i], predictNmax=True).getS()
            NmaxPred.append(NmaxPred_i)
            SPred.append(SPred_i)
        NmaxPred = np.asarray(NmaxPred)
        SPred = np.asarray(SPred)
        axis_min = 0
        axis_max = 2 * max(NmaxObs)
        ax = fig.add_subplot(2, 2, count+1)
        #ax.set_title(r"$\mathbf{N_{max}}$", y=1.03)
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
        print max(NmaxPred)
        r2_all = macroecotools.obs_pred_rsquare(np.log10(NmaxObs), np.log10(NmaxPred))
        r2text = r"${:.{p}f} $".format(r2_all , p=2)
        plt.text(0.18, 0.91, r2text,  fontsize=13,
            horizontalalignment='center',
            verticalalignment='center',transform = ax.transAxes)
        plt.tick_params(axis='both', which='major', labelsize=12)
        plt.subplots_adjust(wspace=0.00001, hspace=0.3)

        #axins = inset_axes(ax, width="30%", height="30%", loc=4)

        ax.set(adjustable='box-forced', aspect='equal')
        #plt.setp(axins, xticks=[], yticks=[])

        count += 1
    fig.text(0.50, 0.03, 'Predicted, ' +r'$log_{10}(N_{max})$', ha='center', va='center', fontsize = 16)
    fig.text(0.09, 0.5, 'Observed, ' +r'$log_{10}(N_{max})$', ha='center', va='center', rotation='vertical',\
        fontsize = 16)
    fig_name = str(mydir + 'figures/' + figname + '.png')
    plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()

def figS1(n, figname = 'FigS1', data_dir=mydir, radius=2, zipfType = 'mle', lognormType = 'pln'):
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

            #if method == 'geom':
            #    ax.set_title("Broken-stick")
            #elif method == 'lognorm':
            #    ax.set_title("Lognormal")
            #elif method == 'mete':
            #    ax.set_title("Log-series")
            #elif method == 'zipf':
            #    ax.set_title("Zipf")
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
            # insert r2 of all data
            r2_all = macroecotools.obs_pred_rsquare(np.log10(obs_all), np.log10(pred_all))
            r2text = r"${:.{p}f} $".format(r2_all , p=2)
            if method == 'geom':
                plt.text(0.18, 0.90, r2text,  fontsize=12,
                    horizontalalignment='center',
                    verticalalignment='center',transform = ax.transAxes)
            else:
                plt.text(0.15, 0.90, r2text,  fontsize=12,
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

    #plt.tight_layout(pad=0.4, w_pad=0.8, h_pad=0.5)
    plt.tight_layout(pad=1.5, w_pad=0.8, h_pad=0.8)
    fig.subplots_adjust(left=0.1)
    #plt.subplots_adjust(wspace=0.2, hspace=0.1)
    fig.text(0.50, 0.02, 'Predicted rank-abundance', ha='center', va='center', fontsize=14)
    fig.text(0.03, 0.5, 'Observed rank-abundance', ha='center', va='center', rotation='vertical', fontsize=14)
    fig_name = str(mydir + 'figures/' + figname + '.png')
    plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()




#352899
#fig2(352899, figname = 'Fig2', data_dir=mydir, \
#    stratify = True, radius=2, remove = 0, zipfType = 'mle', RGF = False, lognormType = 'pln')
#fig3()
figS1(n=352899)
