
import imp, os, signal, datetime, random
import importData
import macroecotools
import numpy as np
import  matplotlib.pyplot as plt
from mpl_toolkits.axes_grid.inset_locator import inset_axes
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
            ax.set_title(r"$\mathbf{Broken-stick}$")
        elif model == 'lognorm':
            ax.set_title(r"$\mathbf{Lognormal}$")
        elif model == 'mete':
            ax.set_title(r"$\mathbf{Log-series}$")
        elif model == 'zipf':
            ax.set_title(r"$\mathbf{Zipf}$")
        macroecotools.plot_color_by_pt_dens(NmaxPred, NmaxObs, radius, loglog=1,
                        plot_obj=plt.subplot(2,2,count+1))
        plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
        plt.xlim(axis_min, axis_max)
        plt.ylim(0, axis_max)
        print max(NmaxPred)
        r2_all = macroecotools.obs_pred_rsquare(np.log10(NmaxObs), np.log10(NmaxPred))
        r2text = r"${}^{{2}}_{{m}} = {:.{p}f} $".format('r',r2_all , p=2)
        plt.text(0.22, 0.91, r2text,  fontsize=13,
            horizontalalignment='center',
            verticalalignment='center',transform = ax.transAxes)
        plt.tick_params(axis='both', which='major', labelsize=7)
        plt.subplots_adjust(wspace=0.5, hspace=0.3)

        #axins = inset_axes(ax, width="30%", height="30%", loc=4)

        ax.set(adjustable='box-forced', aspect='equal')
        #plt.setp(axins, xticks=[], yticks=[])

        count += 1
    fig.text(0.50, 0.04, r'$Predicted\, log_{10}(N_{max})$', ha='center', va='center', fontsize = 16)
    fig.text(0.04, 0.5, r'$Observed\,log_{10}(N_{max})$', ha='center', va='center', rotation='vertical',\
        fontsize = 16)
    fig_name = str(mydir + 'figures/' + figname + '.png')
    plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
    plt.close()

fig4()
