from __future__ import division
import os, sys, signal, collections, argparse, optparse, math, datetime, imp
import importData as importData
import numpy as np
from operator import itemgetter
import models as mo
import metrics as me
import macroecotools

mydir = os.path.expanduser("~/github/MicroMETE/")
importS = imp.load_source('predictS', mydir + 'lognormal/predictS.py')
simLogNormFile = imp.load_source('sim_lognormal', mydir + 'lognormal/sim_lognormal.py')

'''It works!'''
#S = importS.predictS(1000, 900, predictNmax=True).getS()

def generate_obs_pred_data(datasets, methods, size = 0, remove = 0, zipfType = 'mle', lognormType = 'pln'):
    remove = int(remove)
    #if remove != 0:
    #    newpath1 = mydir + "ObsPred/Remove_" + str(remove) + 's/'
    #    if not os.path.exists(newpath1):
    #        os.makedirs(newpath1)
    #    newpath2 = mydir + "NSR2/Remove_" + str(remove) + 's/'
    #    if not os.path.exists(newpath2):
    #        os.makedirs(newpath2)
    for method in methods:
        for dataset in datasets:

            signal.signal(signal.SIGALRM, mo.timeout_handler)

            if (method != 'zipf' and dataset != 'MGRAST'):

                if dataset == 'EMPclosed' or dataset == 'EMPopen' :
                    IN = mydir +"data/" + dataset + '-Data' + '/' + dataset +'-SADs.txt'
                    if method  == 'lognorm':
                        if remove == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method + '_' + lognormType +'_'+dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method + '_' + lognormType +'_'+dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method + '_' + lognormType +'_'+dataset +'_obs_pred_' \
                                + str(remove) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + lognormType +'_'+dataset+'_NSR2_' \
                                + str(remove) + '.txt','w+')
                    else:
                        if remove == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method +'_'+dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method +'_'+dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method +'_'+dataset +'_obs_pred_' \
                                + str(remove) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method +'_'+dataset+'_NSR2_' \
                                + str(remove) + '.txt','w+')
                elif dataset == "HMP":
                    IN = mydir  + "data/" + dataset + '-Data' + '/' + dataset +'-SADs_NAP.txt'
                    if method == 'lognorm':
                        if remove == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method + '_' + lognormType+'_'+dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method + '_' + lognormType+'_'+dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method + '_' + lognormType +'_'+dataset+'_obs_pred_' \
                                + str(remove) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + lognormType +'_'+dataset+'_NSR2_' \
                                + str(remove) + '.txt','w+')
                    else:
                        if remove == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method +'_'+dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method +'_'+dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method +'_'+dataset+'_obs_pred_' \
                                + str(remove) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method +'_'+dataset+'_NSR2_' \
                                + str(remove) + '.txt','w+')
                else:
                    IN = mydir + 'data/MGRAST-Data/' + dataset +  '/' + 'MGRAST-' + dataset + '-SADs.txt'
                    if method == 'lognorm':
                        if remove == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method + '_' + lognormType +'_'+ 'MGRAST' + dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method + '_' + lognormType +'_'+ 'MGRAST' + dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method+ '_' + lognormType +'_'+ 'MGRAST' + dataset+'_obs_pred_' \
                                + str(remove) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + lognormType+'_'+ 'MGRAST' + dataset+'_NSR2_'  \
                                + str(remove) + '.txt','w+')
                    else:
                        if remove == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method +'_'+ 'MGRAST' + dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method +'_'+ 'MGRAST' + dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method +'_'+ 'MGRAST' + dataset+'_obs_pred_' \
                                + str(remove) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method +'_'+ 'MGRAST' + dataset+'_NSR2_'  \
                                + str(remove) + '.txt','w+')
            elif (method == 'zipf' and dataset != 'MGRAST'):

                if dataset == 'EMPclosed' or dataset == 'EMPopen':
                    IN = mydir + "data/" + dataset + '-Data' + '/' + dataset +'-SADs.txt'
                    if remove == 0:
                        OUT1 = open(mydir + "data/ObsPred/" + method + '_' + zipfType + '_'+dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/" + method + '_' + zipfType +'_'+dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_'+ dataset+'_NSR2_'  \
                            + str(remove) + '.txt','w+')
                elif dataset == 'HMP':
                    IN = mydir + "data/" + dataset + '-Data' + '/' + dataset +'-SADs_NAP.txt'
                    if remove == 0:
                        OUT1 = open(mydir + "data/ObsPred/" + method + '_' + zipfType + '_'+dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/" + method + '_' + zipfType +'_'+dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_'+ dataset+'_NSR2_'  \
                            + str(remove) + '.txt','w+')
                else:
                    IN = mydir  + "data/"+ 'MGRAST-Data/' + dataset +  '/' + 'MGRAST-' + dataset + '-SADs.txt'
                    if remove == 0:
                        OUT1 = open(mydir + "data/ObsPred/" + method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/" + method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_NSR2_' \
                            + str(remove) + '.txt','w+')

            elif dataset == 'MGRAST':
                IN = mydir + "data/"+ 'MGRAST-Data/MGRAST/MGRAST-SADs.txt'
                if remove == 0:
                    if method == 'zipf':
                        OUT1 = open(mydir + "data/ObsPred/" + method + '_' + zipfType +'_'+ 'MGRAST_obs_pred.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/" + method + '_' + zipfType +'_'+ 'MGRAST_NSR2.txt','w+')
                    elif method == 'lognorm':
                        OUT1 = open(mydir + "data/ObsPred/" + method + '_' + lognormType +'_'+ 'MGRAST_obs_pred.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/" + method + '_' + lognormType +'_'+ 'MGRAST_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "data/ObsPred/" + method +'_'+ 'MGRAST_obs_pred.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/" + method +'_'+ 'MGRAST_NSR2.txt','w+')
                else:
                    if method == 'zipf':
                        OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_NSR2_' \
                            + str(remove) + '.txt','w+')
                    elif method == 'lognorm':
                        OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method + '_' + lognormType +'_' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + lognormType +'_' + dataset+'_NSR2_' \
                            + str(remove) + '.txt','w+')
                    else:
                        OUT1 = open(mydir + "data/ObsPred/Remove_" + str(remove)  + 's/'+ method +'_' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/Remove_" + str(remove)  + 's/'+ method +'_' + dataset+'_NSR2_' \
                            + str(remove) + '.txt','w+')

            num_lines = sum(1 for line in open(IN))


            random_sites = np.random.randint(num_lines, size=size)

            line_count = 0

            for j,line in enumerate(open(IN)):
                if dataset == "HMP":
                    #line = line.split()
                    line = line.strip().split(',')
                    line = [x.strip(' ') for x in line]
                    line = [x.strip('[]') for x in line]
                    site_name = line[0]
                    line.pop(0)
                elif size == 0:
                    line = eval(line)
                else:
                    line = eval(line)
                    if j not in random_sites:
                        continue
                obs = map(int, line)
                if remove != 0:
                    if int(remove) in obs:
                        obs = obs.remove(int(remove))
                if obs is None:
                    continue
                if len(obs) == 0:
                    continue
                N = sum(obs)
                S = len(obs)
                NmaxObs = np.amax(obs)
                evennessObs = me.e_simpson(obs)
                skewnessObs = me.skewness(obs)


                if S <= 10 or N <= S:
                    num_lines -= 1
                    continue

                num_lines -= 1

                obs.sort()
                obs.reverse()
                if dataset == "HMP":
                    print method, dataset, N, S, site_name, ' countdown: ', num_lines
                else:
                    print method, dataset, N, S, ' countdown: ', num_lines

                if method == 'geom': # Predicted geometric series
                    pred = mo.get_GeomSeries(N, S, False) # False mean no zeros allowed

                elif method == 'mete': # Predicted log-series
                    logSeries = mete.get_mete_rad(S, N)
                    pred = logSeries[0]
                elif method == 'lognorm' and lognormType == 'glm':
                    lognorm_pred = mo.lognorm(obs, lognormType)
                    pred = lognorm_pred.lognorm_glm()
                    pred = np.ceil(pred)
                    pred.astype(int)
                elif method == 'lognorm' and lognormType == 'pln':
                    signal.alarm(20)
                    a = datetime.datetime.now()
                    print "testing"
                    try:
                        # Whatever your function that might hang
                        # use S
                        lognorm_pred = mo.lognorm(obs, lognormType)
                        pred = lognorm_pred.get_rad_from_obs()
                        b = datetime.datetime.now()
                        c = b - a
                        print str(c.seconds) + " seconds"
                    except mo.TimeoutException:
                        continue # continue the for loop if function takes more than x seconds
                    else:
                        # Reset the alarm
                        signal.alarm(0)

                elif method == 'zipf' and zipfType == 'glm':
                    zipf_pred = mo.zipf(obs, 'fmin')
                    zipf_glm = zipf_pred.from_glm()
                    pred = np.ceil(zipf_glm)
                    pred.astype(int)
                    print line_count
                    #numpy.ceil
                elif method == 'zipf' and zipfType == 'rgf':
                    a = datetime.datetime.now()

                    #line = map(int, line)
                    # Start the timer. Once 1 second is over, a SIGALRM signal is sent.
                    signal.alarm(3)
                    # This try/except loop ensures that
                    #   you'll catch TimeoutException when it's sent.
                    try:
                        # Whatever your function that might hang
                        zipf_class_ = mo.zipf(obs, 'fmin')
                        pred = zipf_class_.zipf_rgf()
                        b = datetime.datetime.now()
                        c = b - a
                        print str(c.seconds) + " seconds"
                        print obs
                    except mo.TimeoutException:
                        continue # continue the for loop if function takes more than x seconds
                    else:
                        # Reset the alarm
                        signal.alarm(0)

                elif method == 'zipf' and zipfType == 'mle':

                        #line = map(int, line)
                        # Start the timer. Once 1 second is over, a SIGALRM signal is sent.
                        signal.alarm(4)
                        # This try/except loop ensures that
                        #   you'll catch TimeoutException when it's sent.
                        try:
                            # Whatever your function that might hang
                            # use S
                            #rv = stats.zipf(Zipf_solve_line)
                            zipf_class = mo.zipf(obs, 'fmin')
                            pred_tuple = zipf_class.from_cdf()
                            pred = pred_tuple[0]
                            gamma = pred_tuple[1]
                        except mo.TimeoutException:
                            continue # continue the for loop if function takes more than x seconds
                        else:
                            # Reset the alarm
                            signal.alarm(0)
                print pred
                NmaxPred = np.amax(pred)
                evennessPred = me.e_simpson(pred)
                skewnessPred = me.skewness(pred)
                r2 = macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))
                print " r2:", r2
                if r2 == -float('inf') or r2 == float('inf') or r2 == float('Nan'):
                    print r2 + " is Nan or inf, removing..."
                    continue
                if method == 'zipf':
                    if dataset == 'EMPclosed' or dataset == 'EMPopen' or dataset == 'HMP':
                        OUT2 = open(mydir +"data/" + "NSR2/" + method + '_' + zipfType +'_'+dataset+'_NSR2.txt','a+')
                        if zipfType == 'glm' or zipfType == 'rgf':
                            if dataset == 'HMP':
                                print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                                    evennessPred, skewnessObs, skewnessPred, r2, site_name
                            else:
                                print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                                    evennessPred, skewnessObs, skewnessPred, r2
                        else:
                            if dataset == 'HMP':
                                print>> OUT2, j, N, S, NmaxObs, NmaxPred,evennessObs, \
                                    evennessPred, skewnessObs, skewnessPred, gamma, r2, site_name
                            else:
                                print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                                    evennessPred, skewnessObs, skewnessPred, gamma, r2
                        OUT2.close()
                    else:
                        if zipfType == 'glm' or zipfType == 'rgf':
                            print>> OUT2, j, N, S,NmaxObs, NmaxPred, evennessObs, \
                                evennessPred, skewnessObs, skewnessPred, r2
                        else:
                            print>> OUT2, j, N, S,NmaxObs, NmaxPred, evennessObs, \
                                evennessPred, skewnessObs, skewnessPred, gamma, r2
                else:
                    if dataset == 'HMP':
                        print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                            evennessPred, skewnessObs, skewnessPred, r2, site_name
                    else:
                        print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                            evennessPred, skewnessObs, skewnessPred, r2

                for i, sp in enumerate(pred):
                    print>> OUT1, j, obs[i], pred[i]
                    line_count += 1
            print line_count
            OUT1.close()
            OUT2.close()
        print dataset



def getLogNormSim(testNumber = 100, sample_size = 1000):
    '''This function randomly samples a number of sites set by testNumber from
    'Stratified' dataset using the 75-25 simulation for the lognormal.
    Because some SADs take a very long time to generate (> 2 minutes), the function
    runs on a timer with a timeout function that moves to the next randomly chosed
    SAD (sampled without replacement), stopping when after testNumber of successful
    runs.
    '''
    IN_NSR2 = importData.import_NSR2_data(mydir + 'data/NSR2/Stratified/lognorm_pln_NSR2_stratify.txt')
    IN_Obs_Pred = importData.import_obs_pred_data(mydir + 'data/ObsPred/Stratified/lognorm_pln_obs_pred_stratify.txt')
    OUT = open(mydir + 'data/ObsPred/Stratified/lognorm_75_25_obs_pred_stratify_test.txt', 'w+')
    siteNSR2 = np.asarray(list(((IN_NSR2["site"]))))
    N = np.asarray(list(((IN_NSR2["N"]))))
    S = np.asarray(list(((IN_NSR2["S"]))))
    siteObsPred = np.asarray(list(((IN_Obs_Pred["site"]))))
    obs = np.asarray(list(((IN_Obs_Pred["obs"]))))
    pred = np.asarray(list(((IN_Obs_Pred["pred"]))))

    uniqueSites = np.unique(siteNSR2)
    #randomSites = np.random.choice(uniqueSites, size=testNumber, replace=False)
    obs7525 = []
    pred7525 = []
    sites7525 = []

    pred_pln = []
    sites_pln =[]
    signal.signal(signal.SIGALRM, mo.timeout_handler)

    count = testNumber
    while count > 0:
        #randomSite = np.random.choice(uniqueSites, size=1, replace=False)
        index = random.randint(0,len(uniqueSites)-1)
        randomSite = uniqueSites[index]
        uniqueSites = np.delete(uniqueSites,index)
        for i, site in enumerate(siteNSR2):
            if site == randomSite:
                N_i = N[i]
                S_i = S[i]
                a = datetime.datetime.now()
                #siteNSR2_i = siteNSR2[i]
                SAD =  simLogNormFile.simLogNorm(N_i, S_i, sample_size).SimLogNormInt()
                if len(SAD) != S_i:
                    continue
                print 'countdown: ' + str(count)
                print site, S_i, len(SAD)
                count -= 1
                for j in SAD:
                    pred7525.append(j)
                    sites7525.append(site)

    zipSitesPred7527 = zip(sites7525, pred7525)
    #print zipSitesPred7527
    indexes = np.unique(sites7525, return_index=True)[1]
    uniqueSites7525 = [sites7525[index] for index in sorted(indexes)]
    zipOsPredPln = zip(siteObsPred, obs, pred)
    zipOsPredPlnFilter = [x for x in zipOsPredPln if x[0] in uniqueSites7525]
    #print zipOsPredPlnFilter
    #print len(zipOsPredPlnFilter)
    #zipSitesPred7527Sort = sorted(L, key=itemgetter(0))

    countTest = 0
    for spot, uniqueSite7525 in enumerate(uniqueSites7525):
        for r, s in enumerate(siteObsPred):
            if int(s) == uniqueSite7525:
                print>> OUT, int(zipSitesPred7527[spot][0]),int(obs[r]), int(pred[r]), int(zipSitesPred7527[countTest][1])
                countTest += 1
                #obs7525.append(obs[r])
                #pred_pln.append(pred[r])
                #sites_pln.append(s)
    #print "pred sites obs752527 pred7525    "
    #print len(pred_pln), len(sites_pln), len(obs7525), len(pred7525)
    #for x, site_x in enumerate(sites7525):
    #    print>> OUT, int(site_x), int(sites_pln[x]),int(obs7525[x]), int(pred7525[x]), int(pred_pln[x])
    OUT.close()


datasets = ['EMPclosed']
methods = ['lognorm']

generate_obs_pred_data(datasets, methods, remove = 1)
