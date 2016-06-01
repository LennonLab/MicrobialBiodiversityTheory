from __future__ import division
import os, sys, signal, collections, argparse, optparse, math, datetime, imp
import importData as importData
import numpy as np
from operator import itemgetter
import models as mo
import metrics as me
import macroecotools
import mete as mete

mydir = os.path.expanduser("~/github/MicroMETE/")
importS = imp.load_source('predictS', mydir + 'lognormal/predictS.py')
simLogNormFile = imp.load_source('sim_lognormal', mydir + 'lognormal/sim_lognormal.py')

'''It works!'''
#S = importS.predictS(1000, 900, predictNmax=True).getS()

def generate_obs_pred_data(datasets, methods, size = 0, remove_obs = 0, zipfType = 'mle', lognormType = 'pln'):
    remove_obs = int(remove_obs)
    for method in methods:
        for dataset in datasets:

            signal.signal(signal.SIGALRM, mo.timeout_handler)

            if (method != 'zipf' and dataset != 'MGRAST'):

                if dataset == 'EMPclosed' or dataset == 'EMPopen' :
                    IN = mydir +"data/" + dataset + '-Data' + '/' + dataset +'-SADs.txt'
                    if method  == 'lognorm':
                        if remove_obs == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method + '_' + lognormType +'_'+dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method + '_' + lognormType +'_'+dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method + '_' + lognormType +'_'+dataset +'_obs_pred_' \
                                + str(remove_obs) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method + '_' + lognormType +'_'+dataset+'_NSR2_' \
                                + str(remove_obs) + '.txt','w+')
                    else:
                        if remove_obs == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method +'_'+dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method +'_'+dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method +'_'+dataset +'_obs_pred_' \
                                + str(remove_obs) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method +'_'+dataset+'_NSR2_' \
                                + str(remove_obs) + '.txt','w+')
                elif dataset == "HMP":
                    IN = mydir  + "data/" + dataset + '-Data' + '/' + dataset +'-SADs_NAP.txt'
                    if method == 'lognorm':
                        if remove_obs == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method + '_' + lognormType+'_'+dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method + '_' + lognormType+'_'+dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method + '_' + lognormType +'_'+dataset+'_obs_pred_' \
                                + str(remove_obs) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method + '_' + lognormType +'_'+dataset+'_NSR2_' \
                                + str(remove_obs) + '.txt','w+')
                    else:
                        if remove_obs == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method +'_'+dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method +'_'+dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method +'_'+dataset+'_obs_pred_' \
                                + str(remove_obs) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method +'_'+dataset+'_NSR2_' \
                                + str(remove_obs) + '.txt','w+')
                else:
                    IN = mydir + 'data/MGRAST-Data/' + dataset +  '/' + 'MGRAST-' + dataset + '-SADs.txt'
                    if method == 'lognorm':
                        if remove_obs == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method + '_' + lognormType +'_'+ 'MGRAST' + dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method + '_' + lognormType +'_'+ 'MGRAST' + dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method+ '_' + lognormType +'_'+ 'MGRAST' + dataset+'_obs_pred_' \
                                + str(remove_obs) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method + '_' + lognormType+'_'+ 'MGRAST' + dataset+'_NSR2_'  \
                                + str(remove_obs) + '.txt','w+')
                    else:
                        if remove_obs == 0:
                            OUT1 = open(mydir + "data/ObsPred/" + method +'_'+ 'MGRAST' + dataset+'_obs_pred.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/" + method +'_'+ 'MGRAST' + dataset+'_NSR2.txt','w+')
                        else:
                            OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method +'_'+ 'MGRAST' + dataset+'_obs_pred_' \
                                + str(remove_obs) + '.txt','w+')
                            OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method +'_'+ 'MGRAST' + dataset+'_NSR2_'  \
                                + str(remove_obs) + '.txt','w+')
            elif (method == 'zipf' and dataset != 'MGRAST'):

                if dataset == 'EMPclosed' or dataset == 'EMPopen':
                    IN = mydir + "data/" + dataset + '-Data' + '/' + dataset +'-SADs.txt'
                    if remove_obs == 0:
                        OUT1 = open(mydir + "data/ObsPred/" + method + '_' + zipfType + '_'+dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/" + method + '_' + zipfType +'_'+dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_obs_pred_' \
                            + str(remove_obs) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method + '_' + zipfType +'_'+ dataset+'_NSR2_'  \
                            + str(remove_obs) + '.txt','w+')
                elif dataset == 'HMP':
                    IN = mydir + "data/" + dataset + '-Data' + '/' + dataset +'-SADs_NAP.txt'
                    if remove_obs == 0:
                        OUT1 = open(mydir + "data/ObsPred/" + method + '_' + zipfType + '_'+dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/" + method + '_' + zipfType +'_'+dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_obs_pred_' \
                            + str(remove_obs) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method + '_' + zipfType +'_'+ dataset+'_NSR2_'  \
                            + str(remove_obs) + '.txt','w+')
                else:
                    IN = mydir  + "data/"+ 'MGRAST-Data/' + dataset +  '/' + 'MGRAST-' + dataset + '-SADs.txt'
                    if remove_obs == 0:
                        OUT1 = open(mydir + "data/ObsPred/" + method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/" + method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_obs_pred_' \
                            + str(remove_obs) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_NSR2_' \
                            + str(remove_obs) + '.txt','w+')

            elif dataset == 'MGRAST':
                IN = mydir + "data/"+ 'MGRAST-Data/MGRAST/MGRAST-SADs.txt'
                if remove_obs == 0:
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
                        OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_obs_pred_' \
                            + str(remove_obs) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_NSR2_' \
                            + str(remove_obs) + '.txt','w+')
                    elif method == 'lognorm':
                        OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method + '_' + lognormType +'_' + dataset+'_obs_pred_' \
                            + str(remove_obs) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method + '_' + lognormType +'_' + dataset+'_NSR2_' \
                            + str(remove_obs) + '.txt','w+')
                    else:
                        OUT1 = open(mydir + "data/ObsPred/remove_" + str(remove_obs)  + 's/'+ method +'_' + dataset+'_obs_pred_' \
                            + str(remove_obs) + '.txt','w+')
                        OUT2 = open(mydir + "data/NSR2/remove_" + str(remove_obs)  + 's/'+ method +'_' + dataset+'_NSR2_' \
                            + str(remove_obs) + '.txt','w+')

            num_lines = sum(1 for line in open(IN))


            random_sites = np.random.randint(num_lines, size=size)

            line_count = 0
            print OUT1, OUT2
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
                if remove_obs != 0:
                    obs = [e for e in obs if int(e) > remove_obs]
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
                        signal.alarm(6)
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
                print method
                if method == 'zipf':

                    if dataset == 'EMPclosed' or dataset == 'EMPopen' or dataset == 'HMP':

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
                                print NmaxObs, NmaxPred, r2
                                print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                                    evennessPred, skewnessObs, skewnessPred, gamma, r2
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

def stratifyData(totalSADs = 500, zipfType = 'mle', \
    lognormType = 'pln', remove = True, data_dir= mydir, remove_obs = 0):
    datasets = ['EMPclosed','HMP', 'MGRAST']
    methods = ['geom', 'mete', 'zipf', 'lognorm']
    # Number of lines in each file
    MGRAST_sites = 1174
    HMP_sites = 4504
    EMPclosed_sites = 14979
    Total = MGRAST_sites + HMP_sites + EMPclosed_sites
    for i, method in enumerate(methods):
        if method == 'zipf':
            if remove_obs == 0:
                OUT1 = open(data_dir + 'data/ObsPred/Stratified/'+ method + '_' + zipfType +'_obs_pred_stratify.txt', 'wr')
                OUT2 = open(data_dir + 'data/NSR2/Stratified/'+ method  + '_' + zipfType +'_NSR2_stratify.txt', 'wr')
            else:
                OUT1 = open(data_dir + 'data/ObsPred/Remove_' + str(remove_obs) + \
                    's/Stratified/'+ method + '_' + zipfType +'_obs_pred_' + str(remove_obs) + '_stratify.txt', 'wr')
                OUT2 = open(data_dir + 'data/NSR2/Remove_' + str(remove_obs) + \
                    's/Stratified/'+ method  + '_' + zipfType +'_NSR2_stratify_' +  str(remove_obs) + '_stratify.txt', 'wr')
        elif method == 'lognorm':
            if remove_obs == 0:
                OUT1 = open(data_dir + 'data/ObsPred/Stratified/'+ method + '_' + lognormType +'_obs_pred_stratify.txt', 'wr')
                OUT2 = open(data_dir + 'data/NSR2/Stratified/'+ method  + '_' + lognormType +'_NSR2_stratify.txt', 'wr')
            else:
                OUT1 = open(data_dir + 'data/ObsPred/Remove_' + str(remove_obs) + \
                    's/Stratified/'+ method + '_' + lognormType +'_obs_pred_' + str(remove_obs) + '_stratify.txt', 'wr')
                OUT2 = open(data_dir + 'data/NSR2/Remove_' + str(remove_obs) + \
                    's/Stratified/'+ method  + '_' + lognormType +'_NSR2_stratify_' +  str(remove_obs) + '_stratify.txt', 'wr')
        else:
            if remove_obs == 0:
                OUT1 = open(data_dir + 'data/ObsPred/Stratified/'+ method +'_obs_pred_stratify.txt', 'wr')
                OUT2 = open(data_dir + 'NSR2/Stratified/'+ method  +'_NSR2_stratify.txt', 'wr')
            else:
                OUT1 = open(data_dir + 'data/ObsPred/Remove_' + str(remove_obs) + \
                    's/Stratified/'+ method + '_obs_pred_' + str(remove_obs) + '_stratify.txt', 'wr')
                OUT2 = open(data_dir + 'data/NSR2/Remove_' + str(remove_obs) + \
                    's/Stratified/'+ method  + '_NSR2_stratify_' +  str(remove_obs) + '_stratify.txt', 'wr')
        count2 = 0
        count1 = 0
        for j, dataset in enumerate(datasets):
            lineCount = 0
            removeSADs = []
            print datasets
            if remove_obs == 0:
                get_bad_zipfs = importData.import_NSR2_data(data_dir + 'data/NSR2/' + 'zipf' + '_'+ 'mle' +'_'+dataset +'_NSR2.txt')
            else:
                get_bad_zipfs = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_' +str(remove_obs) + 's/zipf_mle_'+dataset  +'_NSR2_' +  str(remove_obs) +'.txt')
            site = np.asarray(list(((get_bad_zipfs["site"]))))
            N = np.asarray(list(((get_bad_zipfs["N"]))))
            S = np.asarray(list(((get_bad_zipfs["S"]))))
            r2s = np.asarray(list(((get_bad_zipfs["R2"]))))
            zipNsite = zip(site, N, S, r2s)
            for x in zipNsite:
                if float(x[3]) < 0.2:
                    removeSADs.append(int(x[0]))
            removeSADs = np.asarray(removeSADs)
            if dataset == 'MGRAST':
                n = 239 - len(removeSADs)
            else:
                n = totalSADs
            print "sites to remove selected"
            if remove_obs == 0:
                if method == 'zipf':
                    if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                        ( str(dataset) == 'EMPopen'):
                        obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/' + method + '_'+ zipfType +'_'+dataset+'_obs_pred.txt')
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method + '_'+ zipfType +'_'+ dataset+'_NSR2.txt')
                    elif str(dataset) == 'MGRAST':
                        obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/' + method + '_'+ zipfType + '_' + 'MGRAST_obs_pred.txt')
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+ '_'+ zipfType  + '_MGRAST_NSR2.txt')
                    else:
                        obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/' + method + '_'+ zipfType + '_' + 'MGRAST' + dataset +'_obs_pred.txt')
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+ '_'+ zipfType  + '_MGRAST'+dataset+'_NSR2.txt')
                elif method == 'lognorm':
                    if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                        ( str(dataset) == 'EMPopen'):
                        obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/' + method + '_'+ lognormType +'_'+dataset+'_obs_pred.txt')
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method + '_'+ lognormType +'_'+ dataset+'_NSR2.txt')
                    elif str(dataset) == 'MGRAST':
                        obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/' + method + '_'+ lognormType + '_' + 'MGRAST_obs_pred.txt')
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+ '_'+ lognormType  + '_MGRAST_NSR2.txt')
                    else:
                        obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/' + method + '_'+ lognormType + '_' + 'MGRAST' + dataset +'_obs_pred.txt')
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+ '_'+ lognormType  + '_MGRAST'+dataset+'_NSR2.txt')
                else:
                    if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                        ( str(dataset) == 'EMPopen'):
                        obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/' + method+'_'+dataset+'_obs_pred.txt')
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+'_'+dataset+'_NSR2.txt')
                    elif str(dataset) == 'MGRAST':
                        obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/' + method + '_' + 'MGRAST_obs_pred.txt')
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+'_MGRAST_NSR2.txt')
                    else:
                        obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/' + method + '_' + 'MGRAST' + dataset +'_obs_pred.txt')
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+'_MGRAST'+dataset+'_NSR2.txt')
            else:
                if method == 'zipf':
                    obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Remove_' \
                        + str(remove_obs) + 's/' + method + '_'+ zipfType +'_'+dataset+'_obs_pred_' +  str(remove_obs) +'.txt')
                    nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_' \
                        + str(remove_obs) + 's/' + method + '_'+ zipfType +'_'+dataset+'_NSR2_' +  str(remove_obs) +'.txt')
                elif method == 'lognorm':
                    obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Remove_' \
                        + str(remove_obs) + 's/' + method + '_'+ lognormType +'_'+dataset+'_obs_pred_' +  str(remove_obs) +'.txt')
                    nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_' \
                        + str(remove_obs) + 's/' + method + '_'+ lognormType +'_'+dataset+'_NSR2_' +  str(remove_obs) +'.txt')
                else:
                    obs_pred_data = importData.import_obs_pred_data(data_dir + 'data/ObsPred/Remove_' \
                        + str(remove_obs) + 's/' + method +'_'+dataset+'_obs_pred_' +  str(remove_obs) +'.txt')
                    nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_' \
                        + str(remove_obs) + 's/' + method +'_'+dataset+'_NSR2_' +  str(remove_obs) +'.txt')
            siteNSR2 = np.asarray( map(float, list(((nsr2_data["site"])))))
            N = np.asarray(list(((nsr2_data["N"]))))
            S = np.asarray(list(((nsr2_data["S"]))))
            NmaxObs = np.asarray(list(((nsr2_data["NmaxObs"]))))
            NmaxPred = np.asarray(list(((nsr2_data["NmaxPred"]))))
            evennessObs = np.asarray(list(((nsr2_data["evennessObs"]))))
            evennessPred = np.asarray(list(((nsr2_data["evennessPred"]))))
            skewnessObs = np.asarray(list(((nsr2_data["skewnessObs"]))))
            skewnessPred = np.asarray(list(((nsr2_data["skewnessPred"]))))
            R2 = np.asarray(list(((nsr2_data["R2"]))))

            site = ((obs_pred_data["site"]))
            obs = list(((obs_pred_data["obs"])))
            pred = list(((obs_pred_data["pred"])))
            pred
            obs2 = []
            pred2 = []
            site2 = []
            if remove == True:
                siteNSR2_cleaned = np.setdiff1d(siteNSR2, removeSADs)
                uniqueSites = np.unique(siteNSR2_cleaned)
            else:
                uniqueSites = np.unique(siteNSR2)

            randomSites = np.random.choice(uniqueSites, size=n, replace=False)

            #for enumSite, randomSite in enumerate(randomSites):

            print len(np.unique(site)),  len(siteNSR2)
            for enumSite, randomSite in enumerate(randomSites):
                for p, q in enumerate(siteNSR2):
                    if q == randomSite:
                        print>> OUT2, count1, N[p], S[p], NmaxObs[p], NmaxPred[p], \
                            evennessObs[p], evennessPred[p], skewnessObs[p], skewnessPred[p], R2[p]
                for r, s in enumerate(site):
                    if s == randomSite:
                        obs2.append(obs[r])
                        pred2.append(pred[r])
                        site2.append(s)
                        print>> OUT1, count1, obs[r], pred[r]
                count1 += 1
                    #if (r == 0):
                    #    print>> OUT1, count_sites, obs[r], pred[r]
                    #elif r != 0 and obs[r] > obs[r-1]:
                    #    count_sites += 1
                    #    print>> OUT1, count_sites, obs[r], pred[r]
                    #else:
                    #    print>> OUT1, count_sites, obs[r], pred[r]

            print method, dataset

            obs = np.asarray(obs2)
            pred = np.asarray(pred2)
            zippedSiteObsPred = zip(site2,obs2,pred2)
            #site = np.asarray(site2)
            k_minus1 = obs2[0]

        OUT1.close()
        OUT2.close()


#stratifyData(remove_obs = 1)
#datasets = ['EMPclosed','HMP', 'MGRAST']
datasets = ['EMPclosed', 'HMP']
methods = ['zipf']
generate_obs_pred_data(datasets, methods, remove_obs = 1)
