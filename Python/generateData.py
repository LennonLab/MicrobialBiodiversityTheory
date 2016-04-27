from __future__ import division
import imp, os, signal, datetime, random
import importData
import models
import numpy as np
from operator import itemgetter

mydir = os.path.expanduser("~/github/MicroMETE/")
importS = imp.load_source('predictS', mydir + 'lognormal/predictS.py')
simLogNormFile = imp.load_source('sim_lognormal', mydir + 'lognormal/sim_lognormal.py')

'''It works!'''
#S = importS.predictS(1000, 900, predictNmax=True).getS()

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
    signal.signal(signal.SIGALRM, models.timeout_handler)

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


getLogNormSim(testNumber = 100, sample_size = 10000)
