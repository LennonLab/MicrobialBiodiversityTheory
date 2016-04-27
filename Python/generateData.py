from __future__ import division
import imp, os, signal, datetime, random
import importData
import models
import numpy as np

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
    OUT = open(mydir + 'data/ObsPred/Stratified/lognorm_75_25_obs_pred_stratify.txt', 'w+')
    siteNSR2 = np.asarray(list(((IN_NSR2["site"]))))
    N = np.asarray(list(((IN_NSR2["N"]))))
    S = np.asarray(list(((IN_NSR2["S"]))))
    siteObsPred = np.asarray(list(((IN_Obs_Pred["site"]))))
    obs = np.asarray(list(((IN_Obs_Pred["obs"]))))

    uniqueSites = np.unique(siteNSR2)
    #randomSites = np.random.choice(uniqueSites, size=testNumber, replace=False)
    obs7525 = []
    pred7525 = []
    sites7525 = []

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
                signal.alarm(10)
                a = datetime.datetime.now()
                siteNSR2_i = siteNSR2[i]
                try:
                    # Whatever your function that might hang
                    SAD =  simLogNormFile.simLogNorm(N_i, S_i, sample_size).SimLogNormInt()
                    print 'countdown: ' + str(count)
                    b = datetime.datetime.now()
                    c = b - a
                    print str(c.seconds) + " seconds"
                except models.TimeoutException:
                    continue # continue the for loop if function takes more than x seconds
                else:
                    # Reset the alarm
                    signal.alarm(0)
                    count -= 1
                    for j in SAD:
                        pred7525.append(j)
                        sites7525.append(siteNSR2_i)
    for r, s in enumerate(siteObsPred):
        if s in sites7525:
            obs7525.append(obs[r])
            #sites7525.append(s)

    for x, site_x in enumerate(sites7525):
        print>> OUT, int(site_x), int(obs7525[x]), int(pred7525[x])
    OUT.close()


#getLogNormSim(testNumber = 100, sample_size = 1000)
