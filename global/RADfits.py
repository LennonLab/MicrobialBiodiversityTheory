from __future__ import division
import sys

sys.path.append("GenAnalysis/tools/")
import macroeco_distributions as md
import macroecotools as mt
import feasible_functions as ff
import predRADs
import mete
import pln

import cloud
import numpy as np


def getPredRADs(N, S, Nmax):
    
    PRED = []
    
    # Predicted geometric series
    predRAD = predRADs.get_GeomSeries(N, S, False) # False mean no zeros allowed
    PRED.append(predRAD)
    
    # Predicted log-series
    logSeries = mete.get_mete_rad(S, N)
    predRAD = logSeries[0]
    PRED.append(predRAD)
    
    # Predicted PLN
    predRAD = pln.get_rad_from_obs(RAD, 'pln') 
    PRED.append(predRAD)
    
    sample_size = 10
    # Predicted from compositions (basically geometric series)
    predRAD = getPred(N, S, maxn, 'compositions', sample_size)
    PRED.append(predRAD)
    
    # Predicted from Fraction 1: Power Fraction
    predRAD = getPred(N, S, maxn, 'power fraction', sample_size)
    PRED.append(predRAD)
    
    # Predicted from Fraction 2: Random non-preemption
    predRAD = getPred(N, S, maxn, 'random fraction non-preemption', sample_size)
    PRED.append(predRAD)
    
    # Predicted from Fraction 3: Random preemption 
    predRAD = getPred(N, S, maxn, 'random fraction', sample_size)
    PRED.append(predRAD)
    
    return PRED





NScombos = [[10**10, 10**5], [10**20, 10**5], [10**30, 10**5], [10**40, 10**5],
            [10**10, 10**7], [10**20, 10**7], [10**30, 10**7], [10**40, 10**7],
            [10**10, 10**9], [10**20, 10**9], [10**30, 10**9], [10**40, 10**9]]
        

NScombos = [[1000, 50]]            
                                    
models = ['geometric series', 'log series', 'log normal', 'compositions',
          'random fraction preemption', 'random fraction']

# uncomment the following if wanting to overwrite exiting output file
OUT1 = open('/Users/lisalocey/Desktop/global/NS_SimData.txt', 'w+')
print>>OUT1, 'N', 'S', 'geometric series', 'log series', 'log normal',
print>>OUT1, 'compositions', 'power fraction',
print>>OUT1, 'random fraction preemption', 'random fraction'
OUT1.close()

                                          
for combo in NScombos:
    N, S = combo
    Nmax = int(1.00*N)
                                                   #! GETTING READY TO RUN SOME CODE !
    for model in models:
        print N, S, model 
        
        # using PiCloud
        #job_id = cloud.call(getPredRADs, N, S, Nmax, _type='m1') # use picloud
        #PRED = cloud.result(job_id)
        
        # or using the native env
        PRED = getPredRADs(N, S, Nmax)
        
        # Each element of PRED is a list, i.e. PRED[0] is [predicted RAD, r2 of obs vs. pred]
        GS, LS, PLN, FS, COMPS, PowerFrac, RandFracPre, RandFrac = PRED
        
        OUT1 = open('NS_SimData.txt', 'a')
        for i, ab in enumerate(RAD): # each species in each site at each date
            # gets its own record of predicted abundances
            print>>OUT1, N, S, int(GS_rad[i]), LS_rad[i], PLN_rad[i], FS_rad[i]
        
        OUT1.close()