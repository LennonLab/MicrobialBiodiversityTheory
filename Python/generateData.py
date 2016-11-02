from __future__ import division
import os, sys, signal, collections, argparse, optparse, math, datetime, imp
import importData as importData
import numpy as np
import pandas as pd
from operator import itemgetter
import macroeco_distributions as md
import modelsAndMetrics as mo
import macroecotools
import mete as mete
from scipy import stats
import random
mydir = os.path.expanduser("~/github/MicroMETE/")


"""This code was written using MIT liscenced code from the following Weecology
repos: METE (https://github.com/weecology/METE) and macroecotools
(https://github.com/weecology/macroecotools). """


def HMP_OTU_to_sparese_SbyS(path=mydir):
    otu_count = pd.read_table(path + 'data/HMP-Data/hmp1.v35.hq.otu.counts.txt', sep='\t', index_col=False)
    sparsesbys = pd.melt(otu_count, id_vars=['collection'])
    sparsesbys = sparsesbys[sparsesbys['value'] > 0]
    sparsesbys.columns = ['Sample', 'OTU', 'Count']
    sparsesbys.to_csv(path + "data/HMP-Data/HMPsparseSbyS.txt", sep='\t', index=False)

def Match_NAP_to_sparse_SbyS(path=mydir, timeseries=True):
    IN = path + 'data/HMP-Data/ppAll_V35_map.txt'
    metadata = pd.read_table(IN, sep='\t', index_col=False)
    timeseries = bool(timeseries)
    metadata = metadata[np.isfinite(metadata['NAP'])]
    metadata[['NAP']] = metadata[['NAP']].astype(str)
    metadata.drop_duplicates(cols='NAP', take_last=True)
    if timeseries == True:
        pass
    else:
        metadata[['VisitNo']] = metadata[['VisitNo']].astype(float)
        metadata = metadata.loc[metadata['VisitNo'] == float(1)]
        metadata.to_csv(path + "data/HMP-Data/ppAll_V35_map_noTimeseries.txt", sep='\t', index=False)
    print metadata.shape
    # Pull the RSID number & convert to numpy array
    NAPs = metadata[['NAP']].values
    NAPs = NAPs.flatten()
    IN_OTU = path + 'data/HMP-Data/HMPsparseSbyS.txt'
    if timeseries == True:
        OUT = open(path + 'data/HMP-Data/HMPsparseSbyS_NAP.txt','w+')
    else:
        OUT = open(path + 'data/HMP-Data/HMPsparseSbyS_NAP_noTimeseries.txt','w+')
    for j, line in enumerate(open(IN_OTU)):
        if '.PPS' in line:
            continue
        else:
            line_split = line.split()
            OTU = line_split[1]
            count = line_split[2]
            NAP = str(line_split[0].split('.')[0])
            try:
                float(NAP)
                NAP_test = str(NAP) + '.0'
                if NAP_test  in NAPs:
                    print>> OUT, line_split[0], OTU, count
            except ValueError:
                print "Not a float"

def get_SADs(path, name, closedref=True):
    SADdict = {}
    DATA = path + name + '-data.txt'

    with open(DATA) as f:

        for d in f:
            if d.strip():
                d = d.split()
                if name == 'GENTRY':
                    site = d[0]
                    #species = d[1] # Dataset name plus species identifier
                    abundance = float(d[-1])

                else:
                    site = d[0]
                    if closedref == True:
                        for i in d:
                            if 'unclassified' in i:
                                continue
                            elif 'unidentified' in i:
                                continue

                    abundance = float(d[-1])

                if abundance > 0:
                    if site in SADdict:
                        SADdict[site].append(abundance)
                    else:
                        SADdict[site] = [abundance]

    SADs = SADdict.values()
    filteredSADs = []
    for sad in SADs:
        if len(sad) >= 10:
            filteredSADs.append(sad)
    return filteredSADs


def get_SADs_HMP(path= mydir, timeseries=True):
    timeseries = bool(timeseries)
    if timeseries == True:
        IN = path + 'data/HMP-Data/HMPsparseSbyS_NAP.txt'
        OUT =  open(path+'HMP-Data/' + 'HMP-SADs_NAP.txt', 'w+')
    else:
        IN = path + 'data/HMP-Data/HMPsparseSbyS_NAP_noTimeseries.txt'
        OUT =  open(path+'data/HMP-Data/' + 'HMP-SADs_NAP_noTimeseries.txt', 'w+')
    SADdict = {}
    with open(IN) as f:
        for d in f:
            if d.strip():
                d = d.split()
                site = d[0]
                abundance = int(d[-1])
                if abundance > 0:
                    if site not in SADdict:
                        #SADdict[site] = [abundance]
                        SADdict[site] = [abundance]
                    else:
                        SADdict[site].append(abundance)

    filtered_SADdict = {}

    for key, value in SADdict.iteritems():
        if len(value) >= 10:
            filtered_SADdict[key] = value

    print "You have " + str(len(SADdict)) + " sites"

    # first value of the line is the site name, rest is SAD
    # site name filtered out in generate obs pred data
    for key, value in filtered_SADdict.iteritems():
        output = value.insert(0,key)
        print>> OUT, value


def get_SADs_mgrast(path = mydir, thresholds=1):

    datasets = ['BOVINE', 'CATLIN', 'CHU', 'HYDRO', 'LAUB']

    for t in thresholds:
        SADdict = {}

        for d in datasets:
            name = d+t
            print name

            filepath  = path + 'data/MGRAST-Data/'+t+'/'+name+'/'+name+'-data.txt'
            with open(filepath) as f:
                for d in f:
                    if d.strip():
                        d = d.split()
                        site = d[0]
                        abundance = int(d[-1])

                        if abundance > 0:
                            if site in SADdict:
                                SADdict[site].append(abundance)
                            else:
                                SADdict[site] = [abundance]
        SADs = SADdict.values()
        filteredSADs = []

        for sad in SADs:
            if len(sad) >= 10:
                filteredSADs.append(sad)

        print t, len(filteredSADs)

        OUT =  open(path+'MGRAST-Data/'+t+'/MGRAST-'+t+'-SADs.txt', 'w')

        for sad in SADs:
            print>> OUT, sad

def merge_SADs_mgrast(path = mydir):
    SADdict = {}
    fungi_list = map(str,np.arange(4484945.3, 4485075.3,1))
    filepath  = path + 'data/MGRAST-Data/MGRAST/MGRAST-data.txt'
    with open(filepath, 'r') as f:
        for d in f:
            if d.strip():
                d = d.split()
                site = d[0]
                abundance = d[-1]
                if ('-' in str(site))  or (str(site) in fungi_list):
                     continue
                abundance = int(abundance)
                if abundance > 0:
                    if site in SADdict:
                        SADdict[site].append(abundance)
                    else:
                        SADdict[site] = [abundance]
        SADs = SADdict.values()
        filteredSADs = []

        for sad in SADs:
            if len(sad) >= 10:
                filteredSADs.append(sad)

        print  len(filteredSADs)

        OUT =  open(path+'MGRAST-Data/MGRAST/MGRAST-SADs.txt', 'w')

        for sad in SADs:
            print>> OUT, sad


def EMP_SADs(name, mgrast, path= mydir):
    minS = 10
    IN = path + 'data' + name + '-SSADdata.txt'
    n = sum(1 for line in open(IN))

    SiteDict = {}

    with open(IN) as f:
        for d in f:
            n -= 1
            if d.strip():

                d = d.split()
                #species = d[0]
                sample = d[1]
                abundance = float(d[2])

                if abundance > 0:

                    if sample not in SiteDict:
                        SiteDict[sample] = [abundance]

                    else:
                        SiteDict[sample].append(abundance)
    SADs = SiteDict.values()
    filteredSADs = []
    for sad in SADs:
        if len(sad) >= minS:
            filteredSADs.append(sad)
    return filteredSADs



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
            for j,line in enumerate(open(IN)):
                if dataset == "HMP":
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
                evennessObs = mo.e_simpson(obs)
                skewnessObs = mo.skewness(obs)

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
                    pred = mo.get_Geom(N, S, False) # False mean no zeros allowed

                elif method == 'mete': # Predicted log-series
                    logSeries = mete.get_mete_rad(S, N)
                    pred = logSeries[0]
                    p = logSeries[1]
                    beta = math.log(p) * -1
                    ll = sum(np.log(mete.trunc_logser_pmf(range(1, int(N) + 1), p, N)))
                elif method == 'lognorm' and lognormType == 'glm':
                    lognorm_pred = mo.lognorm(obs, lognormType)
                    pred = lognorm_pred.lognorm_glm()
                    pred = np.ceil(pred)
                    pred.astype(int)
                elif method == 'lognorm' and lognormType == 'pln':
                    signal.alarm(21)
                    a = datetime.datetime.now()
                    try:
                        # Whatever your function that might hang
                        # use S
                        lognorm_pred = mo.lognorm(obs, lognormType)
                        get_rad = lognorm_pred.get_rad_from_obs()
                        pred = get_rad[0]
                        mu = get_rad[1]
                        sigma = get_rad[2]
                        ll = md.pln_ll(obs, mu, sigma)
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
                elif method == 'zipf' and zipfType == 'rgf':
                    a = datetime.datetime.now()

                    # Start the timer. Once 1 second is over, a SIGALRM signal is sent.
                    signal.alarm(5)
                    # This try/except loop ensures that
                    #   you'll catch TimeoutException when it's sent.
                    try:
                        # Whatever your function that might hang
                        zipf_class_ = mo.zipf(obs, 'fmin')
                        zipf_return = zipf_class_.zipf_rgf()
                        pred = zipf_return[0]
                        gamma = zipf_return[1]
                        b = datetime.datetime.now()
                        c = b - a
                        print str(c.seconds) + " seconds"
                    except mo.TimeoutException:
                        continue # continue the for loop if function takes more than x seconds
                    else:
                        # Reset the alarm
                        signal.alarm(0)

                elif method == 'zipf' and zipfType == 'mle':

                        # Start the timer. Once 1 second is over, a SIGALRM signal is sent.
                        signal.alarm(7)
                        # This try/except loop ensures that
                        #   you'll catch TimeoutException when it's sent.
                        try:
                            # Whatever your function that might hang
                            # use S
                            zipf_class = mo.zipf(obs, 'fmin')
                            pred_tuple = zipf_class.from_cdf()
                            pred = pred_tuple[0]
                            gamma = pred_tuple[1]
                            ll = md.zipf_ll(obs, gamma)
                        except mo.TimeoutException:
                            continue # continue the for loop if function takes more than x seconds
                        else:
                            # Reset the alarm
                            signal.alarm(0)
                NmaxPred = np.amax(pred)
                evennessPred = mo.e_simpson(pred)
                skewnessPred = mo.skewness(pred)
                r2 = macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))



                print " r2:", r2
                if r2 == -float('inf') or r2 == float('inf') or r2 == float('Nan'):
                    print r2 + " is Nan or inf, removing..."
                    continue
                if dataset == 'EMPclosed' and r2 < 0.2:
                    continue
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
                                    evennessPred, skewnessObs, skewnessPred, gamma, \
                                    ll, r2, site_name
                            else:
                                print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                                    evennessPred, skewnessObs, skewnessPred, gamma, ll, r2
                    else:
                        if zipfType == 'glm' or zipfType == 'rgf':
                            print>> OUT2, j, N, S,NmaxObs, NmaxPred, evennessObs, \
                                evennessPred, skewnessObs, skewnessPred, r2
                        else:
                            print>> OUT2, j, N, S,NmaxObs, NmaxPred, evennessObs, \
                                evennessPred, skewnessObs, skewnessPred, gamma, ll, r2
                elif method == 'mete':
                    if dataset == 'HMP':
                        print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                            evennessPred, skewnessObs, skewnessPred, p, ll, r2, site_name
                    else:
                        print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                            evennessPred, skewnessObs, skewnessPred, p, ll, r2
                elif method == 'lognorm':
                    if dataset == 'HMP':
                        print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                            evennessPred, skewnessObs, skewnessPred, mu, sigma, \
                            ll, r2, site_name
                    else:
                        print>> OUT2, j, N, S, NmaxObs, NmaxPred, evennessObs, \
                            evennessPred, skewnessObs, skewnessPred, mu, sigma, \
                            ll, r2
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


def stratifyDataOnce(datasets, zipfType = 'mle', \
    lognormType = 'pln', remove = True, data_dir= mydir, remove_obs = 0):
    if remove_obs == 0:
        totalSADs = 239
    else:
        totalSADs = 108
    methods = ['geom', 'mete', 'zipf', 'lognorm']
    #methods = ['lognorm','zipf', 'mete']

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
                    's/Stratified/'+ method  + '_' + zipfType +'_NSR2_' +  str(remove_obs) + '_stratify.txt', 'wr')
        elif method == 'lognorm':
            if remove_obs == 0:
                OUT1 = open(data_dir + 'data/ObsPred/Stratified/'+ method + '_' + lognormType +'_obs_pred_stratify.txt', 'wr')
                OUT2 = open(data_dir + 'data/NSR2/Stratified/'+ method  + '_' + lognormType +'_NSR2_stratify.txt', 'wr')
            else:
                OUT1 = open(data_dir + 'data/ObsPred/Remove_' + str(remove_obs) + \
                    's/Stratified/'+ method + '_' + lognormType +'_obs_pred_' + str(remove_obs) + '_stratify.txt', 'wr')
                OUT2 = open(data_dir + 'data/NSR2/Remove_' + str(remove_obs) + \
                    's/Stratified/'+ method  + '_' + lognormType +'_NSR2_' +  str(remove_obs) + '_stratify.txt', 'wr')
        else:
            if remove_obs == 0:
                OUT1 = open(data_dir + 'data/ObsPred/Stratified/'+ method +'_obs_pred_stratify.txt', 'wr')
                OUT2 = open(data_dir + 'data/NSR2/Stratified/'+ method  +'_NSR2_stratify.txt', 'wr')
            else:
                OUT1 = open(data_dir + 'data/ObsPred/Remove_' + str(remove_obs) + \
                    's/Stratified/'+ method + '_obs_pred_' + str(remove_obs) + '_stratify.txt', 'wr')
                OUT2 = open(data_dir + 'data/NSR2/Remove_' + str(remove_obs) + \
                    's/Stratified/'+ method  + '_NSR2_' +  str(remove_obs) + '_stratify.txt', 'wr')
        count2 = 0
        count1 = 0
        for j, dataset in enumerate(datasets):
            lineCount = 0
            removeSADs = []
            print method, dataset
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
                if remove_obs != 0:
                    continue
                else:
                    if float(x[3]) < 0.2:
                        removeSADs.append(int(x[0]))
            removeSADs = np.asarray(removeSADs)
            if dataset == 'MGRAST' and remove_obs == 0:
                n = totalSADs - len(removeSADs)
            elif dataset == 'MGRAST' and remove_obs == 1:
                n = totalSADs - len(removeSADs)
            elif dataset == 'EMPclosed' and method == 'zipf' and remove_obs == 1:
                n = 40
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
            #slope_Nmax_iter, intercept_Nmax_iter, r_value_Nmax_iter, p_value_Nmax_iter, std_err_Nmax_iter = \
            #    stats.linregress(np.log10(N),np.log10(NmaxPred))
            if method != 'geom':
                ll = np.asarray(list(((nsr2_data["ll"]))))
            if method == 'zipf':
                gamma = np.asarray(list(((nsr2_data["gamma"]))))
            if method == 'mete':
                p_mete = np.asarray(list(((nsr2_data["p"]))))
            if method == 'lognorm':
                mu = np.asarray(list(((nsr2_data["mu"]))))
                sigma = np.asarray(list(((nsr2_data["sigma"]))))
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
            print "N", str(n), "Unique", str(len(uniqueSites))
            randomSites = np.random.choice(uniqueSites, size=n, replace=False)

            #print len(np.unique(site)),  len(siteNSR2)
            for enumSite, randomSite in enumerate(randomSites):
                for p, q in enumerate(siteNSR2):
                    if q == randomSite:
                        if method == 'geom':
                            print>> OUT2, count1, N[p], S[p], NmaxObs[p], NmaxPred[p], \
                                evennessObs[p], evennessPred[p], skewnessObs[p], skewnessPred[p], R2[p]
                        elif method == 'lognorm':
                            print>> OUT2, count1, N[p], S[p], NmaxObs[p], NmaxPred[p], \
                                evennessObs[p], evennessPred[p], skewnessObs[p], \
                                skewnessPred[p], mu[p], sigma[p], ll[p], R2[p]
                        elif method == 'mete':
                            print>> OUT2, count1, N[p], S[p], NmaxObs[p], NmaxPred[p], \
                                evennessObs[p], evennessPred[p], skewnessObs[p], \
                                skewnessPred[p], p_mete[p], ll[p], R2[p]
                        elif method == 'zipf':
                            print>> OUT2, count1, N[p], S[p], NmaxObs[p], NmaxPred[p], \
                                evennessObs[p], evennessPred[p], skewnessObs[p], \
                                skewnessPred[p], gamma[p], ll[p], R2[p]
                #print site
                #print randomSites
                for r, s in enumerate(site):
                    if s == randomSite:
                        obs2.append(obs[r])
                        pred2.append(pred[r])
                        site2.append(s)
                        print>> OUT1, count1, obs[r], pred[r]
                count1 += 1
            print method, dataset

            obs = np.asarray(obs2)
            pred = np.asarray(pred2)
            zippedSiteObsPred = zip(site2,obs2,pred2)
            #k_minus1 = obs2[0]

        OUT1.close()
        OUT2.close()

def stratifyDataBootstrap(zipfType = 'mle', iterations = 10000,  \
    lognormType = 'pln', remove = True, data_dir= mydir, remove_obs = 0, seqSim = False):
    # do this for the NSR2 only
    # size of each dataset
    modelSlopes = [0.69744476332, 0.851083946706, 1.04353269385, 1.07279674834]
    modelInterepts = [0.386418541599, 0.295532342646, 1.05321973517, 0.139862729426]
    if seqSim == False:
        datasets = ['EMPclosed','HMP', 'MGRAST']
    else:
        datasets = ['SeqSim']
    methods = ['geom', 'mete', 'zipf', 'lognorm']
    #methods = ['lognorm']
    # Number of lines in each file

    MGRAST_sites = 1174
    HMP_sites = 4504
    EMPclosed_sites = 14979
    Total = MGRAST_sites + HMP_sites + EMPclosed_sites
    for i, method in enumerate(methods):
        if remove_obs == 0 and seqSim == False:
            totalSADs = 200
        elif method == 'zipf':
            totalSADs = 40
        else:
            totalSADs = 100
        print method
        count_iter = iterations
        if seqSim == False:
            if method == 'zipf':
                if remove_obs == 0:
                    OUT2 = open(data_dir + 'data/NSR2/Stratified_Test/'+ method  + '_' + zipfType +'_NSR2_stratify.txt', 'wr')
                else:
                    OUT2 = open(data_dir + 'data/NSR2/Stratified_Test/Remove_' + str(remove_obs) + \
                        's/'+ method  + '_' + zipfType +'_NSR2_' +  str(remove_obs) + '_stratify.txt', 'wr')
            elif method == 'lognorm':
                if remove_obs == 0:
                    OUT2 = open(data_dir + 'data/NSR2/Stratified_Test/'+ method  + '_' + lognormType +'_NSR2_stratify.txt', 'wr')
                else:
                    OUT2 = open(data_dir + 'data/NSR2/Stratified_Test/Remove_' + str(remove_obs) + \
                        's/'+ method  + '_' + lognormType +'_NSR2_' +  str(remove_obs) + '_stratify.txt', 'wr')
            else:
                if remove_obs == 0:
                    OUT2 = open(data_dir + 'data/NSR2/Stratified_Test/'+ method  +'_NSR2_stratify.txt', 'wr')
                else:
                    OUT2 = open(data_dir + 'data/NSR2/Stratified_Test/Remove_' + str(remove_obs) + \
                        's/'+ method  + '_NSR2_' +  str(remove_obs) + '_stratify.txt', 'wr')
        else:
            if method == 'zipf':
                OUT2 = open(data_dir + 'data/NSR2/Stratified_Test/SequenceSimilarity/'+ method  + '_' + zipfType + '_' +seqSim +'_NSR2_stratify.txt', 'wr')
            elif method == 'lognorm':
                OUT2 = open(data_dir + 'data/NSR2/Stratified_Test/SequenceSimilarity/'+ method  + '_' + lognormType + '_' +seqSim  +'_NSR2_stratify.txt', 'wr')
            else:
                OUT2 = open(data_dir + 'data/NSR2/Stratified_Test/SequenceSimilarity/'+ method + '_' +seqSim +'_NSR2_stratify.txt', 'wr')

        count2 = 0
        count1 = 0
        N_iter = []
        S_iter = []
        NmaxObs_iter = []
        NmaxPred_iter = []
        evennessObs_iter = []
        evennessPred_iter = []
        skewnessObs_iter = []
        skewnessPred_iter = []
        R2_iter = []
        R2_Nmax_iter = []
        if method != 'geom':
            ll_iter = []
        if method  == 'mete':
            p_mete_iter = []
        if method == 'zipf':
            gamma_iter = []
        if method == 'lognorm':
            mu_iter = []
            sigma_iter = []


        slope1 = []
        intercept1 = []
        r_value1 = []
        p_value1  = []
        std_err1 = []
        slope2 = []
        intercept2 = []
        r_value2 = []
        p_value2  = []
        std_err2 = []
        slope3 = []
        intercept3 = []
        r_value3 = []
        p_value3  = []
        std_err3 = []
        R2_std_iter = []
        for iteration in range(0, iterations):
            for j, dataset in enumerate(datasets):
                lineCount = 0
                removeSADs = []
                if remove_obs == 0 and seqSim == False:
                    get_bad_zipfs = importData.import_NSR2_data(data_dir + 'data/NSR2/' + 'zipf' + '_'+ 'mle' +'_'+dataset +'_NSR2.txt')
                elif seqSim != False:
                    get_bad_zipfs = importData.import_NSR2_data(data_dir + 'data/NSR2/' + 'zipf_mle_MGRAST' + seqSim +'_NSR2.txt')
                else:
                    get_bad_zipfs = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_' +str(remove_obs) + 's/zipf_mle_'+dataset  +'_NSR2_' +  str(remove_obs) +'.txt')
                site = np.asarray(list(((get_bad_zipfs["site"]))))
                N = np.asarray(list(((get_bad_zipfs["N"]))))
                S = np.asarray(list(((get_bad_zipfs["S"]))))
                r2s = np.asarray(list(((get_bad_zipfs["R2"]))))
                zipNsite = zip(site, N, S, r2s)
                for x in zipNsite:
                    if remove_obs != 0:
                        continue
                    else:
                        if float(x[3]) < 0.2:
                            removeSADs.append(int(x[0]))
                removeSADs = np.asarray(removeSADs)

                n = totalSADs
                if remove_obs == 0:
                    if method == 'zipf':
                        if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                            ( str(dataset) == 'EMPopen'):
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method + '_'+ zipfType +'_'+ dataset+'_NSR2.txt')
                        elif str(dataset) == 'MGRAST':
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+ '_'+ zipfType  + '_MGRAST_NSR2.txt')
                        elif seqSim != False:
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+ '_'+ zipfType  + '_MGRAST'+ seqSim +'_NSR2.txt')
                        else:
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+ '_'+ zipfType  + '_MGRAST'+dataset+'_NSR2.txt')
                    elif method == 'lognorm':
                        if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                            ( str(dataset) == 'EMPopen'):
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method + '_'+ lognormType +'_'+ dataset+'_NSR2.txt')
                        elif str(dataset) == 'MGRAST':
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+ '_'+ lognormType  + '_MGRAST_NSR2.txt')
                        elif seqSim != False:
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+ '_'+ lognormType  + '_MGRAST'+ seqSim +'_NSR2.txt')
                        else:
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+ '_'+ lognormType  + '_MGRAST'+dataset+'_NSR2.txt')
                    else:
                        if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                            ( str(dataset) == 'EMPopen'):
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+'_'+dataset+'_NSR2.txt')
                        elif str(dataset) == 'MGRAST':
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+'_MGRAST_NSR2.txt')
                        elif seqSim != False:
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method + '_MGRAST'+ seqSim +'_NSR2.txt')
                        else:
                            nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/' + method+'_MGRAST'+dataset+'_NSR2.txt')
                else:
                    if method == 'zipf':
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_' \
                            + str(remove_obs) + 's/' + method + '_'+ zipfType +'_'+dataset+'_NSR2_' +  str(remove_obs) +'.txt')
                    elif method == 'lognorm':
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_' \
                            + str(remove_obs) + 's/' + method + '_'+ lognormType +'_'+dataset+'_NSR2_' +  str(remove_obs) +'.txt')
                    else:
                        nsr2_data = importData.import_NSR2_data(data_dir + 'data/NSR2/Remove_' \
                            + str(remove_obs) + 's/' + method +'_'+dataset+'_NSR2_' +  str(remove_obs) +'.txt')

                siteNSR2 = np.asarray( map(float, list(((nsr2_data["site"])))))
                N = np.asarray(list(((nsr2_data["N"]))))
                S = np.asarray(list(((nsr2_data["S"]))))
                NmaxObs = np.asarray(list(((nsr2_data["NmaxObs"]))))
                NmaxPred = np.asarray(list(((nsr2_data["NmaxPred"]))))
                if seqSim == False:
                    evennessObs = np.asarray(list(((nsr2_data["evennessObs"]))))
                    evennessPred = np.asarray(list(((nsr2_data["evennessPred"]))))
                    skewnessObs = np.asarray(list(((nsr2_data["skewnessObs"]))))
                    skewnessPred = np.asarray(list(((nsr2_data["skewnessPred"]))))
                    R2 = np.asarray(list(((nsr2_data["R2"]))))

                    if method == 'zipf':
                        gamma = np.asarray(list(((nsr2_data["gamma"]))))
                        ll = np.asarray(list(((nsr2_data["ll"]))))
                    elif method == 'mete':
                        p_mete = np.asarray(list(((nsr2_data["p"]))))
                        ll = np.asarray(list(((nsr2_data["ll"]))))
                    elif method == 'lognorm':
                        mu = np.asarray(list(((nsr2_data["mu"]))))
                        sigma = np.asarray(list(((nsr2_data["sigma"]))))
                        ll = np.asarray(list(((nsr2_data["ll"]))))

                else:
                    R2 = np.asarray(list(((nsr2_data["R2"]))))

                obs2 = []
                pred2 = []
                site2 = []
                if remove == True:
                    siteNSR2_cleaned = np.setdiff1d(siteNSR2, removeSADs)
                    uniqueSites = np.unique(siteNSR2_cleaned)

                else:
                    uniqueSites = np.unique(siteNSR2)
                randomSites = np.random.choice(uniqueSites, size=n, replace=False)

                N_iter_sample = []
                S_iter_sample = []
                NmaxObs_iter_sample = []
                NmaxPred_iter_sample = []
                evennessObs_iter_sample = []
                evennessPred_iter_sample = []
                skewnessObs_iter_sample = []
                skewnessPred_iter_sample = []
                R2_iter_sample = []

                if method != 'geom':
                    ll_iter_sample = []
                if method == 'zipf':
                    gamma_iter_sample = []
                elif method == 'mete':
                    p_mete_iter_sample = []
                elif method  == 'lognorm':
                    mu_iter_sample = []
                    sigma_iter_sample = []

                for enumSite, randomSite in enumerate(randomSites):
                    for p, q in enumerate(siteNSR2):

                        if q == randomSite :
                            if seqSim == False:
                                N_iter_sample.append( N[p])
                                S_iter_sample.append(S[p])
                                NmaxObs_iter_sample.append(NmaxObs[p])
                                NmaxPred_iter_sample.append( NmaxPred[p])
                                evennessObs_iter_sample.append(evennessObs[p])
                                evennessPred_iter_sample.append(evennessPred[p])
                                skewnessObs_iter_sample.append(skewnessObs[p])
                                skewnessPred_iter_sample.append(skewnessPred[p])
                                R2_iter_sample.append(R2[p])
                                if method != 'geom':
                                    ll_iter_sample.append(ll[p])
                                if method =='zipf':
                                    gamma_iter_sample.append(gamma[p])
                                elif method == 'mete':
                                    p_mete_iter_sample.append(p_mete[p])
                                elif method == 'lognorm':
                                    mu_iter_sample.append(mu[p])
                                    sigma_iter_sample.append(sigma[p])

                            else:
                                N_iter_sample.append( N[p])
                                S_iter_sample.append(S[p])
                                R2_iter_sample.append(R2[p])

                if seqSim == False:
                    N_iter.append(np.mean(N_iter_sample))
                    S_iter.append(np.mean(S_iter_sample))
                    NmaxObs_iter.append(np.mean(NmaxObs_iter_sample))
                    NmaxPred_iter.append(np.mean(NmaxPred_iter_sample))
                    evennessObs_iter.append(np.mean(evennessObs_iter_sample))
                    evennessPred_iter.append(np.mean(evennessPred_iter_sample))
                    skewnessObs_iter.append(np.mean(skewnessObs_iter_sample))
                    skewnessPred_iter.append(np.mean(skewnessPred_iter_sample))
                    R2_iter.append(np.mean(R2_iter_sample))
                    R2_std_iter.append(np.std(R2_iter_sample))
                    if method == 'zipf':
                        ll_iter.append(np.mean(ll_iter_sample))
                        gamma_iter.append(np.mean(gamma_iter_sample))
                    elif method == 'mete':
                        ll_iter.append(np.mean(ll_iter_sample))
                        p_mete_iter.append(np.mean(p_mete_iter_sample))
                    elif method == 'lognorm':
                        ll_iter.append(np.mean(ll_iter_sample))
                        mu_iter.append(np.mean(mu_iter_sample))
                        sigma_iter.append(np.mean(sigma_iter_sample))

                    slope1_iter, intercept1_iter, r_value1_iter, p_value1_iter, std_err1_iter = \
                        stats.linregress(np.log10(N_iter_sample),np.log10(NmaxPred_iter_sample))
                    slope2_iter, intercept2_iter, r_value2_iter, p_value2_iter, std_err2_iter = \
                        stats.linregress(np.log10(N_iter_sample),np.log10(evennessPred_iter_sample))
                    slope3_iter, intercept3_iter, r_value3_iter, p_value3_iter, std_err3_iter = \
                        stats.linregress(np.log10(N_iter_sample),np.log10(skewnessPred_iter_sample))

                    slope1.append(slope1_iter)
                    intercept1.append(intercept1_iter)
                    r_value1.append(r_value1_iter)
                    p_value1.append(p_value1_iter)
                    std_err1.append(std_err2_iter)
                    slope2.append(slope2_iter)
                    intercept2.append(intercept2_iter)
                    r_value2.append(r_value2_iter)
                    p_value2.append(p_value2_iter)
                    std_err2.append(std_err2_iter)
                    slope3.append(slope3_iter)
                    intercept3.append(intercept3_iter)
                    r_value3.append(r_value3_iter)
                    p_value3.append(p_value3_iter)
                    std_err3.append(std_err3_iter)

                    scaling_NmaxPred_iter = []
                    for m in range(len(N_iter_sample)):
                        scaling_NmaxPred_iter_i = mo.predictNmax(N_iter_sample[m]).getNmax(b = (10 ** modelInterepts[i]), slope = modelSlopes[i])
                        scaling_NmaxPred_iter.append(scaling_NmaxPred_iter_i)
                    scaling_NmaxPred_iter = np.asarray(scaling_NmaxPred_iter)

                    scaling_NmaxPred_obs_iter = [k for k in zip(NmaxObs_iter_sample, scaling_NmaxPred_iter) if k[0] < 200000 ]
                    NmaxObs_iter_sample_clean = np.asarray([k[0] for k in scaling_NmaxPred_obs_iter])
                    scaling_NmaxPred_iter_clean = np.asarray([k[1] for k in scaling_NmaxPred_obs_iter])
                    R2_Nmax_iter_sample = macroecotools.obs_pred_rsquare(np.log10(NmaxObs_iter_sample_clean), np.log10(scaling_NmaxPred_iter_clean))
                    R2_Nmax_iter.append(R2_Nmax_iter_sample)


                else:
                    N_iter.append(np.mean(N_iter_sample))
                    S_iter.append(np.mean(S_iter_sample))
                    R2_iter.append(np.mean(R2_iter_sample))
                    R2_std_iter.append(np.std(R2_iter_sample))

                count1 += 1

            print str(count_iter)  + " iteration(s) to go!"
            count_iter -=1

        for k in range(0, iterations):
            if seqSim == False:
                if method == 'geom':
                    print>> OUT2, k, int(np.mean([N_iter[k], N_iter[k+iterations], N_iter[k + (2* iterations)]] )), \
                                    int(np.mean([S_iter[k], S_iter[k+iterations], S_iter[k + (2* iterations)] ])), \
                                    int(np.mean([NmaxObs_iter[k], NmaxObs_iter[k+iterations], NmaxObs_iter[k + (2* iterations)] ])), \
                                    int(np.mean([NmaxPred_iter[k], NmaxPred_iter[k+iterations], NmaxPred_iter[k + (2* iterations)] ])),\
                                    np.mean([evennessObs_iter[k], evennessObs_iter[k+iterations], evennessObs_iter[k + (2* iterations)] ]), \
                                    np.mean([evennessPred_iter[k], evennessPred_iter[k+iterations], evennessPred_iter[k + (2* iterations)] ]), \
                                    np.mean([skewnessObs_iter[k], skewnessObs_iter[k+iterations], skewnessObs_iter[k + (2* iterations)] ]), \
                                    np.mean([skewnessPred_iter[k], skewnessPred_iter[k+iterations], skewnessPred_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_iter[k], R2_iter[k+iterations], R2_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_std_iter[k], R2_std_iter[k+iterations], R2_std_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_Nmax_iter[k], R2_Nmax_iter[k+iterations], R2_Nmax_iter[k + (2* iterations)] ]), \
                                    np.mean([slope1[k], slope1[k+iterations], slope1[k + (2* iterations)] ]), \
                                    np.mean([intercept1[k], intercept1[k+iterations], intercept1[k + (2* iterations)] ]),\
                                    np.mean([r_value1[k], r_value1[k+iterations], r_value1[k + (2* iterations)] ]),\
                                    np.mean([p_value1[k], p_value1[k+iterations], p_value1[k + (2* iterations)] ]),\
                                    np.mean([std_err1[k], std_err1[k+iterations], std_err1[k + (2* iterations)] ]), \
                                    np.mean([slope2[k], slope2[k+iterations], slope2[k + (2* iterations)] ]),\
                                    np.mean([intercept2[k], intercept2[k+iterations], intercept2[k + (2* iterations)] ]),\
                                    np.mean([r_value2[k], r_value2[k+iterations], r_value2[k + (2* iterations)] ]),\
                                    np.mean([p_value2[k], p_value2[k+iterations], p_value2[k + (2* iterations)] ]),\
                                    np.mean([std_err2[k], std_err2[k+iterations], std_err2[k + (2* iterations)] ]),\
                                    np.mean([slope3[k], slope3[k+iterations], slope3[k + (2* iterations)] ]), \
                                    np.mean([intercept3[k], intercept3[k+iterations], intercept3[k + (2* iterations)] ]),\
                                    np.mean([r_value3[k], r_value3[k+iterations], r_value3[k + (2* iterations)] ]), \
                                    np.mean([p_value3[k], p_value3[k+iterations], p_value3[k + (2* iterations)] ]), \
                                    np.mean([std_err3[k], std_err3[k+iterations], std_err3[k + (2* iterations)] ])
                if method == 'zipf':
                    print>> OUT2, k, int(np.mean([N_iter[k], N_iter[k+iterations], N_iter[k + (2* iterations)]] )), \
                                    int(np.mean([S_iter[k], S_iter[k+iterations], S_iter[k + (2* iterations)] ])), \
                                    int(np.mean([NmaxObs_iter[k], NmaxObs_iter[k+iterations], NmaxObs_iter[k + (2* iterations)] ])), \
                                    int(np.mean([NmaxPred_iter[k], NmaxPred_iter[k+iterations], NmaxPred_iter[k + (2* iterations)] ])),\
                                    np.mean([evennessObs_iter[k], evennessObs_iter[k+iterations], evennessObs_iter[k + (2* iterations)] ]), \
                                    np.mean([evennessPred_iter[k], evennessPred_iter[k+iterations], evennessPred_iter[k + (2* iterations)] ]), \
                                    np.mean([skewnessObs_iter[k], skewnessObs_iter[k+iterations], skewnessObs_iter[k + (2* iterations)] ]), \
                                    np.mean([skewnessPred_iter[k], skewnessPred_iter[k+iterations], skewnessPred_iter[k + (2* iterations)] ]), \
                                    np.mean([gamma_iter[k], gamma_iter[k+iterations], gamma_iter[k + (2* iterations)] ]), \
                                    np.mean([ll_iter[k], ll_iter[k+iterations], ll_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_iter[k], R2_iter[k+iterations], R2_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_std_iter[k], R2_std_iter[k+iterations], R2_std_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_Nmax_iter[k], R2_Nmax_iter[k+iterations], R2_Nmax_iter[k + (2* iterations)] ]), \
                                    np.mean([slope1[k], slope1[k+iterations], slope1[k + (2* iterations)] ]), \
                                    np.mean([intercept1[k], intercept1[k+iterations], intercept1[k + (2* iterations)] ]),\
                                    np.mean([r_value1[k], r_value1[k+iterations], r_value1[k + (2* iterations)] ]),\
                                    np.mean([p_value1[k], p_value1[k+iterations], p_value1[k + (2* iterations)] ]),\
                                    np.mean([std_err1[k], std_err1[k+iterations], std_err1[k + (2* iterations)] ]), \
                                    np.mean([slope2[k], slope2[k+iterations], slope2[k + (2* iterations)] ]),\
                                    np.mean([intercept2[k], intercept2[k+iterations], intercept2[k + (2* iterations)] ]),\
                                    np.mean([r_value2[k], r_value2[k+iterations], r_value2[k + (2* iterations)] ]),\
                                    np.mean([p_value2[k], p_value2[k+iterations], p_value2[k + (2* iterations)] ]),\
                                    np.mean([std_err2[k], std_err2[k+iterations], std_err2[k + (2* iterations)] ]),\
                                    np.mean([slope3[k], slope3[k+iterations], slope3[k + (2* iterations)] ]), \
                                    np.mean([intercept3[k], intercept3[k+iterations], intercept3[k + (2* iterations)] ]),\
                                    np.mean([r_value3[k], r_value3[k+iterations], r_value3[k + (2* iterations)] ]), \
                                    np.mean([p_value3[k], p_value3[k+iterations], p_value3[k + (2* iterations)] ]), \
                                    np.mean([std_err3[k], std_err3[k+iterations], std_err3[k + (2* iterations)] ])
                elif method == 'mete':
                    print>> OUT2, k, int(np.mean([N_iter[k], N_iter[k+iterations], N_iter[k + (2* iterations)]] )), \
                                    int(np.mean([S_iter[k], S_iter[k+iterations], S_iter[k + (2* iterations)] ])), \
                                    int(np.mean([NmaxObs_iter[k], NmaxObs_iter[k+iterations], NmaxObs_iter[k + (2* iterations)] ])), \
                                    int(np.mean([NmaxPred_iter[k], NmaxPred_iter[k+iterations], NmaxPred_iter[k + (2* iterations)] ])),\
                                    np.mean([evennessObs_iter[k], evennessObs_iter[k+iterations], evennessObs_iter[k + (2* iterations)] ]), \
                                    np.mean([evennessPred_iter[k], evennessPred_iter[k+iterations], evennessPred_iter[k + (2* iterations)] ]), \
                                    np.mean([skewnessObs_iter[k], skewnessObs_iter[k+iterations], skewnessObs_iter[k + (2* iterations)] ]), \
                                    np.mean([skewnessPred_iter[k], skewnessPred_iter[k+iterations], skewnessPred_iter[k + (2* iterations)] ]), \
                                    np.mean([p_mete_iter[k], p_mete_iter[k+iterations], p_mete_iter[k + (2* iterations)] ]), \
                                    np.mean([ll_iter[k], ll_iter[k+iterations], ll_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_iter[k], R2_iter[k+iterations], R2_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_std_iter[k], R2_std_iter[k+iterations], R2_std_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_Nmax_iter[k], R2_Nmax_iter[k+iterations], R2_Nmax_iter[k + (2* iterations)] ]), \
                                    np.mean([slope1[k], slope1[k+iterations], slope1[k + (2* iterations)] ]), \
                                    np.mean([intercept1[k], intercept1[k+iterations], intercept1[k + (2* iterations)] ]),\
                                    np.mean([r_value1[k], r_value1[k+iterations], r_value1[k + (2* iterations)] ]),\
                                    np.mean([p_value1[k], p_value1[k+iterations], p_value1[k + (2* iterations)] ]),\
                                    np.mean([std_err1[k], std_err1[k+iterations], std_err1[k + (2* iterations)] ]), \
                                    np.mean([slope2[k], slope2[k+iterations], slope2[k + (2* iterations)] ]),\
                                    np.mean([intercept2[k], intercept2[k+iterations], intercept2[k + (2* iterations)] ]),\
                                    np.mean([r_value2[k], r_value2[k+iterations], r_value2[k + (2* iterations)] ]),\
                                    np.mean([p_value2[k], p_value2[k+iterations], p_value2[k + (2* iterations)] ]),\
                                    np.mean([std_err2[k], std_err2[k+iterations], std_err2[k + (2* iterations)] ]),\
                                    np.mean([slope3[k], slope3[k+iterations], slope3[k + (2* iterations)] ]), \
                                    np.mean([intercept3[k], intercept3[k+iterations], intercept3[k + (2* iterations)] ]),\
                                    np.mean([r_value3[k], r_value3[k+iterations], r_value3[k + (2* iterations)] ]), \
                                    np.mean([p_value3[k], p_value3[k+iterations], p_value3[k + (2* iterations)] ]), \
                                    np.mean([std_err3[k], std_err3[k+iterations], std_err3[k + (2* iterations)] ])
                elif method == 'lognorm':
                    print>> OUT2, k, int(np.mean([N_iter[k], N_iter[k+iterations], N_iter[k + (2* iterations)]] )), \
                                    int(np.mean([S_iter[k], S_iter[k+iterations], S_iter[k + (2* iterations)] ])), \
                                    int(np.mean([NmaxObs_iter[k], NmaxObs_iter[k+iterations], NmaxObs_iter[k + (2* iterations)] ])), \
                                    int(np.mean([NmaxPred_iter[k], NmaxPred_iter[k+iterations], NmaxPred_iter[k + (2* iterations)] ])),\
                                    np.mean([evennessObs_iter[k], evennessObs_iter[k+iterations], evennessObs_iter[k + (2* iterations)] ]), \
                                    np.mean([evennessPred_iter[k], evennessPred_iter[k+iterations], evennessPred_iter[k + (2* iterations)] ]), \
                                    np.mean([skewnessObs_iter[k], skewnessObs_iter[k+iterations], skewnessObs_iter[k + (2* iterations)] ]), \
                                    np.mean([skewnessPred_iter[k], skewnessPred_iter[k+iterations], skewnessPred_iter[k + (2* iterations)] ]), \
                                    np.mean([mu_iter[k], mu_iter[k+iterations], mu_iter[k + (2* iterations)] ]), \
                                    np.mean([sigma_iter[k], sigma_iter[k+iterations], sigma_iter[k + (2* iterations)] ]), \
                                    np.mean([ll_iter[k], ll_iter[k+iterations], ll_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_iter[k], R2_iter[k+iterations], R2_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_std_iter[k], R2_std_iter[k+iterations], R2_std_iter[k + (2* iterations)] ]), \
                                    np.mean([R2_Nmax_iter[k], R2_Nmax_iter[k+iterations], R2_Nmax_iter[k + (2* iterations)] ]), \
                                    np.mean([slope1[k], slope1[k+iterations], slope1[k + (2* iterations)] ]), \
                                    np.mean([intercept1[k], intercept1[k+iterations], intercept1[k + (2* iterations)] ]),\
                                    np.mean([r_value1[k], r_value1[k+iterations], r_value1[k + (2* iterations)] ]),\
                                    np.mean([p_value1[k], p_value1[k+iterations], p_value1[k + (2* iterations)] ]),\
                                    np.mean([std_err1[k], std_err1[k+iterations], std_err1[k + (2* iterations)] ]), \
                                    np.mean([slope2[k], slope2[k+iterations], slope2[k + (2* iterations)] ]),\
                                    np.mean([intercept2[k], intercept2[k+iterations], intercept2[k + (2* iterations)] ]),\
                                    np.mean([r_value2[k], r_value2[k+iterations], r_value2[k + (2* iterations)] ]),\
                                    np.mean([p_value2[k], p_value2[k+iterations], p_value2[k + (2* iterations)] ]),\
                                    np.mean([std_err2[k], std_err2[k+iterations], std_err2[k + (2* iterations)] ]),\
                                    np.mean([slope3[k], slope3[k+iterations], slope3[k + (2* iterations)] ]), \
                                    np.mean([intercept3[k], intercept3[k+iterations], intercept3[k + (2* iterations)] ]),\
                                    np.mean([r_value3[k], r_value3[k+iterations], r_value3[k + (2* iterations)] ]), \
                                    np.mean([p_value3[k], p_value3[k+iterations], p_value3[k + (2* iterations)] ]), \
                                    np.mean([std_err3[k], std_err3[k+iterations], std_err3[k + (2* iterations)] ])
            else:
                print>> OUT2, k, int(N_iter[k]), \
                                int(S_iter[k]), \
                                R2_iter[k], \
                                R2_std_iter[k]
        OUT2.close()


def subsample_N(data_dir= mydir, SAD_numbers=10, iterations=100):
    percents = [0.500000, 0.250000, 0.125000, 0.062500, 0.031250, 0.015625]
    models = ['geom', 'mete', 'lognorm', 'zipf']
    for model in models:
        OUT = open(data_dir + 'data/Subsample_N/'+ model + '_SubSampled_Data_2.txt', 'w+')
        if model == 'geom':
            IN = pd.DataFrame(importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/geom_obs_pred_stratify.txt'))
        elif model == 'mete':
            IN = pd.DataFrame(importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/mete_obs_pred_stratify.txt'))
        elif model == 'lognorm':
            IN = pd.DataFrame(importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/lognorm_pln_obs_pred_stratify.txt'))
        elif model == 'zipf':
            IN = pd.DataFrame(importData.import_obs_pred_data(data_dir + 'data/ObsPred/Stratified/zipf_mle_obs_pred_stratify.txt'))
        else:
            continue
        sites = set(IN['site'].tolist())
        S = IN.groupby('site')['obs'].count()
        N = IN.groupby('site')['obs'].sum()
        S_and_N = pd.concat([S, N], axis=1)
        S_and_N.columns = ['S', 'N']
        S_and_N = S_and_N[(S_and_N['S'] *  0.015625 > 10) & (S_and_N['N'] > S_and_N['S'] )]
        sites =  S_and_N.index
        sites = map(int, sites)
        iteration_count = iterations
        while (iteration_count > 0):
            sites_sample = map(int, random.sample(sites, SAD_numbers))
            N_0_iter = []
            S_0_iter = []
            N_max_obs_iter = []
            N_max_pred_iter = []
            r2_iter = []

            for site_sample in sites_sample:
                SAD = list(IN.loc[IN.site == site_sample].obs)
                N_0 = sum(SAD)
                S_0 = len(SAD)
                N_max = max(SAD)
                line_ra = map(lambda x: x/sum(SAD), SAD)

                gm_lines = SAD_numbers
                means = [N_0, S_0, N_max]
                failed_percents = 0

                N_0_list_percent = []
                S_0_list_percent = []
                N_max_obs_list_percent = []
                N_max_pred_list_percent = []
                r2_list_percent = []
                for k, percent in enumerate(percents):
                    sample_size = round(percent * N_0)
                    sample_k = np.random.multinomial(sample_size, line_ra, size = None)
                    sample_k_sorted = -np.sort( -sample_k[sample_k != 0] )
                    N_k = sum(sample_k_sorted)
                    S_k = sample_k_sorted.size
                    print model, N_k, S_k, ' countdown: ', percent, iteration_count
                    Nmax_k_obs = max(sample_k_sorted)
                    result_k = mete.get_mete_rad(S_k, N_k)
                    pred_mete = result_k[0]
                    if model == 'geom':
                        pred = mo.get_Geom(N_k, S_k, False) # False mean no zeros allowed
                    elif model == 'mete':
                        result = mete.get_mete_rad(S_k, N_k)
                        pred = result[0]

                    elif model == 'lognorm':
                        lognorm_pred = mo.lognorm(sample_k_sorted, 'pln')
                        pred = lognorm_pred.get_rad_from_obs()[0]

                    elif model == 'zipf':
                        # Start the timer. Once 1 second is over, a SIGALRM signal is sent.
                        signal.alarm(6)
                        # This try/except loop ensures that
                        #   you'll catch TimeoutException when it's sent.
                        try:
                            # Whatever your function that might hang
                            # use S
                            zipf_class = mo.zipf(sample_k_sorted, 'fmin')
                            pred_tuple = zipf_class.from_cdf()
                            pred = pred_tuple[0]
                            gamma = pred_tuple[1]
                        except mo.TimeoutException:
                            continue # continue the for loop if function takes more than x seconds
                        else:
                            # Reset the alarm
                            signal.alarm(0)
                    else:
                        pass

                    Nmax_k_pred = max(pred)
                    r2_k = macroecotools.obs_pred_rsquare(np.log10(sample_k_sorted), np.log10(pred))
                    if any( (r2 == -float('inf') ) or (r2 == float('inf') ) or (r2 == float('Nan') ) for r2 in r2_list_percent):
                        iteration_count += 1
                        continue
                    N_0_list_percent.append(N_k)
                    S_0_list_percent.append(S_k)
                    N_max_obs_list_percent.append(Nmax_k_obs)
                    N_max_pred_list_percent.append(Nmax_k_pred)
                    r2_list_percent.append(r2_k)

                N_0_iter.append(N_0_list_percent)
                S_0_iter.append(S_0_list_percent)
                N_max_obs_iter.append(N_max_obs_list_percent)
                N_max_pred_iter.append(N_max_pred_list_percent)
                r2_iter.append(r2_list_percent)
            iteration_count -= 1
            N_0_zip = zip(*N_0_iter)
            S_0_zip = zip(*S_0_iter)
            N_max_obs_zip = zip(*N_max_obs_iter)
            N_max_pred_zip = zip(*N_max_pred_iter)
            r2_zip = zip(*r2_iter)
            for i, percent in enumerate(percents):
                print>> OUT, percent, np.mean(N_0_zip[i]),np.mean(S_0_zip[i]), \
                    np.mean(N_max_obs_zip[i]),np.mean(N_max_pred_zip[i]), \
                    np.mean(r2_zip[i])
