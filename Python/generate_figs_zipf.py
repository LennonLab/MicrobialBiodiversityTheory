from __future__ import division
import os, sys, signal, collections, argparse, optparse, math
import scipy as sp
import  matplotlib.pyplot as plt
import numpy as np
from scipy import stats, optimize
from scipy.stats import gaussian_kde, norm
from mpl_toolkits.axes_grid.inset_locator import inset_axes
from sys import argv
import pandas as pd
from sklearn.grid_search import GridSearchCV
from sklearn.neighbors import KernelDensity
import macroeco_distributions as md
import macroecotools
import mete
import statsmodels.formula.api as smf
import statsmodels.api as sm


#mydir = os.path.expanduser("~/github/MicroMETE/data/")
mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'


"""This code was written using MIT liscenced code from the following Weecology
repos: METE (https://github.com/weecology/METE) and macroecotools
(https://github.com/weecology/macroecotools). """


class TimeoutException(Exception):   # Custom exception class
    pass

def timeout_handler(signum, frame):   # Custom signal handler
    raise TimeoutException

class zipf:

    """ A class to obtain a zipf object with inherited mle shape parameter,
    mle form of the rank-abundance distribution, and a rank-abundance curve
    based on fitting the zipf to the observed via a generalized linear model."""

    def __init__(self, obs, estimator):
        self.obs = obs
        self.estimator = estimator

    def zipf_ll(self, ab, a):
        """Log-likelihood of the Zipf distribution with x_min = 1."""
        return sum(stats.zipf.logpmf(ab, a))

    def zipf_solver(self, ab):
        #ab = self.obs
        """Obtain the MLE parameter for a Zipf distribution with x_min = 1."""
        par0 = 1 + len(ab) / (sum(np.log(2 * np.array(ab))))
        def zipf_func(x):
            return -self.zipf_ll(ab, x)
        #par = optimize.fmin(zipf_func, x0 = par0, disp=False)[0]
        estimator = str(self.estimator)
        par = getattr(optimize, estimator)(zipf_func, x0 = par0, disp=False)[0]
        return par

    def from_cdf(self):
        """ Obtain the maximum likelihood form of the Zipf distribution, given
        the mle value for the Zipf shape parameter (a). Using a, this code
        generates a rank-abundance distribution (RAD) from the cumulative
        density function (cdf) using the percent point function (ppf) also known
        as the quantile function.
        see: http://www.esapubs.org/archive/ecol/E093/155/appendix-B.htm

        This is an actual form of the Zipf distribution, obtained from getting
        the mle for the shape parameter.
        """
        p = self.zipf_solver(self.obs)
        S = len(self.obs)
        rv = stats.zipf(a=p)
        rad = []
        for i in range(1, S+1):
            val = (S - i + 0.5)/S
            x = rv.ppf(val)
            rad.append(int(x))
        point = collections.namedtuple('Rad_and_p', ['x', 'y'])
        point_return = point(rad, y = p)
        #print point_return.x, point_return.y
        return point_return

    def from_glm(self):

        """ Fit the Zipf distribution to the observed vector of integer values
        using a generalized linear model.

        Note: This is a fitted curve; not an actual form of the Zipf distribution

        This method was inspired by the vegan
        package's open source code on vegan's public GitHub repository:
        https://github.com/vegandevs/vegan/blob/master/R/rad.zipf.R
        on Thursday, 19 Marth 2015 """

        ranks = np.log(range(1, len(self.obs)+1))
        off = [np.log(sum(self.obs))] * len(self.obs)

        d = pd.DataFrame({'ranks': ranks, 'off': off, 'x':self.obs})

        lm = smf.glm(formula='x ~ ranks', data = d, family = sm.families.Poisson()).fit()
        pred = lm.predict()

        return pred

class lognorm:

    def __init__(self, obs):
        self.obs = obs

    def ppoints(self, n):
        """ numpy analogue or `R`'s `ppoints` function
            see details at http://stat.ethz.ch/R-manual/R-patched/library/stats/html/ppoints.html
            :param n: array type or number

            Obtained from: http://stackoverflow.com/questions/20292216/imitating-ppoints-r-function-in-python
            on 5 April 2016
            """
        if n < 10:
            a = 3/8
        else:
            a = 1/2

        try:
            n = np.float(len(n))
        except TypeError:
            n = np.float(n)
        return (np.arange(n) + 1 - a)/(n + 1 - 2*a)


    def lognorm_glm(self):

        """ Fit the lognormal distribution to the observed vector of integer values
        using a generalized linear model.

        Note: This is a fitted curve; not an actual form of the lognormal distribution

        This method was inspired by the vegan package's open source code on vegan's public
        GitHub repository: https://github.com/vegandevs/vegan/R/rad.lognormal.R
        on Thursday, 5 April 2016
        """

        ranks = np.log(range(1, len(self.obs)+1))
        ranks = -norm.ppf(self.ppoints(len(ranks)))

        d = pd.DataFrame({'rnks': ranks, 'x': self.obs})
        lm = smf.glm(formula='x ~ rnks', data = d, family = sm.genmod.families.family.Poisson(link=sm.genmod.families.links.log)).fit()
        pred = lm.predict()

        return pred

def e_simpson(SAD): # based on 1/D, not 1 - D

    " Simpson's evenness "
    SAD = filter(lambda a: a != 0, SAD)

    D = 0.0
    N = sum(SAD)
    S = len(SAD)

    for x in SAD:
        D += (x*x) / (N*N)

    E = round((1.0/D)/S, 4)

    if E < 0.0 or E > 1.0:
        print 'Simpsons Evenness =',E
    return E

def skewness(RAD):
    skew = stats.skew(RAD)
    # log-modulo skewness
    lms = np.log10(np.abs(float(skew)) + 1)
    if skew < 0:
        lms = lms * -1
    return lms

def get_SADs_mgrast(path, thresholds):

    datasets = ['BOVINE', 'CATLIN', 'CHU', 'HYDRO', 'LAUB']

    for t in thresholds:
        SADdict = {}

        for d in datasets:
            name = d+t
            print name

            filepath  = path + 'MGRAST-Data/'+t+'/'+name+'/'+name+'-data.txt'
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

def get_SADs_mgrast_test(path):
    SADdict = {}
    fungi_list = map(str,np.arange(4484945.3, 4485075.3,1))
    filepath  = path + 'MGRAST-Data/MGRAST/MGRAST-data.txt'
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

def get_SADs_HMP(path, timeseries):
    timeseries = bool(timeseries)
    if timeseries == True:
        IN = path + 'HMP-Data/HMPsparseSbyS_NAP.txt'
        OUT =  open(path+'HMP-Data/' + 'HMP-SADs_NAP.txt', 'w+')
    else:
        IN = path + 'HMP-Data/HMPsparseSbyS_NAP_noTimeseries.txt'
        OUT =  open(path+'HMP-Data/' + 'HMP-SADs_NAP_noTimeseries.txt', 'w+')
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
                                #print 'unclassified'
                                continue
                            elif 'unidentified' in i:
                                #print 'unidentified'
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

def EMP_SADs(path, name, mgrast):
    minS = 10
    IN = path + '/' + name + '-SSADdata.txt'
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

def get_GeomSeries(N,S,zeros):

    rank = range(1,S+1)
    cdf = [(S-i+0.5)/S for i in rank]
    SNratio = S/N
    if zeros == False:
        abd = md.trunc_geom.ppf(np.array(cdf), SNratio, N)
    return abd

def generate_obs_pred_data(datasets, methods, size = 0, remove = 0, zipfType = 'glm'):
    remove = int(remove)
    if remove != 0:
        newpath1 = mydir + "ObsPred/Remove_" + str(remove) + 's/'
        if not os.path.exists(newpath1):
            os.makedirs(newpath1)
        newpath2 = mydir + "NSR2/Remove_" + str(remove) + 's/'
        if not os.path.exists(newpath2):
            os.makedirs(newpath2)
    for method in methods:
        for dataset in datasets:

            signal.signal(signal.SIGALRM, timeout_handler)

            if (method != 'zipf' and dataset != 'MGRAST'):

                if dataset == 'EMPclosed' or dataset == 'EMPopen' :
                    IN = mydir  + dataset + '-Data' + '/' + dataset +'-SADs.txt'
                    if remove == 0:
                        OUT1 = open(mydir + "ObsPred/" + method +'_'+dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "NSR2/" + method +'_'+dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "ObsPred/Remove_" + str(remove)  + 's/'+ method +'_'+dataset +'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "NSR2/Remove_" + str(remove)  + 's/'+ method +'_'+dataset+'_NSR2_' \
                            + str(remove) + '.txt','w+')
                elif dataset == "HMP":
                    IN = mydir  + dataset + '-Data' + '/' + dataset +'-SADs_NAP.txt'
                    if remove == 0:
                        OUT1 = open(mydir + "ObsPred/" + method +'_'+dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "NSR2/" + method +'_'+dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "ObsPred/Remove_" + str(remove)  + 's/'+ method +'_'+dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "NSR2/Remove_" + str(remove)  + 's/'+ method +'_'+dataset+'_NSR2_' \
                            + str(remove) + '.txt','w+')
                else:
                    IN = mydir + 'MGRAST-Data/' + dataset +  '/' + 'MGRAST-' + dataset + '-SADs.txt'
                    if remove == 0:
                        OUT1 = open(mydir + "ObsPred/" + method +'_'+ 'MGRAST' + dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "NSR2/" + method +'_'+ 'MGRAST' + dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "ObsPred/Remove_" + str(remove)  + 's/'+ method +'_'+ 'MGRAST' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "NSR2/Remove_" + str(remove)  + 's/'+ method +'_'+ 'MGRAST' + dataset+'_NSR2_'  \
                            + str(remove) + '.txt','w+')
            elif (method == 'zipf' and dataset != 'MGRAST'):

                if dataset == 'EMPclosed' or dataset == 'EMPopen' or dataset == 'HMP':
                    IN = mydir  + dataset + '-Data' + '/' + dataset +'-SADs_NAP.txt'
                    if remove == 0:
                        OUT1 = open(mydir + "ObsPred/" + method + '_' + zipfType + '_'+dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "NSR2/" + method + '_' + zipfType +'_'+dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "ObsPred/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_'+ dataset+'_NSR2_'  \
                            + str(remove) + '.txt','w+')
                else:
                    IN = mydir + 'MGRAST-Data/' + dataset +  '/' + 'MGRAST-' + dataset + '-SADs.txt'
                    if remove == 0:
                        OUT1 = open(mydir + "ObsPred/" + method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_obs_pred.txt','w+')
                        OUT2 = open(mydir + "NSR2/" + method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "ObsPred/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_'+ 'MGRAST' + dataset+'_NSR2_' \
                            + str(remove) + '.txt','w+')

            elif dataset == 'MGRAST':
                IN = mydir + 'MGRAST-Data/MGRAST/MGRAST-SADs.txt'
                if remove == 0:
                    if method == 'zipf':
                        OUT1 = open(mydir + "ObsPred/" + method + '_' + zipfType +'_'+ 'MGRAST_obs_pred.txt','w+')
                        OUT2 = open(mydir + "NSR2/" + method + '_' + zipfType +'_'+ 'MGRAST_NSR2.txt','w+')
                    else:
                        OUT1 = open(mydir + "ObsPred/" + method +'_'+ 'MGRAST_obs_pred.txt','w+')
                        OUT2 = open(mydir + "NSR2/" + method +'_'+ 'MGRAST_NSR2.txt','w+')
                else:
                    if method == 'zipf':
                        OUT1 = open(mydir + "ObsPred/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "NSR2/Remove_" + str(remove)  + 's/'+ method + '_' + zipfType +'_' + dataset+'_NSR2_' \
                            + str(remove) + '.txt','w+')
                    else:
                        OUT1 = open(mydir + "ObsPred/Remove_" + str(remove)  + 's/'+ method +'_' + dataset+'_obs_pred_' \
                            + str(remove) + '.txt','w+')
                        OUT2 = open(mydir + "NSR2/Remove_" + str(remove)  + 's/'+ method +'_' + dataset+'_NSR2_' \
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
                evennessObs = e_simpson(obs)
                skewnessObs = skewness(obs)


                if S < 10 or N <= S:
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
                    pred = get_GeomSeries(N, S, False) # False mean no zeros allowed

                elif method == 'mete': # Predicted log-series
                    logSeries = mete.get_mete_rad(S, N)
                    pred = logSeries[0]
                elif method == 'lognorm':
                    lognorm_pred = lognorm(obs)
                    pred = lognorm_pred.lognorm_glm()
                    pred = np.ceil(pred)
                    pred.astype(int)

                elif method == 'zipf' and zipfType == 'glm':
                        zipf_pred = zipf(obs, 'fmin')
                        zipf_glm = zipf_pred.from_glm()
                        pred = np.ceil(zipf_glm)
                        pred.astype(int)
                        print line_count
                        #numpy.ceil
                elif method == 'zipf' and zipfType != 'glm':
                        #line = map(int, line)
                        # Start the timer. Once 1 second is over, a SIGALRM signal is sent.
                        signal.alarm(3)
                        # This try/except loop ensures that
                        #   you'll catch TimeoutException when it's sent.
                        try:
                            # Whatever your function that might hang
                            Zipf_solve_line = md.zipf_solver(obs)
                            # use S
                            rv = stats.zipf(Zipf_solve_line)
                            zipf_class = zipf(obs, 'fmin')
                            pred_tuple = zipf_class.from_cdf()
                            pred = pred_tuple[0]
                            gamma = pred_tuple[1]
                        except TimeoutException:
                            continue # continue the for loop if function takes more than x seconds
                        else:
                            # Reset the alarm
                            signal.alarm(0)
                NmaxPred = np.amax(pred)
                evennessPred = e_simpson(pred)
                skewnessPred = skewness(pred)
                r2 = macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))
                print " r2:", r2
                if r2 == -float('inf') or r2 == float('inf') or r2 == float('Nan'):
                    print r2 + " is Nan or inf, removing..."
                    continue
                if method == 'zipf':
                    if dataset == 'EMPclosed' or dataset == 'EMPopen' or dataset == 'HMP':
                        OUT2 = open(mydir + "NSR2/" + method + '_' + zipfType +'_'+dataset+'_NSR2.txt','a+')
                        if zipfType == 'glm':
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
                        if zipfType == 'glm':
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




def import_obs_pred_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8", names = ['site','obs','pred'], delimiter = " ")
    #test = data[0:10000]
    #return test
    return data

def import_subsampled_data(input_filename):
    if ('zipf' in input_filename):
        # 33 for zipf
        # this needs to be fixesd, I put the file name twice in old code
        data = np.genfromtxt(input_filename, \
        dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
        names = ['site','site2','N0','S0','Nmax', \
        'N0_05','S0_05','Nmax_05', 'r2_05', 'gamma_05', \
        'N0_025','S0_025','Nmax_025', 'r2_025', 'gamma_025', \
        'N0_0125','S0_0125','Nmax_0125', 'r2_0125', 'gamma_0125', \
        'N0_00625','S0_00625','Nmax_00625', 'r2_00625', 'gamma_00625', \
        'N0_003125','S0_003125','Nmax_003125', 'r2_003125', 'gamma_003125',
        'N0_0015625','S0_0015625','Nmax_0015625','r2_0015625', 'gamma_0015625'], \
        delimiter = " ")
    else:
        # 27 columns for mete and geom
        data = np.genfromtxt(input_filename, \
        dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
        names = ['site','N0','S0','Nmax', \
        'N0_05','S0_05','Nmax_05','r2_05', \
        'N0_025','S0_025','Nmax_025','r2_025', \
        'N0_0125','S0_0125','Nmax_0125','r2_0125', \
        'N0_00625','S0_00625','Nmax_00625','r2_00625', \
        'N0_003125','S0_003125','Nmax_003125','r2_003125', \
        'N0_0015625','S0_0015625','Nmax_0015625','r2_0015625'], \
        delimiter = " ")
    return data

def import_subsampled_data_pandas(input_filename):
    if ('zipf' in input_filename):
        names = ['site','site2','N0','S0','Nmax', \
        'N0_05','S0_05','Nmax_05', 'r2_05', 'gamma_05', \
        'N0_025','S0_025','Nmax_025', 'r2_025', 'gamma_025', \
        'N0_0125','S0_0125','Nmax_0125', 'r2_0125', 'gamma_0125', \
        'N0_00625','S0_00625','Nmax_00625', 'r2_00625', 'gamma_00625', \
        'N0_003125','S0_003125','Nmax_003125', 'r2_003125', 'gamma_003125',
        'N0_0015625','S0_0015625','Nmax_0015625','r2_0015625', 'gamma_0015625']
        #data_table = pd.read_table(input_filename, names = names, header = None, sep=' ')

    else:
        names = ['site','N0','S0','Nmax', \
        'N0_05','S0_05','Nmax_05','r2_05', \
        'N0_025','S0_025','Nmax_025','r2_025', \
        'N0_0125','S0_0125','Nmax_0125','r2_0125', \
        'N0_00625','S0_00625','Nmax_00625','r2_00625', \
        'N0_003125','S0_003125','Nmax_003125','r2_003125', \
        'N0_0015625','S0_0015625','Nmax_0015625','r2_0015625']
    data_table = pd.read_table(input_filename, names = names, header = None, sep=' ')
    return data_table

def hist_mete_r2(sites, obs, pred):  # TAKEN FROM Macroecotools or the mete_sads.py script used for White et al. (2012)
    """Generate a kernel density estimate of the r^2 values for obs-pred plots"""
    r2s = []
    for site in np.unique(sites):
        if int(site) >= 100:
            break
        else:
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



def obs_pred_r2_multi(methods, datasets, data_dir= mydir): # TAKEN FROM THE mete_sads.py script
    print 'generating 1:1 line R-square values for dataset(s)'

    for i, dataset in enumerate(datasets):
        for j, method in enumerate(methods):
            obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + "_" + dataset + '_obs_pred.txt')
            #obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + "_" + dataset + '_obs_pred.txt')
            obs = ((obs_pred_data["obs"]))
            pred = ((obs_pred_data["pred"]))
            print method, dataset,' ',macroecotools.obs_pred_rsquare(np.log10(obs), np.log10(pred))


def import_NSR2_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    input_filename_str = str(input_filename)
    #NSR2_method = input_filename_split[-4]
    #method = str(NSR2_method.split('/')[1])
    if 'Stratified' in input_filename_str:
        data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
            names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                'evennessPred', 'skewnessObs', 'skewnessPred','R2'], delimiter = " ")
    else:
        if 'HMP' in input_filename_str:
            if ('zipf' in input_filename_str) :
                if ('glm' in input_filename_str) :


                    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                        names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                            'evennessPred', 'skewnessObs', 'skewnessPred','R2', 'NAP'], delimiter = " ")
                else:
                    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                        names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                            'evennessPred', 'skewnessObs', 'skewnessPred', 'gamma','R2', 'NAP'], delimiter = " ")
            else:
                data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                    'evennessPred', 'skewnessObs', 'skewnessPred', 'R2','NAP'], delimiter = " ")
        else:
            if 'zipf' in input_filename_str:
                if ('glm' in input_filename_str) :
                    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                    names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                        'evennessPred', 'skewnessObs', 'skewnessPred', 'R2'], delimiter = " ")
                else:
                    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                    names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                        'evennessPred', 'skewnessObs', 'skewnessPred', 'gamma','R2'], delimiter = " ")
            else:
                data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                names = ['site','N','S','NmaxObs', 'NmaxPred', 'evennessObs', \
                    'evennessPred', 'skewnessObs', 'skewnessPred', 'R2'], delimiter = " ")

    return data


def plot_obs_pred_sad(methods, datasets, n, figname = 'Fig1', data_dir=mydir, \
    stratify = True, radius=2, remove = 0, zipfType = 'glm'):
    # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    # Used for Figure 3 Locey and White (2013)
    """Multiple obs-predicted plotter"""
    fig = plt.figure()
    count = 0
    if stratify == True:
        plot_dim = 2
        for i, method in enumerate(methods):
            obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/Stratified/'+ method +'_obs_pred_stratify_t.txt')
            obs = np.asarray(list(((obs_pred_data["obs"]))))
            pred = np.asarray(list(((obs_pred_data["pred"]))))
            site = np.asarray(list(((obs_pred_data["site"]))))

            obs2 = []
            pred2 = []
            site2 = []
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
                axis_max = 10  * max(pred)
            else:
                axis_min = 0
                axis_max = 2 * max(obs)
            ax = fig.add_subplot(plot_dim, plot_dim, count+1)
            ax.set(adjustable='box-forced', aspect='equal')
            if method == 'geom':
                ax.set_title(r"$\mathbf{Broken-stick}$")
            elif method == 'lognorm':
                ax.set_title(r"$\mathbf{Lognormal}$")
            elif method == 'mete':
                ax.set_title(r"$\mathbf{Log-series}$")
            elif method == 'zipf':
                ax.set_title(r"$\mathbf{Zipf}$")

            macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=1,
                            plot_obj=plt.subplot(plot_dim,plot_dim,count+1))

            plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
            plt.xlim(0, axis_max)
            plt.ylim(0, axis_max)

            plt.tick_params(axis='both', which='major', labelsize=7)
            plt.subplots_adjust(wspace=0.5, hspace=0.3)

            axins = inset_axes(ax, width="30%", height="30%", loc=4)
            #if str(dataset) == 'EMPclosed' or str(dataset) == 'EMPopen':
            INh2 = import_NSR2_data(data_dir + 'NSR2/Stratified/' + method + '_NSR2_stratify_t.txt')
            r2s = ((INh2["R2"]))
            hist_r2 = np.histogram(r2s, range=(0, 1))
            xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
            xvals = xvals[0:len(xvals)-1]
            yvals = hist_r2[0]
            plt.plot(xvals, yvals, 'k-', linewidth=2)
            plt.axis([0, 1, 0, 1.1 * max(yvals)])
            #else:
            #    hist_mete_r2(site, np.log10(obs), np.log10(pred))
            plt.setp(axins, xticks=[], yticks=[])


            count += 1
        fig.text(0.50, 0.04, 'Predicted rank-abundance', ha='center', va='center')
        fig.text(0.05, 0.5, 'Observed rank-abundance', ha='center', va='center', rotation='vertical')
        fig_name = str(mydir[:-6]) + '/figures/' + figname + '.png'
        plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
        plt.close()

    else:
        plot_dim = len(datasets)
        for i, dataset in enumerate(datasets):
            for j, method in enumerate(methods):
                if remove == 0:
                    if method == 'zipf':
                        if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                            ( str(dataset) == 'EMPopen'):
                            obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_'+ zipfType +'_'+dataset+'_obs_pred.txt')
                        elif str(dataset) == 'MGRAST':
                            obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_'+ zipfType + '_' + 'MGRAST_obs_pred.txt')
                        else:
                            obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_'+ zipfType + '_' + 'MGRAST' + dataset +'_obs_pred.txt')
                    else:
                        if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                            ( str(dataset) == 'EMPopen'):
                            obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method+'_'+dataset+'_obs_pred.txt')
                        elif str(dataset) == 'MGRAST':
                            obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_' + 'MGRAST_obs_pred.txt')
                        else:
                            obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_' + 'MGRAST' + dataset +'_obs_pred.txt')
                else:
                    if method == 'zipf':
                        if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                            ( str(dataset) == 'EMPopen') or ( str(dataset) == 'MGRAST'):
                            obs_pred_data = data_dir + "ObsPred/Remove_" + str(remove)  + 's/'+ method + '_'+ zipfType +'_'+dataset+'_obs_pred_' \
                                + str(remove) + '.txt'
                        else:
                            obs_pred_data = import_obs_pred_data(data_dir + "ObsPred/Remove_" + str(remove)  + 's/'+ method + '_'+ zipfType +'_'+ 'MGRAST' + dataset+'_obs_pred_' \
                                + str(remove) + '.txt')
                    else:
                        if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                            ( str(dataset) == 'EMPopen') or ( str(dataset) == 'MGRAST'):
                            obs_pred_data = data_dir + "ObsPred/Remove_" + str(remove)  + 's/'+ method +'_'+dataset+'_obs_pred_' \
                                + str(remove) + '.txt'
                        else:
                            obs_pred_data = import_obs_pred_data(data_dir + "ObsPred/Remove_" + str(remove)  + 's/'+ method +'_'+ 'MGRAST' + dataset+'_obs_pred_' \
                                + str(remove) + '.txt')
                print method, dataset

                site = ((obs_pred_data["site"]))
                print "Dataset " + str(dataset) + " for method " + str(method) + " has "+ str(len(set(site))) + " sites"

                obs = list(((obs_pred_data["obs"])))
                pred = list(((obs_pred_data["pred"])))
                obs2 = []
                pred2 = []
                site2 = []

                if n == 'all' or len(obs) <= n:
                    obs2 = list(obs)
                    pred2 = list(pred)
                    site2 = list(site)

                else:
                    if len(obs) > n:
                        inds = np.random.choice(range(len(obs)), size=n, replace=False)
                        for ind in inds:
                            obs2.append(obs[ind])
    	                pred2.append(pred[ind])
    	                site2.append(site[ind])

                obs = np.asarray(obs2)
                pred = np.asarray(pred2)
                site =  np.asarray(site2)
                if method == 'zipf':
                    axis_min = 0
                    axis_max = 10  * max(pred)
                else:
                    axis_min = 0
                    axis_max = 2 * max(obs)

                ax = fig.add_subplot(plot_dim, plot_dim, count+1)
                ax.set(adjustable='box-forced', aspect='equal')

                if j == 0:
                    if all(x in datasets for x in ['HMP','EMPclosed','EMPopen', 'MGRAST']) == True:
                        if dataset == 'HMP':
                            ax.set_ylabel("HMP", rotation=90, size=12)
                        elif dataset == 'EMPclosed':
                            ax.set_ylabel("EMP (closed)", rotation=90, size=12)
                        elif  dataset == 'EMPopen':
                            ax.set_ylabel("EMP (open)", rotation=90, size=12)
                        elif dataset == 'MGRAST':
                            ax.set_ylabel("MG-RAST", rotation=90, size=12)
                    if all(x in datasets for x in ['95', '97', '99']) == True:
                        if dataset == '95':
                            ax.set_ylabel("MG-RAST 95%", rotation=90, size=12)
                        elif dataset == '97':
                            ax.set_ylabel("MG-RAST 97%", rotation=90, size=12)
                        elif dataset == '99':
                            ax.set_ylabel("MG-RAST 99%", rotation=90, size=12)
                    else:
                        if dataset == 'HMP':
                            ax.set_ylabel("HMP", rotation=90, size=12)
                        elif dataset == 'EMPclosed' or method == 'EMPopen':
                            ax.set_ylabel("EMP", rotation=90, size=12)
                        elif dataset == '97':
                            ax.set_ylabel("MG-RAST 97%", rotation=90, size=12)
                        elif dataset == 'MGRAST':
                            ax.set_ylabel("MG-RAST", rotation=90, size=12)
                if i == 0 and j == 0:
                    ax.set_title("Broken-stick")
                elif i == 0 and j == 1:
                    ax.set_title("METE")
                elif i == 0 and j == 2:
                    ax.set_title("Zipf")

                print len(obs), len(pred)
                macroecotools.plot_color_by_pt_dens(pred, obs, radius, loglog=1,
                                plot_obj=plt.subplot(plot_dim,plot_dim,count+1))

                plt.plot([axis_min, axis_max],[axis_min, axis_max], 'k-')
                plt.xlim(0, axis_max)
                plt.ylim(0, axis_max)

                plt.tick_params(axis='both', which='major', labelsize=7)
                plt.subplots_adjust(wspace=0.5, hspace=0.3)

                axins = inset_axes(ax, width="30%", height="30%", loc=4)
                #if str(dataset) == 'EMPclosed' or str(dataset) == 'EMPopen':
                if method == 'zipf':
                    INh2 = import_NSR2_data(data_dir + 'NSR2/' + method+ '_'+ zipfType + '_'+dataset+'_NSR2.txt')
                else:
                    INh2 = import_NSR2_data(data_dir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')
                r2s = ((INh2["R2"]))
                hist_r2 = np.histogram(r2s, range=(0, 1))
                xvals = hist_r2[1] + (hist_r2[1][1] - hist_r2[1][0])
                xvals = xvals[0:len(xvals)-1]
                yvals = hist_r2[0]
                plt.plot(xvals, yvals, 'k-', linewidth=2)
                plt.axis([0, 1, 0, 1.1 * max(yvals)])
                #else:
                #    hist_mete_r2(site, np.log10(obs), np.log10(pred))
                plt.setp(axins, xticks=[], yticks=[])
                count += 1
            #count += (len(datasets) - len(methods))
            #count += 1

        #ax.set_ylabel(-8.5,500,'Observed rank-abundance',rotation='90',fontsize=10)
        fig.subplots_adjust(wspace=0.0001, left=0.1)

        if plot_dim == 3:
            fig.text(0.5, 0.04, 'Predicted rank-abundance', ha='center', va='center')
        elif plot_dim == 4:
            fig.text(0.30, 0.04, 'Predicted rank-abundance', ha='center', va='center')
        else:
            fig.text(0.37, 0.04, 'Predicted rank-abundance', ha='center', va='center')

        fig.text(0.05, 0.5, 'Observed rank-abundance', ha='center', va='center', rotation='vertical')
        #fig.text(0.35, 0.04, 'Predicted rank-abundance', ha='center', va='center')
        #ax.set_xlabel('Observed rank-abundance',fontsize=10)
        if n != 'all':
            fig_name = str(mydir[:-6]) + '/figures/' + 'obs_pred_' + zipfType + '_' + str(n) + '_plots.png'
        else:
            fig_name = str(mydir[:-6]) + '/figures/' + 'obs_pred_' + zipfType + '_plots.png'
        plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
        plt.close()



# Make a function to generate the histogram.
def NSR2_regression(methods, datasets, figname = 'Fig4', Stratified = True, data_dir= mydir):
    fig = plt.figure()
    count  = 0
    if Stratified == True:
        params = ['N','S']
        for i, param in enumerate(params):
            for j, method in enumerate(methods):
                obs_pred_data = import_NSR2_data(data_dir + 'NSR2/Stratified/'+ method +'_NSR2_stratify_t.txt')
                y = np.asarray(list(((obs_pred_data["R2"]))))
                #if i == 0:
                x = np.log10(np.asarray(list(((obs_pred_data[param])))))
                mean_x = np.mean(x)
                mean_y = np.mean(y)
                std_error = sp.stats.sem(y)
                print method, param
                print "mean modified r2 = " + str(mean_y)
                print "modified r2 standard error = " + str(std_error)
                print "mean " + param  + " is " + str(mean_x)
                ax = fig.add_subplot(len(params), len(methods), count+1)
                print len(x), len(y)
                macroecotools.plot_color_by_pt_dens(x, y, 0.1, loglog=0,
                                plot_obj=plt.subplot(len(params), len(methods), count+1))
                slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
                print "r-value is " + str(r_value)
                print "p-value is " + str(p_value)

                plt.xlim(np.amin(x), np.amax(x))
                plt.ylim(-1,1)

                predict_y = intercept + slope * x
                pred_error = y - predict_y
                degrees_of_freedom = len(x) - 2
                residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
                plt.plot(x, predict_y, 'k-')
                plt.axhline(linewidth=2, color='darkgrey',ls='--')

                #plt.hline(0, xmin, xmax, color="0.3", ls='--')
                if i == 0:
                    plt.xlabel(r'$log_{10}(N)$', fontsize = 'large')
                if i == 1:
                    plt.xlabel(r'$log_{10}(S)$', fontsize = 'large')
                plt.subplots_adjust(wspace=0.2, hspace=0.3)
                if i == 0 and j == 0:
                    ax.set_title(r"$\mathbf{Broken-stick}$")
                if i == 0 and j == 1:
                    ax.set_title(r"$\mathbf{Lognormal}$")
                if i == 0 and j == 2:
                    ax.set_title(r"$\mathbf{Log-series}$")
                if i == 0 and j == 3:
                    ax.set_title(r"$\mathbf{Zipf}$")
                ax.tick_params(axis='x', labelsize=8)
                ax.tick_params(axis='y', labelsize=8)

                count += 1

        fig.subplots_adjust(wspace = 0.2, hspace = 0.2, top=0.70)
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
        fig.text(0.001, 0.5, r'$r^{2}$', ha='center', va='center', size = 'large',rotation='vertical')
        fig_name = str(mydir[:-6]) + '/figures/' + figname  + '.png'
        plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        #plt.xscale()
        plt.close()


    else:
        params = ['N','S', 'N/S']
        for i, dataset in enumerate(datasets):
            for k, param in enumerate(params):
                for j, method in enumerate(methods):

                    if (dataset == 'EMPopen' or dataset == 'EMPclosed'):
                        nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')

                    elif dataset == 'HMP':
                        nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')
                    elif dataset == 'MGRAST':
                        nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method+'_MGRAST_NSR2.txt')
                    else:
                        nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method+'_MGRAST'+dataset+'_NSR2.txt')


                    y = ((nsr2_data["R2"]))
                    mean = np.mean(y)
                    std_error = sp.stats.sem(y)
                    print method, param, dataset
                    print "mean modified r2 = " + str(mean)
                    print "modified r2 standard error = " + str(std_error)


                    #print method, dataset, param
                    ax = fig.add_subplot(3, 3, count+1)
                    #ax.set(adjustable='box-forced', aspect='equal')
                    if param == "N" or param == "S":
                        x = np.log10(((nsr2_data[param])))
                    else:
                        N_count = ((nsr2_data["N"]))
                        S_count = ((nsr2_data["S"]))
                        print dataset, method
                        print "mean N is " + str(np.mean(N_count))
                        print "mean S is " + str(np.mean(S_count))
                        x = np.divide(N_count, S_count)
                        x = np.log10(x)

                    macroecotools.plot_color_by_pt_dens(x, y, 0.1, loglog=0,
                                    plot_obj=plt.subplot(3, 3, count+1))
                    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
                    print "r-value is " + str(r_value)
                    print "p-value is " + str(p_value)

                    plt.xlim(np.amin(x), np.amax(x))
                    plt.ylim(-1,1)

                    predict_y = intercept + slope * x
                    pred_error = y - predict_y
                    degrees_of_freedom = len(x) - 2
                    residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
                    plt.plot(x, predict_y, 'k-')
                    plt.axhline(linewidth=2, color='darkgrey',ls='--')

                    #plt.hline(0, xmin, xmax, color="0.3", ls='--')
                    plt.subplots_adjust(wspace=0.2, hspace=0.3)

                    # Plotting
                    if k == 0:
                        plt.xlabel(r'$N_{0}$', fontsize = 12)
                    elif k == 1:
                        plt.xlabel(r'$S_{0}$', fontsize = 12)
                    elif k == 2:
                        plt.xlabel(r'$N_{0}/S_{0}$', fontsize = 12)

                    if k == 1 and j ==0:
                        plt.ylabel(r'$r^{2}_{m}$', fontsize = 'xx-large')

                    #if k == 0 and i ==0 :
                    #    plt.title('HMP', fontsize = 'large')
                    #elif k == 0 and i ==1:
                    #    plt.title('EMP Closed', fontsize = 'large')
                    #elif k == 0 and i ==2:
                    #    plt.title('EMP Open', fontsize = 'large')
                    #elif k == 0 and i ==3:
                    #    plt.title('MG-RAST, 97%', fontsize = 'large')

                    if j == 0 and k == 0:
                        #plt.title(r'\textbf{Broken-stick}', fontsize = 'large')
                        plt.title('Broken-stick', fontsize = 13)
                    elif j == 1 and k == 0:
                        plt.title('METE', fontsize = 13)
                    elif j == 2 and k == 0:
                        plt.title('Zipf', fontsize = 13)
                    #ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=14, verticalalignment='top', bbox=props)
                    leg = plt.legend(loc=1,prop={'size':10})
                    ax.tick_params(axis='x', labelsize=6)
                    ax.tick_params(axis='y', labelsize=6)

                    #tick.label.set_fontsize(14)
                    #leg.draw_frame(False)
                    #plt.legend(loc='upper left')
                    count += 1
        if '97' in datasets:
            title = "MG-RAST 97%"
            fig.suptitle(title, x = 0.54, y = 1.05, fontsize=16)
        elif '95' in datasets:
            title = "MG-RAST 95%"
            fig.suptitle(title, x = 0.54, y = 1.05, fontsize=16)
        elif '99' in datasets:
            title = "MG-RAST 99%"
            fig.suptitle(title, x = 0.54, y = 1.05, fontsize=16)
        elif 'MGRAST' in datasets:
            title = "MG-RAST"
            fig.suptitle(title, x = 0.54, y = 1.05, fontsize=16)
        elif 'HMP' in datasets:
            fig.suptitle("HMP", x = 0.54, y = 1.05, fontsize=16)
        elif 'EMPclosed' in datasets:
            fig.suptitle("EMP (closed)", x = 0.54, y = 1.05, fontsize=16)
        elif 'EMPopen' in datasets:
            fig.suptitle("EMP (open)", x = 0.54, y = 1.05, fontsize=16)

        fig.subplots_adjust(wspace = 0.2, hspace = 0.2, top=0.70)

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)

        #fig.text(0.02, 0.5, r'$r^{2}$', ha='center', va='center', rotation='vertical', size = 'large')
        #st.set_y(0.95)
        #fig.subplots_adjust(top=0.85)
        #ax.set_ylabel('common ylabel')
        #fig.text(-8,-80,'Rank-abundance at the centre of the feasible set',fontsize=10)
        #plt.suptitle(-8.5,500,r'$r^{2}$',rotation='90',fontsize=10)
        fig_name = 'NSR2_GeomMete' + str(dataset) + '.png'
        plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        #plt.xscale()
        plt.close()


def zipf_mle_plots(datasets, data_dir):
        print "what"
        fig = plt.figure()

        count  = 0
        plot_dim = len(datasets)
        for i, dataset in enumerate(datasets):
            if  (dataset == 'EMPopen' or dataset == 'EMPclosed'):
                zipf_data = import_NSR2_data(data_dir + 'NSR2/zipf_'+ dataset+'_NSR2.txt')
            elif dataset == 'HMP':
                zipf_data = import_NSR2_data(data_dir + 'NSR2/zipf_' + dataset+'_NSR2.txt')
            elif dataset == 'MGRAST':
                zipf_data = import_NSR2_data(data_dir + 'NSR2/zipf_' + dataset+'_NSR2.txt')
            else:
                zipf_data = import_NSR2_data(data_dir + 'NSR2/zipf_MGRAST' + dataset+'_NSR2.txt')

            N =  np.log10(((zipf_data["N"])))
            y = ((zipf_data["gamma"])) * -1
            x = np.log10(((zipf_data["NmaxObs"])))
            ax = fig.add_subplot(plot_dim, plot_dim, count+1)
            mean = np.mean(x)
            std_error = sp.stats.sem(x)

            print dataset
            print "mean gamma = " + str(mean)
            print "gamma standard error = " + str(std_error)

            # Plot 1
            macroecotools.plot_color_by_pt_dens(N, y, 0.1, loglog=0,
                        plot_obj=plt.subplot(plot_dim, plot_dim, count+1))
            slope, intercept, r_value, p_value, std_err = stats.linregress(N,y)
            plt.xlim(np.amin(N), np.amax(N))
            plt.ylim(np.amin(y),0)
            plt.xticks(fontsize = 6) # work on current fig
            plt.yticks(fontsize = 6)
            if dataset == 'HMP':
                ax.set_ylabel(r'HMP' '\n' r'$\alpha$', rotation=90, size=12)
            elif dataset == 'EMPopen':
                ax.set_ylabel(r'EMP (open)' '\n' r'$\alpha$', rotation=90, size=12)
            elif dataset == 'EMPclosed':
                ax.set_ylabel(r'EMP (closed)' '\n' r'$\alpha$', rotation=90, size=12)
            elif dataset == '95':
                ax.set_ylabel(r'MG-RAST 95%' '\n' r'$\alpha$', rotation=90, size=12)
            elif dataset == '97':
                ax.set_ylabel(r'MG-RAST 97%' '\n' r'$\alpha$', rotation=90, size=12)
            elif dataset == '99':
                ax.set_ylabel(r'MG-RAST 99%' '\n' r'$\alpha$', rotation=90, size=12)
            elif dataset == 'MGRAST':
                ax.set_ylabel(r'MG-RAST' '\n' r'$\alpha$', rotation=90, size=12)

            if i == plot_dim-1:
                plt.xlabel(r'$N_{0}$', fontsize = 12)
            predict_y = intercept + slope * N
            pred_error = y - predict_y
            degrees_of_freedom = len(x) - 2
            residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)
            plt.plot(N, predict_y, 'k-')
            plt.axhline(linewidth=2, color='darkgrey',ls='--')
            plt.subplots_adjust(wspace=0.2, hspace=0.3)

            # Plot 2
            macroecotools.plot_color_by_pt_dens(x, y, 0.1, loglog=0,
                        plot_obj=plt.subplot(plot_dim, plot_dim, count+2))

            slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
            plt.xlim(np.amin(x), np.amax(x))
            plt.ylim(np.amin(y),0)
            plt.xticks(fontsize = 6) # work on current fig
            plt.yticks(fontsize = 6)
            predict_y = intercept + slope * x
            pred_error = y - predict_y
            degrees_of_freedom = len(x) - 2
            #residual_std_error = np.sqrt(np.sum(pred_error**2) / degrees_of_freedom)

            plt.plot(x, predict_y, 'k-')
            plt.axhline(linewidth=2, color='darkgrey',ls='--')

            if i == plot_dim-1:
                plt.xlabel(r"$N_{max}$", fontsize = 12)
            #plt.hline(0, xmin, xmax, color="0.3", ls='--')
            plt.ylabel(r'$\alpha$', fontsize = 12)

            plt.subplots_adjust(wspace=0.2, hspace=0.3)
            # Plotting
            # Plot 3
            density = gaussian_kde(y)
            ys = np.linspace(-6,0,200)
            density.covariance_factor = lambda : .25
            density._compute_covariance()
            plt.subplot(plot_dim, plot_dim, count+3)
            plt.plot(ys,density(ys))

            plt.xticks(fontsize = 6) # work on current fig
            plt.yticks(fontsize = 6)
            if i == plot_dim-1:
                plt.xlabel(r'$\alpha$', fontsize = 12)
            plt.ylabel('Probability Density', fontsize = 7)
            #plt.legend(loc='best')
            count += plot_dim

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

        fig.subplots_adjust(wspace = 0.4, hspace = 0.35, top=0.85)
        fig_name = str(mydir[:-6]) + '/figures/' + '_'.join(datasets) + '_zipf_mle_plots.png'
        #ax.set_aspect('equal')
        plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        #plt.xscale()
        plt.close()

def get_slope_row(row):
    #slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
    #return row['a'] % row['c']
    x_vals = [row['N0_05'], row['N0_025'], \
    row['N0_0125'], row['N0_00625'], row['N0_003125'], row['N0_0015625']]
    x_vals = [math.log(x, 10) for x in x_vals]
    y_vals = [row['r2_05'], row['r2_025'], \
    row['r2_0125'], row['r2_00625'], row['r2_003125'], row['r2_0015625']]
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_vals,y_vals)
    return slope*-1

def CV_KDE(x):
    # remove +/- inf
    x = x[np.logical_not(np.isnan(x))]
    grid = GridSearchCV(KernelDensity(),
                    {'bandwidth': np.logspace(0.1, 1.0, 30)},
                    cv=20) # 20-fold cross-validation
    grid.fit(x[:, None])
    x_grid = np.linspace(np.min(x), np.max(x), len(x))
    kde = grid.best_estimator_
    print "bandwidth is " + str(kde.bandwidth)
    pdf = np.exp(kde.score_samples(x_grid[:, None]))
    # returns grod for x-axis,  pdf, and bandwidth
    return_tuple = (x_grid, pdf, kde.bandwidth)
    return return_tuple

def get_kdens_choose_kernel(xlist,expand, kernel = 0.5):
    """ Finds the kernel density function across a vector of values """
    xlist = xlist[np.logical_not(np.isnan(xlist))]
    density = gaussian_kde(xlist)
    n = len(xlist)
    if expand == False:
        xs = np.linspace(min(xlist),max(xlist),n)
    else:
        xs = np.linspace(min(xlist - expand),max(xlist + expand),n)
    #xs = np.linspace(0.0,1.0,n)
    density.covariance_factor = lambda : kernel
    density._compute_covariance()
    D = [xs,density(xs)]
    return D

def plot_subsampled_data(methods, datasets, data_dir= mydir):
    for i, dataset in enumerate(datasets):
        fig = plt.figure()
        count  = 0
        plot_dim = len(methods)
        for j, method in enumerate(methods):
            input_file = data_dir + \
            'SubSampled-Data/' + dataset+ '_' + method + '_SubSampled_Data.txt'
            subsam_data = import_subsampled_data_pandas(input_file)
            subsam_data['slope'] = subsam_data.apply(lambda row: get_slope_row(row), axis=1)
            # Get scatter plot of subsampled data
            Ns = ['N0_05', 'N0_025', 'N0_0125', 'N0_00625', 'N0_003125', 'N0_0015625']
            r2s = ['r2_05', 'r2_025', 'r2_0125', 'r2_00625', 'r2_003125', 'r2_0015625']
            Ns_r2s = zip(Ns, r2s)
            all_N_tuples = []
            all_r2_tuples = []
            for k in Ns_r2s:
                all_N_tuples.append(subsam_data[k[0]])
                all_r2_tuples.append(subsam_data[k[1]])
            all_N_tuples_flat = [item for sublist in all_N_tuples \
                for item in sublist]
            all_r2_tuples_flat = [item for sublist in all_r2_tuples \
                for item in sublist]
            x = np.asarray(np.log10(all_N_tuples_flat))
            y = np.asarray(all_r2_tuples_flat)
            for l in range(2):
                if l == 0:
                    ax = fig.add_subplot(2, plot_dim, count + 1)
                    plt.scatter(x, y, alpha=0.05)
                    slope, intercept, r_value, p_value, std_err = \
                        stats.linregress(x,y)
                    print "r-value is " + str(r_value)
                    print "p-value is " + str(p_value)
                    plt.xlim(np.amin(x), np.amax(x))
                    plt.ylim(np.amin(y), np.amax(y))
                    predict_y = intercept + slope * x
                    pred_error = y - predict_y
                    plt.plot(x, predict_y, 'k-')
                    degrees_of_freedom = len(x) - 2
                    residual_std_error = np.sqrt(np.sum(pred_error**2) / \
                        degrees_of_freedom)
                    plt.xticks(fontsize = 6) # work on current fig
                    plt.yticks(fontsize = 6)
                    plt.ylabel(r'$\overline{r_m^2}$', fontsize = 10, rotation = 0)
                    plt.xlabel(r'$\overline{N_0}$', fontsize = 10)

                elif l == 1:
                    ax = fig.add_subplot(2, plot_dim, count + 4)
                    # Get the KDE
                    slope_data = subsam_data['slope']
                    kde_data = get_kdens_choose_kernel(slope_data, 0.4, kernel = 0.4)
                    ax.plot(kde_data[0], kde_data[1])
                    plt.xlim(np.amin(x), np.amax(x))
                    plt.xticks(fontsize = 6) # work on current fig
                    plt.yticks(fontsize = 6)
                    plt.ylabel('Probability Density', fontsize = 7)
                    plt.xlabel(r'$m$', fontsize = 10)

            count += 1
        if dataset == 'HMP':
            name = str(mydir[:-6]) + '/figures' + '/Fig3.png'
            print name
        elif dataset == 'EMPclosed':
            name = str(mydir[:-6]) + '/figures' + '/Fig4.png'
        elif dataset == 'MGRAST':
            name = str(mydir[:-6]) + '/figures' + '/Fig5.png'
        plt.tight_layout()
        plt.savefig(name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        #plt.xscale()
        plt.close()


def plot_color_by_pt_dens1(x, y, radius, loglog=0, plot_obj=None):
    """Plot bivariate relationships with large n using color for point density

    Inputs:
    x & y -- variables to be plotted
    radius -- the linear distance within which to count points as neighbors
    loglog -- a flag to indicate the use of a loglog plot (loglog = 1)

    The color of each point in the plot is determined by the logarithm (base 10)
    of the number of points that occur with a given radius of the focal point,
    with hotter colors indicating more points. The number of neighboring points
    is determined in linear space regardless of whether a loglog plot is
    presented.
    """
    plot_data = macroecotools.count_pts_within_radius(x, y, radius, loglog)
    sorted_plot_data = np.array(sorted(plot_data, key=lambda point: point[2]))

    if plot_obj == None:
        plot_obj = plt.axes()

    if loglog == 1:
        plot_obj.set_xscale('log')
        #plot_obj.set_yscale('log')
        plot_obj.scatter(sorted_plot_data[:, 0], sorted_plot_data[:, 1],
                         c = np.sqrt(sorted_plot_data[:, 2]), edgecolors='none')
        plot_obj.set_xlim(min(x) * 0.5, max(x) * 2)
        plot_obj.set_ylim(min(y) * 0.5, max(y) * 2)
    else:
        plot_obj.scatter(sorted_plot_data[:, 0], sorted_plot_data[:, 1],
                    c = log10(sorted_plot_data[:, 2]), edgecolors='none')
    return plot_obj


def obs_pred_Nmax_plot(methods,datasets, figname = 'Fig5', zipfType = 'glm', radius=2, stratify = True, data_dir= mydir):
    """Multiple obs-predicted plotter"""
    fig = plt.figure()
    count = 0
    if stratify == True:
        params = ['Nmax', 'evenness', 'skewness']
        #### Add evenness and skewness
        for i, param in enumerate(params):
            for j, method in enumerate(methods):
                obs_pred_data = import_NSR2_data(data_dir + 'NSR2/Stratified/'+ method +'_NSR2_stratify.txt')
                N = np.asarray(list(((obs_pred_data["N"]))))
                if i == 0:
                    pred = list(((obs_pred_data["NmaxPred"])))
                elif i == 1:
                    pred = list(((obs_pred_data["evennessPred"])))
                elif i == 2:
                    pred = list(((obs_pred_data["skewnessPred"])))

                y = np.asarray(pred)

                ax = fig.add_subplot(3, len(methods), count+1)

                x_axis_min = min(N)
                x_axis_max = max(N)
                y_axis_min = min(y)
                y_axis_max = max(y)

                plt.subplots_adjust(wspace=1, hspace=1)

                if i == 2 and (j == 3 or j ==1):
                    merge_filter = [k for k in zip(N, y) if k[1] > 0.000001]
                    N = np.asarray([l[0] for l in merge_filter])
                    y = np.asarray([m[1] for m in merge_filter])

                print param, method

                if i ==0 and j == 0:
                    ax.set_title(r"$\mathbf{Broken-stick}$")
                    plt.ylabel(r'$log_{10}(N_{max})$', fontsize = 'large')
                elif i ==0 and j == 1:
                    ax.set_title(r"$\mathbf{Lognormal}$")
                elif i ==0 and j == 2:
                    ax.set_title(r"$\mathbf{Log-series}$")
                elif i ==0 and j ==3:
                    ax.set_title(r"$\mathbf{Zipf}$")
                elif i == 1 and j == 0:
                    plt.ylabel(r'$Evenness, \, log_{10}$', fontsize = 'large')
                elif i == 2 and j == 0:
                    plt.ylabel(r'$Skewness, \, log_{10}$', fontsize = 'large')

                macroecotools.plot_color_by_pt_dens(N, y, radius, loglog=1,
                                plot_obj=plt.subplot(3,len(methods),count+1))
                slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(N),np.log10(y))
                print "r-value is " + str(r_value)
                print "p-value is " + str(p_value)
                print "slope is " + str(slope)
                if i == 0:
                    plt.plot([0, x_axis_max],[0, x_axis_max], 'k-')
                if (i ==1  or i ==2) and j ==0:
                    plt.xlim(0, x_axis_max)
                    plt.ylim(0, 1)
                elif i == 2 and (j == 2):
                    plt.xlim(0, x_axis_max)
                    plt.ylim(0, 1)
                else:
                    plt.xlim(0, x_axis_max)
                    plt.ylim(0, y_axis_max)


                plt.yscale('log')
                plt.tick_params(axis='both', which='major', labelsize=7)
                plt.subplots_adjust(wspace=0.5, hspace=0.3)
                count += 1

        plt.tight_layout(pad=2, w_pad=2, h_pad=2)
        #fig.text(0.02, 0.5, r'$N_{max}$', ha='center', va='center', size = 'x-large', rotation='vertical')
        fig.text(0.527, 0.05, r'$log_{10}(N)$', ha='center', va='center', size = 'large',rotation='horizontal')
        fig_name = str(mydir[:-6]) + '/figures/' + figname + '.png'
        plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
        plt.close()

    else:
        plot_dim = len(datasets)
        for i, dataset in enumerate(datasets):
            for j, method in enumerate(methods):
                if method == 'zipf':
                    obs_pred_data = import_NSR2_data(data_dir + 'NSR2/' + \
                        method + '_' + zipfType + '_'+ dataset+'_NSR2.txt')
                else:
                    obs_pred_data = import_NSR2_data(data_dir + 'NSR2/' + \
                        method + '_'+ dataset+'_NSR2.txt')
                obs = list(((obs_pred_data["NmaxObs"])))
                pred = list(((obs_pred_data["NmaxPred"])))
                obs = np.asarray(obs)
                pred = np.asarray(pred)
                ax = fig.add_subplot(plot_dim, plot_dim, count+1)
                N = np.asarray(list(((obs_pred_data["N"]))))

                y = pred

                x_axis_min = min(N)
                x_axis_max = max(N)
                y_axis_min = min(y)
                y_axis_max = max(y)

                plt.subplots_adjust(wspace=1, hspace=1)

                if j == 0:
                    if all(x in datasets for x in ['HMP','EMPclosed','EMPopen', 'MGRAST']) == True:
                        if dataset == 'HMP':
                            ax.set_ylabel("HMP", rotation=90, size=12)
                        elif dataset == 'EMPclosed':
                            ax.set_ylabel("EMP (closed)", rotation=90, size=12)
                        elif  dataset == 'EMPopen':
                            ax.set_ylabel("EMP (open)", rotation=90, size=12)
                        elif dataset == 'MGRAST':
                            ax.set_ylabel("MG-RAST", rotation=90, size=12)
                    if all(x in datasets for x in ['95', '97', '99']) == True:
                        if dataset == '95':
                            ax.set_ylabel("MG-RAST 95%", rotation=90, size=12)
                        elif dataset == '97':
                            ax.set_ylabel("MG-RAST 97%", rotation=90, size=12)
                        elif dataset == '99':
                            ax.set_ylabel("MG-RAST 99%", rotation=90, size=12)
                    else:
                        if dataset == 'HMP':
                            ax.set_ylabel("HMP", rotation=90, size=12)
                        elif dataset == 'EMPclosed' or method == 'EMPopen':
                            ax.set_ylabel("EMP", rotation=90, size=12)
                        elif dataset == '97':
                            ax.set_ylabel("MG-RAST 97%", rotation=90, size=12)
                        elif dataset == 'MGRAST':
                            ax.set_ylabel("MG-RAST", rotation=90, size=12)


                if i == 0 and j == 0:
                    ax.set_title("Broken-stick")
                elif i == 0 and j == 1:
                    ax.set_title("METE")
                elif i == 0 and j == 2:
                    ax.set_title("Zipf")

                macroecotools.plot_color_by_pt_dens(N, y, radius, loglog=1,
                                plot_obj=plt.subplot(plot_dim,plot_dim,count+1))
                slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(N),np.log10(y))
                print "r-value is " + str(r_value)
                print "p-value is " + str(p_value)
                print "slope is " + str(slope)

                plt.plot([x_axis_min, x_axis_max],[x_axis_min, x_axis_max], 'k-')
                plt.xlim(0, x_axis_max)
                plt.ylim(0, y_axis_max)

                plt.tick_params(axis='both', which='major', labelsize=7)
                plt.subplots_adjust(wspace=0.5, hspace=0.3)



                count += 1

        #ax.set_ylabel(-8.5,500,'Observed rank-abundance',rotation='90',fontsize=10)
        fig.subplots_adjust(wspace = 0.2, hspace = 0.2, top=0.70)

        if plot_dim == 3:
            fig.text(0.5, 0.04, r'$N$', ha='center', va='center')
        elif plot_dim == 4:
            fig.text(0.30, 0.04, r'$N$', ha='center', va='center')
        else:
            fig.text(0.37, 0.04, r'$N$', ha='center', va='center')

        fig.text(0.02, 0.5, r'$N_{max}$', ha='center', va='center', rotation='vertical')
        #fig.ylabel('ylabel')
        #fig.text(0.35, 0.04, 'Predicted rank-abundance', ha='center', va='center')
        #ax.set_xlabel('Observed rank-abundance',fontsize=10)
        plt.tight_layout(pad=2, w_pad=2, h_pad=2)
        fig_name = str(mydir[:-6]) + '/figures/' + '_'.join(datasets) + '_obs_pred_Nmax_error_' \
            + zipfType + '_plots.png'
        plt.savefig(fig_name, dpi=600)#, bbox_inches = 'tight')#, pad_inches=0)
        plt.close()


def stratifyData(methods,datasets, totalSADs = 1000, zipfType = 'glm', data_dir= mydir):
    # Number of lines in each file
    MGRAST_sites = 1174
    HMP_sites = 4504
    EMPclosed_sites = 14979
    Total = MGRAST_sites + HMP_sites + EMPclosed_sites
    for i, method in enumerate(methods):
        OUT1 = open(data_dir + 'ObsPred/Stratified/'+ method +'_obs_pred_stratify_t.txt', 'w')
        OUT2 = open(data_dir + 'NSR2/Stratified/'+ method +'_NSR2_stratify_t.txt', 'w')
        for j, dataset in enumerate(datasets):
            print "what"
            n = totalSADs
            #if dataset == 'HMP':
            #    n = round((HMP_sites / Total) * totalSADs)
            #elif dataset == 'EMPclosed':
            #    n = round((EMPclosed_sites / Total) * totalSADs)
            #elif dataset == 'MGRAST':
            #    n = round((MGRAST_sites / Total) * totalSADs)

            if method == 'zipf':
                if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                    ( str(dataset) == 'EMPopen'):
                    obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_'+ zipfType +'_'+dataset+'_obs_pred.txt')
                    nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method + '_'+ zipfType +'_'+ dataset+'_NSR2.txt')
                elif str(dataset) == 'MGRAST':
                    obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_'+ zipfType + '_' + 'MGRAST_obs_pred.txt')
                    nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method+ '_'+ zipfType  + '_MGRAST_NSR2.txt')
                else:
                    obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_'+ zipfType + '_' + 'MGRAST' + dataset +'_obs_pred.txt')
                    nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method+ '_'+ zipfType  + '_MGRAST'+dataset+'_NSR2.txt')
            else:
                if ( str(dataset) == 'EMPclosed') or ( str(dataset) == 'HMP') or \
                    ( str(dataset) == 'EMPopen'):
                    obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method+'_'+dataset+'_obs_pred.txt')
                    nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method+'_'+dataset+'_NSR2.txt')
                elif str(dataset) == 'MGRAST':
                    obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_' + 'MGRAST_obs_pred.txt')
                    nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method+'_MGRAST_NSR2.txt')
                else:
                    obs_pred_data = import_obs_pred_data(data_dir + 'ObsPred/' + method + '_' + 'MGRAST' + dataset +'_obs_pred.txt')
                    nsr2_data = import_NSR2_data(data_dir + 'NSR2/' + method+'_MGRAST'+dataset+'_NSR2.txt')
            print "whatt"
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
            obs2 = []
            pred2 = []
            site2 = []
            uniqueSites = np.unique(site)
            print len(uniqueSites)
            print len(N)
            randomSites = np.random.choice(range(len(uniqueSites)), size=n, replace=False)
            count = 0
            for p, q in enumerate(uniqueSites):
                if q in randomSites:

                    print>> OUT2, count, N[p], S[p], NmaxObs[p], NmaxPred[p], \
                        evennessObs[p], evennessPred[p], skewnessObs[p], skewnessPred[p], R2[p]
                    count += 1
            for r, s in enumerate(site):
                if s in randomSites:
                    obs2.append(obs[r])
                    pred2.append(pred[r])
                    site2.append(site[r])
            print method, dataset


            #if n == 'all' or len(obs) <= n:
            #    obs2 = list(obs)
            #    pred2 = list(pred)
            #    site2 = list(site)

            #else:
            #    if len(obs) > n:
            #        #inds = np.random.choice(range(len(obs)), size=n, replace=False)
            #        for randomSite in randomSites:
            #            obs2.append(obs[randomSite])
            #            pred2.append(pred[randomSite])
            #            site2.append(site[randomSite])
            #print len(obs2), len(pred2)
            obs = np.asarray(obs2)
            pred = np.asarray(pred2)

            #site = np.asarray(site2)
            k_minus1 = obs2[0]
            count_sites = 0
            for k, sp in enumerate(obs2):
                print>> OUT1, k, obs2[k], pred2[k]
            #site =  np.asarray(site2)
        OUT1.close()
        OUT2.close()

def plot_example_fig(figname = 'Fig1', data_dir= mydir):
    SAD = [10000, 8000, 6000, 5000, 1000, 200, 100,  20, 18, 16, 14, 12, 10, 4,5,
        4, 4, 3, 3, 2, 2, 2, 2, 2,2, 1, 1, 1, 1, 1,1,1,1, 1, 1, 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
    SAD.sort()
    SAD.reverse()
    x = range(1, len(SAD) +1)
    N = sum(SAD)
    S = len(SAD)
    geom = get_GeomSeries(N, S, False)

    logSeries = mete.get_mete_rad(S, N)[0]

    lognorm_pred = lognorm(SAD)
    lognorm_SAD = lognorm_pred.lognorm_glm()
    lognorm_SAD = np.ceil(lognorm_SAD)
    lognorm_SAD.astype(int)

    zipf_pred = zipf(SAD, 'fmin')
    zipf_glm = zipf_pred.from_glm()
    zipf_SAD = np.ceil(zipf_glm)
    zipf_SAD.astype(int)

    fig = plt.figure()
    plt.plot()
    # try with mle zipf
    zipf_class = zipf(SAD, 'fmin')
    pred_tuple = zipf_class.from_cdf()
    pred = pred_tuple[0]

    plt.plot(x, SAD,color = 'r', linestyle = '-', linewidth=2, label="Emperical distribution")
    plt.plot(x, geom,color = 'b', linestyle = '-', linewidth=2, label="Broken-stick")
    plt.plot(x, lognorm_SAD, color = 'b',linestyle = '--', linewidth=2, label="Lognormal")
    plt.plot(x, logSeries, color = 'b',linestyle = '-.', linewidth=2, label="Log-series")
    plt.plot(x, pred, color = 'b',linestyle = ':',linewidth=2,  label="Zipf")

    plt.tight_layout()
    plt.xlabel(r'$Rank \; Abundance$', fontsize = 18)
    plt.ylabel(r'$Abundance, \, log_{10}$', fontsize = 18)
    output = "dorm_fix_prob.png"
    plt.legend(loc='upper right')
    plt.yscale('log')
    plt.xlim(1, len(SAD))
    plt.ylim(1, max(zipf_SAD))
    #fig.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    fig_name = str(mydir[:-6]) + '/figures/' + figname + '.png'
    plt.savefig(fig_name, bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    #plt.xscale()
    plt.close()

#datasets = ['EMPclosed', 'EMPopen', 'MGRAST', 'HMP', '97', '95', '99']
#datasets = ['HMP','EMPclosed','MGRAST']
datasets = ['EMPclosed']
#methods = ['zipf']
#methods = ['geom', 'lognorm', 'mete']
methods = ['geom', 'lognorm', 'mete', 'zipf']
#plot_subsampled_data(methods, datasets)

#stratifyData(methods,datasets, 352899)
#obs_pred_Nmax_plot(methods, datasets, stratify = True, zipfType = 'glm')

#plot_obs_pred_sad(methods, datasets, 352899, zipfType = 'glm')
#plot_obs_pred_sad(methods, datasets, 352899, zipfType = 'mle')
#plot_obs_pred_sad(methods, datasets, 100, zipfType = 'mle')
#stratifyData(methods,datasets, totalSADs = 100)
#NSR2_regression(methods, datasets)
#obs_pred_Nmax_plot(methods, datasets, stratify = True, zipfType = 'glm')
plot_obs_pred_sad(methods, datasets, n = 'all')
#plot_example_fig()
