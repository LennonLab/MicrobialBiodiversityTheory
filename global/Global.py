# -*- coding: utf-8 -*-
from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import random
import time
import sys

sys.path.append("/Users/lisalocey/Documents/repos/RareBio/global/GenAnalysis/tools/")
import macroeco_distributions as md
import mete
import pln
import datetime


""" There are upwards 1030 organisms on Earth belonging to upwards 108 species.
The vast majority of both species and individuals belong to Bacteria and Archaea.
The cyanobacteria Prochlorococcus and Synechococcus are among the most abundant
microbes, with abundances estimated at 2.9 ± 0.1 × 1027 and 7.0 ± 0.3 × 1026 cells,
respectively. If the abundance of the most abundant microbe is within this range,
then its relative abundance should range from 0.02 to 0.33% of global microbial
abundance.
"""

def truncate(RAD, maxn):

    """ Enforce a maximum abundance """

    while max(RAD) > maxn:
        m = max(RAD)
        i = RAD.index(m)
        n1 = RAD[i]
        n2 = random.randint(1,n1-1)
        RAD[i] -= n2

        i = random.randint(0, len(RAD)-1)
        RAD[i] += n2

    return RAD


def get_GeomSeries(N,S,zeros):

    """ MLE for the geometric series based on N and S """

    rank = range(1,S+1)
    cdf = [(S-i+0.5)/S for i in rank]
    SNratio = S/N
    if zeros == False:
        abd = md.trunc_geom.ppf(np.array(cdf), SNratio, N)

    return abd



def RandomComposition_xx(q, n):
    indices = sorted(np.random.randint(1, q, n-1))
    parts = [(indices + [q])[i] - ([0] + indices)[i] for i in range(len(indices)+1)]
    return parts



def RandCompFast(Q, N):

    """ Random Integer Composition based on N and S """

    composition = []

    indices = []
    while len(indices) < N-1:
        index = random.randint(1, Q-1)
        if index in indices: continue
        else: indices.append(index)

    indices.sort()
    indices.append(Q)

    nsum = 0
    composition = []
    for i in indices:
        i -= sum(composition)
        composition.append(i)
        nsum += i

    composition.sort()
    composition.reverse()

    return composition




def RandBin(N, S, maxn = 0, trunc=False):

    """ Random binomial based on N and S """

    RAD = [0]*S
    while N >= 1:
        i = random.randint(0, S-1)
        RAD[i] += 1
        N -= 1

    if trunc is True:
        RAD = truncate(RAD, maxn)

    RAD.sort()
    RAD.reverse()

    return RAD




def DomDecay(N, mu):

    """ Random fraction Dominance Decay, based on N and mu """

    n1 = random.randint(1, N-1)
    RAD = [n1, N - n1]

    BP  = 1.0
    while BP > mu:

        n1 = max(RAD)
        RAD.remove(n1)

        n2 = random.randint(1, n1-1)
        RAD += [n2, n1 - n2]
        BP = max(RAD)/sum(RAD)

    RAD.sort(reverse=True)
    #RAD.reverse()

    return RAD



def RandFrac(N, mu):

    """ Random fraction based on N and S """

    RAD = [N]
    BP  = 1.0
    while BP > mu:

        i = random.randint(0, len(RAD)-1)
        n1 = RAD[i]
        n2 = random.randint(1, n1-1)
        RAD[i] -= n2
        RAD.append(n2)
        BP = max(RAD)/sum(RAD)

    RAD.sort()
    RAD.reverse()

    return RAD



def Zip_Frac(N, mu):

    """ Random fraction based on N and S """

    RAD = [N]
    while min(RAD) >= 1.0:

        n1 = RAD[-1]
        RAD[-1] = mu * n1
        n1 = n1 - (mu * n1)
        RAD.append(n1)

    RAD.sort()
    RAD.reverse()

    return RAD



def plot_RADs_canonical(N, S):

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    # Predicted geometric series
    print 'generating geometric series'
    t0 = time.time()
    predRAD = get_GeomSeries(N, S, False) # False mean no zeros allowed
    t = time.time() - t0
    print 'time for geometric series:',t
    ranks = range(1, S+1)
    plt.plot(ranks, predRAD, lw = 1, c='m')

    # Predicted log-series
    print 'generating log-series'
    t0 = time.time()
    logSeries = mete.get_mete_rad(S, N)
    t = time.time() - t0
    print 'time for log-series:',t
    predRAD = logSeries[0]
    ranks = range(1, S+1)
    plt.plot(ranks, predRAD, lw = 1, c='c')

    # Predicted PLN
    print 'generating Poisson log-normal'
    t0 = time.time()
    predRAD = pln.get_rad_from_obs(predRAD, 'pln')
    t = time.time() - t0
    print 'time for log-normal:',t
    ranks = range(1, len(predRAD)+1)
    plt.plot(ranks, predRAD, lw = 1, c='gray')

    plt.yscale('log')
    plt.savefig('/Users/lisalocey/Desktop/RareBio/figs/GlobalRADs_N='+str(int(N))+'_S='+str(int(S))+'.png',dpi=600)

    plt.show()

    return



def plot_fracRADs(N, mu):

    fig = plt.figure()
    fig.add_subplot(1,1,1)

    plt.title('Total abundance = '+str('%.1e' % N)+', BP = '+str('%.1e' % mu))

    """
    # Predicted Zipf random fraction
    print 'generating Zipf fraction'
    t0 = time.time()
    predRAD = Zip_Frac(N, mu)
    t = time.time() - t0
    print 'time for Zipf fraction:',t
    BP = round(100*max(predRAD)/sum(predRAD), 3)
    print 'Berger Parker as % = '+str(BP)+'  S = '+str(len(predRAD))+'\n'
    ranks = range(1, len(predRAD)+1)
    plt.plot(ranks, predRAD, lw = 1, c='c', label='Zipf fraction; S = '+str(len(predRAD))+'%')


    # Predicted random fraction
    print 'generating random fraction'
    t0 = time.time()
    predRAD = RandFrac(N, mu)
    t = time.time() - t0
    print 'time for Random fraction:',t
    BP = round(100*max(predRAD)/sum(predRAD), 3)
    print 'Berger Parker as % = '+str(BP)+'  S = '+str(len(predRAD))+'\n'
    ranks = range(1, len(predRAD)+1)
    plt.plot(ranks, predRAD, lw = 1, c='m', label='Random fraction; S = '+str(len(predRAD))+'%')
    """

    # Predicted Dominance Decay
    print 'generating DomDecay'
    t0 = time.time()
    predRAD = DomDecay(N, mu)
    t = time.time() - t0
    print 'time for Dominance Decay:',t
    S = len(predRAD)
    BP = round(100 * max(predRAD)/sum(predRAD), 3)
    print 'Berger Parker as % = '+str(BP)+'  S = '+str(S)+'\n'
    ranks = range(1, len(predRAD)+1)
    plt.scatter(ranks, predRAD, lw = 1, c='m', linewidth=0.0, label='Dominance Decay; S = '+str(S))

    # Predicted geometric series, Ken
    print 'generating geometric series'
    t0 = time.time()
    predRAD = RandCompFast(N, S)
    t = time.time() - t0
    print 'time for Ken\'s geometric series:',t
    BP = round(100*max(predRAD)/sum(predRAD), 3)
    print 'Berger Parker as % = '+str(BP)+'\n'
    ranks = range(1, len(predRAD)+1)
    plt.scatter(ranks, predRAD, lw = 1, c='c', linewidth=0.0, label='Geometric Series for N & S; BP = '+str(BP)+'%')


    leg = plt.legend(loc=3,prop={'size':16})
    leg.draw_frame(False)

    plt.xlabel('Rank in abundance', fontsize = 16)
    plt.ylabel('Abundance', fontsize = 16)
    plt.yscale('log')
    #plt.xscale('log')

    plt.savefig('/Users/lisalocey/Desktop/RareBio/global/figs/Global_N='+str(int(N))+'_S='+str(len(predRAD))+'.png',dpi=600)
    plt.close()

    return



def GeomSeries(N, S):

    fig = plt.figure()
    fig.add_subplot(1,1,1)

    plt.title('Total abundance = '+str('%.1e' % N)+', Richness = '+str('%.1e' % S))

    # Predicted geometric series, Ken
    print 'generating geometric series'
    t0 = time.time()
    predRAD = RandCompFast(N, S)
    t = time.time() - t0
    print 'time for Ken\'s geometric series:',t
    BP = round(100*max(predRAD)/sum(predRAD), 3)
    print 'Berger Parker as % = '+str(BP)+'\n'
    ranks = range(1, len(predRAD)+1)
    plt.plot(ranks, predRAD, lw = 1, c='m', label='Ken\'s Integer Compositions')

    """
    # Predicted geometric series, Xiao Xiao
    print 'generating geometric series'
    t0 = time.time()
    predRAD = RandomComposition_xx(N, S)
    predRAD.sort()
    predRAD.reverse()
    t = time.time() - t0
    print 'time for Xiao\'s geometric series:',t
    BP = 100*max(predRAD)/sum(predRAD)
    print 'Berger Parker as % = '+str(BP)+'\n'
    ranks = range(1, len(predRAD)+1)
    plt.yscale('log')
    plt.plot(ranks, predRAD, lw = 1, c='c', label='Xiao\'s Geom. Series; BP = '+str(BP)+'%')
    """

    # Predicted geometric series, Nathan
    print 'generating geometric series'
    t0 = time.time()
    predRAD = RandomComposition_xx(N, S)
    predRAD.sort()
    predRAD.reverse()
    t = time.time() - t0
    print 'time for Xiao\'s geometric series:',t
    BP = 100*max(predRAD)/sum(predRAD)
    print 'Berger Parker as % = '+str(BP)+'\n'
    ranks = range(1, len(predRAD)+1)
    plt.yscale('log')
    plt.plot(ranks, predRAD, lw = 1, c='c', label='Nathan\'s Broken Stick')

    # Predicted geometric series, from METE
    print 'generating geometric distribution'
    t0 = time.time()
    predRAD = get_GeomSeries(N, S, False) # False mean no zeros allowed
    t = time.time() - t0
    print 'time for geometric distribution:',t
    ranks = range(1, S+1)
    plt.plot(ranks, predRAD, lw = 1, c='gray', label='MaxEnt, Geometric Distribution')

    leg = plt.legend(loc=3,prop={'size':16})
    leg.draw_frame(False)

    plt.xlabel('Rank in abundance', fontsize = 16)
    plt.ylabel('Abundance', fontsize = 16)
    plt.yscale('log')
    #plt.show()
    plt.savefig('/Users/lisalocey/Desktop/GeomDist_N='+str(int(N))+'_S='+str(int(S))+'.png',dpi=600)
    plt.close()

    return


N = 10**15
mu = 0.0005
S = 10**3

#plot_comps(N, S)

t1 = time.time()
GeomSeries(N, S)
t2 = time.time()
print t2-t1,"seconds to complete Geometric Series:"

#t1 = time.time()
#plot_fracRADs(N, mu)
#t2 = time.time()
#print t2-t1,"seconds to complete:"
