from __future__ import division
import sys
import random
import matplotlib.pyplot as plt
import os

mydir = os.path.expanduser("~/Desktop")
sys.path.append(mydir + "/Repos/RareBio/global/GenAnalysis/tools/") # You'll need to change this
import predRADs
import mete
import pln


def RandCompFast(Q, N):
    
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





def plot_comps(N, S):

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)

    for j in range(100):
        
        # Random composition
        RAD = RandCompFast(N, S)
        ranks = range(1, len(RAD)+1)
        plt.plot(ranks, RAD, lw = 1, c='gray')
        
    # Predicted geometric series
    predRAD = predRADs.get_GeomSeries(N, S, False)
        
    ranks = range(1, S+1)
    plt.plot(ranks, predRAD, lw = 2, c='Lime')    
        
    plt.yscale('log')
    plt.ylabel('Abundance')
    plt.xlabel('Rank in abundance')
    plt.title('Smallest abundance of MaxEnt prediction for SAD constrained\nby N and S (Geometric Series):'+str(min(predRAD)))
    plt.show()
    
    return


def plot_RADs_canonical(N, S):

    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    # Predicted geometric series
    predRAD = predRADs.get_GeomSeries(N, S, False) # False mean no zeros allowed
    ranks = range(1, S+1)
    plt.plot(ranks, predRAD, lw = 1, c='m')
    
    # Predicted log-series
    logSeries = mete.get_mete_rad(S, N)
    predRAD = logSeries[0]
    ranks = range(1, S+1)
    plt.plot(ranks, predRAD, lw = 1, c='c')
        
    # Predicted PLN
    #predRAD = pln.get_rad_from_obs(predRAD, 'pln') 
    #ranks = range(1, len(predRAD)+1)
    #plt.plot(ranks, predRAD, lw = 1, c='gray')   
    
    plt.yscale('log')
    plt.show()
    
    return



N = 10**4
S = 10**2
#plot_comps(N, S)

plot_RADs_canonical(N, S)