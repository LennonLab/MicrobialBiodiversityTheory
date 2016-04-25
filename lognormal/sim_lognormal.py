from __future__ import division
import sys
import numpy as np
import math
import random
import matplotlib.pyplot as plt
import scipy.stats
from random import randrange, randint, uniform, seed, choice
import time


def GetRelAbs(RACsample): # a function to transform abundance into relative abundances for many RACs

    RACs = []
    for RAC in RACsample:
        RAC = (np.array(RAC) / sum(RAC)) * 100 # operate on a numpy array (vector multiplication) get relative abundance of each species
        RACs.append(RAC.tolist()) # convert the array back to a list and append to RACs

    return RACs



"""This script codes Broken Stick models for species abundance distribution.
    These are generally conceived as models for how species partition a niche
    dimension or shared and limited resource. -KL """


def SimBrokenStick(N, S, sample_size, rel=False):

    """
    A function to generate random samples of the simultaneous Broken Stick Model
    of MacArthur 1957 or 1961. We need the actual citation.

    N  :  Total abundance, i.e., number of individuals
    S  :  Species richness, i.e., number of species
    sample_size  :  number of rank-abundance vectors to generate

    How this model works.  Create list of range of N, then make S - 1 random
    divisions in range of N, repeat sample_size times. So, first we have to get
    our indices. Which you've figured out. Then, we can split a list of N 1's at
    those indices points. A few rules apply, e.g. 0 can't be one of the indices.
    Say the simultanesouly drawn indices are [2, 9, 5]. Here, simultaneous means
    those numbers were drawn without replacement, i.e., were not allowed to
    can't draw any number twice. Once you sort those numbers then you can...

    Hint: Let N = 15, S = 4, and SortedIndices = [2, 5, 9]
    So: oo|ooo|oooo|oooooo  = [2, 3, 4, 6] -> N = 15 and S = 4
    2-0 = 2  :  5-2 = 3  :  9-5 = 4  :  15-9 = 6

    """

    RACs = []
    RAC = []

    while len(RACs) < sample_size:

        RAC = []

        indices = []
        while len(indices) < S-1:
            index = random.randint(1, N-1)
            if index in indices: continue
            else: indices.append(index)

        indices.sort()
        indices.append(N)

        nsum = 0
        RAC = []
        for i in indices:
            i -= sum(RAC)
            RAC.append(i)
            nsum += i

        RAC.sort()
        RAC.reverse()

        RACs.append(RAC)

        """
        for _list in RACs:
            if sum(RAC) !=N or len(RAC) != S:
                print 'Incorrect N and S: N=',sum(RAC),' S=', len(RAC)
                sys.exit()
        """
    if rel == True: RACs = GetRelAbs(RACs)
    return RACs




def DomPreInt(N, S, sample_size, rel=False): # Works only with positive integers
    '''This script codes Tokeshi's Dominance Preemption Model
    this code does not work well with small N or high S'''
    sample = [] # A list of RACs
    fails = 0

    while len(sample) < sample_size: # number of samples loop

        if fails > 1000:
            print 'Too many attempts'
            return sample
            break

        RAC = [] #RAC is a list
        sp1 = randrange(int(round(N *.5)), N) #Rand select from N to N(.5)
        ab1 = N - sp1
        RAC.extend([sp1, ab1])

        while len(RAC) < S:
            ab2 = RAC.pop()

            if ab2 < S - len(RAC) or ab2 < 2:
                fails += 1
                break

            sp2 = randrange(int(round(ab2*.5)), ab2)
            RAC.extend([sp2, ab2-sp2])

        if len(RAC) == S and sum(RAC) == N:
            sample.append(RAC)

    if rel == True: sample = GetRelAbs(sample)
    return sample



def DomPreFloat(N, S, sample_size, rel=False):#Works with decimal values
    sample = [] # A list of RACs

    while len(sample) < sample_size:
        RAC = []
        sp1 = uniform((N *.5), N)
        ab1 = N - sp1
        RAC.extend([sp1, ab1])

        while len(RAC) < S:
            ab2 = RAC.pop()
            sp2 = uniform((ab2*.5), ab2)
            RAC.extend([sp2, ab2-sp2])

        sample.append(RAC)

    if rel == True: sample = GetRelAbs(sample)
    return sample



def SimLogNormInt(N, S, sample_size, rel=False):
    '''This script codes the Lognormal Model'''
    sample = []
    fails = 0

    while len(sample) < sample_size:

        n = int(round(0.75 * N))
        RAC = [n, N - n]

        if fails > 1000:
                print 'Too many attempts'
                return sample
                break

        while len(RAC) < S:

            ind = randrange(len(RAC))
            v = RAC.pop(ind)
            v1 = int(round(0.75 * v))
            v2 = v - v1   # forcing all abundance values to be integers

            if v1 < 1 or v2 < 1:
                fails += 1
                break  # forcing smallest abundance to be
                                        # greater than one
            RAC.extend([v1, v2])

        if len(RAC) == S and sum(RAC) == N:
            RAC.sort()
            RAC.reverse()
            sample.append(RAC)

    if rel == True: sample = GetRelAbs(sample)
    return sample


def SimLogNormFloat(N, S, sample_size, rel=False):
    '''This script codes the Lognormal Model'''
    sample = []
    while len(sample) < sample_size:

        n = 0.75 * N
        RAC = [n, N - n]

        while len(RAC) < S:

            ind = randrange(len(RAC))
            v = RAC.pop(ind)
            v1 = 0.75 * v
            v2 = v - v1   # forcing all abundance values to be integers

            RAC.extend([v1, v2])

        if len(RAC) == S and int(round(sum(RAC))) == N:
            RAC.sort()
            RAC.reverse()
            sample.append(RAC)

    if rel == True: sample = GetRelAbs(sample)
    return sample


def SimParetoFloat(N, S, sample_size, rel=False):
    '''This script codes the Pareto Model.  The divisions of N follow the 80:20
    rule.  The sections of N are selected at random for division. Returns
    decimal values.'''

    sample = []

    while len(sample) < sample_size:
        RAC = [0.8*N, 0.2*N]

        while len(RAC) < S:
            ind = randrange(len(RAC))
            v = RAC.pop(ind)
            v1 = int(round(0.8 * v))
            v2 = v - v1  # forcing all abundance values to be integers

            RAC.extend([v1, v2])

        RAC.sort(reverse = True)
        sample.append(RAC)

    if rel == True: sample = GetRelAbs(sample)
    return sample



def SimParetoInt(N, S, sample_size, rel=False):
    '''This script codes the Pareto Model.  The divisions of N follow the 80:20
    rule.  The sections of N are selected at random for division. Returns only
    integer values.'''

    sample = []
    fails = 0

    while len(sample) < sample_size:

        n = int(round(0.8 * N))
        RAC = [n, N - n]

        if fails > 1000:
                print 'Too many attempts'
                return sample
                break

        while len(RAC) < S:

            ind = randrange(len(RAC))
            v = RAC.pop(ind)
            v1 = int(round(0.8 * v))
            v2 = v - v1   # forcing all abundance values to be integers

            if v1 < 1 or v2 < 1:
                fails += 1
                break  # forcing smallest abundance to be
                                        # greater than one
            RAC.extend([v1, v2])

        if len(RAC) == S and sum(RAC) == N:
            RAC.sort()
            RAC.reverse()
            sample.append(RAC)

    if rel == True: sample = GetRelAbs(sample)
    return sample


def SimpleRandomFraction(N, S, sample_size, rel=False):

    """
    This function randomly and sequently splits N into two integers by
    starting with a vector where N is the only value, i.e. N, and then splitting
    N into two positive integers at random, [N] -> [x1, x2], where x1+x2 = N.
    The function then draws one of the integers at random from the new vector
    (having 2 values), and then splits the integer into two other integers.
    At this point, the vector has three values [x1, x2, x3], where x1+x2+x3 = N.
    This process keeps on until there are a number of integers equal to S.

    N  :  total abundance; number of individuals in the community
    S  :  species richness; number of species in the community
    sample_size  :  number of random rank-abundance vectors to generate
    """

    sample = []
    for i in range(sample_size):
        RAC = [N]

        while len(RAC) < S:

            sp1_ab = choice(RAC) # pick a species (i.e. vector value) abundance at random
            if sp1_ab == 1:
                continue # this is a control statement to prevent the program
                # from encountering an impossible conditions, i.e., dividing 1
                # at random into two positive integers

            sp2_ab = randrange(1, sp1_ab) # pick a random number (abundance) between 1 and sp1_ab - 1, inclusive
            sp1_index = RAC.index(sp1_ab) # index in the RAC of the species we picked

            RAC[sp1_index] = sp1_ab - sp2_ab # decrease its abundance according to sp_ab
            RAC.append(sp2_ab)

        RAC.sort(reverse = True)  #sorts RAC's in decending order to aid in graphing.
        sample.append(RAC) # appending a list (i.e. an RAC) to another list

    if rel == True: sample = GetRelAbs(sample)
    return sample



def DomDecayFloat(N, S, sample_size, rel=False):
    '''This model randomly divides N into 2 sections, then the largest section of N is then
    divided at random.  This continues until there are S divisions. This model
    returns decimals values.'''

    sample = [] # A list of RACs

    while len(sample) != sample_size: # number of samples loop

        RAC = [] #RAC is a list
        sp1 = uniform(0,  N * .5) #Rand select from N to N(.5)
        ab1 = N - sp1
        RAC.extend([sp1, ab1])

        while len(RAC) < S:
            ab2 = RAC.pop()
            sp2 = uniform(0, ab2 *.5)
            RAC.extend([sp2, ab2-sp2])

        if len(RAC) == S and int(round(sum(RAC))) == N:
            RAC.sort(reverse = True)
            sample.append(RAC)

    if rel == True: sample = GetRelAbs(sample)
    return sample



def DomDecayInt(N, S, sample_size, rel=False): # Works only with positive integers    ...if you think about it, there can never be two 1's
    '''This model randomly divides N into 2 sections, then the largest section of N is then
    divided at random.  This continues until there are S divisions. This model
    returns only integer values.'''

    sample = [] # A list of RACs
    fails = 0

    while len(sample) < sample_size: # number of samples loop

        RAC = [] #RAC is a list
        sp1 = randint(1, int(round(N*.5)))
        ab1 = N - sp1

        if fails > 1000:
            print 'too many attempts'
            return sample
            break

        RAC.extend([sp1, ab1])

        while len(RAC) < S:
            if min(RAC) < 2: break

            ab2 = RAC.pop()
            sp2 = randint(1, int(round(ab2 * .5)))
            RAC.extend([sp2, ab2 - sp2])

        if len(RAC) == S and sum(RAC) == N:
            RAC.sort(reverse = True)
            sample.append(RAC)

    if rel == True: sample = GetRelAbs(sample)
    return sample
