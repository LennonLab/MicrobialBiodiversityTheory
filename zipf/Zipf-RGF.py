from __future__ import division
import numpy as np
import sys
import math
import os


def get_pred_iterative(pmf, S):
    """Function to get predicted abundances (reverse-sorted) for distributions with no analytical ppf."""

    cdf = [(S - i + 0.5) / S for i in range(1, S + 1)]
    cdf = np.sort(cdf)
    rad = []

    j = 0
    i = 1
    cdf_cum = 0

    # print statements below were used as checks in writing the code
    while j < len(cdf):
        #print cdf_cum, len(pmf), i
        cdf_cum += pmf[i-1]

        while cdf_cum >= cdf[j]:
            #print cdf_cum, cdf[j], rad.count(0)

            rad.append(i)
            #print 'N:', N, 'S:', S, 'species:', len(rad), 'abundance:',i, 'cum_abs:', sum(rad)#, 'kmin:', min(rad), 'kmax:', max(rad)

            j += 1
            if j == len(cdf):
                rad.reverse()
                return np.array(rad)
        i += 1


def get_emp_cdf(pmf):
    """Compute the empirical cdf given a list or an array"""
    pmf = np.array(pmf)
    cdf = []

    for point in pmf:
        point_cdf = len(pmf[pmf <= point]) / len(pmf)
        cdf.append(point_cdf)
    return np.array(cdf)



""" Note: In their paper, the authors M is our N. Their N is our S. Their kmax is our Nmax."""

def zipf_rgf_params(obs_rad):

    N = sum(obs_rad)
    S = len(obs_rad)

    kmin = min(obs_rad)
    Nmax = max(obs_rad)
    avg_ab = N/S

    #Set gamma and upper bound of b.
    gamma = 1.99999
    b = 1

    _sum = 0
    for k in range(kmin, N):
        pk = math.exp(-b*k)/k**gamma
        if pk > 0:
            _sum += pk
        else:
            break

    A = 1/_sum

    #Find the best b.
    Nmaxtmp = N # initialize Nmaxtmp as N
    b0 = 2*b
    b1 = 0

    while(abs(Nmaxtmp - Nmax) > 1):
        b = (b0 + b1)/2
        sum1 = 0
        kc = 0

        for k in range(N, kmin,-1):
            sum1 += A * math.exp(-b*k)/k**gamma
            if sum1 > 1/S:
                kc = k
                break

        sum1 = 0
        sum2 = 0

        for k in range(kc, N):
            s1=k*math.exp(-b*k)/k**gamma
            s2=math.exp(-b*k)/k**gamma

            if s1 <= 0 or s2 <= 0:
                break

            sum1 += s1
            sum2 += s2

        Nmaxtmp = sum1/sum2

        #print Nmaxtmp, Nmax,'b =', b

        if Nmaxtmp > Nmax:
            b1 = b

        else:
            b0 = b

    sum1 = 0
    sum2 = 0

    for k in range(kmin, N):
        sum1 += math.exp(-b*k)/k**gamma
        sum2 += k * math.exp(-b*k)/k**gamma

    A = 1/sum1
    #kavm = sum2/sum1

    #print gamma,'\t', b,'\t',A,'\t',kavm,'\t',avg_ab  #gamma, b, A, modeling
    #Compare modelling <k> and real <k>. If they are different, guess another gamma (increasing gamma will decrease <k>).

    return [gamma, b, A, N]




def zipf_pmf(gamma, b, A, N):

    pmf = []
    k = 1
    while k <= N:
        pK = A * np.exp(-b*k) / (k**gamma)
        pmf.append(pK)
        k += 1

    return pmf



def zipf_rgf(obs_rad):

    S = len(obs_rad)
    gamma, b, A, N = zipf_rgf_params(obs_rad)
    pmf = zipf_pmf(gamma, b, A, N) # pmf for non-built-in function
    #print len(pmf),'\t', pmf[0],'\t', pmf[-1]

    rad = get_pred_iterative(pmf, S)
    #print 'predN', sum(rad), 'predS', len(rad), 'Nmax', max(rad), rad[0], 'Nmin', min(rad), rad[-1]

    return rad


# Uncomment these lines to informally test the code

#rad = [500, 400, 350, 320, 310, 300, 250, 240, 230, 220, 200, 180, 160, 140, 130, 120, 100, 80, 70, 60, 40, 20, 10, 9,9,9,9,9,9,9,9,9,9,9,9, 5, 5,5,5,5,5,2,2,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
#print 'N', sum(rad), 'S', len(rad), 'Nmax', max(rad), rad[0], 'Nmin', min(rad), rad[-1]

#rad = zipf_rgf(rad)
#print 'N', sum(rad), 'S', len(rad), 'Nmax', max(rad), rad[0], 'Nmin', min(rad), rad[-1]
