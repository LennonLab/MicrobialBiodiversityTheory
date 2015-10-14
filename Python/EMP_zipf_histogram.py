#!/usr/bin/env sage -python
from __future__ import division
import sys
import os
import pickle
import multiprocessing
import time
#sys.path.append(mydir)
mydir = os.path.expanduser("~/github/MicroMETE/data/")
import matplotlib.gridspec as gridspec
import signal


import scipy as sp
import  matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid.inset_locator import inset_axes
import macroecotools
import mete
import macroeco_distributions as md
from sklearn.preprocessing import PolynomialFeatures
from sklearn import linear_model
from scipy import stats


def import_obs_pred_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8", names = ['site','obs','pred'], delimiter = " ")
    #test = data[0:10000]
    #return test
    return data

obs_pred_data = import_obs_pred_data(mydir + 'ObsPred/zipf_EMPopen_obs_pred_all.txt')

site = ((obs_pred_data["site"]))
obs = ((obs_pred_data["obs"]))
pred = ((obs_pred_data["pred"]))

#axis_min = -0.5 * min(pred)
#axis_max = 100  * max(pred)


#bins = np.linspace(0, 1000000, 100)

#plt.hist(obs, alpha=0.5)
print pred
print max(pred)
print max(obs)
plt.hist(pred,bins=np.arange(min(pred), max(pred) + 10, 100), alpha=0.5)


plt.savefig("hist.png")
plt.close()
