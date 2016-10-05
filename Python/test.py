from __future__ import division
import mete as mete
import os, math
import numpy as np
import modelsAndMetrics as mo
import macroeco_distributions as md
from scipy import stats
import generateData as gd

mydir = os.path.expanduser("~/github/MicroMETE/")

import importData as im


gd.stratifyData(['EMPclosed','HMP', 'MGRAST'], remove_obs = 1)

#gd.stratifyData(['EMPclosed'], remove_obs = 1)
