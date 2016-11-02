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


gd.stratifyDataBootstrap(remove_obs = 0)
#gd.stratifyDataBootstrap(remove_obs = 1)
#gd.stratifyDataBootstrap(remove_obs = 0, seqSim = '95')
#gd.stratifyDataBootstrap(remove_obs = 0, seqSim = '97')
#gd.stratifyDataBootstrap(remove_obs = 0, seqSim = '99')
#gd.stratifyDataOnce(['EMPclosed', 'MGRAST', 'HMP'], remove_obs = 0)
