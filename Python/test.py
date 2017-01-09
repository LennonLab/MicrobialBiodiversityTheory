from __future__ import division
import mete as mete
import os, math
import numpy as np
import modelsAndMetrics as mo
import macroeco_distributions as md
from scipy import stats
import generateData as gd
import generateFigures as gf
import mete

mydir = os.path.expanduser("~/github/MicroMETE/")

import importData as im


#methods = ['geom', 'mete']
#gd.generate_obs_pred_data(['EMPclosed','HMP', 'MGRAST'], methods, remove_obs = 0)
#gd.generate_obs_pred_data(['EMPclosed','HMP', 'MGRAST'], methods, remove_obs = 1)
#gd.generate_obs_pred_data(['95','97', '99'], methods, remove_obs = 0)
#gd.stratifyDataBootstrap(remove_obs = 0)
#gd.stratifyDataBootstrap(remove_obs = 1)
#gd.stratifyDataBootstrap(remove_obs = 0, seqSim = '95')
#gd.stratifyDataBootstrap(remove_obs = 0, seqSim = '97')
#gd.stratifyDataBootstrap(remove_obs = 0, seqSim = '99')
#gd.stratifyDataOnce(['EMPclosed', 'MGRAST', 'HMP'], remove_obs = 0)
#gf.tableS4()

#gd.get_AICc()
#gd.subsample_AICc()
#gf.figS6()
#gf.
