from __future__ import division
import mete as mete
import os, math
import numpy as np
import modelsAndMetrics as mo
import macroeco_distributions as md
from scipy import stats
import generateData as gd
import generateFigures as gf

mydir = os.path.expanduser("~/github/MicroMETE/")

import importData as im

gf.figS5()
#gd.stratifyData1000(remove_obs = 0)
