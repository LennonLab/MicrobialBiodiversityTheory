from __future__ import division
import os
import scipy as sp
import numpy as np
import pandas as pd
import sys

### Regenerate SbyS matrix with EMP site IDs as the row names


# Pull in metadata

# keep these columns: ENV_FEATURE, STUDY_ID, COUNTRY, COMMON_NAME, Description_duplicate, TITLE, ENV_BIOME, ENV_MATTER, TAXON_ID, Description



# Pull in SbyS





# Pull in NSR2 data
mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/NSR2/'
#print mydir
#sys.exit()

def import_NSR2_data(input_filename):

    input_filename_str = str(input_filename)
    if 'zipf' in input_filename_str:

        data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8", \
        names = ['site','N','S', 'Nmax','gamma','R2'], delimiter = " ")

    return data

filename = mydir+"zipf_EMPclosed_NSR2.txt"
dat = import_NSR2_data(filename)

dat = dat[dat['R2'] < 0.2]


# Pull in new Site
