import pandas as pd
import numpy as np
import os
import generate_figs_zipf as gf
# Get sample labels of interest.

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'


#metadata = pd.read_table('~/github/MicroMETE/data/HMP-Data/ppAll_V35_map.txt', sep='\t', index_col=False)
#metadata[['VisitNo']] = metadata[['VisitNo']].astype(float)
#metadata = metadata.loc[metadata['VisitNo'] == float(1)]
#metadata = metadata[np.isfinite(metadata['RSID'])]
#metadata.to_csv("../data/HMP-Data/ppAll_V35_map_noTimeseries.txt", sep='\t', index=False)
# Pull the RSID number
#RSID = metadata[['RSID']]
# Import the count data
def HMP_OTU_to_sparese_SbyS():
    otu_count = pd.read_table('~/github/MicroMETE/data/HMP-Data/hmp1.v35.hq.otu.counts.txt', sep='\t', index_col=False)
    sparsesbys = pd.melt(otu_count, id_vars=['collection'])
    sparsesbys = sparsesbys[sparsesbys['value'] > 0]
    sparsesbys.columns = ['Sample', 'OTU', 'Count']
    sparsesbys.to_csv("../data/HMP-Data/HMPsparseSbyS.txt", sep='\t', index=False)

def Match_NAP_to_sparse_SbyS(timeseries):
    IN = mydir + 'HMP-Data/ppAll_V35_map.txt'
    metadata = pd.read_table(IN, sep='\t', index_col=False)
    timeseries = bool(timeseries)
    metadata = metadata[np.isfinite(metadata['NAP'])]
    metadata[['NAP']] = metadata[['NAP']].astype(str)
    metadata.drop_duplicates(cols='NAP', take_last=True)
    if timeseries == True:
        pass
    else:
        metadata[['VisitNo']] = metadata[['VisitNo']].astype(float)
        metadata = metadata.loc[metadata['VisitNo'] == float(1)]
        metadata.to_csv("../data/HMP-Data/ppAll_V35_map_noTimeseries.txt", sep='\t', index=False)
    print metadata.shape
    # Pull the RSID number & convert to numpy array
    NAPs = metadata[['NAP']].values
    NAPs = NAPs.flatten()
    #print NAPs
    #RSIDs = map(str, RSIDs)
    IN_OTU = mydir + 'HMP-Data/HMPsparseSbyS.txt'
    if timeseries == True:
        OUT = open(mydir + 'HMP-Data/HMPsparseSbyS_NAP.txt','w+')
    else:
        OUT = open(mydir + 'HMP-Data/HMPsparseSbyS_NAP_noTimeseries.txt','w+')
    for j, line in enumerate(open(IN_OTU)):
        if '.PPS' in line:
            continue
        else:
            line_split = line.split()
            OTU = line_split[1]
            count = line_split[2]
            NAP = str(line_split[0].split('.')[0])
            try:
                float(NAP)
                NAP_test = str(NAP) + '.0'
                if NAP_test  in NAPs:
                    print>> OUT, line_split[0], OTU, count
            except ValueError:
                print "Not a float"


        #if ('control' not in RSID) and ('water_blank' not in RSID):
        #    print float(RSID)
def get_SADs_HMP(path, timeseries):
    timeseries = bool(timeseries)
    if timeseries == True:
        IN = path + 'HMP-Data/HMPsparseSbyS_NAP.txt'
        OUT =  open(path+'HMP-Data/' + 'HMP-SADs_NAP.txt', 'w+')
    else:
        IN = path + 'HMP-Data/HMPsparseSbyS_NAP_noTimeseries.txt'
        OUT =  open(path+'HMP-Data/' + 'HMP-SADs_NAP_noTimeseries.txt', 'w+')
    SADdict = {}
    with open(IN) as f:
        for d in f:
            if d.strip():
                d = d.split()
                site = d[0]
                abundance = int(d[-1])
                if abundance > 0:
                    if site not in SADdict:
                        #SADdict[site] = [abundance]
                        SADdict[site] = [abundance]
                    else:
                        SADdict[site].append(abundance)

    #SADs = SADdict.values()
    filtered_SADdict = {}
    #filteredSiteNames = []

    for key, value in SADdict.iteritems():
        if len(value) >= 10:
            filtered_SADdict[key] = value

    print "You have " + str(len(SADdict)) + " sites"

    # first value of the line is the site name, rest is SAD
    # site name filtered out in generate obs pred data
    for key, value in filtered_SADdict.iteritems():
        output = value.insert(0,key)
        #value = ", ".join(value)
        print>> OUT, value

#Match_NAP_to_sparse_SbyS(False)
#get_SADs_HMP(mydir)
#datasets = ['HMP']
#methods = ['geom', 'mete','zipf']
get_SADs_HMP(mydir, False)
get_SADs_HMP(mydir, True)
