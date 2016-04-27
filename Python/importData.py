import numpy as np

def import_obs_pred_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    if '75' in input_filename:
        data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8", names = ['site','obs', 'predPln', 'pred7525'], delimiter = " ")
    else:
        data = np.genfromtxt(input_filename, dtype = "f8,f8,f8", names = ['site','obs','pred'], delimiter = " ")
    #test = data[0:10000]
    #return test
    return data

def import_subsampled_data(input_filename):
    if ('zipf' in input_filename):
        # 33 for zipf
        # this needs to be fixesd, I put the file name twice in old code
        data = np.genfromtxt(input_filename, \
        dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
        names = ['site','site2','N0','S0','Nmax', \
        'N0_05','S0_05','Nmax_05', 'r2_05', 'gamma_05', \
        'N0_025','S0_025','Nmax_025', 'r2_025', 'gamma_025', \
        'N0_0125','S0_0125','Nmax_0125', 'r2_0125', 'gamma_0125', \
        'N0_00625','S0_00625','Nmax_00625', 'r2_00625', 'gamma_00625', \
        'N0_003125','S0_003125','Nmax_003125', 'r2_003125', 'gamma_003125',
        'N0_0015625','S0_0015625','Nmax_0015625','r2_0015625', 'gamma_0015625'], \
        delimiter = " ")
    else:
        # 27 columns for mete and geom
        data = np.genfromtxt(input_filename, \
        dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
        names = ['site','N0','S0','Nmax', \
        'N0_05','S0_05','Nmax_05','r2_05', \
        'N0_025','S0_025','Nmax_025','r2_025', \
        'N0_0125','S0_0125','Nmax_0125','r2_0125', \
        'N0_00625','S0_00625','Nmax_00625','r2_00625', \
        'N0_003125','S0_003125','Nmax_003125','r2_003125', \
        'N0_0015625','S0_0015625','Nmax_0015625','r2_0015625'], \
        delimiter = " ")
    return data

def import_subsampled_data_pandas(input_filename):
    if ('zipf' in input_filename):
        names = ['site','site2','N0','S0','Nmax', \
        'N0_05','S0_05','Nmax_05', 'r2_05', 'gamma_05', \
        'N0_025','S0_025','Nmax_025', 'r2_025', 'gamma_025', \
        'N0_0125','S0_0125','Nmax_0125', 'r2_0125', 'gamma_0125', \
        'N0_00625','S0_00625','Nmax_00625', 'r2_00625', 'gamma_00625', \
        'N0_003125','S0_003125','Nmax_003125', 'r2_003125', 'gamma_003125',
        'N0_0015625','S0_0015625','Nmax_0015625','r2_0015625', 'gamma_0015625']
        #data_table = pd.read_table(input_filename, names = names, header = None, sep=' ')

    else:
        names = ['site','N0','S0','Nmax', \
        'N0_05','S0_05','Nmax_05','r2_05', \
        'N0_025','S0_025','Nmax_025','r2_025', \
        'N0_0125','S0_0125','Nmax_0125','r2_0125', \
        'N0_00625','S0_00625','Nmax_00625','r2_00625', \
        'N0_003125','S0_003125','Nmax_003125','r2_003125', \
        'N0_0015625','S0_0015625','Nmax_0015625','r2_0015625']
    data_table = pd.read_table(input_filename, names = names, header = None, sep=' ')
    return data_table


def import_NSR2_data(input_filename):   # TAKEN FROM THE mete_sads.py script used for White et al. (2012)
    input_filename_str = str(input_filename)
    #NSR2_method = input_filename_split[-4]
    #method = str(NSR2_method.split('/')[1])
    if 'Stratified' in input_filename_str:
        data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
            names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                'evennessPred', 'skewnessObs', 'skewnessPred','R2'], delimiter = " ")
    else:
        if 'HMP' in input_filename_str:
            if ('zipf' in input_filename_str) :
                if ('glm' in input_filename_str) :


                    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                        names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                            'evennessPred', 'skewnessObs', 'skewnessPred','R2', 'NAP'], delimiter = " ")
                else:
                    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                        names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                            'evennessPred', 'skewnessObs', 'skewnessPred','R2', 'NAP'], delimiter = " ")
            else:
                data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                    'evennessPred', 'skewnessObs', 'skewnessPred', 'R2','NAP'], delimiter = " ")
        else:
            if 'zipf' in input_filename_str:
                if ('glm' in input_filename_str) :
                    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                    names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                        'evennessPred', 'skewnessObs', 'skewnessPred', 'R2'], delimiter = " ")

                else:
                    data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                    names = ['site','N','S', 'NmaxObs', 'NmaxPred', 'evennessObs', \
                        'evennessPred', 'skewnessObs', 'skewnessPred','gamma', 'R2'], delimiter = " ")
                        # 'gammma'
            else:
                data = np.genfromtxt(input_filename, dtype = "f8,f8,f8,f8,f8,f8,f8,f8,f8,f8", \
                names = ['site','N','S','NmaxObs', 'NmaxPred', 'evennessObs', \
                    'evennessPred', 'skewnessObs', 'skewnessPred', 'R2'], delimiter = " ")

    return data
