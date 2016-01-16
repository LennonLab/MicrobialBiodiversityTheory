import generate_figs_zipf as gf
import macroeco_distributions as md
from scipy import stats
import signal
from scipy import stats, optimize
import os
import macroeco_distributions as md
import macroecotools
import numpy as np
import mete
import time

mydir = os.path.dirname(os.path.realpath(__file__))
mydir = str(mydir[:-6]) + 'data/'

estimators = ['fmin']

def test_zipf_num_est(datasets, estimators, SAD_number, iterations, fail_threshold):
    percents = [0.500000, 0.250000, 0.125000, 0.062500, 0.031250, 0.015625]
    for dataset in datasets:
        signal.signal(signal.SIGALRM, gf.timeout_handler)
        if dataset == 'MGRAST':
            # fix subset l8r
            IN = mydir  + dataset + '-Data' + '/MGRAST/MGRAST-SADs.txt'
            nsr2_data_zipf = gf.import_NSR2_data(mydir + 'NSR2/' + 'zipf_MGRAST_NSR2.txt')
        elif dataset == '95' or dataset == '97' or dataset == '99':
            IN = mydir  + dataset + '-Data/' + str(dataset) + '/MGRAST-' + str(dataset) + '-SADs.txt'
            nsr2_data_zipf = gf.import_NSR2_data(mydir + 'NSR2/' +'zipf_MGRAST'+dataset+'_NSR2.txt')
        elif dataset == 'HMP':
            IN = mydir  + dataset + '-Data' + '/' + dataset +'-SADs_NAP.txt'
            nsr2_data_zipf = gf.import_NSR2_data(mydir + 'NSR2/' +  'zipf_'+dataset+'_NSR2.txt')
        else:
            IN = mydir  + dataset + '-Data' + '/' + dataset +'-SADs.txt'
            nsr2_data_zipf = gf.import_NSR2_data(mydir + 'NSR2/' +  'zipf_'+dataset+'_NSR2.txt')

        nsr2_data_zipf_N_site = np.column_stack((nsr2_data_zipf["site"], nsr2_data_zipf["N"]))
        # Sort these arrays
        nsr2_data_zipf_sorted = nsr2_data_zipf_N_site[nsr2_data_zipf_N_site[:,1].argsort()[::-1]]
        nsr2_data_zipf_top100 = nsr2_data_zipf_sorted[:SAD_number,]
        # Get the SAD numbers
        zipf_numbers = nsr2_data_zipf_top100[:,0]
        zipf_numbers = zipf_numbers.astype(int)
        successful_SADs_samplings = SAD_number
        for estimator in estimators:
            OUT = open(mydir + 'SubSampled-Data' + '/' + dataset + '_zipf_' + \
                str(estimator) + '_SubSampled_Data.txt', 'w+')
            num_lines = sum(1 for line in open(IN))
            test_lines = 0
            succeess_lines = SAD_number
            while succeess_lines > 0:
                site = nsr2_data_zipf_sorted[test_lines,0]
                for j,line in enumerate(open(IN)):
                    if (j != site):
                        continue
                    else:
                        if dataset == "HMP":
                            line = line.strip().split(',')
                            line = [x.strip(' ') for x in line]
                            line = [x.strip('[]') for x in line]
                            site_name = line[0]
                            line.pop(0)
                        else:
                            line = eval(line)
                    obs = map(int, line)
                    # Calculate relative abundance of each OTU
                    # Use that as weights
                    N_0 = float(sum(obs))
                    S_0 = len(obs)
                    N_max = max(obs)
                    if S_0 < 10 or N_0 <= S_0:
                        test_lines += 1
                        continue
                    line_ra = map(lambda x: x/N_0, obs)
                    sample_sizes = map(lambda x: round(x*N_0), percents)
                    if any(sample_size <= 10 for sample_size in sample_sizes)  == True:
                        test_lines += 1
                        continue
                    zipf_means = [N_0, S_0, N_max]
                    failed_percents = 0
                    for k, percent in enumerate(percents):
                        if failed_percents > 0:
                            continue
                        N_max_list_zipf = []
                        N_0_list_zipf = []
                        S_0_list_zipf = []
                        r2_list_zipf = []
                        gamma_list = []
                        iter_count_current = 0
                        iter_count = iterations
                        iter_failed = 0
                        while iter_count > 0 and iter_failed < fail_threshold:
                            sample_size_k = sample_sizes[0]
                            sample_k = np.random.multinomial(sample_size_k, line_ra, size = None)
                            sample_k_sorted = -np.sort( -sample_k[sample_k != 0] )
                            N_0_k = sum(sample_k_sorted)
                            S_0_k = sample_k_sorted.size
                            if S_k < 10 or N_k <= S_k:
                                continue
                            N_max_k = max(sample_k_sorted)
                            iter_count_current += 1
                            # Start the timer. Once 1 second is over, a SIGALRM signal is sent.
                            signal.alarm(2)
                            # This try/except loop ensures that
                            #   you'll catch TimeoutException when it's sent.
                            #start_time = time.time()
                            try:
                                # Whatever your function that might hang
                                zipf_class = gf.zipf(sample_k_sorted, estimator)
                                pred_tuple = zipf_class.from_cdf()
                                Zipf_solve_line = zipf_class.zipf_solver(sample_k_sorted)
                                rv = stats.zipf(Zipf_solve_line)
                                pred_zipf = pred_tuple[0]
                                gamma = pred_tuple[1]
                                r2_zipf = macroecotools.obs_pred_rsquare(np.log10(sample_k_sorted), np.log10(pred_zipf))
                                if (r2_zipf == -float('inf') ) or (r2_zipf == float('inf') ) or (r2_zipf == float('Nan') ):
                                    continue
                                else:
                                    r2_list_zipf.append(r2_zipf)
                                    gamma_list.append(gamma)
                                    N_max_list_zipf.append(N_max_k)
                                    N_0_list_zipf.append(N_0_k)
                                    S_0_list_zipf.append(S_0_k)

                            except gf.TimeoutException:
                                print "Line " + str(j) + ": " + str(estimator) + " timed out"
                                iter_count -= 1
                                if iter_failed >= fail_threshold:
                                    failed_percents += 1
                                iter_failed += 1
                                continue # continue the for loop if function takes more than x seconds
                            else:
                                iter_count -= 1
                                #print("--- %s seconds ---" % (time.time() - start_time))
                                # Reset the alarm
                                signal.alarm(0)


                        if len(N_0_list_zipf) != iterations:
                            test_lines += 1
                            continue
                        N_0_zipf_mean = np.mean(N_0_list_zipf)
                        zipf_means.append(N_0_zipf_mean)

                        S_0_zipf_mean = np.mean(S_0_list_zipf)
                        zipf_means.append(S_0_zipf_mean)

                        N_max_zipf_mean = np.mean(N_max_list_zipf)
                        zipf_means.append(N_max_zipf_mean)

                        r2_zipf_mean = np.mean(r2_list_zipf)
                        zipf_means.append(r2_zipf_mean)

                        gamma_zipf_mean = np.mean(gamma_list)
                        zipf_means.append(gamma_zipf_mean)

                    '''Now we check if the lists are the right length
                    there are 6 iterations for the percentage
                    mete/ geom, append four items each iteration.
                    4*6 = 24, add three original = 27
                    likewise, for zipf, (5*6) + 3 = 33 '''
                    if len(zipf_means) == 33:
                        test_lines += 1
                        succeess_lines -= 1
                        zipf_means_str = ' '.join(map(str, zipf_means))
                        #OUT1.write(','.join(map(repr, geom_means_str[i]))
                        print>> OUT, j, zipf_means_str
                        print "Line " + str(j) + ": " + str(succeess_lines) + " SADs to go!"
                    else:
                        test_lines += 1
                #print estimator
            print dataset


#datasets = ['MGRAST']
datasets = ['MGRAST']
#test_zipf_num_est(datasets, estimators, 100, 100, 20)
