# Run script.
### This file will contain all the commands to automatically rerun all analyses.

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Run analyses and generate figures for Broken-stick, METE, and Zipf.")
    parser.add_argument('-a', type = str, default = "no", help = "Do you want to run the full analysis for the figure?")
    parser.add_argument('-f', type = str, default = "1", help = "Which figure do you wnat to generate?")
    parser.add_argument('-r', type = str, default = "no", help = "Do you want to rerun all the analyses and generate all the figures?")
    input_options = parser.parse_args()
    generate_data = input_options.a.lower()
    generate_fig = input_options.f.upper()
    generate_all = input_options.r.lower()
    data_points = 352899 # number of obs_pred datapoints to plot ( HMP has ~352899 )
    size = 0
    methods = ['geom', 'mete','zipf']
    params = ['N','S', 'N/S']

    #mydir = os.path.expanduser("~/github/MicroMETE/data/")
    #mydir = os.path.dirname(os.path.realpath(__file__))
    if generate_all == 'yes':
            datasets = [ 'HMP','EMPclosed','EMPopen','MGRAST','95', '97','99']
            generate_obs_pred_data(datasets_all, methods, size)
            datasets1 = [ 'HMP','EMPclosed','MGRAST']
            plot_obs_pred_sad(methods, datasets1, data_points)
            datasets2 = ['HMP']
            NSR2_regression(methods, datasets2, data_dir= mydir)
            datasetsS1 = [ 'HMP','EMPclosed','EMPopen']
            plot_obs_pred_sad(methods, datasetsS1, data_points)
            datasetsS2 = ['95', '97','99']
            plot_obs_pred_sad(methods, datasetsS2, data_points)
            datasetsS3 = ['EMPopen']
            NSR2_regression(methods, datasetsS3, data_dir= mydir)
            datasetsS4 = ['HMP']
            NSR2_regression(methods, datasetsS4, data_dir= mydir)
            datasetsS5 = ['MGRAST']
            NSR2_regression(methods, datasetsS5, data_dir= mydir)
            datasetsS6 = ['95']
            NSR2_regression(methods, datasetsS6, data_dir= mydir)
            datasetsS7 = ['97']
            NSR2_regression(methods, datasetsS7, data_dir= mydir)
            datasetsS8 = ['99']
            NSR2_regression(methods, datasetsS8, data_dir= mydir)
            datasetsS9 = [ 'HMP','EMPclosed','EMPopen', 'MGRAST']
            zipf_mle_plots(datasetsS9, data_dir= mydir)
            datasetsS10 = ['95', '97','99']
            zipf_mle_plots(datasetsS10, data_dir= mydir)

    else:
        if generate_data == 'yes':
            if generate_fig == '1':
                datasets = [ 'HMP','EMPclosed','MGRAST']
                generate_obs_pred_data(datasets, methods, size)
                plot_obs_pred_sad(methods, datasets, data_points)
            elif generate_fig == '2':
                datasets = ['HMP']
                generate_obs_pred_data(datasets, methods, size)
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S1':
                datasets = [ 'HMP','EMPclosed','EMPopen']
                generate_obs_pred_data(datasets, methods, size)
                plot_obs_pred_sad(methods, datasets, data_points)
            elif generate_fig == 'S2':
                datasets = ['95', '97','99']
                generate_obs_pred_data(datasets, methods, size)
                plot_obs_pred_sad(methods, datasets, data_points)
            elif generate_fig == 'S3':
                datasets = ['EMPopen']
                generate_obs_pred_data(datasets, methods, size)
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S4':
                datasets = ['HMP']
                generate_obs_pred_data(datasets, methods, size)
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S5':
                datasets = ['MGRAST']
                generate_obs_pred_data(datasets, methods, size)
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S6':
                datasets = ['95']
                generate_obs_pred_data(datasets, methods, size)
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S7':
                datasets = ['97']
                generate_obs_pred_data(datasets, methods, size)
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S8':
                datasets = ['99']
                generate_obs_pred_data(datasets, methods, size)
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S9':
                datasets = [ 'HMP','EMPclosed','EMPopen', 'MGRAST']
                generate_obs_pred_data(datasets, methods, size)
                zipf_mle_plots(datasets, data_dir= mydir)
            elif generate_fig == 'S10':
                datasets = ['95', '97','99']
                generate_obs_pred_data(datasets, methods, size)
                zipf_mle_plots(datasets, data_dir= mydir)
            else:
                print "Input not valid."

        elif generate_data == 'no':
            if generate_fig == '1':
                datasets = [ 'HMP','EMPclosed','MGRAST']
                plot_obs_pred_sad(methods, datasets, data_points)
            elif generate_fig == '2':
                datasets = ['HMP']
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S1':
                datasets = [ 'HMP','EMPclosed','EMPopen']
                plot_obs_pred_sad(methods, datasets, data_points)
            elif generate_fig == 'S2':
                datasets = ['95', '97','99']
                plot_obs_pred_sad(methods, datasets, data_points)
            elif generate_fig == 'S3':
                datasets = ['EMPopen']
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S4':
                datasets = ['HMP']
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S5':
                datasets = ['MGRAST']
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S6':
                datasets = ['95']
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S7':
                datasets = ['97']
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S8':
                datasets = ['99']
                NSR2_regression(methods, datasets, data_dir= mydir)
            elif generate_fig == 'S9':
                datasets = [ 'HMP','EMPclosed','EMPopen', 'MGRAST']
                zipf_mle_plots(datasets, data_dir= mydir)
            elif generate_fig == 'S10':
                datasets = ['95', '97','99']
                zipf_mle_plots(datasets, data_dir= mydir)
            else:
                print "Input not valid."
