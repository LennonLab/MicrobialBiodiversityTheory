import argparse
import generateFigures as gf
import generateData as gd
### This file will contain all the commands to automatically rerun all analyses.

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Run analyses and/or generate figures.")
    parser.add_argument('-a', '--analysis', default = False, action='store_true',\
        help = 'Do you want to re-run the analysis')
    parser.add_argument('-f', type = str, default = '0', help = "Which figure (T) or table (T) do you wnat to generate?")

    args = parser.parse_args()
    analysis = args.analysis
    figure = args.f.upper()

    if analysis == False:
        if figure == 'T0':
            pass
        elif figure == 'F0':
            pass
        elif figure == 'F1':
            gf.fig1()
        elif figure == 'F2':
            gf.fig2()
        elif figure == 'F3':
            gf.fig3()
        elif figure == 'F4':
            gf.fig4()
        elif figure == 'FS1':
            gf.figS1()
        elif figure == 'FS2':
            gf.figS2()
        elif figure == 'FS3':
            gf.figS3()
        elif figure == 'FS4':
            gf.figS4()
        elif figure == 'FS5':
            gf.figS5()
        elif figure == 'FS6':
            gf.figS6()
        elif figure == 'T1':
            gf.table1()
        elif figure == 'T2':
            gf.table2()
        elif figure == 'TS1':
            gf.tableS1()
        elif figure == 'TS2':
            gf.tableS2()
        elif figure == 'TS3':
            gf.tableS3()
        elif figure == 'TS4':
            gf.tableS4()
        elif figure == 'TS5':
            gf.tableS5()
        elif figure == 'TS6':
            gf.tableS6()
        else:
            print "Command not recognized"
    else:
        gd.HMP_OTU_to_sparese_SbyS()
        gd.Match_NAP_to_sparse_SbyS()
        gd.get_SADs_HMP()
        gd.get_SADs_mgrast()
        gd.merge_SADs_mgrast()

        methods = ['geom','lognorm', 'mete','zipf']
        gd.generate_obs_pred_data(['EMPclosed','HMP', 'MGRAST'], methods, remove_obs = 0)
        gd.generate_obs_pred_data(['EMPclosed','HMP', 'MGRAST'], methods, remove_obs = 1)
        gd.generate_obs_pred_data(['95','97', '99'], methods, remove_obs = 0)
        gd.stratifyDataOnce(['EMPclosed','HMP', 'MGRAST'], remove_obs = 0)
        gd.stratifyDataOnce(['EMPclosed','HMP', 'MGRAST'], remove_obs = 1)
        gd.stratifyDataBootstrap(remove_obs = 0)
        gd.stratifyDataBootstrap(remove_obs = 1)
        gd.stratifyDataBootstrap(remove_obs = 0, seqSim = '95')
        gd.stratifyDataBootstrap(remove_obs = 0, seqSim = '97')
        gd.stratifyDataBootstrap(remove_obs = 0, seqSim = '99')
