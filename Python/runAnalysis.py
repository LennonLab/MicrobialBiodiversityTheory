import argparse
import generateFigures as gf
import generateData as gd
### This file will contain all the commands to automatically rerun all analyses.

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = "Run analyses and/or generate figures.")
    parser.add_argument('-a', '--analysis', default = False, action='store_true',\
        help = 'Do you want to re-run the analysis')
    parser.add_argument('-f', type = str, default = '0', help = "Which figure do you wnat to generate?")

    args = parser.parse_args()
    analysis = args.analysis
    figure = args.f.upper()

    if analysis == False:
        if figure == '0':
            pass
        elif figure == '1':
            gf.fig1()
        elif figure == '2':
            gf.fig2()
        elif figure == '3':
            gf.fig3()
        elif figure == '4':
            gf.fig4()
        elif figure == 'S1':
            gf.figS1()
        elif figure == 'S2':
            gf.figS2()
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
        gd.stratifyData(['EMPclosed','HMP', 'MGRAST'], remove_obs = 0)
        gd.stratifyData(['EMPclosed','HMP', 'MGRAST'], remove_obs = 1)
        gd.stratifyData(['95','97', '99'], remove_obs = 0)
        gd.stratifyData1000(remove_obs = 0)
        gd.stratifyData1000(remove_obs = 1)
        gd.stratifyData1000(remove_obs = 0, seqSim = '95')
        gd.stratifyData1000(remove_obs = 0, seqSim = '97')
        gd.stratifyData1000(remove_obs = 0, seqSim = '99')
