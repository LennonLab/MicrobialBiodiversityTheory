import pandas as pd

otu_count = pd.read_table('/Users/WRShoemaker/Box Sync/METE/HMP/Mothur/V35/hmp1.v35.hq.otu.counts.txt', sep='\t', index_col=False)
sparsesbys = pd.melt(otu_count, id_vars=['collection'])
sparsesbys = sparsesbys[sparsesbys['value'] > 0]
sparsesbys.columns = ['Sample', 'OTU', 'Count']
sparsesbys.to_csv("HMPsparseSbyS.txt", sep='\t', index=False)
