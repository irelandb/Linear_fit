#!/usr/bin/env python
from __future__ import division
import pandas as pd
import sortseq.profile_counts as profile_counts
import numpy as np
import scipy as sp
import sortseq.utils as utils

#create_dfs
names = ['O1','O2']
dfs = []
wts = []
models = []
seq_dict,inv_dict = utils.choose_dict('dna')
val_col = ['val_A','val_C','val_G','val_T']
for i,n in enumerate(names):
    models.append(pd.io.parsers.read_csv('/media/bill/9A1409B614099683/Nathan_data/Lac102016/results/results/' + n + '_wt_operator_MCMC_3',delim_whitespace=True))
    models[i] = np.transpose(np.array(models[i][val_col]))
    #dfs.append(pd.io.parsers.read_csv('../' + n + '_Op_only.txt',delim_whitespace=True))
    #wts.append(''.join(profile_counts.main(dfs[i],'dna',return_wtseq=True,start=54,end=74,bin_k=1)))
    #print wts[i]
    #wts[i] = utils.seq2mat(wts[i],seq_dict)
wts = ['AATTGTGAGCGGATAACAATT','AAATGTGAGCGAGTAACAACC','GGCAGTGAGCGCAACGCAATT']
wts = [utils.seq2mat(wts[i],seq_dict) for i in range(len(wts))]
wt_E = [-15.3,-13.9,-10.0]
#now for each dataset, fix constants by fitting to own data set and one additional
shift = []
for i,n in enumerate(names):
    #fix gauge of each model so that wt_E = 0 
    models[i] = models[i] - np.sum(models[i]*wts[i],axis=0)
    model_df = pd.DataFrame(np.transpose(models[i]))
    model_df.columns = val_col
    pd.set_option('max_colwidth',int(1e8))
    model_df.to_string(
        open(n + '1147_fixed_linear_all','w'), index=False,col_space=10,float_format=utils.format_string)
    shift.append(np.sum(models[i]*wts[i]) - wt_E[i])
    
    #now fix multiplicative shift by using each of the other data sets
    for z,n2 in enumerate(names):
        print ''.join(['model ', n, ' fit using model ', n2])
        if z != i:
            #calc additive shift
            print (np.sum(models[i]*wts[i]))
            C1 = (wt_E[z] - wt_E[i])/((np.sum(models[i]*wts[z]) - (np.sum(models[i]*wts[i]))))
            print C1
            print shift[i]
            for q,n3 in enumerate(names):
                if q != z and q != i:
                    print n3 + ' energy = ' + str(np.sum(wts[q]*models[i])*C1 + wt_E[i]) 



