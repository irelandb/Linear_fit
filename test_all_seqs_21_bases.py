#!/usr/bin/env python
from __future__ import division
import pandas as pd
import sortseq.profile_counts as profile_counts
import numpy as np
import scipy as sp
import sortseq.utils as utils
import sortseq.simulate_library as simulate_library
import sortseq.simulate_evaluate as evaluate_model
import sortseq.Models as Models

operators_df = pd.io.parsers.read_csv('Operators_to_test',delim_whitespace=True)

def calculate_distance(s,wt=None):
    number_mutations = np.sum(np.array(list(s)) != np.array(list(wt)))
    return number_mutations
#create_dfs
names = ['O1','O2']
dfs = []
wts = []
models = []
nbr_model = []
seq_dict,inv_dict = utils.choose_dict('dna')
#model_mult = [-10.53396702,3.25400925,2.5959971]

#model_mult = [3.7509,5.50533,2.5959971]
#old wt energies
# O2 1147 model_mult = [-19.48026649,17.2614184283,2.5959971]
# 446 model_mult = [15.4650335591,-14.5317858596,2.5959971]
# wt model_mult = [26.8446767006,6.86213064582,2.5959971]
#new wt energies
# wt model_mult = [25.0550315872,6.40465526943,2.59]
# 446 model_mult = [14.4340313218,-13.5630001356,2.59]
model_mult = [-18.1815820574,16.1106571997,2.59]
nbr_mult = [-.07002259,0.0334153,0.04766945]
nbr_add = [-169.226309,1,-19.49848]
val_col = ['val_A','val_C','val_G','val_T']
wt_for_sim = ['AATTGTGAGCGGATAACAATT','AAATGTGAGCGAGTAACAACC','GGCAGTGAGCGCAACGCAATT','AATTGTGAGCGAGTAACAATT','AACAGTGAGCGCATCGCAATT']
for i,n in enumerate(names):
    models.append(pd.io.parsers.read_csv('/media/bill/9A1409B614099683/Nathan_data/Lac102016/results/results/' + n + '_wt_operator_MCMC_3',delim_whitespace=True))
    #nbr_model.append(pd.io.parsers.read_csv('~/Documents/energymatrix/Lac3282016/MCMC/' + n + '_operator_MCMC_NBR_2',delim_whitespace=True))
    models[i] = np.transpose(np.array(models[i][val_col]))
    
    #if i == 0:
    #    dfs.append(pd.io.parsers.read_csv('../O1_Op_only.txt',delim_whitespace=True))
    #else:
    #    dfs.append(pd.io.parsers.read_csv('../' + n + '_Op_only.txt',delim_whitespace=True))
    #wts.append(''.join(profile_counts.main(dfs[i],'dna',return_wtseq=True,start=54,end=75,bin_k=1)))
    #print wts[i]
    #wts[i] = utils.seq2mat(wts[i],seq_dict)
wts = ['AATTGTGAGCGGATAACAATT','AAATGTGAGCGAGTAACAACC','GGCAGTGAGCGCAACGCAATT']
wts = [utils.seq2mat(wts[i],seq_dict) for i in range(len(wts))]
models = [models[i] - np.sum(models[i]*wts[i],axis=0) for i in range(len(models))]
wt_E = [-15.3,-13.9,-10.0,-17.3]
mylib = simulate_library.main(wtseq=wt_for_sim[0],numseq=1000000)
for wt in wt_for_sim[1:]:
    temp_lib = simulate_library.main(wtseq=wt,numseq=1000000)
    mylib = mylib.append(temp_lib,ignore_index=True)

mylib['O1_dist'] = mylib['seq'].apply(calculate_distance,args=(wt_for_sim[0],))
mylib['O2_dist'] = mylib['seq'].apply(calculate_distance,args=(wt_for_sim[1],))
mylib['O3_dist'] = mylib['seq'].apply(calculate_distance,args=(wt_for_sim[2],))

#now for each dataset, fix constants by fitting to own data set and one additional
for i,n in enumerate(names):
    #fix gauge of each model so that wt_E = 0 
    model = Models.LinearModel(models[i],'dna')
    operators_df[n] = model.genexp(operators_df['seq'])*model_mult[i] + wt_E[i]
    #write fixed model
    pd.set_option('max_colwidth',int(1e8))
    operators_df[n].to_string(
        open(n + '_wt_operator_MCMC_fixed','w'), index=False,float_format=utils.format_string)
    #model_nbr = Models.NeighborModel(nbr_model[i],'dna',is_df=True)
    col_name = n +'_val'
    col_name_nbr = n+'_nbr_val'
    mylib[col_name] = model.genexp(mylib['seq'])*model_mult[i] + wt_E[i]
    #mylib[col_name_nbr] = (model_nbr.genexp(mylib['seq'])-model_nbr.genexp([wt_for_sim[i]])[0])*nbr_mult[i] +wt_E[i]
    #print model_nbr.genexp([wt_for_sim[i]])[0]
    
    #now fix multiplicative shift by using each of the other data sets
pd.set_option('max_colwidth',int(1e8))
operators_df.to_string(
        open('operators_eval_21_bases_Old_model','w'), index=False,col_space=10,float_format=utils.format_string)

pd.set_option('max_colwidth',int(1e8))
mylib.to_string(
        open('all_seqs_eval_21_bases_Old_model','w'), index=False,col_space=10,float_format=utils.format_string)     



