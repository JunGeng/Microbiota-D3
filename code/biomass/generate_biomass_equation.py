#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/4/20

"""generate_biomass_quation.py
:description : get a models contain biomass reactions from a table
:param : 
:returns: model
:rtype: txt and json
"""

import os
import sys

import cobra
import pandas as pd

sys.path.append('../Functions/')
import biomass

os.chdir('../../data/biomass/')

# %% <IO>
input_file = 'biomass_constituents_id_map_.tsv'
output_file_n = 'biomass_negative_reactions.txt'
output_file_p_iYO844 = 'biomass_positive_reactions_iYO844.txt'
output_file_p_iSMU = 'biomass_positive_reactions_iSMU.txt'
output_file_a = 'biomass_archaea_reactions.txt'

# %% <add biomass to models>

biomass_df = pd.read_csv(input_file, sep='\t', header=0, index_col=0)
observations_list = ['dna', 'rna', 'protein', 'lipid', 'carbohydrate', 'cofactor', 'ion', 'other']

# %%<gram iYO844 positive>
gram_p_df = biomass_df[~biomass_df['metacyc_id'].isna() & ~biomass_df['iYO844'].isna()]
gram_p_df = gram_p_df[['metacyc_id', 'iYO844', 'observations_modified']].sort_values(by='observations_modified')
gram_p_dic = gram_p_df.set_index('metacyc_id').T.to_dict('list')

model_p_iYO844 = biomass.update_biomass_from_dic(gram_p_dic, observations_list)

# %%<gram iSMU positive>
gram_p_df = biomass_df[~biomass_df['metacyc_id'].isna() & ~biomass_df['iSMU'].isna()]
gram_p_df = gram_p_df[['metacyc_id', 'iSMU', 'observations_modified']].sort_values(by='observations_modified')
gram_p_dic = gram_p_df.set_index('metacyc_id').T.to_dict('list')

model_p_iSMU = biomass.changeout_biomass_from_dic(gram_p_dic, observations_list)

# %%<gram negative>
gram_n_df = biomass_df[~biomass_df['metacyc_id'].isna() & ~biomass_df['iML1515'].isna()]
gram_n_df = gram_n_df[['metacyc_id', 'iML1515', 'observations_modified']].sort_values(by='observations_modified')
gram_n_dic = gram_n_df.set_index('metacyc_id').T.to_dict('list')

model_n = biomass.update_biomass_from_dic(gram_n_dic, observations_list)

# %%<archaea>
archaea_df = biomass_df[~biomass_df['metacyc_id'].isna() & ~biomass_df['iAF692'].isna()]
archaea_df = archaea_df[['metacyc_id', 'iAF692', 'observations_modified']].sort_values(by='observations_modified')
archaea_dic = archaea_df.set_index('metacyc_id').T.to_dict('list')

model_a = biomass.update_biomass_from_dic(archaea_dic, observations_list)

# %% <write model>
import gemstool

gemstool.io.gem2txt(model_n, output_file_n)  # write corba model to txt file
gemstool.io.gem2txt(model_p_iYO844, output_file_p_iYO844)
gemstool.io.gem2txt(model_p_iSMU, output_file_p_iSMU)
gemstool.io.gem2txt(model_a, output_file_a)

cobra.io.save_json_model(model_n, output_file_n.replace('txt', 'json'))  # write json model, easy import by python to
cobra.io.save_json_model(model_p_iYO844, output_file_p_iYO844.replace('txt', 'json'))
cobra.io.save_json_model(model_p_iSMU, output_file_p_iSMU.replace('txt', 'json'))
cobra.io.save_json_model(model_a, output_file_a.replace('txt', 'json'))

# #
# # gram_n_df = gram_n_df[(gram_n_df['metacyc']!='')]
# # gram_p_df = gram_p_df[(gram_p_df['metacyc']!='')]
# #
# # gram_n_dic = dict(gram_n_df[['metacyc','coeff']].values.tolist())
# # gram_p_dic = dict(gram_p_df[['metacyc','@iYO844']].values.tolist())
# biomass_n = cobra.Reaction('BIOMASS')
# model = cobra.io.load_matlab_model('draft_GEMs/'+model_name_list_mat[0])
# model.add_reaction(biomass_n)
# # reaction_s = []
# # reaction_p = []
# # for k,v in gram_n_dic.items():
# #     if v > 0:
# #         reaction_p.append(str(v) + ' ' + k)
# #     else:
# #         reaction_s.append(str(-v) + ' ' + k)
# # reaction = ' + '.join(reaction_s) + ' --> '+ ' + '.join(reaction_p)
# # biomass_n.reaction = reaction
