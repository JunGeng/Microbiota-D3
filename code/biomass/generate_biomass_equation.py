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

import cobra
import pandas as pd

os.chdir('../../data/biomass/')

input_file = 'biomass_constituents_id_map.tsv'
output_file_n = 'biomass_negative_reactions.txt'
output_file_p = 'biomass_positive_reactions.txt'
output_file_a = 'biomass_archaea_reactions.txt'


def get_biomass_model(dic_i, observations_list=False):
    '''
    :description : function to get a model contain biomass reactions from a dic
    :param :
        dic_i(dic) : {'met_id': [coefficient, 'pool]}
                    eg. {'10-FORMYL-THF': [-0.000367, 'cofactor'],
                     'AMP': [-0.0046700000000000005, 'cofactor'],
                     'CPD-12125': [-0.000266, 'cofactor'],
                     'CDP': [-0.000251, 'cofactor'],}
        observations_list like:
            eg. ['dna', 'rna', 'protein', 'lipid', 'carbohydrate', 'cofactor', 'ion', 'other']

    :returns: cobra model
    '''

    if observations_list == False:
        observations_list = ['dna', 'rna', 'protein', 'lipid', 'carbohydrate', 'cofactor', 'ion', 'other']
    observations_list = [i + '_c' for i in observations_list]

    biomass_model = cobra.Model('model')  # get a Model

    # %% <Biomass pool reactions>
    biomass_reactions = [cobra.Reaction('Biomass_' + i) for i in observations_list]
    [rxn.add_metabolites({cobra.Metabolite(met): 1}) for rxn, met in
     zip(biomass_reactions, observations_list)]

    biomass_model.add_reactions(biomass_reactions)

    for k, v in dic_i.items():
        # print(k)
        try:  # build metabolites
            met = biomass_model.metabolites.get_by_id(k)
        except:
            met = cobra.Metabolite(k)
        try:  # build reactions
            reaction_i = biomass_model.reactions.get_by_id('Biomass_' + v[1] + '_c')
            reaction_i.add_metabolites({met: v[0]})
        except:
            print(k, v, 'Error')

    # %% <Biomass reaction>
    Biomass = cobra.Reaction('Biomass')
    [Biomass.add_metabolites({biomass_model.metabolites.get_by_id(met): -1}) for met in observations_list]
    Biomass.add_metabolites({cobra.Metabolite('Biomass'): 1})
    biomass_model.add_reactions([Biomass])

    return biomass_model


# %% <add biomass to models>

biomass_df = pd.read_csv(input_file, sep='\t', header=0, index_col=0)
observations_list = ['dna', 'rna', 'protein', 'lipid', 'carbohydrate', 'cofactor', 'ion', 'other']

# %%<gram positive>
gram_p_df = biomass_df[~biomass_df['metacyc_id'].isna() & ~biomass_df['iYO844'].isna()]
gram_p_df = gram_p_df[['metacyc_id', 'iYO844', 'observations_modified']].sort_values(by='observations_modified')
gram_p_dic = gram_p_df.set_index('metacyc_id').T.to_dict('list')

model_p = get_biomass_model(gram_p_dic, observations_list)

# %%<gram negative>
gram_n_df = biomass_df[~biomass_df['metacyc_id'].isna() & ~biomass_df['iML1515'].isna()]
gram_n_df = gram_n_df[['metacyc_id', 'iML1515', 'observations_modified']].sort_values(by='observations_modified')
gram_n_dic = gram_n_df.set_index('metacyc_id').T.to_dict('list')

model_n = get_biomass_model(gram_n_dic, observations_list)

# %%<archaea>
archaea_df = biomass_df[~biomass_df['metacyc_id'].isna() & ~biomass_df['iAF692'].isna()]
archaea_df = archaea_df[['metacyc_id', 'iAF692', 'observations_modified']].sort_values(by='observations_modified')
archaea_dic = archaea_df.set_index('metacyc_id').T.to_dict('list')

model_a = get_biomass_model(archaea_dic, observations_list)

# %% <write model>
import gemstool

gemstool.io.gem2txt(model_n, output_file_n)  # write corba model to txt file
gemstool.io.gem2txt(model_p, output_file_p)
gemstool.io.gem2txt(model_a, output_file_a)

cobra.io.save_json_model(model_n, output_file_n.replace('txt', 'json'))  # write json model, easy import by python to
cobra.io.save_json_model(model_p, output_file_p.replace('txt', 'json'))
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
