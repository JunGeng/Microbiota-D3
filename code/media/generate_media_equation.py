#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/4/20

"""generate_biomass_quation.py
:description : get a models contain media reactions from a table
:param :
:returns: model
:rtype: txt and json
"""

import os

import cobra
import pandas as pd

os.chdir('../../data/media/')

input_file = 'metacyc_growth_media_composition_df.tsv'
output_file = 'media_model_m2.json'


def get_media_dic_form_df(media_i_df, media_name):
    '''
    :param :
        media_i_df (pandas dataframe):
            {'met_id': [coefficient,]}
                    eg. {'chloride ': [51.607],
                            'Ca2+  ': [0.27462000000000003],}
        media_name

    :returns: cobra model
    '''
    media_i_df.loc[:, media_name] = media_i_df[media_name].str.replace(' ', '')  # get dict
    media_i_df.loc[:, 'Constituents'] = media_i_df['Constituents'].str.replace(' ', '')
    media_i_df = media_i_df[media_i_df[media_name] != '']
    media_i_df = media_i_df.dropna(subset=[media_name])
    media_i_df.loc[:, 'coefficient'] = media_i_df[media_name].str.extract(r'([0-9\.]+)', expand=True)[0]
    media_i_df.loc[:, 'unit'] = media_i_df[media_name].str.extract(r'([µMnm]+)', expand=True)[0]

    media_i_df.loc[:, 'coefficient'] = media_i_df.apply(
        lambda x: float(x.coefficient) * 1e-6 if x['unit'] == 'nM' else float(x.coefficient), axis=1)
    media_i_df.loc[:, 'coefficient'] = media_i_df.apply(
        lambda x: float(x.coefficient) * 1e-3 if x.unit == 'µM' else float(x.coefficient), axis=1)

    media_i_dic = media_i_df[['Constituents', 'coefficient']].set_index('Constituents').T.to_dict('list')

    return media_i_dic


def get_media_model_from_dic(media_i_dic):
    '''
    :description : function to get a model contain media exchange and transport reactions from a dic
    :param :
        media_i_dic(dic) : coloums :Constituents,coefficient
                    eg.  Constituents       LAB  coefficient unit
                    0    chloride 51.607mM    51.607000   mM

    :returns: cobra model
    '''

    media_model = cobra.Model('model')  # get a Model

    # %% <exchange and transport reactions>

    for k, v in media_i_dic.items():
        # print(k)
        try:  # build metabolites
            met_c = media_model.metabolites.get_by_id(k + '_c')
            met_e = media_model.metabolites.get_by_id(k + '_e')
        except:
            met_c = cobra.Metabolite(k + '_c')
            met_e = cobra.Metabolite(k + '_e')
        try:  # build reactions
            teansport_reaction_i = cobra.Reaction('TRANS_' + k)
            teansport_reaction_i.add_metabolites({met_c: -1, met_e: 1})

            exchange_reaction_i = cobra.Reaction('Ex_' + k)
            exchange_reaction_i.add_metabolites({met_e: -1})
            exchange_reaction_i.lower_bound = -float(v[0])

            media_model.add_reactions([teansport_reaction_i, exchange_reaction_i])

        except:
            print(k, v, 'Error')

    return media_model


# %% <add process table>

media_df = pd.read_csv(input_file, sep='\t', header=2, index_col=0)
media_m2_df = media_df[['Constituents                      ', 'LAB']]  # select M2 columns
media_m2_df.columns = ['Constituents', 'LAB']
# media_m2_df.loc[:, 'Constituents'] = media_m2_df['Constituents'].str.replace(' ', '')
# media_m2_df.loc[:, 'LAB'] = media_m2_df['LAB'].str.replace(' ', '')

media_i_df = media_m2_df.copy()
media_name = 'LAB'

media_i_dic = get_media_dic_form_df(media_i_df, media_name)
model_i = get_media_model_from_dic(media_i_dic)

# %% <write model>

cobra.io.save_json_model(model_i, output_file)  # write json model, easy import by python to
