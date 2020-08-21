#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/4/20

"""generate_biomass_quation.py
:description : get a models contain media reactions from a table
:description : function to get a model contain media exchange and transport reactions from a dic
:param :
    media_i_dic(dic) : coloums :Constituents,coefficient
                eg.  Constituents       LAB  coefficient unit
                0    chloride 51.607mM    51.607000   mM
:returns: cobra model
:rtype: txt and json

"""

import os

import cobra
import pandas as pd
import sys
sys.path.append('../Functions/')
import media




# %% <IO>
os.chdir('../../data/media/')
input_file = 'metacyc_growth_media_composition_df.tsv'
output_file = 'media_model_m2.json'
ex_and_trans_reactions = pd.read_csv('exchange_and_transport_reactions.tsv', sep='\t', index_col=0, )
media_df = pd.read_csv(input_file, sep='\t', header=0, index_col=0)
model = cobra.io.load_json_model('exchange_and_transport_reactions.json')
# %% <add process table>

media_df.columns = [i.replace(' ', '') for i in media_df.columns[0:4]] + [i for i in media_df.columns[4:]]

media_df = media_df.merge(ex_and_trans_reactions[['id', 'exchange_rxns', 'transport_rxns']], how='left',
                          left_on='metacyc_id', right_on='id')
media_df.to_csv('metacyc_growth_media_composition_df_trimed.tsv', sep='\t')

media_df['Constituents'] = media_df['metacyc_id']
media_name = 'Concentration_' + 'LAB'
media_m2_df = media_df[
    ['Constituents', media_name, 'exchange_rxns', 'transport_rxns']]  # select M2 columns

media_i_df = media_m2_df.copy()
media_i_df.columns = ['Constituents', 'Concentration', 'exchange_rxns', 'transport_rxns']
media_i_dic = media.get_media_dic_form_df(media_i_df, media_name)
model_i = model.copy()
media_model = media.changeout_media_from_dic(model_i, media_i_dic)

# %% <write model>
cobra.io.save_json_model(model_i, output_file)  # write json model, easy import by python to
