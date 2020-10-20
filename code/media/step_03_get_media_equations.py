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
import sys

import cobra
import pandas as pd

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
media_df.loc[media_df['Constituents'] == 'D-glucose', 'metacyc_id'] = 'Glucopyranose'
media_df.loc[media_df['Constituents'] == 'fructose', 'metacyc_id'] = 'fru'
media_df.loc[media_df['Constituents'] == 'lactose', 'metacyc_id'] = 'D-LACTATE'
media_df.loc[media_df['Constituents'] == 'maltose', 'metacyc_id'] = 'MAL'
media_df.loc[media_df['Constituents'] == 'adenine', 'metacyc_id'] = 'ade'

media_df.loc[media_df['Constituents'] == 'dihydrogenphosphate', 'metacyc_id'] = 'DIHYDROXYACETONE'
media_df.loc[media_df['Constituents'] == 'myo-inositol', 'metacyc_id'] = 'inost'
media_df.loc[media_df['Constituents'] == 'nitrate', 'metacyc_id'] = 'no3'
ni_media = ['M1: dGMM', 'M3: dGMM+LAB', 'M4: dGMM+LAB low M/V', 'M5: dGMM+LAB exclude SCFA',
            'M7: dGMM+LAB only monosacharides', 'M8: dGMM+LAB plus Mucin',
            'M9: dGMM+LAB only Mucin', 'M10: dGMM+LAB 10% aminoacids',
            'M11: dGMM+LAB excluding aromatic AA',
            ]
media_df.loc[media_df['Constituents'] == 'Ni2+', ni_media] = '990 nM'
media_df = media_df.merge(ex_and_trans_reactions[['id', 'exchange_rxns', 'transport_rxns']], how='left',
                          left_on='metacyc_id', right_on='id')

# TODO:Ni+

media_df.to_csv('metacyc_growth_media_composition_df_trimed.tsv', sep='\t')

media_df['Constituents'] = media_df['metacyc_id']
media_name = 'Concentration_' + 'LAB'
media_name = 'M2: ' + 'LAB'
media_m2_df = media_df[
    ['Constituents', media_name, 'exchange_rxns', 'transport_rxns']]  # select M2 columns

media_i_df = media_m2_df.copy()
media_i_df.columns = ['Constituents', 'Concentration', 'exchange_rxns', 'transport_rxns']
media_i_dic = media.get_media_dic_form_df(media_i_df, media_name)
model_i = model.copy()
media_model = media.update_media_from_dic(model_i, media_i_dic)

# %% <write model>
cobra.io.save_json_model(model_i, output_file)  # write json model, easy import by python to
import gemstool
gemstool.io.gem2txt(model_i,output_file.replace('json','txt'))