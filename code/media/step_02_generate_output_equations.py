#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/19/20

"""generate_output_equations.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import sys

import cobra
import pandas as pd

sys.path.append('../Functions/')
import media

os.chdir('../../data/template_GEMs/')

# %% <IO>

iML1515 = cobra.io.load_json_model('iML1515_metacyc.json')
# iJO1366 = cobra.io.load_json_model('../archive/iML1515.json')
iYO844 = cobra.io.load_json_model('iYO844_metacyc.json')
iAF692 = cobra.io.load_json_model('iAF692_metacyc.json')
iMSU = cobra.io.read_sbml_model('iMSU_GEMs_supplementary-material-3.xml')
iMSU.id = 'iMSU'


def convert_to_df(exchange_rxns, transport_rxns):
    df_dic = {}
    for [k, v] in exchange_rxns:

        if k not in df_dic.keys():
            df_dic[k] = [[]] * 2
        a = df_dic[k][0].copy()
        a.append(v)
        df_dic[k][0] = a

    for [k, v] in transport_rxns:
        if k not in df_dic.keys():
            df_dic[k] = [[]] * 2
        a = df_dic[k][1].copy()
        a.append(v)
        df_dic[k][1] = a

    for k, v in df_dic.items():
        df_dic[k][0] = '; '.join(set(df_dic[k][0]))
        df_dic[k][1] = '; '.join(set(df_dic[k][1]))
    result_df = pd.DataFrame.from_dict(df_dic, orient='index',
                                       columns=['exchange_rxn', 'transport_rxn'])
    # result_df['met_id'] = result_df.index
    result_df = result_df.reset_index()
    return result_df


# %%
model = iML1515.copy()
exchange_rxns = media.get_exchange_rxns(model)
transport_rxns = media.get_transport_rxns(model)
iML1515_df = convert_to_df(exchange_rxns, transport_rxns)
iML1515_df.columns = ['id', 'exchange_rxns_iML1515', 'transport_rxns_iML1515']

model = iYO844.copy()
exchange_rxns = media.get_exchange_rxns(model)
transport_rxns = media.get_transport_rxns(model)
iYO844_df = convert_to_df(exchange_rxns, transport_rxns)
iYO844_df.columns = ['id', 'exchange_rxns_iYO844', 'transport_rxns_iYO844']

model = iAF692.copy()
exchange_rxns = media.get_exchange_rxns(model)
transport_rxns = media.get_transport_rxns(model)
iAF692_df = convert_to_df(exchange_rxns, transport_rxns)
iAF692_df.columns = ['id', 'exchange_rxns_iAF692', 'transport_rxns_iAF692']

result_df = iML1515_df.merge(iYO844_df, how='outer', on='id', suffixes=['', '_iYO844'])
result_df = result_df.merge(iAF692_df, how='outer', on='id', suffixes=['', '_iAF692'])

# %% <merge exchange>
result_df = result_df.fillna('')
result_df.loc[~(result_df['exchange_rxns_iML1515'] == ''), 'exchange_rxns'] = result_df.loc[
    ~(result_df['exchange_rxns_iML1515'] == ''), 'exchange_rxns_iML1515']
result_df.loc[(result_df['exchange_rxns'] == ''), 'exchange_rxns'] = result_df.loc[
    (result_df['exchange_rxns'] == ''), 'exchange_rxns_iYO844']
result_df.loc[(result_df['exchange_rxns'] == ''), 'exchange_rxns'] = result_df.loc[
    (result_df['exchange_rxns'] == ''), 'exchange_rxns_iAF692']

result_df.loc[~(result_df['transport_rxns_iYO844'] == ''), 'transport_rxns'] = result_df.loc[
    ~(result_df['transport_rxns_iYO844'] == ''), 'transport_rxns_iYO844']
result_df.loc[(result_df['transport_rxns'] == ''), 'transport_rxns'] = result_df.loc[
    (result_df['transport_rxns'] == ''), 'transport_rxns_iAF692']

# <manual checked>
result_df.loc[1, 'exchange_rxns'] = result_df.loc[1, 'exchange_rxns_iYO844']
result_df.loc[81, 'exchange_rxns'] = result_df.loc[81, 'exchange_rxns_iYO844']

result_df.to_csv('../media/exchange_and_transport_reactions.tsv', sep='\t')
# %%
iML1515_rxns = set([i.id for i in iML1515.reactions])
iYO844_rxns = set([i.id for i in iYO844.reactions])
iAF692_rxns = set([i.id for i in iAF692.reactions])

iML1515_mets = set([i.id for i in iML1515.metabolites])
iYO844_mets = set([i.id for i in iYO844.metabolites])
iAF692_mets = set([i.id for i in iAF692.metabolites])
result_df = result_df.fillna('')
add_ex_rxns = []
add_trans_rxns = []
for index_i in result_df.index:

    row_i = result_df.loc[index_i]
    if row_i['exchange_rxns'] != '':
        rxn_id = row_i['exchange_rxns']
        if ';' in rxn_id:
            print(rxn_id)
            rxn_id = rxn_id.split(';')[0]
        if rxn_id in iML1515_rxns:
            add_ex_rxns.append(iML1515.reactions.get_by_id(rxn_id))
        elif rxn_id in iYO844_rxns:
            add_ex_rxns.append(iYO844.reactions.get_by_id(rxn_id))
        elif rxn_id in iAF692_rxns:
            add_ex_rxns.append(iAF692.reactions.get_by_id(rxn_id))

    if row_i['transport_rxns'] != '':
        rxn_ids = row_i['transport_rxns']

        for rxn_id in rxn_ids.split('; '):
            if rxn_id in iYO844_rxns:
                add_trans_rxns.append(iYO844.reactions.get_by_id(rxn_id))
            elif rxn_id in iAF692_rxns:
                add_trans_rxns.append(iAF692.reactions.get_by_id(rxn_id))
            elif rxn_id in iML1515_rxns:
                add_trans_rxns.append(iML1515.reactions.get_by_id(rxn_id))

model = cobra.Model('exchange_and_transport_rxns')
for i in add_ex_rxns:
    i.notes['notes'] = 'exchange_reactions'
    try:
        model.add_reactions([i])
    except:
        print(i)
for i in add_trans_rxns:
    i.notes['notes'] = 'transport_reactions'
    try:
        model.add_reactions([i])
    except:
        print(i)
cobra.io.save_json_model(model, '../media/exchange_and_transport_reactions.json')
