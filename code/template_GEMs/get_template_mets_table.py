#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 7/9/20

"""biomass_constituents_id_map.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import sys

import cobra
import gemstool
import pandas as pd

sys.path.append('../Functions/')
import media

os.chdir('../../data/template_GEMs/')

# %% <IO>

iML1515 = cobra.io.load_json_model('iML1515.json')
# iJO1366 = cobra.io.load_json_model('../archive/iML1515.json')
iYO844 = cobra.io.load_json_model('iYO844.json')
iAF692 = cobra.io.load_json_model('iAF692.json')
iMSU = cobra.io.read_sbml_model('iMSU_GEMs_supplementary-material-3.xml')
iMSU.id = 'iMSU'

biomass_map_df = pd.read_csv('../biomass/biomass_constituents_id_map_.tsv', sep='\t', header=0, index_col=0)

output_mets = 'template_mets_id_map.tsv'
output_rxns = 'template_rxns_id_map.tsv'

exchange_rxns_list = media.get_exchange_rxns(iML1515) + media.get_exchange_rxns(iYO844) + media.get_exchange_rxns(
    iAF692)
exchange_rxns_list = set([i[1] for i in exchange_rxns_list])

# %%< four template model metabolites >
# mets_table = pd.DataFrame(columns=['id'])
# for model in [iML1515, iYO844, iAF692]:
#     # model = iYO844.copy()
#
#     # {'bigg.metabolite': ['2aobut'],
#     #  'biocyc': ['META:AMINO-OXOBUT'],
#     #  'chebi': ['CHEBI:6156',
#     #   'CHEBI:16944'],
#     #  'hmdb': ['HMDB06454'],
#     #  'inchi_key': ['SAUCHDKDCUROAO-VKHMYHEASA-N'],
#     #  'kegg.compound': ['C03508'],
#     #  'lipidmaps': ['LMFA01060172'],
#     #  'metanetx.chemical': ['MNXM114087'],
#     #  'reactome.compound': ['6798652'],
#     #  'sabiork': ['6672'],
#     #  'sbo': 'SBO:0000247',
#     #  'seed.compound': ['cpd02211']}
#
#     mets_table_dic = {'id': [],
#                  'name': [],
#                  'compartment': [],
#                  'metacyc_id': [],
#                  'kegg_id': [],
#                  'seed_id': [],
#                  'metanetx_id': []
#                  }
#
#     for met_i in model.metabolites:
#         annotation_i = met_i.annotation
#
#         mets_table_dic['id'].append(met_i.id)
#         mets_table_dic['name'].append(met_i.name)
#         mets_table_dic['compartment'].append(met_i.compartment)
#
#         biocyc_i = ''
#         keeg_i = ''
#         seed_i = ''
#         metanetx_id = ''
#         if 'biocyc' in annotation_i.keys():
#             biocyc_i = annotation_i['biocyc']
#         if 'kegg.compound' in annotation_i.keys():
#             keeg_i = annotation_i['kegg.compound']
#         if 'seed.compound' in annotation_i.keys():
#             seed_i = annotation_i['seed.compound']
#         if 'metanetx.chemical' in annotation_i.keys():
#             metanetx_id = annotation_i['metanetx.chemical']
#
#         mets_table_dic['metacyc_id'].append(biocyc_i)
#         mets_table_dic['kegg_id'].append(keeg_i)
#         mets_table_dic['seed_id'].append(seed_i)
#         mets_table_dic['metanetx_id'].append(metanetx_id)
#     df = pd.DataFrame.from_dict(mets_table_dic)
#     df[model.id] = model.id
#     mets_table = mets_table.merge(df, how='outer', on=['id'], suffixes=('', '_' + model.id))
#
# mets_table['name_iML1515'] = mets_table['name']
# mets_table['metacyc_id_iML1515'] = mets_table['metacyc_id']
# mets_table['kegg_id_iML1515'] = mets_table['kegg_id']
# mets_table['seed_id_iML1515'] = mets_table['seed_id']
# mets_table['metanetx_id_iML1515'] = mets_table['metanetx_id']
# mets_table = mets_table.fillna('')
# for index_i in mets_table.index:
#     name_i = set(mets_table.loc[index_i, ['name_iML1515', 'name_iYO844', 'name_iAF692']]) - {''}
#     name_i = '; '.join(name_i)
#
#     compartment_i = set(mets_table.loc[index_i, ['compartment', 'compartment_iYO844', 'compartment_iAF692']]) - {''}
#     compartment_i = '; '.join(compartment_i)
#
#     metacyc_i = set(
#         [ii for i in mets_table.loc[index_i, ['metacyc_id_iML1515', 'metacyc_id_iYO844', 'metacyc_id_iAF692']] for ii
#          in i if type(i) == list])
#     if len(metacyc_i) > 1:
#         # print(metacyc_i)
#         # print(mets_table.loc[index_i, ['metacyc_id_iML1515']][0])
#         metacyc_i = metacyc_i & set(mets_table.loc[index_i, ['metacyc_id_iML1515']][0])
#
#     metacyc_i = '; '.join(metacyc_i).replace('META:', '')
#
#     kegg_i = set(
#         [ii for i in mets_table.loc[index_i, ['kegg_id_iML1515', 'kegg_id_iYO844', 'kegg_id_iAF692']] for ii in i if
#          type(i) == list])
#     kegg_i = '; '.join(set(kegg_i))
#
#     seed_i = set(
#         [ii for i in mets_table.loc[index_i, ['seed_id_iML1515', 'seed_id_iYO844', 'seed_id_iAF692']] for ii in i if
#          type(i) == list])
#     seed_i = '; '.join(set(seed_i))
#
#     metanetx_i = set(
#         [ii for i in mets_table.loc[index_i, ['metanetx_id_iML1515', 'metanetx_id_iYO844', 'metanetx_id_iAF692']]
#          for ii in i if type(i) == list])
#     metanetx_i = '; '.join(set(metanetx_i))
#
#     mets_table.loc[index_i, ['metacyc_id']] = metacyc_i
#     mets_table.loc[index_i, ['kegg_id']] = kegg_i
#     mets_table.loc[index_i, ['seed_id']] = seed_i
#     mets_table.loc[index_i, ['metanetx_id']] = metanetx_i
#
# mets_table = mets_table.sort_values(by='id')
# mets_table['id_c'] = mets_table['id']
# mets_table['id'] = mets_table['id'].str.replace('_c', '')
# mets_table['id'] = mets_table['id'].str.replace('_e', '')
# mets_table['id'] = mets_table['id'].str.replace('_p', '')
#
# # %% <map id>
#
# targetlist1, MNX_IDlist = gemstool.mapIDsViaMNXref.mapIDsViaMNXref('mets', mets_table['id'].tolist(), 'biggM',
#                                                                    'metacycM')
# mets_table['metacyc_from_bigg'] = ['; '.join(set(i)) if type(i) == list else i for i in targetlist1]
#
# # merge biomass table
# mets_table = mets_table.merge(biomass_map_df[['metacyc_id','bigg_id']],how = 'left', left_on='id',right_on='bigg_id',suffixes = ['','_biomass_df'])
#
# mets_table['metacyc_id'] = mets_table[['metacyc_id', 'metacyc_from_bigg']].apply(lambda x: x[1] if x[0] == '' else x[0],
#                                                                                  axis=1)
# mets_table[['metacyc_id', 'metacyc_id_biomass_df']] = mets_table[['metacyc_id', 'metacyc_id_biomass_df']].fillna('')
# mets_table['metacyc_id'] = mets_table[['metacyc_id', 'metacyc_id_biomass_df']].apply(lambda x: x[0] if x[1] == '' else x[1],
#                                                                                  axis=1)
#
# columns = ['id_c', 'id', 'compartment', 'name', 'metacyc_id', 'kegg_id', 'seed_id',
#            'metanetx_id', 'metacyc_from_bigg',
#            'iML1515', 'iYO844', 'iAF692',
#            'metacyc_id_iML1515', 'kegg_id_iML1515', 'seed_id_iML1515', 'metanetx_id_iML1515',
#            'metacyc_id_iYO844', 'kegg_id_iYO844', 'seed_id_iYO844', 'metanetx_id_iYO844',
#            'metacyc_id_iAF692', 'kegg_id_iAF692', 'seed_id_iAF692', 'metanetx_id_iAF692',
#            ]
# mets_table = mets_table[columns]
# # a = mets_table[mets_table.duplicated(subset=['id_c'])]
# mets_table = mets_table.drop_duplicates(subset=['id_c', 'id', 'compartment', 'name', 'metacyc_id', ])
# mets_table.to_csv(output_mets, sep='\t', index=True)

# %%< four template model reactions >

rxns_table = pd.DataFrame(columns=['id'])
for model in [iML1515, iYO844, iAF692]:
    # model = iYO844.copy()

    # {'bigg.reaction': ['ASPCT'],
    #  'biocyc': ['META:ASPCARBTRANS-RXN'],
    #  'ec-code': ['2.1.3.2'],
    #  'kegg.reaction': ['R01397'],
    #  'metanetx.reaction': ['MNXR96080'],
    #  'reactome.reaction': ['R-DDI-73573',
    #   'R-XTR-73573'],
    #  'rhea': ['20015', '20014', '20016', '20013'],
    #  'sabiork': ['398'],
    #  'sbo': 'SBO:0000176',
    #  'seed.reaction': ['rxn01018']}

    rxns_table_dic = {'id': [],
                      'name': [],
                      'metacyc_id': [],
                      'kegg_id': [],
                      'seed_id': [],
                      'metanetx_id': []
                      }

    for rea_i in model.reactions:
        annotation_i = rea_i.annotation

        rxns_table_dic['id'].append(rea_i.id)
        rxns_table_dic['name'].append(rea_i.name)

        biocyc_i = ''
        keeg_i = ''
        seed_i = ''
        metanetx_id = ''
        if 'biocyc' in annotation_i.keys():
            biocyc_i = annotation_i['biocyc']
        if 'kegg.reaction' in annotation_i.keys():
            keeg_i = annotation_i['kegg.reaction']
        if 'seed.reaction' in annotation_i.keys():
            seed_i = annotation_i['seed.reaction']
        if 'metanetx.reaction' in annotation_i.keys():
            metanetx_id = annotation_i['metanetx.reaction']

        rxns_table_dic['metacyc_id'].append(biocyc_i)
        rxns_table_dic['kegg_id'].append(keeg_i)
        rxns_table_dic['seed_id'].append(seed_i)
        rxns_table_dic['metanetx_id'].append(metanetx_id)
    df = pd.DataFrame.from_dict(rxns_table_dic)
    df[model.id] = model.id
    rxns_table = rxns_table.merge(df, how='outer', on=['id'], suffixes=('', '_' + model.id))

rxns_table['name_iML1515'] = rxns_table['name']
rxns_table['metacyc_id_iML1515'] = rxns_table['metacyc_id']
rxns_table['kegg_id_iML1515'] = rxns_table['kegg_id']
rxns_table['seed_id_iML1515'] = rxns_table['seed_id']
rxns_table['metanetx_id_iML1515'] = rxns_table['metanetx_id']
rxns_table = rxns_table.fillna('')

for index_i in rxns_table.index:
    name_i = set(rxns_table.loc[index_i, ['name_iML1515', 'name_iYO844', 'name_iAF692']]) - {''}
    name_i = '; '.join(name_i)

    metacyc_i = set(
        [ii for i in rxns_table.loc[index_i, ['metacyc_id_iML1515', 'metacyc_id_iYO844', 'metacyc_id_iAF692']] for ii
         in i if type(i) == list])
    metacyc_i = '; '.join(set(metacyc_i)).replace('META:', '')

    kegg_i = set(
        [ii for i in rxns_table.loc[index_i, ['kegg_id_iML1515', 'kegg_id_iYO844', 'kegg_id_iAF692']] for ii in i if
         type(i) == list])
    kegg_i = '; '.join(set(kegg_i))

    seed_i = set(
        [ii for i in rxns_table.loc[index_i, ['seed_id_iML1515', 'seed_id_iYO844', 'seed_id_iAF692']] for ii in i if
         type(i) == list])
    seed_i = '; '.join(set(seed_i))

    metanetx_i = set(
        [ii for i in rxns_table.loc[index_i, ['metanetx_id_iML1515', 'metanetx_id_iYO844', 'metanetx_id_iAF692']]
         for ii in i if type(i) == list])
    metanetx_i = '; '.join(set(metanetx_i))

    rxns_table.loc[index_i, ['metacyc_id']] = metacyc_i
    rxns_table.loc[index_i, ['kegg_id']] = kegg_i
    rxns_table.loc[index_i, ['seed_id']] = seed_i
    rxns_table.loc[index_i, ['metanetx_id']] = metanetx_i

rxns_table = rxns_table.sort_values(by='id')

# %% <map id>

targetlist1, MNX_IDlist = gemstool.mapIDsViaMNXref.mapIDsViaMNXref('rxns', rxns_table['id'].tolist(), 'biggR',
                                                                   'metacycR')
rxns_table['metacyc_from_bigg'] = ['; '.join(set(i)) if type(i) == list else i for i in targetlist1]

rxns_table['metacyc_id'] = rxns_table[['metacyc_id', 'metacyc_from_bigg']].apply(lambda x: x[1] if x[0] == '' else x[0],
                                                                                 axis=1)

columns = ['id', 'name', 'metacyc_id', 'kegg_id', 'seed_id',
           'metanetx_id', 'metacyc_from_bigg',
           'iML1515', 'iYO844', 'iAF692',
           'metacyc_id_iML1515', 'kegg_id_iML1515', 'seed_id_iML1515', 'metanetx_id_iML1515',
           'metacyc_id_iYO844', 'kegg_id_iYO844', 'seed_id_iYO844', 'metanetx_id_iYO844',
           'metacyc_id_iAF692', 'kegg_id_iAF692', 'seed_id_iAF692', 'metanetx_id_iAF692',
           ]
rxns_table.loc[rxns_table['id'].isin(exchange_rxns_list), 'metacyc_id'] = ''
rxns_table = rxns_table[columns]
rxns_table.to_csv(output_rxns, sep='\t', index=True)
