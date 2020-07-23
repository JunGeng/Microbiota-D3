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

import gemstool
import numpy as np
import pandas as pd

os.chdir('../../data/biomass/')

# %% get template
# dirlist = os.listdir('draft_GEMs/')

# iML1515 = cobra.io.load_json_model('../archive/iML1515.json')
# iJO1366 = cobra.io.load_json_model('../archive/iML1515.json')
# iYO844 = cobra.io.load_json_model('../archive/iML1515.json')
# for i in iML1515.metabolites:
#     print(i.annotation)
# %%< four template model biomass tabble >
carveme_db = pd.read_excel(r'Biomass_compare_summary.xlsx', sheet_name='from carveme')  # iYO844 and iAF692

iML1515 = pd.read_excel(r'Biomass_compare_summary.xlsx', sheet_name='from iML1515 Hao.W ')  # iML1515
iML1515[['bigg_id', 'iML1515']] = iML1515[['bigg.metabolite', 'coeff']]
iML1515['comp'] = iML1515['id'].str[-1]

biomass_summary_df = pd.read_excel(r'Biomass_compare_summary.xlsx', sheet_name='summary', header=1)  # iSMU
biomass_summary_df = biomass_summary_df[['seed_id', 'bigg_id', 'comp', 'metacyc']]
df = pd.read_excel(r'Biomass_Literatures_Templates_Combined_new_correct_upload.xlsx', sheet_name='Biomass References',
                   header=1)
iSMU_df = pd.read_excel(r'iMSU_Biomass_KEGGID.xlsx', sheet_name='Sheet1', header=0)
# biomass_summary_df['iSMU'] = df['iSMU']
# biomass_summary_df = biomass_summary_df.dropna(how='all', subset=['metacyc', 'iSMU'])

# merge tables
biomass_map_df = carveme_db.merge(iML1515, on=['bigg_id', 'comp'], how='outer', suffixes=['', '_y'])
biomass_map_df = biomass_map_df.merge(biomass_summary_df, how='outer', on=['bigg_id', 'comp', 'seed_id'],
                                      suffixes=['', '_y'])
biomass_map_df = biomass_map_df.merge(biomass_summary_df, how='outer', on=['bigg_id', 'comp', 'seed_id'],
                                      suffixes=['', '_y'])
biomass_map_df = biomass_map_df.merge(iSMU_df[['iMSU_Biomass_Name', 'iMSU_Biomass_KEGGID']], how='left', left_on='name',
                                      right_on='iMSU_Biomass_Name',
                                      suffixes=['', '_y'])
keegset = set(biomass_map_df['iMSU_Biomass_KEGGID'])
iML1515_keeg_list = list(biomass_map_df['kegg.compound'])
for index in iSMU_df.index:
    keeg_id = iSMU_df.loc[index,]['iMSU_Biomass_KEGGID']
    if keeg_id not in keegset and keeg_id in iML1515_keeg_list:
        index_2 = iML1515_keeg_list.index(keeg_id)
        biomass_map_df.loc[index, ['iMSU_Biomass_KEGGID']] = keeg_id
        print(keeg_id)

biomass_map_df = biomass_map_df.merge(iSMU_df, how='outer',
                                      left_on='iMSU_Biomass_KEGGID', right_on='iMSU_Biomass_KEGGID',
                                      suffixes=['', '_y'])
biomass_map_df.loc[pd.isnull(biomass_map_df['name']), 'name'] = biomass_map_df[
    pd.isnull(biomass_map_df['name']) ]['iMSU_Biomass_Name_y']

# %%< rename columns >
print(biomass_map_df.columns)
# Index(['bigg_id', 'comp', 'seed_id', '@bacteria', '@grampos', '@gramneg',
#        '@archaea', '@iAF1260', '@iJO1366', '@iYO844', '@iAF692', 'name',
#        'observations', 'id', 'name_y', 'class', 'coeff', 'bigg.metabolite',
#        'biocyc', 'chebi', 'hmdb', 'kegg.compound', 'metanetx.chemical',
#        'seed.compound', 'kegg.drug', 'lipidmaps', 'kegg.glycan', 'iML1515',
#        'metacyc', 'metacyc_y', 'iMSU_Biomass_Name', 'iMSU_Biomass_KEGGID',
#        'iMSU_Biomass_Name_y', 'Coef'],
#       dtype='object')

biomass_map_df = biomass_map_df[
    ['metacyc', 'bigg_id', 'iMSU_Biomass_KEGGID', 'comp', 'name', 'observations', 'seed_id', 'kegg.compound', 'biocyc',
     'iML1515', '@iYO844',
     '@iAF692', 'Coef']]
biomass_map_df.columns = [
    'metacyc_id', 'bigg_id', 'iMSU_Biomass_KEGGID', 'comp', 'name', 'observations', 'seed_id', 'kegg_id_from_iML1515',
    'metacyc_id_from_iML1515', 'iML1515', 'iYO844',
    'iAF692', 'iSMU']

biomass_map_df = biomass_map_df.dropna(how='all', subset=['iML1515', 'iYO844', 'iAF692', 'iSMU'])
biomass_map_df = biomass_map_df.fillna('')
biomass_map_df = biomass_map_df.drop_duplicates(keep='first')
# %% <map id>

targetlist1, MNX_IDlist = gemstool.mapIDsViaMNXref.mapIDsViaMNXref('mets', biomass_map_df['bigg_id'].tolist(), 'biggM',
                                                                   'metacycM')
biomass_map_df['metacyc_from_bigg'] = targetlist1
targetlist2, MNX_IDlist = gemstool.mapIDsViaMNXref.mapIDsViaMNXref('mets', biomass_map_df['seed_id'].tolist(), 'seedM',
                                                                   'metacycM')
biomass_map_df['metacyc_from_seed'] = targetlist2

biomass_map_df = biomass_map_df.fillna('')
biomass_map_df['metacyc_id_from_iML1515'] = biomass_map_df['metacyc_id_from_iML1515'].apply(
    lambda x: x.replace('META:', ''))

# biomass_map_df['metacyc_id'][biomass_map_df['metacyc_id']== '']
# biomass_map_df['metacyc'][(biomass_map_df['metacyc_id']== '')|
#     (biomass_map_df["metacyc_from_bigg"] in biomass_map_df["metacyc_id_from_iML1515"].str.contains(biomass_map_df["metacyc_from_bigg"])) |
#                           (biomass_df_map["metacyc_from_seed"] == biomass_df_map["biocyc"])] = biomass_df_map["biocyc"]

index = biomass_map_df.apply(lambda x: (
        (str(x.metacyc_from_seed) in str(x.metacyc_id_from_iML1515)) and (x.metacyc_id == '') and (
        x.metacyc_from_seed != '')), axis=1)
biomass_map_df['metacyc_id'][index] = biomass_map_df['metacyc_from_seed'][index]

index = biomass_map_df.apply(lambda x: (
        (str(x.metacyc_from_bigg) in str(x.metacyc_id_from_iML1515)) and (x.metacyc_id == '') and (
        x.metacyc_from_bigg != '')), axis=1)
biomass_map_df['metacyc_id'][index] = biomass_map_df['metacyc_from_bigg'][index]

biomass_map_df.to_csv('biomass_constituents_id_map_initial.tsv', sep='\t')

# %% <get reaction>
gram_n_df = biomass_map_df[['metacyc_id', 'bigg_id', 'comp', 'seed_id', 'iML1515']][biomass_map_df['iML1515'] != '']
gram_p_df = biomass_map_df[['metacyc_id', 'bigg_id', 'comp', 'seed_id', 'iYO844']][biomass_map_df['iYO844'] != '']
archaea_df = biomass_map_df[['metacyc_id', 'bigg_id', 'comp', 'seed_id', 'iAF692']][biomass_map_df['iAF692'] != '']
# TODO: some metabolites have no metacyc id
'''
print(gram_n_df[(gram_n_df['metacyc_id']=='')])
1     2fe2s         -0.000026
53      fe2         -0.006715
54      fe3         -0.007808
94      ni2         -0.000323
97    pe161         -0.075214
113     so4         -0.004338

print(gram_p_df[(gram_p_df['metacyc_id']=='')])
24       cdlp_BS         -0.000005
35      d12dg_BS         -0.000110
54           fe3         -0.003450
55           gdp         -0.000180
62   gtca1_45_BS         -0.003624
63   gtca2_45_BS         -0.002347
64   gtca3_45_BS         -0.001819
76   lipo1_24_BS         -0.000007
77   lipo2_24_BS         -0.000006
78   lipo3_24_BS         -0.000018
79   lipo4_24_BS         -0.000015
86          mql7         -0.000266
99    peptido_BS         -0.101817
100      pgly_BS         -0.000176
107    psetha_BS         -0.000560
115     t12dg_BS         -0.000066
116      tcam_BS         -0.003112

print(archaea_df[(archaea_df['metacyc_id']=='')])
   metacyc_id    bigg_id comp   seed_id    iAF692
3                3hdpgpe    c  cpd15817 -0.002698
4                3hdpgpg    c  cpd15818  -0.01012
6                3hdpgps    c  cpd15820  -0.01147
9              adocblhbi    c  cpd15834  -0.00474
40                 dpgpi    c  cpd15860 -0.002698
45                f420_2    c  cpd00649  -4.3e-05
46                f420_3    c  cpd15868  -2.6e-05
47                f420_4    c  cpd15869 -0.000374
48                f420_5    c  cpd15870 -0.000299
49                f420_6    c  cpd15871  -1.1e-05
50                f420_7    c  cpd15872    -1e-06
51                  f430    c              -0.002
56                gdpgpi    c  cpd15880  -0.01619
60              glycogen    c  cpd00155  -0.01543

'''

# %% <add biomass to models>
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
#
# biomass_n.reaction = ' --> '
# for k,v in gram_n_dic.items():
#     # print(k)
#     try:
#         met = model.metabolites.get_by_id(k)
#     except :
#         print(met)
#         met = cobra.Metabolite(k)
#
#     biomass_n.add_metabolites({met:v})
#
# biomass_p = cobra.Reaction('BIOMASS')
# model = cobra.io.load_matlab_model('draft_GEMs/'+model_name_list_mat[0])
# model.add_reaction(biomass_p)
# # reaction_s = []
# # reaction_p = []
# # for k,v in gram_p_dic.items():
# #     if v > 0:
# #         reaction_p.append(str(v) + ' ' + k)
# #     else:
# #         reaction_s.append(str(-v) + ' ' + k)
# # reaction = ' + '.join(reaction_s) + ' --> '+ ' + '.join(reaction_p)
# # biomass_p.reaction = reaction
#
# biomass_p.reaction = ' --> '
# for k,v in gram_n_dic.items():
#     # print(k)
#     try:
#         met = model.metabolites.get_by_id(k)
#     except :
#         print(met)
#         met = cobra.Metabolite(k)
#     biomass_p.add_metabolites({met:v})
