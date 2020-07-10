#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 7/9/20

"""biomass_reactions.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import cobra
import os
import pandas as pd
import gemstool
import My_def

os.chdir('../../data/')

dirlist = os.listdir('draft_GEMs/')
model_name_list_xml = [i for i  in dirlist if i.endswith('xml')]
model_name_list_mat = [i for i  in dirlist if i.endswith('mat')]

# i = 0
# model = cobra.io.read_sbml_model(model_name_list_xml[1])
# model = cobra.io.load_matlab_model('draft_GEMs/'+model_name_list_mat[i])
# iML1515 = cobra.io.load_json_model('../archive/iML1515.json')
# iJO1366 = cobra.io.load_json_model('../archive/iML1515.json')
# iYO844 = cobra.io.load_json_model('../archive/iML1515.json')
# for i in iML1515.metabolites:
#     print(i.annotation)


# %%< mapping >
# biomass_df = pd.read_excel(r'biomass/Biomass_compare_summary.xlsx', sheet_name='summary',header=1)
# biomass_df_map = biomass_df[['seed_id','bigg_id','biocyc']]
# targetlist1, MNX_IDlist = gemstool.mapIDsViaMNXref.mapIDsViaMNXref('mets',biomass_df_map['bigg_id'].tolist(),'biggM','metacycM')
# biomass_df_map['metacyc_from_bigg'] = targetlist1
# targetlist2, MNX_IDlist = gemstool.mapIDsViaMNXref.mapIDsViaMNXref('mets',biomass_df_map['seed_id'].tolist(),'seedM','metacycM')
# biomass_df_map['metacyc_from_seed'] = targetlist2
# biomass_df_map = biomass_df_map.fillna('')
# biomass_df_map['biocyc'] = biomass_df_map.biocyc.apply(lambda x: x.replace('META:',''))

# biomass_df_map['metacyc'] = ''
# biomass_df_map['metacyc'][(biomass_df_map["metacyc_from_bigg"] == biomass_df_map["biocyc"])|
#                       (biomass_df_map["metacyc_from_seed"] == biomass_df_map["biocyc"])] = biomass_df_map["biocyc"]

# with pd.ExcelWriter('biomass/Biomass_compare_summary.xlsx',engine='openpyxl',
#                     mode='a') as writer:
#     biomass_df_map.to_excel(writer, sheet_name='Sheet_name_3')

# %% < get gram p and n df>
# biomass_df = pd.read_excel(r'biomass/Biomass_compare_summary.xlsx', sheet_name='summary',header=1)
# gram_n_df = biomass_df[['bigg_id','metacyc','coeff','observations']].dropna(subset=['coeff'])
# gram_p_df = biomass_df[['bigg_id','metacyc','@iYO844','observations']].dropna(subset=['@iYO844'])
# gram_n_df = gram_n_df.drop_duplicates(subset=['bigg_id']).fillna('')
# gram_p_df = gram_p_df.drop_duplicates(subset=['bigg_id']).fillna('')

# with pd.ExcelWriter('biomass/Biomass_compare_summary.xlsx',engine='openpyxl',
#                     mode='a') as writer:
#     gram_n_df.to_excel(writer, sheet_name='gram_n')
#     gram_p_df.to_excel(writer, sheet_name='gram_p')

# %% <get reaction>
gram_n_df = pd.read_excel(r'biomass/Biomass_compare_summary.xlsx', sheet_name='gram_n',header=1)
gram_p_df = pd.read_excel(r'biomass/Biomass_compare_summary.xlsx', sheet_name='gram_p',header=1)
# NOTE: some metabolites have no metacyc id
'''
print(gram_n_df[(gram_n_df['metacyc']=='')])
1     2fe2s         -0.000026
53      fe2         -0.006715
54      fe3         -0.007808
94      ni2         -0.000323
97    pe161         -0.075214
113     so4         -0.004338

print(gram_p_df[(gram_p_df['metacyc']=='')])
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

'''

gram_n_df = gram_n_df[(gram_n_df['metacyc']!='')]
gram_p_df = gram_p_df[(gram_p_df['metacyc']!='')]

gram_n_dic = dict(gram_n_df[['metacyc','coeff']].values.tolist())
gram_p_dic = dict(gram_p_df[['metacyc','@iYO844']].values.tolist())

# %% <add biomass to models>
biomass_n = cobra.Reaction('BIOMASS')
model = cobra.io.load_matlab_model('draft_GEMs/'+model_name_list_mat[0])
model.add_reaction(biomass_n)
# reaction_s = []
# reaction_p = []
# for k,v in gram_n_dic.items():
#     if v > 0:
#         reaction_p.append(str(v) + ' ' + k)
#     else:
#         reaction_s.append(str(-v) + ' ' + k)
# reaction = ' + '.join(reaction_s) + ' --> '+ ' + '.join(reaction_p)
# biomass_n.reaction = reaction

biomass_n.reaction = ' --> '
for k,v in gram_n_dic.items():
    # print(k)
    try:
        met = model.metabolites.get_by_id(k)
    except :
        print(met)
        met = cobra.Metabolite(k)

    biomass_n.add_metabolites({met:v})

biomass_p = cobra.Reaction('BIOMASS')
model = cobra.io.load_matlab_model('draft_GEMs/'+model_name_list_mat[0])
model.add_reaction(biomass_p)
# reaction_s = []
# reaction_p = []
# for k,v in gram_p_dic.items():
#     if v > 0:
#         reaction_p.append(str(v) + ' ' + k)
#     else:
#         reaction_s.append(str(-v) + ' ' + k)
# reaction = ' + '.join(reaction_s) + ' --> '+ ' + '.join(reaction_p)
# biomass_p.reaction = reaction

biomass_p.reaction = ' --> '
for k,v in gram_n_dic.items():
    # print(k)
    try:
        met = model.metabolites.get_by_id(k)
    except :
        print(met)
        met = cobra.Metabolite(k)
    biomass_p.add_metabolites({met:v})







