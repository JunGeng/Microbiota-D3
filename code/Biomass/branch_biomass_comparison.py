#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 7/1/20

"""branch_biomass_comparison.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import cobra
import scipy
import os
import pandas as pd
import gemstool

os.chdir('../../data/branch_data')

file_list = os.listdir(".")

mat_list = [i for i in file_list if i.endswith('mat')]
xml_list = [i for i in file_list if i.endswith('xml')]

# mat_models = [cobra.io.load_matlab_model(i) for i in file_list if i.endswith('mat')]
# xml_models = [cobra.io.read_sbml_model(i) for i in file_list if i.endswith('xml')]

# summary = pd.DataFrame(columns=['id','name'])
#
# for i in mat_list:
#     try:
#         model_i = cobra.io.load_matlab_model(i)
#         gemstool.io.gem2txt(model_i,i.replace('mat','txt'))
#     except :
#         print(i)
#     biomass = {'id':[],model_i.id:[]}
#     biomass_ = model_i.reactions.get_by_id('biomass0').metabolites
#     for k,v in biomass_.items():
#         biomass['id'].append(k.id)
#         biomass['name'].append(k.name)
#         biomass[model_i.id].append(v)
#     modeli_i_df = pd.DataFrame(biomass)
#     summary = summary.merge(modeli_i_df,left_on=['id'], right_on=['id'],how='outer')
#
# summary.to_csv('mat_models_2.txt',sep = '\t')
# summary = pd.read_csv('mat_models_2.txt',sep = '\t',header=0,index_col=0)
#
# summary['seed_id'] = summary['id'].str.replace('[c0]','',regex=False)
#
# summary = summary.drop(columns=['name_x', 'name_y'])
#
#
# seed_db= pd.read_excel(r'Biomass_compare_summary.xlsx', sheet_name='from SEED')
# seed_db['seed_id'] = seed_db['Compound ID']
# summary_all = seed_db.merge(summary,left_on='seed_id', right_on='seed_id',how='outer')
#
#
# carveme_db= pd.read_excel(r'Biomass_compare_summary.xlsx', sheet_name='from carveme')
# iML1515 = pd.read_excel(r'Biomass_compare_summary.xlsx', sheet_name='from iML1515 Hao.W ')
# iML1515['bigg_id'] = iML1515['bigg.metabolite']
#
# bigg_db = carveme_db.merge(iML1515,left_on='bigg_id', right_on='bigg_id',how='outer')
#
#
# summary_all = bigg_db.merge(summary_all,left_on='seed_id', right_on='seed_id',how='outer')
#
# summary_all.to_csv('summary_all.txt',sep = '\t',index=False)


summary2 = pd.DataFrame(columns=['name'])

for i in xml_list:
    model_i = cobra.io.read_sbml_model(i)
    # gemstool.io.gem2txt(model_i,i.replace('xml','txt'))

    bio_boj = {'iBif452.V01.00.xml':'Biomass','iFap484.V01.00.xml':'Biomass',
               'iBth801 v1.00.xml':'biomass_red','iMsi385.xml':'r449',
               'iEre400 v1.00.xml':'r48'}

    biomass = {'id':[],'name':[],i:[],i+'notes':[],i+'annotation':[]}
    biomass_ = model_i.reactions.get_by_id(bio_boj[i]).metabolites
    for k,v in biomass_.items():
        biomass['id'].append(k.id)
        biomass['name'].append(k.name.split('_')[0])
        biomass[i+'notes'].append(k.notes)
        biomass[i+'annotation'].append(k.annotation)
        biomass[i].append(v)
    modeli_i_df = pd.DataFrame(biomass)
    summary2 = summary2.merge(modeli_i_df,left_on='name', right_on='name',how='outer')

summary2.to_csv('xml_models.txt',sep = '\t')



# carveme_db= pd.read_excel(r'Biomass_compare_summary.xlsx', sheet_name='model1')

