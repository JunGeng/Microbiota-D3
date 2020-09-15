#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/14/20

"""branch_merge_species_table.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import pandas as pd

os.chdir('../../data/')

#%% <IO>
species_table = 'species.tsv'
addition_species_table = 'sequences_processing/addition_species_final.tsv'

experiment_db = pd.read_excel(r'media/NatureMicrobiology2018.3.514–522.xlsx', sheet_name='S4. Growth matrix')
species_df = pd.read_csv(species_table, sep='\t', header=0)
addition_species_df = pd.read_csv(addition_species_table, sep='\t', header=0)

#%% <merge addition species>
addition_species_df['Index_initial'] = 'ref_' + addition_species_df['Index_initial'].astype('str')
addition_species_df['ref_Index_initial'] = addition_species_df['Index_initial']

new_df = species_df.merge(addition_species_df[['species', 'full_name', 'lineage_species', 'genus', 'family', 'order',
                                               'class', 'phylum', 'assembly', 'agora', 'sort', 'nid', 'taxid',
                                               'nearest_kegg_org', 'ref_Index_initial']], how='left',
                          on=['taxid'])

new_df = new_df.append(addition_species_df[~addition_species_df['taxid'].isin(new_df['taxid'])])

#%% <merge experimet data>
new_df['designation in screen'] = new_df['species']
new_df = new_df.drop(['GMM', 'BHI++', 'WCA', 'mGAM', 'M1', 'M2', 'M3', 'M4', 'M5', 'M7', 'M8', 'M9',
                      'M10', 'M11', 'M13', 'M14', 'M15A', 'M15B', 'M16', ], axis=1)

new_df = new_df.merge(experiment_db, how='left', left_on='species', right_on='designation in screen')

#%% <get file name>
temp = new_df.loc[new_df['file_name'].isna(), 'organism_name']
temp = temp.str.replace(' ', '_')
temp = 'addition_' + temp.str.replace('/', '_')
new_df.loc[new_df['file_name'].isna(), 'file_name'] = temp

#%% <get gram_p/n>
gram_database = pd.read_csv('initial_data/data_from_database/genome_metadata', sep='\t', header=0)
gram_database.columns
a = gram_database[['taxon_id','gram_stain','comments']]

new_df = new_df.merge(a,how='left',left_on='taxid',right_on='taxon_id')
new_df = new_df.drop_duplicates(subset='taxid')

new_df['gram_stain'] = new_df['gram_stain'].str.replace('Positive','GramPositive')
new_df['gram_stain'] = new_df['gram_stain'].str.replace('+','GramPositive')
new_df['gram_stain'] = new_df['gram_stain'].str.replace('Negative','GramNegative')
new_df['gram_stain'] = new_df['gram_stain'].str.replace('-','GramNegative')
new_df['gram_stain==bofTemplateType'] = True
new_df['gram_stain==bofTemplateType'] = new_df['gram_stain']==new_df['bofTemplateType']

#%% <out put>
out_put = 'temp_species_6.tsv'
out_put = 'species_with_addition.tsv'
new_df.to_csv(out_put, sep='\t', index=False)
