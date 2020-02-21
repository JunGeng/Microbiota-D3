#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/21/20

"""step3_select_strain.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import pandas as pd

os.chdir('../../data/sequences_processing/')

# load data
refseq_info_df_initial = pd.read_csv('species_info.txt', sep='\t')
refseq_info_df = refseq_info_df_initial.copy()

# sort:
category_sorter = ['reference genome', 'representative genome']
type_sorter = ['assembly from type material']
assembly_level_sorter = ['Complete Genome', 'Chromosome', 'Scaffold', 'Contig']

refseq_info_df.refseq_category = refseq_info_df.refseq_category.astype("category")
refseq_info_df.refseq_category.cat.set_categories(category_sorter, inplace=True)

refseq_info_df.relation_to_type_material = refseq_info_df.relation_to_type_material.astype("category")
refseq_info_df.relation_to_type_material.cat.set_categories(type_sorter, inplace=True)

refseq_info_df.assembly_level = refseq_info_df.assembly_level.astype("category")
refseq_info_df.assembly_level.cat.set_categories(assembly_level_sorter, inplace=True)

refseq_info_df = refseq_info_df.sort_values(
    ['Index_initial', 'refseq_category', 'relation_to_type_material', 'assembly_level', 'seq_rel_date'],
    ascending=[True, True, True, True, False])

# remove duplicates, only keep the fiest one
refseq_info_df = refseq_info_df.drop_duplicates('Index_initial')

final_df = pd.DataFrame()
final_df[['# assembly_accession','Index_initial']] = refseq_info_df[['# assembly_accession','Index_initial']]
final_df = final_df.merge(refseq_info_df_initial,'inner', left_on=['# assembly_accession','Index_initial'],
                          right_on=['# assembly_accession','Index_initial'])
final_df = final_df[['Index_initial', 'Name_initial', 'Name_trimmed',
                     'Note', 'seq_counts',
                     '# assembly_accession', 'bioproject', 'biosample', 'wgs_master',
                     'refseq_category', 'taxid', 'species_taxid', 'organism_name',
                     'infraspecific_name', 'isolate', 'version_status', 'assembly_level',
                     'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter',
                     'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
                     'excluded_from_refseq', 'relation_to_type_material']]
final_df.to_csv('species_final.txt', sep='\t', index=False)

print('Done')

# %% <dataframe of different selections part>
# all_set = set(refseq_info_df['Index_initial'])
# # select unique:
# unique_df = refseq_info_df[(refseq_info_df['seq_counts']==1)]
# unique_set = set(unique_df['Index_initial'])
#
# # select 'refseq_category': 'reference genome' ,'representative genome
# category_df = refseq_info_df[refseq_info_df['refseq_category'].isin({'reference genome','representative genome'})]
# category_df = category_df.sort_values(by=['refseq_category'])
# category_set = set(category_df['Index_initial'])
#
# # select 'relation_to_type_material': 'assembly from type material '
# type_df = refseq_info_df[refseq_info_df['relation_to_type_material']=='assembly from type material']
# type_set = set(type_df['Index_initial'])
#
#
# # the rest: #5
# ## Complete Genome of rest strains #2
# rest_set = all_set-(unique_set|category_set|type_set)
# print(len(rest_set))
# rest_df = refseq_info_df[refseq_info_df['Index_initial'].isin(rest_set)]
#
# complete_df = rest_df[rest_df.assembly_level == 'Complete Genome']
# complete_set = set(complete_df['Index_initial'])
#
# rest2_set = rest_set - complete_set
# rest2_df =rest_df[rest_df['Index_initial'].isin(rest2_set)]
# ##
