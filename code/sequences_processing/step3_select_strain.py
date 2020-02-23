#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/21/20

"""step3_select_strain.py
:description : script to choose the type strain
:param : input: species_info.txt
:returns: output: species_final.txt
:rtype: 
"""

import os
import re
import wget
import pandas as pd

os.chdir('../../data/sequences_processing/')

# load data
refseq_info_df_initial = pd.read_csv('species_info.txt', sep='\t')
refseq_info_df = refseq_info_df_initial.copy()

# %% < sort >
pass
# make category sorter for sort
category_sorter = ['reference genome', 'representative genome']
type_sorter = ['assembly from type material']
assembly_level_sorter = ['Complete Genome', 'Chromosome', 'Scaffold', 'Contig']

refseq_info_df.refseq_category = refseq_info_df.refseq_category.astype("category")
refseq_info_df.refseq_category.cat.set_categories(category_sorter, inplace=True)

refseq_info_df.relation_to_type_material = refseq_info_df.relation_to_type_material.astype("category")
refseq_info_df.relation_to_type_material.cat.set_categories(type_sorter, inplace=True)

refseq_info_df.assembly_level = refseq_info_df.assembly_level.astype("category")
refseq_info_df.assembly_level.cat.set_categories(assembly_level_sorter, inplace=True)

# sort: : category_sorter > type_sorter > assembly_level_sorter
refseq_info_df = refseq_info_df.sort_values(
    ['Index_initial', 'refseq_category', 'relation_to_type_material', 'assembly_level', 'seq_rel_date'],
    ascending=[True, True, True, True, False])

# remove duplicates, only keep the first one
refseq_info_df = refseq_info_df[refseq_info_df.ftp_path != 'na' ]
refseq_info_df = refseq_info_df.drop_duplicates('Index_initial')

# output file and resort column
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
# # detail data for each selection part
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

# %% < select seqs and copy them to new dir >
# os.chdir('data/sequences_processing/')
final_df_copy = pd.read_csv('species_final.txt', sep='\t')
all_seq_list = set(os.listdir('sequences/'))

for index in range(0, final_df_copy.shape[0]):
    seq_path = final_df_copy.iloc[index]['ftp_path']
    if seq_path == 'na':
        continue
    faa_name = seq_path.split('/')[-1] + '_protein.faa.gz'
    fna_name = seq_path.split('/')[-1] + '_genomic.fna.gz'
    gbff_name = seq_path.split('/')[-1] + '_genomic.gbff.gz'
    organism_name = final_df_copy.iloc[index]['organism_name']
    # name error
    organism_name = organism_name.replace(' ','_')
    organism_name = organism_name.replace('/','_')
    # organism_name_list = organism_name.split(' ')
    # infraspecific_name_list = final_df_copy.iloc[index]['infraspecific_name'].split('=')
    # short_name = organism_name_list[0][re.search('[A-Z,a-z]', organism_name_list[0]).start()] + '_' + '_'.join(
    #     organism_name_list[1:])

    for file_name in [gbff_name,faa_name,fna_name] :
        if file_name not in all_seq_list:
            file_url = seq_path + '/' + file_name
            wget.download(file_url, 'sequences/')

    os.system('cp sequences/' + fna_name + ' genomic_sequences/' + organism_name + '_genomic.fna.gz')
    os.system('cp sequences/' + faa_name + ' protein_sequences/' + organism_name + '_protein.faa.gz')
    print(index)
print('Done')