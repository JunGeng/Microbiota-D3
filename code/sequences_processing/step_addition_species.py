#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/13/20

"""step_addation_species.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import pandas as pd
import time
import wget

os.chdir('../../data/sequences_processing/')

#%% load database summery
refseq_database = pd.read_csv('../initial_data/data_from_database/assembly_summary_refseq.txt', sep='\t', header=1)
refseq_database_2 = pd.read_csv('../initial_data/data_from_database/assembly_summary_refseq_historical.txt', sep='\t',
                                header=1)

new_species_table = '../media/organisms.tab'
new_species_df = pd.read_csv(new_species_table, sep='\t', header=0)
new_species_df['Index_initial'] = new_species_df.index

new_refseq_info_df = new_species_df.merge(refseq_database, how='left', on='taxid', left_index=True)
new_refseq_info_df[new_refseq_info_df['organism_name'].isna()] = new_refseq_info_df[
    new_refseq_info_df['organism_name'].isna()].merge(refseq_database_2, how='left', on='taxid', )

# sort: : category_sorter > type_sorter > assembly_level_sorter
new_refseq_info_df = new_refseq_info_df.sort_values(
    ['Index_initial', 'refseq_category', 'relation_to_type_material', 'assembly_level', 'seq_rel_date'],
    ascending=[True, True, True, True, False])
# remove duplicates, only keep the first one
new_refseq_info_df = new_refseq_info_df[~new_refseq_info_df.ftp_path.isna()]
new_refseq_info_df = new_refseq_info_df.drop_duplicates('Index_initial')

new_refseq_info_df.to_csv('addition_species_final.tsv', sep='\t', index=False)



# %% download:

final_df_copy = pd.read_csv('addition_species_final.tsv', sep='\t')
all_seq_list = set(os.listdir('sequences/'))

for index in range(0, final_df_copy.shape[0]):
    seq_path = final_df_copy.iloc[index]['ftp_path']
    if seq_path == 'na':
        # print(1)
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

    os.system('cp sequences/' + fna_name + ' genomic_sequences/addition_' + organism_name + '_genomic.fna.gz')
    os.system('cp sequences/' + faa_name + ' protein_sequences/addition_' + organism_name + '_protein.faa.gz')
    os.system('cp sequences/' + gbff_name + ' gbff_sequences/addition_' + organism_name + '_genomic.gbff.gz')
    print(index)

os.system('gunzip gbff_sequences/*.gz')
os.system('gunzip protein_sequences/*.gz')
os.system('gunzip genomic_sequences/*.gz')

print('Done')

