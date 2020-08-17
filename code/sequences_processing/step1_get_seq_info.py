#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/13/20

"""step1_get_seq_info.py
:description : script to get the seq info
:param : input file: assembly_summary_refseq and species name list
:returns: output files 'species_info.txt' and 'species_summary.txt'
:rtype: txt
"""

import os

import pandas as pd

os.chdir('../../data/sequences_processing/')

# load database summery
refseq_database = pd.read_csv('../initial_data/data_from_database/assembly_summary_refseq.txt', sep='\t', header=1)
refseq_database_2 = pd.read_csv('../initial_data/data_from_database/assembly_summary_refseq_historical.txt', sep='\t',
                                header=1)
name_list_df = pd.read_csv('species_name_list.csv', sep=',', header=0)

# %% select data in '../initial_data/data_from_database/assembly_summary_refseq.txt'
refseq_info_df = pd.DataFrame()
refseq_info_df_summery = name_list_df.copy()

for index in name_list_df.index:

    name_i = name_list_df.loc[index, 'Name_trimmed']
    print(name_i)
    name_i = name_i.replace('[', r'\[')
    name_i = name_i.replace(']', r'\]')
    refseq_info_i_df = pd.DataFrame()
    refseq_info_i_df = refseq_database[refseq_database.organism_name.str.match(name_i)]

    n_match = refseq_info_i_df.shape[0]

    if n_match == 0.0:
        # if not in assembly_summary_refseq find in assembly_summary_refseq_historical
        refseq_info_i_df = refseq_database_2[refseq_database_2.organism_name.str.match(name_i)]
        n_match = refseq_info_i_df.shape[0]
        print(n_match)
    species_taxid = set(refseq_info_i_df.species_taxid.values)
    if len(species_taxid) == 1:
        species_taxid = species_taxid.pop()

    name_list_df.loc[index, 'seq_counts'] = int(n_match)

    if n_match == 0:
        continue

    refseq_info_i_df = refseq_info_i_df.copy()
    refseq_info_i_df.loc[:, 'Index_initial'] = name_list_df.loc[index, 'Index_initial']
    refseq_info_i_df.loc[:, 'Name_trimmed'] = name_list_df.loc[index, 'Name_trimmed']
    refseq_info_i_df.loc[:, 'Name_initial'] = name_list_df.loc[index, 'Name_initial']
    refseq_info_i_df.loc[:, 'Note'] = name_list_df.loc[index, 'Note']
    refseq_info_i_df.loc[:, 'seq_counts'] = int(n_match)
    refseq_info_df = pd.concat([refseq_info_df, refseq_info_i_df])
    name_list_df.loc[index, 'species_taxid'] = str(species_taxid)

    # check match and contains
    # refseq_info_i_df = refseq_database[refseq_database.organism_name.str.contains(name_i)]
    # n_contains = refseq_info_i_df.shape[0]
    # if n_contains!=n_match:
    #     print(n_match,n_contains)

# %% output file
refseq_info_df = refseq_info_df.drop_duplicates()
refseq_info_df = refseq_info_df[['Index_initial', 'Name_initial', 'Name_trimmed',
                                'Note', 'seq_counts',
                                '# assembly_accession', 'bioproject', 'biosample', 'wgs_master',
                                'refseq_category', 'taxid', 'species_taxid', 'organism_name',
                                'infraspecific_name', 'isolate', 'version_status', 'assembly_level',
                                'release_type', 'genome_rep', 'seq_rel_date', 'asm_name', 'submitter',
                                'gbrs_paired_asm', 'paired_asm_comp', 'ftp_path',
                                'excluded_from_refseq', 'relation_to_type_material']]
refseq_info_df.to_csv('species_info.txt', sep='\t', index=False)
name_list_df.to_csv('species_summary.txt', sep='\t', index=False)
print('Done')
# %% download:
