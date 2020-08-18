#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/18/20

"""branch_match_refence_species.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import re

import pandas as pd

os.chdir('../../data/')

# %% <IO>
species_table = 'sequences_processing/species_final.txt'

reference_species_table = 'media/organisms.tab'

refseq_database = pd.read_csv('initial_data/data_from_database/assembly_summary_refseq.txt', sep='\t', header=1)

species_df = pd.read_csv(species_table, sep='\t', header=0)

ref_species_df = pd.read_csv(reference_species_table, sep='\t', header=0)

out_put = 'species.tsv'
# species_df = species_df.fillna('')
# ref_species_df = ref_species_df.fillna('')

# %% <matched species method 1>

taxid = set(species_df['taxid'].dropna()) - {''}
ref_taxid = set(ref_species_df['taxid'].dropna()) - {''}

print('matched taxid/ref_taxid/our_taxid:', '%d/%d/%d' % (len(taxid & ref_taxid), len(ref_taxid), len(taxid)))

species = set([re.sub('[\[\]]', '', i) for i in species_df['Name_trimmed']])
ref_species = set(ref_species_df['genus'] + ref_species_df['lineage_species'].str.split('.', expand=True)[1])

print('matched species/ref_species/our_species:',
      '%d/%d/%d' % (len(species & ref_species), len(ref_species), len(species)))
strains = set([re.sub('[\[\]]', '', i) for i in species_df['organism_name']]) - {''}
ref_strains = [re.sub('(\s\(.*\))', '', i) for i in ref_species_df['full_name']]
ref_strains = set([ii.split(' No')[0] for ii in ref_strains]) - {''}

print('matched strains/ref_strains/our_strains :',
      '%d/%d/%d' % (len(strains & ref_strains), len(ref_strains), len(strains)))



# %% <matched species method 2>

ref_species_df_extend = ref_species_df.merge(refseq_database, how='left', on=['taxid'])

species_taxid = set(species_df['species_taxid'].dropna()) - {''}
ref_species_taxid = set(ref_species_df_extend['species_taxid'].dropna()) - {''}
print('matched species_taxid/ref_species_taxid/our_species_taxid:',
      '%d/%d/%d' % (len(species_taxid & ref_species_taxid), len(ref_species_taxid), len(species_taxid)))

print('matched species:', species & ref_species)
print('matched strains:', strains & ref_strains)

species_df.loc[(species_df['taxid'].isin(taxid & ref_taxid)), 'taxid_matched'] = True
species_df.loc[(species_df['species_taxid'].isin(species_taxid & ref_species_taxid)), 'species_taxid_matched'] = True

species_df.to_csv(out_put, sep='\t', index=False)
