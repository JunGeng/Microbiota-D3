#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/5/20

"""step2_add_biomass_and_media.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os
import sys

import cobra
import pandas as pd

sys.path.append('../../code/biomass/')
import dna
import rna

os.chdir('../../data/draft_GEMs/')

# %% <io>
# model_list = os.listdir('draft_from_RAVEN_metacyc23_5')
# fna_list = os.listdir('../sequences_processing/genomic_sequences')
# gbff_list = os.listdir('../sequences_processing/gbff_sequences')

species_table = '../sequences_processing/species_final.txt'
gram_table = '../initial_data/mostAbundantSpecies.tsv'

biomass_model_list = ['biomass_negative_reactions.json',
                      'biomass_positive_reactions.json',
                      'biomass_archaea_reactions.json'
                      ]

media_model_list = ['media_model_m2.json']

output_file_dir = 'draft_add_biomass_and_media/'

# % <load species table> # get species name list
species_df = pd.read_csv(species_table, sep='\t')
gram_df = pd.read_csv(gram_table, sep='\t')

species_df['merge_index'] = species_df['Index_initial'].astype(int)

species_df = species_df.merge(gram_df[['id', 'bofTemplateType']], how='left', left_on=['merge_index'], right_on=['id'])

species_df = species_df.drop_duplicates(subset=['organism_name'], keep='first')  # NOTE:butyrate_producing_bacterium dup
species_dic = species_df[['Name_trimmed', 'organism_name', 'bofTemplateType']].set_index('organism_name').T.to_dict(
    'list')

name_list = list(species_dic.keys())

index = 0  # try one model first
name_i = name_list[index].replace(' ', '_')

# %% <load draft model i>
model_path_i = 'draft_from_RAVEN_metacyc23_5/' + name_i + '_Metacyc.mat'

model_i = cobra.io.load_matlab_model(model_path_i)

# %% <load biomass reactions type: model>

biomass_n_model = cobra.io.load_json_model('../biomass/' + biomass_model_list[0])
biomass_p_model = cobra.io.load_json_model('../biomass/' + biomass_model_list[1])
biomass_a_model = cobra.io.load_json_model('../biomass/' + biomass_model_list[2])

gram = species_dic[name_list[index]][1]

if gram == 'GramNegative':
    biomass_model = biomass_n_model
    DNA_WEIGHT_FRACTION = 0.032754653286106
    RNA_WEIGHT_FRACTION = 0.2173669926367441
elif gram == 'GramPositive':
    biomass_model = biomass_p_model
    DNA_WEIGHT_FRACTION = 0.024836038302619
    RNA_WEIGHT_FRACTION = 0.06834193130139973
elif gram == 'Archaea':
    biomass_model = biomass_a_model
    DNA_WEIGHT_FRACTION = 0.03365199989709999
    RNA_WEIGHT_FRACTION = 0.22444954697759986
else:
    print('Error')

# %% <load media reactions type: model>

media_model_path = '../media/' + media_model_list[0]
media_model = cobra.io.load_json_model(media_model_path)

# %% <generate new dna and rna coefficients based on seq>

fna_i = '../sequences_processing/genomic_sequences/' + name_i + '_genomic.fna'
gbff_i = '../sequences_processing/gbff_sequences/' + name_i + '_genomic.gbff'
DNA_coefficients = dna.generate_coefficients(fna_i, DNA_WEIGHT_FRACTION=1)
RNA_coefficients, _ratios = rna.generate_coefficients(gbff_i,
                                                      RNA_WEIGHT_FRACTION=1,
                                                      rRNA_WEIGHT_FRACTION=0.9,
                                                      tRNA_WEIGHT_FRACTION=0.05,
                                                      mRNA_WEIGHT_FRACTION=0.05,
                                                      identifier='locus_tag')
dna_reaction = biomass_model.reactions.get_by_id('Biomass_dna_c')
# dna_reaction = cobra.Reaction('Biomass_dna_c')
dna_reaction.reaction = ' --> '
dna_reaction.add_metabolites({biomass_model.metabolites.get_by_id('DATP'): DNA_coefficients['dATP'],
                              biomass_model.metabolites.get_by_id('TTP'): DNA_coefficients['dTTP'],
                              biomass_model.metabolites.get_by_id('DCTP'): DNA_coefficients['dCTP'],
                              biomass_model.metabolites.get_by_id('DGTP'): DNA_coefficients['dGTP'],
                              biomass_model.metabolites.get_by_id('PPI'): DNA_coefficients['ppi_c'],
                              })

rna_reaction = biomass_model.reactions.get_by_id('Biomass_dna_c')
rna_reaction.reaction = ' --> '
rna_reaction.add_metabolites({biomass_model.metabolites.get_by_id('ATP'): RNA_coefficients['atp_c'],
                              biomass_model.metabolites.get_by_id('UTP'): RNA_coefficients['utp_c'],
                              biomass_model.metabolites.get_by_id('CTP'): RNA_coefficients['ctp_c'],
                              biomass_model.metabolites.get_by_id('GTP'): RNA_coefficients['gtp_c'],
                              biomass_model.metabolites.get_by_id('PPI'): RNA_coefficients['ppi_c'],
                              })

# %% <merge draft model and biomass model>
model_i = model_i.merge(biomass_model)

# %% <add compartments>
model_i.id = name_i
model_i.compartments = {'c': 'cytosol',
                        'e': 'extracellular space'}

for met_i in model_i.metabolites:
    if not met_i.id.endswith('_c'):
        met_i.id = met_i.id + '_c'
        met_i.compartment = 'c'

# %% <merge draft model and biomass model>
model_i = model_i.merge(media_model)

# %% <write model>
cobra.io.save_json_model(model_i, output_file_dir + model_i.id + '.json')
cobra.io.save_matlab_model(model_i, output_file_dir + model_i.id + '.mat')
