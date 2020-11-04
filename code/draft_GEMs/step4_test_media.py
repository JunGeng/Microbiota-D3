#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 10/28/20

"""step4_test_media.py
:description : script
:param : 
:returns: 
:rtype: 
"""

# import seaborn as sns
import os
import sys

import cobra
# import matplotlib.pyplot as plt
import pandas as pd

# sys.path.append('../../code/Functions/')
sys.path.append('../Functions/')
import media

os.chdir('../../data/')

# %% <IO>
species_table = 'species_with_addition.tsv'
# species_table = 'temp_species_with_addition_3.tsv'
species_df = pd.read_csv(species_table, sep='\t')
name_list = list(species_df['file_name'])
input_dir_1 = 'draft_GEMs/model_version_1.0/'
input_dir_2 = 'draft_GEMs/draft_from_RAVEN_metacyc23_5/'

species_df = species_df.drop_duplicates(subset=['file_name'], keep='first')  # NOTE:butyrate_producing_bacterium dup
species_dic = species_df[['file_name', 'bofTemplateType', 'M2', 'species_experiment']].set_index('file_name').T.to_dict(
    'list')

media_file = 'media/metacyc_growth_media_composition_df_trimed.tsv'
media_df = pd.read_csv(media_file, sep='\t', header=0, index_col=0)

essential_rxns = ['EX_h2o_e', 'EX_h_e', ]

media_name_dic = \
    {
        'M1: dGMM': 'M1',
        'M2: LAB': 'M2',
        'M3: dGMM+LAB': 'M3',
        'M4: dGMM+LAB low M/V': 'M4',
        'M5: dGMM+LAB exclude SCFA': 'M5',
        'M7: dGMM+LAB only monosacharides': 'M7',
        'M8: dGMM+LAB plus Mucin': 'M8',
        'M9: dGMM+LAB only Mucin': 'M9',
        'M10: dGMM+LAB 10% aminoacids': 'M10',
        'M11: dGMM+LAB excluding aromatic AA': 'M11',
        'M13: B.thetaiotaomicron MM': 'M13',
        'M14: C.perfiringens MM': 'M14',
        'M15A: E.coli MM1': 'M15A',
        'M15B: E.coli MM2': 'M15B',
        'M16: V.parvula defined medium': 'M16',
        'GMM: GMM': 'GMM',
        'mGAM: mGAM': 'mGAM',
        'BHI++: brain heart infusion': 'BHI++',
        'WCA: WCA': 'WCA'}

media_name_list = list(media_name_dic.keys())


def count_gaps(model):
    n = 0
    for rxn in model.reactions:
        if 'notes' in rxn.notes:
            if 'gap' in rxn.notes['notes']:
                n += 1
    return n


def test_media(model_version_1, media_name_list, media_df):
    result = {}
    for media_name in media_name_list:
        print(media_name)
        model_i = model_version_1.copy()
        # f0 = model_i.optimize().objective_value
        print(model_i.optimize())
        # model_i.reactions.EX_glc__D_e.lower_bound = 0

        media_df_temp = media_df.copy()
        media_df_temp['Constituents'] = media_df_temp['metacyc_id']

        media_i_df = media_df_temp[
            ['Constituents', media_name, 'exchange_rxns', 'transport_rxns']]  # select M2 columns

        media_i_df.columns = ['Constituents', 'Concentration', 'exchange_rxns', 'transport_rxns']
        media_i_dic = media.get_media_dic_form_df(media_i_df, media_name)
        keys = list(media_i_dic.keys())
        for k in keys:
            v = media_i_dic[k]
            if v[1] in essential_rxns:
                # pass
                media_i_dic.pop(k, None)
        try:

            media_model = media.update_media_from_dic(model_i, media_i_dic, reste_media=True,
                                                      essential_rxns=essential_rxns)

            ff = media_model.optimize().objective_value
        except:
            print('Error optimizing')
            ff = 0
        result[media_name] = ff
        print(ff, '\n\n')
    return result


# print(a)

for index in range(160, len(name_list)):
    name_i = name_list[index]
    print(index, name_i)

    model_draft = cobra.io.load_matlab_model(input_dir_2 + name_i + '_Metacyc.mat')

    draft_rxn_count = len(model_draft.reactions)
    draft_met_count = len(model_draft.metabolites)
    draft_gen_count = len(model_draft.genes)

    species_df.loc[species_df['file_name'] == name_i, 'draft_rxn_count'] = draft_rxn_count
    species_df.loc[species_df['file_name'] == name_i, 'draft_met_count'] = draft_met_count
    species_df.loc[species_df['file_name'] == name_i, 'draft_gen_count'] = draft_gen_count

    model_version_1 = cobra.io.load_json_model(input_dir_1 + name_i + '_add_gap.json')
    model_version_1.reactions.get_by_id('EX_o2_e').lower_bound = 0

    version_1_rxn_count = len(model_version_1.reactions)
    version_1_met_count = len(model_version_1.metabolites)
    version_1_gen_count = len(model_version_1.genes)
    version_1_gap_count = count_gaps(model_version_1)
    # version_1_growth_rate = model_version_1.optimize().objective_value

    species_df.loc[species_df['file_name'] == name_i, 'version_1_rxn_count'] = version_1_rxn_count
    species_df.loc[species_df['file_name'] == name_i, 'version_1_met_count'] = version_1_met_count
    species_df.loc[species_df['file_name'] == name_i, 'version_1_gen_count'] = version_1_gen_count
    species_df.loc[species_df['file_name'] == name_i, 'version_1_gap_count'] = version_1_gap_count

    result = test_media(model_version_1, media_name_list, media_df)
    for k, v in result.items():
        v_old = species_df.loc[species_df['file_name'] == name_i, k.split(':')[0] + '_model_growth_rate']
        species_df.loc[species_df['file_name'] == name_i, k.split(':')[0] + '_model_growth_rate'] = v
        try:
            if v == v_old:
                print(name_i, '\t', k, '\t v', v, '\t v_old', v_old)
        except:
            print('error')

    # species_df.loc[species_df['file_name'] == name_i,'version_1_growth_rate'] = version_1_growth_rate
    species_df.to_csv(species_table, sep='\t', index=False)

# temp_df = species_df.loc[(~species_df['M2'].isnull()),['organism_name','M2','version_1_growth_rate']]
# temp_df['M2 experiment'] = 1
# temp_df.loc[(temp_df['M2']=='0'),'M2 experiment'] = 0
# temp_df['Model'] = temp_df['version_1_growth_rate']
# temp_df.loc[(temp_df['version_1_growth_rate']>0.001),'Model'] = 1
#
#
# temp_df.index = temp_df['organism_name']
# temp_df = temp_df[['M2 experiment','Model']]
# # plt.imshow(a, cmap='hot', interpolation='nearest')
# # plt.show()
# import matplotlib.colors as mcolors
#
# cmap, norm = mcolors.from_levels_and_colors([0, 0.5, 1], ['grey','CornflowerBlue' ])
# # plt.pcolor(temp_df)
#
# plt.figure(figsize=(9, 10))
# ax = sns.heatmap(temp_df, cmap=cmap, linewidths=0.5, cbar_kws={"shrink": .30})
# ax.set_ylabel('Organism name', fontsize=14)
#
# plt.yticks(np.arange(0.5, len(temp_df.index), 1), temp_df.index,fontsize=10)
# plt.xticks(np.arange(0.5, len(temp_df.columns), 1), temp_df.columns)
# cbar = ax.collections[0].colorbar
# cbar.set_ticks([0.25,  .75, ])
# cbar.set_ticklabels(['not\ngrowth','growth'])
# plt.tight_layout()
# plt.show()
# plt.close()
