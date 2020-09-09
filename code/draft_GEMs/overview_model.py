#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 9/3/20

"""overview_model.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import cobra
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns


os.chdir('../../data/')

# %% <IO>
species_table = 'species.tsv'
species_df = pd.read_csv(species_table, sep='\t')
name_list = list(species_df['file_name'])
input_dir_1 = 'draft_GEMs/model_version_1.0/'
input_dir_2 = 'draft_GEMs/draft_from_RAVEN_metacyc23_5/'

species_df = species_df.drop_duplicates(subset=['file_name'], keep='first')  # NOTE:butyrate_producing_bacterium dup
species_dic = species_df[['file_name', 'bofTemplateType', 'M2', 'species_experiment']].set_index('file_name').T.to_dict(
    'list')


def count_gaps(model):
    n = 0
    for rxn in model.reactions:
        if 'notes' in rxn.notes:
            if 'gap' in rxn.notes['notes']:
                n += 1
    return n


for index in range(0, len(name_list)):
    name_i = name_list[index]
    print(index, name_i)

    model_draft = cobra.io.load_matlab_model(input_dir_2 + name_i + '_Metacyc.mat')

    draft_rxn_count = len(model_draft.reactions)
    draft_met_count = len(model_draft.metabolites)
    draft_gen_count = len(model_draft.genes)

    species_df.loc[species_df['file_name'] == name_i, 'draft_rxn_count'] = draft_rxn_count
    species_df.loc[species_df['file_name'] == name_i, 'draft_met_count'] = draft_met_count
    species_df.loc[species_df['file_name'] == name_i, 'draft_gen_count'] = draft_gen_count

    model_version_1 = cobra.io.load_json_model(input_dir_1+name_i + '_add_gap.json')

    version_1_rxn_count = len(model_version_1.reactions)
    version_1_met_count = len(model_version_1.metabolites)
    version_1_gen_count = len(model_version_1.genes)
    version_1_gap_count = count_gaps(model_version_1)
    version_1_growth_rate = model_version_1.optimize().objective_value

    species_df.loc[species_df['file_name'] == name_i,'version_1_rxn_count'] = version_1_rxn_count
    species_df.loc[species_df['file_name'] == name_i,'version_1_met_count'] = version_1_met_count
    species_df.loc[species_df['file_name'] == name_i,'version_1_gen_count'] = version_1_gen_count
    species_df.loc[species_df['file_name'] == name_i,'version_1_gap_count'] = version_1_gap_count
    species_df.loc[species_df['file_name'] == name_i,'version_1_growth_rate'] = version_1_growth_rate

species_df.to_csv(species_table, sep='\t',index = False)
# the rest removed into overview_model.ipynd


# import seaborn as sns
# df = sns.load_dataset('tips')
#
# # Grouped violinplot
# sns.violinplot(x="day", y="total_bill", hue="smoker", data=df, palette="Pastel1")
# plt.show()


# %%
# fig, axs = plt.subplots(nrows=1, ncols=1, )
#
# # plot violin plot
# axs.violinplot(dataset=[species_df['draft_rxn_count'].values,
#                         species_df['draft_met_count'].values,
#                         species_df['draft_gen_count'].values, ],
#                showmeans=True,
#                showmedians=False)
# axs.set_title('Draft model overview', fontsize=16)
# # axs.set_xlabel('Four separate samples')
# axs.set_ylabel('Counts', fontsize=14)
# plt.setp(axs, xticks=[y + 1 for y in range(3)],
#          xticklabels=['Reactions', 'Metabolites', 'Genes'], )
# axs.tick_params(axis='x', which='major', labelsize=14)
#
# plt.show()
#
# fig, axs = plt.subplots(nrows=1, ncols=1, )
#
# # plot violin plot
# axs.violinplot(dataset=[species_df['version_1_gap_count'].values,
#                         species_df['version_1_met_count'].values,
#                         species_df['version_1_gen_count'].values,
#                         ],
#                showmeans=True,
#                showmedians=False)
# axs.set_title('Gap_filled model overview', fontsize=16)
# # axs.set_xlabel('Four separate samples')
# axs.set_ylabel('Counts', fontsize=14)
# plt.setp(axs, xticks=[y + 1 for y in range(3)],
#          xticklabels=['Reactions', 'Metabolites', 'Genes','Gap_filled'], )
# axs.tick_params(axis='x', which='major', labelsize=14)
#
# plt.show()

# fig, axs = plt.subplots(nrows=1, ncols=1, )
#
# # plot violin plot
# axs.violinplot(dataset=[species_df.loc[species_df['bofTemplateType'] == 'GramPositive','version_1_gap_count'].values,
#                         species_df.loc[species_df['bofTemplateType'] == 'GramNegative','version_1_gap_count'].values,
#                         species_df.loc[species_df['bofTemplateType'] == 'Archaea','version_1_gap_count'].values,
#                         ],
#                showmeans=True,
#                showmedians=False)
# axs.set_title('Gaps number in models', fontsize=16)
# # axs.set_xlabel('Four separate samples')
# axs.set_ylabel('Counts', fontsize=14)
# plt.setp(axs, xticks=[y + 1 for y in range(3)],
#          xticklabels=['GramPositive', 'GramNegative', 'Archaea','Gap_filled'], )
# axs.tick_params(axis='x', which='major', labelsize=14)
#
# plt.show()


fig, axs = plt.subplots(nrows=1, ncols=1, )

# plot violin plot
axs.violinplot(dataset=[species_df.loc[species_df['bofTemplateType'] == 'GramPositive','version_1_growth_rate'].values,
                        species_df.loc[species_df['bofTemplateType'] == 'GramNegative','version_1_growth_rate'].values,
                        species_df.loc[species_df['bofTemplateType'] == 'Archaea','version_1_growth_rate'].values,
                        ],
               showmeans=True,
               showmedians=False)
axs.set_title('Growth rate in M2 medium', fontsize=16)
# axs.set_xlabel('Four separate samples')
axs.set_ylabel('Counts', fontsize=14)
plt.setp(axs, xticks=[y + 1 for y in range(3)],
         xticklabels=['GramPositive', 'GramNegative', 'Archaea','Gap_filled'], )
axs.tick_params(axis='x', which='major', labelsize=14)

plt.show()

# a = species_df.loc[((species_df['bofTemplateType'] == 'GramPositive') & (species_df['M2'] != 0)),'version_1_growth_rate']
# b = species_df.loc[((species_df['bofTemplateType'] == 'GramNegative') & (species_df['M2'] != 0)),'version_1_growth_rate']
# c = species_df.loc[((species_df['bofTemplateType'] == 'Archaea') & (species_df['M2'] != 0)),'version_1_growth_rate']


temp_df = species_df.loc[(~species_df['M2'].isnull()),['organism_name','M2','version_1_growth_rate']]
temp_df['M2 experiment'] = 1
temp_df.loc[(temp_df['M2']=='0'),'M2 experiment'] = 0

temp_df['Model'] = temp_df['version_1_growth_rate']

temp_df.loc[(temp_df['version_1_growth_rate']>0.001),'Model'] = 1

# temp_df['Model'] = temp_df['M2 experiment']

temp_df.index = temp_df['organism_name']
temp_df = temp_df[['M2 experiment','Model']]
# plt.imshow(a, cmap='hot', interpolation='nearest')
# plt.show()
import matplotlib.colors as mcolors

cmap, norm = mcolors.from_levels_and_colors([0, 0.5, 1], ['grey','CornflowerBlue' ])
# plt.pcolor(temp_df)

plt.figure(figsize=(6.5, 6))
ax = sns.heatmap(temp_df, cmap=cmap, linewidths=0.5, cbar_kws={"shrink": .30})
ax.set_ylabel('Organism name', fontsize=14)

plt.yticks(np.arange(0.5, len(temp_df.index), 1), temp_df.index,fontsize=10)
plt.xticks(np.arange(0.5, len(temp_df.columns), 1), temp_df.columns)
cbar = ax.collections[0].colorbar
cbar.set_ticks([0.25,  .75, ])
cbar.set_ticklabels(['not\ngrowth','growth'])
plt.tight_layout()
plt.show()
plt.close()

