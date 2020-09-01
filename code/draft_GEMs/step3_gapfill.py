#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/26/20

"""step_3_gapfill.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import cobra
import gemstool
import pandas as pd

os.chdir('../../data/')

# %%<IO>
iML1515 = cobra.io.load_json_model('template_GEMs/iML1515_metacyc.json')
iYO844 = cobra.io.load_json_model('template_GEMs/iYO844_metacyc.json')
iAF692 = cobra.io.load_json_model('template_GEMs/iAF692_metacyc.json')

species_table = 'species.tsv'
species_df = pd.read_csv(species_table, sep='\t')
name_list = list(species_df['file_name'])
output_dir = 'draft_GEMs/model_version_1.0/'


species_df = species_df.drop_duplicates(subset=['file_name'], keep='first')  # NOTE:butyrate_producing_bacterium dup
species_dic = species_df[['file_name', 'bofTemplateType']].set_index('file_name').T.to_dict('list')

index = name_list.index('Escherichia_coli_O157:H7_str._Sakai')
name_i = name_list[index]

model_i = cobra.io.load_json_model('draft_GEMs/draft_add_biomass_and_media/' + name_i + '.json')
gram = species_dic[name_list[index]][0]

if gram == 'GramNegative':
    model_template = iML1515.copy()
    biomass_id = 'BIOMASS_Ec_iML1515_core_75p37M'
elif gram == 'GramPositive':
    model_template = iYO844.copy()
    biomass_id = 'BIOMASS_BS_10'

elif gram == 'Archaea':
    model_template = iAF692.copy()
    biomass_id = 'BIOMASS_Mb_30'
else:
    print('Error')


#  <gaps>   biomass_gaps_set = set([i.id for i in solution_part_biomass[0]])
def get_gaps(model_gap, model_template, obj_rea_id, except_mets={'h2o_c'}):
    model = model_gap.copy()
    result_dic = {obj_rea_id: [''] * 2}
    all_gaps = set([])

    model.objective = obj_rea_id
    obj_initial = model_gap.optimize().objective_value

    result_dic[obj_rea_id][0] = obj_initial
    gap_bool = True

    if obj_initial > 1e-10:  # no gap
        gap_bool = False
        print('no gap')
        return result_dic, all_gaps

    if gap_bool:  # gapfill all
        try:
            solution_whole_biomass = cobra.flux_analysis.gapfill(model, model_template, demand_reactions=False)
            gaps_set = set([i.id for i in solution_whole_biomass[0]])
            all_gaps = all_gaps | gaps_set
            result_dic[obj_rea_id][1] = gaps_set
            gap_bool = False
            print('gapfill whole succeed')
            return result_dic, all_gaps
        except:
            print('gapfill whole failed')

    if gap_bool:  # gapfill partly(met by met)
        except_mets = except_mets | {'atp_c', 'coa_c', 'h2o_c', 'nad_c', }

        gap_biomass_mets = [i.id for i in model.reactions.get_by_id(obj_rea_id).reactants]

        gap_biomass_mets = set(gap_biomass_mets) - except_mets

        for gap_biomass_met in gap_biomass_mets:
            result_dic[gap_biomass_met] = [''] * 2

            # model_template = model.copy()
            model_i = model.copy()

            rea_temp = cobra.Reaction('object')
            model_i.add_reactions([rea_temp])

            model_i.reactions.get_by_id('object').reaction = ' --> '
            model_i.reactions.get_by_id('object').add_metabolites(
                {model_i.metabolites.get_by_id(gap_biomass_met): -1, })

            # print('biomass apart optimize: ',gap_biomass_met)
            model_i.objective = "object"
            obj_initial = model_i.optimize().objective_value
            result_dic[gap_biomass_met][0] = obj_initial

            if obj_initial > 1e-10:  # no gap
                print(gap_biomass_met, ': no gap')
                continue
            else:
                try:
                    solution_part_biomass = cobra.flux_analysis.gapfill(model_i, model_template,
                                                                        demand_reactions=False)
                    gaps_set = set([i.id for i in solution_part_biomass[0]])
                    all_gaps = all_gaps | gaps_set
                    result_dic[gap_biomass_met][1] = gaps_set
                    print(gap_biomass_met, ': succeed')
                except:
                    print(gap_biomass_met, ': failed', ' * ' * 20)

        print('biomass gaps number:', len(all_gaps), all_gaps)
        return result_dic, all_gaps


BIOMASS = cobra.Reaction('BIOMASS_ALL')
model_i.add_reactions([BIOMASS])
for k, v in model_template.reactions.get_by_id(biomass_id).metabolites.items():
    try:
        met = model_i.metabolites.get_by_id(k.id)
        BIOMASS.add_metabolites({met: v})
    except:
        print(k.id, ' not in model')

obj_rea_id = 'BIOMASS_ALL'

# %% <find gaps>
result_dic, all_gaps = get_gaps(model_i, model_template, obj_rea_id, )

# %% <add gaps>
for gap in all_gaps:
    rea = model_template.reactions.get_by_id(gap)
    if 'notes' in rea.notes.keys():
        rea.notes['notes'].append('gap')
    else:
        rea.notes['notes'] = ['gap']
    model_i.add_reaction(rea)
# %% <test >
model_i.objective = obj_rea_id
print(model_i.id, ' BIOMASS_ALL:', model_i.optimize())

cobra.io.save_json_model(model_i, output_dir + name_i + '_add_gap.json')
cobra.io.save_matlab_model(model_i, output_dir + name_i + '_add_gap.mat')
gemstool.io.gem2txt(model_i, output_dir + name_i + '_add_gap.txt')
