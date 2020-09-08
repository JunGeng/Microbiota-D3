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
                    print(gap_biomass_met, ': failed', ' * ' * 10)
                    model_template_i = model_template.copy()
                    model_template_i.add_reactions([model_i.reactions.get_by_id('object')])
                    model_template_i.objective = "object"
                    f = model_template_i.optimize()
                    if f.objective_value > 1e-10:
                        potential_gaps = set(f.fluxes[abs(f.fluxes) > 0.001].index) - {'object'}
                        gaps_set = potential_gaps - set([rxn.id for rxn in model_i.reactions])
                        all_gaps = all_gaps | gaps_set
                        result_dic[gap_biomass_met][1] = gaps_set
                        print('potential_gaps:', len(gaps_set))
                    else:
                        print('potential_gaps failed', ' * ' * 20)

        print('biomass gaps number:', len(all_gaps), all_gaps)
        return result_dic, all_gaps


def add_gaps(all_gaps, model_i, model_template):
    for gap in all_gaps:
        rxn = model_template.reactions.get_by_id(gap)
        if 'notes' in rxn.notes.keys():
            rxn.notes['notes'].append('gap')
        else:
            rxn.notes['notes'] = ['gap']
        model_i.add_reaction(rxn)


# %%<IO>
iML1515 = cobra.io.load_json_model('template_GEMs/iML1515_metacyc_updated_biomass.json')
iYO844 = cobra.io.load_json_model('template_GEMs/iYO844_metacyc_updated_biomass.json')
iAF692 = cobra.io.load_json_model('template_GEMs/iAF692_metacyc_updated_biomass.json')

species_table = 'species.tsv'
species_df = pd.read_csv(species_table, sep='\t')
name_list = list(species_df['file_name'])
output_dir = 'draft_GEMs/model_version_1.0/'

# experiment_db = pd.read_excel(r'media/NatureMicrobiology2018.3.514–522.xlsx', sheet_name='S4. Growth matrix')
# species_df = species_df.merge(experiment_db,how='left',left_on='species_experiment',right_on='designation in screen')
# species_df.to_csv('temp_species_4.tsv', sep='\t', index=False)

species_df = species_df.drop_duplicates(subset=['file_name'], keep='first')  # NOTE:butyrate_producing_bacterium dup
species_dic = species_df[['file_name', 'bofTemplateType', 'M2', 'species_experiment']].set_index('file_name').T.to_dict(
    'list')

# index = name_list.index('Escherichia_coli_O157:H7_str._Sakai') #16
index = name_list.index('[Eubacterium]_eligens_ATCC_27750')  # 17

for index in range(0, len(name_list)):  # range(0,len(name_list)):
    name_i = name_list[index]

    model_i = cobra.io.load_json_model('draft_GEMs/draft_add_biomass_and_media/' + name_i + '.json')
    # model_i = cobra.io.load_json_model(output_dir + name_i + '_add_gap.json')
    gram = species_dic[name_list[index]][0]

    gap_fill_bool = True
    if species_dic[name_list[index]][1] == '0':
        gap_fill_bool = False

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
    try:
        model_template.reactions.get_by_id('ATPASE-RXN').bounds = (0, 1000)
    except:
        print('no ATPASE-RXN')

    model_i_rxn_set = set([i.id for i in model_i.reactions])
    template_rxn_set = set([i.id for i in model_template.reactions])

    for id in model_i_rxn_set & template_rxn_set:
        rxn_temp_i = model_template.reactions.get_by_id(id)
        rxn_i = model_i.reactions.get_by_id(id)
        if 'p' in rxn_temp_i.compartments:
            continue
        if id.startswith("EX_"):
            rxn_temp_i.bounds = rxn_i.bounds
        else:
            gpr_i = rxn_i.gene_reaction_rule
            rxn_i.remove_from_model()
            rxn_temp_i.gene_reaction_rule = gpr_i
            model_i.add_reaction(rxn_temp_i)

    #  <gaps>   biomass_gaps_set = set([i.id for i in solution_part_biomass[0]])

    BIOMASS_ALL = cobra.Reaction('BIOMASS_ALL')
    model_i.add_reactions([BIOMASS_ALL])
    for k, v in model_template.reactions.get_by_id(biomass_id).metabolites.items():
        try:
            met = model_i.metabolites.get_by_id(k.id)
            BIOMASS_ALL.add_metabolites({met: v})
        except:
            print(k.id, ' not in model')

    obj_rea_id = 'BIOMASS_ALL'
    except_mets = {'WATER_c', 'ATP_c', 'NADP_c', }

    if gap_fill_bool:
        # %% <find gaps>
        result_dic, all_gaps = get_gaps(model_i, model_template, obj_rea_id, except_mets)

        # %% <add gaps>
        add_gaps(all_gaps, model_i, model_template)

    # %% <test >
    model_i.objective = obj_rea_id
    f = model_i.optimize()
    print(model_i.id, ' BIOMASS_ALL:', model_i.optimize())
    if f.objective_value < 1e-10:
        print('%%%' * 50)

    #     model_template.add_reaction(BIOMASS_ALL)
    #     model_template.objective = 'BIOMASS_ALL'
    #     f_temp = model_template.optimize()
    #     potential_gaps = set(f.fluxes[abs(f.fluxes) > 0.001].index) - {'object'}
    #
    #     for gap in potential_gaps:
    #         rea = model_template.reactions.get_by_id(gap)
    #         try:
    #             model_i.add_reaction(rea)
    #             if 'notes' in rea.notes.keys():
    #                 rea.notes['notes'].append('gap')
    #             else:
    #                 rea.notes['notes'] = ['gap']
    #         except:
    #             model_i.reactions.get_by_id(gap).remove
    #             model_template.reactions.get_by_id(gap)

    cobra.io.save_json_model(model_i, output_dir + name_i + '_add_gap.json')
    # cobra.io.save_matlab_model(model_i, output_dir + name_i + '_add_gap.mat')
    gemstool.io.gem2txt(model_i, output_dir + name_i + '_add_gap.txt')
    df_i = pd.DataFrame.from_dict(result_dic, orient='index', columns=['objective_value', 'gaps'])
    df_i['mets'] = df_i.index
    df_i['model'] = name_i
    df_i = df_i[['model', 'mets', 'objective_value', 'gaps']]
    df_i.to_csv(output_dir + name_i + '_gaps.tsv', sep='\t', index=False)
