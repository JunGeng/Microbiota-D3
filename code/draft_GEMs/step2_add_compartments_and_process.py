#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/31/20

"""step2_add_compartments.py
:description : script to add compartments to metabolites id and process transport reaction empty mets problem
:param : 
:returns: 
:rtype: 
"""
import os

import cobra
import gemstool
import pandas as pd


def _add_compartments(model):
    model.compartments = {'c': 'cytosol',
                          'e': 'extracellular space'}

    for met_i in model.metabolites:
        if not met_i.id.endswith('_c'):
            met_i.id = met_i.id + '_c'
        else:
            print(met_i)
        met_i.compartment = 'c'


def _check_compartments(model):
    if model.compartments == {}:
        model.compartments = {'c': 'cytosol',
                              'e': 'extracellular space'}

    for met_i in model.metabolites:
        if met_i.id.endswith('_c'):
            met_i.compartment = 'c'
        elif met_i.id.endswith('_e'):
            met_i.compartment = 'e'


def _process_mets_id(id):
    if id.endswith('_c[c]'):
        id = id.replace('_c[c]', '_c')
    elif id.endswith('_p[p]'):
        id = id.replace('_p[p]', '_e')
    elif id.endswith('_e[e]'):
        id = id.replace('_e[e]', '_e')
    elif id.endswith('_m[m]'):
        id = id.replace('_m[m]', '_e')

    elif id.endswith('_CCO-IN[CCO-IN]'):
        id = id.replace('_CCO-IN[CCO-IN]', '_c')
    elif id.endswith('_CCO-OUT[CCO-OUT]'):
        id = id.replace('_CCO-OUT[CCO-OUT]', '_e')

    elif id.endswith('_CCO-SIDE-1[CCO-SIDE-1]'):
        id = id.replace('_CCO-SIDE-1[CCO-SIDE-1]', '_c')
    elif id.endswith('_CCO-SIDE-2[CCO-SIDE-2]'):
        id = id.replace('_CCO-SIDE-2[CCO-SIDE-2]', '_e')

    return id


def _replace_equation(rxn_initial, rxn_new, model):
    rxn_initial.reaction = ' --> '
    rxn_initial.bounds = rxn_new.bounds
    for k, v in rxn_new.metabolites.items():
        id = _process_mets_id(k.id)
        try:
            met = model.metabolites.get_by_id(id)
        except:
            met = cobra.Metabolite(id)
        rxn_initial.add_metabolites({met: v})


# %% <IO>
os.chdir('../../data/draft_GEMs/draft_from_RAVEN_metacyc23_5/')

species_table = '../../species.tsv'
species_df = pd.read_csv(species_table, sep='\t')
# species_df['file_name'] = species_df['organism_name'].str.replace(' ', '_')
# species_df['file_name'] = species_df['file_name'].str.replace('/', '_')
# species_df.to_csv(species_table, sep='\t', index=False)

name_list = list(species_df['file_name'])

model_db = cobra.io.load_json_model('../../initial_data/data_from_database/MetaCyc_23.5/MetaCyc_cobra.json')
rxns_db_set = set([i.id for i in model_db.reactions])

removed_list = set([])
replaced_list = set([])

# %% <process>

for index in range(0, len(name_list)):  # len(name_list)
    name_i = name_list[index]
    print(index, name_i)

    model_i = cobra.io.load_matlab_model(name_i + '_Metacyc.mat')
    _add_compartments(model_i)  # _add_compartments

    rxn_list = [rxn for rxn in model_i.reactions]
    for rxn_i in rxn_list:
        id = rxn_i.id
        if id[0].isdigit():
            id = '_' + id

        if id in rxns_db_set:
            rxn_py = model_db.reactions.get_by_id(id)
            if len(rxn_py.compartments) > 1:
                replaced_list.add(rxn_i.id)
                _replace_equation(rxn_i, rxn_py, model_i)
                # print(rxn_py)
        else:
            mets_set = set([k.id for k, v in rxn_i.metabolites.items()])
            mets_set = mets_set - {'PROTON_c', 'WATER_c', 'ATP_c', 'ADP_c', 'Pi_c'}
            if len(mets_set) < 1:
                removed_list.add(rxn_i.id)
                rxn_i.remove_from_model()
                # print(rxn_i)
    _check_compartments(model_i)
    cobra.io.save_json_model(model_i, name_i + '_add_compartments.json')

    gemstool.io.gem2txt(model_i, name_i + '_add_compartments.txt')

with open('removed_list.txt', 'w') as f:
    f.write('; '.join(removed_list))
with open('replaced_list.txt', 'w') as f:
    f.write('; '.join(replaced_list))
