#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/19/20

"""get_maped_templated_model.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import re

import cobra
import pandas as pd


os.chdir('../../data/template_GEMs/')


# %% <def functions>
def merge_metabolitesid(model1, new_id, old_id):
    model = model1.copy()
    try:
        model.metabolites.get_by_id(old_id).id = new_id

    except ValueError:

        for rea in model.metabolites.get_by_id(old_id).reactions:
            # rea.reaction = re.sub(r'(^|\b)'+id_in_tp+'(\b|$)', id_in_bigg, rea.reaction)
            rea_equ = rea.reaction
            model.reactions.get_by_id(rea.id).reaction = re.sub(r'\b%s\b' % old_id, new_id, rea.reaction)

        model.metabolites.get_by_id(old_id).remove_from_model()
        print(new_id + ' already in model!!!')

    except KeyError:
        print(old_id + ' not in model!!! skiped')

    return model


def merge_reactionsid(model1, new_id, old_id):
    model = model1.copy()

    try:
        rea1 = model.reactions.get_by_id(new_id)
    except:
        try:
            rea2 = model.reactions.get_by_id(old_id)
            model.reactions.get_by_id(old_id).id = new_id

            # print(new_id + ' not in model, just replase rea_id')
        except:
            print(old_id + ' not in model')

    try:
        model.reactions.get_by_id(old_id).id = new_id

    except ValueError:

        print('\033[1;34;47m')
        print(model.id + 'reaction id: ' + new_id + " and " + old_id + ' both in model!!! check it by hand!!!')
        print('\033[0;30;48m')

        # rea1 = model.reactions.get_by_id(new_id)
        # rea2 = model.reactions.get_by_id(old_id)
        #
        # # rea1.gene_reaction_rule = merge_gprule(rea1.gene_reaction_rule, rea2.gene_reaction_rule)
        #
        # model.reactions.get_by_id(old_id).remove_from_model()

    except KeyError:
        pass
        # print(old_id + ' not in model!!! skiped')

    return model


# %% <IO>

iML1515 = cobra.io.load_json_model('iML1515.json')
# iJO1366 = cobra.io.load_json_model('../archive/iML1515.json')
iYO844 = cobra.io.load_json_model('iYO844.json')
iAF692 = cobra.io.load_json_model('iAF692.json')

mets_table = pd.read_csv('template_mets_id_map.tsv', sep='\t', index_col=0)
rxns_table = pd.read_csv('template_rxns_id_map.tsv', sep='\t', index_col=0)

# %% <get dic from table>

mets_table = mets_table.dropna(subset=['metacyc_id'])  # remove nan

# add _c _e _p
mets_table.loc[mets_table['id_c'].str.endswith('_c'), 'metacyc_id'] = mets_table.loc[mets_table['id_c'].str.endswith(
    '_c'), 'metacyc_id'].apply(lambda x: x.replace(';', '_c;') + '_c')
mets_table.loc[mets_table['id_c'].str.endswith('_e'), 'metacyc_id'] = mets_table.loc[mets_table['id_c'].str.endswith(
    '_e'), 'metacyc_id'].apply(lambda x: x.replace(';', '_e;') + '_e')
mets_table.loc[mets_table['id_c'].str.endswith('_p'), 'metacyc_id'] = mets_table.loc[mets_table['id_c'].str.endswith(
    '_p'), 'metacyc_id'].apply(lambda x: x.replace(';', '_p;') + '_p')

mets_dic = mets_table[['id_c', 'metacyc_id']].set_index('id_c').T.to_dict('list')

rxn_table = rxns_table.dropna(subset=['metacyc_id'])  # remove nan
rxns_dic = rxn_table[['id', 'metacyc_id']].set_index('id').T.to_dict('list')

# %% <iML1515>
iML1515_metacyc = iML1515.copy()
iML1515_mets_set = set([i.id for i in iML1515.metabolites])
for old_id, new_id in mets_dic.items():
    if old_id not in iML1515_mets_set:
        continue
    new_id = new_id[0]
    if ';' in new_id:  # TODO more id process,skip or chose one
        continue
        new_id = new_id.split('; ')[0]
    iML1515_metacyc = merge_metabolitesid(iML1515_metacyc, new_id, old_id)

iML1515_rxns_set = set([i.id for i in iYO844.reactions])

for old_id, new_id in rxns_dic.items():
    if old_id not in iML1515_rxns_set:
        continue
    new_id = new_id[0]
    if ';' in new_id:  # TODO more id process
        continue
        new_id = new_id.split('; ')[0]
    iML1515_metacyc = merge_reactionsid(iML1515_metacyc, new_id, old_id)
print('iML1515 initial opt',iML1515.optimize().objective_value)
print('iML1515_metacyc  opt',iML1515_metacyc.optimize().objective_value)

# %% <iYO844>
iYO844_metacyc = iYO844.copy()
iYO844_mets_set = set([i.id for i in iYO844.metabolites])
for old_id, new_id in mets_dic.items():
    if old_id not in iYO844_mets_set:
        continue
    new_id = new_id[0]
    if ';' in new_id:  # TODO more id process,skip or chose one
        continue
        new_id = new_id.split('; ')[0]
    iYO844_metacyc = merge_metabolitesid(iYO844_metacyc, new_id, old_id)


iYO844_rxns_set = set([i.id for i in iYO844.reactions])

for old_id, new_id in rxns_dic.items():
    if old_id not in iYO844_rxns_set:
        continue
    new_id = new_id[0]
    if ';' in new_id:  # TODO more id process
        continue
        new_id = new_id.split('; ')[0]
    iYO844_metacyc = merge_reactionsid(iYO844_metacyc, new_id, old_id)
print('iYO844 initial opt',iYO844.optimize().objective_value)
print('iYO844_metacyc  opt',iYO844_metacyc.optimize().objective_value)

# %% <iAF692>
iAF692_metacyc = iAF692.copy()
iAF692_mets_set = set([i.id for i in iAF692.metabolites])
for old_id, new_id in mets_dic.items():
    if old_id not in iAF692_mets_set:
        continue
    new_id = new_id[0]
    if ';' in new_id:  # TODO more id process,skip or chose one
        continue
        new_id = new_id.split('; ')[0]
    iAF692_metacyc = merge_metabolitesid(iAF692_metacyc, new_id, old_id)
print(iAF692_metacyc.optimize().objective_value)

iAF692_rxns_set = set([i.id for i in iAF692.reactions])

for old_id, new_id in rxns_dic.items():
    if old_id not in iAF692_rxns_set:
        continue
    new_id = new_id[0]
    if ';' in new_id:  # TODO more id process
        continue
        new_id = new_id.split('; ')[0]
    iAF692_metacyc = merge_reactionsid(iAF692_metacyc, new_id, old_id)

print('iAF692 initial opt',iAF692.optimize().objective_value)
print('iAF692_metacyc  opt',iAF692_metacyc.optimize().objective_value)

# %% <write model>
cobra.io.save_json_model(iML1515_metacyc,'iML1515_metacyc.json')
cobra.io.save_json_model(iYO844_metacyc,'iYO844_metacyc.json')
cobra.io.save_json_model(iAF692_metacyc,'iAF692_metacyc.json')