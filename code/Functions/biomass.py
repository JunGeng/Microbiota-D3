#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/24/20

"""biomass.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import cobra


def update_biomass_from_dic(dic_i, observations_list=False, model=False):
    '''
    :description : function to get a model contain biomass reactions from a dic
    :param :
        dic_i(dic) : {'met_id': [coefficient, 'pool]}
                    eg. {'10-FORMYL-THF': [-0.000367, 'cofactor'],
                     'AMP': [-0.0046700000000000005, 'cofactor'],
                     'CPD-12125': [-0.000266, 'cofactor'],
                     'CDP': [-0.000251, 'cofactor'],}
        observations_list like:
            eg. ['dna', 'rna', 'protein', 'lipid', 'carbohydrate', 'cofactor', 'ion', 'other']

    :returns: cobra model
    '''
    if not model:
        biomass_model = cobra.Model('model')  # get a Model
    else:
        biomass_model = model.copy()

    if observations_list == False:
        observations_list = ['dna', 'rna', 'protein', 'lipid', 'carbohydrate', 'cofactor', 'ion', 'other']
    observations_list_mets = [i + '_c' for i in observations_list]
    observations_list_rxns = ['Biomass_' + i for i in observations_list]

    # %% <Biomass pool reactions>
    try:
        biomass_reactions = [model.reactions.get_by_id(i) for i in observations_list_rxns]
        for i in biomass_reactions: i.reaction = ' --> '  # reset biomass
        [rxn.add_metabolites({model.metabolites.get_by_id(met): 1}) for rxn, met in
         zip(biomass_reactions, observations_list_mets)]
    except:
        biomass_reactions = [cobra.Reaction(i) for i in observations_list_rxns]
        [rxn.add_metabolites({cobra.Metabolite(met): 1}) for rxn, met in
         zip(biomass_reactions, observations_list_mets)]

        biomass_model.add_reactions(biomass_reactions)

    for k, v in dic_i.items():
        # print(k)
        k = k + '_c'
        try:  # build metabolites
            met = biomass_model.metabolites.get_by_id(k)
        except:
            met = cobra.Metabolite(k)
        try:  # build reactions
            reaction_i = biomass_model.reactions.get_by_id('Biomass_' + v[1])
            reaction_i.add_metabolites({met: v[0]})
        except:
            print(k, v, 'Error')

    # %% <Biomass reaction>
    try:
        Biomass = model.reactions.get_by_id('Biomass')
        Biomass.reaction = ' --> '
        Biomass.add_metabolites({model.metabolites.get_by_id('Biomass_c'): 1})
    except:
        Biomass = cobra.Reaction('Biomass')
        biomass_model.add_reactions([Biomass])
        Biomass.add_metabolites({cobra.Metabolite('Biomass_c'): 1})
    [Biomass.add_metabolites({biomass_model.metabolites.get_by_id(met): -1}) for met in observations_list_mets]
    for met_i in biomass_model.metabolites:
        if met_i.id.endswith('_c'):
            if met_i.compartment != 'c':
                # print(met_i)
                met_i.compartment = 'c'

    return biomass_model
