#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/8/20

"""branch_get_template_con.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os

import cobra
import numpy as np

os.chdir('../../data/')


def get_dna_or_rna_weight_fraction(model, biomass_rxn_id, type='dna'):
    '''
    dna_and_rna_weight_fraction
    dna_weight_fraction = (sum('datp_c', 'dctp_c', 'dgtp_c', 'dttp_c') - ppi ) / 1000 mg
    rna_weight_fraction = (sum('atp_c', 'ctp_c', 'gtp_c', 'utp_c') - ppi ) / 1000 mg
                                                # atp_c coefficient should - pi coefficient

    '''
    biomass_rxn = model.reactions.get_by_id(biomass_rxn_id)

    if type == 'dna':
        dna_or_rna_mets_id = ['datp_c', 'dctp_c', 'dgtp_c', 'dttp_c']
    elif type == 'rna':
        dna_or_rna_mets_id = ['atp_c', 'ctp_c', 'gtp_c', 'utp_c']

    dna_or_rna_mets = [model.metabolites.get_by_id(i) for i in dna_or_rna_mets_id]

    dna_or_rna_mets_coff = [biomass_rxn.metabolites[i] for i in dna_or_rna_mets]  # coefficients
    if type == 'rna':
        dna_or_rna_mets_coff[0] = dna_or_rna_mets_coff[0] + biomass_rxn.metabolites[
            model.metabolites.get_by_id('pi_c')]  # ATP coefficients

    dna_or_rna_mets_molar_mass = [i.formula_weight for i in dna_or_rna_mets]  # molar_mass

    dna_or_rna_mets_weight = np.array(dna_or_rna_mets_coff) * np.array(
        dna_or_rna_mets_molar_mass)  # weight = coefficients * molar_mass (mg)
    ppi_weight = sum(dna_or_rna_mets_coff) * model.metabolites.get_by_id('ppi_c').formula_weight
    dna_or_rna_mets_wight = abs(dna_or_rna_mets_weight.sum()) - abs(ppi_weight)
    dna_or_rna_weight_fraction = dna_or_rna_mets_wight / 1000  # units: 1g biomass = 1000 mg

    return dna_or_rna_weight_fraction


# %% <io>

# input_file = 'biomass/biomass_constituents_id_map.tsv'

iML1515_path = 'template_GEMs/iML1515.json'
iYO844_path = 'template_GEMs/iYO844.json'
iAF692_path = 'template_GEMs/iAF692.json'

# %% <load model>

iML1515 = cobra.io.load_json_model(iML1515_path)
iYO844 = cobra.io.load_json_model(iYO844_path)
iAF692 = cobra.io.load_json_model(iAF692_path)

# %% <calculate weight_fraction>
##  iML1515
model = iML1515
biomass_rxn_id = 'BIOMASS_Ec_iML1515_core_75p37M'
dna_weight_fraction = get_dna_or_rna_weight_fraction(model, biomass_rxn_id, type='dna')
rna_weight_fraction = get_dna_or_rna_weight_fraction(model, biomass_rxn_id, type='rna')
print('iML1515 weight_fractions:')
print('\t\t\t\tdna_weight_fraction', dna_weight_fraction)
print('\t\t\t\trna_weight_fraction', rna_weight_fraction)
print('\t\t\t\treference value : dna_weight_fraction == 0.031 , rna_weight_fraction  == 0.205')

##  iYO844
model = iYO844
biomass_rxn_id = 'BIOMASS_BS_10'
dna_weight_fraction = get_dna_or_rna_weight_fraction(model, biomass_rxn_id, type='dna')
rna_weight_fraction = get_dna_or_rna_weight_fraction(model, biomass_rxn_id, type='rna')
print('iYO844 weight_fractions:')
print('\t\t\t\tdna_weight_fraction', dna_weight_fraction)
print('\t\t\t\trna_weight_fraction', rna_weight_fraction)

##  iAF692
model = iAF692
biomass_rxn_id = 'BIOMASS_Mb_30'
dna_weight_fraction = get_dna_or_rna_weight_fraction(model, biomass_rxn_id, type='dna')
rna_weight_fraction = get_dna_or_rna_weight_fraction(model, biomass_rxn_id, type='rna')
print('iAF692 weight_fractions:')
print('\t\t\t\tdna_weight_fraction', dna_weight_fraction)
print('\t\t\t\trna_weight_fraction', rna_weight_fraction)


