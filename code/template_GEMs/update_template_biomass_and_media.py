#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/24/20

"""checkout_biomass_and_media.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os

import cobra

os.chdir('../../data/template_GEMs/')

# %% <IO>
iML1515_initial = cobra.io.load_json_model('iML1515.json')
iYO844_initial = cobra.io.load_json_model('iYO844.json')
iAF692_initial = cobra.io.load_json_model('iAF692.json')

iML1515 = cobra.io.load_json_model('iML1515_metacyc.json')
iYO844 = cobra.io.load_json_model('iYO844_metacyc.json')
iAF692 = cobra.io.load_json_model('iAF692_metacyc.json')

biomass_model_list = ['biomass_negative_reactions.json',
                      'biomass_positive_reactions_iYO844.json',
                      'biomass_archaea_reactions.json'
                      ]

biomass_n_model = cobra.io.load_json_model('../biomass/' + biomass_model_list[0])
biomass_p_model = cobra.io.load_json_model('../biomass/' + biomass_model_list[1])
biomass_a_model = cobra.io.load_json_model('../biomass/' + biomass_model_list[2])

media_model_list = ['media_model_m2.json']
media_model_path = '../media/' + media_model_list[0]
media_model = cobra.io.load_json_model(media_model_path)

# %% <update biomass> # get species name list
iML1515_ = iML1515.copy()
iYO844_ = iYO844.copy()
iAF692_ = iAF692.copy()
for model_i in [iML1515_, iYO844_, iAF692_]:
    # %% <load biomass reactions type: model>

    if model_i.id == 'iML1515':
        biomass_model = biomass_n_model

    elif model_i.id == 'iYO844':
        biomass_model = biomass_p_model

    elif model_i.id == 'iAF692':
        biomass_model = biomass_a_model

    else:
        print('Error')

    # %% <merge draft model and biomass model>
    model_i = model_i.merge(biomass_model)

    # %% <merge draft model and biomass model>
    # model_i = model_i.merge(media_model)

    # %% <check model>
    for met_i in model_i.metabolites:
        if met_i.id.endswith('_c'):
            if met_i.compartment != 'c':
                print(met_i)
                met_i.compartment = 'c'

        elif met_i.id.endswith('_e'):
            if met_i.compartment != 'e':
                print(met_i)
                met_i.compartment = 'e'
        else:
            pass
            # print(met_i)
    model_i.objects = 'Biomass'
    f = model_i.optimize().objective_value
    print(model_i.id, ' opt:', ': ', f)

    # %% <write model>
    # cobra.io.save_json_model(model_i, output_file_dir + model_i.id + '.json')
    # cobra.io.save_matlab_model(model_i, output_file_dir + model_i.id + '.mat')

# %% <check growth rate>
print('iML1515_initial \topt:', iML1515_initial.optimize().objective_value)
print('iML1515_metacyc \topt:', iML1515.optimize().objective_value)
print('iML1515_new \topt:', iML1515_.optimize().objective_value)

print('iYO844_initial \topt:', iYO844_initial.optimize().objective_value)
print('iYO844_metacyc \topt:', iYO844.optimize().objective_value)
print('iYO844_new \topt:', iYO844_.optimize().objective_value)

print('iAF692_initial \topt:', iAF692_initial.optimize().objective_value)
print('iAF692_metacyc \topt:', iAF692.optimize().objective_value)
print('iAF692_new \topt:', iAF692_.optimize().objective_value)

# %% <write model>
cobra.io.save_json_model(iML1515_, 'iML1515_metacyc.json')
cobra.io.save_json_model(iYO844_, 'iYO844_metacyc.json')
cobra.io.save_json_model(iAF692_, 'iAF692_metacyc.json')
# import gemstool
# gemstool.io.solution2txt(iAF692_.optimize(),iAF692_,'temp.txt')
