#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 8/21/20

"""media.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import re

import cobra


def get_exchange_rxns(model):
    exchange_rxns = []
    for rxn_i in model.boundary:
        rxn_i_metabolites = rxn_i.metabolites

        rxn_mets = [i.id for i in rxn_i_metabolites.keys()]

        exchange_rxns.append([re.sub('_[cep]$', '', rxn_mets[0]), rxn_i.id])
    return exchange_rxns


def get_transport_rxns(model):
    transport_rxns = []
    transport_rxns_candidates = (set(model.reactions) - set(model.boundary))

    transport_rxns_candidates = set([rxn for rxn in transport_rxns_candidates if len(rxn.compartments - {'p'}) >= 2])

    for rxn_i in transport_rxns_candidates:
        rxn_reactants = set([met.id for met in rxn_i.reactants])
        rxn_products = set([met.id for met in rxn_i.products])

        rxn_reactants_ = set([re.sub('_[cep]$', '', i) for i in rxn_reactants])
        rxn_products_ = set([re.sub('_[cep]$', '', i) for i in rxn_products])

        transported_mets = rxn_reactants_ & rxn_products_
        add_bool = False

        if len(rxn_reactants) == len(rxn_products_) == 1:
            add_bool = True
            # print(rxn_i)
        else:
            transported_mets = transported_mets - {'PROTON', 'WATER', 'Pi'}

            if len(transported_mets) != 0:
                add_bool = True
                if len(transported_mets) > 1:
                    print(rxn_i)
            elif rxn_reactants_ == rxn_products_:
                add_bool = True
                # print(rxn_i)
        if add_bool:
            for met in transported_mets:
                transport_rxns.append([met, rxn_i.id])
    return transport_rxns


def get_media_dic_form_df(media_i_df, media_name):
    '''
    :param :
            input:
                media_i_df (pandas dataframe):
                media_i_df.columns = ['Constituents', 'Concentration','exchange_rxns','transport_rxns']
            output:

                {'met_id': [coefficient,]}
                        eg. {'CL-': [51.607, 'EX_cl_e', '0mM'],
                                'CA+2': [0.27462000000000003, 'EX_ca2_e', 'CAt4; CA2abc1; CITt14'],
        media_name

    :returns: cobra model
    '''
    media_i_df = media_i_df.fillna('')
    media_i_df['Constituents'] = media_i_df['Constituents'].str.replace(' ', '')

    media_i_df['Concentration'] = media_i_df['Concentration'].str.replace(' ', '')  # get dict
    media_i_df.loc[media_i_df['Concentration'] == '', 'Concentration'] = '0mM'

    media_i_df.loc[:, 'coefficient'] = media_i_df['Concentration'].str.extract(r'([0-9\.]+)', expand=True)[0]
    media_i_df.loc[:, 'unit'] = media_i_df['Concentration'].str.extract(r'([µMnm]+)', expand=True)[0]

    media_i_df.loc[:, 'coefficient'] = media_i_df.apply(
        lambda x: float(x.coefficient) * 1e-6 if x['unit'] == 'nM' else float(x.coefficient), axis=1)
    media_i_df.loc[:, 'coefficient'] = media_i_df.apply(
        lambda x: float(x.coefficient) * 1e-3 if x.unit == 'µM' else float(x.coefficient), axis=1)
    if 'exchange_rxns' not in media_i_df.columns:
        media_i_df['exchange_rxns'] = ''
    if 'transport_rxns' not in media_i_df.columns:
        media_i_df['transport_rxns'] = ''

    media_i_dic = media_i_df[['Constituents', 'coefficient', 'exchange_rxns', 'transport_rxns']].set_index(
        'Constituents').T.to_dict('list')

    return media_i_dic


def check_meida_mets_in_model(model, met_id):
    '''get met_c, met_e from model
        if not in model, creat a new metabolite
    '''

    try:
        met_c = model.metabolites.get_by_id(met_id + '_c')
        met_e = model.metabolites.get_by_id(met_id + '_e')
    except:  # build metabolites
        met_c = cobra.Metabolite(met_id + '_c')
        met_e = cobra.Metabolite(met_id + '_e')
    met_c.compartment = 'c'
    met_e.compartment = 'e'
    return met_c, met_e


def check_exchange_transport_rxns_in_model(model, met_id, lower_bound, met_c, met_e,
                                           exchange_rxns_id='', transport_rxns_id=''):
    '''get teansport_reaction_i, exchange_reaction_i from model
        if not in model, creat a new reactions
        set lower_bound
    '''
    if exchange_rxns_id == '':
        exchange_rxns_id = 'EX_' + met_id + '_e'
    if transport_rxns_id == '':
        transport_rxns_id = 'TRANS_' + met_id
        teansport_reaction_i = cobra.Reaction(transport_rxns_id)
        teansport_reaction_i.add_metabolites({met_c: -1, met_e: 1})
        print(teansport_reaction_i)
        model.add_reactions([teansport_reaction_i])

    try:
        exchange_reaction_i = model.reactions.get_by_id(exchange_rxns_id)
    except:  # build reactions
        exchange_reaction_i = cobra.Reaction(exchange_rxns_id)
        exchange_reaction_i.add_metabolites({met_e: -1})

        # print(exchange_reaction_i)
        model.add_reactions([exchange_reaction_i])
    exchange_reaction_i.lower_bound = lower_bound


def changeout_media_from_dic(model_, media_i_dic):
    model = model_.copy()
    exchange_rxns_list = [i.id for i in model.reactions if i.id.startswith('Ex_')]
    # transport_rxns_list = [i.id for i in model.reactions if i.id.startswith('TRANS_')]
    # rest media to nothing
    for rxn_i in exchange_rxns_list:
        model.reactions.get_by_id(rxn_i).lower_bound = 0
        print('media rested')

    for met_id, v in media_i_dic.items():
        # print(k)
        # build metabolites
        met_c, met_e = check_meida_mets_in_model(model, met_id)
        lower_bound = -float(v[0])
        exchange_rxns_id = v[1]
        transport_rxns_id = v[2]
        check_exchange_transport_rxns_in_model(model, met_id, lower_bound, met_c, met_e, exchange_rxns_id,
                                               transport_rxns_id)

    return model
