#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/25/20

"""step2_draftGEMs.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os
import re

import cobra
import pandas as pd

os.chdir('../../data/draft_GEMs/')


def select_blast(result_file, best_match=True, evalue=10 ** -10, pident=40, length=200, bitscore=0, ppos=0,
                 qcovs=0):
    '''
    #    find BBH(Bidirectional Best Hits) from blast results
    #     balstcmd = 'blastp -db ' + db + ' -query ' + seq + ' -out ' + outfile +  ' -evalue 10e-5  -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"'

    selsect the result hits from blast result
    :param result1: blast result file -outfmt "6 qseqid sseqid evalue pident length bitscore ppos"
    :param best_match: BBH(Bidirectional Best Hits) from blast results, best mached and  unique
    :param evalue:  evalue = 10**-10
    :param pident:  ident
    :param length:  matched length
    :param bitscore:
    :param ppos:
    :return:    a dataframe meet the params
    examlple:
            result1 = 'Lreuteri_refseq_v01_in_Lreuteri_refseq_v02.csv'
            result_df = select_blast(result1, best_match=True,evalue = 10**-10, pident = 40, length = 0, bitscore = 0, ppos = 0)
    '''
    names = ['qseqid', 'sseqid', 'evalue', 'pident', 'length', 'bitscore', 'ppos', 'qcovs']

    result_df = pd.read_csv(result_file, sep='\t', names=names)

    dfi = result_df.copy()
    dfi = dfi[(dfi["evalue"] <= evalue) & (dfi["pident"] >= pident) & (dfi["length"] >= length) & (
            dfi["bitscore"] >= bitscore) & (dfi["ppos"] >= ppos) & (dfi["qcovs"] >= qcovs)]
    dfi = dfi.copy()
    dfi['sseqid'] = dfi.sseqid.apply(lambda x: re.sub('(\|$)|(^.*?\|)', '', x))
    dfi = dfi.sort_values(['qseqid', "evalue", 'pident', 'bitscore'], ascending=[True, True, False, False])

    if best_match:
        dfi = dfi.drop_duplicates(subset='qseqid', keep='first')
    assert isinstance(dfi, object)
    return_df = dfi

    return return_df


def gpr2log(gpr_i, othergeneset, inset=True, notinset=False, emptyset=''):
    torf = ''
    gpr_i = ' ' + gpr_i + ' '

    if not gpr_i:
        if emptyset == '':
            new_gpr_i = gpr_i
        else:
            new_gpr_i = emptyset
    else:
        _, geneset = cobra.core.gene.parse_gpr(gpr_i)
        new_gpr_i = gpr_i
        for gen_i in geneset:
            if gen_i in othergeneset:
                temp_gen_i = inset
            else:
                temp_gen_i = notinset
            new_gpr_i = new_gpr_i.replace(' ' + gen_i + ' ', ' ' + str(temp_gen_i) + ' ')

        if (inset == True and notinset == False):
            torf = eval(new_gpr_i)
    return new_gpr_i, torf


def get_draft_from_template(tp_model1, blast_result_df, remove_missing_genes=True):
    '''
    build a model based on template model
    :param tp_model: template model
    :param blase_result_df: balst result
    :return:  draft model
    '''
    tp_model = tp_model1.copy()
    tp_gene_list = blast_result_df['sseqid'].tolist()
    # my_gene_list = blast_result_df['qseqid'].tolist()
    tp_model_gene_list = [i.id for i in tp_model.genes]

    union_gene_set = set(tp_gene_list) & set(tp_model_gene_list)

    if len(union_gene_set) == 0:
        print('no match genes')
        return None

    temp_df = pd.DataFrame()
    temp_df['sseqid'] = list(union_gene_set)
    temp_df = temp_df.merge(blast_result_df, left_on='sseqid', right_on='sseqid')
    tp_gene_list = temp_df['sseqid'].tolist()
    my_gene_list = temp_df['qseqid'].tolist()

    model = cobra.Model()
    model.description = 'GEM from template' + tp_model.id

    all_reactions = set()
    for i in union_gene_set:
        for ii in tp_model.genes.get_by_id(i).reactions:
            all_reactions.add(ii)

    for rea in all_reactions:
        new_gpr_i, torf = gpr2log(rea.gene_reaction_rule, tp_gene_list)

        if not remove_missing_genes:
            if 'True' in new_gpr_i:
                torf = True

        if torf:
            rea.notes['from'] = [tp_model.id]
            # New gene_reaction_rule
            new_gene_reaction_rule = rea.gene_reaction_rule
            if not new_gene_reaction_rule.startswith(' '):
                new_gene_reaction_rule = ' ' + new_gene_reaction_rule
            if not new_gene_reaction_rule.endswith(' '):
                new_gene_reaction_rule = new_gene_reaction_rule + ' '

            for ii in rea.genes:
                if ii.id in tp_gene_list:
                    k = tp_gene_list.index(ii.id)
                    new_gene_reaction_rule = new_gene_reaction_rule.replace(' ' + ii.id + ' ',
                                                                            ' ' + my_gene_list[k] + ' ')

            rea.gene_reaction_rule = new_gene_reaction_rule
            model.add_reactions([rea])

    removegeneslist = []
    for i in model.genes:
        if i.id not in my_gene_list:
            # i.id=i.id+'_missing'
            removegeneslist.append(i.id)

    if remove_missing_genes:
        cobra.manipulation.remove_genes(model, removegeneslist)

    else:
        for i in removegeneslist:
            reas = model.genes.get_by_id(i).reactions
            for rea in reas:
                rea.gene_reaction_rule = re.sub(i + '(?!_missing)', i + '_missing', rea.gene_reaction_rule)
        list2 = set()
        for gene in model.genes:
            if len(gene.reactions) == 0:
                list2.add(gene.id)

        cobra.manipulation.remove_genes(model, list2)

    return model


# %% <process gene_map_df>
# gene_map_df = pd.DataFrame()
# lind_dic = {}
# for line in open('../initial_data/data_from_database/MetaCyc_23.5/genes.dat'):
#     if line.startswith('UNIQUE-ID'):
#         lind_dic['UNIQUE-ID'] = line.replace('\n','').replace('UNIQUE-ID - ','').replace('-','_')
#     elif line.startswith('COMMON-NAME'):
#         lind_dic['COMMON-NAME'] = line.replace('\n','').replace('COMMON-NAME - ','')
#     elif line.startswith('PRODUCT'):
#         lind_dic['PRODUCT'] = line.replace('\n','').replace('PRODUCT - ','')
#     elif line.startswith('//'):
#         gene_map_df = gene_map_df.append(lind_dic, ignore_index=True)
#         lind_dic = {}
# gene_map_df.to_csv('gene_map.txt', sep='\t', index=False)
gene_map_df = pd.read_csv('gene_map.txt', sep='\t')

# %% <select blust result by cutoff>

file_list = os.listdir('blast/')
blast_result_q_in_s_list = [i for i in file_list if i.startswith('blast_result_q_in_MetaCyc_')]
blast_result_s_in_q_list = [i.replace('blast_result_q_in_MetaCyc_', 'blast_result_q_in_MetaCyc_') for i in
                            blast_result_q_in_s_list]

species_name_list_temp = [i.replace('blast_result_q_in_MetaCyc_', '') for i in blast_result_q_in_s_list]
species_name_list = [i.replace('.csv', '') for i in species_name_list_temp]

MetaCyc_model = cobra.io.read_sbml_model('../initial_data/data_from_database/MetaCyc_23.5/MetaCyc_cobra.xml')
for i in MetaCyc_model.reactions:
    gpr_i = i.gene_reaction_rule
    if not gpr_i.startswith(' '):
        i.gene_reaction_rule = ' ' + i.gene_reaction_rule
    if not gpr_i.endswith(' '):
        i.gene_reaction_rule = i.gene_reaction_rule + ''

# blast_result_q_in_s = 'blast/' + blast_result_q_in_s_list[0]
# blast_result_s_in_q = 'blast/' + blast_result_s_in_q_list[0]
# species_name = species_name_list[0]

for index in range(0, len(species_name_list)):
    blast_result_q_in_s = 'blast/' + blast_result_q_in_s_list[index]
    species_name = species_name_list[index]
    print(species_name)
    blast_result_df = select_blast(blast_result_q_in_s, best_match=True,
                                   evalue=10 ** -10, pident=0, length=0,
                                   bitscore=100, ppos=45, qcovs=0)

    blast_result_df['PRODUCT'] = blast_result_df.sseqid.str.replace('META\|', '')
    blast_result_df = blast_result_df.merge(gene_map_df, left_on='PRODUCT', right_on='PRODUCT')
    blast_result_df.sseqid = blast_result_df['UNIQUE-ID']

    draft_model = get_draft_from_template(MetaCyc_model, blast_result_df, remove_missing_genes=False)
    draft_model.id = 'draft model of %s from MetaCyc' % species_name
    cobra.io.save_json_model(draft_model, 'python/draft_%s_from_py.json' % species_name)

print('Done')
