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

os.chdir('../../data/draft_GEMs/')


def blastp_pairwise(qseq_file, sseq_file, out_dir='', diamond=True):
    '''

    :param qseq_file:
    :param sseq_file:
    :param out_dir:
    :return:  two files
    '''

    os.system('mkdir blast_tmps')

    if diamond:
        mk_db = 'diamond makedb --in ' + qseq_file + ' -d blast_tmps/qseq_db \n' \
                                                     'diamond makedb --in ' + sseq_file + ' -d blast_tmps/sseq_db '
        os.system(mk_db)
        print('diamond blasting...')
        options = ' --top 10 --more-sensitive '
        # options = ''
        diamond_blastp_cmd = 'diamond blastp -d blast_tmps/qseq_db -q ' + sseq_file + options + ' -o ' + out_dir + \
                             'blast_result_s_in_q.csv --evalue 1 --outfmt 6 qseqid sseqid evalue pident length bitscore ppos qcovhsp \n' \
                             'diamond blastp -d blast_tmps/sseq_db -q ' + qseq_file + options + ' -o ' + out_dir + \
                             'blast_result_q_in_s.csv --evalue 1 --outfmt 6 qseqid sseqid evalue pident length bitscore ppos qcovhsp'
        os.system(diamond_blastp_cmd)

    else:
        mk_db = 'makeblastdb -in ' + qseq_file + ' -dbtype prot -out blast_tmps/qseq_db -parse_seqids\n' \
                                                 'makeblastdb -in ' + sseq_file + ' -dbtype prot -out blast_tmps/sseq_db -parse_seqids'
        os.system(mk_db)
        print('blasting...')
        blastp_cmd = 'blastp -db blast_tmps/qseq_db -query ' + sseq_file + ' -out ' + out_dir + 'blast_result_s_in_q.csv -evalue 1 -outfmt "6 qseqid sseqid evalue pident length bitscore ppos qcovs" \n' \
                                                                                                'blastp -db blast_tmps/sseq_db -query ' + qseq_file + ' -out ' + out_dir + 'blast_result_q_in_s.csv -evalue 1 -outfmt "6 qseqid sseqid evalue pident length bitscore ppos qcovs"'
        os.system(blastp_cmd)

    os.system('rm -r blast_tmps')
    print('out put files are ' + out_dir + 'blast_result_s_in_q.csv and blast_result_q_in_s.csv')
    return out_dir + 'blast_result_s_in_q.csv', out_dir + 'blast_result_q_in_s.csv'


# %% <diamond makedb> our species query sequence and subject sequence in the MetaCyc.
print('diamond makedb...')
# make dimon database for MetaCyc
sseq_file = '../initial_data/data_from_database/MetaCyc_23.5/protseq.fsa'
MetaCyc_db_name = 'blast/MetaCyc'
# mk_db = 'diamond makedb --in '+ sseq_file +' -d '+ MetaCyc_db_name
# os.system(mk_db)

# make dimon database for each species
# qseq_list = os.listdir('../sequences_processing/protein_sequences/')
qseq_list = [i for i in os.listdir('../sequences_processing/protein_sequences/') if i.endswith('.faa')]
species_name_list = [i.replace('_protein.faa', '') for i in qseq_list]

# for qseq_file,species_name in zip(qseq_list,species_name_list):
#     mk_db = 'diamond makedb --in ../sequences_processing/protein_sequences/' + qseq_file +\
#             ' -d blast/'+ species_name
#     os.system(mk_db)

# %% <diamond blasting>
print('diamond blasting...')
options = ' --top 10 --more-sensitive '
for qseq_file, species_name in zip(qseq_list, species_name_list):
    diamond_blastp_cmd_1 = 'diamond blastp -d blast/MetaCyc -q ../sequences_processing/protein_sequences/%s ' \
                           '%s -o blast/blast_result_q_in_MetaCyc_%s.csv' \
                           ' --evalue 1 --outfmt 6 qseqid sseqid evalue pident length bitscore ppos qcovhsp' \
                           % (qseq_file, options, species_name)
    os.system(diamond_blastp_cmd_1)

    diamond_blastp_cmd_2 = 'diamond blastp -d blast/%s -q ../initial_data/data_from_database/MetaCyc_23.5/protseq.fsa ' \
                           '%s -o blast/blast_result_MetaCyc_in_q_%s.csv' \
                           ' --evalue 1 --outfmt 6 qseqid sseqid evalue pident length bitscore ppos qcovhsp' \
                           % (species_name, options, species_name)
    os.system(diamond_blastp_cmd_2)

# %% <select blust result by cutoff>
