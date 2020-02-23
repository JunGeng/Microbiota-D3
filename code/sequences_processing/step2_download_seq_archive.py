#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 2/13/20

"""step1_get_seq_info.py
:description : script to get the seq info
:param : input file: 'species_info.txt'
:returns: sequences: .faa and .fna
:rtype: 
"""

import os
import time
import pandas as pd
import wget

os.chdir('../../data/sequences_processing/')

# load data
refseq_info_df = pd.read_csv('species_info.txt', sep='\t')

# %% download:
# list_old = set(os.listdir('sequences2/'))
list_old = []
list_final = set(os.listdir('sequences/'))


f = open('Error_final_1.txt', 'w')
for index in range(0, refseq_info_df.shape[0]):  #
    # refseq_info_df.columns
    # refseq_info_df['# assembly_accession']
    # refseq_info_df.ftp_path

    # ignore the species that have too many seqs/strains: Escherichia coli (18000+); Klebsiella pneumoniae (8000+)
    if ('Escherichia coli' in refseq_info_df.iloc[index]['organism_name'] )or 'Klebsiella pneumoniae' in refseq_info_df.iloc[index]['organism_name']:
        print(index, end=' passed \t')
        if index % 5 == 0:
            print('')
        continue

    ftp_path = refseq_info_df.iloc[index]['ftp_path']
    faa_name = ftp_path.split('/')[-1] + '_protein.faa.gz'
    fna_name = ftp_path.split('/')[-1] + '_genomic.fna.gz'

    if faa_name in list_final and fna_name in list_final:
        print(index, end=' downloaded \t')
        if index % 5 == 0:
            print('')
        continue

    # strain_name = refseq_info_df.iloc[index]['organism_name']
    # strain_name = refseq_info_df.iloc[index]['infraspecific_name'].split('=')[-1]
    try:
        if faa_name not in list_old:
            file_url_faa = ftp_path + '/' + faa_name
            wget.download(file_url_faa, 'sequences/')
        else:
            os.system('mv sequences_process2/' + faa_name + ' sequences_process/' + faa_name)

        file_url_fna = ftp_path + '/' + fna_name
        wget.download(file_url_fna, 'sequences/')
        # os.rename('sequience/'+ftp_path.split('/')[-1] + '_genomic.gbff', os.path.join('all_Lreu_strains_gbff/', strain_name+ '.gbff'))
        # os.remove(my_file)
        print(index, end='\t')
        if index % 5 == 0:
            print('')
    except:
        #raise
        f.write('ERROR!!!\t' + str(index) + '\t' +ftp_path)
        f.flush()
        print('ERROR!!!\t' + str(index) + '\t' +ftp_path)

    if (index ) % 50 == 0:
        print('sleeping')
        time.sleep(0)
        # if (index ) % 100 == 0:
        #     print('sleeping')
        #     time.sleep(129)

f.close()
print('Done')
# all: 32466
# 818
# 819
# 19114