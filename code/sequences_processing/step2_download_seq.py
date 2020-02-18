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

import pandas as pd
import wget

os.chdir('../../data/sequences_processing/')

# load data
refseq_info_df = pd.read_csv('species_info.txt', sep='\t')

# %% download:
# list_1 = os.listdir('sequences_process2/')
list_1 = []

for index in range(0, refseq_info_df.shape[0]):  #
    # refseq_info_df.columns
    # refseq_info_df['# assembly_accession']
    # refseq_info_df.ftp_path

    ftp_path = refseq_info_df.iloc[index]['ftp_path']
    faa_name = ftp_path.split('/')[-1] + '_protein.faa.gz'
    fna_name = ftp_path.split('/')[-1] + '_genomic.fna.gz'

    # strain_name = refseq_info_df.iloc[index]['organism_name']
    # strain_name = refseq_info_df.iloc[index]['infraspecific_name'].split('=')[-1]
    try:
        if faa_name not in list_1:
            file_url_faa = ftp_path + '/' + faa_name
            wget.download(file_url_faa, 'sequences/')
        else:
            os.system('mv sequences_process2/' + faa_name + ' sequences_process/' + faa_name)

        file_url_fna = ftp_path + '/' + fna_name
        wget.download(file_url_fna, 'sequences/')
        # os.rename('sequience/'+ftp_path.split('/')[-1] + '_genomic.gbff', os.path.join('all_Lreu_strains_gbff/', strain_name+ '.gbff'))
        # os.remove(my_file)
        print(index)
    except:
        raise
        print('ERROR!!!   ' + file_url)
