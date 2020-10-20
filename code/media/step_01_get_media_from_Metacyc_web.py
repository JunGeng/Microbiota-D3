#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 7/9/20

"""get_media_from_Metacyc_web.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import os
import time

import pandas as pd
import requests
from bs4 import BeautifulSoup
from fake_useragent import UserAgent

os.chdir('../../data/media/')


def get_records(table):
    records = []
    for tr in table.findAll("tr"):
        trs = tr.findAll("td")
        record = []
        for each in trs:
            text = each.text.replace('\n', '')
            record.append(text)
            try:
                link = each.find('a')['href']
                record.append(link)
                # print(each)
            except:
                pass
        records.append(record)
    return records


# %% summary page
# url = 'https://metacyc.org/META/NEW-IMAGE?object=Growth-Media'
#
# tables = pd.read_html(url)
# summary_df = tables[7]          # find the table index is 7
# new_header = summary_df.iloc[0] #grab the first row for the header
# summary_df = summary_df[1:] #take the data less the header row
# summary_df.columns = new_header #set the header row as the df header
#
# response = requests.get(url)
# soup = BeautifulSoup(response.text, 'html.parser')
# # print(soup.prettify())
# table = soup.findAll('table')
#
# Link = []
# for tr in table[7].findAll("tr"):
#     trs = tr.findAll("td")
#     for each in trs:
#         try:
#             link = each.find('a')['href']
#             print(link)
#             if '/META/NEW-IMAGE?type=Growth-Media&object'in link:
#                 Link.append(link)
#         except:
#             pass
# summary_df['Link'] = Link
# # summary_df.to_csv('metacyc_growth_media_list.tsv',sep = '\t',)

# %%  media_i page
summary_df = pd.read_csv('metacyc_growth_media_list.tsv', sep='\t', index_col=0)

sub_columns = ['Substances', 'Link', 'Concentration', 'Role']
substances_df = pd.DataFrame(columns=sub_columns)

com_columns = ['Constituents', 'Link', 'Concentration']
composition_df = pd.DataFrame(columns=com_columns)  # 'Constituents'
headers = {
    'User-Agent': 'Mozilla/5.0 (Macintosh; Intel Mac OS X 10_10_1) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/39.0.2171.95 Safari/537.36'}
ua = UserAgent()
print(ua.chrome)
host_link = 'https://metacyc.org'
headers = {'User-Agent': str(ua.chrome)}
download = False
if download:

    for index_i in summary_df.index[1:]:  # summary_df.index:
        # index_i = 1
        print(index_i)
        name_i = summary_df['Medium Name'].loc[index_i]
        url_i = host_link + summary_df['Link'].loc[index_i]

        response = requests.get(url_i, headers=headers)
        soup = BeautifulSoup(response.text, 'html.parser')
        tables = soup.findAll('table')

        if 'Substances' in str(tables[5]):
            records = get_records(tables[5])
            columns = ['Substances', 'Link', 'Concentration_' + name_i, 'Role_' + name_i]
            for i in range(len(records)):
                record = records[i]
                if len(record) == 4:
                    continue
                elif len(record) < 4:
                    records[i] = record + [''] * (4 - len(record))
                elif len(record) > 4:
                    print(i)
            substances_df_i = pd.DataFrame(data=records[1:], columns=columns)
            substances_df_i.to_csv('substances_df_' + str(index_i) + '_.tsv', sep='\t', )
            substances_df = substances_df.merge(substances_df_i, how='outer', on=['Substances', 'Link'])
        else:
            print(name_i, index_i, 'Substances', 'failure')

        if 'Constituents' in str(tables[6]):
            records = get_records(tables[6])
            columns = ['Constituents', 'Link', 'Concentration_' + name_i, 'Role_' + name_i]
            for i in range(len(records)):
                record = records[i]
                if len(record) == 3:
                    continue
                elif len(record) < 3:
                    records[i] = record + [''] * (4 - len(record))
                elif len(record) > 3:
                    print(i)
            composition_df_i = pd.DataFrame(data=records[1:], columns=columns)
            composition_df_i.to_csv('composition_df_' + str(index_i) + '_.tsv', sep='\t', )
            composition_df = composition_df.merge(composition_df_i, how='outer', on=['Constituents', 'Link'])

        else:
            print(name_i, index_i, 'Constituents', 'failure')
        time.sleep(200)

# %% due to the anti-reptile of metacyc, the scrips is unstable
sub_columns = ['Substances', 'Link', ]
substances_df = pd.DataFrame(columns=sub_columns)

com_columns = ['Constituents', 'Link', ]
composition_df = pd.DataFrame(columns=com_columns)  # 'Constituents'
media_name_dic = \
    {
     'Concentration_dGMM': 'M1',
     'Concentration_LAB': 'M2',
     'Concentration_dGMM+LAB': 'M3',
     'Concentration_dGMM+LAB low M/V': 'M4',
     'Concentration_dGMM+LAB exclude SCFA': 'M5',
     'Concentration_dGMM+LAB only monosacharides': 'M7',
     'Concentration_dGMM+LAB plus Mucin': 'M8',
     'Concentration_dGMM+LAB only Mucin': 'M9',
     'Concentration_dGMM+LAB 10% aminoacids': 'M10',
     'Concentration_dGMM+LAB excluding aromatic AA': 'M11',
     'Concentration_B.thetaiotaomicron MM': 'M13',
     'Concentration_C.perfiringens MM': 'M14',
     'Concentration_E.coli MM1': 'M15A',
     'Concentration_E.coli MM2': 'M15B',
     'Concentration_V.parvula defined medium': 'M16',
     'Concentration_GMM': 'GMM',
     'Concentration_mGAM': 'mGAM',
     'Concentration_brain heart infusion': 'BHI++',
     'Concentration_WCA': 'WCA'}
media_name_list = []
for index_i in summary_df.index:
    substances_df_i = pd.read_csv('substances_df_' + str(index_i) + '_.tsv', sep='\t', index_col=0)
    substances_df_i.columns = [i.replace(' ', '') for i in substances_df_i.columns[0:2]] + [i for i in
                                                                                            substances_df_i.columns[2:]]
    substances_df_i = substances_df_i.dropna(how='all', axis=0)
    substances_df_i = substances_df_i.dropna(how='all', axis=1)
    substances_df_i['Substances'] = substances_df_i['Substances'].str.replace(' ', '')
    substances_df_i['Link'] = substances_df_i['Link'].str.replace(' ', '')
    media_columns = substances_df_i.columns[2]
    media_columns = media_name_dic[media_columns] + ': ' + media_columns.replace('Concentration_', '')
    # media_name_list.append(media_columns)
    substances_df_i[media_columns] = substances_df_i[substances_df_i.columns[2]]
    substances_df = substances_df.merge(substances_df_i, how='outer', on=['Substances', 'Link'])

    composition_df_i = pd.read_csv('composition_df_' + str(index_i) + '_.tsv', sep='\t', index_col=0)
    composition_df_i.columns = [i.replace(' ', '') for i in composition_df_i.columns[0:2]] + [i for i in
                                                                                              composition_df_i.columns[
                                                                                              2:]]

    composition_df_i = composition_df_i.dropna(how='all', axis=0)
    composition_df_i = composition_df_i.dropna(how='all', axis=1)
    composition_df_i['Constituents'] = composition_df_i['Constituents'].str.replace(' ', '')
    composition_df_i['Link'] = composition_df_i['Link'].str.replace(' ', '')
    media_columns = composition_df_i.columns[2]
    media_columns = media_name_dic[media_columns] + ': ' + media_columns.replace('Concentration_', '')
    media_name_list.append(media_columns)
    composition_df_i[media_columns] = composition_df_i[composition_df_i.columns[2]]
    composition_df_i.loc[composition_df_i['Constituents'] == 'fructose', media_columns] = '5.5 mM'  # (1/180*1000)
    composition_df_i.loc[composition_df_i['Constituents'] == 'lactose', media_columns] = '2.92 mM'  # (1/180*1000)
    composition_df_i.loc[composition_df_i['Constituents'] == 'maltose', media_columns] = '2.92 mM'  # (1/180*1000)

    composition_df = composition_df.merge(composition_df_i, how='outer', on=['Constituents', 'Link'])

substances_df['metacyc_id'] = substances_df['Link'].str.split('=', expand=True)[2]
substances_df['class'] = substances_df['Link'].str.split('=', expand=True)[1]
substances_df['Link'] = host_link + substances_df['Link']
composition_df['metacyc_id'] = composition_df['Link'].str.split('=', expand=True)[2]
composition_df['metacyc_id'] = composition_df['metacyc_id'].str.replace('%2b', '+')
composition_df['class'] = composition_df['Link'].str.split('=', expand=True)[1]
composition_df['Link'] = host_link + composition_df['Link']

# columns = list(composition_df.columns[[0, 1, -1, -2]]) + list(composition_df.columns[2:-2])
columns = ['Constituents', 'Link', 'class', 'metacyc_id'] + list(set(media_name_list))
composition_df = composition_df[columns]
composition_df = composition_df.sort_values(by=['class', 'Constituents'])

# columns = list(substances_df.columns[[0, 1, -1, -2]]) + list(substances_df.columns[2:-2])
columns = ['Substances', 'Link', 'class', 'metacyc_id'] + list(set(media_name_list))
substances_df = substances_df[columns]

substances_df.to_csv('metacyc_growth_media_substances_df.tsv', sep='\t', )
composition_df.to_csv('metacyc_growth_media_composition_df.tsv', sep='\t', )
