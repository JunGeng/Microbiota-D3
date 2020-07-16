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
headers = {'User-Agent': str(ua.chrome)}

host_link = 'https://metacyc.org'
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

for index_i in summary_df.index:
    substances_df_i = pd.read_csv('substances_df_' + str(index_i) + '_.tsv', sep='\t', index_col=0)
    substances_df_i = substances_df_i.dropna(how='all', axis=0)
    substances_df_i = substances_df_i.dropna(how='all', axis=1)

    substances_df = substances_df.merge(substances_df_i, how='outer', on=['Substances', 'Link', ])

    composition_df_i = pd.read_csv('composition_df_' + str(index_i) + '_.tsv', sep='\t', index_col=0)
    composition_df_i = composition_df_i.dropna(how='all', axis=0)
    composition_df_i = composition_df_i.dropna(how='all', axis=1)

    composition_df = composition_df.merge(composition_df_i, how='outer', on=['Constituents', 'Link'])

substances_df['metacyc_id'] = substances_df['Link'].str.split('=', expand=True)[2]
substances_df['class'] = substances_df['Link'].str.split('=', expand=True)[1]
substances_df['Link'] = host_link + substances_df['Link']
composition_df['metacyc_id'] = composition_df['Link'].str.split('=', expand=True)[2]
composition_df['class'] = composition_df['Link'].str.split('=', expand=True)[1]
composition_df['Link'] = host_link + composition_df['Link']

columns = list(composition_df.columns[[0, 1, -1, -2]]) + list(composition_df.columns[2:-2])
composition_df = composition_df[columns]

columns = list(substances_df.columns[[0, 1, -1, -2]]) + list(substances_df.columns[2:-2])
substances_df = substances_df[columns]

substances_df.to_csv('metacyc_growth_media_substances_df.tsv', sep='\t', )
composition_df.to_csv('metacyc_growth_media_composition_df.tsv', sep='\t', )

#
# url = 'https://metacyc.org/META/NEW-IMAGE?type=Growth-Media&object=MIX-13'
#
# response = requests.get(url)
# soup = BeautifulSoup(response.text, 'html.parser')
# tables = soup.findAll('table')
#
# if 'Substances' in str(tables[5]) :
#     records = get_records(tables[5])
#     columns = ['Substances','Link','Concentration' ,'Role' ]
#     recipe_substances_df_ = pd.DataFrame(data=records[1:],columns = columns )
# else:
#     for i in range(len(tables)):
#         print(i)
#
# if 'Constituents' in str(tables[6]) :
#     records = get_records(tables[6])
#     columns = ['Constituents','Link','Concentration'  ]
#     composition_df_i = pd.DataFrame(data=records[1:],columns = columns )
# else:
#     for i in range(len(tables)):
#         print(i)
