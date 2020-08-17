# issues and notes

## sequences_prossessing
1. species name list: data/initial_data/MostAbundant_SignificantSpecies.txt
- issue 1: wrong species name: eg. Bacteroides_dorei_vulgatus, that is the combination of two species name
    - current process: splitted to two species name:  Bacteroides dorei and Bacteroides vulgatus
    - comment:
- issue 2: some species with strain name(Clostridium_sp__HGF2) some without (Clostridium_spiroforme)
    - current process: if no strain name, representative strain selected by scripts; if with strain name, it is used as representative strain
    - comment:
- issue 3: duplicate eg. 'butyrate_producing_bacterium_SS3_4' and 'butyrate_producing_bacterium
    - current process: SS3_4 as representative strain
    - comment:

## draft GEMs
1. 




