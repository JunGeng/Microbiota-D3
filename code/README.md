# code/


## sequences processing/
download  sequences
- step1_get_seq_info.py
  - get info, get summary

- step2_download_seq.py
  - download seqs
  - not necessary! replaced by step3_select_stran

- step3_select_stran.py
  - sort by: 'refseq_category' > 'relation_to_type_material' > 'assembly_level' > 'seq_rel_date'
    - 'refseq_category': 'reference genome' > 'representative genome'
    - 'assembly from type material'
    - 'Complete Genome' > 'Chromosome' > 'Scaffold' > 'Contig'
    - 'seq_rel_date': ascendingbool = False

  - after select the 'refseq_category' and 'relation_to_type_material', 115 species 'type' strain were selected
  - organize seq files. copy selected files to genomic_sequences/ and protein_sequences/


## draft GEMs/
- step1_draftGEMs_from_metacyc.m
  - get draft GEM from MetaCyc by RAVEN


## biomass/
- dna.py and rna.py functions to generate nda and rna coefficients adopt from: [BOFdat: Generating biomass objective functions for genome-scale metabolic models from experimental data
](https://github.com/jclachance/BOFdat)
- 

