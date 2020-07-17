# data/
data

## sequences processing/

download  sequences
- species_name_list.csv: manualy curated species name list
- species_info.txt: information
- species_summary:
- species_final: the final list (selected strain list)

## initial data/
- MostAbundant_SignificantSpecies.txt
  - From Jun (Box)
- data_from_database/assembly_summary_refseq.txt and assembly_summary_refseq_historical.txt
  - From NCBI FTP [ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/](ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/)
- data_from_database/MetaCyc_23.5
  - [MetaCyc](https://metacyc.org/) database: http://bioinformatics.ai.sri.com/ecocyc/dist/flatfiles-52983746/
  - MetaCyc_cobra.xml modifed by Cobra toolbox

### genomic sequences/
- genomic sequences files

### protein sequences/
- protein sequences files

### sequences/
- temp dir to download data
- .faa and .fna files from NCBI

## draft GEMs/

### ./
- draft models from RAVEN
### blast
- diamond blast db(.dmnd) and results

### python
- draft models from python scrips


## biomass/
- Biomass_compare_summary.xlsx the summary.xlsx
 - including ID mapping information
- Biomass_Literatures_Templates_Combined_new_correct_upload.xlsx
  -  Record the comparison results of multiple GEMs biomass equation from different sources of gut microbiota specific species
  - mat GEMs from: [Gut microbiota dysbiosis is associated with malnutrition and reduced plasma amino acid levels: Lessons from genome-scale metabolic modeling](https://doi.org/10.1016/j.ymben.2018.07.018)
  - xml GEMs from: [Quantifying Diet-Induced Metabolic Changes of the Human Gut Microbiome](http://dx.doi.org/10.1016/j.cmet.2015.07.001)
  - iMSU GEM for Strephococcus mutans: DOI: 10.1128/mSystems.00529-19
  - iAF692 is selected as the biomass template for Archea: doi: 10.1155/2017/9763848

## media/
- mediaComposition.tsv from [Nutritional preferences of human gut bacteria reveal their metabolic idiosyncrasies](https://www.nature.com/articles/s41564-018-0123-9),[github](https://github.com/sandrejev/growth_curves/blob/master/data/compounds.xlsx)
- mediaComposition.xlsx sheet M2 from [MetaCyc Growth Medium: dGMM+LAB](https://metacyc.org/META/NEW-IMAGE?type=Growth-Media&object=MIX-13)

