clear all
% files = gunzip('data/sequences_processing/protein_sequences/*.gz')
seqs_list = dir('data/sequences_processing/protein_sequences/*.faa');

for seq_index = 1:length(seqs_list)
    seq_name = seqs_list(seq_index).name;
    species_name = strsplit(seq_name,'_protein.faa');
    species_name = species_name{1};
    model_name = strcat(species_name,'_Metacyc');
    seq_index
    model_name
    
    MetaCycDraftModel_i=getMetaCycModelForOrganism(model_name,strcat(...
             'data/sequences_processing/protein_sequences/',seq_name),1);
         
    exportModel(MetaCycDraftModel_i,strcat(...
        'data/draft_GEMs/',model_name),true,false);
        
end
