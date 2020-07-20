function [reduceModel,biomassEqu]=ConstructBiomass(model,BM_COMP)
%% this function is to reconstruct Biomass equation based on multiple synthesis reactions of biomass components, i.e., protein, DNA, RNA, Lipids, CellWall and Pool
% eg,  BM_COMP={'Protein' 'DNA' 'RNA' 'Lipid' 'Cellwall' 'Pool'};
% output: biomass equation, and new model with one lumped biomass equation
%% Jun Geng (gejun@chalmers.se), 2020.07.10

% BM_COMP=strcat('1g',BM_COMP);
BM_COMP_index=getIndexes(model,BM_COMP,'mets');
% model.mets(BM_COMP_index)
BM_index=getIndexes(model,'Biomass','mets');
model.metNames(BM_index)
model.mets(BM_COMP_index)
BMRxn_index=find(model.S(BM_index,:)==1);
model.rxns(BMRxn_index)
constructEquations(model,BMRxn_index)
Syn_rxn=cellfun(@(x) find(model.S(x,:)==1),num2cell(BM_COMP_index),'un',0) %% look for the coefficient of protein component
model.rxnNames(cell2mat(Syn_rxn))
  biomass_coef=zeros(length(model.mets),1);
for i=1:length(Syn_rxn)
  COMP_coef=model.S(BM_COMP_index(i),BMRxn_index);
  COMP_synRxn=Syn_rxn{i};
  model.rxnNames(COMP_synRxn)
  constructEquations(model,COMP_synRxn)
  COMP_metIndex=BM_COMP_index(i)
  COMP_COMPIndex=find(model.S(:,COMP_synRxn)<0)
  model.mets(COMP_COMPIndex)
  COMP_SynCoef=model.S(COMP_COMPIndex,COMP_synRxn).*COMP_coef.*-1
  biomass_coef_i=zeros(length(model.mets),1);
  biomass_coef_i(COMP_COMPIndex)=COMP_SynCoef;
  biomass_coef=biomass_coef+biomass_coef_i;
  model.S(COMP_metIndex,:)=zeros(1,length(model.rxns));
  model.S(:,COMP_synRxn)=zeros(length(model.mets),1);
end
  model.S(:,BMRxn_index)=biomass_coef;
  model.S(BM_index,BMRxn_index)=1;

%   model.mets(find(model.S(:,BMRxn_index)<0))
%% remove original biomass component synthesis reactions
  constructEquations(model,BMRxn_index)
  constructEquations(iMH_Raven,BMRxn_index)
  reduceModel=removeReactionsFull(model,model.rxns(cell2mat(Syn_rxn)));
  reduceModel=removeMets(reduceModel,model.mets(BM_COMP_index));
  exportToExcelFormat(reduceModel,[iMH_path,'iMH551_constructBiomass.xlsx']);
  BM_index=getIndexes(reduceModel,'Biomass','mets');
reduceModel.metNames(BM_index)
BMRxn_index=find(reduceModel.S(BM_index,:)==1);
   biomassEqu=constructEquations(reduceModel,BMRxn_index);