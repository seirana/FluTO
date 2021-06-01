function mdl_ = original_model(adrs, file_name)
%%% read the mdl_
mdl_ = readCbModel(strcat(adrs,file_name, '.mat'));

%%% remove non-used fields from the mdl_CVN
mdl_ = rmfield(mdl_,'C');
mdl_ = rmfield(mdl_,'ctrs');
mdl_ = rmfield(mdl_,'d');
mdl_ = rmfield(mdl_,'dsense');

%%% add rxn & metabolite numbers and rxn type to the mdl_
mdl_.rxnNumber = zeros(size(mdl_.rxns,1),1);
for i = 1:size(mdl_.rxns,1)
    mdl_.rxnNumber(i,1) = i;
end
mdl_.metNumber = zeros(size(mdl_.mets,1),1);
for i = 1:size(mdl_.mets,1)
    mdl_.metNumber(i,1) = i;
end

mdl_.rxnType = cell(size(mdl_.rxns,1),1);

%%% remove 'core' biomass from network
f = find(string(mdl_.rxns) == 'BIOMASS_Ec_iJO1366_core_53p95M');
mdl_.lb(f) = 0;
mdl_.ub(f) = 0;
mdl_.c(f) = 0;

%%% set biomass
f = string(mdl_.rxns) == 'BIOMASS_Ec_iJO1366_WT_53p95M';
mdl_.c(f) = 1;
end
