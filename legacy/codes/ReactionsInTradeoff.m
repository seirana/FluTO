clc;clear;

addpath(genpath('./cobratoolbox'));
addpath(genpath('./Trade-offs/functions'));

changeCobraSolver('glpk', 'LP');
solver = 'glpk';

adrs = './Trade-offs/';
file_name = 'Ecoli_iJO1366';

%%% read carbon source list and read the model
[~,~,cb_rxns] = xlsread(strcat(adrs, file_name, '_CarbonSourceList.xlsx'));
cb_rxns = cb_rxns(3:end,:);

%%%  read the original model(the mrtabolic network)
model = original_model(adrs, file_name);

%%% define the variables
mdl_file =  cell(size(model.rxns,1),3+size(cb_rxns,1)*3);
%%% define a cell to keep results for flux ratio under different conditions
for i = 1:size(model.rxns,1)
    mdl_file(i,1:3) = {i model.rxns(i) model.subSystems(i,1)};
end

%%% do it for different biomasses
for cb = 1:size(cb_rxns,1)
    rxns_flux = cell(size(model.rxns,1),3);    
    mdl0 = model;
    mdl1 = change_the_model(mdl0, cb_rxns, cb); %%% change the model and remove blocked rxns and mets
    mdl2 = FluxCone(mdl1); %%% flux cone and flux classification
    [mdl, rxns_flux] = FluxClasification(mdl2, rxns_flux); %%% flux classification
    fctable = makeFCMatrix(mdl2); %%% make F2C2Matrix

    mdl_file(:,1+cb*3:3+cb*3) = rxns_flux(:,1:3);   
 
    tradeoffs(mdl, fctable, cb, adrs, file_name); %%% find_tradeoffs
end

T = cell2table(mdl_file);
writetable(T, strcat(adrs, file_name, '_flux_', string(cb), '.xlsx'));
