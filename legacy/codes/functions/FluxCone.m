%%% flux in entire cone
function mdl_ = FluxCone(mdl_)
tmp = zeros(size(mdl_.rxns,1),2);
for i = 1:size(mdl_.rxns,1)
    objective = mdl_.rxns(i);
    objectiveCoeff = 1;
    mdl_ = changeObjective(mdl_, objective, objectiveCoeff); % change objective of the model
    maxi = optimizeCbModel(mdl_, 'max'); % maximize the model
    mini = optimizeCbModel(mdl_, 'min'); % minimize the model
    tmp(i,1) = round(mini.f,5);
    tmp(i,2) = round(maxi.f,5);
end

for i = 1:size(mdl_.rxns,1)
    mdl_.lb(i) = tmp(i,1);
    mdl_.ub(i) = tmp(i,2);
end

%%% correct model, remove blocked reactions and metaolites, zero columns
mdl_ = removeBlockedReactionsAndMetabolites(mdl_);
end
