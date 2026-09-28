%%% remove blocked reactions and metabolites, zero columns and metabolites
function mdl_ = removeBlockedReactionsAndMetabolites(mdl_)

%%% change negative irreversible rxns to positive
for i = 1:size(mdl_.rxns,1)
    if mdl_.lb(i) < 0 && mdl_.ub(i) <= 0
        tmp = mdl_.lb(i);
        mdl_.lb(i) = -1*mdl_.ub(i);
        mdl_.ub(i) = -1*tmp;
        mdl_.S(:,i) = -1*mdl_.S(:,i);
    end
end

%%% remove blocked reactions from model
sz = size(mdl_.rxns,1);
for i = sz:-1:1
    sm = 0;
    for j = 1:size(mdl_.mets,1)
        sm = sm + abs(mdl_.S(j,i));
    end
    if sm == 0 || (mdl_.lb(i) == 0 && mdl_.ub(i) == 0)
        mdl_.S(:,i) = [];
        mdl_.rxns(i,:) = [];
        mdl_.lb(i,:) = [];
        mdl_.ub(i,:) = [];
        mdl_.c(i,:) = [];
        mdl_.rules(i,:) = [];
        mdl_.rxnNames(i,:) = [];
        mdl_.subSystems(i,:) = [];
        mdl_.rxnNumber(i,:) = [];
        mdl_.rxnType(i,:) = [];
    end
end

%%% remove blocked metabolites from model
sz = size(mdl_.mets,1);
for i = sz:-1:1
    sm = 0;
    for j = 1:size(mdl_.rxns,1)
        sm = sm + abs(mdl_.S(i,j));
    end
    if sm == 0
        mdl_.S(i,:) = [];
        mdl_.mets(i,:) = [];
        mdl_.b(i,:) = [];
        mdl_.csense(i,:) = [];
        mdl_.metFormulas(i,:) = [];
        mdl_.metNames(i,:) = [];
        mdl_.metNumber(i,:) = [];
    end
end
end