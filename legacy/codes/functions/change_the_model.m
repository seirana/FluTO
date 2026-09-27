function mdl_  = change_the_model(mdl_, cb_rxns, cb)

%%% design our model
%%% block all carbon sources
for c = 1:size(cb_rxns,1)
    f = mdl_.rxnNumber == str2double(string(cb_rxns(c,1)));
    mdl_.lb(f) = 0;
    mdl_.ub(f) = 0;
end

%%% active one carbon source
f = find(mdl_.rxnNumber == str2double(string(cb_rxns(cb,1))));
mdl_.lb(f) = -10;
mdl_.ub(f) = -10;

%%% mdl_.rxns(200) = 'Phosphate exchange'
f = find(mdl_.rxnNumber == 200);
mdl_.lb(f) = -.9476;
mdl_.ub(f) = -.9476;

%%% llimit ATP synthase
f = find(mdl_.rxnNumber == 716);
mdl_.lb(f) = 3.15;
mdl_.ub(f) = 3.15;

%%% llimit ATP synthase
f = find(mdl_.rxnNumber == 721);
mdl_.lb(f) = .0931;
mdl_.ub(f) = .0931;

%%% llimit ATP synthase
f = find(mdl_.rxnNumber == 723);
mdl_.lb(f) = 12.612;
mdl_.ub(f) = 12.612;
end
