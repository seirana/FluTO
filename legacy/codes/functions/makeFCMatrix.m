function fctable = makeFCMatrix(mdl_)
% For each pair of Ri ←→ Rj , they are fully coupled to each other if and only
% if the optimal values of
% maximize vi
% subject to v E C
% vj = 1,
% and
% minimize vi
% subject to v E C
% vj = 1,
% are equal.

rnd = 6;

fctable = zeros(size(mdl_.rxns,1));

for i = 1:size(mdl_.rxns,1)
    if mdl_.lb(i,1) ~= mdl_.ub(i,1)
        f = find(fctable(i,:) == 1);
        if size(f,2) == 0
            for j = i+1:size(mdl_.rxns,1)
                if mdl_.lb(j,1) ~= mdl_.ub(j,1) && ((mdl_.lb(j,1) == mdl_.lb(i,1) && mdl_.lb(i,1) == 0) || (mdl_.lb(j,1) > 0 && mdl_.lb(i,1) > 0))
                    f = find(fctable(j,:) == 1);
                    if size(f,2) == 0
                        %%% find flux cone
                        model = mdl_;
                        objective = model.rxns(i,1);
                        objectiveCoeff = 1;
                        ave = (model.lb(j,1) + model.ub(j,1))/2;
                        model.lb(j,1) = ave;
                        model.ub(j,1) = ave;
                        model = changeObjective(model, objective, objectiveCoeff); % change objective of the model
                        mini = optimizeCbModel(model, 'min'); % minmodelize the model
                        maxi = optimizeCbModel(model, 'max'); % maxmodelize the model
                        if round(mini.f,rnd) == round(maxi.f,rnd) && round(mini.f,rnd) > 0
                            %%% find flux cone
                            model = mdl_;
                            objective = model.rxns(i,1);
                            objectiveCoeff = 1;
                            ave = (model.lb(j,1) + model.ub(j,1))/2;
                            model.lb(j,1) = ave;
                            model.ub(j,1) = ave;
                            model = changeObjective(model, objective, objectiveCoeff); % change objective of the model
                            mini2 = optimizeCbModel(model, 'min'); % minmodelize the model
                            maxi2 = optimizeCbModel(model, 'max'); % maxmodelize the model
                            if round(mini2.f,rnd) == round(maxi2.f,rnd) && round(mini2.f,rnd) > 0
                                if round(mini2.f,rnd) > round(mini.f,rnd)
                                    fctable(i,j) = 1;
                                    fctable(j,i) = 1;
                                end
                            end
                        end
                    end
                end
                f = find(fctable(i,:) == 1);
                for m = 1:size(f,2)-1
                    for n = m+1:size(f,2)
                        fctable(f(m),f(n)) = 1;
                        fctable(f(n),f(m)) = 1;
                    end
                end
            end
        end
    end
end

end
