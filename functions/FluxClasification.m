%%detect positive, negative and reversible rxns_flux
function [mdl_, rxns_flux] = FluxClasification(mdl_, rxns_flux)
for i = 1:size(mdl_.rxns,1)
    rxns_flux(mdl_.rxnNumber(i),1) = {mdl_.lb(i)};
    rxns_flux(mdl_.rxnNumber(i),2) = {mdl_.ub(i)};
    
    if mdl_.lb(i,1) >= 0 && mdl_.ub(i,1) > 0
        rxns_flux(mdl_.rxnNumber(i),3) = {'variable'};
        mdl_.rxnType(i) = {'variable'};
    end
    
    if  mdl_.lb(i,1) == mdl_.ub(i,1)
        rxns_flux(mdl_.rxnNumber(i),3) = {'fixed'};
        mdl_.rxnType(i) = {'fixed'};
    end
    
    if mdl_.lb(i,1) < 0 && mdl_.ub(i,1) > 0
        rxns_flux(mdl_.rxnNumber(i),3) = {'reversible'};
        mdl_.rxnType(i) = {'reversible'};
    end
end
end
