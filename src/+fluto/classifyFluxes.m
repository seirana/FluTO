function [model, fluxTable] = classifyFluxes(model, tolerance)
%CLASSIFYFLUXES Classify reactions from their feasible flux ranges.
%
%   Categories follow the terminology used by the FluTO implementation:
%     fixed      - lower and upper ranges are equal within tolerance
%     variable   - fixed-sign variable reaction (non-negative orientation)
%     reversible - sign-variable reaction crossing zero
%
%   Call canonicalizeReactionDirections first when negative irreversible
%   reactions are present.

    if nargin < 2 || isempty(tolerance)
        tolerance = 1e-9;
    end
    validateattributes(tolerance, {'numeric'}, {'scalar', 'real', 'nonnegative', 'finite'});

    model = fluto.validateModel(model);

    lower = model.lb;
    upper = model.ub;
    fixedMask = abs(upper - lower) <= tolerance;
    reversibleMask = lower < -tolerance & upper > tolerance;
    variableMask = ~fixedMask & lower >= -tolerance & upper > tolerance;
    negativeVariableMask = ~fixedMask & upper <= tolerance & lower < -tolerance;

    if any(negativeVariableMask)
        error( ...
            "fluto:NegativeIrreversibleReaction", ...
            ["Model contains negative irreversible reactions. " ...
             "Call fluto.canonicalizeReactionDirections before classification."]);
    end

    reactionType = strings(numel(model.rxns), 1);
    reactionType(fixedMask) = "fixed";
    reactionType(reversibleMask) = "reversible";
    reactionType(variableMask) = "variable";

    unclassified = reactionType == "";
    if any(unclassified)
        indices = find(unclassified);
        error( ...
            "fluto:UnclassifiedFluxRange", ...
            "Could not classify reaction indices: %s", ...
            strjoin(string(indices(:)'), ", "));
    end

    model.rxnType = cellstr(reactionType);
    fluxTable = table( ...
        model.rxnNumber, ...
        string(model.rxns), ...
        lower, ...
        upper, ...
        reactionType, ...
        "VariableNames", ...
        {"ReactionNumber", "ReactionID", "LowerBound", "UpperBound", "FluxClass"});
end
