function [model, report] = removeBlockedReactions(model, tolerance)
%REMOVEBLOCKEDREACTIONS Remove structurally empty and zero-range reactions.
%
%   This function delegates model-field bookkeeping to COBRA Toolbox's
%   removeRxns when available, then restores FluTO's stable reaction numbers.
%   Zero-row metabolites are removed with removeMetabolites when that helper
%   is available; otherwise they are retained because they do not change the
%   feasible flux space.

    if nargin < 2 || isempty(tolerance)
        tolerance = 1e-9;
    end
    validateattributes(tolerance, {'numeric'}, {'scalar', 'real', 'nonnegative', 'finite'});

    model = fluto.validateModel(model);
    originalRxns = string(model.rxns);
    originalNumbers = model.rxnNumber;

    zeroColumn = full(sum(abs(model.S), 1))' <= tolerance;
    zeroRange = abs(model.lb) <= tolerance & abs(model.ub) <= tolerance;
    blockedMask = zeroColumn | zeroRange;
    blockedRxns = model.rxns(blockedMask);

    if any(blockedMask)
        if exist("removeRxns", "file") ~= 2
            error( ...
                "fluto:CobraDependencyMissing", ...
                ["COBRA Toolbox function removeRxns was not found. " ...
                 "Initialize COBRA Toolbox before removing blocked reactions."]);
        end
        model = removeRxns(model, blockedRxns);
    end

    model = fluto.validateModel(model);

    [found, oldIndex] = ismember(string(model.rxns), originalRxns);
    if ~all(found)
        error( ...
            "fluto:ReactionIdentityLost", ...
            "Could not map reduced model reactions back to original reaction IDs.");
    end
    model.rxnNumber = originalNumbers(oldIndex);

    zeroMetaboliteMask = full(sum(abs(model.S), 2)) <= tolerance;
    removedMetabolites = strings(0, 1);

    if any(zeroMetaboliteMask) && exist("removeMetabolites", "file") == 2
        removedMetabolites = string(model.mets(zeroMetaboliteMask));
        model = removeMetabolites( ...
            model, ...
            model.mets(zeroMetaboliteMask), ...
            false);
        model = fluto.validateModel(model);
    end

    report = struct( ...
        "blockedReactionIDs", string(blockedRxns), ...
        "nBlockedReactions", numel(blockedRxns), ...
        "removedMetaboliteIDs", removedMetabolites, ...
        "nRemovedMetabolites", numel(removedMetabolites));
end
