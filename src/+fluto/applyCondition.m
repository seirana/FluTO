function model = applyCondition(model, condition)
%APPLYCONDITION Apply an explicit FluTO condition specification.
%
%   CONDITION is a struct with optional fields:
%     blockReactionNumbers : reaction numbers to set to [0, 0]
%     activeReactionNumber : one reaction number to activate
%     activeBounds         : [lower upper] for the active reaction
%     fixedReactionNumbers : reaction numbers with fixed flux values
%     fixedFluxValues      : fixed values matching fixedReactionNumbers
%
%   This replaces machine- and model-specific index manipulation with a
%   validated configuration object.

    model = fluto.validateModel(model);
    if nargin < 2 || isempty(condition)
        return;
    end
    if ~isstruct(condition)
        error("fluto:InvalidCondition", "condition must be a struct.");
    end

    if isfield(condition, "blockReactionNumbers") && ~isempty(condition.blockReactionNumbers)
        blockNumbers = double(condition.blockReactionNumbers(:));
        indices = localReactionIndices(model, blockNumbers, "blockReactionNumbers");
        model.lb(indices) = 0;
        model.ub(indices) = 0;
    end

    if isfield(condition, "activeReactionNumber") && ~isempty(condition.activeReactionNumber)
        if ~isfield(condition, "activeBounds") || numel(condition.activeBounds) ~= 2
            error( ...
                "fluto:InvalidCondition", ...
                "activeBounds must be [lower upper] when activeReactionNumber is supplied.");
        end
        activeIndex = localReactionIndices( ...
            model, ...
            double(condition.activeReactionNumber), ...
            "activeReactionNumber");
        bounds = double(condition.activeBounds(:));
        if bounds(1) > bounds(2)
            error("fluto:InvalidCondition", "activeBounds lower value must be <= upper value.");
        end
        model.lb(activeIndex) = bounds(1);
        model.ub(activeIndex) = bounds(2);
    end

    hasFixedNumbers = isfield(condition, "fixedReactionNumbers") ...
        && ~isempty(condition.fixedReactionNumbers);
    hasFixedValues = isfield(condition, "fixedFluxValues") ...
        && ~isempty(condition.fixedFluxValues);

    if xor(hasFixedNumbers, hasFixedValues)
        error( ...
            "fluto:InvalidCondition", ...
            "fixedReactionNumbers and fixedFluxValues must be supplied together.");
    end

    if hasFixedNumbers
        numbers = double(condition.fixedReactionNumbers(:));
        values = double(condition.fixedFluxValues(:));
        if numel(numbers) ~= numel(values)
            error( ...
                "fluto:InvalidCondition", ...
                "fixedReactionNumbers and fixedFluxValues must have equal length.");
        end
        indices = localReactionIndices(model, numbers, "fixedReactionNumbers");
        model.lb(indices) = values;
        model.ub(indices) = values;
    end
end

function indices = localReactionIndices(model, numbers, fieldName)
    [isFound, indices] = ismember(numbers, model.rxnNumber);
    if any(~isFound)
        missing = numbers(~isFound);
        error( ...
            "fluto:UnknownReactionNumber", ...
            "%s contains unknown reaction numbers: %s", ...
            fieldName, ...
            strjoin(string(missing(:)'), ", "));
    end
end
