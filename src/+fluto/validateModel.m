function model = validateModel(model)
%VALIDATEMODEL Validate and normalize the COBRA-model fields used by FluTO.
%
%   MODEL = fluto.validateModel(MODEL) checks the dimensions and numerical
%   validity of the stoichiometric model and adds stable reaction/metabolite
%   indices when they are absent.
%
%   FluTO requires a COBRA-compatible model with fields S, rxns, mets, lb,
%   ub, and c. This function does not mutate the biochemical constraints.

    if ~isstruct(model)
        error("fluto:InvalidModel", "Model must be a struct.");
    end

    requiredFields = ["S", "rxns", "mets", "lb", "ub", "c"];
    missingFields = requiredFields(~isfield(model, cellstr(requiredFields)));
    if ~isempty(missingFields)
        error( ...
            "fluto:MissingModelFields", ...
            "Model is missing required fields: %s", ...
            strjoin(missingFields, ", "));
    end

    if ~isnumeric(model.S) && ~issparse(model.S)
        error("fluto:InvalidStoichiometry", "model.S must be numeric or sparse.");
    end
    if ndims(model.S) ~= 2
        error("fluto:InvalidStoichiometry", "model.S must be two-dimensional.");
    end
    if any(~isfinite(nonzeros(model.S)))
        error("fluto:InvalidStoichiometry", "model.S contains non-finite values.");
    end

    nMetabolites = size(model.S, 1);
    nReactions = size(model.S, 2);

    if numel(model.rxns) ~= nReactions
        error( ...
            "fluto:ReactionDimensionMismatch", ...
            "numel(model.rxns) must equal size(model.S, 2).");
    end
    if numel(model.mets) ~= nMetabolites
        error( ...
            "fluto:MetaboliteDimensionMismatch", ...
            "numel(model.mets) must equal size(model.S, 1).");
    end

    model.rxns = cellstr(string(model.rxns(:)));
    model.mets = cellstr(string(model.mets(:)));
    model.lb = double(model.lb(:));
    model.ub = double(model.ub(:));
    model.c = double(model.c(:));

    vectorFields = {"lb", "ub", "c"};
    for i = 1:numel(vectorFields)
        fieldName = vectorFields{i};
        if numel(model.(fieldName)) ~= nReactions
            error( ...
                "fluto:ReactionDimensionMismatch", ...
                "model.%s must contain one value per reaction.", ...
                fieldName);
        end
        if any(~isfinite(model.(fieldName)))
            error( ...
                "fluto:NonFiniteBounds", ...
                "model.%s contains non-finite values.", ...
                fieldName);
        end
    end

    if any(model.lb > model.ub)
        error( ...
            "fluto:InvalidBounds", ...
            "Every lower bound must be <= its upper bound.");
    end

    if ~isfield(model, "rxnNumber") || numel(model.rxnNumber) ~= nReactions
        model.rxnNumber = (1:nReactions)';
    else
        model.rxnNumber = double(model.rxnNumber(:));
    end

    if ~isfield(model, "metNumber") || numel(model.metNumber) ~= nMetabolites
        model.metNumber = (1:nMetabolites)';
    else
        model.metNumber = double(model.metNumber(:));
    end

    if ~isfield(model, "rxnType") || numel(model.rxnType) ~= nReactions
        model.rxnType = repmat({""}, nReactions, 1);
    else
        model.rxnType = cellstr(string(model.rxnType(:)));
    end

    if ~isfield(model, "subSystems") || numel(model.subSystems) ~= nReactions
        model.subSystems = repmat({""}, nReactions, 1);
    end
end
