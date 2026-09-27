function result = enumerateTradeoffs(model, couplingMatrix, options)
%ENUMERATETRADEOFFS Enumerate FluTO trade-off supports with MILP.
%
%   RESULT = fluto.enumerateTradeoffs(MODEL, COUPLINGMATRIX, OPTIONS)
%   implements the published FluTO mixed-integer formulation with explicit
%   solver status handling and deterministic support de-duplication.
%
%   Required product: Optimization Toolbox (intlinprog).
%
%   OPTIONS fields:
%     minDegree     first trade-off degree to search (default 2)
%     maxDegree     largest trade-off degree to search (default 9)
%     bigM          historical MILP big-M value (default 1001)
%     variableBound absolute bound for k/alpha/u variables (default 100)
%     display       intlinprog display mode (default "off")
%     maxSolutions  optional global cap, Inf by default

    if nargin < 3 || isempty(options)
        options = struct();
    end
    options = localOptions(options);

    if exist("intlinprog", "file") ~= 2
        error( ...
            "fluto:OptimizationToolboxMissing", ...
            "Optimization Toolbox function intlinprog is required.");
    end

    model = fluto.validateModel(model);
    nReactions = numel(model.rxns);

    if ~isequal(size(couplingMatrix), [nReactions, nReactions])
        error( ...
            "fluto:InvalidCouplingMatrix", ...
            "couplingMatrix must be nReactions-by-nReactions.");
    end

    reactionType = string(model.rxnType(:));
    validTypes = ["fixed", "variable", "reversible"];
    if any(~ismember(reactionType, validTypes))
        error( ...
            "fluto:MissingFluxClassification", ...
            "Run fluto.classifyFluxes before trade-off enumeration.");
    end

    fixedIndices = find(reactionType == "fixed");
    variableIndices = find(reactionType == "variable");
    reversibleIndices = find(reactionType == "reversible");

    if isempty(fixedIndices)
        error( ...
            "fluto:NoFixedReactions", ...
            "Trade-off formulation requires at least one fixed reaction.");
    end
    if numel(variableIndices) < options.minDegree
        result = localEmptyResult(model, variableIndices, options);
        return;
    end

    nFixed = numel(fixedIndices);
    nVariable = numel(variableIndices);

    Nfixed = model.S(:, fixedIndices);
    Nvariable = model.S(:, variableIndices);
    Nreversible = model.S(:, reversibleIndices);
    N = [Nfixed, Nvariable, Nreversible];

    nMetabolites = size(N, 1);
    nFluxColumns = size(N, 2);

    % Decision vector:
    % [k(nMetabolites); alpha(nFixed+nVariable); u(nVariable); s(nVariable)]
    alphaCount = nFixed + nVariable;
    totalVariables = nMetabolites + alphaCount + nVariable + nVariable;

    alphaStart = nMetabolites + 1;
    variableAlphaStart = nMetabolites + nFixed + 1;
    uStart = nMetabolites + alphaCount + 1;
    sStart = uStart + nVariable;

    alphaBlock = sparse(nFluxColumns, alphaCount);
    alphaBlock(1:alphaCount, 1:alphaCount) = -speye(alphaCount);

    AeqBase = [ ...
        sparse(N'), ...
        alphaBlock, ...
        sparse(nFluxColumns, nVariable), ...
        sparse(nFluxColumns, nVariable)];
    beqBase = zeros(nFluxColumns, 1);

    Ain = sparse(0, totalVariables);
    bin = zeros(0, 1);

    % At least one fixed alpha coefficient contributes.
    row = sparse(1, totalVariables);
    row(alphaStart:(alphaStart + nFixed - 1)) = -1;
    Ain(end + 1, :) = row;
    bin(end + 1, 1) = -1;

    % |alpha_variable| <= u, alpha_variable <= 0 in the historical model.
    for i = 1:nVariable
        alphaIndex = variableAlphaStart + i - 1;
        uIndex = uStart + i - 1;
        sIndex = sStart + i - 1;

        row = sparse(1, totalVariables);
        row(alphaIndex) = 1;
        row(uIndex) = -1;
        Ain(end + 1, :) = row;
        bin(end + 1, 1) = 0;

        row = sparse(1, totalVariables);
        row(alphaIndex) = -1;
        row(uIndex) = -1;
        Ain(end + 1, :) = row;
        bin(end + 1, 1) = 0;

        row = sparse(1, totalVariables);
        row(alphaIndex) = 1;
        row(sIndex) = options.bigM;
        Ain(end + 1, :) = row;
        bin(end + 1, 1) = options.bigM - 1;

        row = sparse(1, totalVariables);
        row(alphaIndex) = -1;
        row(sIndex) = -options.bigM;
        Ain(end + 1, :) = row;
        bin(end + 1, 1) = 0;
    end

    variableCoupling = logical(couplingMatrix(variableIndices, variableIndices));

    % Fully-coupled alternatives cannot be selected simultaneously.
    [coupledLeft, coupledRight] = find(triu(variableCoupling, 1));
    for pairIndex = 1:numel(coupledLeft)
        row = sparse(1, totalVariables);
        row(sStart + coupledLeft(pairIndex) - 1) = 1;
        row(sStart + coupledRight(pairIndex) - 1) = 1;
        Ain(end + 1, :) = row;
        bin(end + 1, 1) = 1;
    end

    lb = [ ...
        -options.variableBound * ones(nMetabolites, 1); ...
        -options.variableBound * ones(nFixed, 1); ...
        -options.variableBound * ones(nVariable, 1); ...
        zeros(nVariable, 1); ...
        zeros(nVariable, 1)];
    ub = [ ...
        options.variableBound * ones(nMetabolites, 1); ...
        options.variableBound * ones(nFixed, 1); ...
        zeros(nVariable, 1); ...
        options.variableBound * ones(nVariable, 1); ...
        ones(nVariable, 1)];

    objective = zeros(totalVariables, 1);
    objective(uStart:(uStart + nVariable - 1)) = 1;

    % Preserve the historical integer formulation.
    integerVariables = 1:totalVariables;
    solverOptions = optimoptions( ...
        "intlinprog", ...
        "Display", options.display, ...
        "IntegerPreprocess", "none", ...
        "LPPreprocess", "none");

    supports = zeros(0, nVariable);
    coefficients = cell(0, 1);
    solverRuns = 0;
    searchedDegrees = zeros(0, 1);

    for degree = options.minDegree:options.maxDegree
        searchedDegrees(end + 1, 1) = degree;

        while size(supports, 1) < options.maxSolutions
            degreeRow = sparse(1, totalVariables);
            degreeRow(sStart:(sStart + nVariable - 1)) = 1;
            Aeq = [AeqBase; degreeRow];
            beq = [beqBase; degree];

            [solution, ~, exitFlag, solverOutput] = intlinprog( ...
                objective, ...
                integerVariables, ...
                Ain, ...
                bin, ...
                Aeq, ...
                beq, ...
                lb, ...
                ub, ...
                solverOptions);
            solverRuns = solverRuns + 1;

            if exitFlag == -2
                break;
            end
            if exitFlag <= 0 || isempty(solution)
                error( ...
                    "fluto:TradeoffSolverFailed", ...
                    "intlinprog failed at degree %d (exitFlag=%d): %s", ...
                    degree, ...
                    exitFlag, ...
                    localSolverMessage(solverOutput));
            end

            solution = round(solution);
            alphaVariable = solution( ...
                variableAlphaStart:(variableAlphaStart + nVariable - 1));
            selected = find(alphaVariable < 0);

            if numel(selected) ~= degree
                error( ...
                    "fluto:UnexpectedMILPSupport", ...
                    ["MILP selected %d negative variable coefficients at degree %d. " ...
                     "Check solver tolerances and formulation assumptions."], ...
                    numel(selected), ...
                    degree);
            end

            alternatives = fluto.expandCoupledAlternatives( ...
                selected, ...
                variableCoupling);
            if isempty(alternatives)
                alternatives = selected(:)';
            end

            for rowIndex = 1:size(alternatives, 1)
                supportMask = false(1, nVariable);
                supportMask(alternatives(rowIndex, :)) = true;

                if ~any(all(supports == supportMask, 2))
                    supports(end + 1, :) = supportMask;
                    coefficients{end + 1, 1} = alphaVariable(:);

                    exclusion = sparse(1, totalVariables);
                    chosen = find(supportMask);
                    exclusion(sStart + chosen - 1) = 1;
                    Ain(end + 1, :) = exclusion;
                    bin(end + 1, 1) = numel(chosen) - 1;
                end
            end

            % Always exclude the exact MILP support even if expansion mapped
            % to a previously known equivalent support.
            exactExclusion = sparse(1, totalVariables);
            exactExclusion(sStart + selected - 1) = 1;
            Ain(end + 1, :) = exactExclusion;
            bin(end + 1, 1) = numel(selected) - 1;
        end

        if size(supports, 1) >= options.maxSolutions
            break;
        end
    end

    tradeoffs = localTradeoffTable( ...
        model, ...
        variableIndices, ...
        supports, ...
        coefficients);

    result = struct( ...
        "tradeoffs", tradeoffs, ...
        "variableReactionIndices", variableIndices, ...
        "variableReactionNumbers", model.rxnNumber(variableIndices), ...
        "variableReactionIDs", string(model.rxns(variableIndices)), ...
        "supportMatrix", logical(supports), ...
        "solverRuns", solverRuns, ...
        "searchedDegrees", searchedDegrees, ...
        "options", options);
end

function options = localOptions(options)
    defaults = struct( ...
        "minDegree", 2, ...
        "maxDegree", 9, ...
        "bigM", 1001, ...
        "variableBound", 100, ...
        "display", "off", ...
        "maxSolutions", Inf);

    names = fieldnames(defaults);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(options, name)
            options.(name) = defaults.(name);
        end
    end

    validateattributes(options.minDegree, {'numeric'}, {'scalar', 'integer', 'positive'});
    validateattributes(options.maxDegree, {'numeric'}, {'scalar', 'integer', '>=', options.minDegree});
    validateattributes(options.bigM, {'numeric'}, {'scalar', 'real', 'positive', 'finite'});
    validateattributes(options.variableBound, {'numeric'}, {'scalar', 'real', 'positive', 'finite'});
    validateattributes(options.maxSolutions, {'numeric'}, {'scalar', 'real', 'positive'});
    options.display = char(string(options.display));
end

function message = localSolverMessage(output)
    message = "No solver message available.";
    if isstruct(output) && isfield(output, "message") && ~isempty(output.message)
        message = string(output.message);
    end
end

function tableOut = localTradeoffTable(model, variableIndices, supports, coefficients)
    nTradeoffs = size(supports, 1);
    tradeoffID = (1:nTradeoffs)';
    degree = zeros(nTradeoffs, 1);
    reactionNumbers = strings(nTradeoffs, 1);
    reactionIDs = strings(nTradeoffs, 1);
    coefficientText = strings(nTradeoffs, 1);

    for i = 1:nTradeoffs
        localIndices = find(supports(i, :));
        modelIndices = variableIndices(localIndices);
        degree(i) = numel(localIndices);
        reactionNumbers(i) = strjoin(string(model.rxnNumber(modelIndices)), ";");
        reactionIDs(i) = strjoin(string(model.rxns(modelIndices)), ";");

        alpha = coefficients{i};
        coefficientText(i) = strjoin(string(alpha(localIndices)'), ";");
    end

    tableOut = table( ...
        tradeoffID, ...
        degree, ...
        reactionNumbers, ...
        reactionIDs, ...
        coefficientText, ...
        "VariableNames", ...
        {"TradeoffID", "Degree", "ReactionNumbers", "ReactionIDs", "VariableCoefficients"});
end

function result = localEmptyResult(model, variableIndices, options)
    result = struct( ...
        "tradeoffs", table( ...
            zeros(0, 1), ...
            zeros(0, 1), ...
            strings(0, 1), ...
            strings(0, 1), ...
            strings(0, 1), ...
            "VariableNames", ...
            {"TradeoffID", "Degree", "ReactionNumbers", "ReactionIDs", "VariableCoefficients"}), ...
        "variableReactionIndices", variableIndices, ...
        "variableReactionNumbers", model.rxnNumber(variableIndices), ...
        "variableReactionIDs", string(model.rxns(variableIndices)), ...
        "supportMatrix", false(0, numel(variableIndices)), ...
        "solverRuns", 0, ...
        "searchedDegrees", zeros(0, 1), ...
        "options", options);
end
