function [couplingMatrix, diagnostics] = buildFullyCoupledMatrix(model, options)
%BUILDFULLYCOUPLEDMATRIX Identify pairwise fully coupled variable reactions.
%
%   The historical FluTO implementation contained a copy/paste defect in the
%   second coupling check: it solved the same constrained optimization twice
%   and then compared the duplicate results. This implementation performs the
%   intended symmetric check:
%
%     1. fix reaction j at an interior feasible value and test whether i is fixed;
%     2. fix reaction i at an interior feasible value and test whether j is fixed.
%
%   Connected fully-coupled components are completed transitively.

    if nargin < 2 || isempty(options)
        options = struct();
    end
    options = localOptions(options);

    if exist("changeObjective", "file") ~= 2 || exist("optimizeCbModel", "file") ~= 2
        error( ...
            "fluto:CobraDependencyMissing", ...
            "COBRA Toolbox functions changeObjective and optimizeCbModel are required.");
    end

    model = fluto.validateModel(model);
    nReactions = numel(model.rxns);
    couplingMatrix = false(nReactions);

    variableMask = ...
        abs(model.ub - model.lb) > options.tolerance ...
        & model.lb >= -options.tolerance ...
        & model.ub > options.tolerance;
    variableIndices = find(variableMask);

    testedPairs = 0;
    coupledPairs = 0;

    for leftPosition = 1:numel(variableIndices) - 1
        i = variableIndices(leftPosition);

        for rightPosition = leftPosition + 1:numel(variableIndices)
            j = variableIndices(rightPosition);

            if ~localCompatibleSignClass(model, i, j, options.tolerance)
                continue;
            end

            testedPairs = testedPairs + 1;

            referenceJ = localInteriorFlux(model.lb(j), model.ub(j), options.tolerance);
            [minimumI, maximumI] = localRangeWithFixedReaction( ...
                model, ...
                i, ...
                j, ...
                referenceJ);

            if ~localApproximatelyEqual(minimumI, maximumI, options.tolerance) ...
                    || abs(0.5 * (minimumI + maximumI)) <= options.tolerance
                continue;
            end

            referenceI = localInteriorFlux(model.lb(i), model.ub(i), options.tolerance);
            [minimumJ, maximumJ] = localRangeWithFixedReaction( ...
                model, ...
                j, ...
                i, ...
                referenceI);

            if localApproximatelyEqual(minimumJ, maximumJ, options.tolerance) ...
                    && abs(0.5 * (minimumJ + maximumJ)) > options.tolerance
                couplingMatrix(i, j) = true;
                couplingMatrix(j, i) = true;
                coupledPairs = coupledPairs + 1;
            end
        end
    end

    couplingMatrix = localTransitiveClosure(couplingMatrix);
    couplingMatrix(1:nReactions + 1:end) = false;

    diagnostics = struct( ...
        "nReactions", nReactions, ...
        "nEligibleVariableReactions", numel(variableIndices), ...
        "nPairsTested", testedPairs, ...
        "nDirectPairsCoupled", coupledPairs, ...
        "nPairsAfterTransitiveClosure", nnz(triu(couplingMatrix, 1)), ...
        "tolerance", options.tolerance);
end

function options = localOptions(options)
    if ~isfield(options, "tolerance")
        options.tolerance = 1e-7;
    end
    validateattributes(options.tolerance, {"numeric"}, {"scalar", "real", "positive", "finite"});
end

function tf = localCompatibleSignClass(model, i, j, tolerance)
    bothTouchZero = ...
        abs(model.lb(i)) <= tolerance ...
        && abs(model.lb(j)) <= tolerance;
    bothStrictPositive = ...
        model.lb(i) > tolerance ...
        && model.lb(j) > tolerance;
    tf = bothTouchZero || bothStrictPositive;
end

function value = localInteriorFlux(lower, upper, tolerance)
    value = 0.5 * (lower + upper);
    if abs(value) <= tolerance
        if upper > tolerance
            value = 0.5 * upper;
        elseif lower < -tolerance
            value = 0.5 * lower;
        else
            error("fluto:NoInteriorFlux", "Reaction has no non-zero interior flux.");
        end
    end
end

function [minimumValue, maximumValue] = localRangeWithFixedReaction( ...
        model, objectiveIndex, fixedIndex, fixedValue)

    constrained = model;
    constrained.lb(fixedIndex) = fixedValue;
    constrained.ub(fixedIndex) = fixedValue;
    constrained = changeObjective(constrained, constrained.rxns{objectiveIndex}, 1);

    minimumSolution = optimizeCbModel(constrained, "min");
    maximumSolution = optimizeCbModel(constrained, "max");

    minimumValue = localObjectiveValue(minimumSolution);
    maximumValue = localObjectiveValue(maximumSolution);
end

function value = localObjectiveValue(solution)
    if ~isstruct(solution) || ~isfield(solution, "f") || isempty(solution.f)
        error("fluto:OptimizationFailed", "COBRA optimization returned no objective value.");
    end
    if isfield(solution, "stat") && ~isempty(solution.stat) && solution.stat ~= 1
        error( ...
            "fluto:OptimizationFailed", ...
            "COBRA optimization returned status %g.", ...
            double(solution.stat));
    end
    value = double(solution.f);
    if ~isscalar(value) || ~isfinite(value)
        error("fluto:OptimizationFailed", "COBRA optimization returned a non-finite objective.");
    end
end

function tf = localApproximatelyEqual(first, second, tolerance)
    scale = max([1, abs(first), abs(second)]);
    tf = abs(first - second) <= tolerance * scale;
end

function closure = localTransitiveClosure(adjacency)
    closure = logical(adjacency);
    n = size(closure, 1);

    for k = 1:n
        closure = closure | (closure(:, k) & closure(k, :));
    end

    closure = closure | closure';
end
