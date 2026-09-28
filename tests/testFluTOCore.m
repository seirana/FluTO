function tests = testFluTOCore
%TESTFLUTOCORE Unit tests for solver-independent maintained FluTO utilities.

    tests = functiontests(localfunctions);
end

function testValidateModelAddsStableIndices(testCase)
    model = localToyModel();

    validated = fluto.validateModel(model);

    verifyEqual(testCase, validated.rxnNumber, (1:4)');
    verifyEqual(testCase, validated.metNumber, (1:2)');
    verifySize(testCase, validated.rxnType, [4, 1]);
end

function testValidateModelRejectsInvalidBounds(testCase)
    model = localToyModel();
    model.lb(2) = 2;
    model.ub(2) = 1;

    verifyError( ...
        testCase, ...
        @() fluto.validateModel(model), ...
        "fluto:InvalidBounds");
end

function testCanonicalizeNegativeIrreversibleReaction(testCase)
    model = localToyModel();
    originalColumn = model.S(:, 2);

    [canonical, flipped] = ...
        fluto.canonicalizeReactionDirections(model, 1e-9);

    verifyTrue(testCase, flipped(2));
    verifyEqual(testCase, canonical.lb(2), 0);
    verifyEqual(testCase, canonical.ub(2), 10);
    verifyEqual(testCase, canonical.S(:, 2), -originalColumn);
end

function testClassifyFluxesUsesTolerance(testCase)
    model = localToyModel();
    [model, ~] = fluto.canonicalizeReactionDirections(model, 1e-9);

    model.lb = [1; 0; -2; 0];
    model.ub = [1 + 1e-10; 4; 3; 0];

    [classified, tableOut] = ...
        fluto.classifyFluxes(model, 1e-8);

    verifyEqual( ...
        testCase, ...
        string(classified.rxnType), ...
        ["fixed"; "variable"; "reversible"; "fixed"]);
    verifyEqual( ...
        testCase, ...
        tableOut.FluxClass, ...
        ["fixed"; "variable"; "reversible"; "fixed"]);
end

function testApplyConditionUsesReactionNumbers(testCase)
    model = fluto.validateModel(localToyModel());
    model.rxnNumber = [10; 20; 30; 40];

    condition = struct( ...
        "blockReactionNumbers", [10, 20], ...
        "activeReactionNumber", 20, ...
        "activeBounds", [-7, -7], ...
        "fixedReactionNumbers", [30, 40], ...
        "fixedFluxValues", [1.5, 2.5]);

    changed = fluto.applyCondition(model, condition);

    verifyEqual(testCase, changed.lb, [0; -7; 1.5; 2.5]);
    verifyEqual(testCase, changed.ub, [0; -7; 1.5; 2.5]);
end

function testApplyConditionRejectsUnknownReactionNumber(testCase)
    model = fluto.validateModel(localToyModel());

    condition = struct( ...
        "fixedReactionNumbers", 999, ...
        "fixedFluxValues", 1);

    verifyError( ...
        testCase, ...
        @() fluto.applyCondition(model, condition), ...
        "fluto:UnknownReactionNumber");
end

function testExpandCoupledAlternativesIncludesSelf(testCase)
    coupling = false(4);
    coupling(1, 3) = true;
    coupling(3, 1) = true;
    coupling(2, 4) = true;
    coupling(4, 2) = true;

    combinations = fluto.expandCoupledAlternatives( ...
        [1, 2], ...
        coupling);

    expected = [ ...
        1, 2; ...
        1, 4; ...
        2, 3; ...
        3, 4];

    verifyEqual(testCase, combinations, expected);
end

function testExpandCoupledAlternativesAvoidsDuplicateMembers(testCase)
    coupling = false(3);
    coupling(1, 2) = true;
    coupling(2, 1) = true;

    combinations = fluto.expandCoupledAlternatives( ...
        [1, 2], ...
        coupling);

    verifyEqual(testCase, combinations, [1, 2]);
end

function model = localToyModel()
    model = struct();
    model.S = [ ...
        1, -1, 0, 0; ...
        0, 1, -1, 0];
    model.rxns = {"R1"; "R2"; "R3"; "R4"};
    model.mets = {"M1"; "M2"};
    model.lb = [0; -10; -2; 0];
    model.ub = [10; 0; 3; 0];
    model.c = zeros(4, 1);
    model.subSystems = {"A"; "A"; "B"; "B"};
end
