function output = runFluTO(modelFile, conditionFile, outputDirectory)
%RUNFLUTO Run the maintained FluTO pipeline from files.
%
%   OUTPUT = runFluTO(MODELFILE, CONDITIONFILE, OUTPUTDIRECTORY)
%   loads a COBRA-compatible model, reads a JSON condition specification,
%   runs FluTO, and writes CSV/JSON artifacts.
%
%   Example:
%       addpath("src")
%       addpath("scripts")
%       output = runFluTO( ...
%           "models/Ecoli_iJO1366.mat", ...
%           "examples/condition.template.json", ...
%           "artifacts/example");

    arguments
        modelFile (1, 1) string
        conditionFile (1, 1) string
        outputDirectory (1, 1) string = "artifacts"
    end

    repositoryRoot = fileparts(fileparts(mfilename("fullpath")));
    addpath(fullfile(repositoryRoot, "src"));

    if exist("readCbModel", "file") ~= 2
        error( ...
            "fluto:CobraDependencyMissing", ...
            ["COBRA Toolbox function readCbModel was not found. " ...
             "Initialize COBRA Toolbox before running FluTO."]);
    end

    if ~isfile(modelFile)
        error("fluto:ModelFileNotFound", "Model file not found: %s", modelFile);
    end
    if ~isfile(conditionFile)
        error( ...
            "fluto:ConditionFileNotFound", ...
            "Condition file not found: %s", ...
            conditionFile);
    end

    condition = jsondecode(fileread(conditionFile));
    model = readCbModel(modelFile);

    output = fluto.runCondition(model, condition, struct());

    metadata = struct( ...
        "modelFile", char(modelFile), ...
        "conditionFile", char(conditionFile), ...
        "matlabRelease", version("-release"), ...
        "matlabVersion", version, ...
        "inputSHA256", struct( ...
            "model", localSha256(modelFile), ...
            "condition", localSha256(conditionFile)));

    if isfield(condition, "name")
        conditionName = string(condition.name);
    else
        [~, conditionName] = fileparts(conditionFile);
        conditionName = string(conditionName);
    end

    output.paths = fluto.writeResults( ...
        outputDirectory, ...
        conditionName, ...
        output.fluxClassification, ...
        output.tradeoffResult, ...
        metadata);
end

function digest = localSha256(path)
    bytes = filereadBytes(path);
    engine = java.security.MessageDigest.getInstance("SHA-256");
    engine.update(bytes);
    digest = lower(string(reshape(dec2hex(typecast(engine.digest(), "uint8"), 2).', 1, [])));
end

function bytes = filereadBytes(path)
    fid = fopen(path, "rb");
    if fid < 0
        error("fluto:FileReadFailed", "Could not open %s.", path);
    end
    cleanup = onCleanup(@() fclose(fid));
    bytes = fread(fid, Inf, "*uint8");
    clear cleanup;
end
