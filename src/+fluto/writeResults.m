function paths = writeResults(outputDirectory, conditionName, fluxTable, result, metadata)
%WRITERESULTS Write FluTO results as machine-readable files.

    outputDirectory = string(outputDirectory);
    conditionName = string(conditionName);
    if strlength(conditionName) == 0
        error("fluto:InvalidConditionName", "conditionName must not be empty.");
    end

    if ~isfolder(outputDirectory)
        mkdir(outputDirectory);
    end

    stem = regexprep(conditionName, "[^A-Za-z0-9._-]+", "_");

    fluxPath = fullfile(outputDirectory, stem + "_flux_ranges.csv");
    tradeoffPath = fullfile(outputDirectory, stem + "_tradeoffs.csv");
    metadataPath = fullfile(outputDirectory, stem + "_metadata.json");

    writetable(fluxTable, fluxPath);
    writetable(result.tradeoffs, tradeoffPath);

    metadata.solverRuns = result.solverRuns;
    metadata.searchedDegrees = result.searchedDegrees;
    metadata.nTradeoffs = height(result.tradeoffs);
    metadata.tradeoffOptions = result.options;

    fid = fopen(metadataPath, "w");
    if fid < 0
        error("fluto:OutputWriteFailed", "Could not open %s for writing.", metadataPath);
    end
    cleanup = onCleanup(@() fclose(fid));
    fprintf(fid, "%s", jsonencode(metadata, "PrettyPrint", true));
    clear cleanup;

    paths = struct( ...
        "fluxRanges", string(fluxPath), ...
        "tradeoffs", string(tradeoffPath), ...
        "metadata", string(metadataPath));
end
