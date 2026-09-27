function combinations = expandCoupledAlternatives(selectedIndices, couplingMatrix)
%EXPANDCOUPLEDALTERNATIVES Expand a support through fully-coupled equivalents.
%
%   Each selected reaction may be replaced by itself or by a reaction in its
%   fully-coupled component. Returned rows are sorted unique reaction-index
%   sets and never collapse two support positions onto the same reaction.

    selected = double(selectedIndices(:)');
    n = size(couplingMatrix, 1);

    if size(couplingMatrix, 2) ~= n
        error("fluto:InvalidCouplingMatrix", "couplingMatrix must be square.");
    end
    if any(selected < 1) || any(selected > n) || any(mod(selected, 1) ~= 0)
        error("fluto:InvalidReactionIndex", "selectedIndices contains invalid indices.");
    end

    if isempty(selected)
        combinations = zeros(0, 0);
        return;
    end

    choices = cell(1, numel(selected));
    for i = 1:numel(selected)
        equivalents = find(couplingMatrix(selected(i), :));
        choices{i} = unique([selected(i), equivalents], "stable");
    end

    grids = cell(1, numel(choices));
    [grids{:}] = ndgrid(choices{:});
    combinations = zeros(numel(grids{1}), numel(choices));

    for i = 1:numel(grids)
        combinations(:, i) = grids{i}(:);
    end

    valid = false(size(combinations, 1), 1);
    for row = 1:size(combinations, 1)
        combinations(row, :) = sort(combinations(row, :));
        valid(row) = numel(unique(combinations(row, :))) == numel(selected);
    end

    combinations = unique(combinations(valid, :), "rows", "stable");
end
