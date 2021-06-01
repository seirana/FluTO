function Wtite_down_results(whole_tradeoffs, lst_variable, adrs_mdl_folder, file_name, cb)
%%% write down the results
title = strings(1,size(whole_tradeoffs,2)+4);
title(1,1) = 'Tradeoff-degree';
blnk = {'#' '#' '#' '#'};
solution_info = [blnk; string(lst_variable)];
solution_ = [solution_info whole_tradeoffs];
solution_ = [title; solution_];

for i = 5:size(solution_,2)
    si = find(str2double(string(solution_(3:size(lst_variable,1)+2,i))) ~= 0);
    solution_(1,i) = {size(si,1)};
end

nonzero_rows = {};
for i = 2:size(solution_,1)
    a = find(str2double(string(solution_(i,5:end)))~= 0);
    if size(a) > 0
        nonzero_rows(end+1,:) = {solution_(i,:)};
    end
end

T = cell2table(nonzero_rows);
writetable(T, strcat(adrs_mdl_folder, file_name, '_tradeoffs_', string(cb), '.xlsx'));
clear T;
clear solution_;
clear title;
end