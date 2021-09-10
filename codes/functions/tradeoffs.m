function tradeoffs(mdl_, fctable, cb, adrs_mdl_folder, file_name)
%%% 'reversible' reaction is called 'sign variable' reaction in the paper
%%% 'variable is called 'fixed-sign variable' reaction in the paper
for i = size(mdl_.rxnType,1):-1:1
    if strcmp(string(mdl_.rxnType(i)), 'fixed') == 1 || strcmp(string(mdl_.rxnType(i)), 'reversible') == 1  
        fctable(i,:) = [];
        fctable(:,i) = [];
    end
end

%%% N_fixed
N_fixed = [];
j = 1;
lst_fixed = {};
for i = 1:size(mdl_.rxns,1)
    if strcmp(mdl_.rxnType(i),{'fixed'})== 1
        N_fixed(:,j) = mdl_.S(:,i);
        lst_fixed(j,1)  = {'fixed'};
        lst_fixed(j,2) = mdl_.rxns(i);
        lst_fixed(j,3) = {mdl_.rxnNumber(i)};
        if size(mdl_.subSystems{i,1},2) <= 1
            lst_fixed(j,4) = mdl_.subSystems{i,1};
        else
            if size(mdl_.subSystems{i,1},2) > 1
                lst_fixed(j,4) = {mdl_.subSystems{i,1}{1,1}};
            end
        end
        j = j+1;
    end
end

%%% N_variable
N_variable = [];
j = 1;
lst_variable = {};
for i = 1:size(mdl_.rxns,1)
    if strcmp(mdl_.rxnType(i),{'variable'}) == 1
        N_variable(:,j) = mdl_.S(:,i);
        lst_variable(j,1) = {'variable'};
        lst_variable(j,2) = mdl_.rxns(i);
        lst_variable(j,3) = {mdl_.rxnNumber(i)};
        if size(mdl_.subSystems{i,1},2) <= 1
            lst_variable(j,4) = mdl_.subSystems{i,1};
        else
            if size(mdl_.subSystems{i,1},2) > 1
                lst_variable(j,4) = {mdl_.subSystems{i,1}{1,1}};
            end
        end
        j = j+1;
    end
end

if size(fctable,1) == size(N_fixed,1)
    disp('checked');
    disp('');
end

%%% N_reversible
N_reversible = [];
j = 1;
for i = 1:size(mdl_.rxns,1)
    if strcmp(mdl_.rxnType(i),{'reversible'}) == 1
        N_reversible(:,j) = mdl_.S(:,i);
        j = j+1;
    end
end

lp_tradeoffs = [];
whole_tradeoffs = [];
tr_fc_lst = zeros(size(lst_variable,1)+2,1);
rep = 0;
M = 1001;

N = [N_fixed, N_variable, N_reversible];

%%% Aeq = [k, alpha, u, s]
r = size(N,2);
m = size(N,1);
v2 = size(N_variable,2);
f2 = size(N_fixed,2);

alpha = zeros(r,f2+v2);
for xi = 1:f2+v2
    alpha(xi,xi) = -1;
end
Aeq1 = [N', alpha, zeros(r,v2), zeros(r,v2)]; % k*N = alpha
beq1 = zeros(r,1);

tr_size = 2;
Aeq2 = [zeros(1,m), zeros(1,f2+v2), zeros(1,v2), ones(1,v2)]; % tr_size <= sum s %%% = trade-off length
beq2 = tr_size;

Aeq = [Aeq1;Aeq2];
beq = [beq1;beq2];

Ain1 = [zeros(1,m), -1*ones(1,f2), zeros(1,v2), zeros(1,v2), zeros(1,v2)]; % 1 <= sum alpha_fixed
bin1 = -1;

Ain2 = [zeros(v2,m), zeros(v2,f2), eye(v2), -1*eye(v2), zeros(v2)]; % alpha_variable <= u
bin2 = zeros(v2,1);

Ain3 = [zeros(v2,m), zeros(v2,f2), -1*eye(v2), -1*eye(v2), zeros(v2)]; % -u <= alpha_variable
bin3 = zeros(v2,1);

Ain4 = [zeros(v2,m), zeros(v2,f2), eye(v2), zeros(v2), M*eye(v2)]; % alpha_variable + M*s <= M-1
bin4 = (M-1)*ones(v2,1);

Ain5 = [zeros(v2,m), zeros(v2,f2), -1*eye(v2), zeros(v2), -1*M*eye(v2)]; % 0 <= alpha_variable_i + M*s
bin5= zeros(v2,1);

Ain = [Ain1;Ain2;Ain3;Ain4;Ain5];
bin = [bin1;bin2;bin3;bin4;bin5];

%%% fully coupled reactions can not
for i = 1:v2-1
    for j = i+1:v2
        if fctable(i,j) == 1
            s_in = zeros(1,v2);
            s_in(1,i) = 1;
            s_in(1,j) = 1;
            new_bin = [zeros(1,m), zeros(1,f2+v2), zeros(1,v2), s_in];
            Ain = [Ain; new_bin];
            bin = [bin; 1];
        end
    end
end

lb = [-100*ones(m,1); -100*ones(f2,1); -100*ones(v2,1); zeros(v2,1); zeros(v2,1)];
ub = [100*ones(m,1); 100*ones(f2,1); zeros(v2,1); 100*ones(v2,1); ones(v2,1)];


intcon = 1:size(lb,1); % integer variables
fk = [zeros(m,1); zeros(f2+v2,1); ones(v2,1); zeros(v2,1)]; % factors for objective function
options = optimoptions('intlinprog', 'IntegerPreprocess','none','LPPreprocess','none');

while 1
    sol = intlinprog(fk,intcon,Ain,bin,Aeq,beq,lb,ub,options);
    sol = round(sol);
    %%% add recently found trade-off to the trade-off list;
    if ~isempty(sol)
        
        %%% update LP trade-off list
        lp_tradeoffs = [lp_tradeoffs sol(m+1:m+f2+v2)];
        
        %%% find the trade-off
        new_tradeoff_set = find(sol(m+f2+1:m+f2+v2) < 0);
        
        %%% if there is any fully coupled reactions, find them
        elements = {};
        for i = 1:size(new_tradeoff_set,1)
            fc = find(fctable(new_tradeoff_set(i),:) == 1);
            elements = [elements {fc}];
        end
        
        %%% find all combination of tradeoff among fully coupled rxns
        combinations = cell(1, numel(elements));
        [combinations{:}] = ndgrid(elements{:});
        combinations = cellfun(@(x) x(:), combinations,'uniformoutput',false);
        result = [combinations{:}]; % list of new trade-offs
        disp('result');
        disp(result);
        disp('-------');
        
        %%% add new trade-offs
        rep = rep+1;
        for i = 1:size(result,1)
            % Whole_tradeoffs
            sol = zeros(1+v2,1);
            sol(1) = rep;
            for j = 1:size(result,2)
                sol(result(i,j)+1,1) = 1;
            end
            
            %%% update Whole_tradeoff list
            whole_tradeoffs = [whole_tradeoffs sol];
            
            %%% add new conditions to the program
            s_in = zeros(1,v2);
            for j = 1:size(result,2)
                s_in(1,result(i,j)) = 1;
            end
            
            new_bin = [zeros(1,m), zeros(1,f2+v2), zeros(1,v2), s_in];
            Ain = [Ain; new_bin];
            bin = [bin; size(result,2)-1];
        end
        
        tr_fc_lst(1,rep) = rep;
        tr_fc_lst(2,rep) = tr_size;
        for j = 1:size(result,2)
            % update tr_fc_lst
            for i = 1:size(result,1)
                tr_fc_lst(result(i,j)+2,rep) = j;
            end
        end
        
        [Ain, bin, tr_fc_lst, fctable] = add_new_tradeoffs(fctable, rep, Ain, bin, tr_fc_lst, m, f2, v2);
    else
        if tr_size < 9
            tr_size = tr_size+1;
            beq(end,1) = tr_size;
        else
            break;
        end
    end
    %%% write down the results
    Wtite_down_results(whole_tradeoffs, lst_variable, adrs_mdl_folder, file_name, cb);
end
end
