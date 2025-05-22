%written by Kyle Adams

%this file creates the supplementary sensitivity index tables
%you need to have all csv files which should have the total indices and
%first order indices separate. this is automatically saved when SobolMain
%is run.


%creating an array of csv files for total sensitivity indices
csv_files_ST = { %rename to actual file names in environment
    'sensitivity_resultsST_20250514_124510.csv', ...
    'sensitivity_resultsST_20250514_134514.csv', ...
    'sensitivity_resultsST_20250514_151721.csv', ...
    'sensitivity_resultsST_20250514_140844.csv',
};

%creating an array of csv files for first order sensitivity indices
csv_files_S1 = {
    'sensitivity_resultsS1_20250514_124510.csv', ...
    'sensitivity_resultsS1_20250514_134514.csv', ...
    'sensitivity_resultsS1_20250514_151721.csv', ...
    'sensitivity_resultsS1_20250514_140844.csv',

};

%naming columns in descending number of base samples
column_names = {'ST_175k', 'ST_150k', 'ST_125k', 'ST_100k'};
column_names_s1 = {'S1_175k', 'S1_150k', 'S1_125k', 'ST_100k'};

%loading paramnames, setting num_runs if you have different number of runs
p = parameters();
param_names = fieldnames(p);
num_params = length(param_names);
num_runs = length(csv_files_ST);

%populating ST table
ST_matrix = zeros(num_params, num_runs);
for j = 1:num_runs
    T = readtable(csv_files_ST{j});  % columns: parameter, ST
    [~, idx] = ismember(param_names, T.Parameter);
    ST_matrix(:, j) = T.ST(idx);
end

%populating S1 table
S1_matrix = zeros(num_params, num_runs);
for j = 1:num_runs
    T = readtable(csv_files_S1{j});  % columns: parameter, S1
    [~, idx] = ismember(param_names, T.Parameter);
    S1_matrix(:, j) = T.S1(idx);
end

%sorting by descending ST values from run with largest base samples (loaded as
%column 1 for us)
[~, sort_idx_ST] = sort(ST_matrix(:, 1), 'descend');
ST_sorted = ST_matrix(sort_idx_ST, :);
S1_sorted = S1_matrix(sort_idx_ST, :);
param_names_sorted_ST = param_names(sort_idx_ST);

%making tables with sorted sensitivity indices and matching parameter names
T_ST = array2table(ST_sorted, ...
    'VariableNames', column_names, ...
    'RowNames', param_names_sorted_ST);

T_S1 = array2table(S1_sorted, ...
    'VariableNames', column_names_s1, ...
    'RowNames', param_names_sorted_ST);

%saving tables to environment
writetable(T_ST, 'summary_sensitivity_ST_table.csv', 'WriteRowNames', true);
writetable(T_S1, 'summary_sensitivity_S1_table.csv', 'WriteRowNames', true);

