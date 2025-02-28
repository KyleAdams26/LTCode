% Plots for LT

% Necessary files: Sobol indices from sensitivity analysis;
% QOI distribution from varying all parameters, varying most influential 7 parameters, and
% varying least influential 28 parameters


% Reading in data from the MATLAB SA results (equivalent of pd.read_csv in Python)
SA_values = readtable('3milSAValues.txt', 'Delimiter', ',', 'ReadVariableNames', false);
SA_values.Properties.VariableNames = {'S1', 'ST'};

% Adding a column of parameter titles
parameters = {'lL', 'dA', 'sR', 'dR', 'aIR', 'bIR', 'aCI', 'bCI', 'aHI', 'bHI', 'lC', 'gC', 'KC', ...
    'aIC', 'bIC', 'lH', 'gH', 'KH', 'aIH', 'bIH', 'lR', 'dI', 'aHC', 'bHC', 'dC', 'aAH', 'bAH', 'aRA', 'bRA', ...
    'aIRA', 'bIRA', 'dH', 'dL', 'aCL', 'bCL'};

% Check if the number of rows in SA_values matches the number of parameters
% (used for debugging)
num_rows = height(SA_values);
if num_rows ~= length(parameters)
    error('The number of rows in SA_values (%d) does not match the number of parameters (%d).', num_rows, length(parameters));
end

% Add the parameter list as a new column
SA_values.param = parameters(:); % Ensure it's a column vector

% Sort by 'ST' in descending order and select the top 10 influential parameters
SA_values = sortrows(SA_values, 'ST', 'descend');
SA_values = SA_values(1:10, :);

% Extract ST and S1 for plotting
ST = SA_values.ST;
S1 = SA_values.S1;

% Set up plotting
x = SA_values.param;
x_pos = 1:length(x); % For bar placement

% Bar width for grouped bars
bar_width = 0.35;

% Create the figure
figure;

% Plot ST and S1 as bar charts
hold on;
rectsST = bar(x_pos - bar_width / 2, ST, bar_width, 'FaceColor', [0 0.188 0.69]);
rectsS1 = bar(x_pos + bar_width / 2, S1, bar_width, 'FaceColor', [1 0.79 0.63]);

% Add labels and title
set(gca, 'XTick', x_pos, 'XTickLabel', x, 'XTickLabelRotation', 45);
xlabel('Parameters', 'FontSize', 12, 'FontName', 'serif');
ylabel('Values', 'FontSize', 12, 'FontName', 'serif');
title('Sobol Indices', 'FontSize', 12, 'FontName', 'serif');
legend({'ST', 'S1'}, 'FontSize', 12, 'FontName', 'serif');
hold off;

% Show the plot
grid on;

% Load and plot PDFs of QOIs with all parameters varying vs 7 varying (and
% then all vs 28)

% Loading data from text files
ALLCELLS = readmatrix('3milALLCELLS.txt', 'Delimiter', ',');


SEVENCELLS = readmatrix('3milSEVENCELLS.txt'); 
%SEVENCELLS = SEVENCELLS(~isnan(SEVENCELLS));  % Remove NaN values
TWENTYEIGHTCELLS = readmatrix('3milTWENTYEIGHTCELLS.txt');

% Use the subset function to get a random subset of ALLCELLS with the same size as SEVENCELLS
subset_allcells_seven = subset(ALLCELLS, length(SEVENCELLS));

% Define common bin edges for the histograms
%num_bins = 100; % Specify the number of bins
%bin_edges = linspace(0, 2e11, num_bins+1); % Create bin edges from 0 to 2e11

% Plot histograms of the random subset against SEVENCELLS
figure;
hold on;
histogram(subset_allcells_seven, 'EdgeColor', 'red', 'DisplayStyle', 'stairs'); % Random subset of ALLCELLS
histogram(SEVENCELLS,'EdgeColor', 'black', 'DisplayStyle', 'stairs'); %seven cells

% Create custom legend handles using 'line'
all_params_handle = plot(nan, nan, 'Color', 'red', 'LineWidth', 1.5);
seven_params_handle = plot(nan, nan, 'Color', 'black', 'LineWidth', 1.5);

% Add legend
legend([all_params_handle, seven_params_handle], ...
    {'All Parameters Varied', '7 Parameters Varied'}, ...
    'FontSize', 12, 'FontName', 'serif', 'Location', 'northwest');
xlabel('QOI Value', 'FontSize', 12, 'FontName', 'serif');
ylabel('Frequency', 'FontSize', 12, 'FontName', 'serif');
title('Histogram Comparison', 'FontSize', 12, 'FontName', 'serif');
xlim([0 2e11]); % Set x-axis limits from 0 to 2x10^11
hold off;

% Use the subset function to get a random subset of ALLCELLS with the same size as TWENTYEIGHTCELLS
subset_allcells_twenty_eight = subset(ALLCELLS, length(TWENTYEIGHTCELLS));

% Plot histograms of the random subset against 3milTWENTYEIGHTCELLS
figure;
hold on;
histogram(subset_allcells_twenty_eight, 'EdgeColor', 'red', 'DisplayStyle', 'stairs'); % Random subset of ALLCELLS
histogram(TWENTYEIGHTCELLS, 100, 'EdgeColor', 'blue', 'DisplayStyle', 'stairs');

% Create custom legend handles using 'line'
all_params_handle = plot(nan, nan, 'Color', 'red', 'LineWidth', 1.5);
twentyeight_params_handle = plot(nan, nan, 'Color', 'blue', 'LineWidth', 1.5);

% Add legend
legend([all_params_handle, twentyeight_params_handle], ...
    {'All Parameters Varied', '28 Parameters Varied'}, ...
    'FontSize', 12, 'FontName', 'serif', 'Location', 'northwest');
xlabel('QOI Value', 'FontSize', 12, 'FontName', 'serif');
ylabel('Frequency', 'FontSize', 12, 'FontName', 'serif');
title('Histogram Comparison', 'FontSize', 12, 'FontName', 'serif');
xlim([0 2e11]); % Set x-axis limits from 0 to 2x10^11
hold off;

%this is me debugging
%fprintf('Length of ALLCELLS: %d\n', length(ALLCELLS));
%fprintf('Length of SEVENCELLS: %d\n', length(SEVENCELLS));
%fprintf('Length of subset_allcells_seven: %d\n', length(subset_allcells_seven));

%fprintf('Length of TWENTYEIGHTCELLS: %d\n', length(TWENTYEIGHTCELLS));
%fprintf('Length of subset_allcells_twenty_eight: %d\n', length(subset_allcells_twenty_eight));