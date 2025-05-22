% sobol file to run the Sobol sensitivity analysis
% on the QSP ODE model. 
% This code originally stems from Jaimit Parikh's SobolMain.m script.
% Kyle Adams added comments and adjusted this code to compare QOI
% distributions resulting from fixing some parameters.

random_seed = 1;
rng(random_seed);  %Setting seed 

% ###Step 1
%setting up bounds and distributions for sobol sampling
lower_percentage = 0.5;
upper_percentage = 1.5;
base_samples = 175000;
param_dist = {'Uniform'};
p = parameters();

% ###Step 2
%defining parameter sets to compare QOI distributions
param_sets = {
    'All', fieldnames(p)', ...
    'Top 6 Influential Parameters Varied', {'aCL', 'dL', 'bCL', 'KC', 'gC', 'aIC'}, ...
    'Least 29 Influential Parameters Varied', {'aHI', 'dC', 'dI', 'aIH', 'gH', 'dH', 'bHI', 'bIH', 'bIC', ...
                  'KH', 'aHC', 'aCI', 'bHC', 'bCI', 'aIC', 'dA', 'lL', 'sR', ...
                  'bIR', 'aRA', 'aIRA', 'lR', 'aIR', 'bRA', 'bIRA', 'bAH', ...
                  'dR', 'aAH', 'lH'} ...
    };

num_param_sets = length(param_sets) / 2;
QOIs = cell(1, num_param_sets);
labels = cell(1, num_param_sets);

%performing a sensitivity analysis for each parameter set
for i = 1:num_param_sets
    set_name = param_sets{2*i - 1};
    varying_params = param_sets{2*i};
    labels{i} = set_name;

    %using the name of the parameters from param_sets, extract all
    %parameter info from parameters.m
    p_subset = struct();
    for j = 1:length(varying_params)
        field = varying_params{j};
        if isfield(p, field)
            p_subset.(field) = p.(field);
        end
    end

    param_names = fieldnames(p_subset);
    num_params = length(param_names);
    lowBounds = zeros(1, num_params);
    upBounds = zeros(1, num_params);

    for j = 1:num_params
        param = param_names{j};
        lowBounds(j) = p_subset.(param) * lower_percentage;
        upBounds(j) = p_subset.(param) * upper_percentage;
    end

    parsObj.name = param_names';
    parsObj.lb = num2cell(repmat(-inf, 1, length(parsObj.name)));
    parsObj.ub = num2cell(inf(1, length(parsObj.name)));
    parsObj.dist = repmat(param_dist, 1, num_params);
    parsObj.N = base_samples;
    parsObj.parameters = arrayfun(@(i) {'lower', lowBounds(i), 'upper', upBounds(i)}, ...
                                  1:num_params, 'UniformOutput', false);

    %get sobol samples
    samples = getSamplesSobol(parsObj, false);

    %simulate the QOI
    QOI = zeros(1, length(samples));
    pN = cell(1,length(samples));
    parsName = parsObj.name;
    parfor ii = 1:length(samples)
        pN{ii} = updatePars(p, parsName, samples(ii, :));
        QOI(ii) = qoi(pN{ii});
    end
    QOIs{i} = QOI;
end
smallest_num_model_evals = min(cellfun(@length, QOIs)); %used to take subsets of bigger QOI distributions to compare fairly

% ###Step 3
%plotting histogram comparisons of QOI distributions
figure;
hold on;
colors = {'red', 'black', 'blue'};
for i = 1:num_param_sets
    histogram(subset(QOIs{i}, smallest_num_model_evals, random_seed), ...
        'EdgeColor', colors{i}, ...
        'DisplayStyle', 'stairs', ...
        'LineWidth', 2);
end
legend(labels, 'FontSize', 12, 'FontName', 'serif', 'Location', 'northwest');
xlabel('QOI Value', 'FontSize', 12, 'FontName', 'serif');
ylabel('Frequency', 'FontSize', 12, 'FontName', 'serif');
title('Histogram Comparison of QOI Distributions', 'FontSize', 14, 'FontName', 'serif');
xlim([0 2e11]);
hold off;

%update parameters function
function p = updatePars(p, parsName, parsValue)
    for ii = 1:length(parsName)
        p.(parsName{ii}) = parsValue(ii);
    end
end
