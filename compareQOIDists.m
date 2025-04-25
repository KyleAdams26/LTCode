% sobol file to run the Sobol sensitivity analysis
% on the QSP ODE model. 
% This code originally stems from Dr. Jaimit Parikh's SobolMain.m script.
% Kyle Adams added comments and adjusted this code to compare QOI
% distributions resulting from fixing some parameters.

%%changeable here is how you create your bounds. default is
%%lower bounds are 50% of the nominal value, and upper bounds are 150%
%% change me %%
lower_percentage = 0.5;
upper_percentage = 1.5;
base_samples = 175000;
param_dist = {'Uniform'};
%varying_params = {'aHI', 'dC', 'dI', 'aIH', 'gH', 'dH', 'bHI', 'bIH', 'bIC', 'KH', 'aHC', 'aCI', 'aIC', 'bHC', 'bCI', 'dA', 'lL', 'sR', 'bIR', 'aRA', 'aIRA', 'lR', 'aIR', 'bRA', 'bIRA', 'bAH', 'dR', 'aAH', 'lH'};
%numVarying = 29;
varying_params = {'aCL', 'dL', 'bCL', 'KC', 'gC', 'aIC'};
numVarying = 6;
%% ---------- %%

% 1. Get subset varying_params of parameters of the model from parameters.m
p = parameters();
p_subset = struct();
for i = 1:length(varying_params)
    field = varying_params{i};
    if isfield(p,field)
        p_subset.(field) = p.(field);
    end
end
disp(p_subset)
paramNames = fieldnames(p_subset);
numParam = length(paramNames);

lowBounds = zeros(1, numParam); %initializes vector of length numParam with 0s
upBounds = zeros(1, numParam);


%2a. Set up bounds
for i = 1:numParam
    param = paramNames{i};
    lowBounds(i) = p_subset.(param)*lower_percentage;
    upBounds(i) = p_subset.(param)*upper_percentage;
end

%storing this info in structure called parsObj, used in getSamplesSobol
parsObj.name = paramNames'; %transpose so the dimensions of each field match
parsObj.lb = num2cell(repmat(-inf, 1, length(parsObj.name)));
parsObj.ub = num2cell(inf(1, length(parsObj.name)));
parsObj.dist = repmat(param_dist, 1, numParam);
parsObj.N = base_samples; %how many sobol samples you want
%below is a field that stores each parameter's bounds
parsObj.parameters = arrayfun(@(i) {'lower', lowBounds(i), 'upper', upBounds(i)}, 1:numParam, 'UniformOutput', false);

samples = getSamplesSobol(parsObj, false);

%%
% 3. Simulating model for the desired parameter set and extracting 
% QOI cells at the end of 30 days

parsName = parsObj.name;
QOI = zeros(1, length(samples));
pN = cell(1,length(samples));
parfor ii = 1:length(samples)
    pN{ii} = updatePars(p, parsName, samples(ii, :));
    QOI(ii) = qoi(pN{ii});
end

%saving QOI with some fixed parameters for histogram comparison
%QOI_with_fixed_params = QOI;
QOI_6 = QOI;

%plot histograms of comparing QOI distributions and perform KS test
plotHistogram(numVarying, QOI_all_varying, QOI_with_fixed_params);
[h,p, ks2stat] = kstest2(QOI_all_varying, QOI_with_fixed_params);


%cdf plots
figure;
hold on;
c1 = cdfplot(QOI_with_fixed_params);
c2 = cdfplot(subset(QOI_all_varying, length(QOI_with_fixed_params)));
c1.LineWidth = 2;
c2.LineWidth = 2;
legend({'All Parameters Varied', sprintf('%d Parameters Varied', numVarying)}, ...
    'FontSize', 12, 'FontName', 'serif', 'Location', 'northwest');
hold off;

function p = updatePars(p, parsName, parsValue)

for ii = 1:length(parsName)
    p.(parsName{ii}) = parsValue(ii);
end

end

function plotHistogram(numVarying, QOI_all, QOI_some)
figure;
hold on;
subset_QOI_all = subset(QOI_all, 1400000);
histogram(subset_QOI_all,'EdgeColor', 'red', 'DisplayStyle', 'stairs', 'LineWidth', 2);
histogram(QOI_some,'EdgeColor', 'blue', 'DisplayStyle', 'stairs', 'LineWidth', 2);
legend({'All Parameters Varied', sprintf('%d Parameters Varied', numVarying)}, ...
    'FontSize', 12, 'FontName', 'serif', 'Location', 'northwest');
xlabel('QOI Value', 'FontSize', 12, 'FontName', 'serif');
ylabel('Frequency', 'FontSize', 12, 'FontName', 'serif');
title('Histogram Comparison', 'FontSize', 12, 'FontName', 'serif');
xlim([0 2e11]);
hold off;

end
%saving data

%writematrix(QOI, 'QOIvalues.txt') %saves QOI values
%writetable(mytable, "SAValues.txt") %saves sensitivity indices


