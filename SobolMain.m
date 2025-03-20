% sobol file to run the Sobol sensitivity analysis
% on the QSP ODE model. 
% This code originally stems from Dr. Jaimit Parikh, and we modified it to
% be proper for our purpose. Adjustments and comments were added by Kyle Adams.

%%changeable here is how you create your bounds. default is
%%lower bounds are 50% of the nominal value, and upper bounds are 150%
%% change me %%
lower_percentage = 0.5;
upper_percentage = 1.5;
base_samples = 100000;
param_dist = {'Uniform'};
%% ---------- %%

% 1. Get parameters of the model from parameters.m
p = parameters();
paramNames = fieldnames(p);
numParam = length(paramNames);

lowBounds = zeros(1, numParam); %initiliazes vector of length numParam with 0s
upBounds = zeros(1, numParam);
% 2. Set up bounds for parameters 
for i = 1:numParam %populates lowBounds and upBounds with bounds listed in parametersBetter
    param = paramNames{i};
    lowBounds(i) = p.(param)*lower_percentage; %can adjust percent change of parameter here
    upBounds(i) = p.(param)*upper_percentage;
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

%stores QOI distribution for histogram comparison later
QOI_all_varying = QOI;

%%
% 4. Scatter Plot of QOI vs selected patrameters
 plotScatter(samples, QOI, parsName, ...
    'scatterSobol.png');
ylabel('QOI')

%%
% 5. Estimating the Sobol Index values
Y = QOI';
S = sobolAnalysis(Y, length(parsName), parsObj.N, 1, 0.95);

Si = sobolAnalysis(Y, length(parsObj.name), parsObj.N);

 mytable = cell2table([Si.S1, Si.ST], 'VariableNames', ["S1", "ST"],...
    'RowNames',parsObj.name); 
sortTable = sortrows(mytable, 1, 'descend'); 
disp(sortTable)

writetable(sortTable, 'sortedSobol2.csv', 'WriteRowNames',true)


%preparing for figure in order of descending total sensitivity indices
S1 = cell2mat(S.S1);
ST = cell2mat(S.ST);
[ST_sorted, idx] = sort(ST, 'descend');
S1_sorted = S1(idx);
paramNames_sorted = paramNames(idx);

%plotting sensitivity indices
hold on;
figure('DefaultAxesFontSize', 16);
b = bar([S1_sorted, ST_sorted], 'grouped');
b(1).FaceColor = [1 0.79 0.63];
b(2).FaceColor =  [0 0.188 0.69];
set(gca, 'FontName', 'Times New Roman')
xticks(1:length(paramNames_sorted));
xticklabels(paramNames_sorted)
xtickangle(45);
ylabel('Sensitivity Index', 'FontSize', 16, 'FontName', 'serif');
legend({'S1', 'ST'}, 'FontSize', 16, 'FontName', 'serif');
hold off;


%gives a scatterplot of each parameter's value vs what QOI value it
%produced
function plotScatter(samples, QOI,  parsName, fname) 
f = figure('DefaultAxesFontSize', 14);
tiledlayout('flow')
for ii = 1:size(samples, 2)
    nexttile;
    plot(samples(:, ii), QOI, 'ko'); xlabel(parsName(ii));
    ylabel('QOI');
end
exportgraphics(f, fname, 'resolution', 300);
end

function p = updatePars(p, parsName, parsValue)

for ii = 1:length(parsName)
    p.(parsName{ii}) = parsValue(ii);
end

end

%saving data

%writematrix(QOI, 'QOIvalues.txt') %saves QOI values
%writetable(mytable, "SAValues.txt") %saves sensitivity indices


