% sobol file to run the Sobol sensitivity analysis
% on the QSP ODE model. 
% This code originally stems from Dr. Jaimit Parikh, and we modified it to
% be proper for our purpose.

%%
% 1. Get parameters of the model
parameters = parameters();


%%
% 2. Sample using Sobol sampling the parameter values and running
% simulations for the samples

% Sampling for the desired parameters

parsObj.name = ["aLA", "dA", ...
                "sR", "dR", "aIR", "bIR", ...
                "aCI", "bCI", "aHI", "bHI", "lC", "gC", "KC", "bIC", "lH", "gH", "KH", "bIH", "lR", "dI" ...
                "aHC", "bHC", "dC", ...
                "aAH", "bAH", "aRA", "bRA", "aIRA", "bIRA", "dH", ...
                "dL", "aCL", "bCL"];
            

parsObj.dist = {'Uniform', 'Uniform' , 'Uniform',...
    'Uniform', 'Uniform', 'Uniform',...
    'Uniform','Uniform', 'Uniform',...
    'Uniform', 'Uniform', 'Uniform',...
    'Uniform','Uniform', 'Uniform',...
    'Uniform', 'Uniform', 'Uniform',...
    'Uniform','Uniform', 'Uniform',...
    'Uniform', 'Uniform', 'Uniform',...
    'Uniform','Uniform', 'Uniform',...
    'Uniform', 'Uniform', 'Uniform',...
    'Uniform','Uniform', 'Uniform'};

nParameters = length(parsObj.name);

parsObj.parameters = {
    {'lower', 0.5,      'upper', 1.5},... 1%aLA
    {'lower', 0.0416665,   'upper', 0.1249995},... 2%dA;
    {'lower', 0.025,     'upper', 0.075},... 3%sR
    {'lower', 0.03290,       'upper', 0.0987},... 4%dR
    {'lower', 0.304,       'upper', 0.912},...    5%aIR
    {'lower', 0.004166665,      'upper', 0.012499995},... 6%bIR
    {'lower', 0.18,     'upper', 0.54},... 7%aCI
    {'lower', 175,   'upper', 525},... 8%bCI
    {'lower', 35.35,       'upper', 106.5},... 9%aHI
    {'lower', 50,       'upper', 150},... 10%bHI
    {'lower', 0.0025,   'upper', 0.0075},... 11%IC
    {'lower', 1.0395,     'upper', 3.1185},... 12%gC
    {'lower', 299.1,       'upper', 897.3},... 13%KC
    {'lower', 0.089,     'upper', 0.267},... 14%bIC
    {'lower', 0.000002999994,   'upper', 0.000008999982},... 15%IH
    {'lower', 0.756,     'upper', 2.268},... 16%gH
    {'lower', 211.2,       'upper', 633.6},... 17%KH
    {'lower', 0.089,       'upper', .267},...    18%bIH
    {'lower', 0.000001666666,      'upper', 0.000004999998},... 19%IR
    {'lower', 83.1775,        'upper', 249.5325},... %20dI
    {'lower', 0.5,     'upper', 1.5},... %21aHC
    {'lower', 176,       'upper', 528},...    %22bHC
    {'lower', 0.2926829269,       'upper', 0.8780487806},...%23dC
    {'lower', 0.00001305,     'upper', 0.00003915},...%24aAH
    {'lower', 2,     'upper', 6},... %25bAH   
    {'lower', 0.2,     'upper', 0.6},...%26aRA
    {'lower', 10,     'upper', 30},... %27bRA
    {'lower', 1,       'upper', 3},... %28aIRA   
    {'lower', 178,      'upper', 534},... %29bIRA
    {'lower', 0.16665,   'upper', 0.49995},... %30dH
    {'lower', 0.025,      'upper', 0.075},... %31dL
    {'lower', 5,     'upper', 15},... %32aCL
    {'lower', 100,       'upper', 300}};  %33bCL



size(parsObj.parameters)
parsObj.lb = num2cell(repmat(-inf, 1, length(parsObj.name)));
parsObj.ub = num2cell(inf(1, length(parsObj.name)));
parsObj.N=1000; 
samples = getSamplesSobol(parsObj, false);

%%
% 3. Simulating model for the desired parameter set and extracting 
% QOI cells at the end of 30 days

parsName = parsObj.name;
t0 = 0; tfinal = 30;
QOI = zeros(1, length(samples));
pN = cell(1,length(samples));
parfor ii = 1:length(samples)
    pN{ii} = updatePars(parameters, parsName, samples(ii, :));
    QOI(ii) = qoi(pN{ii});
end
QOICells = QOI ;

size(QOI)
size(samples)
%%
% 4. Scatter Plot of QOI vs selected parameters
 plotScatter(samples, QOICells, parsName, ...
    'scatterSobol.png');
ylabel('QOI')

%%
% 5. Estimating the finite difference based
Y = QOI';
S = sobolAnalysis(Y, length(parsName), parsObj.N, 1000, 0.95);

Si = sobolAnalysis(Y, length(parsObj.name), parsObj.N);

 mytable = cell2table([Si.S1, Si.ST], 'VariableNames', ["S1", "ST"],...
    'RowNames',parsObj.name); 
sortTable = sortrows(mytable, 1, 'descend'); 
disp(sortTable)

writetable(sortTable, 'sortedSobol2.csv', 'WriteRowNames',true)




f=figure('DefaultAxesFontSize', 16);
subplot(1,2,1)
bar(cell2mat(S.S1), 'FaceColor', 'k');
ylabel('S1', 'FontSize', 16);

subplot(1,2,2)
bar(cell2mat(S.ST), 'FaceColor', 'k');
ylabel('ST', 'FontSize', 16);


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
