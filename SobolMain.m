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


parsObj.name = ["lL", "dA", ...
                "sR", "dR", "aIR", "bIR", ...
                "aCI", "bCI", "aHI", "bHI", "lC", "gC", "KC", "aIC", "bIC", "lH", "gH", "KH", "aIH", "bIH", "lR", "dI" ...
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
    'Uniform','Uniform', 'Uniform',...
    'Uniform', 'Uniform'};

nParameters = length(parsObj.name);

parsObj.parameters = {
    {'lower', 0.5,      'upper', 1.5},... %1lL
    {'lower', 0.0416665,   'upper', 0.1249995},... %2dA;
    {'lower', 0.025,     'upper', 0.075},... %3sR
    {'lower', 0.03290,       'upper', 0.0987},... %4dR
    {'lower', 0.304,       'upper', 0.912},...    %5aIR
    {'lower', 0.004166665,      'upper', 0.012499995},... %6bIR
    {'lower', 0.18,     'upper', 0.54},... %7aCI
    {'lower', 175,   'upper', 525},... %8bCI
    {'lower', 35.35,       'upper', 106.5},... %9aHI
    {'lower', 50,       'upper', 150},... %10bHI
    {'lower', 0.0025,   'upper', 0.0075},... %11lC
    {'lower', 1.0395,     'upper', 3.1185},... %12gC
    {'lower', 299.1,       'upper', 897.3},... %13KC
    {'lower', 1,       'upper', 3},... %14aIC
    {'lower', 0.089,     'upper', 0.267},... %15bIC
    {'lower', 0.000002999994,   'upper', 0.000008999982},... %16lH
    {'lower', 0.756,     'upper', 2.268},... %17gH
    {'lower', 211.2,       'upper', 633.6},... %18KH
    {'lower', 1,       'upper', 3},... %19aIH
    {'lower', 0.089,       'upper', .267},...    %20bIH
    {'lower', 0.000001666666,      'upper', 0.000004999998},... %21IR
    {'lower', 83.1775,        'upper', 249.5325},... %22dI
    {'lower', 0.5,     'upper', 1.5},... %23aHC
    {'lower', 176,       'upper', 528},...    %24bHC
    {'lower', 0.2926829269,       'upper', 0.8780487806},...%25dC
    {'lower', 0.00001305,     'upper', 0.00003915},...%26aAH
    {'lower', 2,     'upper', 6},... %27bAH   
    {'lower', 0.2,     'upper', 0.6},...%28aRA
    {'lower', 10,     'upper', 30},... %29bRA
    {'lower', 1,       'upper', 3},... %30aIRA   
    {'lower', 178,      'upper', 534},... %31bIRA
    {'lower', 0.16665,   'upper', 0.49995},... %32dH
    {'lower', 0.025,      'upper', 0.075},... %33dL
    {'lower', 5,     'upper', 15},... %34aCL
    {'lower', 100,       'upper', 300}};  %35bCL



size(parsObj.parameters)
parsObj.lb = num2cell(repmat(-inf, 1, length(parsObj.name)));
parsObj.ub = num2cell(inf(1, length(parsObj.name)));
parsObj.N=100; 
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
S = sobolAnalysis(Y, length(parsName), parsObj.N, 100, 0.95);

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
