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
    {'lower', 0.0000226,      'upper', 0.0000678},... %1lL
    {'lower', 0.04165,   'upper', 0.12495},... %2dA;
    {'lower', 0.0535,     'upper', 0.1605},... %3sR
    {'lower', 0.0329,       'upper', 0.0987},... %4dR
    {'lower', 0.3125,       'upper', 0.9375},...    %5aIR
    {'lower', 0.004165,      'upper', 0.012495},... %6bIR
    {'lower', 0.18,     'upper', 0.54},... %7aCI
    {'lower', 176,   'upper', 528},... %8bCI
    {'lower', 35.35,       'upper', 106.5},... %9aHI
    {'lower', 49.85,       'upper', 149.55},... %10bHI
    {'lower', 0.0005,   'upper', 0.0015},... %11lC
    {'lower', 1.04,     'upper', 3.12},... %12gC
    {'lower', 299,       'upper', 897},... %13KC
    {'lower', 1,       'upper', 3},... %14aIC
    {'lower', 0.089,     'upper', 0.267},... %15bIC
    {'lower', 0.00000315,   'upper', 0.00000945},... %16lH
    {'lower', 0.755,     'upper', 2.265},... %17gH
    {'lower', 211,       'upper', 633},... %18KH
    {'lower', 1,       'upper', 3},... %19aIH
    {'lower', 0.089,       'upper', .267},...    %20bIH
    {'lower', 0.0000003335,      'upper', 0.0000010005},... %21lR
    {'lower', 83,        'upper', 249},... %22dI
    {'lower', 0.5,     'upper', 1.5},... %23aHC
    {'lower', 17.5,       'upper', 52.5},...    %24bHC
    {'lower', 0.2925,       'upper', 0.8775},...%25dC
    {'lower', 0.00001305,     'upper', 0.00003915},...%26aAH
    {'lower', 2,     'upper', 6},... %27bAH   
    {'lower', 0.2,     'upper', 0.6},...%28aRA
    {'lower', 10,     'upper', 30},... %29bRA
    {'lower', 1,       'upper', 3},... %30aIRA   
    {'lower', .178,      'upper', .534},... %31bIRA
    {'lower', 0.16665,   'upper', 0.49995},... %32dH
    {'lower', 0.0025,      'upper', 0.0075},... %33dL
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

%saving data to send

writematrix(QOICells, 'ALLCELLS.txt')
writetable(mytable, "SAValues.txt")
