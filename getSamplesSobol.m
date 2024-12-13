function [samples, sampleDist] = getSamplesSobol(paramObj, ...
    calcSecondOrder)
arguments
    paramObj;%note: this is called parsObj in SobolMain
    calcSecondOrder logical; % 
end
% getSamplesSobol provides the sobol samples using structured sampling
% Inputs:
%   paramObj; structure with parameter information
% Outputs:
%   samples: the number of returned samples = N(k+2) if only calculating
%   first order index or 2N(k+1) if also calculating second order index
% Author: Jaimit Parikh
% Date Last modified: 07-10-2023 (by Dr. Parikh)
% see also sobolset sgenerator makedist


N = paramObj.N; % Number of base samples - this will get multiplied by (K+2) for first-order S1
k = length(paramObj.name);
sgenerator = sobolset(k * 2); % doubling the number of samples to create 2 matrices (A & B) below
%the line above creates a higher-dimensional object of dimension k*2
samples = sgenerator(1:N, :); %Taking out the desired number of samples from the higher-dimensional object



samples = mat2cell(samples, N, [k, k]);%separates the matrix into a row of N kxk cells
samples = cell2mat(samples');%Transforms the cells into a column of N kxk matrices

%calls the function below beginning on line 44 to take the samples and generate the matrices and their products
[A, B, AB, BA] = createMatricesForSobolIndices(samples);

if ~calcSecondOrder
    samples = [A; B; AB];%puts matrices into desired object for second order
else
    samples = [A;B;AB;BA];%puts matrices into desired object for first order
end

sampleDist = cell(1, k);%makes an object of type cell with dimension 1xk
[sampleDist{:}] = deal(makedist('Uniform'));%makes a standard uniform distribution [0,1] and fills each cell with that
%Calls a function to take the parameter object along with the sobol matrices and the
%standard uniform distributions and return sobol matrices that are in the parameter space
samples = getSamplesDesiredDist(paramObj, samples, sampleDist);
end
%This function takes the NxKxK object and transforms it into sobol matrices
function [A, B, AB, BA] = createMatricesForSobolIndices(samples)
arguments
    samples = reshape(1:8, 2, [])'; % default matrix for test
end


N = length(samples) / 2;
nInputs = size(samples, 2);

% randomly permute the rows of the samples
% samples = samples(randperm(2*N), :);

% create matrix A, B and AB by splitting the samples
A = samples(1:N, :);
B = samples(N + 1:end, :);

Arepeat = repmat(A, nInputs, 1); % Matrix A is duplicated nInputs times

Bcell = num2cell(B, 1);
Bdiag = blkdiag(Bcell{:});
AB = Bdiag + Arepeat .* not(Bdiag);


Brepeat = repmat(B, nInputs, 1); % Matrix  B is duplicated nInputs times
Acell = num2cell(A, 1);
Adiag = blkdiag(Acell{:});
BA = Adiag + Brepeat .*not(Adiag);


end
