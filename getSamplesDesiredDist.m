function samples = getSamplesDesiredDist(paramObj, samples, sampleDist)
%Author: Jaimit Parikh
%Comments by Skylar Grey

%Separate out parameter distribution types from parameter object
paramDistTypes = paramObj.dist;
%separate out lower bounds for truncation function
lb = paramObj.lb;
%separate out upper bounds for truncation function
ub = paramObj.ub;
%separate out upper and lower bounds for parameter distributions
distParamNamesValues = paramObj.parameters;

%cellfun applies the function to each cell of the array one at a time
%in other words it loops through each of the parameters without a loop
%makedist makes the distributions of the specified types with the specified
%properties, while truncate cuts those distributions off at the specified
%upper and lower bounds (not needed for all uniform distributions, however
%this code is general enough to work for other distribution types)
paramDists = cellfun(@(x, y, l, u)truncate(makedist(x, y{:}), l, u), ...
    paramDistTypes, ...
    distParamNamesValues, lb, ub, ...
    'UniformOutput',false);
%splitapply splits each object into groups (in this case along columns)
%and applies the function to each group one at a time
%In this case, it applies the user-written function beginning on line 38
%which transforms the sobol matrices into the parameter space
samples = splitapply(...
    @(a,b,c)inverseTransformSampling(a, b{:}, c{:}), ...
    samples, sampleDist, paramDists,...
    1:size(samples, 2));
end

%inverseTransformSampling is passing elements of samples, and unpacking sampleDist and paramDists
% for each column of samples, inverseTransform sampling is transforming the distribution
% of the column from a uniform distribution to a desired distribution, specified by paramDists.

%this function is called by the code in line 25, see line 41 for description of function
function transformedSamples = inverseTransformSampling(samples, ...
    sampleDist, desiredDist)
arguments
    samples;
    sampleDist;
    desiredDist;
end
%% transformedSamples = inverseTransformSampling(samples, sampleDist,
% desiredDist) returns the transformedSamples with the desired
% distribution given orginal samples sampled from sampleDist
%samples containes the sobol matrices and their products
%this function takes a uniform distribution with min 0 and max 1
%and applies its cdf to the sobol matrices, then it does an inverse
%cdf with the parameter distributions. this transforms the samples into
%the parameter space.
transformedSamples = desiredDist.icdf(sampleDist.cdf(samples));

end
