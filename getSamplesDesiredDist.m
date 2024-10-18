function samples = getSamplesDesiredDist(paramObj, samples, sampleDist)
paramDistTypes = paramObj.dist;
lb = paramObj.lb;
ub = paramObj.ub;
distParamNamesValues = paramObj.parameters;

paramDists = cellfun(@(x, y, l, u)truncate(makedist(x, y{:}), l, u), ...
    paramDistTypes, ...
    distParamNamesValues, lb, ub, ...
    'UniformOutput',false);

samples = splitapply(...
    @(a,b,c)inverseTransformSampling(a, b{:}, c{:}), ...
    samples, sampleDist, paramDists,...
    1:size(samples, 2));
end


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

transformedSamples = desiredDist.icdf(sampleDist.cdf(samples));


end