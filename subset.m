% this function was made to take a random subset of any bigger set. it was helpful in plotting QOI distributions
% as the Sobol main file ends with a sample of (2N+1)*base samples QOI values, where N is # of parameters
% so it was necessary to restrict the QOI distribution to a random subset in order to compare variability
% produced by a different sets of parameters 


% subset.m (with validation check)
function random_subset = subset(data, subset_size)
    rng(59); % Set random seed
    if subset_size > length(data)
        error('Subset size is greater than the data size');
    end
    random_indices = randperm(length(data), subset_size);  % Get random indices
    random_subset = data(random_indices);  % Return the random subset
end
