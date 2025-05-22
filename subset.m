% this function was made to take a random subset of any bigger set. it was helpful in plotting QOI distributions
% as the Sobol main file ends with a sample of (2N+1)*base samples QOI values, where N is # of parameters
% so it was necessary to restrict the QOI distribution to a random subset in order to compare variability
% produced by a different sets of parameters. written by Kyle Adams

function random_subset = subset(data, subset_size, seed)
    rng(seed); %set random seed
    random_indices = randperm(length(data), subset_size);  %get random indices
    random_subset = data(random_indices);  %return the random subset
end