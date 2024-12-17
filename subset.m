% subset.m (with validation check)
function random_subset = subset(data, subset_size)
    rng(59); % Set random seed
    if subset_size > length(data)
        error('Subset size is greater than the data size');
    end
    random_indices = randperm(length(data), subset_size);  % Get random indices
    random_subset = data(random_indices);  % Return the random subset
end
