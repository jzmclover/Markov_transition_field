function [X_mtf] = markovTransitionField(X)
% Calculate Markov Transition Field (MTF) matrix for a univariate time serie X.
% Input 
%       X: 1 * time_length
% Output 
%       X_mtf: time_length * time_length

    % Quantify time series
    len = length(X);
    n_bins = 8;
    quantizedX = uniformBinDivision(X, n_bins);
    
    % Build transition matrix
    X_mtm = markov_transition_matrix(quantizedX, len, n_bins);
    
    % Normalize transition matrix
    sum_mtm = sum(X_mtm,2);
    sum_mtm(sum_mtm==0) = 1;
    X_mtm = X_mtm ./ sum_mtm;
    
    % Generate MTF matrix
    X_mtf = markov_transition_field(quantizedX, X_mtm, len);
    
end

function [X_mtm] = markov_transition_matrix(X_binned, n_timestamps, n_bins)
    X_mtm = zeros(n_bins, n_bins);
    for j = 1:n_timestamps - 1
        X_mtm(X_binned(j), X_binned(j + 1)) = X_mtm(X_binned(j), X_binned(j + 1)) + 1;
    end
end

function [X_mtf] = markov_transition_field(X_binned, X_mtm, n_timestamps)
    X_mtf = zeros(n_timestamps, n_timestamps);
    for j = 1:n_timestamps
        for k = 1:n_timestamps
            X_mtf(j, k) = X_mtm(X_binned(j), X_binned(k));
        end
    end
end

function binIndices = uniformBinDivision(array, N)
    % Sort the array
    sortedArray = sort(array);
    
    % Calculate the number of elements that each bin should contain
    elementsPerBin = ceil(length(sortedArray) / N);
    
    % Create edge values for each bin
    binEdges = sortedArray(1:elementsPerBin:end);
    if length(binEdges) < N + 1
        binEdges = [binEdges; sortedArray(end)];
    end
    
    % Determine which bin each element belongs to via discretize function
    binIndices = discretize(array, binEdges);
    
    % If an element falls outside the last bin, it will be included in the last bin
    binIndices(binIndices == N+1) = N;
end

