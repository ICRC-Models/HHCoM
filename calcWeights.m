function [weights] = calcWeights(paramSetMatrix , alphaSets_prev , alphaWeights_prev , pIdx)

%% Specify dimensions
numParam = size(alphaSets_prev , 1); % number parameters
numPart_prev = size(alphaSets_prev , 2); % number particles (t-1)
numPart = size(paramSetMatrix , 2); % number particles (t)

%% Calculate numerator
[paramsAll] = genParamStruct();
numerator = zeros(numPart,1);
for i = 1 : numPart
    d_prior = 1;
    for s = 1 : length(pIdx)
        d_prior = d_prior * prod(1 ./ (paramsAll{pIdx(s)}.ub - paramsAll{pIdx(s)}.lb));
    end
    numerator(i) = d_prior;
end

%% Calculate denominator
% Specify constant multiplier in Gaussian distribution
const = (1/sqrt(2 * pi));
  
% Create empty vector in which to store denominator values for particles in iteration t
denominator = zeros(numPart,1);
  
% For each particle simulated in iteration t, calculate the probability of moving to that particle from each accepted particle in iteration t-1
parfor i = 1 : numPart
    denominator_temp = 0.0;
    for j = 1 : numPart_prev
        kernel = 1;
        for p = 1 : numParam
            v = 2 .* var(alphaSets_prev(p,:)); % 2x variance of previous accepted particles (normalized by numPart-1)
            kernel = kernel * (const * (1 / sqrt(v)) * exp(-0.5 * ((1/v)*(paramSetMatrix(p,i) - alphaSets_prev(p,j))) ^ 2));
        end
        denominator_temp = denominator_temp + alphaWeights_prev(j) * kernel;
    end
    denominator(i) = denominator_temp;
end
  
%% Calculate weights
weights = (numerator./denominator);
