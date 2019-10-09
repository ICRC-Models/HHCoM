% Calculates distribution of individuals who go on ART by disease and viral
% load at time of therapy initiation. Obtained by averaging the
% distribution over the past stepsPerYear*2 time steps. When less than stepsPerYear*2 time steps
% have elapsed, the distribution is obtained by averaging over all the 
% elapsed time steps.
% Accepts artDistList, a linked list of matrices describing the distributon of individuals who
% go on ART in the past stepsPerYear*2 time steps.
% Returns artDist, a matrix describing the distribution of individuals who
% went on ART averaged over the past stepsPerYear*2 time steps.
function[artDist] = calcDist(artDistList , disease , viral , gender , age , ...
    risk , sumall)
s = zeros(prod([disease , viral , gender , age , risk]) , 1); % initialize sum variable as a row vector
for i = 0 : artDistList.size() - 1 % Note: Java uses 0-based indexing
    dist = reshape(artDistList.get(i) , size(s)); % reshape artDistList at time i to match the shape of s
    if sumall(dist) ~= 0
        s = s + dist ./ sumall(dist(:)); % sum distribution going on ART over all time steps
    end
end
artDist = s ./ artDistList.size(); % divide by number of time steps: size(artDistList) <= 20
