% Calculates distribution of individuals who go on ART by disease and viral
% load at time of therapy initiation. Obtained by averaging the
% distribution over the past 20 time steps. When less than 20 time steps
% have elapsed, the distribution is obtained by averaging over all the 
% elapsed time steps.
% Accepts artTreat, a matrix describing the distrubiton of individuals who
% go on ART in the current time step.
% Returns artDist, a matrix describing the distribution of individuals who
% went on ART averaged over the past 20 time steps.
function[artDist] = calcDist(artDistList , disease , viral , gender , age , ...
    risk)
s = zeros(prod([disease , viral , gender , age , risk]) , 1); % initialize sum variable
sumall = @(x) sum(x(:)); % helper function
for i = 0 : artDistList.size() - 1 % Note: Java uses 0-based indexing
    dist = reshape(artDistList.get(i) , size(s));
    if sumall(dist) ~= 0
        s = s + dist ./ sumall(dist(:));
    end
end
artDist = s ./ artDistList.size(); % size(artDistList) <= 20
% timeSteps = size(artTreat , 6);
% artDist = sum(artTreat , timeSteps) / timeSteps;
% artDist = reshape(artTreat , [prod(dim) , 1]);
  


  