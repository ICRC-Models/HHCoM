% ART discontinuation and initiation to maintain coverage limits
% Calculate initiation and discontinuation matrices to maintain
% ART coverage minimums and maximums by age

function [artOut , treat , maxAges , excMaxAges , minAges , excMinAges] = ...
    artMinMax(artOut , treat, minCoverLim , maxCoverLim , ageFracART , ageVec , ...
    ageHIVallSubTots , ageHIVeligSubTots , g , risk , ageSexDebut , dRange)

% Discontinue persons from age groups with coverage > maxCover
maxInds = ageFracART > maxCoverLim; % find inds of ages above max coverage
if any(maxInds)
    maxAges = ageVec(maxInds == 1); % ages above max coverage
    excMaxAges = ageVec(~maxInds); % ages that are good, ie below max coverage
    excMaxAges = excMaxAges(excMaxAges > 2); % exclude ages ineligible for ART
    coverMax = (ageFracART - maxCoverLim) ./ ageFracART; % discontinuation fraction needed to meet max coverage
    formatMax = ones(1 , length(maxAges) , 1); 
    formatMax(:) = coverMax(maxAges); % set up matrix values by age
    artOut(g , maxAges , :) = artOut(g , maxAges , :) + ...
        bsxfun(@times , ones(1 , length(maxAges) , risk) , formatMax); % discontinuation matrix by gender, age, risk
else
    maxAges = [];
    excMaxAges = ageVec(ageSexDebut:end);
end
% Initiate persons in age groups with coverage < minCover
minInds = ageFracART < minCoverLim; % find inds of ages below min coverage
if sum(minInds) > 2
    minAges = ageVec(minInds == 1); % ages below min coverage
    minAges = minAges(minAges > 2); % exclude ages ineligible for ART
    excMinAges = ageVec(~minInds); % ages that are good, ie above min coverage
    coverMin = ((minCoverLim - ageFracART) ./ (1 - ageFracART)) .* ... % initiation needed to meet min coverage
        (ageHIVallSubTots ./ ageHIVeligSubTots); % adjust coverage b/c only initiating persons with eligible CD4 = ((coverMin*HIVtot)/HIVelig)
    formatMin = ones(1 , 1 , 1 , length(minAges) , 1);
    formatMin(:) = coverMin(minAges); % set up matrix values by age
    treat(dRange , 1 : 5 , g , minAges , :) = treat(dRange , 1 : 5 , g , minAges , :) + ...
        bsxfun(@times , ones(length(dRange) , 5 , 1 , length(minAges) , risk) , formatMin); % initiation matrix by disease, VL, gender, age, risk
else
    minAges = [];
    excMinAges = ageVec(ageSexDebut:end);
end




