% ART discontinuation and initiation to match population-level coverage
% Calculate initiation and discontinuation matrices to maintain
% population-level ART coverage

function [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
    excMinAges , popCoverInd , g , risk , ...
    ageHIVSubTots , ageARTSubTots , maxAges , minAges , ...
    fracARTAge , minCoverLim , maxCoverLim , ageSexDebut , ...
    agePopSubTots , dRange)

sumall = @(x) sum(x(:));

% Calculate population-level ART coverage
fracARTAge_adj = fracARTAge; % baseline matrix
fracARTAge_adj(1 , minAges) = minCoverLim; % adjust for ages brought up to min coverage
fracARTAge_adj(1 , maxAges) = maxCoverLim; % adjust for ages brought down to max coverage
fracART = sum(fracARTAge_adj(1 , ageSexDebut:end) .* ...
    (agePopSubTots(1 , ageSexDebut:end) ./ sum(agePopSubTots(1 , ageSexDebut:end)))); % calculate weighted average of adjusted age-level coverages

fracExcMinAges = ageHIVSubTots(excMinAges) ./ sum(ageHIVSubTots(excMinAges)); % population HIV proportion by age, for ages > MIN coverage
toRedisMinAges = sum(ageHIVSubTots(minAges)); % summed HIV-positives of ages < MIN to redistribute
% Discontinue persons to match population coverage
if fracART > popCoverInd
    cover = (fracART - popCoverInd) ./ (fracART); % discontinuation fraction needed to meet population-level coverage
    totToRedis = fracExcMinAges .* sumall(toRedisMinAges .* max(cover , 0)); % calculate discontinuation #s, distributing persons that would have been from ages < MIN
    fracToRedis = totToRedis ./ ageHIVSubTots(excMinAges); % calculate discontinuation proportions
    formatToRedis = ones(1 , length(excMinAges) , 1);
    formatToRedis(:) = (max(cover , 0) + fracToRedis'); % set up matrix values by age
    formatModifier = ones(1 , length(excMinAges) , 1);
    formatModifier(:) = ageARTSubTots ./ ageHIVSubTots;
    % Discontinuation matrix by gender, age, risk, accounting for min/max adjustment
    % adjustment factor = ((HIVpop+#discont)/HIVpop) = (1 + #discont/HIVpop) 
    % = (1 + ((artOut*ARTpop)/HIVpop)) = (1 + artOut * (ARTpop/HIVpop))
    artOut(g , excMinAges , :) = artOut(g , excMinAges , :) + ...
        bsxfun(@times , ones(1 , length(excMinAges) , risk) , formatToRedis) .* ...
        (1 + artOut(g , excMinAges , :) .* ...
        bsxfun(@times , ones(1 , length(excMinAges) , risk) , formatModifier));
end

fracExcMaxAges = ageHIVSubTots(excMaxAges) ./ sum(ageHIVSubTots(excMaxAges)); % population HIV proportion by age, for ages < MAX coverage
toRedisMaxAges = sum(ageHIVSubTots(maxAges)); % summed HIV-positives of ages > MAX coverage to redistribute
% Initiate persons to match population coverage
if fracART < popCoverInd
    cover = (popCoverInd - fracART) ./ (1 - fracART); % initiation fraction needed to meet population-level coverage
    totToRedis = fracExcMaxAges .* sumall(toRedisMaxAges .* max(cover , 0)); % calculate initiation #s, distributing persons that would have been from ages > MAX
    fracToRedis = totToRedis ./ ageHIVSubTots(excMaxAges); % calculate initiation proportions
    formatToRedis = ones(1 , 1 , 1 , length(excMaxAges) , 1);
    formatToRedis(:) = (max(cover , 0) + fracToRedis'); % set up matrix values by age
    % Initiation matrix by disease, VL, gender, age, risk, accounting for min/max adjustment
    % adjustment factor = ((HIVpop-#init)/HIVpop) = = (1 - #init/HIVpop) = (1 - treat)
    treat(dRange , 1 : 5 , g , excMaxAges , :) = treat(dRange , 1 : 5 , g , excMaxAges , :) + ...
        bsxfun(@times , ones(length(dRange) , 5 , 1 , length(excMaxAges) , risk) , formatToRedis) .* ...
        (1 - treat(dRange , 1 : 5 , g , excMaxAges , :)); 
end
