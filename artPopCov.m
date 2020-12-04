% ART discontinuation and initiation to match population-level coverage
% Calculate initiation and discontinuation matrices to maintain
% population-level ART coverage

function [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
    excMinAges , popCoverInd , g , risk , ...
    ageHIVallSubTots , ageHIVeligSubTots , ageARTSubTots , maxAges , minAges , ...
    fracARTAge , minCoverLim , maxCoverLim , ageSexDebut , ...
    agePopSubTots , dRange)

% Calculate population-level ART coverage
fracARTAge_adj = fracARTAge; % baseline matrix
fracARTAge_adj(1 , minAges) = minCoverLim; % adjust for ages brought up to min coverage
fracARTAge_adj(1 , maxAges) = maxCoverLim; % adjust for ages brought down to max coverage
fracART = sum(fracARTAge_adj(1 , ageSexDebut+1:end) .* ...
    (agePopSubTots(1 , ageSexDebut+1:end) ./ sum(agePopSubTots(1 , ageSexDebut+1:end)))); % calculate weighted average of adjusted age-level coverages

fracExcMinAges = ageARTSubTots(excMinAges) ./ sum(ageARTSubTots(excMinAges)); % population ART proportion by age, for ages > MIN coverage
toRedisMinAges = sum(ageARTSubTots(minAges)); % summed HIV-positives on ART of ages < MIN to redistribute
% Discontinue persons to match population coverage
% if fracART > popCoverInd
%     cover = (fracART - popCoverInd) ./ (fracART); % discontinuation fraction needed to meet population-level coverage
%     totToRedis = fracExcMinAges .* sumall(toRedisMinAges .* max(cover , 0)); % calculate discontinuation #s, distributing persons that would have been from ages < MIN
%     fracToRedis = totToRedis ./ ageARTSubTots(excMinAges); % calculate discontinuation proportions
%     formatToRedis = ones(1 , length(excMinAges) , 1);
%     formatToRedis(:) = (max(cover , 0) + fracToRedis'); % set up matrix values by age
%     % Discontinuation matrix by gender, age, risk, accounting for min/max adjustment
%     % adjustment factor = ((ARTpop-#disc)/ARTpop) = (1 - #disc/ARTpop) = (1 - artOut)
%     artOut(g , excMinAges , :) = artOut(g , excMinAges , :) + ...
%         bsxfun(@times , ones(1 , length(excMinAges) , risk) , formatToRedis) .* ...
%         (1 - artOut(g , excMinAges , :));
% end

fracExcMaxAges = ageHIVeligSubTots(excMaxAges) ./ sum(ageHIVeligSubTots(excMaxAges)); % population HIV proportion by age, for ages < MAX coverage
toRedisMaxAges = ageHIVeligSubTots(maxAges); % HIV-positives of ages > MAX coverage to redistribute
% Initiate persons to match population coverage
if fracART < popCoverInd
    cover = ((popCoverInd - fracART) ./ (1 - fracART)) .* ... % initiation fraction needed to meet population-level coverage
        (ageHIVallSubTots ./ ageHIVeligSubTots); % adjust coverage b/c only initiating persons with eligible CD4 = ((cover*HIVtot)/HIVelig)
    totToRedis = fracExcMaxAges .* sumall(toRedisMaxAges .* max(cover(maxAges) , 0)); % calculate initiation #s, distributing persons that would have been from ages > MAX
    fracToRedis = totToRedis ./ ageHIVeligSubTots(excMaxAges); % calculate initiation proportions
    formatToRedis = ones(1 , 1 , 1 , length(excMaxAges) , 1);
    formatToRedis(:) = (max(cover(excMaxAges) , 0)' + fracToRedis'); % set up matrix values by age
    % Initiation matrix by disease, VL, gender, age, risk, accounting for min/max adjustment
    % adjustment factor = ((HIVpop-#init)/HIVpop) = (1 - #init/HIVpop) = (1 - treat)
    treat(dRange , 1 : 5 , g , excMaxAges , :) = treat(dRange , 1 : 5 , g , excMaxAges , :) + ...
        bsxfun(@times , ones(length(dRange) , 5 , 1 , length(excMaxAges) , risk) , formatToRedis) .* ...
        (1 - treat(dRange , 1 : 5 , g , excMaxAges , :)); 
end
