% HPV vaccination in a limited vaccine scenario
function[dPop , vaxRemain] = hpvVaxLmtd(pop , year , vaxLimitPerYr , ...
    disease , viral , risk , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , vaxCoverL , vaxRemain , vaxGL , toInd)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
hpvVaxd = 0;

%% If within first vaxLimitYrs-many vaccine-limited years
if rem(year,1) == 0.0    % reset vaxRemain at start of each new year to the number of available vaccines per year
    vaxRemain = vaxLimitPerYr;
end

fromNonVSusL_noScrn = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ...
    vaxGL , vaxAgeL , 1 : risk));
fromNonVSusL_scrn = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 : 3 , 3 , ...
    vaxGL , vaxAgeL , 1 : risk));
fromNonVImmL_noScrn = toInd(allcomb(1 : disease , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ...
    vaxGL , vaxAgeL , 1 : risk));
fromNonVImmL_scrn = toInd(allcomb(1 : disease , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 : 3 , 3 ,...
    vaxGL , vaxAgeL , 1 : risk));
toVSusL_noScrn = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 : 3 , 2 , ...
    vaxGL , vaxAgeL , 1 : risk));
toVSusL_scrn = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 : 3 , 4 , ...
    vaxGL , vaxAgeL , 1 : risk));
toVImmL_noScrn = toInd(allcomb(1 : disease , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 : 3 , 2 , ...
    vaxGL , vaxAgeL , 1 : risk));
toVImmL_scrn = toInd(allcomb(1 : disease , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 : 3 , 4 , ...
    vaxGL , vaxAgeL , 1 : risk));
otherVL = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , [2,4] , vaxGL , vaxAgeL , 1 : risk));
allVNonVL = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , vaxGL , vaxAgeL , 1 : risk));

fracVaxd = sumall(pop(otherVL)) / ...
    sumall(pop(allVNonVL)); % find proportion of population that is currently vaccinated
if vaxCoverL - fracVaxd > 10 ^ -6 % when proportion vaccinated is below target vaccination level
    vaxCover = max(0 , (vaxCoverL - fracVaxd) ./ (1 - fracVaxd)); % vaccinate enough people in age group to reach target
    vaxdGroupSum = (sumall(pop(fromNonVSusL_noScrn)) + ...
        sumall(pop(fromNonVSusL_scrn)) + sumall(pop(fromNonVImmL_noScrn)) + ...
        sumall(pop(fromNonVImmL_scrn)));
    if vaxRemain >= 1    % when vaccines remain
        if (vaxRemain >= vaxdGroupSum)    % when remaining vaccines >= targeted coverage
            usedVax = vaxdGroupSum;
        else    % when remaining vaccines < targeted coverage
            usedVax = min(vaxdGroupSum, vaxRemain);    
            vaxdGroupSus_noScrn = (usedVax .* (pop(fromNonVSusL_noScrn) ./ vaxdGroupSum))...
                .* (pop(fromNonVSusL_noScrn) > 0);    % proportionately divide available vaccines across compartments
            vaxdGroupSus_scrn = (usedVax .* (pop(fromNonVSusL_scrn) ./ vaxdGroupSum))...
                .* (pop(fromNonVSusL_scrn) > 0);
            vaxdGroupImm_noScrn = (usedVax .* (pop(fromNonVImmL_noScrn) ./ vaxdGroupSum))...
                .* (pop(fromNonVImmL_noScrn) > 0);
            vaxdGroupImm_scrn = (usedVax .* (pop(fromNonVImmL_scrn) ./ vaxdGroupSum))...
                .* (pop(fromNonVImmL_scrn) > 0);
        end
        vaxRemain = vaxRemain - usedVax;
    else 
        vaxdGroupSus_noScrn = pop(fromNonVSusL_noScrn) .* 0;    % when vaccines are depleted
        vaxdGroupSus_scrn = pop(fromNonVSusL_scrn) .* 0;
        vaxdGroupImm_noScrn = pop(fromNonVImmL_noScrn) .* 0;
        vaxdGroupImm_scrn = pop(fromNonVImmL_scrn) .* 0;
    end
    dPop(fromNonVSusL_noScrn) = dPop(fromNonVSusL_noScrn) - vaxdGroupSus_noScrn;
    dPop(fromNonVSusL_scrn) = dPop(fromNonVSusL_scrn) - vaxdGroupSus_scrn;
    dPop(fromNonVImmL_noScrn) = dPop(fromNonVImmL_noScrn) - vaxdGroupImm_noScrn;
    dPop(fromNonVImmL_scrn) = dPop(fromNonVImmL_scrn) - vaxdGroupImm_scrn;
    dPop(toVSusL_noScrn) = dPop(toVSusL_noScrn) + vaxdGroupSus_noScrn;
    dPop(toVSusL_scrn) = dPop(toVSusL_scrn) + vaxdGroupSus_scrn;
    dPop(toVImmL_noScrn) = dPop(toVImmL_noScrn) + vaxdGroupImm_noScrn;
    dPop(toVImmL_scrn) = dPop(toVImmL_scrn) + vaxdGroupImm_scrn;
    hpvVaxd = hpvVaxd + sumall(vaxdGroupSus_noScrn) + ...
        sumall(vaxdGroupSus_scrn) + sumall(vaxdGroupImm_noScrn) + ...
        sumall(vaxdGroupImm_scrn); % count number of people vaccinated at current time step
end

%% Convert dPop to a column vector for output to ODE solver
dPop = sparse(dPop);
