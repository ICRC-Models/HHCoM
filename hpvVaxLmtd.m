% HPV vaccination in a limited vaccine scenario
function[dPop , hpvVaxd , vaxRemain] = hpvVaxLmtd(pop , k , year , vaxLimitPerYr , ...
    disease , viral , risk , hpvTypes , hpvStates , periods , vaxCoverL , ...
    vaxRemain , vaxGL)

%% Set constants and initialize vectors
toInd = @(x) (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + (x(: , 6) - 1) * k(5) ...
    + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + (x(: , 3) - 1) * k(2) ...
    + (x(: , 2) - 1) * k(1) + x(: , 1);
sumall = @(x) sum(x(:));
hpvVaxd = 0;

dPop = zeros(size(pop));

% If within first vaxLimitYrs-many vaccine-limited years
if rem(year,1) == 0.0    % reset vaxRemain at start of each new year to the number of available vaccines per year
    vaxRemain = vaxLimitPerYr;
end

fromNonVSusL_noScrn = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 , ...
    vaxGL , vaxAgeL , 1 : risk));
fromNonVSusL_scrn = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 6 , ...
    vaxGL , vaxAgeL , 1 : risk));
fromNonVImmL_noScrn = toInd(allcomb(1 : disease , 1 : viral , 2 , 10 , 1 , ...
    vaxGL , vaxAgeL , 1 : risk));
fromNonVImmL_scrn = toInd(allcomb(1 : disease , 1 : viral , 2 , 10 , 6 , ...
    vaxGL , vaxAgeL , 1 : risk));
fromNonVImmLNonV_noScrn = toInd(allcomb(1 : disease , 1 : viral , 3 , 10 , 1 , ...
    vaxGL , vaxAgeL , 1 : risk));
fromNonVImmLNonV_scrn = toInd(allcomb(1 : disease , 1 : viral , 3 , 10 , 6 , ...
    vaxGL , vaxAgeL , 1 : risk));  
toVL_noScrn = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 1 , ...
    vaxGL , vaxAgeL , 1 : risk));
toVL_scrn = toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , 6 , ...
    vaxGL , vaxAgeL , 1 : risk));
otherVL = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    [2,4] , vaxGL , vaxAgeL , 1 : risk));
allVNonVL = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , vaxGL , vaxAgeL , 1 : risk));

fracVaxd = (sumall(pop(toVL_noScrn)) + sumall(pop(toVL_scrn)) + ...
    sumall(pop(otherVL))) / (sumall(pop(allVNonVL))); % find proportion of population that is currently vaccinated
if vaxCoverL - fracVaxd > 10 ^ -6 % when proportion vaccinated is below target vaccination level
    vaxCover = max(0 , (vaxCoverL - fracVaxd) ./ (1 - fracVaxd)); % vaccinate enough people in age group to reach target
    vaxdGroupSum = (sumall(pop(fromNonVSusL_noScrn)) + ...
        sumall(pop(fromNonVSusL_scrn)) + sumall(pop(fromNonVImmL_noScrn)) + ...
        sumall(pop(fromNonVImmL_scrn)) + sumall(pop(fromNonVImmLNonV_noScrn)) + sumall(pop(fromNonVImmLNonV_scrn)));
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
            vaxdGroupImmNonV_noScrn = (usedVax .* (pop(fromNonVImmLNonV_noScrn) ./ vaxdGroupSum))...
                .* (pop(fromNonVImmLNonV_noScrn) > 0);
            vaxdGroupImmNonV_scrn = (usedVax .* (pop(fromNonVImmLNonV_scrn) ./ vaxdGroupSum))...
                .* (pop(fromNonVImmLNonV_scrn) > 0);
        end
        vaxRemain = vaxRemain - usedVax;
    else 
        vaxdGroupSus_noScrn = pop(fromNonVSusL_noScrn) .* 0;    % when vaccines are depleted
        vaxdGroupSus_scrn = pop(fromNonVSusL_scrn) .* 0;
        vaxdGroupImm_noScrn = pop(fromNonVImmL_noScrn) .* 0;
        vaxdGroupImm_scrn = pop(fromNonVImmL_scrn) .* 0;
        vaxdGroupImmNonV_noScrn = pop(fromNonVImmLNonV_noScrn) .* 0;
        vaxdGroupImmNonV_scrn = pop(fromNonVImmLNonV_scrn) .* 0;
    end
    dPop(fromNonVSusL_noScrn) = dPop(fromNonVSusL_noScrn) - vaxdGroupSus_noScrn;
    dPop(fromNonVSusL_scrn) = dPop(fromNonVSusL_scrn) - vaxdGroupSus_scrn;
    dPop(fromNonVImmL_noScrn) = dPop(fromNonVImmL_noScrn) - vaxdGroupImm_noScrn;
    dPop(fromNonVImmL_scrn) = dPop(fromNonVImmL_scrn) - vaxdGroupImm_scrn;
    dPop(fromNonVImmLNonV_noScrn) = dPop(fromNonVImmLNonV_noScrn) - vaxdGroupImmNonV_noScrn;
    dPop(fromNonVImmLNonV_scrn) = dPop(fromNonVImmLNonV_scrn) - vaxdGroupImmNonV_scrn;
    dPop(toVL_noScrn) = dPop(toVL_noScrn) + vaxdGroupSus_noScrn + vaxdGroupImm_noScrn + vaxdGroupImmNonV_noScrn;
    dPop(toVL_scrn) = dPop(toVL_scrn) + vaxdGroupSus_scrn + vaxdGroupImm_scrn + vaxdGroupImmNonV_scrn;
    hpvVaxd = hpvVaxd + sumall(vaxdGroupSus_noScrn) + ...
        sumall(vaxdGroupSus_scrn) + sumall(vaxdGroupImm_noScrn) + ...
        sumall(vaxdGroupImm_scrn) + sumall(vaxdGroupImmNonV_noScrn) + sumall(vaxdGroupImmNonV_scrn); % count number of people vaccinated at current time step
end
                        
dPop = sparse(dPop);
