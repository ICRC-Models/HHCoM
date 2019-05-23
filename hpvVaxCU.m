% HPV vaccination
function[dPop , hpvVaxd] = hpvVaxCU(pop , k , disease , viral , risk , ...
    hpvTypes , hpvStates , periods , vaxAgeCU , vaxCoverCU , vaxGCU)

%% Set constants and initialize vectors
toInd = @(x) (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + (x(: , 6) - 1) * k(5) ...
    + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + (x(: , 3) - 1) * k(2) ...
    + (x(: , 2) - 1) * k(1) + x(: , 1);
sumall = @(x) sum(x(:));
hpvVaxd = 0;

dPop = zeros(size(pop));

% Apply catch-up vaccination regimen
for d = 1 : disease
    for v = 1 : viral
        for g = min(vaxGCU) : max(vaxGCU) 
            for r = 1 : risk
                for aV = 1:length(vaxAgeCU) % apply appropriate coverage rate for each age group catch-up vaccinated    
                    a = vaxAgeCU(aV);
                    fromNonVSusCU_noScrn(:,aV) = toInd(allcomb(d , v , 1 , 1 , 1 , ...
                        g , a , r)); 
                    fromNonVSusCU_scrn(:,aV) = toInd(allcomb(d , v , 1 , 1 , 6 , ...
                        g , a , r)); 
                    fromNonVImmCU_noScrn(:,aV) = toInd(allcomb(d , v , 2 , 10 , 1 , ...
                        g , a , r)); 
                    fromNonVImmCU_scrn(:,aV) = toInd(allcomb(d , v , 2 , 10 , 6 , ...
                        g , a , r)); 
                    toVCU_noScrn(:,aV) = toInd(allcomb(d , v , 1 , 9 , 1 , ...
                       g , a , r));
                    toVCU_scrn(:,aV) = toInd(allcomb(d , v , 1 , 9 , 6 , ...
                       g , a , r));
                    otherVCU(:,aV) = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , ...
                        [2,4] , g , a , r));
                    allVNonVCU(:,aV) = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , ...
                        1 : periods , g , a , r));
                
                    fracVaxd = (sumall(pop(toVCU_noScrn(:,aV))) + sumall(pop(toVCU_scrn(:,aV))) ...
                        + sumall(pop(otherVCU(:,aV)))) / (sumall(pop(allVNonVCU(:,aV)))); % find proportion of population that is currently vaccinated
                    if vaxCoverCU(aV) - fracVaxd > 10 ^ -6 % when proportion vaccinated is below target vaccination level
                        vaxCover = max(0 , (vaxCoverCU(aV) - fracVaxd) ./ (1 - fracVaxd)); % vaccinate enough people in age group to reach target
                        vaxdGroupSus_noScrn = vaxCover .* pop(fromNonVSusCU_noScrn(:,aV));
                        vaxdGroupSus_scrn = vaxCover .* pop(fromNonVSusCU_scrn(:,aV));
                        vaxdGroupImm_noScrn = vaxCover .* pop(fromNonVImmCU_noScrn(:,aV));
                        vaxdGroupImm_scrn = vaxCover .* pop(fromNonVImmCU_scrn(:,aV));
                        dPop(fromNonVSusCU_noScrn(:,aV)) = dPop(fromNonVSusCU_noScrn(:,aV)) - vaxdGroupSus_noScrn;
                        dPop(fromNonVSusCU_scrn(:,aV)) = dPop(fromNonVSusCU_scrn(:,aV)) - vaxdGroupSus_scrn;
                        dPop(fromNonVImmCU_noScrn(:,aV)) = dPop(fromNonVImmCU_noScrn(:,aV)) - vaxdGroupImm_noScrn; 
                        dPop(fromNonVImmCU_scrn(:,aV)) = dPop(fromNonVImmCU_scrn(:,aV)) - vaxdGroupImm_scrn; 
                        dPop(toVCU_noScrn(:,aV)) = dPop(toVCU_noScrn(:,aV)) + vaxdGroupSus_noScrn + vaxdGroupImm_noScrn;
                        dPop(toVCU_scrn(:,aV)) = dPop(toVCU_scrn(:,aV)) + vaxdGroupSus_scrn + vaxdGroupImm_scrn;
                        hpvVaxd = hpvVaxd + sumall(vaxdGroupSus_noScrn) + ...
                            sumall(vaxdGroupSus_scrn) + sumall(vaxdGroupImm_noScrn) + ...
                            sumall(vaxdGroupImm_scrn); % count number of people vaccinated at current time step
                    end

                end
            end
        end
    end
end

dPop = sparse(dPop);
