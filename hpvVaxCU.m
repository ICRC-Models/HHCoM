% HPV catch-up vaccination
function[dPop , hpvVaxd] = hpvVaxCU(pop , viral , risk , ...
    hpvVaxStates , hpvNonVaxStates , endpoints , intervens , vaxAgeCU , ...
    vaxCoverCU , vaxGCU , vaxDiseaseIndsCU , toInd , sumall)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
hpvVaxd = 0;

% Apply catch-up vaccination regimen
for dS = 1 : length(vaxDiseaseIndsCU)
    d = vaxDiseaseIndsCU(dS);
    for v = 1 : viral
        for g = min(vaxGCU) : max(vaxGCU) 
            for r = 1 : risk
                for aV = 1:length(vaxAgeCU) % apply appropriate coverage rate for each age group catch-up vaccinated    
                    a = vaxAgeCU(aV);
                    fromNonVSusCU_noScrn = toInd(allcomb(d , v , 1 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ...
                        g , a , r)); 
                    fromNonVSusCU_scrn = toInd(allcomb(d , v , 1 , 1 : hpvNonVaxStates , 1 : 3 , 3 , ...
                        g , a , r)); 
                    fromNonVImmCU_noScrn = toInd(allcomb(d , v , 7 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ...
                        g , a , r)); 
                    fromNonVImmCU_scrn = toInd(allcomb(d , v , 7 , 1 : hpvNonVaxStates , 1 : 3 , 3 , ...
                        g , a , r)); 
                    toVSusCU_noScrn = toInd(allcomb(d , v , 1 , 1 : hpvNonVaxStates , 1 : 3 , 2 , ...
                       g , a , r));
                    toVSusCU_scrn = toInd(allcomb(d , v , 1 , 1 : hpvNonVaxStates , 1 : 3 , 4 , ...
                       g , a , r));
                    toVImmCU_noScrn = toInd(allcomb(d , v , 7 , 1 : hpvNonVaxStates , 1 : 3 , 2 , ...
                       g , a , r));
                    toVImmCU_scrn = toInd(allcomb(d , v , 7 , 1 : hpvNonVaxStates , 1 : 3 , 4 , ...
                       g , a , r));
                    otherVCU = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                        1 : endpoints , [2,4] , g , a , r));
                    allVNonVCU = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                        1 : endpoints , 1 : intervens , g , a , r));
                
                    fracVaxd = sumall(pop(otherVCU)) / ...
                        sumall(pop(allVNonVCU)); % find proportion of population that is currently vaccinated
                    if vaxCoverCU(aV) - fracVaxd > 10 ^ -6 % when proportion vaccinated is below target vaccination level
                        vaxCover = max(0 , (vaxCoverCU(aV) - fracVaxd) ./ (1 - fracVaxd)); % vaccinate enough people in age group to reach target
                        vaxdGroupSus_noScrn = vaxCover .* pop(fromNonVSusCU_noScrn);
                        vaxdGroupSus_scrn = vaxCover .* pop(fromNonVSusCU_scrn);
                        vaxdGroupImm_noScrn = vaxCover .* pop(fromNonVImmCU_noScrn);
                        vaxdGroupImm_scrn = vaxCover .* pop(fromNonVImmCU_scrn);
                        dPop(fromNonVSusCU_noScrn) = dPop(fromNonVSusCU_noScrn) - vaxdGroupSus_noScrn;
                        dPop(fromNonVSusCU_scrn) = dPop(fromNonVSusCU_scrn) - vaxdGroupSus_scrn;
                        dPop(fromNonVImmCU_noScrn) = dPop(fromNonVImmCU_noScrn) - vaxdGroupImm_noScrn; 
                        dPop(fromNonVImmCU_scrn) = dPop(fromNonVImmCU_scrn) - vaxdGroupImm_scrn;
                        dPop(toVSusCU_noScrn) = dPop(toVSusCU_noScrn) + vaxdGroupSus_noScrn;
                        dPop(toVSusCU_scrn) = dPop(toVSusCU_scrn) + vaxdGroupSus_scrn;
                        dPop(toVImmCU_noScrn) = dPop(toVImmCU_noScrn) + vaxdGroupImm_noScrn;
                        dPop(toVImmCU_scrn) = dPop(toVImmCU_scrn) + vaxdGroupImm_scrn;
                        hpvVaxd = hpvVaxd + sumall(vaxdGroupSus_noScrn) + ...
                            sumall(vaxdGroupSus_scrn) + sumall(vaxdGroupImm_noScrn) + ...
                            sumall(vaxdGroupImm_scrn); % count number of people vaccinated at current time step
                    end

                end
            end
        end
    end
end

%% Convert dPop to a column vector for output to ODE solver
dPop = sparse(dPop);
