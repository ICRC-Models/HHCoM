% HPV catch-up vaccination
function[dPop , hpvVaxd, hpvVaxd1, hpvVaxd2] = hpvVaxCU(pop , viral , risk , ...
    hpvVaxStates , hpvNonVaxStates , endpoints , intervens , vaxAgeCU , ...
    vaxCoverCU , vaxGCU , vaxDiseaseIndsCU , toInd)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
hpvVaxd = 0;
% *** Uncomment when want CU vax with different vax eff for HIV -/+; This
% sets up tracking number of vaccines given to each group
%hpvVaxd1 = 0; %HIV neg vaccination 
%hpvVaxd2 = 0; % HIV positive


%% Apply catch-up vaccination regimen
for dS = 1 : length(vaxDiseaseIndsCU)
    d = vaxDiseaseIndsCU(dS);

%Classify disease state * UNCOMMENT WHEN RUNNING DIFFERENT VAX EFF FOR HIV -/+ 
%TRACKED FOR CU
%isHIVneg = (d == 1 || d == 2); 
%isHIVpos = (d >= 3 && d <= 8);

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
                           sumall(vaxdGroupImm_scrn); % count number of people vaccinated at current time step **Comment out when tracking hiv-/+ seperate vacciantion


% Uncomment when tracking hiv -/+ vaccination, groups vaccines given to
% those + or - when cu is on
%numVaxd = sumall(vaxdGroupSus_noScrn) + ...
         % sumall(vaxdGroupSus_scrn) + ...
         % sumall(vaxdGroupImm_noScrn) + ...
          %sumall(vaxdGroupImm_scrn);

% Update overall vaccination count
%hpvVaxd = hpvVaxd + numVaxd;

 %Track totals by HIV status
%if isHIVneg
   % hpvVaxd1 = hpvVaxd1 + numVaxd;   % disease states 1–2
%elseif isHIVpos
  %hpvVaxd2 = hpvVaxd2 + numVaxd;   % disease states 3–8
%end


                    end

                end
            end
        end
    end
end

%% Convert dPop to a column vector for output to ODE solver
dPop = sparse(dPop);
