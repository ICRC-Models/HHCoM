% HPV screening and treatment
function[dPop , ccScreen , ccTreatImm , ccTreatHpv , ccTreatHyst] = hpvScreen(pop , ...
    disease , viral , hpvTypes , hpvStates , risk , screenYrs , screenAlgs , ...
    year , stepsPerYear , screenAgeAll , screenAgeS , noVaxNoScreen , ...
    noVaxToScreen , vaxNoScreen , vaxToScreen , noVaxToScreenTreatImm , ...
    vaxToScreenTreatImm , noVaxToScreenTreatHpv , vaxToScreenTreatHpv , ...
    noVaxToScreenHyst , vaxToScreenHyst , screenAlgorithm , numScreenAge)

%% Set constants and initialize vectors
sumall = @(x) sum(x(:));
ccScreen = zeros(disease , viral , hpvTypes , hpvStates , numScreenAge , risk , 2);
ccTreatImm = ccScreen;
ccTreatHpv = ccScreen;
ccTreatHyst = ccScreen;

dPop = zeros(size(pop));

for i = 1 : length(screenAlgs)
    prevAL = 0;
    if i == 2
        prevAL = length(screenAlgs{1}.screenAge);
    end
    % Screening level
    dataYr1 = screenYrs(1);
    dataYrLast = screenYrs(size(screenYrs , 1));
    baseYrInd = max(find(year >= screenYrs , 1, 'last') , 1); % get index of first year <= current year
    baseYr = screenYrs(baseYrInd);
    screenRate = screenAlgs{i}.screenCover_vec{1}(1); % screening coverage up to 1st year
    if year < dataYrLast && year > dataYr1 % screening coverage between 1st and last year
        screenRate = screenAlgs{i}.screenCover_vec{baseYrInd}(round((year - baseYr) * stepsPerYear) + 1);
    elseif year >= dataYrLast % screening coverage last year and after
        lastInd = size(screenAlgs{i}.screenCover_vec , 1);
        screenRate = screenAlgs{i}.screenCover_vec{lastInd}(size(screenAlgs{i}.screenCover_vec{lastInd} , 2));
    end

    for aS = (prevAL + 1) : (prevAL + length(screenAlgs{i}.screenAge))
        for dS = 1 : length(screenAlgs{i}.diseaseInds)
            d = screenAlgs{i}.diseaseInds(dS);
            for v = 1 : viral
                for h = 1 : 2
                    for s = 1 : hpvStates
                        for r = 1 : risk
                            fracScreend = (sumall(pop(screenAgeS(d,v,h,s,:,aS,r))) / sumall(pop(screenAgeAll(d,v,h,s,:,aS,r)))); % find proportion of population that is currently screened
                            if screenRate - fracScreend > 10 ^ -6 % when proportion screened is below target screening level
                                screenCover = max(0 , (screenRate - fracScreend) ./ (1 - fracScreend)); % screen enough people in each compartment to reach target

                                % Baseline screening or CISNET or WHO screening algorithm
                                if any(s == [1 : 2 , 8 : 10])
                                    toScreenMult = 1.0;
                                    toScreenTreatImmMult = 0.0;
                                    toScreenTreatHpvMult = 0.0;
                                    toScreenTreatHystMult = 0.0;
                                elseif any(s == [3 : 4]) 
                                    toScreenMult = ((1-screenAlgs{i}.testSens(s)) + (screenAlgs{i}.testSens(s) * (1 - screenAlgs{i}.colpoRetain)) + ...
                                        (screenAlgs{i}.testSens(s) * screenAlgs{i}.colpoRetain * (1 - screenAlgs{i}.cinTreatRetain)) + ...
                                        (screenAlgs{i}.testSens(s) * screenAlgs{i}.colpoRetain * screenAlgs{i}.cinTreatRetain * (1-screenAlgs{i}.cinTreatEff(d))));
                                    toScreenTreatImmMult = screenAlgs{i}.testSens(s) * screenAlgs{i}.colpoRetain * screenAlgs{i}.cinTreatRetain * screenAlgs{i}.cinTreatEff(d) * ...
                                        (1.0-((screenAlgs{i}.cinTreatHpvPersist - (1-screenAlgs{i}.cinTreatEff(d)))/screenAlgs{i}.cinTreatEff(d)));
                                    toScreenTreatHpvMult = screenAlgs{i}.testSens(s) * screenAlgs{i}.colpoRetain * screenAlgs{i}.cinTreatRetain * screenAlgs{i}.cinTreatEff(d) * ...
                                        ((screenAlgs{i}.cinTreatHpvPersist - (1-screenAlgs{i}.cinTreatEff(d)))/screenAlgs{i}.cinTreatEff(d));
                                    toScreenTreatHystMult = 0.0;
                                elseif any(s == [5 : 7]) 
                                    toScreenMult = ((1-screenAlgs{i}.testSens(s)) + (screenAlgs{i}.testSens(s) * (1 - screenAlgs{i}.colpoRetain)) + ...
                                        (screenAlgs{i}.testSens(s) * screenAlgs{i}.colpoRetain * (1 - screenAlgs{i}.ccTreatRetain)));
                                    toScreenTreatImmMult = 0.0;
                                    toScreenTreatHpvMult = 0.0;
                                    toScreenTreatHystMult = screenAlgs{i}.testSens(s) * screenAlgs{i}.colpoRetain * screenAlgs{i}.ccTreatRetain;
                                end

                                noVaxScreend = screenCover .* pop(noVaxNoScreen(d,v,h,s,aS,r));
                                dPop(noVaxNoScreen(d,v,h,s,aS,r)) = dPop(noVaxNoScreen(d,v,h,s,aS,r)) - noVaxScreend;
                                dPop(noVaxToScreen(d,v,h,s,aS,r)) = dPop(noVaxToScreen(d,v,h,s,aS,r)) + toScreenMult .* noVaxScreend;
                                dPop(noVaxToScreenTreatImm(d,v,aS,r)) = dPop(noVaxToScreenTreatImm(d,v,aS,r)) + toScreenTreatImmMult .* noVaxScreend;
                                dPop(noVaxToScreenTreatHpv(d,v,aS,r)) = dPop(noVaxToScreenTreatHpv(d,v,aS,r)) + toScreenTreatHpvMult .* noVaxScreend;
                                dPop(noVaxToScreenHyst(d,v,aS,r)) = dPop(noVaxToScreenHyst(d,v,aS,r)) + toScreenTreatHystMult .* noVaxScreend;
                                ccScreen(d , v , h , s , aS , r , 1) = sumall(noVaxScreend);
                                ccTreatImm(d , v , h , s , aS , r , 1) = sumall(toScreenTreatImmMult .* noVaxScreend);
                                ccTreatHpv(d , v , h , s , aS , r , 1) = sumall(toScreenTreatHpvMult .* noVaxScreend);
                                ccTreatHyst(d , v , h , s , aS , r , 1) = sumall(toScreenTreatHystMult .* noVaxScreend);

                                vaxScreend = screenCover .* pop(vaxNoScreen(d,v,h,s,aS,r));
                                dPop(vaxNoScreen(d,v,h,s,aS,r)) = dPop(vaxNoScreen(d,v,h,s,aS,r)) - vaxScreend;
                                dPop(vaxToScreen(d,v,h,s,aS,r)) = dPop(vaxToScreen(d,v,h,s,aS,r)) + toScreenMult .* vaxScreend;
                                dPop(vaxToScreenTreatImm(d,v,aS,r)) = dPop(vaxToScreenTreatImm(d,v,aS,r)) + toScreenTreatImmMult .* vaxScreend;
                                dPop(vaxToScreenTreatHpv(d,v,aS,r)) = dPop(vaxToScreenTreatHpv(d,v,aS,r)) + toScreenTreatHpvMult .* vaxScreend;
                                dPop(vaxToScreenHyst(d,v,aS,r)) = dPop(vaxToScreenHyst(d,v,aS,r)) + toScreenTreatHystMult .* vaxScreend;
                                ccScreen(d , v , h , s , aS , r , 2) = sumall(vaxScreend);
                                ccTreatImm(d , v , h , s , aS , r , 2) = sumall(toScreenTreatImmMult .* vaxScreend);
                                ccTreatHpv(d , v , h , s , aS , r , 2) = sumall(toScreenTreatHpvMult .* vaxScreend);
                                ccTreatHyst(d , v , h , s , aS , r , 2) = sumall(toScreenTreatHystMult .* vaxScreend);
                            end
                        end
                    end
                end
            end
        end
    end
end
dPop = sparse(dPop);
