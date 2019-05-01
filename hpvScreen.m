% HPV screening and treatment
function[dPop , ccScreen , ccTreatImm , ccTreatHpv , ccTreatHyst] = hpvScreen(pop , disease , viral , hpvTypes , ...
    hpvStates , risk , screenYrs , screenCover_vec , testSens , screenAge , ...
    year , stepsPerYear , screenAgeAll , screenAgeS , noVaxNoScreen , ...
    noVaxToScreen , vaxNoScreen , vaxToScreen , noVaxToScreenTreatImm , ...
    vaxToScreenTreatImm , noVaxToScreenTreatHpv , vaxToScreenTreatHpv , ...
    noVaxToScreenHyst , vaxToScreenHyst , colpoRetain , cinTreatRetain , ...
    cinTreatEff , cinTreatHpvPersist , ccTreatRetain , nhScreenParams)

%% Set constants and initialize vectors
sumall = @(x) sum(x(:));
ccScreen = zeros(disease , viral , hpvTypes , hpvStates , risk , 2);
ccTreatImm = ccScreen;
ccTreatHpv = ccScreen;
ccTreatHyst = ccScreen;

dPop = zeros(size(pop));

% Screening level
dataYr1 = screenYrs(1);
dataYrLast = screenYrs(size(screenYrs , 1));
baseYrInd = max(find(year >= screenYrs , 1, 'last') , 1); % get index of first year <= current year
baseYr = screenYrs(baseYrInd);
screenRate = screenCover_vec{1}(1); % screening coverage up to 1st year
if year < dataYrLast && year > dataYr1 % screening coverage between 1st and last year
    screenRate = screenCover_vec{baseYrInd}(round((year - baseYr) * stepsPerYear) + 1);
elseif year >= dataYrLast % screening coverage last year and after
    lastInd = size(screenCover_vec , 1);
    screenRate = screenCover_vec{lastInd}(size(screenCover_vec{lastInd} , 2));
end
screenRate = screenRate * 0.20; % find 1/5 of age group (represents 35 year olds, only)

for aS = 1 : length(screenAge)
    for d = 1 : disease
        for v = 1 : viral
            for h = 1 : 2
                for s = 1 : hpvStates
                    for r = 1 : risk
                        fracScreend = (sumall(pop(screenAgeS(d,v,h,s,:,aS,r))) / sumall(pop(screenAgeAll(d,v,h,s,:,aS,r)))); % find proportion of population that is currently screened
                        if screenRate - fracScreend > 10 ^ -6 % when proportion screened is below target screening level
                            screenCover = max(0 , (screenRate - fracScreend) ./ (1 - fracScreend)); % screen enough people in each compartment to reach target
                            
                            if nhScreenParams
                                % NH screening algorithm
                                if any(s == [1 : 2 , 8 : 10])
                                    toScreenMult = 1.0;
                                    toScreenTreatImmMult = 0.0;
                                    toScreenTreatHpvMult = 0.0;
                                    toScreenTreatHystMult = 0.0;
                                elseif any(s == [3 : 4]) 
                                    toScreenMult = ((1-testSens) + (testSens * (1 - colpoRetain)) + ...
                                        (testSens * colpoRetain * (1 - cinTreatRetain)) + ...
                                        (testSens * colpoRetain * cinTreatRetain * (1-cinTreatEff(d))));
                                    toScreenTreatImmMult = testSens * colpoRetain * cinTreatRetain * cinTreatEff(d) * ...
                                        (1.0-((cinTreatHpvPersist - (1-cinTreatEff(d)))/cinTreatEff(d)));
                                    toScreenTreatHpvMult = testSens * colpoRetain * cinTreatRetain * cinTreatEff(d) * ...
                                        ((cinTreatHpvPersist - (1-cinTreatEff(d)))/cinTreatEff(d));
                                    toScreenTreatHystMult = 0.0;
                                elseif any(s == [5 : 7]) 
                                    toScreenMult = ((1-testSens) + (testSens * (1 - colpoRetain)) + ...
                                        (testSens * colpoRetain * (1 - ccTreatRetain)));
                                    toScreenTreatImmMult = 0.0;
                                    toScreenTreatHpvMult = 0.0;
                                    toScreenTreatHystMult = testSens * colpoRetain * ccTreatRetain;
                                %elseif any(s == [6 : 7])
                                %    toScreenMult = 1.0;
                                %    toScreenTreatImmMult = 0.0;
                                %    toScreenTreatHpvMult = 0.0;
                                %    toScreenTreatHystMult = 0.0;
                                end
                            else
                                % WHO screening algorithm
                                if any(s == [5 : 10]) || (s==1 && h==1)
                                    toScreenMult = 1.0;
                                    toScreenTreatImmMult = 0.0;
                                    toScreenTreatHpvMult = 0.0;
                                    toScreenTreatHystMult = 0.0;
                                elseif any(s == [2 : 4]) 
                                    toScreenMult = ((1-testSens) + (testSens * (1 - colpoRetain)) + ...
                                        (testSens * colpoRetain * (1 - cinTreatRetain)) + ...
                                        (testSens * colpoRetain * cinTreatRetain * (1-cinTreatEff(d))));
                                    toScreenTreatImmMult = testSens * colpoRetain * cinTreatRetain * cinTreatEff(d) * ...
                                        (1.0-((cinTreatHpvPersist - (1-cinTreatEff(d)))/cinTreatEff(d)));
                                    toScreenTreatHpvMult = testSens * colpoRetain * cinTreatRetain * cinTreatEff(d) * ...
                                        ((cinTreatHpvPersist - (1-cinTreatEff(d)))/cinTreatEff(d));
                                    toScreenTreatHystMult = 0.0;
                                elseif (s == 1 && h == 2)
                                    toScreenMult = ((1-testSens) + (testSens * (1 - colpoRetain)) + ...
                                        (testSens * colpoRetain * (1 - cinTreatRetain)));
                                    toScreenTreatImmMult = testSens * colpoRetain * cinTreatRetain * (1.0-cinTreatHpvPersist);
                                    toScreenTreatHpvMult = testSens * colpoRetain * cinTreatRetain * cinTreatHpvPersist;
                                    toScreenTreatHystMult = 0.0;
                                end
                            end
                            
                            noVaxScreend = screenCover .* pop(noVaxNoScreen(d,v,h,s,aS,r));
                            dPop(noVaxNoScreen(d,v,h,s,aS,r)) = dPop(noVaxNoScreen(d,v,h,s,aS,r)) - noVaxScreend;
                            dPop(noVaxToScreen(d,v,h,s,aS,r)) = dPop(noVaxToScreen(d,v,h,s,aS,r)) + toScreenMult .* noVaxScreend;
                            dPop(noVaxToScreenTreatImm(d,v,aS,r)) = dPop(noVaxToScreenTreatImm(d,v,aS,r)) + toScreenTreatImmMult .* noVaxScreend;
                            dPop(noVaxToScreenTreatHpv(d,v,aS,r)) = dPop(noVaxToScreenTreatHpv(d,v,aS,r)) + toScreenTreatHpvMult .* noVaxScreend;
                            dPop(noVaxToScreenHyst(d,v,aS,r)) = dPop(noVaxToScreenHyst(d,v,aS,r)) + toScreenTreatHystMult .* noVaxScreend;
                            ccScreen(d , v , h , s , r , 1) = sumall(noVaxScreend);
                            ccTreatImm(d , v , h , s , r , 1) = sumall(toScreenTreatImmMult .* noVaxScreend);
                            ccTreatHpv(d , v , h , s , r , 1) = sumall(toScreenTreatHpvMult .* noVaxScreend);
                            ccTreatHyst(d , v , h , s , r , 1) = sumall(toScreenTreatHystMult .* noVaxScreend);
                            
                            vaxScreend = screenCover .* pop(vaxNoScreen(d,v,h,s,aS,r));
                            dPop(vaxNoScreen(d,v,h,s,aS,r)) = dPop(vaxNoScreen(d,v,h,s,aS,r)) - vaxScreend;
                            dPop(vaxToScreen(d,v,h,s,aS,r)) = dPop(vaxToScreen(d,v,h,s,aS,r)) + toScreenMult .* vaxScreend;
                            dPop(vaxToScreenTreatImm(d,v,aS,r)) = dPop(vaxToScreenTreatImm(d,v,aS,r)) + toScreenTreatImmMult .* vaxScreend;
                            dPop(vaxToScreenTreatHpv(d,v,aS,r)) = dPop(vaxToScreenTreatHpv(d,v,aS,r)) + toScreenTreatHpvMult .* vaxScreend;
                            dPop(vaxToScreenHyst(d,v,aS,r)) = dPop(vaxToScreenHyst(d,v,aS,r)) + toScreenTreatHystMult .* vaxScreend;
                            ccScreen(d , v , h , s , r , 2) = sumall(vaxScreend);
                            ccTreatImm(d , v , h , s , r , 1) = sumall(toScreenTreatImmMult .* noVaxScreend);
                            ccTreatHpv(d , v , h , s , r , 1) = sumall(toScreenTreatHpvMult .* noVaxScreend);
                            ccTreatHyst(d , v , h , s , r , 1) = sumall(toScreenTreatHystMult .* noVaxScreend);
                        end
                    end
                end
            end
        end
    end
end

dPop = sparse(dPop);
