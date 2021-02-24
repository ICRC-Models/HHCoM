% HPV screening and treatment
function[dPop , ccScreen] = hpvScreen(pop , ...
    disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , risk , ...
    screenYrs , screenAlgs , year , stepsPerYear , screenAgeAll , screenAgeS , ...
    noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , noVaxToScreenTreatImm , ...
    vaxToScreenTreatImm , noVaxToScreenTreatHpv , vaxToScreenTreatHpv , ...
    noVaxToScreenTreatVaxHpv , vaxToScreenTreatVaxHpv , noVaxToScreenTreatNonVaxHpv , ...
    vaxToScreenTreatNonVaxHpv , noVaxToScreenHyst , vaxToScreenHyst , numScreenAge , ageMultsComb)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
ccScreen = zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , numScreenAge , risk , 2);
% ccTreatImm = ccScreen;
% ccTreatHpv = ccScreen;
% ccTreatHyst = ccScreen;

%% Run screening algorithm
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
        screenRateAge = screenRate * ageMultsComb(aS); % apply screening to fraction of age group (1/5 represents 35 year olds only with 5-year age groups)
        for dS = 1 : length(screenAlgs{i}.diseaseInds)
            d = screenAlgs{i}.diseaseInds(dS);
            for v = 1 : viral
                for h = 1 : hpvVaxStates
                    for s = 1 : hpvNonVaxStates
                        for x = 1 : endpoints
                            for r = 1 : risk
                                fracScreend = (sumall(pop(screenAgeS(d,v,h,s,x,:,aS,r))) / sumall(pop(screenAgeAll(d,v,h,s,x,:,aS,r)))); % find proportion of population that is currently screened
                                if screenRateAge - fracScreend > 10 ^ -6 % when proportion screened is below target screening level
                                    screenCover = max(0 , (screenRateAge - fracScreend) ./ (1 - fracScreend)); % screen enough people in each compartment to reach target

                                    % Apply selected screening algorithm(s)
                                    % if you're susceptible/infected/CIN1/immune to both HPV types or have had a hysterectomy
                                    if ( ((h<=3) || (h==7)) && ((s<=3) || (s==7)) && (x==1) ) || (x==4)
                                        toScreenMult = 1.0;
                                        toScreenTreatImmMult = 0.0;
                                        toScreenTreatHpvMult = 0.0;
                                        toScreenTreatHystMult = 0.0;
                                    % if you have CIN2+ of either HPV type
                                    elseif ( (((h==4) || (h==5)) && ((s<=5) || (s==7))) || (((h<=5) || (h==7)) && ((s==4) || (s==5))) ) && (x==1) 
                                        toScreenMult = ((1-screenAlgs{i}.testSens(2)) + (screenAlgs{i}.testSens(2) * (1 - screenAlgs{i}.colpoRetain)) + ...
                                            (screenAlgs{i}.testSens(2) * screenAlgs{i}.colpoRetain * (1 - screenAlgs{i}.cinTreatRetain)) + ...
                                            (screenAlgs{i}.testSens(2) * screenAlgs{i}.colpoRetain * screenAlgs{i}.cinTreatRetain * (1-screenAlgs{i}.cinTreatEff(d))));
                                        toScreenTreatImmMult = screenAlgs{i}.testSens(2) * screenAlgs{i}.colpoRetain * screenAlgs{i}.cinTreatRetain * screenAlgs{i}.cinTreatEff(d) * ...
                                            (1.0-(screenAlgs{i}.cinTreatHpvPersistHivNeg/screenAlgs{i}.cinTreatEff(d)));
                                        toScreenTreatHpvMult = screenAlgs{i}.testSens(2) * screenAlgs{i}.colpoRetain * screenAlgs{i}.cinTreatRetain * screenAlgs{i}.cinTreatEff(d) * ...
                                            (screenAlgs{i}.cinTreatHpvPersistHivNeg/screenAlgs{i}.cinTreatEff(d));
                                        toScreenTreatHystMult = 0.0;
                                    % if you have cervical cancer
                                    elseif ( (x==1) && ((h==6) || (s==6)) ) || (x==2) || (x==3)
                                        toScreenMult = ((1-screenAlgs{i}.testSens(3)) + (screenAlgs{i}.testSens(3) * (1 - screenAlgs{i}.colpoRetain)) + ...
                                            (screenAlgs{i}.testSens(3) * screenAlgs{i}.colpoRetain * (1 - screenAlgs{i}.ccTreatRetain)));
                                        toScreenTreatImmMult = 0.0;
                                        toScreenTreatHpvMult = 0.0;
                                        toScreenTreatHystMult = screenAlgs{i}.testSens(3) * screenAlgs{i}.colpoRetain * screenAlgs{i}.ccTreatRetain;
                                    end
                                    
                                    noVaxScreend = screenCover .* pop(noVaxNoScreen(d,v,h,s,x,aS,r));
                                    dPop(noVaxNoScreen(d,v,h,s,x,aS,r)) = dPop(noVaxNoScreen(d,v,h,s,x,aS,r)) - noVaxScreend;
                                    dPop(noVaxToScreen(d,v,h,s,x,aS,r)) = dPop(noVaxToScreen(d,v,h,s,x,aS,r)) + toScreenMult .* noVaxScreend;
                                    dPop(noVaxToScreenTreatImm(d,v,aS,r)) = dPop(noVaxToScreenTreatImm(d,v,aS,r)) + toScreenTreatImmMult .* noVaxScreend;
                                    dPop(noVaxToScreenHyst(d,v,aS,r)) = dPop(noVaxToScreenHyst(d,v,aS,r)) + toScreenTreatHystMult .* noVaxScreend;
                                    
                                    vaxScreend = screenCover .* pop(vaxNoScreen(d,v,h,s,x,aS,r));
                                    dPop(vaxNoScreen(d,v,h,s,x,aS,r)) = dPop(vaxNoScreen(d,v,h,s,x,aS,r)) - vaxScreend;
                                    dPop(vaxToScreen(d,v,h,s,x,aS,r)) = dPop(vaxToScreen(d,v,h,s,x,aS,r)) + toScreenMult .* vaxScreend;
                                    dPop(vaxToScreenTreatImm(d,v,aS,r)) = dPop(vaxToScreenTreatImm(d,v,aS,r)) + toScreenTreatImmMult .* vaxScreend;
                                    dPop(vaxToScreenHyst(d,v,aS,r)) = dPop(vaxToScreenHyst(d,v,aS,r)) + toScreenTreatHystMult .* vaxScreend;
                                            
                                    % if you have CIN2+ of either HPV type and are susceptible or immune to the other HPV type
                                    % Note: don't want to move individuals susceptible/immune to the other HPV type falsely into an HPV infected compartment
                                    if ( (((h==4) || (h==5)) && (((s==1) || (s==7)))) || ((((h==1) || (h==7))) && ((s==4) || (s==5))) ) && (x==1)
                                        if ((h==4) || (h==5))
                                            dPop(noVaxToScreenTreatVaxHpv(d,v,s,aS,r)) = dPop(noVaxToScreenTreatVaxHpv(d,v,s,aS,r)) + toScreenTreatHpvMult .* noVaxScreend;
                                            dPop(vaxToScreenTreatVaxHpv(d,v,s,aS,r)) = dPop(vaxToScreenTreatVaxHpv(d,v,s,aS,r)) + toScreenTreatHpvMult .* vaxScreend;
                                        elseif ((s==4) || (s==5))    
                                            dPop(noVaxToScreenTreatNonVaxHpv(d,v,h,aS,r)) = dPop(noVaxToScreenTreatNonVaxHpv(d,v,h,aS,r)) + toScreenTreatHpvMult .* noVaxScreend;
                                            dPop(vaxToScreenTreatNonVaxHpv(d,v,h,aS,r)) = dPop(vaxToScreenTreatNonVaxHpv(d,v,h,aS,r)) + toScreenTreatHpvMult .* vaxScreend;
                                        end
                                    % if you have any other combination of HPV types and states
                                    else     
                                        dPop(noVaxToScreenTreatHpv(d,v,aS,r)) = dPop(noVaxToScreenTreatHpv(d,v,aS,r)) + toScreenTreatHpvMult .* noVaxScreend;
                                        dPop(vaxToScreenTreatHpv(d,v,aS,r)) = dPop(vaxToScreenTreatHpv(d,v,aS,r)) + toScreenTreatHpvMult .* vaxScreend;
                                    end    
                                        
                                    ccScreen(d , v , h , s , x , aS , r , 1) = sumall(noVaxScreend);
                                    % ccTreatImm(d , v , h , s , x , aS , r , 1) = sumall(toScreenTreatImmMult .* noVaxScreend);
                                    % ccTreatHpv(d , v , h , s , x , aS , r , 1) = sumall(toScreenTreatHpvMult .* noVaxScreend);
                                    % ccTreatHyst(d , v , h , s , x , aS , r , 1) = sumall(toScreenTreatHystMult .* noVaxScreend);

                                    ccScreen(d , v , h , s , x , aS , r , 2) = sumall(vaxScreend);
                                    % ccTreatImm(d , v , h , s , x , aS , r , 2) = sumall(toScreenTreatImmMult .* vaxScreend);
                                    % ccTreatHpv(d , v , h , s , x , aS , r , 2) = sumall(toScreenTreatHpvMult .* vaxScreend);
                                    % ccTreatHyst(d , v , h , s , x , aS , r , 2) = sumall(toScreenTreatHystMult .* vaxScreend);
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

%% Convert dPop to a column vector for output to ODE solver
dPop = sparse(dPop);
