% HPV screening and treatment
function[dPop , extraOut] = hpvScreen(t , pop , disease , hpvTypes , hpvStates , ...
    screenYrs , screenCover_vec , testSens , screenAge , ...
    year , stepsPerYear , screenAgeAll , screenAgeS , ...
    noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , ...
    noVaxToScreenTreatImm , vaxToScreenTreatImm , noVaxToScreenTreatHpv , ...
    vaxToScreenTreatHpv , cryoEff , cryoElig , cryoHpvClear)

%% Set constants and initialize vectors
sumall = @(x) sum(x(:));
ccScreen = zeros(hpvTypes , hpvStates);

dPop = zeros(size(pop));

% Screening level
dataYr1 = screenYrs(1);
dataYrLast = screenYrs(size(screenYrs , 1));
baseYrInd = max(find(year >= screenYrs , 1, 'last') , 1); % get index of first year <= current year
baseYr = screenYrs(baseYrInd);
screenRate = screenCover_vec{1}(1); % screening coverage up to 1st year
if year < dataYrLast && year > dataYr1 % screening coverage between 1st and last year
    screenRate = screenCover_vec{baseYrInd}(round((year - baseYr) * stepsPerYear) + 1);
elseif year > dataYrLast % assortativity in last year
    lastInd = size(screenCover_vec , 1);
    screenRate = screenCover_vec{lastInd}(size(screenCover_vec{lastInd} , 2));
end
screenRate = screenRate * 0.20; % find 1/5 of age group (represents 35 year olds, only)

for aS = 1 : length(screenAge)
    fracScreend = (sumall(pop(screenAgeS(:,aS))) / sumall(pop(screenAgeAll(:,aS)))); % find proportion of population that is currently screened
    if screenRate - fracScreend > 10 ^ -6 % when proportion screened is below target screening level
        screenCover = max(0 , (screenRate - fracScreend) ./ (1 - fracScreend)); % screen enough people in each compartment to reach target

        for d = 1 : disease
            for h = 1 : 2
                for s = [1 , 8 : 10]
                    noVaxScreend = screenCover .* pop(noVaxNoScreen(d,h,s,:,aS));
                    dPop(noVaxNoScreen(d,h,s,:,aS)) = dPop(noVaxNoScreen(d,h,s,:,aS)) - noVaxScreend;
                    dPop(noVaxToScreen(d,h,s,:,aS)) = dPop(noVaxToScreen(d,h,s,:,aS)) + noVaxScreend;

                    vaxScreend = screenCover .* pop(vaxNoScreen(d,h,s,:,aS));
                    dPop(vaxNoScreen(d,h,s,:,aS)) = dPop(vaxNoScreen(d,h,s,:,aS)) - vaxScreend;
                    dPop(vaxToScreen(d,h,s,:,aS)) = dPop(vaxToScreen(d,h,s,:,aS)) + vaxScreend;
                end
                for s = 2 : 7
                    noVaxScreend = screenCover .* pop(noVaxNoScreen(d,h,s,:,aS));
                    dPop(noVaxNoScreen(d,h,s,:,aS)) = dPop(noVaxNoScreen(d,h,s,:,aS)) - noVaxScreend;
                    dPop(noVaxToScreen(d,h,s,:,aS)) = dPop(noVaxToScreen(d,h,s,:,aS)) + ...
                        ((1-testSens) + (testSens * (1 - cryoElig(s-1))) + (testSens * cryoElig(s-1) * (1-cryoEff(d)))) .* noVaxScreend;
                    dPop(noVaxToScreenTreatImm(d,:,aS)) = dPop(noVaxToScreenTreatImm(d,:,aS)) + ...
                        testSens * cryoElig(s-1) * cryoEff(d) * (1.0-((cryoHpvClear - (1-cryoEff(d)))/cryoEff(d))) .* noVaxScreend;
                    dPop(noVaxToScreenTreatHpv(d,:,aS)) = dPop(noVaxToScreenTreatHpv(d,:,aS)) + ...
                        testSens * cryoElig(s-1) * cryoEff(d) * ((cryoHpvClear - (1-cryoEff(d)))/cryoEff(d)) .* noVaxScreend;

                    vaxScreend = screenCover .* pop(vaxNoScreen(d,h,s,:,aS));
                    dPop(vaxNoScreen(d,h,s,:,aS)) = dPop(vaxNoScreen(d,h,s,:,aS)) - vaxScreend;
                    dPop(vaxToScreen(d,h,s,:,aS)) = dPop(vaxToScreen(d,h,s,:,aS)) + ...
                        ((1-testSens) + (testSens * (1 - cryoElig(s-1))) + (testSens * cryoElig(s-1) * (1-cryoEff(d)))) .* vaxScreend;
                    dPop(vaxToScreenTreatImm(d,:,aS)) = dPop(vaxToScreenTreatImm(d,:,aS)) + ...
                        testSens * cryoElig(s-1) * cryoEff(d) * (1.0-((cryoHpvClear - (1-cryoEff(d)))/cryoEff(d))) .* vaxScreend;
                    dPop(vaxToScreenTreatHpv(d,:,aS)) = dPop(vaxToScreenTreatHpv(d,:,aS)) + ...
                        testSens * cryoElig(s-1) * cryoEff(d) * ((cryoHpvClear - (1-cryoEff(d)))/cryoEff(d)) .* vaxScreend;
                end
            end
        end
    end
end

extraOut{1} = ccScreen;
dPop = sparse(dPop);