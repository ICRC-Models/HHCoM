baseCirc = load('H:\HHCoM_Results\baseCirc.mat');
circHigh = load('H:\HHCoM_Results\circHigh.mat');
c = fix(clock);
currYear = c(1); % get the current year
yearNow = round((currYear - startYear) * stepsPerYear);
%% Age-aggregated Disease Incidence
wVec = zeros(age , 1);
wVec(5 : age) = [0.188 , 0.18 , 0.159 , 0.121 , 0.088 , 0.067 , 0.054 , ...
    0.046 , 0.038 , 0.029 , 0.017 , 0.013]; 
figure()
newHiv_Arr = {baseCirc.newHiv , circHigh.newHiv};
popVec_Arr = {baseCirc.popVec , circHigh.popVec};
incMat = zeros(age , size(baseCirc.popVec , 1) / stepsPerYear);

inc = {incMat , incMat , incMat , incMat};
tVec = baseCirc.tVec;
for i = 1 : length(newHiv_Arr)
    for a = 4 : age
        newHiv = sum(sum(sum(newHiv_Arr{i}(1 : end , 1 : gender , a , :)...
            ,2),3),4);
        popVec = popVec_Arr{i};
        hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 1 : gender , a , 1 : risk)); ...
            toInd(allcomb(7 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 1 : gender , a , 1 : risk))];
        hivSus = sum(popVec(1 : end , hivSusInds) , 2);
        hivSus_Year = sum(reshape(hivSus , stepsPerYear , size(hivSus , 1) ...
            / stepsPerYear)) ./ stepsPerYear; % average susceptible population size per year
        newHiv_Year = sum(reshape(newHiv , stepsPerYear , size(newHiv , 1) ...
            /stepsPerYear)); % total new HIV infections per year
        inc{i}(a , :) = newHiv_Year ./ hivSus_Year .* 100;
    end
end

for i = 1 : length(newHiv_Arr)
    incAS = sum(bsxfun(@times , inc{i} , wVec));
    plot(tVec(1 : stepsPerYear : end) , incAS)
    xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
    hold on
end
xlabel('Year'); ylabel('Incidence per 100'); title('HIV Incidence')
% Women aged 16-29:  35% on treatment (30% suppressed)
% Women aged 30+: 60% on treatment (55% suppressed)
% Men aged 16-29: 25% on treatment (20% suppressed)
% Men aged 30+: 50% on treatment (40% suppressed).
legend('Base (40% Circumcision)' , '90% 16-29 YO Males circumcised' , 'northeastoutside')

%% ART uptake
figure()
hivInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 1 : age , 1 : risk)); toInd(allcomb(10 , 6 , ...
    1 : hpvTypes , 1 : hpvStates, 1 : periods , 1 : 2 , 1 : age , 1 : risk))];
artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 1 : age , 1 : risk));
plot(tVec , sum(baseCirc.popVec(: , artInds) , 2) ./ sum(baseCirc.popVec(: , hivInds) , 2) * 100 ,...
    tVec , sum(circHigh.popVec(: , artInds) , 2) ./ sum(circHigh.popVec(: , hivInds) , 2) * 100)

legend('Base', '90% 16-29 YO Males circumcised' , ...
    'Location' , 'northeastoutside')
title('Overall ART Coverage')
xlabel('Year'); ylabel('Coverage (%)')

%% HIV prevalence
hivInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk)); toInd(allcomb(10 , 6 , ...
    1 : hpvTypes , 1 : hpvStates, 1 : periods , 1 : 2 , 4 : 10 , 1 : risk))];
allInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk))];
artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk));
figure()
plot(tVec , sum(baseCirc.popVec(: , hivInds) , 2) ./ sum(baseCirc.popVec(: , allInds) , 2) * 100 ,...
    tVec , sum(circHigh.popVec(: , hivInds) , 2) ./ sum(circHigh.popVec(: , allInds) , 2) * 100)
legend('Base', '90% 16-29 YO Males circumcised' , ...
    'Location' , 'northeastoutside')
title('HIV Prevalence')
xlabel('Year'); ylabel('Prevalence (%)')

%% Circumcised
circInds = [toInd(allcomb(7, 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 , 1 : age, 1 : risk))];
allMInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 , 1 :  age , 1 : risk))];
figure()
plot(tVec , sum(baseCirc.popVec(: , circInds) , 2) ./ sum(baseCirc.popVec(: , allMInds) , 2) * 100 ,...
    tVec , sum(circHigh.popVec(: , circInds) , 2) ./ sum(circHigh.popVec(: , allMInds) , 2) * 100)
legend('Base', '90% 16-29 YO Males circumcised' , ...
    'Location' , 'northeastoutside')
title('Circumcision')
xlabel('Year'); ylabel('Coverage (%)')