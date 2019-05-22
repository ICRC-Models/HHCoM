unequalAll_Age = load('H:\HHCoM_Results\unequalAll_Age.mat');
equalMale_Age = load('H:\HHCoM_Results\equalMale_Age.mat');
equalFemale_Age = load('H:\HHCoM_Results\equalFemale_Age.mat');
equalAll_Age = load('H:\HHCoM_Results\equalAll_Age.mat');
c = fix(clock);
currYear = c(1); % get the current year
yearNow = round((currYear - startYear) * stepsPerYear);
%% Age-aggregated Disease Incidence
wVec = zeros(age , 1);
wVec(5 : age) = [0.188 , 0.18 , 0.159 , 0.121 , 0.088 , 0.067 , 0.054 , ...
    0.046 , 0.038 , 0.029 , 0.017 , 0.013]; 
figure()
newHiv_Arr = {unequalAll_Age.newHiv , equalMale_Age.newHiv , equalFemale_Age.newHiv , ...
    equalAll_Age.newHiv};
popVec_Arr = {unequalAll_Age.popVec , equalMale_Age.popVec , equalFemale_Age.popVec , ...
    equalAll_Age.popVec};
incMat = zeros(age , size(unequalAll_Age.popVec , 1) / stepsPerYear);

inc = {incMat , incMat , incMat , incMat};
tVec = unequalAll_Age.tVec;
for i = 1 : length(newHiv_Arr)
    for a = 4 : age
        newHiv = sum(sum(sum(newHiv_Arr{i}(1 : end , 1 : gender , a , :)...
            ,2),3),4);
        popVec = popVec_Arr{i};
        hivSusInds = toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 1 : gender , a , 1 : risk));
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
legend('Base' , 'Equal in Males' , 'Equal in Females' , 'Equal in All', ...
    'Location' , 'northeastoutside')

