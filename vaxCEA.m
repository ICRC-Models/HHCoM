function vaxCEA()

%% load results
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'])
nSims = size(dir([pwd , '\HHCoM_Results\Vaccine\*.mat']) , 1);
curr = load([pwd , '\HHCoM_Results\toNow.mat']); % Population up to 2018

% helper functions
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

% time
c = fix(clock); % get time
currYear = c(1); % get the current year from time

vaxResult = cell(nSims , 1);

parfor n = 1 : nSims
    % load results from vaccine run into cell array
    vaxResult{n} = load([pwd , '\HHCoM_Results\Vaccine\vaxSimResult' ,...
        num2str(n), '.mat']);
    % concatenate vectors/matrices of population up to 2018 to population
    % matrices for years past 2018
    vaxResult{n}.popVec = [curr.popVec ; vaxResult{n}.popVec];
    vaxResult{n}.newHpv= [curr.newHpv ; vaxResult{n}.newHpv];
    vaxResult{n}.newCC = [curr.newCC ; vaxResult{n}.newCC];
    vaxResult{n}.newHiv = [curr.newHiv ; vaxResult{n}.newHiv];
    vaxResult{n}.tVec = [curr.tVec , vaxResult{n}.tVec];
end

% Find no vaccine scenario
noVaxInd = -1;
for n = 1 : nSims
    if vaxResult{n}.vaxEff == 0
        noVaxInd = n;
    end
end
noV = vaxResult{noVaxInd};
tVec = noV.tVec;

%%
reset(0)
set(0 , 'defaultlinelinewidth' , 2)
%% Calculate life years saved

yrIntStart = 2018;

% CC Costs
ccCost = [2617 , 8533 , 8570]; % local, regional, distant

% all scenarios
for n = 1 : nSims
    vaxResult{n}.popVecAnnAvg = annAvg(sum(vaxResult{n}.popVec , 2)); % average population in the year
    vaxResult{n}.ccTreatedAnn = annlz(vaxResult{n}.ccTreated); % annual number of people treated for cervical cancer
end
% no vaccine scenario
noV.popVecAnnAvg = annAvg(noV.popVec); % average population in the year
noV.ccTreatedAnn = annlz(noV.ccTreated); % annual cost of cervical cancer treatment

for i = currYear : stepsPerYear: tVec(end)
    for n = 1 : length(vaxResult)
        % If y = current year, count benefits and CC treatment costs for women aged
        % >= y - B, where B = last year eligible for inclusion
        % Since 5 year age groups are being used, at ea ch year y, count benefits
        % for women in age group (round((y-B)/5)) and above.
        a = min(max(round((tVec(i) - yrIntStart) / 5) , 1) , age);
        ageCounted = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 2 , a : age , 1 : risk));
        
        % CC Indices
        localInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 5 , 1 : periods , ...
            1 : gender , a : age , 1 : risk));
        regionalInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 6 , 1 : periods , ...
            1 : gender , a : age , 1 : risk));
        distantInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 7 , 1 : periods , ...
            1 : gender , a : age , 1 : risk));
        
        % Count life years
        vaxResult{n}.ly(i) = sum(vaxResult{n}.popVecAnnAvg(i , ageCounted) , 2);
        
        % Count DALYs
        % Adjust life years for CC by region according to disability
        % weights
        % calculate CC DALYs for each time step
        cc_dalys = ...
            sum(vaxResult{n}.popVecAnnAvg(i , localInds) , 2) .* ccDalyWeight(1) ...
            + sum(vaxResult{n}.popVecAnnAvg(i , regionalInds) , 2) .* ccDalyWeight(2)...
            + sum(vaxResult{n}.popVecAnnAvg(i , distantInds) , 2) .* ccDalyWeight(3); 
        
        % DALYs obtained by subtracting full life years corresponding to CC
        % then adding DALYs corresponding to CC
        vaxResult{n}.daly(i) = vaxResult{n}.ly(i) ...
            - (sum(vaxResult{n}.popVecAnnAvg(i , localInds) , 2) ... % subtract full LY corresponding to CC 
            + sum(vaxResult{n}.popVecAnnAvg(i , regionalInds) , 2) ...
            + sum(vaxResult{n}.popVecAnnAvg(i , distantInds) , 2)) ...
            + cc_dalys; % Add CC DALYs


        % Cervical cancer costs
        vaxResult{n}.ccTreatedAnn(i) = ...
            sum(sum(sum(vaxResult{n}.ccTreatedAnn(i , : , : , a : age , 1) , 2) , 3) , 4) .* ccCost(1) + ...
            sum(sum(sum(vaxResult{n}.ccTreatedAnn(i , : , : , a : age , 2) , 2) , 3) , 4) .* ccCost(2) + ...
            sum(sum(sum(vaxResult{n}.ccTreatedAnn(i , : , : , a : age , 3) , 2) , 3) , 4) .* ccCost(3);
    end
    
    % no vaccine scenario
    % Count life years
    noV.ly(i) = sum(noV.popVecAnn(i , ageCounted) , 2);
    
    % Cervical cancer costs
    noV.ccCosts(i) = sum(sum(sum(noV.ccTreated(i , : , : , a : age , 1) , 2) , 3) , 4) .* ccCost(1) + ...
        sum(sum(sum(noV.ccTreated(i , : , : , a : age , 2) , 2) , 3) , 4) .* ccCost(2) + ...
        sum(sum(sum(noV.ccTreated(i , : , : , a : age , 3) , 2) , 3) , 4) .* ccCost(3);
end

for n = 1 : length(vaxResult)
    vaxResult{n}.lys = vaxResult{n}.ly - noV.ly;
end
%% Calculate annual costs
% HIV Costs (Not used in ICER calculation. Included for completeness)
hospCost = [117 , 56 , 38 , 38]; % <200 , 200-350, >350 , on ART
artCost = 260; 

% CC Costs (Incurred once per person at time of cervical cancer detection)
ccCost = [2617 , 8533 ,8570]; % local, regional, distant


% HIV Indices (Not used)
% above350Inds = toInd(allcomb(4 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%     1 : gender , 1 : age , 1 : risk));
% cd200_350Inds = toInd(allcomb(5 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%     1 : gender , 1 : age , 1 : risk));
% under200Inds = toInd(allcomb(6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%     1 : gender , 1 : age , 1 : risk));
% artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%     1 : gender , 1 : age , 1 : risk));

% HIV costs
% c90_2vFull.hivCosts = sum(c90_2vFull.popVec(: , above350Inds) , 2) .* hospCost(1) + ...
% sum(c90_2vFull.popVec(: , cd200_350Inds) , 2) .* hospCost(2) + ...
% sum(c90_2vFull.popVec(: , under200Inds) , 2) .* hospCost(3) + ...
% sum(c90_2vFull.popVec(: , artInds) , 2) .* (hospCost(4) + artCost);
% 
% c90_9vFull.hivCosts = sum(c90_9vFull.popVec(: , above350Inds) , 2) .* hospCost(1) + ...
% sum(c90_9vFull.popVec(: , cd200_350Inds) , 2) .* hospCost(2) + ...
% sum(c90_9vFull.popVec(: , under200Inds) , 2) .* hospCost(3) + ...
% sum(c90_9vFull.popVec(: , artInds) , 2) .* (hospCost(4) + artCost);

% CC Costs

% c90_2vFull.ccCosts = sum(sum(sum(c90_2vFull.ccTreated(: , : , : , : , 1) , 2) , 3) , 4) .* ccCost(1) + ...
% sum(sum(sum(c90_2vFull.ccTreated(: , : , : , : , 2) , 2) , 3) , 4) .* ccCost(2) + ...
% sum(sum(sum(c90_2vFull.ccTreated(: , : , : , : , 3) , 2) , 3) , 4) .* ccCost(3);
% 
% c90_9vFull.ccCosts = sum(sum(sum(c90_9vFull.ccTreated(: , : , : , : , 1) , 2) , 3) , 4) .* ccCost(1) + ...
% sum(sum(sum(c90_9vFull.ccTreated(: , : , : , : , 2) , 2) , 3) , 4) .* ccCost(2) + ...
% sum(sum(sum(c90_9vFull.ccTreated(: , : , : , : , 3) , 2) , 3) , 4) .* ccCost(3);

% Vaccination
cost2v = 27; % cost of 2 doses of bivalent vaccine
c90_2vFull.vaxCost = annlz(c90_2vFull.vaxd) * cost2v;

% NPV of vaccination cost
discountRate = 0.03; % discount rate of 3% per annum
c90_2vFull.vaxCostNPV = pvvar(c90_2vFull.vaxCost , discountRate);
%% Find price threshold

% ceThreshold = 1540; % USD per LYS
% priceGuess = 100;
% ce9v = @(x) abs(pvvar(annlz(c90_9vFull.vaxd) * x - annlz(c90_2vFull.vaxd) .* cost2v ...
%     + annAvg((c90_9vFull.hivCosts + c90_9vFull.ccCosts)) - annAvg((c90_2vFull.hivCosts + c90_2vFull.ccCosts)) , discountRate) ...
%     / pvvar(annAvg(c90_9vFull.lys) - annAvg(c90_2vFull.lys) , discountRate) - ceThreshold);
% priceThreshold_9v = fminsearch(ce9v , priceGuess);
% disp(['With a cost-effectiveness threshold of ' , num2str(ceThreshold) , ' USD, ' ,...
%     'the unit cost of 9v vaccine must be less than or equal to ' , ...
%     num2str(round(priceThreshold_9v , 2)),' USD.'])

%% CC only price threshold
% 3 thresholds: 0.5x GDP , 1x GDP , 500 USD per LYS
ceThreshold = 1540; % USD per LYS
ceThresholds = [0.5 * ceThreshold , ceThreshold , 500];
for i = 1 : length(ceThresholds)
    priceGuess = 100; % Enter a price guess for 9v to seed the search process
    % ce9v is an anonymous function that finds the vaccine price that
    % places the 9v vaccine right at the cost-effectiveness threshold
    % specified by ceThresholds(i)
    ce9v = @(x) abs(pvvar(annlz(c90_9vFull.vaxd) * x - annlz(c90_2vFull.vaxd) .* cost2v ...
        + annlz((c90_9vFull.ccCosts)) - annlz((c90_2vFull.ccCosts)) , discountRate) ...
        / pvvar(annlz(c90_9vFull.lys) - annlz(c90_2vFull.lys) , discountRate) - ceThresholds(i));
    priceThreshold_9v = fminsearch(ce9v , priceGuess);
    disp(['9v vs 2v: Considering only CC costs, with a cost-effectiveness threshold of ' , num2str(ceThresholds(i)) , ' USD, ' ,...
        'the unit cost of 9v vaccine must be less than or equal to ' , ...
        num2str(round(priceThreshold_9v , 2)),' USD.']) 
end

% 3 thresholds: 0.5x GDP , 1x GDP , 500 USD per LYS
ceThreshold = 1540; % USD per LYS
ceThresholds = [0.5 * ceThreshold , ceThreshold , 500];
for i = 1 : length(ceThresholds)
    priceGuess = 27; % Enter a price guess for 2v to seed the search process
    % ce9v is an anonymous function that finds the vaccine price that
    % places the 9v vaccine right at the cost-effectiveness threshold
    % specified by ceThresholds(i)
    ce9v = @(x) abs(pvvar(annlz(c90_2vFull.vaxd) * x ...
        + annlz((c90_2vFull.ccCosts)) - annlz((noV.ccCosts)) , discountRate) ...
        / pvvar(annlz(c90_2vFull.lys) , discountRate) - ceThresholds(i));
    priceThreshold_9v = fminsearch(ce9v , priceGuess);
    disp(['2v vs No Vaccine: Considering only CC costs, with a cost-effectiveness threshold of ' , num2str(ceThresholds(i)) , ' USD, ' ,...
        'the unit cost of 9v vaccine must be less than or equal to ' , ...
        num2str(round(priceThreshold_9v , 2)),' USD.']) 
end
%% YLS
figure()
plot(tVec(1 : stepsPerYear : end) , annAvg(c90_9vFull.lys) - annAvg(c90_2vFull.lys))
title('Years of Life Saved'); xlabel('Year'); ylabel('Years of life Saved')

figure()
plot(tVec(1 : stepsPerYear : end) , annlz(c90_9vFull.vaxd))
title('Vaccinated with 9v'); xlabel('Year'); ylabel('Number vaccinated')
%% CC incidence reduction

inds = {':' , [2 : 6 , 10] , [2 : 6] , 1 , 10};
files = {'CEA CC_General_Hpv_VaxCover' , 'CEA CC_HivAll_Hpv_VaxCover' , ...
    'CEA CC_HivNoART_Hpv_VaxCover' , 'CEA CC_HivNeg_Hpv_VaxCover' ,...
    'CEA CC_ART_HPV_VaxCover'};
plotTits = {'General' , 'HIV-Positive' , 'HIV-Positive (No ART)' , ....
    'HIV-Negative' , 'HIV-Positive on ART'};
fac = 10 ^ 5;
noV_Hpv = zeros(1 , length(tVec) / stepsPerYear);
% c70_2vPartial_Inc = noV_Hpv;
% c90_9vFullInc = noV_HpvAge;
% c90_2vFullInc = noV_HpvAge;
% v90_2vFullInc = noV_HpvAge;
for i = 1 : length(inds)
        % general
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 4 : age , 1 : risk))];
        % All HIV-positive women
        allHivF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 4 : age , 1 : risk));...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 4 : age , 1 : risk))];
        % HIV-positive women not on ART
        hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 4 : age , 1 : risk))];
        % All HIV-negative women
        hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
            2 , 4 : age , 1 : risk))];
        % Women on ART
        artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 4 : age , 1 : risk))];
        
        genArray = {allF , allHivF , hivNoARTF , hivNeg , artF};
        
        noV_Hpv = ...
            annlz(sum(sum(sum(noV.newCC(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
            (annlz(sum(noV.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac;
       
        c30_2vFullInc = ...
            annlz(sum(sum(sum(c30_2vFull.newCC(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
           	(annlz(sum(c30_2vFull.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
        c60_2vFullInc = ...
            annlz(sum(sum(sum(c60_2vFull.newCC(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
            (annlz(sum(c60_2vFull.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
        c90_2vFullInc = ...
            annlz(sum(sum(sum(c90_2vFull.newCC(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
            (annlz(sum(c90_2vFull.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
        c30_9vFullInc = ...
            annlz(sum(sum(sum(c30_9vFull.newCC(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
            (annlz(sum(c30_9vFull.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
        c60_9vFullInc = ...
            annlz(sum(sum(sum(c60_9vFull.newCC(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
            (annlz(sum(c60_9vFull.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
        c90_9vFullInc = ...
            annlz(sum(sum(sum(c90_9vFull.newCC(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
            (annlz(sum(c90_9vFull.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
        figure()
        plot(tVec(1 : stepsPerYear : end) , noV_Hpv , ...
            tVec(1 : stepsPerYear : end) , c30_2vFullInc , ...
            tVec(1 : stepsPerYear : end) , c60_2vFullInc , ...
            tVec(1 : stepsPerYear : end) , c90_2vFullInc , ...
            tVec(1 : stepsPerYear : end) , c30_9vFullInc , ...
            tVec(1 : stepsPerYear : end) , c60_9vFullInc , ...
            tVec(1 : stepsPerYear : end) , c90_9vFullInc)
    title([plotTits{i} , ' Cervical Cancer Incidence'])
    xlabel('Year'); ylabel('Incidence per 100,000')
    legend('No vaccination' , '30% Coverage (Full 2v)' , ...
        '60% coverage (Full 2v)' , '90% coverage (Full 2v)' , ...
        '30% Coverage (Full 9v)' , '60% coverage (Full 9v)' , ...
        '90% coverage (Full 9v)')
    % Reduction
    c30_2vFull_Red = (c30_2vFullInc - noV_Hpv) ./ noV_Hpv * 100;
    c60_2vFull_Red = (c60_2vFullInc - noV_Hpv) ./ noV_Hpv * 100;
    c90_2vFull_Red = (c90_2vFullInc - noV_Hpv) ./ noV_Hpv * 100;
    c30_9vFull_Red = (c30_9vFullInc - noV_Hpv) ./ noV_Hpv * 100;
    c60_9vFull_Red = (c60_9vFullInc - noV_Hpv) ./ noV_Hpv * 100;
    c90_9vFull_Red = (c90_9vFullInc - noV_Hpv) ./ noV_Hpv * 100;
    
    figure()
    plot(tVec(1 : stepsPerYear : end) , c30_2vFull_Red , ...
        tVec(1 : stepsPerYear : end) , c60_2vFull_Red , ...
        tVec(1 : stepsPerYear : end) , c90_2vFull_Red , ...
        tVec(1 : stepsPerYear : end) , c30_9vFull_Red , ...
        tVec(1 : stepsPerYear : end) , c60_9vFull_Red , ...
        tVec(1 : stepsPerYear : end) , c90_9vFull_Red)
    title([plotTits{i} , ' Cervical Cancer Incidence Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('30% Coverage (Full 2v)' , '60% coverage (Full 2v)' , ...
        '90% coverage (Full 2v)' , '30% Coverage (Full 9v)' , ...
        '60% coverage (Full 9v)' , '90% coverage (Full 9v)')
    axis([tVec(2) tVec(end) -100 0])
%     
%     T = table(tVec(1 : stepsPerYear : end)' , c30_2vFull_Inc' , ...
%         c90_9vFullInc' , c60_2vFullInc' , ...
%         c90_2vFull_Red' , c60_9vFull_Red' , c70_2vPartial_Red');
%     writetable(T , [files{i} , '.csv'] , 'Delimiter' , ',')
end
%%
figure()
for g = 1 : 2
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : 10 , 1 : risk));
    artPop = sum(noV.popVec(: , artInds) , 2);
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hivPop = sum(noV.popVec(: , hivInds) , 2);
    plot(tVec , 100 * artPop ./ (hivPop + artPop))
    hold on
end
xlabel('Year')
ylabel('Proportion of HIV Population')
title('Proportion on ART')
legend('Model (Male)' , 'Model (Female)')

%%
figure()
for g = 1 : 2
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : 10 , 1 : risk));
    artPop = sum(noV.popVec(: , artInds) , 2);
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk)); 
    hivPop = sum(noV.popVec(: , hivInds) , 2);
    allPop = sum(noV.popVec(: , allInds) , 2);
    plot(tVec , 100 * (hivPop + artPop) ./ allPop)
    hold on
end
xlabel('Year')
ylabel('Prevalence')
title('HIV Prevalence')
%%
hold on
for g = 1 : 2
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : 10 , 1 : risk));
    artPop = sum(c90_2vFull.popVec(: , artInds) , 2);
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk)); 
    hivPop = sum(c90_2vFull.popVec(: , hivInds) , 2);
    allPop = sum(c90_2vFull.popVec(: , allInds) , 2);
    plot(tVec , 100 * (hivPop + artPop) ./ allPop)
    hold on
end
legend('Male' , 'Female' , 'Male Vax' , 'Female Vax')

%%
figure()    
for g = 1 : 2
    hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : 10 , 1 : risk)); ...
        toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : 10 , 1 : risk))];
    hivSus = annlz(sum(c90_2vFull.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;    
    plot(tVec(1 : stepsPerYear : end) , ...
        annlz(sum(sum(c90_2vFull.newHiv(: , g , 4 : 10 , :) ...
        , 3) , 4)) ./ hivSus * 100)
    hold on
end

xlabel('Year'); ylabel('Rate Per 100'); title('HIV Incidence')
hold on

for g = 1 : 2
    hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : 10 , 1 : risk)); ...
        toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : 10 , 1 : risk))];
    hivSusNo = annlz(sum(noV.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
    plot(tVec(1 : stepsPerYear : end) , ...
        annlz(sum(sum(noV.newHiv(: , g , 4 : 10 , :) ...
        , 3) , 4)) ./ hivSusNo * 100 )
end
legend('Male' , 'Female' , 'Male No Vax' , 'Female No vax')
