% function vaxCEA()

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
    vaxResult{n}.popVec = [curr.popVec(1 : end  , :) ; vaxResult{n}.popVec];
    vaxResult{n}.newHpv= [curr.newHpv(1 : end , : , : , : , :) ; vaxResult{n}.newHpv];
    vaxResult{n}.newImmHpv= [curr.newImmHpv(1 : end , : , : , : , :) ; vaxResult{n}.newImmHpv];
%     vaxResult{n}.ccDeath = [curr.ccDeath(1 : end - 1 , : , : , :) ; vaxResult{n}.ccDeath];
    vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :) ; vaxResult{n}.newCC];
    vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , :) ; vaxResult{n}.newHiv];
    vaxResult{n}.tVec = [curr.tVec(1 : end) , vaxResult{n}.tVec];
%     vaxResult{n}.ccTreated = [curr.ccTreated(1 : end - 1) , vaxResult{n}.ccTreated];
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
tVecYr = tVec(1 : stepsPerYear : end);
%%
reset(0)
set(0 , 'defaultlinelinewidth' , 2)
%% Calculate life years saved

yrIntStart = 2018;
for n = 1 : length(vaxResult)
    vaxResult{n}.ly = zeros((length(tVec) - length(curr.tVec)) , 1);
    vaxResult{n}.daly = zeros((length(tVec) - length(curr.tVec)) , 1);
end
noV.ly = zeros((length(tVec) - length(curr.tVec)) , 1);
noV.daly = zeros((length(tVec) - length(curr.tVec)) , 1);
%% CC Costs
ccCost = [2617 , 8533 , 8570]; % local, regional, distant
ccDalyWeight = 1 - [0.288 , 0.288 , 0.288]; % corresponds to local, regional, distant CC

for i = 1 : (length(tVec) - length(curr.tVec))
    % If y = current year, count benefits and CC treatment costs for women aged
        % >= y - B, where B = last year eligible for inclusion
        % Since 5 year age groups are being used, at each year y, count benefits
        % for women in age group (round((y-B)/5)) and above.
        a = min(max(round((tVec(i + length(curr.tVec) - 1) - yrIntStart) / 5) , 1) , age);
        ageCounted = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 2 , a : age , 1 : risk));
    for n = 1 : length(vaxResult)    
        % CC Indices
        localInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 5 , 1 : periods , ...
            1 : gender , a : age , 1 : risk));
        regionalInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 6 , 1 : periods , ...
            1 : gender , a : age , 1 : risk));
        distantInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 7 , 1 : periods , ...
            1 : gender , a : age , 1 : risk));
        
        % Count life years
        vaxResult{n}.ly(i) = sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1  , ageCounted) , 2);
        
        % Count DALYs
        % Adjust life years for CC by region according to disability
        % weights
        % calculate CC DALYs for each time step
        cc_dalys = ...
            sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , localInds) , 2) .* ccDalyWeight(1) ...
            + sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , regionalInds) , 2) .* ccDalyWeight(2)...
            + sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , distantInds) , 2) .* ccDalyWeight(3); 
        
        % DALYs obtained by subtracting full life years corresponding to CC
        % then adding DALYs corresponding to CC (done since LY already
        % calculated above)
        vaxResult{n}.daly(i) = vaxResult{n}.ly(i) ...
            - (sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , localInds) , 2) ... % subtract full LY corresponding to CC 
            + sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , regionalInds) , 2) ...
            + sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , distantInds) , 2)) ...
            + cc_dalys; % Add CC DALYs


        % Cervical cancer costs
        vaxResult{n}.ccCosts(i) = ...
            sum(sum(sum(vaxResult{n}.ccTreated(i , : , : , a : age , 1) , 2) , 3) , 4) .* ccCost(1) + ...
            sum(sum(sum(vaxResult{n}.ccTreated(i , : , : , a : age , 2) , 2) , 3) , 4) .* ccCost(2) + ...
            sum(sum(sum(vaxResult{n}.ccTreated(i , : , : , a : age , 3) , 2) , 3) , 4) .* ccCost(3);
    end
    
    % no vaccine scenario
    % Count life years
    noV.ly(i) = sum(noV.popVec(i + length(curr.tVec) - 1 , ageCounted) , 2);
    
    cc_dalys = ...
        sum(noV.popVec(i + length(curr.tVec) - 1 , localInds) , 2) .* ccDalyWeight(1) ...
        + sum(noV.popVec(i + length(curr.tVec) - 1 , regionalInds) , 2) .* ccDalyWeight(2)...
        + sum(noV.popVec(i + length(curr.tVec) - 1 , distantInds) , 2) .* ccDalyWeight(3);
    
    % DALYs obtained by subtracting full life years corresponding to CC
    % then adding DALYs corresponding to CC (done since LY already
    % calculated above)
    noV.daly(i) = noV.ly(i) ...
        - (sum(noV.popVec(i + length(curr.tVec) - 1 , localInds) , 2) ... % subtract full LY corresponding to CC
        + sum(noV.popVec(i + length(curr.tVec) - 1 , regionalInds) , 2) ...
        + sum(noV.popVec(i + length(curr.tVec) - 1 , distantInds) , 2)) ...
        + cc_dalys; % Add CC DALYs
    
    % Cervical cancer costs
    noV.ccCosts(i) = sum(sum(sum(noV.ccTreated(i, : , : , a : age , 1) , 2) , 3) , 4) .* ccCost(1) + ...
        sum(sum(sum(noV.ccTreated(i , : , : , a : age , 2) , 2) , 3) , 4) .* ccCost(2) + ...
        sum(sum(sum(noV.ccTreated(i , : , : , a : age , 3) , 2) , 3) , 4) .* ccCost(3);
end
%% 
for n = 1 : length(vaxResult)
    vaxResult{n}.lys = vaxResult{n}.ly - noV.ly;
end

figure()
for n = 1 : length(vaxResult)
    plot(tVec(length(curr.tVec) + 1 : end) , vaxResult{n}.lys , ...
        'DisplayName' , ['Vaccine Efficacy: ' , ...
        num2str(round(vaxResult{n}.vaxEff * 100)) ,'%, ' , ...
        'Vaccine Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) ,'%'])
    hold on
    legend('-DynamicLegend' , 'Location' , 'NorthwestOutside')
end

figure()
for n = 1 : length(vaxResult)
    plot(tVec(length(curr.tVec) + 1 : end) , sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end , :) , 2)...
        -sum(noV.popVec(length(curr.tVec) + 1 : end , :),2), ...
        'DisplayName' , ['Vaccine Efficacy: ' , ...
        num2str(round(vaxResult{n}.vaxEff * 100)) ,'%, ' , ...
        'Vaccine Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) ,'%'])
    hold on
    legend('-DynamicLegend' , 'Location' , 'NorthwestOutside')
end

for n = 1 : length(vaxResult)
    vaxResult{n}.dalyPlus = vaxResult{n}.daly - noV.daly;
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

% Vaccination
discountRate = 0.03; % discount rate of 3% per annum
cost2v = 27; % cost of 2 doses of bivalent vaccine
import java.util.LinkedList
sims2v = LinkedList();
full9v = LinkedList();
full2v = LinkedList();
maxRate9vSim = -1;
maxRate2vSim = -1;
maxRate9v = -1;
maxRate2v = -1;
for n = 1 : length(vaxResult)
    if vaxResult{n}.vaxEff <= 0.7 && vaxResult{n}.vaxEff > 0 % Vaccine efficacy <= 70% corresponds to 2v
        sims2v.add(n); % save index of 2v scenarios
        vaxResult{n}.vax2vCost = annlz(vaxResult{n}.vaxd) * cost2v;
        vaxResult{n}.vax2vCostNPV = pvvar(vaxResult{n}.vax2vCost , discountRate); % NPV of vaccination cost
    end
    if vaxResult{n}.vaxEff == 0.85
        full9v.add(n);
        if vaxResult{n}.vaxRate > maxRate9v
            maxRate9vSim = n;
            maxRate9v = vaxResult{n}.vaxRate;
        end
    elseif vaxResult{n}.vaxEff == 0.65
        full2v.add(n);
        if vaxResult{n}.vaxRate > maxRate2v
            maxRate2vSim = n;
            maxRate2v = vaxResult{n}.vaxRate; 
        end
    end
end

%% Find price threshold for 9v (CC costs only)
% 3 thresholds: 0.5x GDP , 1x GDP , 500 USD per LYS (BMGF)
ceThreshold = 1540; % USD per LYS
ceThresholds = [0.5 * ceThreshold , ceThreshold , 500];

% Using Life years
% High coverage scenario (9v vs 2v)
for i = 1 : length(ceThresholds)
    priceGuess = 100; % Enter a price guess for 9v to seed the search process
    % ce9v is an anonymous function that finds the vaccine price that
    % places the 9v vaccine right at the cost-effectiveness threshold
    % specified by ceThresholds(i)
    ce9v = @(x) abs(pvvar(annlz(vaxResult{maxRate9vSim}.vaxd) * x - annlz(vaxResult{maxRate2vSim}.vaxd) .* cost2v ... % difference in vaccine cost for 9v vs 2v 
        + annlz(vaxResult{maxRate9vSim}.ccCosts') - annlz(vaxResult{maxRate2vSim}.ccCosts') , discountRate) ... % difference in CC cost for 9v vs 2v scenario
        / pvvar(annAvg(vaxResult{maxRate9vSim}.lys) - annAvg(vaxResult{maxRate2vSim}.lys) , discountRate) - ceThresholds(i)); % difference in LYS for 9v vs 2v scenario
    priceThreshold_9v = fminsearch(ce9v , priceGuess);
    fprintf(['\n 9v vs 2v: Considering only CC costs, with a cost-effectiveness \n' , ...
        ' threshold of ' , num2str(ceThresholds(i)) , ' USD per LYS, ' ,...
        ' the unit cost of 9v vaccine must be less than or equal to \n ' , ...
        num2str(round(priceThreshold_9v , 2)),' USD. \n']) 
end

disp(' ')
% Using DALYs
% High coverage scenario (9v vs 2v)
for i = 1 : length(ceThresholds)
    priceGuess = 100; % Enter a price guess for 9v to seed the search process
    % ce9v is an anonymous function that finds the vaccine price that
    % places the 9v vaccine right at the cost-effectiveness threshold
    % specified by ceThresholds(i)
    ce9v = @(x) abs(pvvar(annlz(vaxResult{maxRate9vSim}.vaxd) * x - annlz(vaxResult{maxRate2vSim}.vaxd) .* cost2v ... % difference in vaccine cost for 9v vs 2v 
        + annlz(vaxResult{maxRate9vSim}.ccCosts') - annlz(vaxResult{maxRate2vSim}.ccCosts') , discountRate) ... % difference in CC cost for 9v vs 2v scenario
        / pvvar(annAvg(vaxResult{maxRate9vSim}.daly) - annAvg(vaxResult{maxRate2vSim}.daly) , discountRate) - ceThresholds(i)); % difference in DALYs for 9v vs 2v scenario
    priceThreshold_9v = fminsearch(ce9v , priceGuess);
    fprintf(['\n 9v vs 2v: Considering only CC costs, with a cost-effectiveness \n' , ...
        ' threshold of ' , num2str(ceThresholds(i)) , ' USD per DALY, ' ,...
        ' the unit cost of 9v vaccine must be less than or equal to \n ' , ...
        num2str(round(priceThreshold_9v , 2)),' USD.\n']) 
end

%% YLS

% figure()
% plot(tVec(1 : stepsPerYear : end) , annlz(c90_9vFull.vaxd))
% title('Vaccinated with 9v'); xlabel('Year'); ylabel('Number vaccinated')


%% CC incidence reduction
 
inds = {':' , [2 : 6 , 10] , [2 : 6] , 1 , 10};
files = {'CEA CC_General_Hpv_VaxCover' , 'CEA CC_HivAll_Hpv_VaxCover' , ...
     'CEA CC_HivNoART_Hpv_VaxCover' , 'CEA CC_HivNeg_Hpv_VaxCover' ,...
     'CEA CC_ART_HPV_VaxCover'};
 plotTits = {'General' , 'HIV-Positive' , 'HIV-Positive (No ART)' , ....
     'HIV-Negative' , 'HIV-Positive on ART'};
 fac = 10 ^ 5;
% noV = zeros(1 , (length(tVec) - stepsPerYear + 1) / stepsPerYear);
% % c70_2vPartial_Inc = noV_Hpv;
% % c90_9vFullInc = noV_HpvAge;
% % c90_2vFullInc = noV_HpvAge;
% % v90_2vFullInc = noV_HpvAge;
figure()
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
        
        noV.ccInc = ...
            annlz(sum(sum(sum(noV.newCC(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
            (annlz(sum(noV.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac;
        
        figure()
        
        plot(tVec(1 : stepsPerYear : end) , noV.ccInc ,'DisplayName' , ...
            ['Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
            'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
        legend('-DynamicLegend')
        
        hold on
        
        for n = 1 : length(vaxResult)
            vaxResult{n}.ccInc = ...
            annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
            plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccInc , 'DisplayName' , ...
                ['Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
                'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
            fname = ['Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.csv'];
            hold on
            if exist(fname , 'file') == 2
                M = csvread(fname);
                M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , M);
                csvwrite(fname , M)
%                 dlmwrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , ...
%                     '-append' , 'delimiter' , ',' , 'coffset' , 1)
            else
                csvwrite(fname , [tVec(1 : stepsPerYear : end)' , ...
                    vaxResult{n}.ccInc'] , 0 , i)
            end
        end
        title([plotTits{i} , ' Cervical Cancer Incidence'])
        xlabel('Year'); ylabel('Incidence per 100,000')
        hold off
        % Reduction
        figure()
        for n = 1 : length(vaxResult)
            vaxResult{n}.ccRed = (vaxResult{n}.ccInc - noV.ccInc) ./ noV.ccInc * 100;
            plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccRed , 'DisplayName' , ...
                ['Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
                'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
            legend('-DynamicLegend')
            hold on
            fname = ['Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.csv'];
            if exist(fname , 'file') == 2
                M = csvread(fname);
                M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccRed'] , M);
                csvwrite(fname , M)
%                 dlmwrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccRed'] , ...
%                     '-append' , 'delimiter' , ',' , 'coffset' , 1)
            else
                csvwrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccRed']...
                    , 0 , length(inds) + i)
            end
            title([plotTits{i} , ' Cervical Cancer Incidence Reduction'])
            xlabel('Year'); ylabel('Reduction (%)')
        end
        
        hold off
    
end
% %%

%% CC Mortality 
%% CC mortality reduction
 
inds = {':' , [2 : 6 , 10] , [2 : 6] , 1 , 10};
files = {'CEA CC_General_Hpv_VaxCover' , 'CEA CC_HivAll_Hpv_VaxCover' , ...
     'CEA CC_HivNoART_Hpv_VaxCover' , 'CEA CC_HivNeg_Hpv_VaxCover' ,...
     'CEA CC_ART_HPV_VaxCover'};
 plotTits = {'General' , 'HIV-Positive' , 'HIV-Positive (No ART)' , ....
     'HIV-Negative' , 'HIV-Positive on ART'};
 fac = 10 ^ 5;
% noV = zeros(1 , (length(tVec) - stepsPerYear + 1) / stepsPerYear);
% % c70_2vPartial_Inc = noV_Hpv;
% % c90_9vFullInc = noV_HpvAge;
% % c90_2vFullInc = noV_HpvAge;
% % v90_2vFullInc = noV_HpvAge;
figure()
for i = 1 : length(inds)
    % general
    allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 2 , 4 : age , 1 : risk));
    % All HIV-positive women
    allHivF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 2 , 4 : age , 1 : risk)); ...
        toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 2 , 4 : age , 1 : risk))];
    % HIV-positive women not on ART
    hivNoARTF = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 2 , 4 : age , 1 : risk));
    % All HIV-negative women
    hivNeg = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , 4 : age , 1 : risk));
    % Women on ART
    artF = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 2 , 4 : age , 1 : risk));
    
    genArray = {allF , allHivF , hivNoARTF , hivNeg , artF};
    
    noV.ccMort = ...
        annlz(sum(sum(sum(noV.ccDeath(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
        (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , genArray{i}) , 2) ./ stepsPerYear))* fac;
    
    figure()
    
    plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , noV.ccMort ,'DisplayName' , ...
        ['Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
        'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
    legend('-DynamicLegend')
    
    hold on
    
    for n = 1 : length(vaxResult)
        vaxResult{n}.ccMort = ...
            annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , 4 : age),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , vaxResult{n}.ccMort , 'DisplayName' , ...
            ['Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
            'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
        hold on
        fname = ['Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
            'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.csv'];
        if exist(fname , 'file') == 2
            M = csvread(fname);
            M = catpad(2 , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , vaxResult{n}.ccMort'] , M);
            csvwrite(fname , M)
%             dlmwrite(fname , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , vaxResult{n}.ccMort'] , ...
%                 '-append' , 'delimiter' , ',' , 'coffset' , 1)
        else
            csvwrite(fname , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , vaxResult{n}.ccMort'] , 0 , ...
                2 * length(inds) + i)
        end
    end
    title([plotTits{i} , ' Cervical Cancer Mortality'])
    xlabel('Year'); ylabel('Mortality per 100,000')
    hold off
    % Reduction
    figure()
    for n = 1 : length(vaxResult)
        vaxResult{n}.ccMortRed = (vaxResult{n}.ccMort - noV.ccMort) ./ noV.ccMort * 100;
        plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end)  , vaxResult{n}.ccMortRed , 'DisplayName' , ...
            ['Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
            'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
        legend('-DynamicLegend')
        hold on
        title([plotTits{i} , ' Cervical Cancer Mortality Reduction'])
        xlabel('Year'); ylabel('Reduction (%)')
        fname = ['Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
            'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.csv'];
        if exist(fname , 'file') == 2
            M = csvread(fname);
            M = catpad(2 , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , vaxResult{n}.ccMortRed'] , M);
            csvwrite(fname , M)
%             dlmwrite(fname , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , vaxResult{n}.ccMortRed'] , ...
%                 '-append' , 'delimiter' , ',' ,'coffset' , 1)
        else
            csvwrite(fname , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , vaxResult{n}.ccMortRed'] , ...
                0 , 3 * length(inds) + i)
        end
    end
    hold off
end

%% Population Size
figure()
for n = 1 : length(vaxResult)
    plot(tVec , sum(vaxResult{n}.popVec , 2) , 'DisplayName' , ...
        ['Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
        'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
    legend('-DynamicLegend')
    hold on
end
title('Population Size')
xlabel('Year'); ylabel('Individuals')
hold off
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
% 
% %%
% figure()
% for g = 1 : 2
%     artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
%         1 : periods , g , 4 : 10 , 1 : risk));
%     artPop = sum(noV.popVec(: , artInds) , 2);
%     hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
%         1 : periods , g , 4 : 10 , 1 : risk));
%     allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
%         1 : periods , g , 4 : 10 , 1 : risk)); 
%     hivPop = sum(noV.popVec(: , hivInds) , 2);
%     allPop = sum(noV.popVec(: , allInds) , 2);
%     plot(tVec , 100 * (hivPop + artPop) ./ allPop)
%     hold on
% end
% xlabel('Year')
% ylabel('Prevalence')
% title('HIV Prevalence')
% %%
% hold on
% for g = 1 : 2
%     artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
%         1 : periods , g , 4 : 10 , 1 : risk));
%     artPop = sum(c90_2vFull.popVec(: , artInds) , 2);
%     hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
%         1 : periods , g , 4 : 10 , 1 : risk));
%     allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
%         1 : periods , g , 4 : 10 , 1 : risk)); 
%     hivPop = sum(c90_2vFull.popVec(: , hivInds) , 2);
%     allPop = sum(c90_2vFull.popVec(: , allInds) , 2);
%     plot(tVec , 100 * (hivPop + artPop) ./ allPop)
%     hold on
% end
% legend('Male' , 'Female' , 'Male Vax' , 'Female Vax')
% 
% %%
% figure()    
% for g = 1 : 2
%     hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
%         1 : periods , g , 4 : 10 , 1 : risk)); ...
%         toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
%         1 : periods , g , 4 : 10 , 1 : risk))];
%     hivSus = annlz(sum(c90_2vFull.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;    
%     plot(tVec(1 : stepsPerYear : end) , ...
%         annlz(sum(sum(c90_2vFull.newHiv(: , g , 4 : 10 , :) ...
%         , 3) , 4)) ./ hivSus * 100)
%     hold on
% end
% 
% xlabel('Year'); ylabel('Rate Per 100'); title('HIV Incidence')
% hold on
% 
% for g = 1 : 2
%     hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
%         1 : periods , g , 4 : 10 , 1 : risk)); ...
%         toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
%         1 : periods , g , 4 : 10 , 1 : risk))];
%     hivSusNo = annlz(sum(noV.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
%     plot(tVec(1 : stepsPerYear : end) , ...
%         annlz(sum(sum(noV.newHiv(: , g , 4 : 10 , :) ...
%         , 3) , 4)) ./ hivSusNo * 100 )
% end
% legend('Male' , 'Female' , 'Male No Vax' , 'Female No vax')
