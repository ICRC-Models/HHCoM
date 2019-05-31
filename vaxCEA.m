function vaxCEA(pathModifier)

waning = 0;    % turn waning on or off

%% Load parameters
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'],'stepsPerYear','circ','condUse','disease','viral',...
    'hpvTypes','hpvStates','periods','gender','age','risk','dim','k','toInd','sumall','modelYr1')

sumall = @(x) sum(x(:));

% Load results
nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
curr = load([pwd , '\HHCoM_Results\toNow_053119_noVax']); % Population up to current year

% Helper functions
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

% Time
c = fix(clock); % get time
currYear = c(1); % get the current year from time

vaxResult = cell(nSims , 1);
resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxSimResult'];
if waning
    resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxWaneSimResult'];
end
parfor n = 1 : nSims
    % load results from vaccine run into cell array
    vaxResult{n} = load([resultFileName , num2str(n), '.mat']);
    % concatenate vectors/matrices of population up to current year to population
    % matrices for years past current year
    vaxResult{n}.popVec = [curr.popVec(1 : end  , :) ; vaxResult{n}.popVec];
    vaxResult{n}.newHpv= [curr.newHpv(1 : end , : , : , : , :) ; vaxResult{n}.newHpv];
    vaxResult{n}.newImmHpv= [curr.newImmHpv(1 : end , : , : , : , :) ; vaxResult{n}.newImmHpv];
    vaxResult{n}.newVaxHpv= [curr.newVaxHpv(1 : end , : , : , : , :) ; vaxResult{n}.newVaxHpv];
    %vaxResult{n}.ccDeath = [curr.ccDeath(1 : end - 1 , : , : , :) ; vaxResult{n}.ccDeath];
    vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :) ; vaxResult{n}.newCC];
    vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , :) ; vaxResult{n}.newHiv];
    vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :) ; vaxResult{n}.artTreatTracker];
    vaxResult{n}.tVec = [curr.tVec(1 : end) , vaxResult{n}.tVec];
%     vaxResult{n}.ccTreated = [curr.ccTreated(1 : end - 1) , vaxResult{n}.ccTreated];
end

% Find no vaccine scenario
% noVaxInd = -1;
% for n = 1 : nSims
%     if vaxResult{n}.vaxEff == 0
%         noVaxInd = n;
%     end
% end
noVaxInd = nSims;

noV = vaxResult{noVaxInd};
tVec = noV.tVec;
tVecYr = tVec(1 : stepsPerYear : end);

%% Plot settings
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

%% Vaccination
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
% % 3 thresholds: 0.5x GDP , 1x GDP , 500 USD per LYS (BMGF)
% ceThreshold = 1540; % USD per LYS
% ceThresholds = [0.5 * ceThreshold , ceThreshold , 500];
% 
% % Using Life years
% % High coverage scenario (9v vs 2v)
% for i = 1 : length(ceThresholds)
%     priceGuess = 100; % Enter a price guess for 9v to seed the search process
%     % ce9v is an anonymous function that finds the vaccine price that
%     % places the 9v vaccine right at the cost-effectiveness threshold
%     % specified by ceThresholds(i)
%     ce9v = @(x) abs(pvvar(annlz(vaxResult{maxRate9vSim}.vaxd) * x - annlz(vaxResult{maxRate2vSim}.vaxd) .* cost2v ... % difference in vaccine cost for 9v vs 2v 
%         + annlz(vaxResult{maxRate9vSim}.ccCosts') - annlz(vaxResult{maxRate2vSim}.ccCosts') , discountRate) ... % difference in CC cost for 9v vs 2v scenario
%         / pvvar(annAvg(vaxResult{maxRate9vSim}.lys) - annAvg(vaxResult{maxRate2vSim}.lys) , discountRate) - ceThresholds(i)); % difference in LYS for 9v vs 2v scenario
%     priceThreshold_9v = fminsearch(ce9v , priceGuess);
%     fprintf(['\n 9v vs 2v: Considering only CC costs, with a cost-effectiveness \n' , ...
%         ' threshold of ' , num2str(ceThresholds(i)) , ' USD per LYS, ' ,...
%         ' the unit cost of 9v vaccine must be less than or equal to \n ' , ...
%         num2str(round(priceThreshold_9v , 2)),' USD. \n']) 
% end
% 
% disp(' ')
% % Using DALYs
% % High coverage scenario (9v vs 2v)
% for i = 1 : length(ceThresholds)
%     priceGuess = 100; % Enter a price guess for 9v to seed the search process
%     % ce9v is an anonymous function that finds the vaccine price that
%     % places the 9v vaccine right at the cost-effectiveness threshold
%     % specified by ceThresholds(i)
%     ce9v = @(x) abs(pvvar(annlz(vaxResult{maxRate9vSim}.vaxd) * x - annlz(vaxResult{maxRate2vSim}.vaxd) .* cost2v ... % difference in vaccine cost for 9v vs 2v 
%         + annlz(vaxResult{maxRate9vSim}.ccCosts') - annlz(vaxResult{maxRate2vSim}.ccCosts') , discountRate) ... % difference in CC cost for 9v vs 2v scenario
%         / pvvar(annAvg(vaxResult{maxRate9vSim}.daly) - annAvg(vaxResult{maxRate2vSim}.daly) , discountRate) - ceThresholds(i)); % difference in DALYs for 9v vs 2v scenario
%     priceThreshold_9v = fminsearch(ce9v , priceGuess);
%     fprintf(['\n 9v vs 2v: Considering only CC costs, with a cost-effectiveness \n' , ...
%         ' threshold of ' , num2str(ceThresholds(i)) , ' USD per DALY, ' ,...
%         ' the unit cost of 9v vaccine must be less than or equal to \n ' , ...
%         num2str(round(priceThreshold_9v , 2)),' USD.\n']) 
%end
% end

%% YLS

% figure()
% plot(tVec(1 : stepsPerYear : end) , annlz(c90_9vFull.vaxd))
% title('Vaccinated with 9v'); xlabel('Year'); ylabel('Number vaccinated')

%% CC incidence and incidence reduction
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]};
files = {'CC_General_Hpv_VaxCover' , ...
     'CC_HivNeg_Hpv_VaxCover' , 'CC_HivNoART_Hpv_VaxCover' ,...
     'CC_ART_HPV_VaxCover' , 'CC_HivAll_HPV_VaxCover'};
plotTits1 = {'General' , 'HIV-' , ...
    'HIV+ no ART' , 'HIV+ ART' , 'HIV all'};
plotTits2 = {'80% coverage: Total female population' , '90% coverage'  , ...
    '80% coverage: HIV-Negative' , '90% coverage' , ...
    '80% coverage: HIV-Positive no ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive all' , '90% coverage'};
fac = 10 ^ 5;
linStyle = {'-' , '--'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g'};

figure();

for i = 1 : length(inds)-1
%     plotTits = {plotTits2{(i*2-1):(i*2)}};
% %     figure();
% %     noV.ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
% %     for n = 1 : length(vaxResult)-1
% %         vaxResult{n}.ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
% %     end
%     
% %     % General, all ages
% %     allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
% %         1 : periods , 2 , 3 : age , 1 : risk)); ...
% %         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
% %         1 : periods , 2 , 3 : age , 1 : risk))];
% %     allhivNegFAge = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
% %             2 , 3 : age , 1 : risk)); ...
% %             toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
% %             2 , 3 : age , 1 : risk))];
    
% %     for a = 3 : age
        % General
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % All HIV-negative women
        hivNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , 3 : age , 1 : risk)); ...
            toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
            2 , 3 : age , 1 : risk))];
        % HIV-positive women not on ART
        hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % Women on ART
        artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % All HIV-positive women
        hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};

        ccIncRef = ...
            (annlz(sum(sum(sum(noV.newCC(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(noV.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac);
% %             .* (annlz(sum(noV.popVec(: , genArray{3}) , 2) ./ stepsPerYear));
% %         noV.ccIncRef = noV.ccIncRef + ccIncRef; 
        noV.ccIncRef = ccIncRef; 
                
        for n = 1 : length(vaxResult)-1
            ccIncRef = ...
                (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , 3 : age),2),3),4)) ./ ...
                (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
% %                 .* (annlz(sum(vaxResult{n}.popVec(: , genArray{3}) , 2) ./ stepsPerYear));
% %             vaxResult{n}.ccIncRef = vaxResult{n}.ccIncRef + ccIncRef;
            vaxResult{n}.ccIncRef = ccIncRef;
        end
        
% %     end
    
% %     noV.ccInc = noV.ccIncRef ./ (annlz(sum(noV.popVec(: , allhivNegFAge) , 2) ./ stepsPerYear));
    noV.ccInc = noV.ccIncRef;
    plot(tVec(1 : stepsPerYear : end) , noV.ccInc ,'DisplayName' , ...
         [plotTits1{i} , ': Efficacy: ' , num2str(round(noV.vaxEff * 100)) '% ,', ...
         'Coverage: ' , num2str(round(noV.vaxRate * 100)) , '%']);
    legend('-DynamicLegend');
    hold all;
    for n = 1 : length(vaxResult)-1
% %         vaxResult{n}.ccInc = vaxResult{n}.ccIncRef ./ (annlz(sum(vaxResult{n}.popVec(: , allhivNegFAge) , 2) ./ stepsPerYear));
        vaxResult{n}.ccInc = vaxResult{n}.ccIncRef;
        plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccInc , 'DisplayName' , ...
            [plotTits1{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
            'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
        
        % Reduction
%         vaxResult{n}.ccRed = (vaxResult{n}.ccInc - noV.ccInc) ./ noV.ccInc * 100;
%         plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccRed , 'LineStyle' , linStyle{n} , 'Color' , linColor{i} , 'DisplayName' , plotTits{n})
% %             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
% %             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
%         grid on
%         legend('-DynamicLegend')
%         xlim([2019 2099]);
%         xticks([2019 : 10 : 2099]);
%         hold all
              
        % Save reduction results
%         fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
%             'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
%         sname = [plotTits1{i} , '_IncRed'];
%         if exist(fname , 'file') == 2
%             M = xlsread(fname);
%             M = catpad(2 , [tVec(1 : stepsPerYear : end)' , noV.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , M);
%             xlswrite(fname , M , sname)
%         else
%             xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , sname)
%         end
        
    end
    title(' Cervical Cancer Incidence')
    xlabel('Year'); ylabel('Incidence per 100,000')
    hold all;

%     title('Figure 1: Percent reduction in cervical cancer cases, by HIV status')
%     xlabel('Year'); ylabel('Percent change')
end        
   
%             fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
%                 'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_' , pathModifier , '.xlsx'];
%             hold on
%             if exist(fname , 'file') == 2
%                 M = xlsread(fname);
%                 M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , M);
%                 xlswrite(fname , M , 'CC Incidence')
% %                 dlmwrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , ...
% %                     '-append' , 'delimiter' , ',' , 'coffset' , 1)
%             else
%                 xlswrite(fname , [tVec(1 : stepsPerYear : end)' , ...
%                     vaxResult{n}.ccInc'] , 'CC Incidence')
%             end

%% CC incidence by risk
% inds = {':' , [2 : 6] , 1 , 10};
% files = {'CC_General_Hpv_VaxCover' , ...
%      'CC_HivNoART_Hpv_VaxCover' , 'CC_HivNeg_Hpv_VaxCover' ,...
%      'CC_ART_HPV_VaxCover'};
% plotTits = {'General' , 'HIV-Positive (No ART)' , ....
%      'HIV-Negative' , 'HIV-Positive on ART'};
% fac = 10 ^ 5;
% 
% figure();
% for i = 2 : length(inds)
% %     figure();
%     % General, all risk groups
%     allFRisk = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
%         1 : periods , 2 , 8 : age , 1 : risk)); ...
%         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%         1 : periods , 2 , 8 : age , 1 : risk))];
%     % HIV-positive women not on ART
%     hivNoARTFRisk = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
%         1 : periods , 2 , 8 : age , 1 : risk)); ...
%         toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%         1 : periods , 2 , 8 : age , 1 : risk))];
%     % All HIV-negative women
%     hivNegRisk = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
%         2 , 8 : age , 1 : risk)); ...
%         toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
%         2 , 8 : age , 1 : risk))];
%     % Women on ART
%     artFRisk = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
%         1 : periods , 2 , 8 : age , 1 : risk)); ...
%         toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
%         1 : periods , 2 , 8 : age , 1 : risk))];
%     genArrayRisk = {allFRisk , hivNoARTFRisk , hivNegRisk , artFRisk};
%     
%     for r = 1 : risk
%         % General
%         allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
%             1 : periods , 2 , 8 : age , r)); ...
%             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 8 : age , r))];
%         % HIV-positive women not on ART
%         hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
%             1 : periods , 2 , 8 : age , r)); ...
%             toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 8 : age , r))];
%         % All HIV-negative women
%         hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
%             2 , 8 : age , r)); ...
%             toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
%             2 , 8 : age , r))];
%         % Women on ART
%         artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
%             1 : periods , 2 , 8 : age , r)); ...
%             toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 8 : age , r))];
%         genArray = {allF , hivNoARTF , hivNeg , artF};
% 
%         prop = ((annlz(sum(noV.popVec(: , genArray{i}) , 2))) ./ (annlz(sum(noV.popVec(: , genArrayRisk{i}) , 2)))) .* 100;
%             
%         
%         plot(tVec(1 : stepsPerYear : end) , prop ,'DisplayName' , ...
%              ['Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%              'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
%         ylim([0 100]);
%         hold all;
%         legend('Low risk' , 'Medium risk' , 'High risk');
%         
%     end
%     title([plotTits{i} , ' Proportion by Risk']);
%     xlabel('Year'); ylabel('Percent');
% %     hold off;
% end        

%% CC mortality and mortality reduction
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]};
files = {'CC_General_Hpv_VaxCover' , ...
     'CC_HivNeg_Hpv_VaxCover' , 'CC_HivNoART_Hpv_VaxCover' ,...
     'CC_ART_HPV_VaxCover' , 'CC_HivAll_HPV_VaxCover'};
plotTits1 = {'General' , 'HIV-' , ...
    'HIV+ no ART' , 'HIV+ ART' , 'HIV all'};
plotTits2 = {'80% coverage: Total female population' , '90% coverage'  , ...
    '80% coverage: HIV-Negative' , '90% coverage' , ...
    '80% coverage: HIV-Positive no ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive all' , '90% coverage'};
fac = 10 ^ 5;
linStyle = {'-' , '--'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g'};

figure();

for i = 1 : length(inds)-1
    plotTits = {plotTits2{(i*2-1):(i*2)}};
% %     figure();
% %     noV.ccMortRef = zeros(length(tVec(length(curr.tVec) + 1 : stepsPerYear : end)),1)';
% %     for n = 1 : length(vaxResult)-1
% %         vaxResult{n}.ccMortRef = zeros(length(tVec(length(curr.tVec) + 1 : stepsPerYear : end)),1)';
% %     end
    
    % General, all ages
% %     allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
% %         1 : periods , 2 , 3 : age , 1 : risk)); ...
% %         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
% %         1 : periods , 2 , 3 : age , 1 : risk))];
% %     allhivNegFAge = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
% %             2 , 3 : age , 1 : risk)); ...
% %             toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
% %             2 , 3 : age , 1 : risk))];
      
% %     for a = 3 : age
        % General
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % All HIV-negative women
        hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
            2 , 3 : age , 1 : risk))];
        % HIV-positive women not on ART
        hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % Women on ART
        artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % All HIV-positive women
        hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};

        ccMortRef = ...
            (annlz(sum(sum(sum(noV.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , genArray{i}) , 2) ./ stepsPerYear))* fac);
% %             .* (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , genArray{1}) , 2) ./ stepsPerYear));
% %         noV.ccMortRef = noV.ccMortRef + ccMortRef;
        noV.ccMortRef = ccMortRef;
                
        for n = 1 : length(vaxResult)-1
            ccMortRef = ...
                (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
                (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end  , genArray{i}) , 2) ./ stepsPerYear)) * fac);
% %                 .* (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end  , genArray{1}) , 2) ./ stepsPerYear));
% %             vaxResult{n}.ccMortRef = vaxResult{n}.ccMortRef + ccMortRef;
            vaxResult{n}.ccMortRef = ccMortRef;
        end
        
% %     end
% %     noV.ccMort = noV.ccMortRef ./ (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , allFAge) , 2) ./ stepsPerYear));
      noV.ccMort = noV.ccMortRef;
%     plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , noV.ccMort ,'DisplayName' , ...
%          [plotTits1{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%          'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
%     legend('-DynamicLegend');
%     hold all;
    for n = 1 : length(vaxResult)-1
% %         vaxResult{n}.ccMort = vaxResult{n}.ccMortRef ./ (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end , allFAge) , 2) ./ stepsPerYear));
        vaxResult{n}.ccMort = vaxResult{n}.ccMortRef;
%         plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , vaxResult{n}.ccMort , 'DisplayName' , ...
%             [plotTits1{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
        
        % Reduction
        vaxResult{n}.ccRed = (vaxResult{n}.ccMort - noV.ccMort) ./ noV.ccMort * 100;
        plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , vaxResult{n}.ccRed , 'LineStyle' , linStyle{n} , 'Color' , linColor{i} , 'DisplayName' , plotTits{n}) 
%             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
        grid on        
        legend('-DynamicLegend')
        xlim([2019 2099]);
        xticks([2019 : 10 : 2099]);
        hold all
        
        % Save reduction results
%         fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
%             'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_Mort' , '.xlsx'];
%         sname = [plotTits1{i} , '_MortRed'];
%         if exist(fname , 'file') == 2
%             M = xlsread(fname);
%             M = catpad(2 , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , noV.ccMort' , vaxResult{n}.ccMort' , vaxResult{n}.ccRed'] , M);
%             xlswrite(fname , M , sname)
%         else
%             xlswrite(fname , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , noV.ccMort' , vaxResult{n}.ccMort' , vaxResult{n}.ccRed'] , sname)
%         end
        
    end
%     title('Cervical Cancer Mortality')
%     xlabel('Year'); ylabel('Incidence per 100,000')
%     hold all;

    title('Figure 2: Percent reduction in cervical cancer mortality, by HIV status')
    xlabel('Year'); ylabel('Percent change')
end

%% Population Size
% HIV-positive women not on ART
hivNoART = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 7 , ...
    1 : periods , 1 : gender , 1 : age , 1 : risk)); ...
    toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
    1 : periods , 1 : gender , 1 : age , 1 : risk))];
% All HIV-negative women
hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 7 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk)); ...
    toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk))];
% Women on ART
art = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 7 , ...
    1 : periods , 1 : gender , 1 : age , 1 : risk)); ...
    toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
    1 : periods , 1 : gender , 1 : age , 1 : risk))];
genArray = {hivNoART , hivNeg , art};

figure()
subplot(1,2,1)
for i = 1 : length(genArray)
    plot(tVec , sum(noV.popVec(: , genArray{i}) , 2))
    hold all;
end
subplot(1,2,2)
plot(tVec , sum(noV.popVec(: , genArray{1}) , 2)+sum(noV.popVec(: , genArray{2}) , 2)+sum(noV.popVec(: , genArray{3}) , 2))
title('Population Size')
xlabel('Year'); ylabel('Individuals')
xlim([1910 2099]);
legend('HIV+ , no ART' , 'HIV-' , 'HIV+ , ART');
hold off

%% Proportion HIV population on ART
figure()
for g = 1 : 2
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 3 : age , 1 : risk));
    artPop = sum(noV.popVec(: , artInds) , 2);
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 3 : age , 1 : risk));
    hivPop = sum(noV.popVec(: , hivInds) , 2);
    plot(tVec , 100 * artPop ./ (hivPop + artPop))
    hold on
end
xlabel('Year')
ylabel('Proportion of HIV Population')
title('Proportion on ART')
legend('Model (Male)' , 'Model (Female)')
% 
%% HIV prevalance
figure()
for g = 1 : 2
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 1 : age , 1 : risk));
    artPop = sum(noV.popVec(: , artInds) , 2);
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 1 : age , 1 : risk));
    allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 1 : age , 1 : risk)); 
    hivPop = sum(noV.popVec(: , hivInds) , 2);
    allPop = sum(noV.popVec(: , allInds) , 2);
    plot(tVec , 100 * (hivPop + artPop) ./ allPop)
    hold on
end
xlabel('Year')
ylabel('Prevalence')
title('HIV Prevalence')

%% HPV prevalence by HIV group
figure()
linStyle = {'--' , '-' , ':'};
for a = 1 : age
    for r = 1 : risk
        % HIV+
        hpvHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 2:4 , 1, ...
             [1:2] , 2 , a , r));
        hpvHivPop = sum(noV.popVec(: , hpvHivInds) , 2);
        popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates ,  1 : periods , ...
            2 , a , r)));
        %ART
        hpvArtInds = toInd(allcomb(10 , 6 , 2:4 , 1, ...
             [1,2] , 2 , a , r));
        hpvArtPop = sum(noV.popVec(: , hpvArtInds) , 2);
        popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            2 , a , r)));
        %HIV-
        hpvHivNegInds = toInd(allcomb([1,7:9] , 1 , 2:4 , 1 , ...
             [1:2] , 2 , a , r));
        hpvHivNegPop = sum(noV.popVec(: , hpvHivNegInds) , 2);
        popHivNegTot = noV.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates ,  1 : periods , ...
            2 , a , r)));

        subplot(4,4,a)
        plot(tVec , 100 * hpvHivNegPop ./ sum(popHivNegTot , 2),linStyle{r})
        hold all
        plot(tVec , 100 * hpvHivPop ./ sum(popHivTot , 2),linStyle{r})
        hold all
        plot(tVec , 100 * hpvArtPop ./ sum(popArtTot , 2),linStyle{r})
        %axis([tVec(1) tVec(end) 0 100])
        xlim([2018 2050])
        xlabel('Year'); ylabel('Prevalence (%)'); title(' HPV Prevalence all ages')
    end
end
legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')

%% HPV/CIN/CC prevalence by HIV group
% figure()
% linStyle = {'--' , '-' , ':'};
% for a = 1 : age
%     for r = 1 : risk
%         % HIV+
%         hpvHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 2 : 4 , 1:7, ...
%              1 : periods , 2 , a , r));
%         hpvHivPop = sum(noV.popVec(: , hpvHivInds) , 2);
%         popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates ,  1 : periods , ...
%             2 , a , r)));
%         %ART
%         hpvArtInds = toInd(allcomb(10 , 6 , 2 : 4 , 1:7, ...
%              1 : periods , 2 , a , r));
%         hpvArtPop = sum(noV.popVec(: , hpvArtInds) , 2);
%         popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%             2 , a , r)));
%         %HIV-
%         hpvHivNegInds = toInd(allcomb([1,7:9] , 1 , 2 : 4 , 1:7 , ...
%              1 : periods , 2 , a , r));
%         hpvHivNegPop = sum(noV.popVec(: , hpvHivNegInds) , 2);
%         popHivNegTot = noV.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates ,  1 : periods , ...
%             2 , a , r)));
% 
%         subplot(4,4,a)
%         plot(tVec , 100 * hpvHivNegPop ./ sum(popHivNegTot , 2),linStyle{r})
%         hold all
%         plot(tVec , 100 * hpvHivPop ./ sum(popHivTot , 2),linStyle{r})
%         hold all
%         plot(tVec , 100 * hpvArtPop ./ sum(popArtTot , 2),linStyle{r})
%         %axis([tVec(1) tVec(end) 0 100])
%         xlim([1980 2100])
%         xlabel('Year'); ylabel('Prevalence (%)'); title(' CC Prevalence all ages')
%     end
% end
% legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')

%% CC prevalence by HIV group
% % HIV+
% ccHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 2 : 4 , 5 : 7, ...
%      [1,6] , 2 , 3 : age , 1 : risk));
% ccHivPop = sum(noV.popVec(: , ccHivInds) , 2);
% popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , [1:8,10] ,  [1,6] , ...
%     2 , 3 : age , 1 : risk)));
% %ART
% ccArtInds = toInd(allcomb(10, 6 , 2 : 4 , 5 : 7, ...
%      [1,6] , 2 , 3 : age , 1 : risk));
% ccArtPop = sum(noV.popVec(: , ccArtInds) , 2);
% popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , [1:8,10] ,  [1,6] , ...
%     2 , 3 : age , 1 : risk)));
% %HIV-
% ccHivNegInds = toInd(allcomb(1 , 1 , 2 : 4 , 5 : 7, ...
%      [1,6] , 2 , 3 : age , 1 : risk));
% ccHivNegPop = sum(noV.popVec(: , ccHivNegInds) , 2);
% popHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvTypes , [1:8,10] ,  [1,6] , ...
%     2 , 3 : age , 1 : risk)));
% 
% figure();
% plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'o')
% hold all
% plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'o')
% hold all
% plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'o')
% %axis([tVec(1) tVec(end) 0 100])
% xlabel('Year'); ylabel('Prevalence (%)'); title(' CC Prevalence')
% legend('HIV-' , 'HIV+ noART' , 'ART')

%% Vaccinated proportion by HIV group
figure();
linStyle = {'--' , '-' , ':'};
for a = 1 : age
    for r = 1 : risk
    % HIV+
    vaxHivInds = [toInd(allcomb(2 : 6 , 1 : 5 , 1 , 9 , [1,6] , 2 , a , r)); ...
        toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , [2,4] , 2 , a , r))];
    vaxHivPop = sum(vaxResult{2}.popVec(: , vaxHivInds) , 2);
    popHivTot = vaxResult{2}.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r)));
    %ART
    vaxArtInds = [toInd(allcomb(10 , 6 , 1 , 9 , [1,6] , 2 , a , r)); ...
        toInd(allcomb(10, 6 , 1 : hpvTypes , 1 : hpvStates , [2,4] , 2 , a , r))];
    vaxArtPop = sum(vaxResult{2}.popVec(: , vaxArtInds) , 2);
    popArtTot = vaxResult{2}.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r)));
    %HIV-
    vaxHivNegInds = [toInd(allcomb([1,7:9] , 1 , 1 , 9 , [1,6] , 2 , a , r)); ...
    toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , [2,4] , 2 , a , r))];
    vaxHivNegPop = sum(vaxResult{2}.popVec(: , vaxHivNegInds) , 2);
    popHivNegTot = vaxResult{2}.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r)));

    subplot(4,4,a)
    plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2) , linStyle{r})
    hold all
    plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2) , linStyle{r})
    hold all
    plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2) , linStyle{r})
    %axis([tVec(1) tVec(end) 0 100])
    xlabel('Year'); ylabel('Proportion (%)'); title('Vaccinated Proportion')
    xlim([2018 2100])
    
    hold all;
    end
end
legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')

%% Vaccinated proportion by CD4 count
figure();
linStyle = {'--' , '-' , ':'};
for c = 2 : 6
% HIV+
vaxHivInds = [toInd(allcomb(c , 1 : 5 , 1 , 9 , [1,6] , 2 , 7 , 1 : risk)); ...
    toInd(allcomb(c , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , [2,4] , 2 , 7 , 1 : risk))];
vaxHivPop = sum(vaxResult{2}.popVec(: , vaxHivInds) , 2);
popHivTot = vaxResult{2}.popVec(: , toInd(allcomb(c , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    2 , 7 , 1 : risk)));
%ART
vaxArtInds = [toInd(allcomb(10 , 6 , 1 , 9 , [1,6] , 2 , 7 , 1 : risk)); ...
   toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , [2,4] , 2 , 7 , 1 : risk))];
vaxArtPop = sum(vaxResult{2}.popVec(: , vaxArtInds) , 2);
popArtTot = vaxResult{2}.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    2 , 7 , 1 : risk)));
%HIV-
vaxHivNegInds = [toInd(allcomb([1,7:9] , 1 , 1 , 9 , [1,6] , 2 , 7 , 1 : risk)); ...
    toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , [2,4] , 2 , 7 , 1 : risk))];
vaxHivNegPop = sum(vaxResult{2}.popVec(: , vaxHivNegInds) , 2);
popHivNegTot = vaxResult{2}.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    2 , 7 , 1 : risk)));

% hold all
% plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2),'k--')
% hold all
plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2),'m:')
hold all
% plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2),':')
% hold all
axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Proportion (%)'); title('Vaccinated Proportion- age group 7')
xlim([2018 2100])
end

% legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')
legend('HIV-','ART','Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
   'CD4 <= 200 cells/uL')

%% Vaccinated proportion by risk
figure()
for r = 1 : risk
    %HIV-
    vaxHivNegInds = [toInd(allcomb([1,7:9] , 1 , 1 , 9 , [1,6] , 2 , 3 , r)); ...
    toInd(allcomb([1,7:9] , 1 , 2 : 4 , 1 : hpvStates , [2,4] , 2 , 3 , r))];
    vaxHivNegPop = sum(vaxResult{1}.popVec(: , vaxHivNegInds) , 2);
    popHivNegTot = vaxResult{1}.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , 3 , r)));
    % HIV+
    vaxHivInds = [toInd(allcomb(2 : 6 , 1 : 5 , 1 , 9 , [1,6] , 2 , 3 , r)); ...
        toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , [2,4] , 2 , 3 , r))];
    vaxHivPop = sum(vaxResult{1}.popVec(: , vaxHivInds) , 2);
    popHivTot = vaxResult{1}.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , 3 , r)));
    %ART
    vaxArtInds = [toInd(allcomb(10, 6 , 1 , 9 , [1,6] , 2 , 3 , r)); ...
        toInd(allcomb(10, 6 , 2 : 4 , 1 : hpvStates , [2,4] , 2 , 3 , r))];
    vaxArtPop = sum(vaxResult{1}.popVec(: , vaxArtInds) , 2);
    popArtTot = vaxResult{1}.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , 3 , r)));

    plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2),'o')
    hold all
    plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2),'o')
    hold all
    plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2),'o')
    %axis([tVec(1) tVec(end) 0 100])
    xlabel('Year'); ylabel('Proportion (%)'); title('Vaccinated Proportion')
    legend('lr HIV-' , 'lr HIV+' , 'lr ART' , 'mr HIV-' , 'mr HIV+' , 'mr ART' , 'hr HIV-' , 'hr HIV+' , 'hr ART')
end

%% Vaccinateable prevalence by HIV group
% figure();
% linStyle = {'--' , '-' , ':'};
% for a = 1 : age
%     for r = 1 : risk
%         % HIV+
%         vaxHivInds = [toInd(allcomb(2 : 6 , 1 : 5 , 1 , 1 , [1,6] , 2 , a , r)); ...
%             toInd(allcomb(2 : 6 , 1 : 5 , 2:4 , 10 , [1,6] , 2 , a , r))];
%         vaxHivPop = sum(noV.popVec(: , vaxHivInds) , 2);
%         popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%             2 , a , r)));
%         %ART
%         vaxArtInds = [toInd(allcomb(10, 6 , 1 , 1 , [1,6] , 2 , a , r)); ...
%             toInd(allcomb(10, 6 , 2:4 , 10 , [1,6] , 2 , a , r))];
%         vaxArtPop = sum(noV.popVec(: , vaxArtInds) , 2);
%         popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%             2 , a , r)));
%         %HIV-
%         vaxHivNegInds = [toInd(allcomb([1,7:9] , 1 , 1 , 1 , [1,6] , 2 , a , r)); ...
%             toInd(allcomb([1,7:9] , 1 , 2:4 , 10 , [1,6] , 2 , a , r))];
%         vaxHivNegPop = sum(noV.popVec(: , vaxHivNegInds) , 2);
%         popHivNegTot = noV.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%             2 , a , r)));
% 
%         subplot(4,4,a)
%         plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2),linStyle{r})
%         hold all
%         plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2),linStyle{r})
%         hold all
%         plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2),linStyle{r})
%         %axis([tVec(1) tVec(end) 0 100])
%         xlabel('Year'); ylabel('Prevalence (%)'); title(' Vaccinateable Prevalence')
%         ylim([0 100])
%     end
% end
% legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')

%% Population by risk and HIV group
figure();
for a = 1 : age
for r = 1 : risk
    % HIV+
    rpopHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r)));
    popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk)));
    %ART
    rpopArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r)));
    popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk)));
    %HIV-
    rpopHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r)));
    popHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk)));

    subplot(4,4,a)
    plot(tVec , 100 * sum(rpopHivNegTot , 2) ./ sum(popHivNegTot , 2))
    hold all
    plot(tVec , 100 * sum(rpopHivTot , 2) ./ sum(popHivTot , 2))
    hold all
    plot(tVec , 100 * sum(rpopArtTot , 2) ./ sum(popArtTot , 2))
    xlabel('Year'); ylabel('Prevalence (%)'); title(' Risk Proportion')
    hold all;
end
end
legend('HIV- : lr' , 'HIV+ noART : lr' , 'ART : lr' , 'HIV- : mr' , 'HIV+ noART : mr' , 'ART : mr' , 'HIV- : hr' , 'HIV+ noART : hr'  , 'ART : hr')

%% Population by age and HIV group
% for a = 3 : 6
%     % HIV+
%     rpopHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , a , 1 : risk)));
%     popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , 3 : age , 1 : risk)));
%     %ART
%     rpopArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , a , 1 : risk)));
%     popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , 3 : age , 1 : risk)));
%     %HIV-
%     rpopHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , a , 1 : risk)));
%     popHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , 3 : age , 1 : risk)));
% 
%     figure();
%     % plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
%     plot(tVec , 100 * sum(rpopHivNegTot , 2) ./ sum(popHivNegTot , 2))
%     hold all
%     plot(tVec , 100 * sum(rpopHivTot , 2) ./ sum(popHivTot , 2))
%     hold all
%     plot(tVec , 100 * sum(rpopArtTot , 2) ./ sum(popArtTot , 2))
%     xlabel('Year'); ylabel('Prevalence (%)'); title([num2str(a) , ': Age Distribution'])
%     legend('HIV-' , 'HIV+ noART' , 'ART' )
%     hold all;
% end
    
%% Population by "p"
figure();
subplot(2,2,1);
for p = 1 : periods
    % General
    inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
         p , 1 : gender , 8 , 1 : risk));
    pop = sum(vaxResult{1}.popVec(: , inds) , 2);
    popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
         1 : periods , 1 : gender , 8 , 1 : risk)));
    plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
    xlabel('Year'); ylabel('Proportion (%)'); title(' p Proportion')
    legend('1' , '2' , '3' , '4' , '5' ,'6')
    hold all;
end

subplot(2,2,3);
inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
     [1,6] , 1 : gender , 8 , 1 : risk));
pop = sum(vaxResult{1}.popVec(: , inds) , 2);
popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
     1 : periods , 1 : gender , 8 , 1 : risk)));
plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
xlabel('Year'); ylabel('Proportion (%)'); title(' Proportion unvaccinated or vaccinated and not reinfected ')

subplot(2,2,4);
inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
     [2,4] , 1 : gender , 8 , 1 : risk));
pop = sum(vaxResult{1}.popVec(: , inds) , 2);
popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
     1 : periods , 1 : gender , 8 , 1 : risk)));
plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
xlabel('Year'); ylabel('Proportion (%)'); title(' Proportion vaccinated and reinfected ')

subplot(2,2,2);
inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
     [4,6] , 1 : gender , 8 , 1 : risk));
pop = sum(vaxResult{1}.popVec(: , inds) , 2);
popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
     1 : periods , 1 : gender , 8 , 1 : risk)));
plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
xlabel('Year'); ylabel('Proportion (%)'); title(' Proportion screened')

%% Screened proportion by HIV group
figure();
linStyle = {'--' , '-' , ':'};
for a = 1 : age
    for r = 1 : risk
    % HIV+
    vaxHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 4 : 6 , 2 , a , r));
    vaxHivPop = sum(noV.popVec(: , vaxHivInds) , 2);
    popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r)));
    %ART
    vaxArtInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 4 : 6 , 2 , a , r));
    vaxArtPop = sum(noV.popVec(: , vaxArtInds) , 2);
    popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r)));
    %HIV-
    vaxHivNegInds = toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , 4 : 6 , 2 , a , r));
    vaxHivNegPop = sum(noV.popVec(: , vaxHivNegInds) , 2);
    popHivNegTot = noV.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r)));

    subplot(4,4,a)
    plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2) , linStyle{r})
    hold all
    plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2) , linStyle{r})
    hold all
    plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2) , linStyle{r})
    %axis([tVec(1) tVec(end) 0 100])
    xlabel('Year'); ylabel('Proportion (%)'); title('Screened Proportion')
%     xlim([2018 2100])
    
    hold all;
    end
end
legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')

%%
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

%%
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

%% HPV incidence
inds = {':' , [2 : 6] , [1,7:9] , 10};
files = {'CC_General_Hpv_VaxCover' , ...
     'CC_HivNoART_Hpv_VaxCover' , 'CC_HivNeg_Hpv_VaxCover' ,...
     'CC_ART_HPV_VaxCover'};
plotTits = {'General' , 'HIV-Positive (No ART)' , ....
     'HIV-Negative' , 'HIV-Positive on ART'};
fac = 10 ^ 5;

figure();

for i = 2 : length(inds)
%     figure();
    noV.hpvIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
    for n = 1 : length(vaxResult)-1
        vaxResult{n}.hpvIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
    end
    
    % General, all ages
    allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , ...
        1 : periods , 2 , 3 : age , 1 : risk)); ...
        toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 10 , ...
        1 : periods , 2 , 3 : age , 1 : risk))];
    allhivNegFAge = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , 3 : age , 1 : risk)); ...
            toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
            2 , 3 : age , 1 : risk))];
    
    for a = 3 : age
        % General
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , ...
            1 : periods , 2 , a , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 10 , ...
            1 : periods , 2 , a , 1 : risk))];
        % HIV-positive women not on ART
        hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 , ...
            1 : periods , 2 , a , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 10 , ...
            1 : periods , 2 , a , 1 : risk))];
        % All HIV-negative women
        hivNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 , 1 , 1 : periods , ...
            2 , a, 1 : risk)); ...
            toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 10 , 1 : periods , ...
            2 , a , 1 : risk))];
        % Women on ART
        artF = [toInd(allcomb(10 , 6 , 1 , 1 , ...
            1 : periods , 2 , a , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 10 , ...
            1 : periods , 2 , a , 1 : risk))];
        genArray = {allF , hivNoARTF , hivNeg , artF};

        hpvIncRef = ...
            ((annlz(sum(sum(noV.newHpv(: , 2 , inds{i} , a , :),3),5))+annlz(sum(sum(noV.newImmHpv(: , 2 , inds{i} , a , :),3),5))+annlz(sum(sum(noV.newVaxHpv(: , 2 , inds{i} , a , :),3),5))) ./ ...
            (annlz(sum(noV.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac) ...
            .* (annlz(sum(noV.popVec(: , genArray{3}) , 2) ./ stepsPerYear));
        noV.hpvIncRef = noV.hpvIncRef + hpvIncRef;
                
%         for n = 1 : length(vaxResult)-1
%             hpvIncRef = ...
%                 ((annlz(sum(sum(vaxResult{n}.newHpv(: , 2 , inds{i} , a , :),3),5))+annlz(sum(sum(vaxResult{n}.newImmHpv(: , 2 , inds{i} , a , :),3),5))+annlz(sum(sum(vaxResult{n}.newVaxHpv(: , 2 , inds{i} , a , :),3),5))) ./ ...
%                 (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac) ...
%                 .* (annlz(sum(vaxResult{n}.popVec(: , genArray{3}) , 2) ./ stepsPerYear));
%             vaxResult{n}.hpvIncRef = vaxResult{n}.hpvIncRef + hpvIncRef;
%         end
        
    end
    noV.hpvInc = noV.hpvIncRef ./ (annlz(sum(noV.popVec(: , allhivNegFAge) , 2) ./ stepsPerYear));
    plot(tVec(1 : stepsPerYear : end) , noV.hpvInc ,'DisplayName' , ...
         [plotTits{i} , ': Efficacy: ' , num2str(round(noV.vaxEff * 100)) '% ,', ...
         'Coverage: ' , num2str(round(noV.vaxRate * 100)) , '%']);
    legend('-DynamicLegend');
    hold all;
%     for n = 1 : length(vaxResult)-1
%         vaxResult{n}.hpvInc = vaxResult{n}.hpvIncRef ./ (annlz(sum(vaxResult{n}.popVec(: , allhivNegFAge) , 2) ./ stepsPerYear));
%         plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.hpvInc , 'DisplayName' , ...
%             [plotTits{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
%     end
    title(' HPV Incidence')
    xlabel('Year'); ylabel('Incidence per 100,000')
    hold all;
end       


%% HIV incidence
fac = 10 ^ 5;

figure();

for a = 1 : age
for r = 1 : risk

    % All HIV-negative women
    hivNeg = toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r));

    noV.hivInc = ...
        (annlz(noV.newHiv(: , 2 , a , r)) ./ (annlz(sum(noV.popVec(: , hivNeg) , 2) ./ stepsPerYear))* fac);

    subplot(4,4,a)
    plot(tVec(1 : stepsPerYear : end) , noV.hivInc);
    title(' HIV Incidence')
    xlabel('Year'); ylabel('Incidence per 100,000')
    xlim([2018 2050])
    hold all;
end
end

legend('lr','mr','hr')

%% ART 'incidence'
fac = 10 ^ 5;

figure();

for a = 1 : age
    for r = 1 : risk

    % All HIV-positive women
    hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , r));

    noV.artInc = ...
        (annlz(sum(sum(noV.artTreatTracker(: , 2 : 6 , 1 : 5 , 2 , a , r),2),3)) ./ (annlz(sum(noV.popVec(: , hivPos) , 2) ./ stepsPerYear))* fac);

    subplot(4,4,a)
    plot(tVec(1 : stepsPerYear : end) , noV.artInc);
    title(' ART "Incidence"')
    xlabel('Year'); ylabel('Incidence per 100,000')
    xlim([2018 2100])
    hold all;
    end
end

legend('lr','mr','hr')

%% HPV incidence by age, HIV+ on ART
inds = {':' , [2 : 6] , 1 , 10};
fac = 10 ^ 5;

figure();
    
for a = 3 : 10
    for r = 1 : risk
    % Women on ART
    artF = [toInd(allcomb(10 , 6 , 1 , 1 , ...
        1 : periods , 2 , a , r)); ...
        toInd(allcomb(10 , 6 , 1 : hpvTypes , 10 , ...
        1 : periods , 2 , a , r))];

    noV.hpvInc = ...
        ((annlz(sum(noV.newHpv(: , 2 , inds{4} , a , r),3))+annlz(sum(noV.newImmHpv(: , 2 , inds{4} , a , r),3))+annlz(sum(noV.newVaxHpv(: , 2 , inds{4} , a , r),3))) ./ ...
        (annlz(sum(noV.popVec(: , artF) , 2) ./ stepsPerYear))* fac);
    subplot(1,3,r)
    plot(tVec(1 : stepsPerYear : end) , noV.hpvInc ,'DisplayName' , ...
     num2str(a));
    ylim([0 12*10^4])
    title(' HPV Incidence: HIV+ on ART')
    xlabel('Year'); ylabel('Incidence per 100,000')    
    hold all;
    end
end
legend('-DynamicLegend');

%% HPV incidence by age, HIV+
inds = {':' , [2 : 6] , 1 , 10};
fac = 10 ^ 5;

figure();
    
for a = 3 : 10
    for r = 1 : risk
    % HIV-positive women not on ART
    hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 , ...
        1 : periods , 2 , a , r)); ...
        toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 10 , ...
        1 : periods , 2 , a , r))];

    noV.hpvInc = ...
        ((annlz(sum(noV.newHpv(: , 2 , inds{2} , a , r),3))+annlz(sum(noV.newImmHpv(: , 2 , inds{2} , a , r),3))+annlz(sum(noV.newVaxHpv(: , 2 , inds{2} , a , r),3)))./ ...
        (annlz(sum(noV.popVec(: , hivNoARTF) , 2) ./ stepsPerYear))* fac);

    subplot(1,3,r)
    plot(tVec(1 : stepsPerYear : end) , noV.hpvInc ,'DisplayName' , ...
     num2str(a));
    ylim([0 12*10^4])
    title(' HPV Incidence: HIV+')
    xlabel('Year'); ylabel('Incidence per 100,000')
    hold all;
    end
end   
legend('-DynamicLegend');

%% HPV incidence by age, HIV-
inds = {':' , [2 : 6] , 1 , 10};
fac = 10 ^ 5;

figure();
    
for a = 3 : 10
    for r = 1 : risk
    % HIV-negative women
    hivNeg = [toInd(allcomb(1 , 1 : viral , 1 , 1 , 1 : periods , ...
            2 , a, r)); ...
            toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 10 , 1 : periods , ...
            2 , a , r))];

    noV.hpvInc = ...
        (annlz(sum(noV.newHpv(: , 2 , inds{3} , a , r),3)) ./ ...
        (annlz(sum(noV.popVec(: , hivNeg) , 2) ./ stepsPerYear))* fac);

    subplot(1,3,r)
    plot(tVec(1 : stepsPerYear : end) , noV.hpvInc ,'DisplayName' , ...
     num2str(a));
    title(' HPV Incidence: HIV-')
    xlabel('Year'); ylabel('Incidence per 100,000')
    ylim([0 5*10^4])
    hold all;
    end
end
legend('-DynamicLegend');

%% HIV by age group
hivAge = zeros(length(tVec) , 12);
for a = 1 : age
    hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , a , 1 : risk));
    hivArt = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , a , 1 : risk));
    hivNeg = toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , a , 1 : risk));
    hivAge(: , a) = sum(noV.popVec(: , hivPos) , 2) + sum(noV.popVec(: , hivArt) , 2);
end
hivPosAllInd = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
hivArtAllInd = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
hivNegAllInd = toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , 1 : risk));
hivPosAll = sum(noV.popVec(: , hivPosAllInd) , 2 ) + sum(noV.popVec(: , hivArtAllInd),2);

figure()
subplot(1 , 2 , 1)
area(tVec , bsxfun(@rdivide , hivAge , hivPosAll));
title('HIV Status by Age Group')
xlabel('Year')
ylabel('Proportion of HIV Positive')
legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')

% HIV by risk group
hivRisk = zeros(length(tVec) , risk);
for r = 1 : risk
    hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , r));
    hivArt = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , r));
    hivRisk(: , r) = sum(noV.popVec(: , hivPos) , 2) + sum(noV.popVec(: , hivArt) , 2);

end
subplot(1 , 2 , 2)
area(tVec , bsxfun(@rdivide , hivRisk , hivPosAll));
title('HIV Status by Risk Group')
xlabel('Year')
ylabel('Proportion of HIV Positive')
legend('Low' , 'Medium' , 'High' , 'Location' , 'NorthEastOutside')

%% ART treatment tracker- CD4/risk/age
cd4ARTFrac = zeros(length(tVec) , risk);
for i = 1 : length(tVec)
    %currTot = sumall(noV.artTreatTracker(i , 2 : 6 , 1 : 5 , 2 , 3 , 1 : risk));
    for c = 2 : 6
        curr = sumall(noV.artTreatTracker(i , c , 1 : viral , 2 , 3 : age , 1 : risk));
        %hivArt = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        %    1 : periods , 2 , 3 , r));
        %curr = sum(noV.popVec(i , hivArt));
%         hivPos = toInd(allcomb(c , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%             1 : periods , 2 , 4 , 1 : risk));
        hivPos = sumall(noV.artTreatTracker(i , 2 : 6 , 1 : viral , 2 , 3 : age , 1 : risk));
        %cd4ARTFrac(i , r) = (curr / (sum(noV.popVec(i , hivPos))+curr)) .* 100;
%         cd4ARTFrac(i , c-1) = (curr / (sum(noV.popVec(i , hivPos)))).* 100;
        cd4ARTFrac(i , c-1) = (curr / hivPos).* 100;
    end
end

figure()
plot(tVec , cd4ARTFrac)
legend('Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
   'CD4 <= 200 cells/uL' , 'Location' , 'NorthEastOutside')
%legend('1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , 'Location' , 'NorthEastOutside')
xlabel('Year')
ylabel('Initiated ART')

%% Population Size by age: vax vs. non-vax
figure()
for a = 1: age  
    subplot(4,4,a)
    % HIV-positive women not on ART
    hivNoARTnoV = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , [1:7,10] , ...
        [1,6] , 2 , a , 1 : risk));
    hivNoARTV = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        [2,4] , 2 , a , 1 : risk)); ...
        toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 , ...
        [1,6] , 2 , a , 1 : risk))];
    % All HIV-negative women
    hivNegnoV = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , [1:7,10] , [1,6] , ...
        2 , a , 1 : risk));
    hivNegV = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , [2,4] , ...
        2 , a , 1 : risk)); ...
        toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 , [1,6] , ...
        2 , a , 1 : risk))];
    % Women on ART
    artnoV = toInd(allcomb(10 , 6 , 1 : hpvTypes , [1:7,10] , ...
        [1,6] , 2 , a , 1 : risk));
    artV = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        [2,4] , 2 , a , 1 : risk)); ...
        toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 , ...
        [1,6] , 2 , a , 1 : risk))];
    genArraynoV = {hivNoARTnoV , hivNegnoV , artnoV};
    genArrayV = {hivNoARTV , hivNegV , artV};

    for i = 1 : length(genArraynoV)
        plot(tVec , sum(vaxResult{2}.popVec(: , genArraynoV{i}) , 2),'-')
        hold all;
    end
    set(gca,'ColorOrderIndex',1)
    for i = 1 : length(genArrayV)
        plot(tVec , sum(vaxResult{2}.popVec(: , genArrayV{i}) , 2),'--')
        hold all;
    end
    title('Population Size')
    xlabel('Year'); ylabel('Individuals')
    xlim([1910 2200]);
end
legend('HIV+ , no ART: no vax' , 'HIV-: no vax' , 'HIV+ , ART: no vax');

%% Population Size by age
figure()
for a = 1: age  
    subplot(4,4,a)
    % HIV-positive women not on ART
    hivNoART = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 2 , a , 1 : risk));
    % All HIV-negative women
    hivNeg = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));
    % Women on ART
    art = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 2 , a , 1 : risk));
    genArray = {hivNoART , hivNeg , art};

    for i = 1 : length(genArray)
        plot(tVec , sum(vaxResult{2}.popVec(: , genArray{i}) , 2),'-')
        hold all;
    end
    title('Population Size')
    xlabel('Year'); ylabel('Individuals')
    xlim([1910 2200]);
end
legend('HIV+ , no ART' , 'HIV-' , 'HIV+ , ART');

%% Population age distribution by vaccination status
% figure()
% linStyle = {'o' , '*'};
% linColor = {'[0, 0.4470, 0.7410]' , '[0.8500, 0.3250, 0.0980]' , '[0.9290, 0.6940, 0.1250]'};
% 
% for a = 3: age  
%     % HIV-positive women not on ART
%     hivNoARTnoV = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , [1:7,10] , ...
%         1 , 2 , a , 1 : risk));
%     hivNoARTV = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%         2 , 2 , a , 1 : risk)); ...
%         toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 , ...
%         1 , 2 , a , 1 : risk))];
%     hivNoARTnoVTot = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , [1:7,10] , ...
%         1 , 2 , 3 : age , 1 : risk));
%     hivNoARTVTot = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%         2 , 2 , 3 : age , 1 : risk)); ...
%         toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 , ...
%         1 , 2 , 3 : age , 1 : risk))];
%     % All HIV-negative women
%     hivNegnoV = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , [1:7,10] , 1 , ...
%         2 , a , 1 : risk));
%     hivNegV = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 2 , ...
%         2 , a , 1 : risk)); ...
%         toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 , 1 , ...
%         2 , a , 1 : risk))];
%     hivNegnoVTot = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , [1:7,10] , 1 , ...
%         2 , 3 : age , 1 : risk));
%     hivNegVTot = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 2 , ...
%         2 , 3 : age , 1 : risk)); ...
%         toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 , 1 , ...
%         2 , 3 : age , 1 : risk))];
%     % Women on ART
%     artnoV = toInd(allcomb(10 , 6 , 1 : hpvTypes , [1:7,10] , ...
%         1 , 2 , a , 1 : risk));
%     artV = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
%         2 , 2 , a , 1 : risk)); ...
%         toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 , ...
%         1 , 2 , a , 1 : risk))];
%     artnoVTot = toInd(allcomb(10 , 6 , 1 : hpvTypes , [1:7,10] , ...
%         1 , 2 , 3 : age , 1 : risk));
%     artVTot = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
%         2 , 2 , 3 : age , 1 : risk)); ...
%         toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 , ...
%         1 , 2 , 3 : age , 1 : risk))];
%     genArraynoV = {hivNoARTnoV , hivNegnoV , artnoV};
%     genArrayV = {hivNoARTV , hivNegV , artV};
%     genArraynoVTot = {hivNoARTnoVTot , hivNegnoVTot , artnoVTot};
%     genArrayVTot = {hivNoARTVTot , hivNegVTot , artVTot};
% 
%     for i = 1 : length(genArraynoV)
%         scatter(a , sum(vaxResult{2}.popVec(end-40*stepsPerYear , genArraynoV{i}) , 2)/sum(vaxResult{2}.popVec(end-70*stepsPerYear , genArraynoVTot{i}) , 2),'MarkerEdgeColor',linColor{i},'MarkerFaceColor','k')
%         hold all;
%     end
%     for i = 1 : length(genArrayV)
%         scatter(a , sum(vaxResult{2}.popVec(end-40*stepsPerYear , genArrayV{i}) , 2)/sum(vaxResult{2}.popVec(end-70*stepsPerYear , genArrayVTot{i}) , 2),'MarkerEdgeColor',linColor{i})
%         hold all;
%     end
% 
%     title('Population Size')
%     xlabel('Age'); ylabel('Age Proportion')
% end
% legend('HIV+ , no ART' , 'HIV-' , 'HIV+ , ART');

%% Screened by risk and HIV group
% figure();
% % for r = 1 : risk
%     for v = 1 : 2
%         for i = 1 : (length(vaxResult{1}.tVec) - length(curr.tVec))
%         %HIV-
%         rpopHivNegTot(i,1) = sumall(vaxResult{1}.newScreen(i , [1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : risk , v ));
% 
%         %HIV+
%         rpopHivTot(i,1) = sumall(vaxResult{1}.newScreen(i , 2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : risk , v));
% 
%         %ART
%         rpopArtTot(i,1) = sumall(vaxResult{1}.newScreen(i , 10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : risk , v));
%         end
% 
%         subplot(1,3,1)
%         plot([1:(length(vaxResult{1}.tVec)-length(curr.tVec))]' , rpopHivNegTot)
%         hold all
%         subplot(1,3,2)
%         plot([1:(length(vaxResult{1}.tVec)-length(curr.tVec))]' , rpopHivTot)
%         hold all
%         subplot(1,3,3)
%         plot([1:(length(vaxResult{1}.tVec)-length(curr.tVec))]' , rpopArtTot)
%         xlabel('Year'); ylabel('Count'); title('Number')
%         hold all;
%     end
% % end
% legend('1','2','3','4','5','6','7','8','9','10')

