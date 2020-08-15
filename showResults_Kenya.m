function[] = showResults(pathModifier)

%% Load parameters and results
paramDir = [pwd , '\Params\'];

[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd, annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , partnersMmult, maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4,...
    mue , mue2 , mue3 , mue4 , epsA_vec , epsR_vec , yr , ...
    hivOn , betaHIV_mod , hiv_hpvMult, muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , waning , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear ,circNatStartYear , vaxStartYear , baseline , cisnet , who , whob , ...
    circProtect , condProtect , MTCTRate , hyst , OMEGA , ...
    ccInc2012_dObs , cc_dist_dObs , cin3_dist_dObs , cin1_dist_dObs , ...
    hpv_dist_dObs , cinPos2007_dObs , cin1_2010_dObs , cin2_2010_dObs, ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_all_dObs , hpv_hiv2009_dObs , ...
    hivPrevF_dObs , hivPrevM_dObs , hivPrevAll_dObs, popAgeDist_dObs , totPopSize_dObs , hivCurr , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ageInd , riskInd , ...
    hivNegNonVMMCinds , hivNegVMMCinds , vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3  , dFertMat3, d_partnersMmult, riskAdj, d_riskAdj, ...
    deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2_Kenya(1 , 0 , [] , [] , []);

% Load results
resultsDir = [pwd , '\HHCoM_Results\'];
toNowName = ['toNow_21Jul_RR_1-5_HIVtrans-00075']
load([resultsDir ,toNowName]) %change from pathModifier to file name
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)

% excel output file 
filename = [pwd, '\Calibration_comparison_Kenya.xlsx']


%% Population size over time vs. validation data
totalPop0_79 = zeros(2, length(tVec));
for g = 1 : gender
% All HIV-negative women
hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , g , 1 : 16 , 1 : risk));
% HIV-positive women not on ART
hivNoART = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , g , 1 : 16 , 1 : risk));
% Women on ART
art = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , g , 1 : 16 , 1 : risk));

genArray = {hivNeg , hivNoART , art};

totalPop0_79(g, :) = sum(popVec(:,genArray{1}),2) + sum(popVec(:,genArray{2}),2) + sum(popVec(:,genArray{3}),2);
end
% Load calibration data from Excel (years, values)
file = [pwd , '/Config/Population_validation_targets_Kenya.xlsx'];
historicalPop0_79 = zeros(15,2);
futurePop = zeros(16,2);
historicalPop0_79(:,1) = xlsread(file , 'Demographics' , 'B80:P80'); % years
historicalPop0_79(:,2) = xlsread(file , 'Demographics' , 'B97:P97') .* 1000; % estimates
futurePop(:,1) = xlsread(file , 'Demographics' , 'B101:Q101'); % years
futurePop(:,2) = xlsread(file , 'Demographics' , 'B102:Q102') .* 1000; % projections

figure()
area(tVec , totalPop0_79' );
hold all;
plot(historicalPop0_79(:,1) , historicalPop0_79(:,2) , 'o');
hold all;
plot(futurePop(:,1) , futurePop(:,2) , 'o');
title('Kenya Population Size Ages 0-79')
xlabel('Year'); ylabel('Individuals')
xlim([1920 2120]);
legend('Model, male' ,'Model, female', 'Kenya historical estimates (UN)' , 'Kenya future projections (UN)', ...
    'Location', 'Northwest')
hold off
%%
totalPopVec = sum(totalPop0_79(1 : 2, :),1); 
%%
sheet = ['pop'];
cols1 = {toNowName};
cols2 = {'Year', 'Model pop'}; %, 'UN pop'};
xlswrite(filename, cols1, sheet, 'D1')
xlswrite(filename, cols2, sheet, 'D2')
xlswrite(filename, [totalPopVec(1 : stepsPerYear * 5 : end)'], sheet, 'E3')
% xlswrite(filename, [historicalPop0_79(:, 2)], sheet,'G8')
% xlswrite(filename, [futurePop(:, 1)], sheet,'E23')
% xlswrite(filename, [futurePop(:, 2)], sheet,'G23')

%% Population size by age vs. validation data

% % Load calibration data from Excel
file = [pwd , '/Config/Kenya_parameters_Feb20.xlsx'];
% years = xlsread(file , 'Demographics' , 'B91:F91');    % years
% kzn_popByage_yrs(: , :) = xlsread(file , 'Demographics' , 'M92:Q107').*1000;    % males and females by age in 1996-2019

years = 1990:2020;
ageGroup = {'0-9', '10-19' , '20-29' , '30-39' , '40-49' , '50-59', '60-79'};
popPropYrs = zeros(length(years),5);
popNumYrs = popPropYrs ;
popPropYrs_obs = zeros(7,8);
popNumYrs_obs = popPropYrs_obs;
popNumYrs_obs = xlsread(file , 'Population' , 'H178:O184'); %Number
popPropYrs_obs = xlsread(file , 'Population' , 'H169:O175'); %proportions

ageVec = {[1:2], [3:4] , [5:6] , [7:8] , [9:10] , [11:12], [13:16]};
for y = 1 : length(years)
    yearCurr = years(y);
    for aInd = 1 : length(ageVec)
        a = ageVec{aInd};
        popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1:2 , a , 1 : risk));
        popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1:2  , 1 : age , 1 : risk));
        popPropYrs(y,aInd) = sum(popVec(((yearCurr - startYear) * stepsPerYear +1) , popAge),2) ...
            ./ sum(popVec(((yearCurr - startYear) * stepsPerYear +1) , popTot),2); 
        popNumYrs(y,aInd) = sum(popVec(((yearCurr - startYear) * stepsPerYear +1) , popAge),2);            
    end
end
%%
figure();
plot(years , popPropYrs);
newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
colororder();
hold on;
plot(popPropYrs_obs(:, 1) , popPropYrs_obs(:, 2) , 'o', ...
    popPropYrs_obs(:, 1) , popPropYrs_obs(:, 3) , 'o', ... 
    popPropYrs_obs(:, 1) , popPropYrs_obs(:, 4) , 'o',...
    popPropYrs_obs(:, 1) , popPropYrs_obs(:, 5) , 'o', ... 
    popPropYrs_obs(:, 1) , popPropYrs_obs(:, 6) , 'o', ...
    popPropYrs_obs(:, 1) , popPropYrs_obs(:, 7) , 'o', ...
    popPropYrs_obs(:, 1) , popPropYrs_obs(:, 8) , 'o');
colororder()
ylabel('Proportions'); xlabel('Year'); title('Age distribution in Kenya'); 
legend('0-9', '10-19' , '20-29' , '30-39' , '40-49' , '50-59', '60-79',...
    '0-9 obs', '10-19 obs', '20-29 obs', '30-39 obs', '40-49 obs', ...
    '50-59 obs', '60-79 obs', 'Location', 'NorthEastOutside');
%%

figure()
plot(years , popNumYrs);
newcolors = {'#F00','#F80','#FF0','#0B0','#00F','#50F','#A0F'};
colororder();
hold on;
plot(popNumYrs_obs(:, 1) , popNumYrs_obs(:, 2) , 'o', ...
    popNumYrs_obs(:, 1) , popNumYrs_obs(:, 3) , 'o', ... 
    popNumYrs_obs(:, 1) , popNumYrs_obs(:, 4) , 'o',...
    popNumYrs_obs(:, 1) , popNumYrs_obs(:, 5) , 'o', ... 
    popNumYrs_obs(:, 1) , popNumYrs_obs(:, 6) , 'o', ...
    popNumYrs_obs(:, 1) , popNumYrs_obs(:, 7) , 'o', ...
    popNumYrs_obs(:, 1) , popNumYrs_obs(:, 8) , 'o');
colororder()
ylabel('Number of individuals'); xlabel('Year'); title('Age distribution in Kenya'); 
legend('0-9', '10-19' , '20-29' , '30-39' , '40-49' , '50-59', '60-79',...
    '0-9 obs', '10-19 obs' , '20-29 obs' , '30-39 obs' , '40-49 obs' , '50-59 obs', '60-79 obs', ...
    'Location', 'NorthEastOutside');
%%
sheet = ['Pop_by_Age'];
cols1 = {toNowName};
cols2 = {'Year','0-9', '10-19' , '20-29' , '30-39' , '40-49' , '50-59', '60-79'}; %,...
    %'0-9 obs', '10-19 obs' , '20-29 obs' , '30-39 obs' , '40-49 obs' , '50-59 obs', '60-79 obs'};
xlswrite(filename, cols1, sheet, 'A12')
xlswrite(filename, cols2, sheet, 'A13')
xlswrite(filename, cols2, sheet, 'J13')
xlswrite(filename, [years(1: 5 : end)', popPropYrs(1 : 5 : end, :)], sheet, 'A14')
xlswrite(filename,[years(1: 5 : end)', popNumYrs(1 : 5 : end, :)], sheet, 'J14')


%% Fertility
% Load validation data from Excel (years, values)
file = [pwd , '/Config/Population_validation_targets_Kenya.xlsx'];
fertilityVal = xlsread(file , 'Demographics' , 'B4:G33');

fertilityVec = [];
for y = 1 : stepsPerYear : length(tVec)
    year = tVec(y);
    fertilityAnl = fertility;
    if year > 1970 && year <= 1990
        dt = (year - 1970) * stepsPerYear;
        dFert = (fertility2 - fertility) ./ ((1990 - 1970) * stepsPerYear);
        fertilityAnl = fertility + dFert .* dt;
%     elseif year > 2000 && year <= 2010
%         fertilityAnl = fertility2;
    elseif year > 1990 && year <=2020
        dt = (year - 1990) * stepsPerYear;
        dFert = (fertility3 - fertility2) ./ ((2020 - 1990) * stepsPerYear);
        fertilityAnl = fertility2 + dFert .* dt;
    elseif year > 2020
        fertilityAnl = fertility3;
    end
    
    diseaseVec = {[1:2,8] , 3 , 4 , 5 , 6 , 7};
    aSum = 0;        
    for a = 4 : 10
        allDinds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : 3 , 1 : intervens , 2 , a , 1 : risk));
        allTot = sumall(popVec(y,allDinds));
        for d = 1 : length(diseaseVec)
            subDinds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : 3 , 1 : intervens , 2 , a , 1 : risk));
            dProp = sumall(popVec(y,subDinds)) / allTot;
            aSum = aSum + fertilityAnl(a,d)*dProp*5;
        end
    end
    
    fertilityVec = [fertilityVec; [year aSum]];
end

figure;
plot(fertilityVec(:,1) , fertilityVec(:,2) , '-');
hold all;
plot(fertilityVal(:,1) , fertilityVal(:,2) , 'o');
hold all;
plot(fertilityVal(:,1) , fertilityVal(:,3) , 'o');
hold all;
plot(fertilityVal(:,1) , fertilityVal(:,4) , 'o');
hold all;
plot(fertilityVal(:,1) , fertilityVal(:,5) , 'o');
hold all;
plot(fertilityVal(:,1) , fertilityVal(:,6) , 'o');
title('Total fertility rate');
legend('Model' , 'UN estimates (Kenya)' , 'Lower 95' , ...
    'Lower 80' , 'Upper 80' , 'Upper 95');
ylim([0 10]);
xlim([1920 2100]);

%%
sheet = ['Fertility'];
cols1 = {toNowName};
cols2 = {'Year', 'Model'} %, 'Year', 'UN', 'Lower 95', 'Lower 80', 'Upper 80', 'Upper 95'};
xlswrite(filename, cols1, sheet, 'H1')
xlswrite(filename, cols2, sheet, 'H2')
xlswrite(filename, [fertilityVec(1:5:end, :)], sheet, 'H3')
%xlswrite(filename, [fertilityVal], sheet, 'C8')

%% Kenya HIV prevalence vs. observed data by year
% Total HIV positive
hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 : 2 , 4 : 10 , 1 : risk));
hivPop = sum(popVec(: , hivInds) , 2);
artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 : 2 , 4 : 10 , 1 : risk));
art = sum(popVec(: , artInds) , 2);
popTot = popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
    1 : intervens , 1 : 2 , 4 : 10 , 1 : risk)));
hivPrev = (hivPop + art) ./ sum(popTot , 2) * 100;

% Compared to national HIV data
file = [pwd , '/Config/Kenya_parameters_Feb20.xlsx'];
HIV_Ken_spectrum = xlsread(file , ['HIV ' 'prevalence'] , 'B184:E212');

figure()

plot(tVec , hivPrev , HIV_Ken_spectrum(: , 1)' , HIV_Ken_spectrum(: , 2)' , '*')
hold on 
% yPosError = abs(upper_prevVal - prevVal);
% yNegError = abs(lower_prevVal - prevVal);
% errorbar(prevValYrs , prevVal , yNegError , yPosError , 'ms')
xlabel('Year'); ylabel('Proportion of Population (%)'); title('HIV Prevalence (Ages 15-49)')
legend('Model' , 'Kenya (Spectrum)')
xlim([1972 2020])
ylim([0, 20])
%%
sheet = ['HIV_prev'];
cols1 = {toNowName};
cols2 = {'Final'} %, 'ANC data (Kisumu)', 'Year', 'Spectrum data (Nyanza)'};
xlswrite(filename, cols1, sheet, 'K1')
xlswrite(filename, cols2, sheet, 'K2')
xlswrite(filename, [hivPrev(331:stepsPerYear:end) ], sheet, 'K3')
% xlswrite(filename, [HIV_ANC_Kisumu'], sheet, 'C3')
% xlswrite(filename, [HIV_Kenya_spectrum'], sheet, 'E3')


%% HIV prevalance, all ages by gender
figure()
hivObsGender = zeros(4,3)
hivObsGender(:,3) = [2003 2007 2009 2012]; 
hivObsGender(:,1) = [8.7 8.4 8.0 6.9]; 
hivObsGender(:,2) = [4.6 5.4 4.6 4.4]; 

% sheet = ['HIV_by_sex'];
% cols1 = {toNowName};
% cols2 = {'M, Final', 'F, Final'} %, 'Year', 'Males, DHS/KAIS', 'Females, DHS/KAIS',};
% xlswrite(filename, cols1, sheet, 'K1')
% xlswrite(filename, cols2, sheet, 'K2')

for g = 1 : 2
    artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , 4 : age , 1 : risk));
    artPop = sum(popVec(: , artInds) , 2);
    hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , g , 4 : age , 1 : risk));
    allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , g , 4 : age , 1 : risk)); 
    hivPop = sum(popVec(: , hivInds) , 2);
    allPop = sum(popVec(: , allInds) , 2);
    hivPrev_sex =  100 * (hivPop + artPop) ./ allPop;
    plot(tVec , hivPrev_sex)
    hold on
    plot(hivObsGender(:, 3), hivObsGender(:, g), 'o')
    hold on
%     cell1 = ['K', 'L'];
%     cell = ([cell1(g) +'3']);
%     xlswrite(filename, [hivPrev_sex(331:stepsPerYear:end)], sheet, cell)
end
xlabel('Year')
xlim([1980 2020])
ylabel('Prevalence')
title('HIV Prevalence (aged 14-59)')
legend('Males, model' , 'Males, DHS/KAIS', 'Females, model', 'Females, DHS/KAIS')

%xlswrite(filename, [hivObsGender(:, 3), hivObsGender(:, 1:2)], sheet, 'E3')


%% HIV prevalance, by sex and age over time 

hivSexAge = zeros(length(tVec), 7);
ageGroup = ["15-19"; "20-24" ; "25-29"; "30-34"; "35-39" ; "40-44" ; "45-49"];
yr = 1980 : 2020;
sheet = ['HIV_sex_age'];
cols1 = {toNowName,"Final" }
%cols2 = ["HPVCalib_2.5xSexAct"]; %, 'F, 005% riskAdj'} %, 'Year', 'Males, DHS/KAIS', 'Females, DHS/KAIS',};
xlswrite(filename, cols1, sheet, 'J11')
%xlswrite(filename, cols2, sheet, 'A12')
%xlswrite(filename, yr', sheet, 'P4')
xlswrite(filename, ['Year', 'Sex', ageGroup' ], sheet, 'J12')
for g = 1  : 2
    for a = 4: 10
        %a = aVec{aInd};
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivSexAge(:, g, a-3) =  sum(popVec(: , hivInds) , 2) ./ sum(popVec(:,ageInds),2)*100;
    end
end

sex_label = repmat(["M" "F"], 41, 1)
    xlswrite(filename, [yr', sex_label(:, 1), squeeze(hivSexAge(331:stepsPerYear:end, 1, :))], sheet, "J13")
    xlswrite(filename, [yr', sex_label(:, 2), squeeze(hivSexAge(331:stepsPerYear:end, 2, :))], sheet, "J54")
    

%xlswrite(filename, [hivObsGender(:, 3), hivObsGender(:, 1:2)], sheet, 'E3')

%% HIV prevalence overall
figure()
artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , 1 : gender , 4 : age , 1 : risk));
artPop = sum(popVec(: , artInds) , 2);
hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 : gender , 4 : age , 1 : risk));
allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 : gender , 4 : age , 1 : risk)); 
hivPop = sum(popVec(: , hivInds) , 2);
allPop = sum(popVec(: , allInds) , 2);
plot(tVec , 100 * (hivPop + artPop) ./ allPop)

xlabel('Year')
ylabel('Prevalence')
xlim([1980 2020])
title('HIV Prevalence (aged 14-59)')
%% Proportion HIV-negative males circumcised by broad age groups over time  - CHANGE vmmcYr
circPropYr_obs = vmmcYr(1:5);
circProp_obs = vmmcRate(1:5, :)' .* 100;
%circProp_obs = [0.0 0.0 0.0 0.0 0.0 0.0 0.0; circProp_obs];

ageVec = {4 , 5 , 6, [7:8] , [9:10], [11:age]}; % Ages: (15-19), (20-24), (25-49), (50+)
circProp = zeros(length(tVec) , length(ageVec));

figure()
for aInd = 1 : length(ageVec)
    a = ageVec{aInd};
    circInds = toInd(allcomb(2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    circPop = sum(popVec(: , circInds) , 2);
    hivNegInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    hivNegPop = sum(popVec(: , hivNegInds) , 2);
    circProp(: , aInd) = 100 * circPop ./ hivNegPop;
end
plot(tVec , circProp);
set(gca,'ColorOrderIndex',1)
hold on;
plot(circPropYr_obs , circProp_obs , 'o--');
xlabel('Year')
ylabel('Proportion of HIV-Negative Males Circumcised by Broad Age Groups (%)')
title('Circumcision Indicator')
xlim([1960 2020]);
grid on;
legend('15-19, model' , '20-24, model ' , '25-29, model' ,...
    '30-39, model', '40-49, model', '50+, model' , ...
    '15-19, observed' , '20-24, observed' , '25-29, observed' ,...
    '30-39, observed', '40-49, observed', '50+, observed' ,...
    'Location' , 'NorthWest');
  
  
%% Kenya HIV by age
genderVec = {'Males (on and off ART)' , 'Females (on and off ART)'};

hiv2007 = zeros(7 , 2);
hiv2003 = hiv2007;
hiv2009 = hiv2007;
hiv2012 = hiv2007;
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};
file = [pwd , '/Config/Kenya_parameters_Feb20.xlsx'];
Kais2007_Ken = zeros(7,2);
Kais2007_Ken(:,1) = xlsread(file , 'HIV prevalence' , 'D18:D24'); % Males estimates
Kais2007_Ken(:,2) = xlsread(file , 'HIV prevalence' , 'C18:C24'); % Female estimates
Kais2012_Ken = Kais2007_Ken;
Kais2012_Ken(:,1) = xlsread(file , 'HIV prevalence' , 'M34:M40'); % Males estimates
Kais2012_Ken(:,2) = xlsread(file , 'HIV prevalence' , 'J34:J40'); % Female estimates
DHS2003_Ken = hiv2007;
DHS2003_Ken(:,1) = xlsread(file , 'HIV prevalence' ,'F4:F10' ); % Males estimates
DHS2003_Ken(:,2) = xlsread(file , 'HIV prevalence' , 'C4:C10' ); % Female estimates
DHS2009_Ken = hiv2007;
DHS2009_Ken(:,1) = xlsread(file , 'HIV prevalence' , 'Q62:Q68'); % Males estimates
DHS2009_Ken(:,2) = xlsread(file , 'HIV prevalence' , 'N62:N68'); % Female estimates

figure;
%2003
for g = 1 : gender
    for a = 4: 10
        %a = aVec{aInd};
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hiv2003(a-3 , g) = (sum(popVec((2003 - startYear)* stepsPerYear,hivInds),2) ... 
            /sum(popVec((2003 - startYear)* stepsPerYear,ageInds),2))*100;
    end
end
    plot(1 : size(hiv2003 , 1) , hiv2003(: , :) , '-', 1 : size(hiv2003 , 1), DHS2003_Ken(:,:), '--o');
    hold all;
    xlabel('Age Group'); ylabel('HIV Prevalence')
    set(gca , 'xtick' , 1 : length(hiv2003) , 'xtickLabel' , ageGroup(4:10));
    title('2003 National HIV prevalence by sex and age')
    ylim([0 30])
    legend({'Males (model)', 'Females (model)', 'DHS 2003 males', 'DHS 2003 females'}, ...
        'Location' , 'northeast')
    grid on;
 
figure()
%2007
for g = 1 : gender
    for a = 4: 10
        %a = aVec{aInd};
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hiv2007(a-3 , g) = (sum(popVec((2007 - startYear)* stepsPerYear,hivInds),2) ... 
            /sum(popVec((2007 - startYear)* stepsPerYear,ageInds),2))*100;
    end
end
    plot(1 : size(hiv2007 , 1) , hiv2007(: , :) , '-', 1 : size(hiv2007 , 1), Kais2007_Ken(:,:).* 100, '--o');
    hold all;
    xlabel('Age Group'); ylabel('HIV Prevalence')
    set(gca , 'xtick' , 1 : length(hiv2007) , 'xtickLabel' , ageGroup(4:10));
    title('2007 National HIV prevalence by sex and age')
    ylim([0 30])
    legend({'Males (model)', 'Females (model)', 'Kais 2007 males', 'Kais 2007 females'}, ...
        'Location' , 'northeast')
    grid on;
    hold off;
%2009
figure()
for g = 1 : gender
    for a = 4: 10
        %a = aVec{aInd};
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hiv2009(a-3 , g) = (sum(popVec((2009 - startYear)* stepsPerYear,hivInds),2) ... 
            /sum(popVec((2009 - startYear)* stepsPerYear,ageInds),2))*100;
    end
end
 plot(1 : size(hiv2009 , 1) , hiv2009(: , :) , '-', 1 : size(hiv2009 , 1), DHS2009_Ken(:,:), '--o');
    hold all;
    xlabel('Age Group'); ylabel('HIV Prevalence')
    set(gca , 'xtick' , 1 : length(hiv2009) , 'xtickLabel' , ageGroup(4:10));
    title('2009 National HIV prevalence by sex and age')
    ylim([0 30])
    legend({'Males (model)', 'Females (model)', 'DHS 2008/09 males', 'DHS 2008/09 females'}, ...
        'Location' , 'northeast')
    grid on; 
    hold off;
% 2012
figure()
for g = 1 : gender
    for a = 4: 10
        %a = aVec{aInd};
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hiv2012(a-3 , g) = (sum(popVec((2012 - startYear)* stepsPerYear,hivInds),2) ... 
            /sum(popVec((2012 - startYear)* stepsPerYear,ageInds),2))*100;
    end
end

    plot(1 : size(hiv2012, 1) , hiv2012(: , :) , '-', 1 : size(hiv2012, 1), Kais2012_Ken(:,:), '--o');
    hold all;
    xlabel('Age Group'); ylabel('HIV Prevalence')
    set(gca , 'xtick' , 1 : length(hiv2012) , 'xtickLabel' , ageGroup(4:10));
    title('2012 National HIV prevalence by sex and age')
    ylim([0 30])
    legend({'Males (model)', 'Females (model)', 'Kais 2012 males', 'Kais 2012 females'}, ...
        'Location' , 'northeast')
    grid on;
     
%%
sheet = ['HIV_by_age'];
cols1 = {toNowName};
cols2 = {'Age group', 'Male (2003)', 'Female (2003)', ...
    'Male (2007)', 'Female (2007)', ...
    'Male (2009)', 'Female (2009)',...
    'Male (2012)', 'Female (2012)'};
xlswrite(filename, cols1, sheet, 'A64')
xlswrite(filename, cols2, sheet, 'A65')
xlswrite(filename, ageGroup(4:10)', sheet, 'A66')
xlswrite(filename, [hiv2003], sheet, 'B66')
xlswrite(filename, [hiv2007], sheet, 'D66')
xlswrite(filename, [hiv2009], sheet, 'F66')
xlswrite(filename, [hiv2012], sheet, 'H66')

%% total mortality by age over time 
bkgdDeathsF = zeros(length(tVec(1:stepsPerYear:end-1)), length(age));
popTotAgeF = bkgdDeathsF;
hivDeathsF = bkgdDeathsF;
ageGroup = {'0-4','5-9' ,'10-14' , '15-19' , '20-24' , '25-29' ,...
     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
     '60-64' , '65-69' , '70-74' , '75-79'};
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear));
% popPropYrs = zeros(length(years),5);
% popPropYrs_obs = zeros(7,8);
% popPropYrs_obs = xlsread(file , 'Population' , 'H178:O184');
% ageVec = {[1:2], [3:4] , [5:6] , [7:8] , [9:10] , [11:12], [13:16]};
%women
    for a = 1 : age
        %a = ageVec{aInd};
        popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
%         popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : endpoints , 1 : intervens , 2  , 1 : age , 1 : risk));
        %hivDeathInds = toInd(allcomb(1, a));
        hivDeathsF(:, a) = annlz(sum(hivDeaths(1:end-1, 2, a), 2)) ;
        bkgdDeathsF(:, a) = annlz(sum(deaths(1:end-1, popAge), 2));  
        popTotAgeF(:, a) = annlz(sum(popVec(1:end-1 , popAge), 2)) ./ stepsPerYear;  
           
    end

totalDeathsF = ((bkgdDeathsF + hivDeathsF)./ popTotAgeF) ;
%%
figure()
%plot(tVec, squeeze(totalDeaths(:, 1, :)))
%plot(tVec, bsxfun(@rdivide ,squeeze(bkgdDeathsF(:, 1: 10)), popTotAgeF(:, 1: 10)) .* 100)
plot(tVec(1:stepsPerYear:end-1), totalDeathsF(:, 1:10))
xlim([1930 2020])
ylim([0 0.12])
ylabel('Mortality rates')
legend(ageGroup(1:10), 'Location', 'NorthEastOutside')
title('All-cause mortality')


%% HIV mortality by age
ageGroup = {'0-4','5-9' ,'10-14' , '15-19' , '20-24' , '25-29' ,...
     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
     '60-64' , '65-69' , '70-74' , '75-79'};
figure()
% subplot(2, 1, 1)
plot(tVec , squeeze(hivDeaths(:, 2, 4:10)) )
xlim([1980 2020])
xlabel('Year'); ylabel('Deaths'); title('Female HIV-associated Deaths, 15-49')
legend('15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , 'Location' , 'NorthEastOutside')
% subplot(2,1, 2)
plot(tVec(1:stepsPerYear:end-1) , (hivDeathsF(:, 4:10))./popTotAgeF(:, 4:10) .* 100)
xlim([1980 2020])
xlabel('Year'); ylabel('Deaths relative to population (%)'); title('Female HIV-associated Deaths, 15-49')
legend(ageGroup(4:10), 'Location' , 'NorthEastOutside')

%%
figure()
subplot(2, 1, 1)
plot(tVec , squeeze(hivDeaths(:,2, 11:16)) )
xlim([1980 2020])
xlabel('Year'); ylabel('Deaths'); title('Female HIV-associated Deaths, 50+')
legend('50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')
subplot(2,1, 2)
plot(tVec , bsxfun(@rdivide , squeeze(hivDeaths(:,2,11:16)) , sum(popVec , 2)) .* 100)
xlim([1980 2020])
xlabel('Year'); ylabel('Deaths relative to population (%)'); title('Female HIV-associated Deaths, 50+')
legend('50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')

%% Relative HIV prevalence for untreated and on ART
figure()
artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , 1 : gender , 4 : 10 , 1 : risk));
artPop = sum(popVec(: , artInds) , 2);
hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 : gender , 4 : 10 , 1 : risk));
hivPop = sum(popVec(: , hivInds) , 2);
popInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 4 : 10 , 1 : risk));
popTot = sum(popVec(: , popInds) , 2);
hiv_art = [100 * hivPop ./ popTot , 100 * artPop ./ popTot];
area(tVec , hiv_art); %art ./ sum(popVec , 2) , tVec , hiv ./ sum(popVec , 2))
xlabel('Year')
ylabel('Proportion of Population ages 15-49 (%)')
title('Relative HIV Prevalence')
legend('Untreated', 'On ART' , 'Location' , 'NorthWest')
xlim([1980 2020]);


%% Proportion of total HIV+ population on ART (CD4-eligible and ineligible)
figure()

artUNAIDS = [31 38 42 45 51 60 67 73 69];
yrsUNAIDS = 2010:2018;
artWB = [10 14 19 27 33 41 46 49 55];
yrsWB = 2006:2014;
artMOH = [73, 67, 80];
yrsMOH = [2013 2015 2017]

artIndsF = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , 2 , 3 : age , 1 : risk));
hivAllIndsF = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 2 , 3 : age , 1 : risk));
artIndsM = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , 1 , 3 : age , 1 : risk));
hivAllIndsM = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 , 3 : age , 1 : risk));

plot(tVec , 100 * sum(popVec(: , artIndsF) , 2) ./ (sum(popVec(: , hivAllIndsF) , 2) + sum(popVec(: , artIndsF) , 2)) , ...
    tVec , 100 * sum(popVec(: , artIndsM) , 2) ./ (sum(popVec(: , hivAllIndsM) , 2) + sum(popVec(: , artIndsM) , 2)) , ...
    (artYr + 1) , maxRateF .* 100 , 'o' , ...
    (artYr + 1) , maxRateM .* 100 , 'o', ...
    yrsUNAIDS , artUNAIDS , 'v' , ...
    2018 , 78.8 , '+' , ...
    2012 , 38.3 , '+' , ...    
    yrsMOH , artMOH , 'v')
xlabel('Year')
xlim([2000 2020])
ylabel('Proportion of HIV Population')
title('Viral suppression coverage')
legend('Model: Females' , ...
    'Model: Males' , ...
    'Observed Kenya: Females' , ...
    'Observed Kenya: Males', ...
    'UNAIDS (Nyanza, ART)' , ...
    'PHIA (Nyanza)' , ...
    'KAIS (Nyanza)' , ...
    'MOH reports (Nyanza, ART)',...
    'Location', 'Northwest')

%% Proportion of HIV+ population on ART, by age
ageGroup = {'0-4','5-9' ,'10-14' , '15-19' , '20-24' , '25-29' ,...
     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
     '60-64' , '65-69' , '70-74' , '75-79'};
gVec = {'men', 'women'};

for g = 1 : 2
    figure()
    
artInds_tot = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , 1:age , 1 : risk));
artPop_tot = sum(popVec(: , artInds_tot) , 2);

hivInds_tot = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , g , 1:age, 1 : risk));
hivPop_tot = sum(popVec(: , hivInds_tot) , 2);
for a = 1:4
    artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , a , 1 : risk));
    artPop = sum(popVec(: , artInds) , 2);
    
    hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , g , a , 1 : risk));
    hivPop = sum(popVec(: , hivInds) , 2);
    
    subplot(2, 3, 1)
    plot(tVec, 100 * hivPop ./ hivPop_tot)
     xlim([2000 2020])
     ylabel('Proportion of HIV Population')
     title(['Age distribution of HIV+', gVec(g)])
    hold on
    subplot(2, 3, 4)
    plot(tVec , 100 * artPop ./ (artPop_tot))
    ylabel('Proportion of VS Population')
    hold on
    legend(ageGroup(1:4) , 'Location' , 'NorthWest')
    title(['Age distribution of VS', gVec(g)])
    xlim([2000 2020]);
end
for a = 5:10
    artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , a , 1 : risk));
    artPop = sum(popVec(: , artInds) , 2);
    
    hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , g , a , 1 : risk));
    hivPop = sum(popVec(: , hivInds) , 2);
    subplot(2, 3, 2)
    plot(tVec, 100 * hivPop ./ hivPop_tot)
     xlim([2000 2020]);
    ylabel('Proportion of HIV Population')
    hold on
    subplot(2, 3, 5)
    plot(tVec , 100 * artPop ./ (artPop_tot))
    ylabel('Proportion of VS Population')
    hold on
    legend(ageGroup(5:10) , 'Location' , 'NorthWest')
    xlim([2000 2020]);
end
for a = 11:16
    artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , a , 1 : risk));
    artPop = sum(popVec(: , artInds) , 2);
    
    hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , g , a , 1 : risk));
    hivPop = sum(popVec(: , hivInds) , 2);
    subplot(2, 3, 3)
    plot(tVec, 100 * hivPop ./ hivPop_tot)
     xlim([2000 2020]);
    ylabel('Proportion of HIV Population')
    hold on
    subplot(2, 3, 6)
    plot(tVec , 100 * artPop ./ (artPop_tot))
    hold on
    ylabel('Proportion of VS Population')
    legend(ageGroup(11:16) , 'Location' , 'NorthWest')
    xlim([2000 2020]);
    hold on
end
end

%% On ART by age
aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80}; %{10:15,16:25,26:35,36:50,51:75};
ageGroup = {'0-4','5-9' ,'10-14' , '15-19' , '20-24' , '25-29' ,...
     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
     '60-64' , '65-69' , '70-74' , '75-79'};
aMatrix = zeros(1 , age); %length(aVec));
for a = 1 : age %length(aVec)
    %a = aVec{aInd};
    artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    artPop = sum(popVec(end , artInds) , 2); %end-605
    hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hivPop = sum(popVec(end , hivInds) , 2);
    hiv_art = [100 * artPop ./ hivPop];
    aMatrix(1 , a) = hiv_art;
end

phia2018_art=zeros(2, age);
phia2018_art(1, :) = 1: age;
phia2018_art(2, :) = [53.5  53.5 53.5 59 59 70.8 70.8 77 77 83 83 85.1 85.1 NaN NaN NaN];
% PHIA estimates = proportion of all PLHIV virally suppressed

figure;
hold all;
plot([1:age] , aMatrix(1,:) , '-', phia2018_art(1,:), phia2018_art(2, :), 'o' )
hold all;
ylabel('Percent virally suppressed');
title('Proportion of VS among HIV+ women in 2020 by age');
ylim([0 100])
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
legend('Model' , 'PHIA 2018', 'Location', 'Southeast')
grid on;
%legend('Without ART dropout' , 'With ART dropout');
% legend('Without ART dropout' , 'With ART dropout: 6.19%' , 'With ART dropout: 11.8%' , 'With ART dropout: 11.8%, HIV mort on ART');

%% HPV prevalence over time by HIV status and gender
genders = {'Male' , 'Female'};
figure()
for g = 1 : gender
    % General
    hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , g , 4 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , g , 4 : age , 1 : risk))]);
    hpvPop = sum(popVec(: , hpvInds) , 2);
    popTot = popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , 4 : age , 1 : risk)));
    genPopHPV = 100 * hpvPop ./ sum(popTot , 2);
    % HIV+
    hpvHivInds = unique([toInd(allcomb(3 : 7 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , g , 4 : age , 1 : risk)); toInd(allcomb(3 : 7 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , g , 4 : age , 1 : risk))]);
    hpvHivPop = sum(popVec(: , hpvHivInds) , 2);
    popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , 4 : age , 1 : risk)));
    hivPosHPV = 100 * hpvHivPop ./ sum(popHivTot , 2);
    %ART
    hpvArtInds = unique([toInd(allcomb(8 , 6 , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , g , 4 : age , 1 : risk)); toInd(allcomb(8 , 6 , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , g , 4 : age , 1 : risk))]);
    hpvArtPop = sum(popVec(: , hpvArtInds) , 2);
    popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , 4 : age , 1 : risk)));
    artHPV =  100 * hpvArtPop ./ sum(popArtTot , 2);
    %HIV-
    hpvHivNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , g , 4 : age , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , g , 4 : age , 1 : risk))]);
    hpvHivNegPop = sum(popVec(: , hpvHivNegInds) , 2);
    popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , 4 : age , 1 : risk)));
    hivNegHPV = 100 * hpvHivNegPop ./ sum(popHivNegTot , 2);

    subplot(2 , 1 , g)
    plot(tVec , genPopHPV)
    hold on
    plot(tVec , hivNegHPV)
    plot(tVec , hivPosHPV)
    plot(tVec , artHPV)
     xlim([1925 2020])
    ylim([0 80])
    xlabel('Year'); ylabel('Prevalence (%)'); title(['Kenya ', genders{g} , ' HPV Prevalence (ages 15+)'])
    legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Location' , 'NorthWest')


% sheet = ['hpvPrev'];
% cols1 = {toNowName, "2.2xC3toCC"};
% cols2 = {[genders{g},' Gen Pop'], [genders{g},' HIV-neg'], [genders{g},' HIV-pos'], [genders{g},' ART']};
% xlswrite(filename, cols1, sheet, 'R')
% cell = ['R', 'V'];

% xlswrite(filename, cols2, sheet, [cell(g) +'2'])
% xlswrite(filename, [genPopHPV(1 : stepsPerYear * 5 : end), ...
%     hivNegHPV(1 : stepsPerYear * 5 : end), hivPosHPV(1 : stepsPerYear * 5 : end), ...
%     artHPV(1 : stepsPerYear * 5 : end)], sheet, [cell(g) +'3'])
end
%% HPV Prevalence by age in 2005 vs. Yamada and Luchter data
ageGroup = {'17 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};

yr = 2005;
year = {yr};
hpv2005 = zeros(9 , 1);
hpvHIV2005 = hpv2005;
hpvNeg2005 = hpv2005;

aVec = {18:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
%for aInd = 1 : 13
for a = 4 : 12
    %a = aVec{aInd};
    hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpv2005(a - 3 , 1) = sum(popVec((yr - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;
    
    % HIV+
    hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpvHIV2005(a - 3 , 1) = sum(popVec((yr - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;
    
    % HIV-
    hpvInds = unique([toInd(allcomb(1 : 2 , 1  , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : 2 , 1  , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpvNeg2005(a - 3 , 1) = sum(popVec((yr - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;
end

% DeVuyst, 2012 (data collected 2009)
% HIV+ Nairobi
hpvHivObs(: , 1) = [NaN
NaN
0.613
0.535
0.529
0.495
0.495
NaN
NaN];

%HIV+ family planning/gyn health center, any HPV prevalence, but >50% types tested for
% were HR: 16,18,31,33,35,39,45,51,52,56,58,59,67,68,82,(6,11,13,44),(26, 30, 53, 66,70)
hpvHivObs(: , 2) = [1
0.58
0.56
0.43
0.46
0.43
0.43
0.43
0.43];
%HIV- family planning/gyn health center
hpvHivObs(: , 3) = [0.44
0.25
0.12
0.15
0.10
0.10
0.1
0.1
0.1];

hpvHivObs = hpvHivObs * 100;
% hpvNegObs = hpvNegObs * 100;
figure()
% plot(1 : length(hpv2002) , hpv2002 , 'o-')
% hold all;
plot(1 : length(hpvHIV2005) , hpvHIV2005 , 'o-');
hold all;
plot(1 : length(hpvNeg2005) , hpvNeg2005 , 'o-')
hold all;
set(gca , 'xtickLabel' , ageGroup);

plot(1 : length(hpvHivObs) , hpvHivObs(: ,1) , '*--')
plot(1 : length(hpvHivObs) , hpvHivObs(: ,2) , '+--')
plot(1 : length(hpvHivObs) , hpvHivObs(: ,3) , '+--')

set(gca , 'xtick' , 1 : length(hpvHivObs) , 'xtickLabel' , ageGroup);
legend('Model HIV-pos' , 'Model HIV-neg' ,...
    'Obs HIV-pos: DeVuyst (CC screening 2009)' ,...
    'Obs HIV-pos: Yamada (Nairobi, clinic)','Obs HIV-neg: Yamada (Nairobi, clinic)')
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)')
title(['Kenya HPV prevalence in women in', year ])
ylim([0 100])
%%
sheet = ['HPV_by_age_2005'];
cols1 = {toNowName};
cols2 = {'HIV+, higherHIVmort 2020', 'HIV-, higherHIVmort 2020'} %, 'DeVuyst HIV+ 2009', 'Yamada HIV+', 'Yamada HIV-'};
xlswrite(filename, cols1, sheet, 'V1')
xlswrite(filename, cols2, sheet, 'V2')
%xlswrite(filename, ageGroup(1:9)', sheet, 'A3')
xlswrite(filename, [hpvHIV2005, hpvNeg2005], sheet, 'V3')

%% Age specific HPV prevalence data 
ageGroup = {'17 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
yr = 2000;
hpv2000 = zeros(9 , 1);
hpvHIV2000 = hpv2000;
hpvNeg2000 = hpv2000;

aVec = {18:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
%for aInd = 1 : 13
for a = 4 : 12
    %a = aVec{aInd};
    hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpv2000(a - 3 , 1) = sum(popVec((yr - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;
%  
end

% DeVuyst, 2003 (data collected 1998-2000)
% Gen pop, prevalence
hpvHivObs2(: , 1) = [NaN
NaN
0.414
0.277
0.409
0.186
0.243
0.084
NaN];

% HIV+ 
hpvHivObs2(: , 2) = [NaN
NaN
0.613
0.535
0.529
0.495
0.495
NaN
NaN];

hpvHivObs2 = hpvHivObs2 * 100;
% hpvNegObs = hpvNegObs * 100;
figure()
plot(1 : length(hpv2000) , hpv2000 , 'o-')
hold all;
set(gca , 'xtickLabel' , ageGroup);

plot(1 : length(hpvHivObs2) , hpvHivObs2(: ,1) , '*--')
%plot(1 : length(hpvHivObs) , hpvHivObs(: ,2) , '+--')

set(gca , 'xtick' , 1 : length(hpvHivObs2) , 'xtickLabel' , ageGroup);
legend('Model General pop' ,...
    'Obs gen pop: DeVuyst (Nairobi, 2000)')
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)')
title('HPV prevalence among women in 2000')
ylim([0 100])

sheet = ['HPV_by_age_2000'];
cols1 = {toNowName,'HPV prevalence among women in 2000'};
%cols2 = {'Age', 'Model Gen Pop'}; %, 'DeVuyst gen pop', 'DeVuyst HIV+'};
xlswrite(filename, cols1, sheet, 'I1')
%xlswrite(filename, cols2, sheet, 'I2')
xlswrite(filename, ageGroup(1:9)', sheet, 'I3')
xlswrite(filename, [hpv2000], sheet, 'J3')


%% HPV prevalence by age and HIV status in 2008 vs. Mbulawa data
yearPrev = 2008;
ageVec = {[4:5],[6:7],[8:9],[10:13]};
ageGroup2 = {'18-24' , '25-34' , '35-44' , '45-66' };
hpv_hivM2008 = zeros(length(ageVec) , 2);

for aV = 1 : length(ageVec)
    a = ageVec{aV};
    hpvHivPosInd = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 , [1 : 2 , 7] , ...
        1 , 1 : intervens , 1 , a  , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
    popHivInd = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , ...
        1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    hpv_hivM2008(aV , 1) = (sum(popVec((yearPrev - startYear) * stepsPerYear +1 , hpvHivPosInd)) ...
        ./ sum(popVec((yearPrev - startYear) * stepsPerYear +1 , popHivInd))) * 100;
    
    hpvHivNegInd = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 , [1 : 2 , 7] , 1 , ...
        1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
    popNegInd = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , ...
        1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    hpv_hivM2008(aV , 2) = (sum(popVec((yearPrev - startYear) * stepsPerYear +1 , hpvHivNegInd)) ...
        ./ sum(popVec((yearPrev - startYear) * stepsPerYear +1 , popNegInd))) * 100;
end

smithHPV = zeros(4,1);
smithHPV(:, 1) = [0.353 NaN NaN NaN ];

figure()
plot([1 : length(ageVec)] , hpv_hivM2008(: , 1)' , 'b-');
hold all;
plot([1 : length(ageVec)] , hpv_hivM2008(: , 2)' , 'r-')
hold all;
plot([1 : length(ageVec)] , hpv_hivM2008_dObs(: , 2)' .* 100 , 'bo');
hold all;
plot([1 : length(ageVec)] , hpv_hivMNeg2008_dObs(: , 2)' .* 100 , 'ro', ...
    [1 : length(ageVec)] , smithHPV(: , 1)' .* 100, '*');
set(gca , 'xtick' , [1 : length(ageVec)] , 'xtickLabel' , ageGroup2);
legend('HIV+ Model', 'HIV- Model', 'HIV+ Mbulawa 2008 (Cape Town)', ...
    'HIV- Mbulawa 2008 (Cape Town)', 'HIV- Smith 2010 (Kisumu)');
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)'); 
title('HPV prevalence in 2008 men by HIV status');
ylim([0 100]);


%%  CIN Prevalence by HIV stat among women
yr = 2020;
year = {yr};
cinHiv_ccScreen = [0.18	0.105 0.089] .* 100; % Observed, HIV+, Nairobi, Memiah, 2013
cinHiv_ccScreen2 = [33.1 10.5 12.8]; %Observed, HIV clinic, De Vuyst 2013
cinHIV_FSW = [19.1 NaN 13.1]; %Observed, HIV+ FSW, Nairobi, Sweet 2020
cinNeg_FSW = [8.0 NaN 1.5]; % Observed, HIV- FSW, Nairobi, Sweet 2020
cinHIV_FSW2 = [17.9 NaN 11.9]; %Observed, HIV+ FSW, Nairobi, Njagi 2013
cinNeg_FSW2 = [9.0 NaN 1.88]; %Observed, HIV- FSW, Nairobi, Njagi 2013
cinHiv2010 = zeros(3 , 1);
cinHivNeg2010 = cinHiv2010;
cinNoART2010 = cinHiv2010;
cinART2010 = cinHiv2010;

for cin = 3 : 5
        % HIV+, all 
        cinPosInds = unique([toInd(allcomb(3 : 8 , 1 : viral , cin , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 5:12 , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , cin , 1 , 1 : intervens , 2 , 5:12 , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 5:12 , 1 : risk));
        cinHiv2010(cin-2) = (sum(popVec((yr - startYear) * stepsPerYear , cinPosInds)))...
        ./ sum(popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;
        
        % HIV+, no ART
        cinPosInds = unique([toInd(allcomb(3 : 7 , 1 : viral , cin , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 5:12 , 1 : risk)); toInd(allcomb(3 : 7 , 1 : viral , ...
        [1 : 5 , 7] , cin , 1 , 1 : intervens , 2 , 5:12 , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 5:12 , 1 : risk));
        cinNoART2010(cin-2) = (sum(popVec((yr - startYear) * stepsPerYear , cinPosInds)))...
        ./ sum(popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;
    
        % ART
        cinArtInds = unique([toInd(allcomb(8 , 1 : viral , cin , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 5:12 , 1 : risk)); toInd(allcomb(8 , 1 : viral , ...
        [1 : 5 , 7] , cin , 1 , 1 : intervens , 2 , 5:12 , 1 : risk))]);
        ageInds = toInd(allcomb(8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 5:12 , 1 : risk));
        cinART2010(cin-2) = (sum(popVec((yr - startYear) * stepsPerYear , cinArtInds)))...
        ./ sum(popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;
    
        % HIV-
        cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , cin , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 5:12 , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , cin , 1 , 1 : intervens , 2 , 5:12 , 1 : risk))]);
        ageInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 5:12 , 1 : risk));
        cinHivNeg2010(cin-2) = (sum(popVec((yr - startYear) * stepsPerYear , cinNegInds)))...
        ./ sum(popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;     
end    

figure()
cinGroup = {'CIN 1' , 'CIN 2' , 'CIN 3'};
plot(1 : length(cinHiv2010) , cinHiv2010 , '-' , 1 : length(cinHiv2010), cinHivNeg2010, '-' ,...
    1 : length(cinHiv2010), cinART2010, '--', 1 : length(cinHiv2010), cinNoART2010, '--')
hold on
p = plot(1: length(cinGroup) , cinHiv_ccScreen , '+',  1: length(cinGroup),cinHiv_ccScreen2, '+',...
    1 : length(cinGroup), cinHIV_FSW, 'g+', 1 : length(cinGroup), cinNeg_FSW, 'go',...
    1 : length(cinGroup), cinHIV_FSW2, 'm+', 1 : length(cinGroup), cinNeg_FSW2, 'mo');
p(1).MarkerSize = 10;
p(2).MarkerSize = 8;
ylabel('Prevalence (%)')
set(gca , 'XTick', 1:3, 'xtickLabel' , cinGroup);
xlabel('CIN Stage')
legend('Model HIV+' , 'Model HIV-', 'Model ART', 'Model, HIV+, no ART', ...
    'HIV+ screening pop (Memiah)', 'HIV+ screening pop (DeVuyst)', ...
    'HIV+ FSW (Sweet)', 'HIV- FSW (Sweet)', 'HIV+ FSW (Njagi)', 'HIV- FSW (Njagi)')
title(['Kenya CIN Prevalence by HIV status in ', year])
text(2.2, -3.25, 'Note: All observed data are from Nairobi')
%annotation('textbox', [0, 0, 0.5, 0.1], 'string', 'All observed data are from Nairobi')

%%
sheet = ['CIN_by_HIV_2010'];
cols1 = {[toNowName], ['CIN Prevalence by HIV status in ', year]};
cols2 = {'HIV+, higherHIVmort 2020', 'HIV-, higherHIVmort 2020'} %, 'HIV+ screening pop (Memiah)',...
    %'HIV+ screening pop (DeVuyst)', ...
    %'HIV+ FSW (Sweet)', 'HIV- FSW (Sweet)', 'HIV+ FSW (Njagi)', 'HIV- FSW (Njagi)'};
xlswrite(filename, cols1, sheet, 'X1')
xlswrite(filename, cols2, sheet, 'X2')
%xlswrite(filename, cinGroup', sheet, 'K3')
xlswrite(filename, [cinHiv2010, cinHivNeg2010], sheet, 'X3')
% xlswrite(filename, [cinHiv2010, cinHivNeg2010, cinHiv_ccScreen',cinHiv_ccScreen2',...
%     cinHIV_FSW', cinNeg_FSW', cinHIV_FSW2', cinNeg_FSW2'], sheet, 'B9')
 
%% CIN2/3 prevalence for All HR HPV types combined by age in 2000 vs. deVuyst 2003 data
cin1Pos2000 = zeros(10 , 1);
cin3Pos2000 = cin1Pos2000;
ageGroup = {'17-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};
year = 2000;
%aVec = {18:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
for a = 4 : 13  %note, age group 4 is 17-19 in the data
    %a = aVec{aInd};
    % CIN 1 prevalence 
    cin1Inds = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 3 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cin1Pos2000(a - 3) = (sum(popVec((year - startYear) * stepsPerYear , cin1Inds)))...
        ./ sum(popVec((year - startYear) * stepsPerYear , ageInds)) * 100;
    % CIN 2/3 prevalence 
    cin3Inds = unique([toInd(allcomb(1 : disease , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cin3Pos2000(a - 3) = (sum(popVec((year - startYear) * stepsPerYear , cin3Inds)))...
        ./ sum(popVec((year - startYear) * stepsPerYear , ageInds)) * 100;
end

% DeVuyst 2003
 % CIN 1
cinPosAct(: , 1) = [NaN
NaN
0.145
0.097
0.0593
0.0142
0.0435
0.0841
NaN
NaN]; 

 %CIN 2/3
cinPosAct(: , 2) = [NaN
NaN
0.0526
0.0842
0.0880
0.0315
0.0427
0.0029
NaN
NaN]; 


figure();
cinPosAct = cinPosAct .* 100; % convert to %

plot(1 : length(cin1Pos2000) , cin1Pos2000 ,'o-')
hold all 
plot(1 : length(cin3Pos2000) , cin3Pos2000 ,'o-')
hold all 
set(gca , 'xtickLabel' , ageGroup);

plot(1 : length(cinPosAct) , cinPosAct(: , 1), '+--')
plot(1 : length(cinPosAct) , cinPosAct(: , 2), '+--')
legend('Model CIN 1', 'Model CIN 2/3' , 'DeVuyst CIN 1', 'DeVuyst CIN 2/3')
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('Age specific CIN prevalence among all women in 2000')
ylim([0 30])
%% 

sheet = ['CIN_by_age_2000'];
cols1 = {toNowName, 'Age specific CIN prevalence among all women in 2000'};
cols2 = {'Age', 'Model CIN 1', 'Model CIN 2/3'}; % , 'DeVuyst CIN 1', 'DeVuyst CIN 2/3'};
xlswrite(filename, cols1, sheet, 'O1')
xlswrite(filename, cols2, sheet, 'O2')
xlswrite(filename, ageGroup', sheet, 'O3')
xlswrite(filename, [cin1Pos2000, cin3Pos2000], sheet, 'P3')
% xlswrite(filename, [cin1Pos2000, cin3Pos2000, cinPosAct], sheet, 'H3')



%% HPV type distribution by state over time (coinfections grouped as 9v-type HPV)
% HPV infected
hpvInds_vax = toInd(allcomb(1 : disease , 1 : viral , 2 , [1 : 2 , 7] , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
hpvInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 , 7] , 2 , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
hpvInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 2 , [1 : 2 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 , 7] , 2 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
hpv_vax = sum(popVec(: , hpvInds_vax) , 2)...
    ./ sum(popVec(: , hpvInds_tot) , 2) * 100;
hpv_nonVax = sum(popVec(: , hpvInds_nonVax) , 2)...
    ./ sum(popVec(: , hpvInds_tot) , 2) * 100;

figure;
subplot(2,3,1)
plot(tVec , hpv_vax , 'k')
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 46.82 , 'ko')
hold all;
plot(tVec , hpv_nonVax);
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 53.18 , 'o');
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('HPV');
ylim([0 100]);
xlim([2000 2015]);
%legend('9v-type HPV' , 'Observed 2011: 9v' , 'Non-9v-type HPV' , 'Observed 2011: non-9v');

% CIN1
cinInds_vax = toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cinInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 2 , 7] , 3 , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cinInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
cin_vax = sum(popVec(: , cinInds_vax) , 2)...
    ./ sum(popVec(: , cinInds_tot) , 2) * 100;
cin_nonVax = sum(popVec(: , cinInds_nonVax) , 2)...
    ./ sum(popVec(: , cinInds_tot) , 2) * 100;

subplot(2,3,2)
plot(tVec , cin_vax , 'k')
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 51.92 , 'ko')
hold all;
plot(tVec , cin_nonVax);
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 48.08 , 'o')
ylim([0 100]);
xlim([2000 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('CIN1')
%legend('9v-type HPV' , 'Observed 2011: 9v' , 'Non-9v-type HPV' , 'Obsersved 2011: non-9v');

% CIN2
cinInds_vax = toInd(allcomb(1 : disease , 1 : viral , 4 , [1 : 4 , 7] , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cinInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 3 , 7] , 4 , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cinInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 4 , [1 : 4 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 3 , 7] , 4 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
cin_vax = sum(popVec(: , cinInds_vax) , 2)...
    ./ sum(popVec(: , cinInds_tot) , 2) * 100;
cin_nonVax = sum(popVec(: , cinInds_nonVax) , 2)...
    ./ sum(popVec(: , cinInds_tot) , 2) * 100;

subplot(2,3,3)
plot(tVec , cin_vax , 'k')
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 62.81 , 'ko')
hold all;
plot(tVec , cin_nonVax);
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 37.19 , 'o')
ylim([0 100]);
xlim([2000 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('CIN2')
%legend('9v-type HPV' , 'Observed 2011: 9v' , 'Non-9v-type HPV' , 'Observed 2011: non-9v');

% CIN3
cinInds_vax = toInd(allcomb(1 : disease , 1 : viral , 5 , [1 : 5 , 7] , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cinInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 4 , 7] , 5 , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cinInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
cin_vax = sum(popVec(: , cinInds_vax) , 2)...
    ./ sum(popVec(: , cinInds_tot) , 2) * 100;
cin_nonVax = sum(popVec(: , cinInds_nonVax) , 2)...
    ./ sum(popVec(: , cinInds_tot) , 2) * 100;

subplot(2,3,4)
plot(tVec , cin_vax , 'k')
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 73.71 , 'ko')
hold all;
plot(tVec , cin_nonVax);
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 26.29 , 'o')
ylim([0 100]);
xlim([2000 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('CIN3')
%legend('9v-type HPV' , 'Observed 2011: 9v' , 'Non-9v-type HPV' , 'Observed 2011: non-9v');

% CC
ccInds_vax = toInd(allcomb(1 : disease , 1 : viral , 6 , [1 : 6 , 7] , ...
    1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
ccInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , 6 , ...
    1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
ccInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 6 , [1 : 6 , 7] , ...
        1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
cc_vax = sum(popVec(: , ccInds_vax) , 2)...
    ./ sum(popVec(: , ccInds_tot) , 2) * 100;
cc_nonVax = sum(popVec(: , ccInds_nonVax) , 2)...
    ./ sum(popVec(: , ccInds_tot) , 2) * 100;

subplot(2,3,5)
plot(tVec , cc_vax , 'k')
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 85.78 , 'ko')
hold all;
plot(tVec , cc_nonVax);
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 14.22 , 'o')
ylim([0 100]);
xlim([2000 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('Cervical Cancer')
legend('9v-type HPV' , 'Observed 2011: 9v' , 'Non-9v-type HPV' , 'Observed 2011: non-9v',...
    'Location', 'SouthEastOutside');

%% Cervical cancer incidence by age
%ccIncYears = [2017 , 2003 , 1994 , 2012];
ccIncYears = [2009, 2012];
ccAgeRel = zeros(age , length(ccIncYears));
ccAgeNegRel = ccAgeRel;
ccAgePosRel = zeros(age , 5 , length(ccIncYears));
ccArtRel = ccAgeRel;
ccNegPosArt = zeros(age , 3 , length(ccIncYears));
fScale = 10^5;
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 
ccYrs = ((ccIncYears - startYear) * stepsPerYear : ...
    (ccIncYears + 1 - startYear) * stepsPerYear);

%aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
%for aInd = 1 : 16
for a = 1 : age
    %a = aVec{aInd};
    for y = 1 : length(ccIncYears)
        % Year
        yr_start = (ccIncYears(y) - 1 - startYear)  .* stepsPerYear;
        yr_end = (ccIncYears(y) - startYear) .* stepsPerYear - 1;
        
        % Total population
        ageInds = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
            1 : intervens , 2 , a , 1 : risk));
        ccAgeRel(a , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end , ...
            1 : disease , a , :) , 2) , 3) , 4)) ...
            ./ (annlz(sum(popVec(yr_start : yr_end , ageInds) , 2)) ...
            ./ stepsPerYear) * fScale;

        % HIV Negative
        ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
            1 : intervens , 2 , a , 1 : risk));
        ccAgeNegRel(a , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end...
            , 1 : 2 , a , :) , 2) , 3) , 4)) ...
            ./ (annlz(sum(popVec(yr_start : yr_end , ageNegInds) , 2)) ...
            ./ stepsPerYear) * fScale;
        
        % Acute and CD4 > 500
        agePosInds = toInd(allcomb(3 : 4 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
            1 : intervens , 2 , a , 1 : risk));
        ccAgePosRel(a , 1 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end ...
            , 3 : 4 , a , :), 2) , 3) , 4)) ...
            ./ (annlz(sum(popVec(yr_start : yr_end , agePosInds) , 2)) ...
            ./ stepsPerYear) * fScale;
        
        % HIV Positive CD4 500-350 -> CD4 < 200
        for d = 5 : 7
            agePosInds = toInd(allcomb(d , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
                1 : intervens , 2 , a , 1 : risk));
            ccAgePosRel(a , d - 3 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end...
                , d , a , :), 2) , 3) , 4)) ...
                ./ (annlz(sum(popVec(yr_start : yr_end , agePosInds) , 2)) ...
                ./ stepsPerYear) * fScale;
        end

        % All HIV+ no ART
         ageAllPosInds = toInd(allcomb(3 : 7 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
            1 : intervens , 2 , a , 1 : risk));
         ccAgePosRel(a , 5 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end ...
            , 3 : 7 , a , :), 2) , 3) , 4)) ...
            ./ (annlz(sum(popVec(yr_start : yr_end , ageAllPosInds) , 2)) ...
            ./ stepsPerYear) * fScale;
        
        % On ART
        ageArtInds = toInd(allcomb(8 , 6 , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
            1 : intervens , 2 , a , 1 : risk));
        ccArtRel(a , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end ...
            , 8 , a , :) , 2) , 3) , 4)) ...
            ./ (annlz(sum(popVec(yr_start : yr_end , ageArtInds) , 2)) ...
            ./ stepsPerYear) * fScale;

        % Proportion of cervical cancers by HIV/ART status and age
        % Total by age
        ageTotal = annlz(sum(popVec(yr_start : yr_end , ageInds), 2 )) ./ stepsPerYear;

        % HIV-
        ccNegPosArt(a , 1 , y) = (annlz(sum(sum(sum(newCC(yr_start : yr_end , ...
            1 , a , :), 2) , 3) , 4)) ...
            ./ ageTotal) .* fScale;
        % HIV+
        ccNegPosArt(a , 2 , y) = (annlz(sum(sum(sum(newCC(yr_start : yr_end , ...
            3 : 7 , a , :) , 2) , 3) , 4)) ...
            ./ ageTotal) .* fScale;
        % ART
        ccNegPosArt(a , 3 , y) = (annlz(sum(sum(sum(newCC(yr_start : yr_end , ...
            8 , a , :), 2) , 3) , 4)) ...
            ./ ageTotal) .* fScale;
    end
end

% Nairobi and Eldoret cancer registion 2007-2011 (N), 2008-2011 (E), and
% 2004-2008 (N) Nairobi 2004-2008	Nairobi 2007-2011	Eldoret 2008-2011	Eldoret, 2012-2016	Globocan 2012
file = [pwd , '/Config/Kenya_parameters_Feb20.xlsx'];
ken_canreg = zeros(13,5);
ken_canreg = xlsread(file , 'CC prevalence' , 'C46:G58'); % years

for y = 1 
    ccIncYear = ccIncYears(y);
    
    figure()
    % Plot model outputs
    plot(1 : size(ccAgeRel , 1) , ccAgeRel(: , y) , '-ko' , 1 : size(ccAgeNegRel(: , y) , 1) , ...
        ccAgeNegRel(: , y) , '-kp' , 1 : size(ccAgePosRel , 1) , ccAgePosRel(: , 5 , y) , '-k+' , ...
        1 : size(ccArtRel , 1) , ccArtRel(: , y) , '-k^');
    hold on
    % Plot observed data
     plot(4 : age ,  ken_canreg(:, 2)  , '--' , 4 : age ,  ken_canreg(:, 3), '--' );
    xlabel('Year'); ylabel('Incidence per 100,000');
    set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
    title(['Cervical Cancer Incidence by age in ' num2str(ccIncYear)]);
       legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Nairobi 2007-11' , 'Eldoret 2008-11' )
end

for y = 2 
    ccIncYear = ccIncYears(y);
    
    figure()
    % Plot model outputs
    plot(1 : size(ccAgeRel , 1) , ccAgeRel(: , y) , '-ko' , 1 : size(ccAgeNegRel(: , y) , 1) , ...
        ccAgeNegRel(: , y) , '-kp' , 1 : size(ccAgePosRel , 1) , ccAgePosRel(: , 5 , y) , '-k+' , ...
        1 : size(ccArtRel , 1) , ccArtRel(: , y) , '-k^');
    hold on
    % Plot observed data
     plot(4 : age ,  ken_canreg(:, 2)  , '--' , 4 : age ,  ken_canreg(:, 4), '--', 4 : age ,  ken_canreg(:, 5), '--' );
    xlabel('Year'); ylabel('Incidence per 100,000');
    set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
    title(['Cervical Cancer Incidence by age in ' num2str(ccIncYear)]);
       legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Nairobi 2007-2011' , 'Eldoret 2012-2016', 'Globocan 2012' )
end

%%
sheet = ['CC_inc_by_age_2009'];
cols1 = {toNowName, 'Cervical Cancer Incidence by age in 2009'};
cols2 = {'Age', 'General' , 'HIV-' , 'HIV+' , 'ART'}; % , 'Nairobi 2007-11' , 'Eldoret 2008-11'};
xlswrite(filename, cols1, sheet, 'V1')
xlswrite(filename, cols2, sheet, 'V2')
xlswrite(filename, ageGroup', sheet, 'V3')
xlswrite(filename, [ccAgeRel, ccAgeNegRel, ccAgePosRel(:, 5), ccArtRel], sheet, 'W3')
% xlswrite(filename, [ken_canreg], sheet, 'F6')

%% Cervical cancer incidence over time
fScale = 10^5;
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 

% Total population
ageInds = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
    1 : intervens , 2 , 1 : age , 1 : risk));
ccAgeRel = annlz(sum(sum(sum(newCC(1:end-1 , ...
    1 : disease , 1 : age , :) , 2) , 3) , 4)) ...
    ./ (annlz(sum(popVec(1:end-1 , ageInds) , 2)) ...
    ./ stepsPerYear) * fScale;

% HIV Negative
ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
    1 : intervens , 2 , 1 : age , 1 : risk));
ccAgeNegRel = annlz(sum(sum(sum(newCC(1:end-1 , ...
    1 : 2 , 1 : age , :) , 2) , 3) , 4)) ...
    ./ (annlz(sum(popVec(1:end-1 , ageNegInds) , 2)) ...
    ./ stepsPerYear) * fScale;

% All HIV+ no ART
 ageAllPosInds = toInd(allcomb(3 : 7 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
    1 : intervens , 2 , 1 : age , 1 : risk));
 ccAgePosRel = annlz(sum(sum(sum(newCC(1:end-1 ...
    , 3 : 7 , 1 : age , :), 2) , 3) , 4)) ...
    ./ (annlz(sum(popVec(1:end-1 , ageAllPosInds) , 2)) ...
    ./ stepsPerYear) * fScale;

% On ART
ageArtInds = toInd(allcomb(8 , 6 , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
    1 : intervens , 2 , 1 : age , 1 : risk));
ccArtRel = annlz(sum(sum(sum(newCC(1:end-1 ...
    , 8 , 1 : age , :) , 2) , 3) , 4)) ...
    ./ (annlz(sum(popVec(1:end-1 , ageArtInds) , 2)) ...
    ./ stepsPerYear) * fScale;

gbd_cc= zeros(28, 4);
gbd_cc = [
1990	28.34	21.28	41.36
1991	28.30	21.13	41.17
1992	28.45	21.24	40.73
1993	28.54	21.39	40.09
1994	28.83	21.67	39.98
1995	29.01	21.67	39.75
1996	28.90	21.50	39.28
1997	28.92	21.33	38.35
1998	28.70	21.15	37.67
1999	28.14	20.79	36.71
2000	27.92	20.64	36.02
2001	27.14	20.04	35.03
2002	27.28	20.14	34.87
2003	27.26	20.09	34.86
2004	26.90	19.95	34.50
2005	26.58	19.81	34.22
2006	26.28	19.90	34.21
2007	25.87	19.77	33.89
2008	25.69	19.74	34.01
2009	25.69	19.87	34.40
2010	26.10	20.28	35.24
2011	25.54	19.89	34.80
2012	24.98	19.62	34.44
2013	24.34	19.41	33.85
2014	23.76	19.23	33.27
2015	23.40	19.10	32.86
2016	23.11	18.87	32.50
2017	22.72	18.55	31.90
];

globocan_EA = zeros(2, 2);
globocan_EA = [2008 34.5
    2012 42.7];
globocan_Ken = zeros(3, 2);
globocan_Ken = [2008 23.4
    2012 40.03
    2018 32.5]

figure()
% Plot model outputs
plot(tVec(1 : stepsPerYear : end-1) , ccAgeRel)
hold all
plot(tVec(1 : stepsPerYear : end-1) , ccAgeNegRel)
hold all
plot(tVec(1 : stepsPerYear : end-1) , ccAgePosRel)
hold all
plot(tVec(1 : stepsPerYear : end-1) , ccArtRel)
hold all
errorbar(gbd_cc(:, 1) , gbd_cc(:, 2), gbd_cc(:, 2) - gbd_cc(: , 3), ...
    gbd_cc(:, 4)-gbd_cc(: , 2))
hold all
plot(globocan_EA(: ,1) , globocan_EA(: ,2), 'o', ...
    globocan_Ken(:, 1), globocan_Ken(:,2), 'o')
hold all
xlabel('Time'); ylabel('Incidence per 100,000');
title('Cervical Cancer Incidence ');
xlim([1970 2020]);
legend('General' , 'HIV-' , 'HIV+, no ART' , 'HIV+, ART', 'GBD Kenya 2018', 'Globocan E. Africa', 'Globocan Kenya',...
    'Location', 'NorthWest')

%%
sheet = ['CC_inc'];
cols1 = {toNowName, 'Cervical Cancer Incidence'};
cols2 = {'Gen, higherHIVmort' , 'HIV-, higherHIVmort' , 'HIV+ ART-, higherHIVmort' , 'HIV+ ART+, higherHIVmort'}; %,...
   % 'Year', 'GBD Kenya 2018', 'GBD LB', 'GBD UB', 'Year', 'Globocan Kenya'};
xlswrite(filename, cols1, sheet, 'S1')
xlswrite(filename, cols2, sheet, 'S2')
xlswrite(filename, [ccAgeRel', ccAgeNegRel', ccAgePosRel',...
    ccArtRel'], sheet, 'S3')
% xlswrite(filename, [gbd_cc], sheet, 'F68')
% xlswrite(filename, globocan_Ken, sheet, 'J86')


%% CC incidence over time (age-standardized)
inds = {':' , 1 : 2 , 3 : 7 , 8 , 3 : 8}; % HIV state inds
fac = 10 ^ 5;
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g'};
% worldStandard_WP2015 = [325428 (311262/5.0) 295693 287187 291738 299655 272348 ...
%         247167 240167 226750 201603 171975 150562 113118 82266 64484];
worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
        247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 23477 9261 2155];
% sheet = ['CC_inc_standardized'];
% cols1 = {toNowName, 'Age-standardized cervical cancer incidence'};
% cols2 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all' }; %,...
   % 'Year', 'GBD Kenya 2018', 'GBD LB', 'GBD UB', 'Year', 'Globocan Kenya'};
% xlswrite(filename, cols1, sheet, 'N1')
% xlswrite(filename, cols2, sheet, 'N2')
% xlswrite(filename, [tVec(1: stepsPerYear :end-1)'] , sheet, 'N3')

figure();
for i = 1 : length(inds)  
    ccIncRefVec = zeros(length(tVec(1 : stepsPerYear : end-1)),1)';          
    for aInd = 1 : age + 4
        if aInd >= age
            a = age;
        else
            a = aInd;
        end
       
        % General
        allF = toInd(allcomb(1 : disease, 1 : viral, [1 : 5, 7], [1 : 5, 7], ...
            1 , 1 : intervens , 2 , a , 1 : risk));
        % All HIV-negative women
        hivNeg = toInd(allcomb(1 : 2, 1 : viral, [1 : 5 , 7], [1 : 5 , 7], 1, 1 : intervens, 2, a, 1 : risk));
        % HIV-positive women not on ART
        hivNoARTF = toInd(allcomb(3 : 7, 1 : viral, [1 : 5, 7], [1 : 5, 7], ...
            1, 1 : intervens, 2, a, 1 : risk));
        % Women on ART
        artF = toInd(allcomb(8 , 6 , [1 : 5 , 7] , [1 : 5 , 7] , 1 , 1 : intervens , 2 , a , 1 : risk));
        % All HIV-positive women
        hivAllF = toInd(allcomb(3 : 8 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk));
        genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};
       
        % Calculate incidence
        if aInd <= age
            ccIncRef = ...
                (annlz(sum(sum(newCC(1:end-1 , inds{i} , a , :),2),4)) ./ ...
                (annlz(sum(popVec(1:end-1 , genArray{i}) , 2) ./ stepsPerYear))* fac) ...
                .* (worldStandard_WP2015(a));
            if (i == 4) && (a < 3) && (max(annlz(sum(sum(newCC(1:end-1 , inds{i} , a , :),2),4))) == 0.0)
                ccIncRef = zeros(length(tVec(1 : stepsPerYear : end-1)),1)';
            end
        elseif aInd > age
            ccIncRef = ...
                (annlz(sum(sum(newCC(1:end-1 , inds{i} , a , :),2),4)) ./ ...
                (annlz(sum(popVec(1:end-1 , genArray{i}) , 2) ./ stepsPerYear)) * fac);
            ccIncRef = [(ones(1,aInd-a).*ccIncRef(1,1)) , ccIncRef(1,1:end-(aInd-a))];
            ccIncRef = ccIncRef .* worldStandard_WP2015(aInd);
        end
        ccIncRefVec = ccIncRefVec + ccIncRef;
    end
   
    ccInc = ccIncRefVec ./ (sum(worldStandard_WP2015(1:age+4)));
 
    plot(tVec(1 : stepsPerYear : end-1) , ccInc ,'DisplayName' , plotTits1{i});
    hold all;
    
    cell = ['N'; 'O'; 'P'; 'Q'; 'R'];
    xlswrite(filename, ccInc', sheet, [cell(i) +'3'])

end
    errorbar(gbd_cc(:, 1) , gbd_cc(:, 2), gbd_cc(:, 2) - gbd_cc(: , 3), ...
    gbd_cc(:, 4)-gbd_cc(: , 2))
    hold on
    plot(globocan_EA(: ,1) , globocan_EA(: ,2), 'o', globocan_Ken(:, 1), globocan_Ken(:,2), 'o')
    hold on;
    legend('General' , 'HIV-' , 'HIV+, no ART' , 'HIV+, ART', 'HIV all', ...
        'GBD Kenya 2018', 'Globocan E. Africa', 'Globocan Kenya',...
    'Location', 'NorthWest');
    title('Age-standardized Cervical Cancer Incidence ');
    xlabel('Year'); ylabel('Incidence per 100,000');
    grid off;
    xlim([1925 2020]);
    %ylim([0 150]);
    
% xlswrite(filename, [gbd_cc], sheet, 'G68')
% xlswrite(filename, globocan_Ken, sheet, 'K86')

%% HIV incidence by age and gender
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
     '60-64' , '65-69' , '70-74' , '75-79'};
sex = {'males ', 'females '};
 for g = 1 : gender
figure()
for a = 1 : age
    subplot(4,4,a)
   % for r = 1 : risk
%     for g = 1 : 2
        hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens, g , a , 1: risk)); ...
            toInd(allcomb(7 : 9 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens, g , a , 1: risk))];
        hivSus = annlz(sum(popVec(1:end-1 , hivSusInds) , 2)) ./ stepsPerYear;
       % plot(tVec(1 : stepsPerYear : end-1) , ...
        %    annlz(sum(newHiv(1:end-1 , 2 , a), 3)) ./ hivSus * 100000)
         plot(tVec(1 : stepsPerYear : end-1) , ...
              annlz(sum(sum(sum(sum(sum(newHiv(1:end-1 , : , : , : , g , a , :), 2), 3), 4), 6), 7)) ./ hivSus * 100)
        %axis([startYear , endYear , 0 , 100])
        hold all
%     end
  %  end
    xlabel('Year'); ylabel('Rate Per 100'); title(['HIV incidence: ', sex(g), ageGroup(a)])
    xlim([1980 2020]);
    ylim([0 10]);
end

 end
%legend('Male' , 'Female')
% title('HIV incidence')
% 
%% HIV incidence by age and risk among women
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
     '60-64' , '65-69' , '70-74' , '75-79'};
riskgroup = {'low ', 'medium ', 'high'};

figure()
for r = 1 : risk
        hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens, 2 , 4:16 , r)); ...
            toInd(allcomb(7 : 9 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens, 2 , 4:16 , r))];
        hivSus = annlz(sum(popVec(1:end-1 , hivSusInds) , 2)) ./ stepsPerYear;
         plot(tVec(1 : stepsPerYear : end-1) , ...
              annlz(sum(sum(sum(sum(sum(newHiv(1:end-1 , : , : , : , 2 , 4:16 , r), 2), 3), 4), 6), 7)) ./ hivSus * 100)
        hold all
end
xlabel('Year'); ylabel('Rate Per 100'); title('HIV incidence')
xlim([1980 2020]);
ylim([0 10]);
legend('Low', 'Medium', 'High')


%% HPV incidence by gender and HIV status 
figure()
for g = 1 : 2
    hpvSusInds = [toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 , ...
        1 : endpoints , 1 : intervens, g , 7 , 1 : risk)); ...
        toInd(allcomb(1 : 2 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens, g , 7, 1 : risk))];
    hpvSus = annlz(sum(popVec(1:end-1 , hpvSusInds) , 2)) ./ stepsPerYear;
    hpvSusIndsPos = [toInd(allcomb(3 : disease , 1 : viral , 1 : hpvVaxStates , 1 , ...
        1 : endpoints , 1 : intervens, g , 7, 1 : risk)); ...
        toInd(allcomb(3 : disease , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens, g , 7 , 1 : risk))];
    hpvSusPos = annlz(sum(popVec(1:end-1 , hpvSusIndsPos) , 2)) ./ stepsPerYear;
    plot(tVec(1 : stepsPerYear : end-1) , annlz(sum(sum(sum(sum(newHpvVax(1:end-1 , g , 1:2 , 7 , :, :)...
        , 3) , 4) , 5), 6)) ./ hpvSus * 100)
    hold on
    plot(tVec(1 : stepsPerYear : end-1) , annlz(sum(sum(sum(sum(newHpvVax(1:end-1 , g , 3 : disease , 7 , :, :)...
        , 3) , 4) , 5), 6)) ./ hpvSusPos * 100)
    axis([startYear , endYear , 0 , 100])    
    hold on
end
legend('HIV- Male' , 'HIV+ Male',  'HIV- Female', 'HIV+ Female')
ylim([0 40])
xlabel('Year'); ylabel('Rate Per 100'); title('HPV Incidence by gender and HIV status')
%  gender , disease , age , risk , intervens


%% HPV incidence by gender and HIV status 
figure()
for r = 1 : risk
    hpvSusInds = [toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 , ...
        1 : endpoints , 1 : intervens, 2 , 4:16 , r)); ...
        toInd(allcomb(1 : 2 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens, 2 , 4: 16, r))];
    hpvSus = annlz(sum(popVec(1:end-1 , hpvSusInds) , 2)) ./ stepsPerYear;
    hpvSusIndsPos = [toInd(allcomb(3 : disease , 1 : viral , 1 : hpvVaxStates , 1 , ...
        1 : endpoints , 1 : intervens, 2 , 4 : 16, r)); ...
        toInd(allcomb(3 : disease , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens, 2 , 4 : 16 , r))];
    hpvSusPos = annlz(sum(popVec(1:end-1 , hpvSusIndsPos) , 2)) ./ stepsPerYear;
    plot(tVec(1 : stepsPerYear : end-1) , annlz(sum(sum(sum(sum(newHpvVax(1:end-1 , 2 , 1:2 , 4 :16 , r, :)...
        , 3) , 4) , 5), 6)) ./ hpvSus * 100)
    hold on
    plot(tVec(1 : stepsPerYear : end-1) , annlz(sum(sum(sum(sum(newHpvVax(1:end-1 , 2 , 3 : disease , 4:16 , r, :)...
        , 3) , 4) , 5), 6)) ./ hpvSusPos * 100)
    axis([startYear , endYear , 0 , 100])    
    hold on
end
legend('HIV- Low' , 'HIV+ Low',  'HIV- Medium', 'HIV+ Medium', 'HIV- high', 'HIV+ high')
%ylim([0 50])
xlabel('Year'); ylabel('Rate Per 100'); title('HPV Incidence by risk group and HIV status')


%% HIV+ HPV incidence no ART by gender
% % figure()
% % for g = 1 : 2
% %     hpvHivSusInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
% %           1 : endpoints , g , 16 : 50 , 1 : risk)); ...
% %         toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvVaxStates , 10 , ...
% %           1 : endpoints , g , 16 : 50 , 1 : risk))];
% %     hpvHivSus = annlz(sum(popVec(1:end-1 , hpvHivSusInds) , 2)) ./ stepsPerYear;
% %     plot(tVec(1 : stepsPerYear : end-1) , ...
% %         annlz(sum(sum(sum(newHpv(1:end-1 , g , 2 : 6 , 16 : 50 , :) ...
% %         + newImmHpv(1:end-1 , g , 2 : 6 , 16 : 50 , :) ...
% %         , 3), 4), 5)) ./ hpvHivSus * 100)
% %     hold on
% %     xlabel('Year'); ylabel('Rate Per 100'); title(['HPV Incidence in HIV+ , no ART'])
% %     %axis([startYear , endYear , 0 , 100])
% % end
% % legend('Male' , 'Female')
% 
%% HIV+ HPV incidence on ART by gender
% % figure()
% % for g = 1 : 2
% %     hpvHivSusInds = [toInd(allcomb(10 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
% %           1 : endpoints , g , 16 : 50 , 1 : risk)); ...
% %         toInd(allcomb(10 , 1 : viral , 2 : hpvVaxStates , 10 , ...
% %           1 : endpoints , g , 16 : 50 , 1 : risk))];
% %     hpvHivSus = annlz(sum(popVec(1:end-1 , hpvHivSusInds) , 2)) ./ stepsPerYear;
% %     plot(tVec(1 : stepsPerYear : end-1) , ...
% %         annlz(sum(sum(sum(newHpv(1:end-1 , g , 10 , 16 : 50 , :) ...
% %         + newImmHpv(1:end-1 , g , 10 , 16 : 50 , :) ...
% %         , 3), 4), 5)) ./ hpvHivSus * 100)
% %     hold on
% %     xlabel('Year'); ylabel('Rate Per 100'); title([' HPV Incidence in HIV+ , on ART'])
% %     axis([startYear , endYear , 0 , 100])
% % end
% % legend('Male' , 'Female')
% 

%% New infections
% % figure()
% % for g = 1 : 2
% %     for a = 1 : age
% %         subplot(4 , 4 , a)
% %         hpvHivSusInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , g , a , 1 : risk)); ...
% %             toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvVaxStates , 10 , ...
% %             1 : endpoints , g , a , 1 : risk)); ...
% %             toInd(allcomb(10 , 6 , 1 , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , g , a , 1 : risk));...
% %             toInd(allcomb(10 , 6 , 2 : hpvVaxStates , 10 , ...
% %             1 : endpoints , g , a , 1 : risk))];
% %         hpvHivSus = annlz(sum(popVec(1:end-1 , hpvHivSusInds) , 2)) ./ stepsPerYear;
% %         plot(tVec(1 : stepsPerYear : end-1) , ...
% %             annlz(sum(sum(newHpv(1:end-1 , g , 2 : 6 , a , :) ...
% %             + newHpv(1:end-1 , g , 10 , a , :) ...
% %             + newImmHpv(1:end-1 , g , 2 : 6 , a , :) ...
% %             + newImmHpv(1:end-1 , g , 10 , a , :), 3), 5)))
% %         hold on
% %         xlabel('Year'); ylabel('New Infections'); title([' Age group ' , ages{a} , ' New HPV in HIV+'])
% % 
% %     end
% % end
% % legend('Male' , 'Female')
% 

%% VL breakdown
hivInf = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
totalHiv = sum(popVec(: , hivInf) , 2);
viralPop = zeros(length(tVec) , 6);
for v = 1 : viral
    viralGroup = toInd(allcomb(3 : 7 , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
    viralPop(: , v) = sum(popVec(: , viralGroup) , 2);
end
figure()
subplot(1,2,1);
area(tVec , bsxfun(@rdivide , viralPop , totalHiv));
legend('Acute Infection' , 'VL < 1000' , 'VL 1,000 - 10,000' , 'VL 10,000 - 50,000' ,...
    'VL > 50,000')
xlabel('Year')
ylabel('Proportion of HIV infected')
xlim([1980 2020]);
ylim([0 1]);

subplot(1,2,2);
bar(tVec , viralPop , 'stacked')
legend('Acute Infection' , 'VL < 1000' , 'VL 1,000 - 10,000' , 'VL 10,000 - 50,000' ,...
    'VL > 50,000')
xlabel('Year')
ylabel('HIV Infected')
xlim([1980 2020]);

%% CD4 breakdown
cd4Pop = zeros(length(tVec) , 5);
for d = 3 : 7
    cd4Group = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
    cd4Pop(: , d - 2) = sum(popVec(: , cd4Group) , 2);
end
figure()
subplot(1,2,1);
area(tVec , bsxfun(@rdivide , cd4Pop , totalHiv));
legend('Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
    'CD4 <= 200 cells/uL')
xlabel('Year');
ylabel('Proportion of HIV infected');
xlim([1980 2020]);
ylim([0 1]);

subplot(1,2,2);
bar(tVec , cd4Pop , 'stacked')
legend('Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
    'CD4 <= 200 cells/uL')
xlabel('Year');
ylabel('HIV Infected');
xlim([1980 2020]);

%% HIV by age group
% % hivAge = zeros(length(tVec) , 12);
% % for a = 1 : age
% %     hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %         1 : endpoints , 1 : gender , a , 1 : risk));
% %     hivArt = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %         1 : endpoints , 1 : gender , a , 1 : risk));
% %     hivNeg = toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %         1 : endpoints , 1 : gender , a , 1 : risk));
% %     hivAge(: , a) = sum(popVec(: , hivPos) , 2) + sum(popVec(: , hivArt) , 2);
% % end
% % hivPosAllInd = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %     1 : gender , 1 : age , 1 : risk));
% % hivArtAllInd = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %     1 : gender , 1 : age , 1 : risk));
% % hivNegAllInd = toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %         1 : endpoints , 1 : gender , 1 : age , 1 : risk));
% % hivPosAll = sum(popVec(: , hivPosAllInd) , 2 ) + sum(popVec(: , hivArtAllInd),2);
% % 
% % figure()
% % subplot(1 , 2 , 1)
% % area(tVec , bsxfun(@rdivide , hivAge , hivPosAll));
% % title('HIV Status by Age Group')
% % xlabel('Year')
% % ylabel('Proportion of HIV Positive')
% % legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
% %     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
% %     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')
% % 
% % % HIV by risk group
% % hivRisk = zeros(length(tVec) , risk);
% % for r = 1 : risk
% %     hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %         1 : endpoints , 1 : gender , 1 : age , r));
% %     hivArt = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %         1 : endpoints , 1 : gender , 1 : age , r));
% %     hivRisk(: , r) = sum(popVec(: , hivPos) , 2) + sum(popVec(: , hivArt) , 2);
% % end
% % subplot(1 , 2 , 2)
% % area(tVec , bsxfun(@rdivide , hivRisk , hivPosAll));
% % title('HIV Status by Risk Group')
% % xlabel('Year')
% % ylabel('Proportion of HIV Positive')
% % legend('Low' , 'Medium' , 'High' , 'Location' , 'NorthEastOutside')
% % 
% % figure()
% % subplot(2,1,1)
% % area(tVec , bsxfun(@rdivide , hivAge , sum(popVec , 2)));
% % title('HIV Prevalence by Age Group')
% % xlabel('Year')
% % ylabel('Proportion of HIV Positive')
% % legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
% %     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
% %     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')
% % 
% % subplot(2 , 1 , 2)
% % area(tVec , bsxfun(@rdivide , hivRisk , sum(popVec , 2)));
% % title('HIV Prevalence by Risk Group')
% % xlabel('Year')
% % ylabel('Proportion of HIV Positive')
% % legend('Low' , 'Medium' , 'High' , 'Location' , 'NorthEastOutside')
% 

%% CC prevalence by HIV group and HPV type
% Vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , 6 , 1 : hpvNonVaxStates , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , 6 , 1 : hpvNonVaxStates , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(popVec(: , ccHivNegInds) , 2);
popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

figure();
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'-')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'-')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'-')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' CC Prevalence')
%legend('HIV-' , 'HIV+ noART' , 'ART')
hold all;

% Non-vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , ...
     [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(popVec(: , ccHivNegInds) , 2);
popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

hold all;
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'--')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'--')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'--')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' CC Prevalence')
legend('HIV-' , 'HIV+ noART' , 'ART' , 'HIV-, nonVax' , 'HIV+ noART, nonVax' , 'ART, nonVax')
hold all;

%% CIN3 prevalence by HIV group and HPV type
% Vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , 5 , [1 : 5 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , 5 , [1 : 5 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , 5 , [1 : 5 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(popVec(: , ccHivNegInds) , 2);
popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

figure();
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'-')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'-')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'-')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' CIN3 Prevalence')
%legend('HIV-' , 'HIV+ noART' , 'ART')
hold all;

% Non-vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , ...
     [1 : 4 , 7] , 5 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     [1 : 4 , 7] , 5 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     [1 : 4 , 7] , 5 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(popVec(: , ccHivNegInds) , 2);
popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

hold all;
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'--')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'--')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'--')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' CIN3 Prevalence')
legend('HIV-' , 'HIV+ noART' , 'ART' , 'HIV-, nonVax' , 'HIV+ noART, nonVax' , 'ART, nonVax')
hold all;

%% CIN2 prevalence by HIV group and HPV type
% Vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , 4 , [1 : 4 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , 4 , [1 : 4 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , 4 , [1 : 4 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(popVec(: , ccHivNegInds) , 2);
popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

figure();
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'-')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'-')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'-')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' CIN2 Prevalence')
%legend('HIV-' , 'HIV+ noART' , 'ART')
hold all;

% Non-vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , ...
     [1 : 3 , 7] , 4 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     [1 : 3 , 7] , 4 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     [1 : 3 , 7] , 4 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(popVec(: , ccHivNegInds) , 2);
popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

hold all;
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'--')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'--')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'--')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' CIN2 Prevalence')
legend('HIV-' , 'HIV+ noART' , 'ART' , 'HIV-, nonVax' , 'HIV+ noART, nonVax' , 'ART, nonVax')
hold all;

%% CIN1 prevalence by HIV group and HPV type
% Vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , 3 , [1 : 3 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , 3 , [1 : 3 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , 3 , [1 : 3 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(popVec(: , ccHivNegInds) , 2);
popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

figure();
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'-')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'-')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'-')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' CIN1 Prevalence')
%legend('HIV-' , 'HIV+ noART' , 'ART')
hold all;

% Non-vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , ...
     [1 : 2 , 7] , 3 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     [1 : 2 , 7] , 3 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     [1 : 2 , 7] , 3 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(popVec(: , ccHivNegInds) , 2);
popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

hold all;
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'--')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'--')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'--')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' CIN1 Prevalence')
legend('HIV-' , 'HIV+ noART' , 'ART' , 'HIV-, nonVax' , 'HIV+ noART, nonVax' , 'ART, nonVax')
hold all;

%% HPV prevalence by HIV group and HPV type
% Vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , 2 , [1 : 2 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , 2 , [1 : 2 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , 2 , [1 : 2 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(popVec(: , ccHivNegInds) , 2);
popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

figure();
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'-')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'-')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'-')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' HPV Prevalence')
%legend('HIV-' , 'HIV+ noART' , 'ART')
hold all;

% Non-vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , ...
     [1 , 7] , 2 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     [1 , 7] , 2 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     [1 , 7] , 2 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(popVec(: , ccHivNegInds) , 2);
popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

hold all;
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'--')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'--')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'--')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' HPV Prevalence')
legend('HIV-' , 'HIV+ noART' , 'ART' , 'HIV-, nonVax' , 'HIV+ noART, nonVax' , 'ART, nonVax')
hold all;

