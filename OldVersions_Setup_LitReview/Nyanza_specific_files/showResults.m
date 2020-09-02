function[] = showResults(pathModifier)

%% Load parameters and results
paramDir = [pwd , '\Params\'];

[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd , ...
    annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , partnersMmult, maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4, ...
    mue , mue2 , mue3 , mue4 , epsA_vec , epsR_vec , ...
    yr , ...
    hivOn , betaHIV_mod , muHIV , kCD4 , ...
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
    ccInc2011_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , ...
    hivPrevM_dObs , hivPrevF_dObs , popAgeDist_dObs , totPopSize_dObs , ...
    hivCurr , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , ...
    hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ageInd , riskInd , ...
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    d_partnersMmult, riskAdj, d_riskAdj, ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3  , dFertMat3,  ...
    deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2(1 , 0 , [] , [] , []);

% Load results
resultsDir = [pwd , '\HHCoM_Results\'];
toNowName = ['toNow_30May20_N_increaseClearHIV_increasekCC_5_muART_final']

load([resultsDir ,toNowName]) %change from pathModifier to file name
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)

% excel output file 
filename = [pwd, '\Calibration_comparison_Nyanza.xlsx']

%% 
%% Population size by gender
% figure()
% mInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%     1 : endpoints , 1 : intervens , 1 , 1 : age , 1 : risk));
% fInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%     1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
% plot(tVec , sum(popVec(: , mInds) , 2))
% hold on
% plot(tVec , sum(popVec(: , fInds) , 2))
% legend('Males' , 'Females')
% xlabel('Year'); ylabel('Population')

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
xlim([1980 2120]);
legend('Model, male' ,'Model, female', 'Kenya historical estimates (UN)' , 'Kenya future projections (UN)', ...
    'Location', 'Northwest')
hold off
%%
sheet = ['pop'];
cols1 = {toNowName};
cols2 = {'Year', 'Model pop'}; %, 'UN pop'};
xlswrite(filename, cols1, sheet, 'L1')
xlswrite(filename, cols2, sheet, 'L2')
xlswrite(filename, [tVec(1 : stepsPerYear * 5  : end)', totalPop0_79(1 : stepsPerYear * 5 : end)], sheet, 'L3')
% xlswrite(filename, [historicalPop0_79(:, 2)], sheet,'G8')
% xlswrite(filename, [futurePop(:, 1)], sheet,'E23')
% xlswrite(filename, [futurePop(:, 2)], sheet,'G23')

%% Population size by age vs. validation data

% % Load calibration data from Excel
file = [pwd , '/Config/Kenya_parameters_Feb20.xlsx'];
years = 1990:2020;
ageGroup = {'0-9', '10-19' , '20-29' , '30-39' , '40-49' , '50-59', '60-79'};
popPropYrs = zeros(length(years),5);
popNumYrs = popPropYrs ;
popPropYrs_obs = zeros(7,8);
popPropYrs_Nyan = popPropYrs_obs;
popNumYrs_obs = popPropYrs_obs;
popNumYrs_obs = xlsread(file , 'Population' , 'H178:O184'); %Number
popPropYrs_obs = xlsread(file , 'Population' , 'H169:O175'); %proportions
popPropYrs_Nyan = xlsread(file, 'Population' , 'H198:O201'); % proportions nyanza

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
    popPropYrs_obs(:, 1) , popPropYrs_obs(:, 8) , 'o', ...
    popPropYrs_Nyan(:, 1), popPropYrs_Nyan(:, 2) , 'x', ...
    popPropYrs_Nyan(:, 1), popPropYrs_Nyan(:, 3) , 'x', ...
    popPropYrs_Nyan(:, 1), popPropYrs_Nyan(:, 4) , 'x', ...
    popPropYrs_Nyan(:, 1), popPropYrs_Nyan(:, 5) , 'x', ...
    popPropYrs_Nyan(:, 1), popPropYrs_Nyan(:, 6) , 'x', ...
    popPropYrs_Nyan(:, 1), popPropYrs_Nyan(:, 7) , 'x', ...
    popPropYrs_Nyan(:, 1), popPropYrs_Nyan(:, 8) , 'x');
colororder()
ylabel('Proportions'); xlabel('Year'); title('Age distribution in broad groups'); 
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
ylabel('Number of individuals'); xlabel('Year'); title('Age distribution in broad groups'); 
legend('0-9', '10-19' , '20-29' , '30-39' , '40-49' , '50-59', '60-79',...
    '0-9 obs', '10-19 obs' , '20-29 obs' , '30-39 obs' , '40-49 obs' , '50-59 obs', '60-79 obs', ...
    'Location', 'NorthEastOutside');
%%
sheet = ['Pop_by_Age'];
cols1 = {toNowName};
cols2 = {'Year','0-9', '10-19' , '20-29' , '30-39' , '40-49' , '50-59', '60-79'}; %,...
    %'0-9 obs', '10-19 obs' , '20-29 obs' , '30-39 obs' , '40-49 obs' , '50-59 obs', '60-79 obs'};
xlswrite(filename, cols1, sheet, 'A75')
xlswrite(filename, cols2, sheet, 'A76')
xlswrite(filename, [years(1: 5 : end)', popPropYrs(1 : 5 : end, :)], sheet, 'A77')
xlswrite(filename,[years(1: 5 : end)', popNumYrs(1 : 5 : end, :)], sheet, 'J77')


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
xlswrite(filename, cols1, sheet, 'P1')
xlswrite(filename, cols2, sheet, 'P2')
xlswrite(filename, [fertilityVec(1:5:end, :)], sheet, 'P3')
%xlswrite(filename, [fertilityVal], sheet, 'C8')


%% HIV prevalence vs. HIV data by year
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

% Compared to ANC data
HIV_ANC_Kisumu(1 , :) = [1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2008 2010 2011 ] ;
HIV_ANC_Kisumu(2 , :) = [18 18 19 19 29 24 26 32 27 25 33 29 26 26 11.2 15.1 18.5 16.9 18.5 15.5];
% 
% Compared to National AIDS Control Council (Kenya)
HIV_Kenya_spectrum(1, :) = [1983 1987 1990 1994 1999 2003 2007 2011 2015 2019];
HIV_Kenya_spectrum(2, :) = [0.002 0.0209 0.1233 0.2431 0.2302 0.171 0.1372 0.1419 0.1456 0.1252 ] .* 100;
% 
% 
figure()

plot(tVec , hivPrev , HIV_ANC_Kisumu(1 , :) , HIV_ANC_Kisumu(2 , :) , '*', ...
    HIV_Kenya_spectrum(1, :), HIV_Kenya_spectrum(2, :), 'o')
hold on 
% yPosError = abs(upper_prevVal - prevVal);
% yNegError = abs(lower_prevVal - prevVal);
% errorbar(prevValYrs , prevVal , yNegError , yPosError , 'ms')
xlabel('Year'); ylabel('Proportion of Population (%)'); title('HIV Prevalence (Ages 15-49)')
legend('Model' , 'ANC data (Kisumu)', 'Spectrum data (Nyanza)')
xlim([1980 2020])
ylim([0, 40])

%%
sheet = ['HIV_prev']
cols1 = {toNowName};
cols2 = {'Final_2'}; %, 'ANC data (Kisumu)', 'Year', 'Spectrum data (Nyanza)'};
xlswrite(filename, cols1, sheet, 'F1')
xlswrite(filename, cols2, sheet, 'F2')
xlswrite(filename, [hivPrev(331:stepsPerYear:end) ], sheet, 'F3')
%xlswrite(filename, [HIV_ANC_Kisumu'], sheet, 'C3')
%xlswrite(filename, [HIV_Kenya_spectrum'], sheet, 'E3')

%% HIV prevalance, all ages by gender
figure()
hivObsGender = zeros(4,3)
hivObsGender(:,3) = [1998 2003 2007 2008]; 
hivObsGender(:,1) = [19.8 12.28 11.0 11.6]; 
hivObsGender(:,2) = [30.1 18.25 18.0 15.97]; 

sheet = ['HIV_by_sex']
cols1 = {toNowName};
cols2 = {'Male, Final_2', 'Female, Final_2'} %, 'Year', 'Males, DHS/KAIS', 'Females, DHS/KAIS',};
xlswrite(filename, cols1, sheet, 'L1')
xlswrite(filename, cols2, sheet, 'L2')

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

    cell1 = ['L', 'M'];
    cell = ([cell1(g) +'3']);
    xlswrite(filename, [hivPrev_sex(331:stepsPerYear:end)], sheet, cell)
end
xlabel('Year')
xlim([1980 2020])
ylabel('Prevalence')
title('HIV Prevalence (aged 14-59)')
legend('Males, model' , 'Males, DHS/KAIS', 'Females, model', 'Females, DHS/KAIS')

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

%% HIV prevalence by age group
figure()

years = 1980:2020;
ageGroup = {'0-14', '15-29' , '30-39' , '40-49' , '50-59', '60-79'};
hivPrevYrs = zeros(length(years),7);
ageVec = {[1:2], [4:6] , [7:8] , [9:10] , [11:12], [13:16]};

for y = 1 : length(years)
    yearCurr = years(y);
    for aInd = 1 : length(ageVec)
        a = ageVec{aInd};
        hivAge = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2   , a , 1 : risk));
        hivPrevYrs(y,aInd) = sum(popVec(((yearCurr - startYear) * stepsPerYear +1) , hivAge),2) ...
            ./ sum(popVec(((yearCurr - startYear) * stepsPerYear +1) , popTot),2);
    end
end
plot(years , hivPrevYrs);
xlabel('Year')
ylabel('Prevalence')
legend('0-14',  '15-29' , '30-39' , '40-49' , '50-59', '60-79')
xlim([1980 2020])
title('HIV Prevalence by age in women')

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
plot(circPropYr_obs , circProp_obs , 'o');
xlabel('Year')
ylabel('Proportion of HIV-Negative Males Circumcised by Broad Age Groups (%)')
title('Circumcision Indicator')
xlim([1960 2020]);
ylim([0 100])
grid on;
legend('15-19, model' , '20-24, model ' , '25-29, model' ,...
    '30-39, model', '40-49, model', '50+, model' , ...
    '15-19, observed' , '20-24, observed' , '25-29, observed' ,...
    '30-39, observed', '40-49, observed', '50+, observed' ,...
    'Location' , 'NorthWest');
   
%% HIV prevalence by age on x-axis
genderVec = {'Males (on and off ART)' , 'Females (on and off ART)'};
hiv2003 = zeros(7 , 2);
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};
file = [pwd , '/Config/Kenya_parameters_Feb20.xlsx'];
DHS2003_Nyanza = zeros(7,6);
DHS2008_Nyanza = zeros(7,6);
DHS2003_Nyanza(:,1:3) = xlsread(file , 'HIV prevalence' , 'R4:T10'); % Males estimates
DHS2003_Nyanza(:,4:6) = xlsread(file , 'HIV prevalence' , 'O4:Q10'); % Female estimates
DHS2008_Nyanza(:,1:3) = xlsread(file , 'HIV prevalence' , 'Q62:S68'); % Males estimates
DHS2008_Nyanza(:,4:6) = xlsread(file , 'HIV prevalence' , 'N62:P68'); % Female estimates

Glynn1998_Kisumu(:, 1) = [3.5 18.3 18.3 33.1 33.1 27.7 27.7];
Glynn1998_Kisumu(:, 2) = [23.2 38.8 37.7 29.8 29.9 NaN NaN];

% DHS2003_Nyanza(:, 1) = [0.1419 5.6208 23.1835 17.7284 19.6896 24.6223 18.6658];
% DHS2003_Nyanza(:, 2) = [4.3396 29.9703 21.8535 16.3208 17.9828 33.1712 15.8049];
% 
% DHS2008_Nyanza(:, 1) = [1.83 5.75 24.33 15.30 22.51 24.73 12.91];
% DHS2008_Nyanza(:, 2) = [10.75 11.91 22.20 25.58 22.25 9.27 17.16 ];

figure;
for g = 1 : gender
    %aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
    for a = 4: 10
        %a = aVec{aInd};
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hiv1998(a-3 , g) = (sum(popVec((1998 - startYear)* stepsPerYear,hivInds),2) ... 
            /sum(popVec((1998 - startYear)* stepsPerYear,ageInds),2))*100;
    end
end
    subplot(3 , 1 , 1)
    plot(1 : size(hiv1998 , 1) , hiv1998(: , :) , '-', 1 : size(hiv1998 , 1), Glynn1998_Kisumu(:,:), '--o');
    hold all;
    xlabel('Age Group'); ylabel('HIV Prevalence')
    set(gca , 'xtick' , 1 : length(hiv1998) , 'xtickLabel' , ageGroup(4:10));
    title('1998 HIV prevalence by sex and age')
    ylim([0 70])
    legend({'Males (model)', 'Females (model)', '1998 Kisumu survey males', '1998 Kisumu survey females'}, ...
        'Location' , 'northeast')
    grid on;
    
for g = 1 : gender
    %aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
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
    subplot(3 , 1 , 2)
    plot(1 : size(hiv2003 , 1) , hiv2003(: , :) , '-');
    hold all;
    errorbar( 1 : size(hiv2003 , 1) , DHS2003_Nyanza(:, 1), DHS2003_Nyanza(:,1) - DHS2003_Nyanza(: , 2), ...
    DHS2003_Nyanza(:, 3)-DHS2003_Nyanza(: , 1));
    hold all;
    errorbar( 1 : size(hiv2003 , 1) , DHS2003_Nyanza(:, 4), DHS2003_Nyanza(:,4) - DHS2003_Nyanza(: , 5), ...
    DHS2003_Nyanza(:, 6)-DHS2003_Nyanza(: , 4));
    hold all;
    xlabel('Age Group'); ylabel('HIV Prevalence')
    set(gca , 'xtick' , 1 : length(hiv2003) , 'xtickLabel' , ageGroup(4:10));
    title('2003 HIV prevalence by sex and age')
    ylim([0 70])
    legend('Males (model)', 'Females (model)', '2003 DHS males', '2003 DHS females', ...
        'Location' , 'northeast')
    grid on;
    
for g = 1 : gender
    for a = 4: 10
        %a = aVec{aInd};
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hiv2008(a-3 , g) = (sum(popVec((2008 - startYear)* stepsPerYear,hivInds),2) ... 
            /sum(popVec((2008 - startYear)* stepsPerYear,ageInds),2))*100;
    end
end
    subplot(3 , 1 , 3)
    plot(1 : size(hiv2008 , 1) , hiv2008(: , :) , '-' ); 
    hold all;
    errorbar( 1 : size(hiv2008 , 1) , DHS2008_Nyanza(:, 1), DHS2008_Nyanza(:,1) - DHS2008_Nyanza(: , 2), ...
    DHS2008_Nyanza(:, 3)-DHS2008_Nyanza(: , 1));
hold all;
    errorbar( 1 : size(hiv2008 , 1) , DHS2008_Nyanza(:, 4), DHS2008_Nyanza(:,4) - DHS2008_Nyanza(: , 5), ...
    DHS2008_Nyanza(:, 6)-DHS2008_Nyanza(: , 4));
hold all;
    xlabel('Age Group'); ylabel('HIV Prevalence')
    set(gca , 'xtick' , 1 : length(hiv2008) , 'xtickLabel' , ageGroup(4:10));
    title('2008 HIV prevalence by sex and age')
    ylim([0 70])
    legend('Males (model)', 'Females (model)', '2008 DHS males', '2008 DHS females', ...
        'Location' , 'northeast')
    %legend('Without ART dropout' , 'With ART dropout');
    grid on;
    
%%
sheet = ['HIV_by_age'];
cols1 = {toNowName};
% cols2 = {'Age group', 'Male (1998)', 'Female (1998)', ' Male (Kisumu 1998)', 'Female (Kisumu 1998)', ...
%     'Male (2003)', 'Female (2003)', 'Female (DHS 2003)', 'Female - LB (DHS 2003)', 'Female - UB (DHS 2003)',...
%     'Male (DHS 2003)', 'Male - LB (DHS 2003)', 'Male - UB (DHS 2003)',...
%     'Male (2008)', 'Female (2008)', 'Female (DHS 2008)', 'Female - LB (DHS 2008)', 'Female - UB (DHS 2008)',...
%     'Male (DHS 2008)', 'Male - LB (DHS 2008)', 'Male - UB (DHS 2008)'};
cols2 = {'Male (1998)', 'Female (1998)', ...
    'Male (2003)', 'Female (2003)', ...
    'Male (2008)', 'Female (2008)'};

xlswrite(filename, cols1, sheet, 'A32')
xlswrite(filename, cols2, sheet, 'A33')
% xlswrite(filename, ageGroup(4:10)', sheet, 'A24')
xlswrite(filename, [hiv1998], sheet, 'A34')
%xlswrite(filename, [Glynn1998_Kisumu], sheet, 'D3')
xlswrite(filename, [hiv2003], sheet, 'C34')
%xlswrite(filename, [DHS2003_Nyanza], sheet, 'H3')
xlswrite(filename, [hiv2008], sheet, 'E34')

%% HIV prevalance, by sex and age over time 
hivSexAge = zeros(length(tVec), 7);
ageGroup = {'15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' };
yr = 1980 : 2020;
sheet = ['HIV_sex_age']
cols1 = {toNowName};
cols2 = {'FutureSim'} %, 'F, 005% riskAdj'} %, 'Year', 'Males, DHS/KAIS', 'Females, DHS/KAIS',};
xlswrite(filename, cols1, sheet, 'A11')
xlswrite(filename, cols2, sheet, 'A12')
xlswrite(filename, ['Year', 'Sex', ageGroup], sheet, 'A13')

for g = 1  : 2
    for a = 4: 10
        %a = aVec{aInd};
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivSexAge(:, g, a-3) =  sum(popVec(: , hivInds) , 2) ... 
            ./ sum(popVec(:,ageInds),2)*100;
    
    end
%     plot(tVec , hivSexAge(g, :))
%     hold on
%     plot(hivObsGender(:, 3), hivObsGender(:, g), 'o')
%     hold on
end

   sex_label = repmat(["M" "F"], 41, 1)
    xlswrite(filename, [yr', sex_label(:, 1), squeeze(hivSexAge(331:stepsPerYear:end, 1, :))], sheet, "A14")
    xlswrite(filename, [yr', sex_label(:, 2), squeeze(hivSexAge(331:stepsPerYear:end, 2, :))], sheet, "A55")
%xlswrite(filename, [hivObsGender(:, 3), hivObsGender(:, 1:2)], sheet, 'E3')

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

%% Proportion of CD4-eligible population on ART
% figure()
% % artActual = [0	0	0	0	1	2	3	6, ...
% %     9	14	19	27	34	40	45	48];
% % yrsArtActual = [2000	2001	2002	2003	2004	2005 ...
% %     2006	2007	2008	2009	2010	2011	2012	2013 ...
% %     2014	2015];
% artUNAIDS = [31 38 42 45 51 60 67 73 69];
% yrsUNAIDS = 2010:2018;
% artWB = [10 14 19 27 33 41 46 49 55];
% yrsWB = 2006:2014;
% artMOH = [73, 67, 80];
% yrsMOH = [2013 2015 2017]
% 
% artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%     1 : endpoints , 1 : intervens , 1 : gender , 3 : age , 1 : risk));
% hiv200Inds = toInd(allcomb(7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%     1 : endpoints , 1 : intervens , 1 : gender , 3 : age , 1 : risk));
% hiv350Inds = toInd(allcomb(6 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%     1 : endpoints , 1 : intervens , 1 : gender , 3 : age , 1 : risk));
% hiv500Inds = toInd(allcomb(5 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%     1 : endpoints , 1 : intervens , 1 : gender , 3 : age , 1 : risk));
% hivAllInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%     1 : endpoints , 1 : intervens , 1 : gender , 3 : age , 1 : risk));
% 
% plot(tVec , 100 * sum(popVec(: , artInds) , 2) ./ (sum(popVec(: , hivAllInds) , 2) + sum(popVec(: , artInds) , 2)) , ...
%     tVec(1:(2010-startYear)*stepsPerYear) , 100 * sum(popVec(1:(2010-startYear)*stepsPerYear , artInds) , 2) ./ (sum(popVec(1:(2010-startYear)*stepsPerYear , hiv200Inds) , 2) + sum(popVec(1:(2010-startYear)*stepsPerYear , artInds) , 2)) , ...
%     tVec((2010-startYear)*stepsPerYear+1 : (2013-startYear)*stepsPerYear) , 100 * sum(popVec((2010-startYear)*stepsPerYear+1 : (2013-startYear)*stepsPerYear , artInds) , 2) ./ (sum(popVec((2010-startYear)*stepsPerYear+1 : (2013-startYear)*stepsPerYear , hiv350Inds) , 2) + sum(popVec((2010-startYear)*stepsPerYear+1 : (2013-startYear)*stepsPerYear , artInds) , 2)) , ...
%     tVec((2013-startYear)*stepsPerYear+1 : (2015-startYear)*stepsPerYear) , 100 * sum(popVec((2013-startYear)*stepsPerYear+1 : (2015-startYear)*stepsPerYear , artInds) , 2) ./ (sum(popVec((2013-startYear)*stepsPerYear+1 : (2015-startYear)*stepsPerYear , hiv500Inds) , 2) + sum(popVec((2013-startYear)*stepsPerYear+1 : (2015-startYear)*stepsPerYear , artInds) , 2)) , ...
%     tVec((2015-startYear)*stepsPerYear+1 : end) , 100 * sum(popVec((2015-startYear)*stepsPerYear+1 : end , artInds) , 2) ./ (sum(popVec((2015-startYear)*stepsPerYear+1 : end , hivAllInds) , 2) + sum(popVec((2015-startYear)*stepsPerYear+1 : end , artInds) , 2)) , ...
%     yrsUNAIDS , artUNAIDS , '*' , ...
%     yrsWB , artWB , '*' , ...
%     2018 , 78.8 , '*' , ...
%     2012 , 38.3 , '*' , ...    
%     yrsMOH , artMOH , 'o')
% xlabel('Year')
% ylabel('Proportion of HIV Population')
% title('Proportion on ART')
% legend('Model all HIV+' , ...
%     'Model CD4 <= 200' , ...
%     'Model CD4 <= 350' , ...
%     'Model CD4 <= 500' , ...
%     'Model any CD4' , ...
%     'UNAIDS (Nyanza)' , ...
%     'WB (Kenya)' , ...
%     'PHIA (Nyanza)' , ...
%     'KAIS (Nyanza)' , ...
%     'MOH reports (Nyanza)', ...
%     'Location', 'Northwest')

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
    xlim([1975 2020])
    ylim([0 80])
    xlabel('Year'); ylabel('Prevalence (%)'); title([genders{g} , ' HPV Prevalence (ages 15+)'])
    legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Location' , 'NorthWest')

%
sheet = ['hpvPrev'];
cols1 = {[toNowName, '2.2xC3toCC']};
cols2 = {[genders{g},' Gen Pop'], [genders{g},' HIV-neg'], [genders{g},' HIV-pos'], [genders{g},' ART']};
xlswrite(filename, cols1, sheet, 'R1')
cell = ['R', 'V'];

xlswrite(filename, cols2, sheet, [cell(g) +'2'])
xlswrite(filename, [genPopHPV(1 : stepsPerYear * 5 : end), ...
    hivNegHPV(1 : stepsPerYear * 5 : end), hivPosHPV(1 : stepsPerYear * 5 : end), ...
    artHPV(1 : stepsPerYear * 5 : end)], sheet, [cell(g) +'3'])
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
    hpvInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
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
title(['Nyanza HPV prevalence in women ', year])
ylim([0 100])
%%
sheet = ['HPV_by_age_2005'];
cols1 = {toNowName};
cols2 = {'HIV+, FutureSim', 'HIV-, FutureSim'} %, 'DeVuyst HIV+ 2009', 'Yamada HIV+', 'Yamada HIV-'};
xlswrite(filename, cols1, sheet, 'V1')
xlswrite(filename, cols2, sheet, 'V2')
%xlswrite(filename, ageGroup(1:9)', sheet, 'F3')
xlswrite(filename, [hpvHIV2005, hpvNeg2005], sheet, 'V3')

%% Age specific HPV prevalence data 
ageGroup = {'17 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
hpv2000 = zeros(9 , 1);
hpvHIV2000 = hpv2000;
hpvNeg2000 = hpv2000;

yr = 2020;
year = {yr};
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
title(['HPV prevalence among women in ', year])
ylim([0 100])

% sheet = ['HPV_by_age_2000'];
% cols1 = {toNowName,'HPV prevalence among women in 2000'};
% %cols2 = {'Age', 'Model Gen Pop'}; %, 'DeVuyst gen pop', 'DeVuyst HIV+'};
% xlswrite(filename, cols1, sheet, 'I1')
% %xlswrite(filename, cols2, sheet, 'I2')
% xlswrite(filename, ageGroup(1:9)', sheet, 'I3')
% xlswrite(filename, [hpv2000], sheet, 'J3')


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

yr = 2010;
year = {yr};
cinHiv_ccScreen = [0.18	0.105 0.089] .* 100; % Observed, HIV+, Nairobi, Memiah, 2013
cinHiv_ccScreen2 = [33.1 10.5 12.8]; %Observed, HIV clinic, De Vuyst 2013
cinHIV_FSW = [19.1 NaN 13.1]; %Observed, HIV+ FSW, Nairobi, Sweet 2020
cinNeg_FSW = [8.0 NaN 1.5]; % Observed, HIV- FSW, Nairobi, Sweet 2020
cinHIV_FSW2 = [17.9 NaN 11.9]; %Observed, HIV+ FSW, Nairobi, Njagi 2013
cinNeg_FSW2 = [9.0 NaN 1.88]; %Observed, HIV- FSW, Nairobi, Njagi 2013
cinHiv2010 = zeros(3 , 1);
cinHivNeg2010 = cinHiv2010;

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
plot(1 : length(cinHiv2010) , cinHiv2010 , 'o-' , 1 : length(cinHiv2010), cinHivNeg2010, 'o-' ,...
    1 : length(cinHiv2010), cinART2010, '--', 1 : length(cinHiv2010), cinNoART2010, '--*')
hold on
p = plot(1: length(cinGroup) , cinHiv_ccScreen , '+',  1: length(cinGroup),cinHiv_ccScreen2, '+',...
    1 : length(cinGroup), cinHIV_FSW, 'g+', 1 : length(cinGroup), cinNeg_FSW, 'go',...
    1 : length(cinGroup), cinHIV_FSW2, 'm+', 1 : length(cinGroup), cinNeg_FSW2, 'mo');
p(1).MarkerSize = 10;
p(2).MarkerSize = 8;
ylabel('Prevalence (%)')
set(gca , 'XTick', 1:3, 'xtickLabel' , cinGroup);
xlabel('CIN Stage')
legend('Model HIV+' , 'Model HIV-', 'Model ART+', 'Model ART-', ...
    'HIV+ screening pop (Memiah)', 'HIV+ screening pop (DeVuyst)', ...
    'HIV+ FSW (Sweet)', 'HIV- FSW (Sweet)', 'HIV+ FSW (Njagi)', 'HIV- FSW (Njagi)')
title(['CIN Prevalence by HIV status in ', year])
text(2.2, -3.25, 'Note: All observed data are from Nairobi')
%annotation('textbox', [0, 0, 0.5, 0.1], 'string', 'All observed data are from Nairobi')

%%
sheet = ['CIN_by_HIV_2010'];
cols1 = {toNowName, 'CIN Prevalence by HIV status in 2010'};
cols2 = {'HIV+, FutureSim', 'HIV-, FutureSim'}; %, 'HIV+ screening pop (Memiah)',...
    %'HIV+ screening pop (DeVuyst)', ...
    %'HIV+ FSW (Sweet)', 'HIV- FSW (Sweet)', 'HIV+ FSW (Njagi)', 'HIV- FSW (Njagi)'};
xlswrite(filename, cols1, sheet, 'X1')
xlswrite(filename, cols2, sheet, 'X2')
%xlswrite(filename, cinGroup', sheet, 'A23')
xlswrite(filename, [cinHiv2010, cinHivNeg2010], sheet, 'X3')
% xlswrite(filename, [cinHiv2010, cinHivNeg2010, cinHiv_ccScreen',cinHiv_ccScreen2',...
%     cinHIV_FSW', cinNeg_FSW', cinHIV_FSW2', cinNeg_FSW2'], sheet, 'B9')
 
%% CIN2/3 prevalence for All HR HPV types combined by age in 2000 vs. deVuyst 2003 data
cin1Pos2000 = zeros(10 , 1);
cin3Pos2000 = cin1Pos2000;
ageGroup = {'17-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};
%aVec = {18:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
for a = 4 : 13  %note, age group 4 is 17-19 in the data
    %a = aVec{aInd};
    % CIN 1 prevalence 
    cin1Inds = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 3 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cin1Pos2000(a - 3) = (sum(popVec((2000 - startYear) * stepsPerYear , cin1Inds)))...
        ./ sum(popVec((2000 - startYear) * stepsPerYear , ageInds)) * 100;
    % CIN 2/3 prevalence 
    cin3Inds = unique([toInd(allcomb(1 : disease , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cin3Pos2000(a - 3) = (sum(popVec((2000 - startYear) * stepsPerYear , cin3Inds)))...
        ./ sum(popVec((2000 - startYear) * stepsPerYear , ageInds)) * 100;
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
% 

sheet = ['CIN_by_age_2000'];
cols1 = {toNowName, 'Age specific CIN prevalence among all women in 2000'};
cols2 = {'Age', 'Model CIN 1', 'Model CIN 2/3'}; % , 'DeVuyst CIN 1', 'DeVuyst CIN 2/3'};
xlswrite(filename, cols1, sheet, 'O1')
xlswrite(filename, cols2, sheet, 'O2')
xlswrite(filename, ageGroup', sheet, 'O3')
xlswrite(filename, [cin1Pos2000, cin3Pos2000], sheet, 'P3')
% xlswrite(filename, [cin1Pos2000, cin3Pos2000, cinPosAct], sheet, 'H3')


%% HPV prevalence by age and gender over time
% % hivAge = zeros(age , length(tVec));
% ageGroup = {'10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
%     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
%     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
% figure()
% for g = 1 : gender
%     %aVec = {11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
%     %for aInd = 1 : 14
%     for a = 3 : age
%         %a = aVec{aInd};
%         hpvAgeInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvVaxStates , 1 : 7 , 1 : endpoints , ...
%             g , a , 1 : risk));
%         ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
%             g , a , 1 : risk));
%         hpvAgeRel = bsxfun(@rdivide , sum(popVec(: , hpvAgeInds) , 2)' , sum(popVec(: , ageInds) , 2)') * 100;
%         subplot(5 , 3 , aInd)
%         hold on
%         plot(tVec , hpvAgeRel);
%         xlabel('Year'); ylabel('% HPV'); title([' Age group ' , ageGroup{aInd} , ' HPV Prevalence'])
%     end
% end
% legend('Male' , 'Female')

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

%% Cervical cancer prevalence by age in 2017
% % ccAgeRel = zeros(age , 1);
% % ccAgeNegRel = ccAgeRel;
% % ccAgePosRel = ccAgeRel;
% % ccNegPosArt = zeros(age , 2);
% % ccArtRel = ccAgeRel;
% % for a = 1 : age
% %     % Total population
% %     ccInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 5 : 7 , 1 : endpoints , ...
% %         2 , a  , 1 : risk));
% %     ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %             2 , a , 1 : risk));
% %     ccAgeRel(a) = sum(popVec((2017 - startYear) * stepsPerYear , ccInds) , 2) ...
% %         / sum(popVec((2017 - startYear) * stepsPerYear , ageInds) , 2) * 100;
% % 
% %     % HIV Negative
% %     ccHivNegInds = toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 5 : 7 , 1 : endpoints , ...
% %         2 , a  , 1 : risk));
% %     ageNegInds = toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %             2 , a , 1 : risk));
% %     ccAgeNegRel(a) = sum(popVec((2017 - startYear) * stepsPerYear , ccHivNegInds) , 2) ...
% %         / (sum(popVec((2017 - startYear) * stepsPerYear , ageNegInds) , 2)) * 100;
% % 
% %     % HIV Positive
% %     ccHivPosInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 5 : 7 , 1 : endpoints , ...
% %         2 , a  , 1 : risk));
% %     agePosInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         2 , a , 1 : risk));
% %     ccAgePosRel(a) = sum(popVec((2017 - startYear) * stepsPerYear , ccHivPosInds) , 2) ...
% %         / (sum(popVec((2017 - startYear) * stepsPerYear , agePosInds) , 2)) * 100;
% % 
% %     % On ART
% %     ccArtInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 5 : 7 , 1 : endpoints , ...
% %         2 , a  , 1 : risk));
% %     ageArtInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         2 , a , 1 : risk));
% %     ccArtRel(a) = sum(popVec((2017 - startYear) * stepsPerYear , ccArtInds) , 2) ...
% %         / (sum(popVec((2017 - startYear) * stepsPerYear , ageArtInds) , 2)) * 100;
% % 
% %     % Each group relative to total female population
% %     % HIV-
% %     ccNegPosArt(a , 1) = sum(popVec((2017 - startYear) * stepsPerYear , ccHivNegInds) , 2) ...
% %         / sum(popVec((2017 - startYear) * stepsPerYear , ageInds) , 2) * 100;
% %     % HIV+
% %     ccNegPosArt(a , 2) = sum(popVec((2017 - startYear) * stepsPerYear , ccHivPosInds) , 2) ...
% %         / sum(popVec((2017 - startYear) * stepsPerYear , ageInds) , 2) * 100;
% %     % ART
% %     ccNegPosArt(a , 3) = sum(popVec((2017 - startYear) * stepsPerYear , ccArtInds) , 2) ...
% %         / sum(popVec((2017 - startYear) * stepsPerYear , ageInds) , 2) * 100;
% % 
% % end
% % 
% % figure()
% % plot(1 : length(ccAgeRel) , ccAgeRel , '-o' , 1 : length(ccAgeNegRel) , ...
% %     ccAgeNegRel , '-o' , 1 : length(ccAgePosRel) , ccAgePosRel , '-o' , ...
% %     1 : length(ccArtRel) , ccArtRel , '-o');
% % xlabel('Age Group'); ylabel('Prevalence (%)')
% % set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
% % title('Cervical Cancer Prevalence in 2017')
% % legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Location' , 'NorthWest')
% % 
% % % figure()
% % % bar(1 : length(ccAgeRel) , ccNegPosArt , 'stacked')
% % % hold on
% % % plot(1 : length(ccAgeRel) , ccAgeRel , '-o')
% % % xlabel('Age Group'); ylabel('Prevalence (%)')
% % % set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
% % % title('Cervical Cancer Prevalence in 2017')
% % % legend('HIV-' , 'HIV+' , 'ART' , 'Location' , 'NorthWest')
% 
%% Cervical cancer incidence by age
%ccIncYears = [2017 , 2003 , 1994 , 2012];
ccIncYears = [2010];
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
% 2004-2008 (N)
ken_canreg =[1		2	0
12		3	1.67
34		7	3.63
86		21	14.37
156		27	41.72
187		35	70.34
163		36	74.19
144		27	127.79
120		21	148.35
90		15	193.66
57		11	199.95
40		9	170.4
48		15	134.78
];

for y = 1 : length(ccIncYears)
    ccIncYear = ccIncYears(y);
    
    figure()
    % Plot model outputs
    plot(1 : size(ccAgeRel , 1) , ccAgeRel(: , y) , '-ko' , 1 : size(ccAgeNegRel(: , y) , 1) , ...
        ccAgeNegRel(: , y) , '-kp' , 1 : size(ccAgePosRel , 1) , ccAgePosRel(: , 5 , y) , '-k+' , ...
        1 : size(ccArtRel , 1) , ccArtRel(: , y) , '-k^');
    hold on
    % Plot observed data
     plot(4 : age ,  ken_canreg(:, 1)  , 'r--' , 4 : age ,  ken_canreg(:, 2), '--' );
    xlabel('Year'); ylabel('Incidence per 100,000');
    set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
    title(['Nyanza cervical Cancer Incidence by age in ' num2str(ccIncYear)]);
       legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Nairobi 2007-11' , 'Eldoret 2008-11' )
%     figure()
%     bar(1 : length(ccNegPosArt(: , : , y)) , ccNegPosArt(: , : , y), 'stacked')
%     xlabel('Age Group'); ylabel('Incidence per 100,000')
%     set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
%     title(['Cervical Cancer Incidence Distribution in ' , num2str(ccIncYear)])
%     legend('HIV-' , 'HIV+' , 'ART')
end
% 
% sheet = ['CC_inc_by_age_2009'];
% cols1 = {toNowName, 'Cervical Cancer Incidence by age in 2009'};
% cols2 = {'Age', 'General' , 'HIV-' , 'HIV+' , 'ART'}; % , 'Nairobi 2007-11' , 'Eldoret 2008-11'};
% xlswrite(filename, cols1, sheet, 'V1')
% xlswrite(filename, cols2, sheet, 'V2')
% xlswrite(filename, ageGroup', sheet, 'V3')
% xlswrite(filename, [ccAgeRel, ccAgeNegRel, ccAgePosRel(:, 5), ccArtRel], sheet, 'W3')
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
title('Cervical Cancer Incidence in Nyanza');
xlim([1980 2020]);
legend('General' , 'HIV-' , 'HIV+, no ART' , 'HIV+, ART', 'GBD Kenya 2018', 'Globocan E. Africa', 'Globocan Kenya',...
    'Location', 'NorthWest')

%%
sheet = ['CC_inc'];
cols1 = {toNowName, 'Cervical Cancer Incidence'};
cols2 = {'General, FutureSim' , 'HIV-, FutureSim' , 'HIV+, no ART, FutureSim' , 'HIV+, ART, FutureSim'}; %,...
   % 'Year', 'GBD Kenya 2018', 'GBD LB', 'GBD UB', 'Year', 'Globocan Kenya'};
xlswrite(filename, cols1, sheet, 'P1')
xlswrite(filename, cols2, sheet, 'P2')
xlswrite(filename, [ ccAgeRel', ccAgeNegRel', ccAgePosRel',...
    ccArtRel'], sheet, 'P3')
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
sheet = ['CC_inc_standardized'];
cols1 = {toNowName, 'Age-standardized cervical cancer incidence'};
cols2 = {'Year', 'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all' }; %,...
   % 'Year', 'GBD Kenya 2018', 'GBD LB', 'GBD UB', 'Year', 'Globocan Kenya'};
xlswrite(filename, cols1, sheet, 'AB1')
xlswrite(filename, cols2, sheet, 'AB2')
xlswrite(filename, [tVec(1: stepsPerYear :end-1)'] , sheet, 'AB3')

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
    
    cell = ['AC'; 'AD'; 'AE'; 'AF'; 'AG'];
    xlswrite(filename, ccInc', sheet, [cell(i, 1:2) +'3'])

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

    
%% CC Cumulative Probability of Incidence- early years
% % ccIncYears = [1980,1990,2000,2010];
% % ccAgeCI = zeros(1 , length(ccIncYears));
% % 
% % fScale = 10^5;
% % ageGroup = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 - 24' , '25 - 29' ,...
% %     '30 - 34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
% %     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
% % annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 
% % ccYrs = ((ccIncYears - startYear) * stepsPerYear :...
% %     (ccIncYears + 1 - startYear) * stepsPerYear);
% % 
% % for y = 1 : length(ccIncYears)
% %     for a = 1 : age
% %         % Year
% %         yr_start = (ccIncYears(y) - 1 - startYear)  .* stepsPerYear;
% %         yr_end = (ccIncYears(y) - startYear) .* stepsPerYear - 1;
% %         % Total population
% %         ageInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %             2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 8 : 10 , 1 : endpoints , ...
% %             2 , a , 1 : risk))];
% %         ccAgeCI(1 , y) = ccAgeCI(1 , y) + (-1) .* annlz(sum(sum(sum(newCC(yr_start : yr_end , ...
% %             1 : disease , 1 : hpvVaxStates , a) , 2) , 3) , 4)) ...
% %             ./ (annlz(sum(popVec(yr_start : yr_end , ageInds) , 2)) ...
% %             ./ stepsPerYear) ;
% %     end
% %     ccAgeCI(1 , y) = 1 - exp(ccAgeCI(1 , y));
% % end
% % 
% % forouzanfar =[4.7
% % 4.3
% % 3.9
% % 3.4] ./ 100;
% % 
% % forouzanfar_ub = [6.1
% % 5.5
% % 5.7
% % 5.3] ./ 100;
% % 
% % forouzanfar_lb = [3.0
% % 2.7
% % 3.0
% % 2.5] ./ 100;
% % 
% % figure()
% % plot(ccIncYears, ccAgeCI(: , y) , '-o');
% % xlabel('Year'); ylabel('Cumulative Probability of Incidence')
% % title(['Cumulative Probability of Incidence'])
% % hold on
% % % globocan data
% % plot(ccIncYears , forouzanfar , '-' , ccIncYears , forouzanfar_ub , 'r--' , ccIncYears , forouzanfar_lb , 'r--')
% 
%% Cervical cancer incidence type distribution
% % newCCTotal = sum(sum(sum(newCC(: , : , : , :) , 2) , 3) , 4);
% % newCCType = zeros(size(newCC , 1) , 3);
% % for h = 2 : hpvVaxStates
% %     newCCType(: , h - 1) = sum(sum(newCC(: , : , h  , :) , 2) , 4) ...
% %         ./ newCCTotal;
% % end
% % figure(); area(tVec , newCCType)
% % legend('HPV 16/18' , 'Non-4v HPV' , 'Non-Vaccine HPV')
% % title('Cervical Cancer Incidence Type Distribution')
% 
%% cervical cancer incidence by age over time
% % ageTotal = zeros(length(tVec) , age);
% % for a = 1 : age
% %     ageTotal(: , a) = sum(popVec(1 : size(popVec , 1) , ...
% %             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %                 2 , a , 1 : risk))), 2) + sum(popVec(1 : size(popVec , 1) , ...
% %             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 8 : 10 , 1 : endpoints , ...
% %                 2 , a , 1 : risk))) , 2);
% % end
% % 
% % newCCAge = annlz((sum(sum(sum(newCC(1 : size(newCC , 1) , 1 : disease ,...
% %     1 : hpvVaxStates , 1 : age) , 2) , 3) , 4)));
% % 
% % ccIncEvo = bsxfun(@rdivide , permute(newCCAge , [1 , 4 , 3 , 2]) * fScale , ageTotal);
% % figure()
% % mesh(1 : age , tVec , ccIncEvo)
% % set(gca , 'yLim' , [tVec(1) tVec(end)]);
% % set(gca , 'xtick' , 1 : age , 'xtickLabel' , ageGroup);
% % xlabel('Age Group'); ylabel('Year'); zlabel('Incidence per 100,000')
% % title('Cervical Cancer Incidence')
% % %%
% % hold on
% % [px , py] = meshgrid(1 : age , 2017 ...
% %     .* ones(age , 1));
% % pz = bsxfun(@times , ones(size(px , 1) , size(py , 1)) , linspace(0 , max(ccIncEvo(:)) * 1.2 , size(px , 1)));
% % m = surf(px , py , pz' , 'edgecolor' , 'r');
% % set(m , 'facecolor' , 'r')
% % alpha(0.4)
% % %% cervical cancer incidence by age over time in HIV-
% % ageTotal = zeros(length(tVec) , age);
% % for a = 1 : age
% %     ageTotal(: , a) = sum(popVec(1 : size(popVec , 1) , ...
% %             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %                 2 , a , 1 : risk))), 2) + sum(popVec(1 : size(popVec , 1) , ...
% %             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 8 : 10 , 1 : endpoints , ...
% %                 2 , a , 1 : risk))) , 2);
% % end
% % 
% % newCCAge = (sum(sum(sum(newCC(1 : size(newCC , 1) , 1 ,...
% %     1 : hpvVaxStates , 1 : age) , 2) , 3) , 4));
% % 
% % ccIncEvo = bsxfun(@rdivide , permute(newCCAge , [1 , 4 , 3 , 2]) * fScale , ageTotal);
% % figure()
% % mesh(1 : age , tVec , ccIncEvo)
% % set(gca , 'yLim' , [tVec(1) tVec(end)]);
% % set(gca , 'xtick' , 1 : age , 'xtickLabel' , ageGroup);
% % xlabel('Age Group'); ylabel('Year'); zlabel('Incidence per 100,000')
% % title('HIV- Cervical Cancer Incidence')
% % %%
% % hold on
% % [px , py] = meshgrid(1 : age , 2017 ...
% %     .* ones(age , 1));
% % pz = bsxfun(@times , ones(size(px , 1) , size(py , 1)) , linspace(0 , max(ccIncEvo(:)) * 1.2 , size(px , 1)));
% % m = surf(px , py , pz' , 'edgecolor' , 'r');
% % set(m , 'facecolor' , 'r')
% % alpha(0.4)
% % % hold on
% % % [px , py] = meshgrid(1 : age , 2003 ...
% % %     .* ones(age , 1));
% % % pz = bsxfun(@times , ones(size(px , 1) , size(py , 1)) , linspace(0 , max(ccIncEvo(:)) * 1.2 , size(px , 1)));
% % % m = surf(px , py , pz' , 'edgecolor' , 'c');
% % % set(m , 'facecolor' , 'c')
% % % alpha(0.4)
% % %% Set up recording parameters (optional), and record
% % % OptionZ.FrameRate=15;OptionZ.Duration=10;OptionZ.Periodic=true;
% % % CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'CCInc',OptionZ)
% 
%% Compare simulation population size to actual population size
% % % Male (0 - 14)
% % m0_14all = zeros(length(tVec) , 1);
% % % Female(0 - 14)
% % f0_14all = zeros(length(tVec) , 1);
% % % Male (15 -64)
% % m15_64all = zeros(length(tVec) , 1);
% % % Female (15- 64)
% % f15_64all = zeros(length(tVec) , 1);
% % actualYrs = [1985 , 1996 , 2001 , 2011];
% % 
% % for i = 1 : length(tVec)
% %     m0_14Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , ...
% %         1 : hpvNonVaxStates , 1 : endpoints , 1 , 1 : 15 , 1 : risk));
% %     m0_14all(i) = sumall(popVec(i , m0_14Inds));
% % 
% %     f0_14Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , ...
% %         1 : hpvNonVaxStates , 1 : endpoints , 2 , 1 : 15 , 1 : risk));
% %     f0_14all(i) = sumall(popVec(i , f0_14Inds));
% % 
% %     m15_64Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , ...
% %         1 : hpvNonVaxStates , 1 :  endpoints , 1 , 16 : 65 , 1 : risk));
% %     m15_64all(i) = sumall(popVec(i , m15_64Inds));
% % 
% %     f15_64Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , ...
% %         1 : hpvNonVaxStates , 1 : endpoints , 2 , 16 : 65 , 1 : risk));
% %     f15_64all(i) = sumall(popVec(i , f15_64Inds));
% % end
% % 
% % % actual values
% % % males 0 - 14
% % actual(1 , :) = [1107298
% %     1526597
% %     1669704
% %     1658047];
% % 
% % % females 0 -14
% % actual(2 , :) = [1108238
% %     1537145
% %     1675104
% %     1621472];
% % 
% % % males 15 - 64
% % actual(3 , :) = [1549446
% %     2292625
% %     2659850
% %     3046456];
% % 
% % % females 15 - 64
% % actual(4 , :) = [1655333
% %     2711917
% %     3133521
% %     3433274];
% % 
% % actualYrs = [1985 , 1996 , 2001 , 2011];
% % popTotal = zeros(1 , size(popVec , 1));
% % for i = 1 : size(popVec , 1)
% %     popTotal(i) = sumall(popVec(i , :));
% % end
% % 
% % figure()
% % plot(tVec , popTotal , actualYrs , sum(actual , 1) , 'o')
% % xlabel('Year') ; ylabel('Population'); title('Population Size')
% % legend('Projected' , 'Actual')
% % figure()
% % plot(tVec , m0_14all , actualYrs , actual(1 , :) , 'o')
% % title('Males (0 - 14)')
% % xlabel('Year'); ylabel('Population Size')
% % legend('Projected' , 'Actual')
% % figure()
% % plot(tVec , f0_14all , actualYrs , actual(2 , :) , 'o')
% % title('Females (0 - 14)')
% % xlabel('Year'); ylabel('Population Size')
% % legend('Projected' , 'Actual')
% % figure()
% % plot(tVec , m15_64all , actualYrs , actual(3 , :) , 'o')
% % title('Males (15 - 64)')
% % xlabel('Year'); ylabel('Population Size')
% % legend('Projected' , 'Actual')
% % figure()
% % plot(tVec , f15_64all , actualYrs , actual(4 , :) , 'o')
% % title('Females (15 - 64)')
% % xlabel('Year'); ylabel('Population Size')
% % legend('Projected' , 'Actual')
% 
%% Population by age group over time
popByAge = zeros(length(tVec) , age);
%for i = 1 : length(tVec)
    for a = 1 : age
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
            1 : endpoints , 1: intervens, 1 : gender , a , 1 : risk));
        popByAge(a) = sum(popVec((2010 - startYear) * stepsPerYear , ageInds), 2 );
    end
%end
%
figure()
area(tVec(571), bsxfun(@rdivide , popByAge , sum(popVec(2020 - startYear) , 2)))
xlabel('Year'); ylabel('Relative Population Size'); title('Population')
legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')

%%
% figure()
% ages = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
%     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
%     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
% surf(1 : age , tVec , popByAge);
% set(gca , 'xtickLabel' , ages);
% set(gca , 'xLim' , [1 age]);
% set(gca , 'yLim' , [tVec(1) tVec(end)]);
% xlabel('Ages'); ylabel('Year'); zlabel('Population Size'); title('Population by Age')
% 
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
    ylim([0 15]);
end
 end
%legend('Male' , 'Female')
% title('HIV incidence')
% 
%% HPV incidence by gender
% % figure()
% % for g = 1 : 2
% %     hpvSusInds = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , ...
% %         1 : endpoints , g , 16 : 50 , 1 : risk));...
% %         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 10 , ...
% %         1 : endpoints , g , 16 : 50 , 1 : risk))];
% %     hpvSus = annlz(sum(popVec(1:end-1 , hpvSusInds) , 2)) ./ stepsPerYear;
% %     plot(tVec(1 : stepsPerYear : end-1) , annlz(sum(sum(sum(newHpv(1:end-1 , g , : , 16 : 50 , :) ...
% %         + newImmHpv(1:end-1 , g , : , 16 : 50 , :)...
% %         , 3) , 4) , 5)) ./ hpvSus * 100)
% %     axis([startYear , endYear , 0 , 100])    
% %     hold on
% % end
% % legend('Male' , 'Female')
% % xlabel('Year'); ylabel('Rate Per 100'); title('HPV Incidence')
% 
%% HIV incidence by age
% genderVec = {'Males' , 'Females'};
% hivAge = zeros(16 , 2);
% fScale = 100;
% ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
%     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
%     '60-64' , '65-69' , '70-74' , '75-79'};
% annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 
% 
% %figure;
% for g = 1 : gender
%     aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
%     for aInd = 1 : 16
%         a = aVec{aInd};
%         hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
%             1 : periods , g , a , 1 : risk)); ...
%             toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
%             1 : periods , g , a , 1 : risk))];
%         hivAge(aInd , g) = annlz(sum(sum(newHiv(end-6:end-1 , g , a , 1 : risk) , 3) , 4)) ...
%             ./ (annlz(sum(popVec(end-6:end-1 , hivSusInds) , 2)) ./ stepsPerYear) * fScale;
%     end
%     hold all;
%     subplot(1,2,g);
%     plot(1 : size(hivAge , 1) , hivAge(: , g) , '--');
%     hold all;
%     xlabel('Age Group'); ylabel('HIV incidence per 100')
%     set(gca , 'xtick' , 1 : length(hivAge) , 'xtickLabel' , ageGroup);
%     title(genderVec{g})
%     legend('Without ART dropout' , 'With ART dropout');
% end
% 
%% HPV incidence from immune only by gender
% % figure()
% % for g = 1 : 2
% %     hpvSusInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 10 , ...
% %         1 : endpoints , g , 16 : 50 , 1 : risk))];
% %     hpvSus = annlz(sum(popVec(1:end-1 , hpvSusInds) , 2)) ./ stepsPerYear;
% %     plot(tVec(1 : stepsPerYear : end-1) , annlz(sum(sum(sum(newImmHpv(1:end-1 , g , : , 16 : 50 , :)...
% %         , 3) , 4) , 5)) ./ hpvSus * 100)
% %     axis([startYear , endYear , 0 , 100])    
% %     hold on
% % end
% % legend('Male' , 'Female')
% % xlabel('Year'); ylabel('Rate Per 100'); title('HPV Incidence from Immune')
% 
%% HPV Incidence Among HIV-Negative
% % figure()
% % for g = 1 : 2
% %     hpvHivNegSusInds = [toInd(allcomb(1 , 1 : viral , 1 , 1 , ...
% %         1 : endpoints , g , 16 : 50 , 1 : risk));...
% %         toInd(allcomb(1 , 1 : viral , 2 : hpvVaxStates , 10 , ...
% %         1 : endpoints , g , 16 : 50 , 1 : risk));...
% %         toInd(allcomb(7 : 9 , 1 : viral , 1 , 1 , ...
% %         1 : endpoints , g , 16 : 50 , 1 : risk)); ...
% %         toInd(allcomb(7 : 9 , 1 : viral , 2 : hpvVaxStates , 10 , ...
% %         1 : endpoints , g , 16 : 50 , 1 : risk))];
% %     hpvSus = annlz(sum(popVec(1:end-1 , hpvHivNegSusInds) , 2)) ./ stepsPerYear;
% %     plot(tVec(1 : stepsPerYear : end-1) , annlz(sum(sum(sum(newHpv(1:end-1 , g , 1 , 16 : 50 , :) ...
% %         + newHpv(1:end-1 , g , 7 : 9 , 16 : 50 , :)...
% %         + newImmHpv(1:end-1 , g , 1 , 16 : 50 , :) ...
% %         + newImmHpv(1:end-1 , g , 7 : 9 , 16 : 50 , :)...
% %         , 3) , 4) , 5)) ./ hpvSus * 100)
% %     axis([startYear , endYear , 0 , 100])    
% %     hold on
% % end
% % legend('Male' , 'Female')
% % xlabel('Year'); ylabel('Rate Per 100'); title('HPV Incidence Among HIV-Negative')
% %  
%% HIV+ HPV incidence by age and gender
% % ages = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
% %     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
% %     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
% % aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
% % figure()
% % for g = 1 : 2
% %     for aInd = 1 : 16
% %         a = aVec{aInd};
% %         subplot(4 , 4 , aInd)
% %         hpvHivSusInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , g , a , 1 : risk)); ...
% %             toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvVaxStates , 10 , ...
% %             1 : endpoints , g , a , 1 : risk)); ...
% %              toInd(allcomb(10 , 6 , 1 , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , g , a , 1 : risk));...
% %             toInd(allcomb(10 , 6 , 2 : hpvVaxStates , 10 , ...
% %             1 : endpoints , g , a , 1 : risk))];
% %         hpvHivSus = annlz(sum(popVec(1:end-1 , hpvHivSusInds) , 2)) ./ stepsPerYear;
% %         plot(tVec(1 : stepsPerYear : end-1)' , ...
% %             annlz(sum(sum(sum(newHpv(1:end-1 , g , 2 : 6 , a , :) ...
% %             + newHpv(1:end-1 , g , 10 , a , :) ...
% %             + newImmHpv(1:end-1 , g , 2 : 6 , a , :) ...
% %             + newImmHpv(1:end-1 , g , 10 , a , :)...
% %             , 3), 4), 5)) ./ hpvHivSus * 100)
% %         hold on
% %         xlabel('Year'); ylabel('Rate Per 100'); title([' Age group ' , ages{aInd} , ' HPV Incidence in HIV+'])
% %         axis([startYear , endYear , 0 , 100])
% %     end
% % end
% % legend('Male' , 'Female')
% 
%% HIV+ HPV incidence by gender
% % figure()
% % for g = 1 : 2
% %     hpvHivSusInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
% %           1 : endpoints , g , 16 : 50 , 1 : risk)); ...
% %         toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvVaxStates , 10 , ...
% %           1 : endpoints , g , 16 : 50 , 1 : risk)); ...
% %         toInd(allcomb(10 , 6 , 1 , 1 : hpvNonVaxStates , ...
% %           1 : endpoints , g , 16 : 50 , 1 : risk));...
% %         toInd(allcomb(10 , 6 , 2 : hpvVaxStates , 10 , ...
% %           1 : endpoints , g , 16 : 50 , 1 : risk))];
% %     hpvHivSus = annlz(sum(popVec(1:end-1 , hpvHivSusInds) , 2)) ./ stepsPerYear;
% %     plot(tVec(1 : stepsPerYear : end-1) , ...
% %         annlz(sum(sum(sum(newHpv(1:end-1 , g , 2 : 6 , 16 : 50 , :) ...
% %         + newHpv(1:end-1 , g , 10 , 16 : 50 , :) ...
% %         + newImmHpv(1:end-1 , g , 2 : 6 , 16 : 50 , :) ...
% %         + newImmHpv(1:end-1 , g , 10 , 16 : 50 , :)...
% %         , 3), 4), 5)) ./ hpvHivSus * 100)
% %     hold on
% %     xlabel('Year'); ylabel('Rate Per 100'); title([' HPV Incidence in HIV+'])
% %     axis([startYear , endYear , 0 , 100])
% % end
% % legend('Male' , 'Female')
% 
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
%% CC incidence for ages 15+
% % inds = {':' , [2 : 6 , 10] , [2 : 6] , 1 , 10};
% % files = {'CEA CC_General_Hpv' , 'CEA CC_HivAll_Hpv' , ...
% %      'CEA CC_HivNoART_Hpv' , 'CEA CC_HivNeg_Hpv' ,...
% %      'CEA CC_ART_HPV'};
% % plotTits = {'General' , 'HIV-Positive' , 'HIV-Positive (No ART)' , ....
% %      'HIV-Negative' , 'HIV-Positive on ART'};
% % fac = 10 ^ 5;
% % 
% % for i = 1 : length(inds)
% %         % general
% %         allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
% %             1 : endpoints , 2 , 16 : age , 1 : risk)); ...
% %             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
% %             1 : endpoints , 2 , 16 : age , 1 : risk))];
% %         % All HIV-positive women
% %         allHivF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
% %             1 : endpoints , 2 , 16 : age , 1 : risk)); ...
% %             toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
% %             1 : endpoints , 2 , 16 : age , 1 : risk));...
% %             toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : 4 , ...
% %             1 : endpoints , 2 , 16 : age , 1 : risk)); ...
% %             toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 : 10 , ...
% %             1 : endpoints , 2 , 16 : age , 1 : risk))];
% %         % HIV-positive women not on ART
% %         hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
% %             1 : endpoints , 2 , 16 : age , 1 : risk)); ...
% %             toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
% %             1 : endpoints , 2 , 16 : age , 1 : risk))];
% %         % All HIV-negative women
% %         hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %             2 , 16 : age , 1 : risk)); ...
% %             toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : endpoints , ...
% %             2 , 16 : age , 1 : risk))];
% %         % Women on ART
% %         artF = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : 4 , ...
% %             1 : endpoints , 2 , 16 : age , 1 : risk)); ...
% %             toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 : 10 , ...
% %             1 : endpoints , 2 , 16 : age , 1 : risk))];
% %         
% %         genArray = {allF , allHivF , hivNoARTF , hivNeg , artF};
% %         
% %         ccInc = ...
% %             annlz(sum(sum(sum(newCC(1:end-1 , inds{i} , : , 16 : age),2),3),4)) ./ ...
% %             (annlz(sum(popVec(1:end-1 , genArray{i}) , 2) ./ stepsPerYear))* fac;
% %         
% %         figure()
% %         
% %         plot(tVec(1 : stepsPerYear : end-1) , ccInc)
% %         
% %         title([plotTits{i} , ' Cervical Cancer Incidence'])
% %         xlabel('Year'); ylabel('Incidence per 100,000')
% % end
% 
%% CC incidence standardized for ages 10+
% % inds = {':' , [2 : 6] , [1,7:9] , 10};
% % files = {'CC_General_Hpv_VaxCover' , ...
% %      'CC_HivNoART_Hpv_VaxCover' , 'CC_HivNeg_Hpv_VaxCover' ,...
% %      'CC_ART_HPV_VaxCover'};
% % plotTits = {'General' , 'HIV-Positive (No ART)' , ....
% %      'HIV-Negative' , 'HIV-Positive on ART'};
% % fac = 10 ^ 5;
% % 
% % figure();
% % for i = 2 : length(inds)
% %     ccIncRef = zeros(length(tVec(1 : stepsPerYear : end-1)),1)';
% %     
% %     % General, all ages
% %     allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 3 , ...
% %         1 : endpoints , 2 , 11 : age , 1 : risk)); ...
% %         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
% %         1 : endpoints , 2 , 11 : age , 1 : risk))];
% %     allhivNegFAge = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : 3 , 1 : endpoints , ...
% %             2 , 11 : age , 1 : risk)); ...
% %             toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : endpoints , ...
% %             2 , 11 : age , 1 : risk))];
% %     
% %     for a = 11 : age
% %         % General
% %         allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
% %             1 : endpoints , 2 , a , 1 : risk)); ...
% %             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
% %             1 : endpoints , 2 , a , 1 : risk))];
% %         % HIV-positive women not on ART
% %         hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
% %             1 : endpoints , 2 , a , 1 : risk)); ...
% %             toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
% %             1 : endpoints , 2 , a , 1 : risk))];
% %         % All HIV-negative women
% %         hivNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %             2 , a, 1 : risk)); ...
% %             toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : endpoints , ...
% %             2 , a , 1 : risk))];
% %         % Women on ART
% %         artF = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : 4 , ...
% %             1 : endpoints , 2 , a , 1 : risk)); ...
% %             toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 : 10 , ...
% %             1 : endpoints , 2 , a , 1 : risk))];
% %         genArray = {allF , hivNoARTF , hivNeg , artF};
% % 
% %         ccIncRefTemp = ...
% %             (annlz(sum(sum(newCC(1:end-1 , inds{i} , : , a),2),3)) ./ ...
% %             (annlz(sum(popVec(1:end-1 , genArray{i}) , 2) ./ stepsPerYear))* fac) ...
% %             .* (annlz(sum(popVec(1:end-1 , genArray{1}) , 2) ./ stepsPerYear));
% % 
% %         ccIncRef = ccIncRef + ccIncRefTemp;
% % 
% %     end
% %     
% %     %annlz(sum(sum(sum(newCC(: , inds{i} , : , 3:age),2),3),4))
% %     %annlz(sum(popVec(: , genArray{i}) , 2) ./ stepsPerYear)
% %     
% %     ccInc = ccIncRef ./ (annlz(sum(popVec(1:end-1 , allFAge) , 2) ./ stepsPerYear));
% %     plot(tVec(1 : stepsPerYear : end-1) , ccInc)
% %     legend('HivNoART' , 'HivNeg' , 'HivART');
% %     hold all;
% %     
% %     title([plotTits{i} , ' Cervical Cancer Incidence'])
% %     xlabel('Year'); ylabel('Incidence per 100,000')
% % end       
% 
%% General CC incidence validation
% % fac = 10 ^ 5;
% % % ccInc = annlz(sum(sum(sum(newCC(: , : , : , 4 : age),2),3),4)) ./ ...
% % %     (annlz(sum(popVec(: , allF) , 2) ./ stepsPerYear))* fac;
% % 
% % worldStandard_Segi1964 = [12000 10000 9000 9000 8000 8000 6000 6000 6000 ...
% %     6000 5000 4000 4000 3000 2000 1000 500 500];
% % 
% % % General, all ages
% % % allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
% % %     1 : endpoints , 2 , 16 : age , 1 : risk)); ...
% % %     toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
% % %     1 : endpoints , 2 , 16 : age , 1 : risk))];
% % ccIncRef = zeros(length(tVec(1 : stepsPerYear : end-1)),1)';
% % 
% % aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
% % for aInd = 4 : 16
% %     a = aVec{aInd};
% %     % General
% %     allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
% %         1 : endpoints , 2 , a , 1 : risk)); ...
% %         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
% %         1 : endpoints , 2 , a , 1 : risk))];
% % 
% %     ccIncRefAge = ...
% %         (annlz(sum(sum(sum(newCC(1:end-1 , : , : , a),2),3),4)) ./ ...
% %         (annlz(sum(popVec(1:end-1 , allF) , 2) ./ stepsPerYear))* fac) ...
% %         .* (worldStandard_Segi1964(aInd));
% %     ccIncRef = ccIncRef + ccIncRefAge; 
% %       
% % end
% % ccInc = ccIncRef ./ (sum(worldStandard_Segi1964(4:16)));
% % 
% % globocan = ... %[1.0	3.4	17.1	20.2	34.0	59.4	50.4	104.1	87.0	96.7	57.3	70.0 70.0];	
% % [0
% % 2.646467154
% % 8.848389036
% % 45.1937379
% % 53.40682334
% % 63.4
% % 68.3
% % 70.7
% % 73
% % 77.4
% % 82.7
% % 88.6
% % 95.2];
% % 
% % ccIncRef = 0;
% % for aInd = 4 : 16
% %     a = aVec{aInd};
% %     ccIncRefAge = globocan(aInd-3) .* (worldStandard_Segi1964(aInd));
% %     ccIncRef = ccIncRef + ccIncRefAge; 
% % end
% % ccIncGlobocan = ccIncRef ./ (sum(worldStandard_Segi1964(4:16)));
% % 
% % olorunfemi = [1994.0648457561042, 22.241027817219138;
% %     1994.4057758035783, 22.48378323297575;
% %     1994.7466755440043, 22.735208485009384;
% %     1995.08512041352, 23.68889047548178;
% %     1995.4213831755596, 25.266800677899745;
% %     1995.7576762446477, 26.83604104404069;
% %     1996.0971212467564, 27.503618437371372;
% %     1996.4385058999546, 27.616326308972656;
% %     1996.7799208602014, 27.720364344296918;
% %     1997.1203659949033, 28.101837140485877;
% %     1997.4603868309293, 28.604687644553138;
% %     1997.800589509245, 29.055519130958274;
% %     1998.1458534646226, 28.058487959100766;
% %     1998.4929964569938, 26.52392693806791;
% %     1998.8399576070751, 25.04138493469718;
% %     1999.1845245003426, 24.243759997211175;
% %     1999.5283034103547, 23.671550802927737;
% %     1999.8720217062705, 23.116681281198343;
% %     2000.2147398695931, 22.847916356610668;
% %     2000.5572458835777, 22.639840285962144;
% %     2000.8996306693693, 22.466443560421705;
% %     2001.2412274719056, 22.518462578083838;
% %     2001.5827333532973, 22.596491104577034;
% %     2001.9247847615577, 22.518462578083838;
% %     2002.2689576631976, 21.833545512199112;
% %     2002.6134336353202, 21.061930083544173;
% %     2002.956818553705, 20.60242876086202;
% %     2003.2965060121999, 21.200647463976523;
% %     2003.6358297861154, 21.90290420241529;
% %     2003.9758203150932, 22.414424542759576;
% %     2004.3175080387743, 22.440434051590643;
% %     2004.6592866836004, 22.440434051590643;
% %     2005.0006410297506, 22.561811759468945;
% %     2005.3409649362593, 22.977963900765992;
% %     2005.681288842768, 23.39411604206304;
% %     2006.0226734959665, 23.506823913664324;
% %     2006.3654522733857, 23.220719316522604;
% %     2006.7083219719495, 22.90860521054982;
% %     2007.0505552224997, 22.77855766639449;
% %     2007.3919398756982, 22.891265537995775;
% %     2007.7333851429933, 22.986633737043014;
% %     2008.0747697961915, 23.0993416086443;
% %     2008.4162150634866, 23.194709807691538;
% %     2008.7576603307816, 23.290078006738778];
% % 
% % figure()
% % plot(tVec(1 : stepsPerYear : end-1) , ccInc)
% % hold on;
% % scatter(olorunfemi(:,1),olorunfemi(:,2))
% % hold all;
% % scatter(2012,ccIncGlobocan , '*')
% % hold all;
% % scatter(2018, 61.9 , '*') % Globocan 2018 South Africa incidence rates were estimated from ...
% %                     % national mortality estimates by modelling, using mortality-to-incidence ratios ...
% %                     % derived from cancer registries in neighbouring countries: ASR = 61.9/100,000 for females 15-79
% % xlim([1995 2020]);
% % title('General Cervical Cancer Incidence')
% % xlabel('Year'); ylabel('Incidence per 100,000')
% % legend('Model' , 'Olorunfemi Validation' , 'Globocan 2012 Validation' , 'Globocan 2018 Validation')
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
%% susceptibles
% % ages = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
% %     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
% %     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
% % 
% % % figure()
% % % for g = 1 : 2
% % %     for a = 1 : age
% % %         subplot(4 , 4 , a)
% % %         hpvHivSusInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
% % %             1 : endpoints , g , a , 1 : risk)); ...
% % %             toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvVaxStates , 10 , ...
% % %             1 : endpoints , g , a , 1 : risk)); ...
% % %             toInd(allcomb(10 , 6 , 1 , 1 : hpvNonVaxStates , ...
% % %             1 : endpoints , g , a , 1 : risk));...
% % %             toInd(allcomb(10 , 6 , 2 : hpvVaxStates , 10 , ...
% % %             1 : endpoints , g , a , 1 : risk))];
% % %         ageInd = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates, 1 : hpvNonVaxStates , ...
% % %             1 : endpoints , g , a , 1 : risk)); ...
% % %             toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% % %             1 : endpoints , g , a , 1 : risk))];
% % %         plot(tVec , sum(popVec(: , hpvHivSusInds) , 2) ./ sum(popVec(: , ageInd) , 2) * 100)
% % %         hold on
% % %         xlabel('Year'); ylabel('%'); title([' Age group ' , ages{a} , ' HPV Susceptible HIV+'])
% % %         axis([startYear , endYear , 0 , 100])
% % %     end
% % % end
% % % legend('Male' , 'Female')
% % 
% % figure()
% % for g = 2
% %     for a = 1 : age
% %         subplot(4 , 4 , a)
% %         hpvHivSusInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , g , a , 1 : risk)); ...
% %             toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvVaxStates , 10 , ...
% %             1 : endpoints , g , a , 1 : risk))];
% %         ageInd = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates, 1 : hpvNonVaxStates , ...
% %             1 : endpoints , g , a , 1 : risk))];
% %         plot(tVec , sum(popVec(: , hpvHivSusInds) , 2) ./ sum(popVec(: , ageInd) , 2) * 100)
% %         hold on
% %         xlabel('Year'); ylabel('%'); title([' Age group ' , ages{a} , ' HPV Susceptible HIV+, no ART: Females'])
% %         axis([startYear , endYear , 0 , 100])
% %     end
% % end
% % % legend('Male' , 'Female')
% % legend('HIV+,noART','HIV+,ART')
% % hold all;
% % 
% % % figure()
% % for g = 2
% %     for a = 1 : age
% %         subplot(4 , 4 , a)
% %         hpvHivSusInds = [toInd(allcomb(10 , 6 , 1 , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , g , a , 1 : risk));...
% %             toInd(allcomb(10 , 6 , 2 : hpvVaxStates , 10 , ...
% %             1 : endpoints , g , a , 1 : risk))];
% %         ageInd = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , g , a , 1 : risk))];
% %         plot(tVec , sum(popVec(: , hpvHivSusInds) , 2) ./ sum(popVec(: , ageInd) , 2) * 100)
% %         hold on
% %         xlabel('Year'); ylabel('%'); title([' Age group ' , ages{a} , ' HPV Susceptible HIV+, ART: Females'])
% %         axis([startYear , endYear , 0 , 100])
% %     end
% % end
% % % legend('Male' , 'Female')
% % 
% % % figure()
% % % for g = 1 : 2
% % %     for a = 1 : age
% % %         subplot(4 , 4 , a)
% % %         hpvHivSusInds = [toInd(allcomb([1,7:9] , 1 , 1 , 1 : hpvNonVaxStates , ...
% % %             1 : endpoints , g , a , 1 : risk)); ...
% % %             toInd(allcomb([1,7:9] , 1 : viral , 2 : hpvVaxStates , 10 , ...
% % %             1 : endpoints , g , a , 1 : risk))];
% % %         ageInd = [toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates, 1 : hpvNonVaxStates , ...
% % %             1 : endpoints , g , a , 1 : risk))];
% % %         plot(tVec , sum(popVec(: , hpvHivSusInds) , 2) ./ sum(popVec(: , ageInd) , 2) * 100)
% % %         hold on
% % %         xlabel('Year'); ylabel('%'); title([' Age group ' , ages{a} , ' HPV Susceptible HIV-'])
% % %         axis([startYear , endYear , 0 , 100])
% % %     end
% % % end
% % % legend('Male' , 'Female')
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
%% Population Size
% % % HIV-positive women not on ART
% % hivNoART = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
% %     1 : endpoints , 1 : gender , 3 : age , 1 : risk)); ...
% %     toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
% %     1 : endpoints , 1 : gender , 3 : age , 1 : risk))];
% % % All HIV-negative women
% % hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %     1 : gender , 3 : age , 1 : risk)); ...
% %     toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : endpoints , ...
% %     1 : gender , 3 : age , 1 : risk))];
% % % Women on ART
% % art = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : 4 , ...
% %     1 : endpoints , 1 : gender , 3 : age , 1 : risk)); ...
% %     toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 : 10 , ...
% %     1 : endpoints , 1 : gender , 3 : age , 1 : risk))];
% % genArray = {hivNoART , hivNeg , art};
% % 
% % figure()
% % for i = 1 : length(genArray)
% %     plot(tVec , sum(popVec(: , genArray{i}) , 2))
% %     hold all;
% % end
% % title('Population Size')
% % xlabel('Year'); ylabel('Individuals')
% % xlim([1910 2099]);
% % legend('HIV+ , no ART' , 'HIV-' , 'HIV+ , ART');
% % hold off
% 
%% Population by risk and HIV group
% % figure();
% % for r = 1 : risk
% %     % HIV+
% %     rpopHivTot = popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         1 , 3 , r)));
% %     popHivTot = popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         1 , 3 , 1 : risk)));
% %     %ART
% %     rpopArtTot = popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         1 , 3 , r)));
% %     popArtTot = popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         1 , 3 , 1 : risk)));
% %     %HIV-
% %     rpopHivNegTot = popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         1 , 3 , r)));
% %     popHivNegTot = popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         1 , 3 , 1 : risk)));
% % 
% %     % plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
% %     plot(tVec , 100 * sum(rpopHivNegTot , 2) ./ sum(popHivNegTot , 2))
% %     hold all
% %     plot(tVec , 100 * sum(rpopHivTot , 2) ./ sum(popHivTot , 2))
% %     hold all
% %     plot(tVec , 100 * sum(rpopArtTot , 2) ./ sum(popArtTot , 2))
% %     xlabel('Year'); ylabel('Prevalence (%)'); title(' Risk Proportion')
% %     legend('HIV- : lr' , 'HIV+ noART : lr' , 'ART : lr' , 'HIV- : mr' , 'HIV+ noART : mr' , 'ART : mr' , 'HIV- : hr' , 'HIV+ noART : hr'  , 'ART : hr')
% %     %legend('HIV+ noART : lr' ,  'HIV+ noART : mr' , 'HIV+ noART : hr')
% %     hold all;
% % end

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

%% HIV by age and risk
% % hivAgeRisk = zeros(length(tVec) , age , risk);
% % for a = 1 : age
% %     for r = 1 : risk
% %         hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , 1 : gender , a , r));
% %         hivArt = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , 1 : gender , a , r));
% %         hivAgeRisk(: , a , r) = sum(popVec(: , hivPos) , 2) + sum(popVec(: , hivArt) , 2);
% %     end
% % end
% 
%% ART treatment tracker- cd4
% figure()
% cd4ARTFrac = zeros(length(tVec) , 5);
% for a = 1 : 5 : age
%     for i = 1 : length(tVec)
%     currTot = sumall(artTreatTracker(i , 2 : 6 , 1 : 5 , 1 : gender , a , 1 :risk));
%         for d = 2 : 6
%             curr = sumall(artTreatTracker(i , d , 1 : 5 , 1 : gender , a , 1 : risk));
%             cd4ARTFrac(i , d - 1) = curr / currTot;
%         end
%     end
%     subplot(4,4,a)
%     area(tVec , cd4ARTFrac)
%     xlabel('Year')
%     ylabel('Initiated ART')
% end
% legend('Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
%         'CD4 <= 200 cells/uL' , 'Location' , 'NorthEastOutside')
% 
%% ART treatment tracker- risk/age
% % aARTFrac = zeros(length(tVec) , 5);
% % for i = 1 : length(tVec)
% %     currTot = sumall(artTreatTracker(i , 2 : 6 , 1 : 5 , 2 , 1 : age , 1 : risk));
% %     for a = 1 : age
% %         curr = sumall(artTreatTracker(i , 2 : 6 , 1 : 5 , 2 , a , 1 : risk));
% %         aARTFrac(i , a) = curr / currTot;
% %     end
% % end
% % 
% % figure()
% % area(tVec , aARTFrac)
% % legend('1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , '11' ,'12' , '13' , '14' , '15' ,'16' , 'Location' , 'NorthEastOutside')
% % % legend('lr','mr','hr');
% % xlabel('Year')
% % ylabel('Initiated ART')
% 
%% CD4 vs VL
% % figure()
% % cd4Vl = zeros(length(tVec) , 6 , viral);
% % dHivPos = [2 : 6 , 10];
% % f = zeros(length(tVec) , 1);
% % set(gca, 'nextplot','replacechildren', 'Visible','off');
% % % mov(1:f) = struct('cdata',[], 'colormap',[]);
% % % set(0,'DefaultFigureVisible','off')
% % vid = VideoWriter([date '_cd4_vs_vl_heatmap']);
% % cd4 = {'Acute infection' , 'CD4 > 500' , 'CD4 500-350' , ...
% %     'CD4 350-200' , 'CD4 <= 200' , 'HIV-positive on ART'};
% % vl = {'Acute infection' , 'VL < 1000' , 'VL 1000-10,000', 'VL 10,000-50,000' , ...
% %     'VL > 50,000' , 'HIV-positive on ART'};
% % vlHeatMap = flip(vl);
% % for i = 1 : length(tVec)
% %     for j = 1 : length(dHivPos)
% %         d = dHivPos(j);
% %         for v = 1: viral
% %             curr = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates ,...
% %                 1 : endpoints , 1 : gender , 1 : age , 1 : risk));
% %             cd4Vl(i , j , v) = sumall(popVec(i , curr));
% %         end
% %     end
% % end
% 
%% HPV incidence standardized
% % inds = {':' , [2 : 6] , [1,7:9] , 10};
% % files = {'CC_General_Hpv_VaxCover' , ...
% %      'CC_HivNoART_Hpv_VaxCover' , 'CC_HivNeg_Hpv_VaxCover' ,...
% %      'CC_ART_HPV_VaxCover'};
% % plotTits = {'General' , 'HIV-Positive (No ART)' , ....
% %      'HIV-Negative' , 'HIV-Positive on ART'};
% % fac = 10 ^ 5;
% % 
% % figure();
% % 
% % for i = 2 : length(inds)
% % %     figure();
% %    hpvIncRef = zeros(length(tVec(1 : stepsPerYear : end-1)),1)';
% %     
% %     % General, all ages
% %     allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , ...
% %         1 : endpoints , 2 , 11 : age , 1 : risk)); ...
% %         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 10 , ...
% %         1 : endpoints , 2 , 11 : age , 1 : risk))];
% %     allhivNegFAge = [toInd(allcomb([1,7:9] , 1 : viral , 1 , 1 , 1 : endpoints , ...
% %             2 , 11 : age , 1 : risk)); ...
% %             toInd(allcomb([1,7:9] , 1 : viral , 1 , 10 , 1 : endpoints , ...
% %             2 , 11 : age , 1 : risk))];
% %     
% %     for a = 11 : age
% %         % General
% %         allF = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , ...
% %             1 : endpoints , 2 , a , 1 : risk)); ...
% %             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 10 , ...
% %             1 : endpoints , 2 , a , 1 : risk))];
% %         % HIV-positive women not on ART
% %         hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 , ...
% %             1 : endpoints , 2 , a , 1 : risk)); ...
% %             toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 10 , ...
% %             1 : endpoints , 2 , a , 1 : risk))];
% %         % All HIV-negative women
% %         hivNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 , 1 , 1 : endpoints , ...
% %             2 , a, 1 : risk)); ...
% %             toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 10 , 1 : endpoints , ...
% %             2 , a , 1 : risk))];
% %         % Women on ART
% %         artF = [toInd(allcomb(10 , 6 , 1 , 1 , ...
% %             1 : endpoints , 2 , a , 1 : risk)); ...
% %             toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 10 , ...
% %             1 : endpoints , 2 , a , 1 : risk))];
% %         genArray = {allF , hivNoARTF , hivNeg , artF};
% % 
% %         hpvIncRef = ...
% %             ((annlz(sum(sum(newHpv(1:end-1 , 2 , inds{i} , a , :),3),5))+annlz(sum(sum(newImmHpv(1:end-1 , 2 , inds{i} , a , :),3),5))) ./ ...
% %             (annlz(sum(popVec(1:end-1 , genArray{i}) , 2) ./ stepsPerYear))* fac) ...
% %             .* (annlz(sum(popVec(1:end-1 , genArray{3}) , 2) ./ stepsPerYear));
% %         hpvIncRef = hpvIncRef + hpvIncRef;
% %         
% %     end
% %     hpvInc = hpvIncRef ./ (annlz(sum(popVec(1:end-1 , allhivNegFAge) , 2) ./ stepsPerYear));
% %     plot(tVec(1 : stepsPerYear : end-1) , hpvInc ,'DisplayName' , ...
% %          plotTits{i});
% %     legend('-DynamicLegend');
% %     hold all;
% %     title(' HPV Incidence')
% %     xlabel('Year'); ylabel('Incidence per 100,000')
% %     hold all;
% % end   
% 
%% Population by "p"
% % figure();
% % for d = 1 : disease
% %     for p = 1 : endpoints
% %         subplot(3,4,d);
% %         % General
% %         inds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %              p , 1 : gender , 36:40 , 1 : risk));
% %         pop = sum(popVec(: , inds) , 2);
% %         popTot = popVec(: , toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %              1 : endpoints , 1 : gender , 36:40 , 1 : risk)));
% %         plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
% %         xlabel('Year'); ylabel('Proportion (%)'); title(' p Proportion')
% %         legend('1' , '2' , '3' , '4' , '5' , '6')
% %         hold all;
% %     end
% % end
% 
%% Screened proportion by HIV group
% % figure();
% % linStyle = {'--' , '-' , ':'};
% % for a = 36
% %     for r = 1 : risk
% %     % HIV+
% %     vaxHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 4 : 6 , 2 , a , r));
% %     vaxHivPop = sum(popVec(: , vaxHivInds) , 2);
% %     popHivTot = popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         2 , a , r)));
% %     %ART
% %     vaxArtInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 4 : 6 , 2 , a , r));
% %     vaxArtPop = sum(popVec(: , vaxArtInds) , 2);
% %     popArtTot = popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         2 , a , r)));
% %     %HIV-
% %     vaxHivNegInds = toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 4 : 6 , 2 , a , r));
% %     vaxHivNegPop = sum(popVec(: , vaxHivNegInds) , 2);
% %     popHivNegTot = popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
% %         2 , a , r)));
% % 
% % %     subplot(4,4,a)
% %     plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2) , linStyle{r})
% %     hold all
% %     plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2) , linStyle{r})
% %     hold all
% %     plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2) , linStyle{r})
% %     xlabel('Year'); ylabel('Proportion (%)'); title('Screened Proportion')
% %     
% %     hold all;
% %     end
% % end
% % legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')
% 
%%
% % reset(gca)
% % figure()
% % open(vid);
% % x = 1 : 6;
% % y = 1 : 6;
% % % [x , y] = meshgrid(x, y);
% % % xq = 1 : 0.1 : 6;
% % % yq = 1 : 0.1 : 6;
% % % [xq , yq] = meshgrid(xq , yq);
% % % for i = 1 : length(tVec)
% % %     zq = griddata(x , y , squeeze(cd4Vl(i , : , :)) , xq , yq , 'nearest');
% % %     mesh(xq , yq , zq)
% % %     set(gca , 'xtick' , 1 : 6);
% % %     set(gca, 'XTick', 1 : 6, 'XTickLabel', cd4)
% % %
% % %     set(gca , 'ytick' , 1 : 6);
% % %     set(gca, 'YTick', 1 : 6, 'YTickLabel', vl)
% % %
% % %     set(gca , 'yTickLabelRotation' , 15)
% % % %     set(gca , 'zLim' , [0 max(max(cd4Vl(i , : , :))) + 100]);
% % %     set(gca , 'zLim' , [0 max(cd4Vl(:))]);
% % %     title(['Year : ' , num2str(floor(tVec(i)))])
% % %     writeVideo(vid , getframe(gcf));
% % % end
% %
% % for i = 1 : length(tVec)
% %     %     %     s = stem3(1 : 6 , 1 : 6 , squeeze(cd4Vl(i , : , :)) , '--');
% %     %     %     s.Color = 'red';
% %     %     %     s.MarkerFaceColor = 'red';
% %     %     %     s.MarkerSize = 8;
% %     %     %     set(gca , 'xtick' , 6 : 1);
% %     %     %     set(gca , 'ytick' , 1 : 6);
% %     %     %     set(gca , 'yTickLabelRotation' , 15)
% %     %     %     set(gca , 'zLim' , [1 max(cd4Vl(:))]);
% %     colormap('hot')
% %     grid on
% %     caxis([1 max(cd4Vl(:))])
% %     colorbar
% %     imagesc(flipud(squeeze(cd4Vl(i , : , :))))
% %     set(gca, 'XTick', 1 : 6, 'XTickLabel', cd4)
% %     set(gca , 'XTickLabelRotation' , 90)
% %     set(gca, 'YTick', 1 : 6, 'YTickLabel', vlHeatMap)
% %     title(['Year : ' , num2str(floor(tVec(i)))])
% %     writeVideo(vid , getframe(gcf));
% % end
% % close(gcf);
% % close(vid);
% % winopen([date '_cd4_vs_vl_heatmap.avi'])
% % set(0,'DefaultFigureVisible','on')
% % reset(gca)
% %%
% % for f = 1 : numel(figHandles)
% %     %   saveas(figHandles(f),sprintf('figure_%d.jpg',f))
% %     baseFileName = sprintf('figure_%d.jpg',f);
% %     % Specify some particular, specific folder:
% %     fullFileName = fullfile('C:\Users\nicktzr\Google Drive\ICRC\CISNET\Model\Figures', baseFileName);
% %     figure(f); % Activate the figure again.
% %     export_fig(fullFileName); % Using export_fig instead of saveas.
% % end
% %%
% % resultOut()
