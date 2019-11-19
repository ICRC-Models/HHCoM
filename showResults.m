function[] = showResults(pathModifier)

%% Load parameters and results
paramDir = [pwd , '\Params\'];

% Load saved parameters
% load([paramDir, 'generalParams'], 'stepsPerYear' , 'timeStep' , ...
%     'disease' , 'viral' , 'hpvVaxStates' , 'hpvNonVaxStates' , 'endpoints' , ...
%     'intervens' , 'gender' , 'age' , 'risk' , 'hpvTypeGroups' , 'dim' , 'k' , 'toInd' , ...
%     'sumall');
% load([paramDir, 'demoBehavParams'], 'ageSexDebut' , 'mInit' , 'fInit' , 'partnersM' , 'partnersF' , ...
%     'maleActs' , 'femaleActs' , 'riskDist' , 'mue' , 'epsA_vec' , 'epsR_vec' , 'yr');
% load([paramDir, 'hivParams'], 'betaHIVM2F' , 'betaHIVF2M' , 'muHIV' , 'kVl' , 'kCD4');
% load([paramDir, 'hpvParams'], 'perPartnerHpv_vax' , 'perPartnerHpv_nonV' , ...
%     'fImm' , 'rImmune' , 'kCin1_Inf' , 'kCin2_Cin1' , 'kCin3_Cin2' , 'kCC_Cin3' , ...
%     'rNormal_Inf' , 'kInf_Cin1' , 'kCin1_Cin2' , 'kCin2_Cin3' , 'lambdaMultImm' , ...
%     'hpv_hivClear' , 'rImmuneHiv' , 'c3c2Mults' , 'c2c1Mults' , 'muCC' , ...
%     'kRL' , 'kDR' , 'artHpvMult' , 'hpv_hivMult');
% load([paramDir, 'intervenParams'], 'circ' , 'condUse' , ...
%     'maxRateM1' , 'maxRateF1' , 'maxRateM2' , 'maxRateF2' , 'hivStartYear' , 'circStartYear' , ...
%     'vaxStartYear' , 'baseline' , 'circProtect' , 'condProtect' , 'MTCTRate');
% load([paramDir , 'calibData'], 'cinPos2002_dObs' , 'cinNeg2002_dObs' , ...
%     'hpv_hiv_dObs' , 'hpv_hivNeg_dObs' , 'hpv_hivM2008_dObs' , 'hpv_hivMNeg2008_dObs' , ...
%     'hivPrevM_dObs' , 'hivPrevF_dObs');
[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd , ...
    sumall , annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , mue , epsA_vec , epsR_vec , ...
    yr , ...
    hivOn , betaHIVM2F , betaHIVF2M , muHIV , kVl , kCD4 , ...
    hpvOn , perPartnerHpv_vax , perPartnerHpv_nonV , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , ...
    circ , condUse , screenYrs , hpvScreenStartYear , waning , ...
    maxRateM1 , maxRateF1 , maxRateM2 , maxRateF2 , hivStartYear , ...
    circStartYear , vaxStartYear , baseline , cisnet , who , whob , ...
    circProtect , condProtect , MTCTRate , hyst , ...
    OMEGA , ...
    ccInc2011_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , ...
    hivPrevM_dObs , hivPrevF_dObs , ...
    mCurr , fCurr , mCurrArt , fCurrArt , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , ...
    hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ...
    ageInd , riskInd , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , deathMat , circMat , circMat2] = loadUp2(1 , 0 , [] , [] , []);

resultsDir = [pwd , '\HHCoM_Results\'];
load([resultsDir , pathModifier])
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 

% Plot Settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)

%% HIV prevalence vs. ANC data
% Total HIV positive
hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 : 2 , 4 : 10 , 1 : risk));
hivPop = sum(popVec(: , hivInds) , 2);
artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 : 2 , 4 : 10 , 1 : risk));
art = sum(popVec(: , artInds) , 2);
% Compared to ANC data
overallHivPrev_KZN_AC(1 , :) = 1990 : 2009;
overallHivPrev_KZN_AC(2 , :) = [0.464072571
    0.985438052
    1.506803533
    5.576509907
    8.126044402
    13.04177608
    12.54905705
    16.61876343
    19.50632609
    19.52064932
    22.2391979
    20.22439723
    22.09787539
    22.78825495
    25.16877536
    26.19622822
    25.36548102
    27.2380043
    27.42134161
    28.44974934];

% Compared to Africa Center (validation years)
prevValYrs = 2010 : 2016;
prevVal = [0.290
0.290
0.293
0.312
0.338
0.338
0.344
] .* 100;

upper_prevVal = [0.30
0.30
0.31
0.32
0.35
0.35
0.36
] .* 100;

lower_prevVal = [0.27
0.27
0.27
0.29
0.32
0.32
0.33] .* 100;

figure()
popTot = popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
    1 : intervens , 1 : 2 , 4 : 10 , 1 : risk)));
plot(tVec , (hivPop + art) ./ sum(popTot , 2) * 100 , overallHivPrev_KZN_AC(1 , :) , overallHivPrev_KZN_AC(2 , :) , '*')
hold on 
% yPosError = abs(upper_prevVal - prevVal);
% yNegError = abs(lower_prevVal - prevVal);
% errorbar(prevValYrs , prevVal , yNegError , yPosError , 'ms')
xlabel('Year'); ylabel('Proportion of Population (%)'); title('HIV Prevalence (Ages 15-49)')
legend('Model' , 'National ANC Data' , 'Validation set: KZN Actual (Africa Center Data)')

%% HIV prevalence by age vs. AC data
hivAge = zeros(length(tVec) , age);
ageGroup = {'15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74'};

% 2010-2016 AC data
hivPrevF_val = [9.29	9.02	10.45	9.33	10.37	11.00	9.35
31.41	31.68	30.64	33.95	34.56	34.12	33.42
53.27	51.72	50.80	51.33	51.94	53.98	52.41
59.18	61.35	58.66	64.90	62.57	64.71	63.09
53.97	54.08	58.77	65.12	65.28	64.66	66.95
42.69	43.27	45.29	49.16	54.25	56.37	61.28
32.34	34.30	39.18	41.47	48.21	49.57	50.23
];

hivPrevM_val = [1.60	1.85	2.75	3.46	2.87	3.95	4.50
9.56	8.02	9.87	9.65	11.86	7.19	8.02
28.99	21.92	24.88	29.84	35.40	27.65	27.31
46.47	44.51	39.49	47.22	46.35	41.64	42.08
52.03	44.30	49.61	63.33	51.41	52.05	51.35
41.73	41.53	51.55	51.64	59.40	52.69	51.18
36.64	37.12	33.01	40.00	40.54	44.52	52.17
];

hivM = hivPrevM_dObs(: , 2) .* 100;
hivF = hivPrevF_dObs(: , 2) .* 100;
prevYears = unique(hivPrevF_dObs(: , 1));
prevYears2 = [2010 : 2016];
gen = {'Male' , 'Female'};
for g = 1 : gender
    hivPrevs = hivM;
    hivPrevs2 = hivPrevM_val;
    if g == 2
        hivPrevs = hivF;
        hivPrevs2 = hivPrevF_val;
    end
    figure()
    aVec = {16:20,21:25,26:30,31:35,36:40,41:45,46:50}; %,51:55,56:60,61:65,66:70,71:75};
    for a = 4 : 10
        %a = aVec{aInd};
        hivAgeInds = [toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , g , a , 1 : risk)); toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , g , a , 1 : risk))];
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , g , a , 1 : risk));
        hivAge(: , a-3) = sum(popVec(: , hivAgeInds) , 2);
        hivAgeRel = bsxfun(@rdivide , hivAge(: , a-3) , sum(popVec(: , ageInds) , 2)) * 100;
        subplot(3 , 3 , a-3) %subplot(4 , 3 , aInd)
        if a <= 11
            plot(tVec' , hivAgeRel , prevYears , hivPrevs((a-3) : 7 : end) ,'o'); % , ...
                %prevYears2 , hivPrevs2((aInd) : 7 : end) , 'bo');
        else
            plot(tVec' , hivAgeRel);
        end
        xlabel('Year'); ylabel('Prevalence (%)'); title([gen{g} , 's ages ' , ageGroup{a-3}]) % , ' HIV Prevalence'])
        xlim([1980 2019])
    end
    legend('Model' , 'Observed KZN: AHRI DHHS' , 'Observed KZN: AHRI DHHS')
end

%% HIV prevalence by gender vs. AC data
figure;
prevYears = unique(hivPrevF_dObs(: , 1));
hivRaw(:,:,1) = hivPrevM_dObs(: , 4:5);
hivRaw(:,:,2) = hivPrevF_dObs(: , 4:5);

hivData(: , : , 1) = zeros(length(prevYears) , 1);
hivData(: , : , 2) = zeros(length(prevYears) , 1);

for i = 1 : length(prevYears)
    for g = 1 : gender
        hivData(i,1,g) = (sumall(hivRaw(((i-1)*7+1):(i*7) , 1 , g)) ./ sumall(hivRaw(((i-1)*7+1):(i*7) , 2 , g))) .* 100;
    end
end

gen = {'Male' , 'Female'};
for g = 1 : gender
    hivInds = [toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , g , 4 : 10 , 1 : risk)); toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , g , 4 : 10 , 1 : risk))];
    totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , g , 4 : 10 , 1 : risk));
    hivPop = sum(popVec(: , hivInds) , 2);
    hivPopPrev = bsxfun(@rdivide , hivPop , sum(popVec(: , totInds) , 2)) * 100;
    subplot(1,2,g)
    plot(tVec' , hivPopPrev);
    hold all;
    plot(prevYears , hivData(:,:,g) , 'ro');
    xlabel('Year'); ylabel('Prevalence (%)'); title(gen{g});
    xlim([1980 2020])
    legend('Model' , 'Africa Center Data (Calibration)')
end

%% HIV prevalance, all ages by gender
figure()
for g = 1 : 2
    artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , 1 : age , 1 : risk));
    artPop = sum(popVec(: , artInds) , 2);
    hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , g , 1 : age , 1 : risk));
    allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , g , 1 : age , 1 : risk)); 
    hivPop = sum(popVec(: , hivInds) , 2);
    allPop = sum(popVec(: , allInds) , 2);
    plot(tVec , 100 * (hivPop + artPop) ./ allPop)
    hold on
end
xlabel('Year')
ylabel('Prevalence')
title('HIV Prevalence')

figure()
artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
artPop = sum(popVec(: , artInds) , 2);
hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk)); 
hivPop = sum(popVec(: , hivInds) , 2);
allPop = sum(popVec(: , allInds) , 2);
plot(tVec , 100 * (hivPop + artPop) ./ allPop)

xlabel('Year')
ylabel('Prevalence')
title('HIV Prevalence All')

%% HIV prevalence by age on x-axis
genderVec = {'Males (on and off ART)' , 'Females (on and off ART)'};
hivAge = zeros(16 , 2);
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

figure;
for g = 1 : gender
    %aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
    for a = 1 : 16
        %a = aVec{aInd};
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , a , 1 : risk));
        hivAge(a , g) = (sum(popVec(end,hivInds),2)/sum(popVec(end,ageInds),2))*100;
    end
    hold all;
    subplot(1,2,g);
    plot(1 : size(hivAge , 1) , hivAge(: , g) , '-');
    hold all;
    xlabel('Age Group'); ylabel('HIV Prevalence')
    set(gca , 'xtick' , 1 : length(hivAge) , 'xtickLabel' , ageGroup);
    title(genderVec{g})
    legend('Without ART dropout' , 'With ART dropout');
    grid on;
end

%% HIV mortality by age
% figure()
% bar(tVec , squeeze(hivDeaths(:,2,:)) , 'stacked')
% xlabel('Year'); ylabel('Deaths'); title('Female HIV-associated Deaths')
% legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
%     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
%     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthWest')
% figure()
% bar(tVec , bsxfun(@rdivide , squeeze(hivDeaths(:,2,:)) , sum(popVec , 2)) , 'stacked')
% xlabel('Year'); ylabel('Deaths relative to population'); title('Female HIV-associated Deaths')
% legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
%     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
%     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthWest')

%% Population size by gender
figure()
mInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , 1 , 1 : age , 1 : risk));
fInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
plot(tVec , sum(popVec(: , mInds) , 2))
hold on
plot(tVec , sum(popVec(: , fInds) , 2))
legend('Males' , 'Females')
xlabel('Year'); ylabel('Population')

%% Relative HIV prevalence for untreated and on ART
figure()
artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : endpoints , 1 : intervens , 1 : gender , 4 : 10 , 1 : risk));
artPop = sum(popVec(: , artInds) , 2);
hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
    1 : endpoints , 1 : intervens , 1 : gender , 4 : 10 , 1 : risk));
hivPop = sum(popVec(: , hivInds) , 2);
hiv_art = [100 * hivPop ./ sum(popTot , 2), 100 * artPop ./ sum(popTot , 2)];
area(tVec , hiv_art); %art ./ sum(popVec , 2) , tVec , hiv ./ sum(popVec , 2))
xlabel('Year')
ylabel('Proportion of Population ages 15-49 (%)')
title('Relative HIV Prevalence')
legend('Untreated', 'On ART' , 'Location' , 'NorthWest')
xlim([1980 2020]);

%% Proportion of HIV+ population on ART
figure()
artActual = [0	0	0	0	1	2	3	6, ...
    9	14	19	27	34	40	45	48];
yrsArtActual = [2000	2001	2002	2003	2004	2005 ...
    2006	2007	2008	2009	2010	2011	2012	2013 ...
    2014	2015];

plot(tVec , 100 * artPop ./ (hivPop + art) , yrsArtActual , artActual , '*')
xlabel('Year')
ylabel('Proportion of HIV Population')
title('Proportion on ART')
legend('Model' , 'Observed')

%% Proportion of HIV+ population on ART, by gender
figure()
artActualM = [2005	0.00483	0.004493	0.005193
2006	0.02004	0.018562	0.021656
2007	0.046485	0.042782	0.050622
2008	0.081635	0.074986	0.089097
2009	0.116299	0.107197	0.126437
2010	0.131834	0.122535	0.142661
2011	0.162606	0.15032	0.177079
2012	0.16853	0.154812	0.184915];

artActualF = [2005	0.009496	0.009117	0.009907
2006	0.044063	0.04216	0.046146
2007	0.085786	0.082091	0.08983
2008	0.143365	0.137734	0.149967
2009	0.200526	0.192269	0.210197
2010	0.222291	0.214523	0.229979
2011	0.271331	0.261925	0.281438
2012	0.306655	0.294358	0.319155];

artActualM(: , 2 : 4) = artActualM(: , 2 : 4) .* 100;
artActualF(: , 2 : 4) = artActualF(: , 2 : 4) .* 100;
artAct = {artActualM , artActualF};

for g = 1 : 2
    artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , 4 : 15 , 1 : risk));
    artPop = sum(popVec(: , artInds) , 2);
    hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , g , 4 : 15 , 1 : risk));
    hivPop = sum(popVec(: , hivInds) , 2);
    plot(tVec , 100 * artPop ./ (hivPop + artPop))
    hold on
    errorbar(artAct{g}(: , 1) , artAct{g}(: , 2) , ...
        artAct{g}(: , 2) - artAct{g}(: , 3) , ...
        artAct{g}(: , 4) - artAct{g}(: , 2))
    hold on
end
xlabel('Year')
ylabel('Proportion of HIV Population')
title('Proportion on ART')
legend('Model (Male)' , 'Observed(Male)' , 'Model (Female)' , 'Observed (Female)')

%% On ART by age
aVec = {1:5,6:1011:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80}; %{10:15,16:25,26:35,36:50,51:75};
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
figure;
hold all;
plot([1:age] , aMatrix(1,:) , '-')
hold all;
ylabel('Percent females on ART in 2020 by age');
ylim([0 110])
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
legend('Without ART dropout' , 'With ART dropout');
% legend('Without ART dropout' , 'With ART dropout: 6.19%' , 'With ART dropout: 11.8%' , 'With ART dropout: 11.8%, HIV mort on ART');

%% HPV prevalence over time by HIV status and gender
genders = {'Male' , 'Female'};
figure()
for g = 1 : gender
    % General
    hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 3 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 3 : age , 1 : risk))]);
    hpvPop = sum(popVec(: , hpvInds) , 2);
    popTot = popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 3 : age , 1 : risk)));
    % HIV+
    hpvHivInds = unique([toInd(allcomb(3 : 7 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 3 : age , 1 : risk)); toInd(allcomb(3 : 7 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 3 : age , 1 : risk))]);
    hpvHivPop = sum(popVec(: , hpvHivInds) , 2);
    popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 3 : age , 1 : risk)));
    %ART
    hpvArtInds = unique([toInd(allcomb(8 , 6 , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 3 : age , 1 : risk)); toInd(allcomb(8 , 6 , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 3 : age , 1 : risk))]);
    hpvArtPop = sum(popVec(: , hpvArtInds) , 2);
    popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , g , 3 : age , 1 : risk)));
    %HIV-
    hpvHivNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 3 : age , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 3 : age , 1 : risk))]);
    hpvHivNegPop = sum(popVec(: , hpvHivNegInds) , 2);
    popHivNegTot = popVec(: , toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 3 : age , 1 : risk)));

    subplot(2 , 1 , g)
    plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
    hold on
    plot(tVec , 100 * hpvHivNegPop ./ sum(popHivNegTot , 2))
    plot(tVec , 100 * hpvHivPop ./ sum(popHivTot , 2))
    plot(tVec , 100 * hpvArtPop ./ sum(popArtTot , 2))
    %axis([tVec(1) tVec(end) 0 100])
    % axis([-inf inf 0 min(max(hpvPop./sum(popTot , 2) * 100) + 10 , 100)])
    xlabel('Year'); ylabel('Prevalence (%)'); title([genders{g} , ' HPV Prevalence (ages 10+)'])
    legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Location' , 'NorthEastOutside')
end

%% HPV Prevalence by age in 2002 vs. McDonald 2014 data
ageGroup = {'17 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
hpv2002 = zeros(9 , 1);
hpvHIV2002 = hpv2002;
hpvNeg2002 = hpv2002;

aVec = {18:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
%for aInd = 1 : 13
for a = 4 : 12
    %a = aVec{aInd};
    hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpv2002(a - 3 , 1) = sum(popVec((2002 - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
    
    % HIV+
    hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpvHIV2002(a - 3 , 1) = sum(popVec((2002 - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
    
    % HIV-
    hpvInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpvNeg2002(a - 3 , 1) = sum(popVec((2002 - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
end

% McDonald 2014
hpvHivObs(: , 1) = [0.75
0.61
0.60
0.55
0.46
0.42
0.43
0.54
0.35];

hpvHivObs(: , 2) = [0.63
0.54
0.54
0.47
0.42
0.34
0.32
0.35
0.16];

hpvHivObs(: ,3) = [0.87
0.67
0.66
0.62
0.51
0.50
0.55
0.72
0.53];

hpvNegObs(: , 1) = [0.60
0.38
0.24
0.20
0.19
0.18
0.13
0.17
0.15];

hpvNegObs(: , 2) = [0.53
0.34
0.21
0.17
0.18
0.16
0.11
0.14
0.12];

hpvNegObs(: , 3) = [0.67
0.41
0.27
0.23
0.21
0.20
0.15
0.19
0.18];

hpvHivObs = hpvHivObs * 100;
hpvNegObs = hpvNegObs * 100;
figure()
% plot(1 : length(hpv2002) , hpv2002 , 'o-')
% hold all;
plot(1 : length(hpvHIV2002) , hpvHIV2002 , 'o-');
hold all;
plot(1 : length(hpvNeg2002) , hpvNeg2002 , 'o-')
hold all;
set(gca , 'xtickLabel' , ageGroup);

% general
% yPosError = abs(hrHpvObs(: , 3) - hrHpvObs(: , 1));
% yNegError = abs(hrHpvObs(: , 2) - hrHpvObs(: , 1));
% errorbar(1 : length(hrHpvObs) , hrHpvObs(: , 1) , yNegError , yPosError , 'rs')
% HIV+
yPosError = abs(hpvHivObs(: , 3) - hpvHivObs(: , 1));
yNegError = abs(hpvHivObs(: , 2) - hpvHivObs(: , 1));
errorbar(1 : length(hpvHivObs) , hpvHivObs(: , 1) , yNegError , yPosError , 'bs')
hold all;
%HIV-
yPosError = abs(hpvNegObs(: , 3) - hpvNegObs(: , 1));
yNegError = abs(hpvNegObs(: , 2) - hpvNegObs(: , 1));
errorbar(1 : length(hpvNegObs) , hpvNegObs(: , 1) , yNegError , yPosError , 'rs')

set(gca , 'xtick' , 1 : length(hpvNegObs) , 'xtickLabel' , ageGroup);
legend('HIV-Positive (year 2002)' , 'HIV-Negative (year 2002)' , 'Observed HIV-Positive: McDonald 2014' , 'Observed HIV-Negative: McDonald 2014')
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)')
%title('Age Specific hrHPV Prevalence in 2002')

%% HPV Prevalence by age and type in 2002 vs. McDonald 2014 data
ageGroup = {'17 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
hpv2002 = zeros(hpvTypeGroups , 9);
hpvHIV2002 = hpv2002;
hpvNeg2002 = hpv2002;

aVec = {18:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
%for aInd = 1 : 13
for a = 4 : 12
    %a = aVec{aInd};
    
    % HIV+
    hpvInds_vax = toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk));
    hpvInds_nonVax = toInd(allcomb(3 : 8 , 1 : viral , [1 , 7] , 2 : 5 , ...
        1 , 1 : intervens , 2 , a , 1 : risk));
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpvHIV2002(1 , a - 3) = sum(popVec((2002 - startYear) * stepsPerYear , hpvInds_vax))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
    hpvHIV2002(2 , a - 3) = sum(popVec((2002 - startYear) * stepsPerYear , hpvInds_nonVax))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
    
    % HIV-
    hpvInds_vax = toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk));
    hpvInds_nonVax = toInd(allcomb(1 : 2 , 1 : viral , [1 , 7] , 2 : 5 , ...
        1 , 1 : intervens , 2 , a , 1 : risk));
    ageInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpvNeg2002(1 , a - 3) = sum(popVec((2002 - startYear) * stepsPerYear , hpvInds_vax))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
    hpvNeg2002(2 , a - 3) = sum(popVec((2002 - startYear) * stepsPerYear , hpvInds_nonVax))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
end

% McDonald 2014
hpvHivObs(: , 1) = [0.75
0.61
0.60
0.55
0.46
0.42
0.43
0.54
0.35];

hpvHivObs(: , 2) = [0.63
0.54
0.54
0.47
0.42
0.34
0.32
0.35
0.16];

hpvHivObs(: ,3) = [0.87
0.67
0.66
0.62
0.51
0.50
0.55
0.72
0.53];

hpvNegObs(: , 1) = [0.60
0.38
0.24
0.20
0.19
0.18
0.13
0.17
0.15];

hpvNegObs(: , 2) = [0.53
0.34
0.21
0.17
0.18
0.16
0.11
0.14
0.12];

hpvNegObs(: , 3) = [0.67
0.41
0.27
0.23
0.21
0.20
0.15
0.19
0.18];

hpvHivObs = hpvHivObs * 100;
hpvNegObs = hpvNegObs * 100;
figure()
% plot(1 : length(hpv2002) , hpv2002 , 'co-')
% hold on
plot(1 : length(hpvHIV2002) , hpvHIV2002 , 'bo-');
hold all;
plot(1 : length(hpvNeg2002) , hpvNeg2002 , 'ro-')
hold all;
set(gca , 'xtickLabel' , ageGroup);

% general
% yPosError = abs(hrHpvObs(: , 3) - hrHpvObs(: , 1));
% yNegError = abs(hrHpvObs(: , 2) - hrHpvObs(: , 1));
% errorbar(1 : length(hrHpvObs) , hrHpvObs(: , 1) , yNegError , yPosError , 'rs')
% HIV+
yPosError = abs(hpvHivObs(: , 3) - hpvHivObs(: , 1));
yNegError = abs(hpvHivObs(: , 2) - hpvHivObs(: , 1));
errorbar(1 : length(hpvHivObs) , hpvHivObs(: , 1) , yNegError , yPosError , 'bs')
%HIV-
hold all;
yPosError = abs(hpvNegObs(: , 3) - hpvNegObs(: , 1));
yNegError = abs(hpvNegObs(: , 2) - hpvNegObs(: , 1));
errorbar(1 : length(hpvNegObs) , hpvNegObs(: , 1) , yNegError , yPosError , 'rs')

set(gca , 'xtick' , 1 : length(hpvNegObs) , 'xtickLabel' , ageGroup);
%legend('General' , 'HIV+' , 'HIV-' , 'McDonald 2014 - HIV+' , 'McDonald 2014 - HIV-')
legend('HIV-Positive, 9v (year 2002)' , 'HIV-Positive, non-9v (year 2002)' , 'HIV-Negative, 9v (year 2002)' , 'HIV-Negative, non-9v (year 2002)' , 'Observed HIV-Positive: McDonald 2014' , 'Observed HIV-Negative: McDonald 2014')
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)')
%title('Age Specific hrHPV Prevalence in 2002')

%% HPV prevalence by age and HIV status in 2008 vs. Mbulawa data, bar graph
% yearPrev = 2008;
% hpvHivPosInd = [toInd(allcomb(3 : 7 , 1 : viral , 2 : 5 , 1 , ...
%     1 , 1 : intervens , 2 , 19 : 67 , 1 : risk)); toInd(allcomb(8 , 6 , ...
%     2 : 5 , 1 , 1 , 1 : intervens , 2 , 19 : 67 , 1 : risk))];
% hpvHivNegInd = toInd(allcomb(1 : 2 , 1 , 2 : 5 , 1 , 1 , 1 : intervens , 2, ...
%     19 : 67 , 1 : risk));
% genHpvInd = toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , 1 , ...
%     1 , 1 : intervens , 2 , 19 : 67 , 1 : risk));
% 
% popHivInd = [toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%     1 : endpoints , 1 : intervens , 2 , 19 : 67 , 1 : risk)); toInd(allcomb(8 , 6 ,...
%     1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 2 , 19 : 67 , 1 : risk))];
% popNegInd = toInd(allcomb(1:2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%     1 : endpoints , 1 : intervens , 2 , 19 : 67 , 1 : risk));
% popGenInd = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%     1 : endpoints , 1 : intervens , 2 , 19 : 67 , 1 : risk));
% 
% popGen = sum(popVec((yearPrev - startYear) * stepsPerYear , popGenInd));
% popHiv = sum(popVec((yearPrev - startYear) * stepsPerYear , popHivInd));
% popNeg = sum(popVec((yearPrev - startYear) * stepsPerYear , popNegInd));
% 
% 
% hpvHivPos = sum(popVec((yearPrev - startYear) * stepsPerYear , hpvHivPosInd))...
%     ./ popHiv * 100;
% hpvHivNeg = sum(popVec((yearPrev - startYear) * stepsPerYear , hpvHivNegInd))...
%     ./ popNeg * 100;
% genHpv = sum(popVec((yearPrev - startYear) * stepsPerYear , genHpvInd))...
%     ./ popGen * 100;
% 
% %Mbulawa KZN
% kzn_HpvData = [0.39 , 0.35 , 0.44; % gen hpv
%     0.52, 0.46 , 0.58; % hpvHivPos
%     0.23 , 0.17 , 0.28] .* 100; %hpvHivNeg
% 
% upper = abs(kzn_HpvData(: , 2) - kzn_HpvData(: , 1));
% lower = abs(kzn_HpvData(: , 3) - kzn_HpvData(: , 1));
% hpvGroup = {'General Population' , 'HIV+' , 'HIV-'};
% figure()
% bar([genHpv , hpvHivPos , hpvHivNeg])
% hold on
% errorbar(1 : 3, kzn_HpvData(: , 1) , lower , upper , 'rs') 
% set(gca , 'xtick' , 1 : 3 , 'xtickLabel' , hpvGroup);
% legend('Model' , 'Mbulawa')
% title(['HPV Prevalence by HIV Status in ' , num2str(yearPrev)])
% ylabel('Prevalence (%)')
% 
%% HPV prevalence by HIV status over time vs. Mbulawa data
% figure()
% popGen = sum(popVec(: , popGenInd) , 2);
% popHiv = sum(popVec(: , popHivInd) , 2);
% popNeg = sum(popVec(: , popNegInd) , 2);
% hpvHivPos = sum(popVec(: , hpvHivPosInd) , 2) ./ popHiv * 100;
% hpvHivNeg = sum(popVec(:, hpvHivNegInd) , 2) ./ popNeg * 100;
% genHpv = sum(popVec(: , genHpvInd) , 2) ./ popGen * 100; 
% plot(tVec , genHpv)
% hold on
% plot(tVec , hpvHivPos)
% hold on
% plot(tVec , hpvHivNeg)
% hold on
% % plot([tVec(1) , tVec(end)] , [kzn_HpvData(1 , 1) , kzn_HpvData(1 , 1)] , 'r-'); 
% plot(tVec , ones(length(tVec) , 1) .* kzn_HpvData(1 , 1) , 'r-')
% hold on
% plot(tVec , ones(length(tVec) , 1) .* kzn_HpvData(2 , 1) , 'b-')
% hold on
% plot(tVec , ones(length(tVec) , 1) .* kzn_HpvData(3 , 1) , 'c-')
% xlabel('Year'); ylabel('HPV Prevalence'); title('HPV Prevalence by HIV Status')
% legend('General Population' , 'HIV+' , 'HIV-' , 'Mbulawa (General)' , ...
%     'Mbulawa (HIV+)' , 'Mbulawa (HIV-)' , 'location' , 'NorthEastOutside')
% 
%% CIN2/3 prevalence by HIV status, HPV type, and age in 2002 vs. McDonald 2014 data and Firnhaber data
cinPos2002 = zeros(hpvTypeGroups - 1 , 10);
cinNeg2002 = cinPos2002;
ageGroup = {'17 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
%aVec = {18:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
%for aInd = 1 : 13
for a = 4 : 13 % note, age group 4 is 17-19 in the data
    %a = aVec{aInd};
    % HIV-positive (on and not on ART)
    cinInds_vax = toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk));
    cinInds_nonVax = toInd(allcomb(3 : 8 , 1 : viral , [1 : 3 , 7] , 4 : 5 , ...
        1 , 1 : intervens , 2 , a , 1 : risk));
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cinPos2002(1 , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , cinInds_vax)))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
    cinPos2002(2 , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , cinInds_nonVax)))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
    % HIV-negative
    cinNegInds_vax = toInd(allcomb(1 : 2 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk));
    cinNegInds_nonVax = toInd(allcomb(1 : 2 , 1 : viral , [1 : 3 , 7] , 4 : 5 , ...
        1 , 1 : intervens , 2 , a , 1 : risk));
    ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cinNeg2002(1 , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , cinNegInds_vax)))...
        ./ (sum(popVec((2002 - startYear) * stepsPerYear , ageNegInds))) * 100;
    cinNeg2002(2 , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , cinNegInds_nonVax)))...
        ./ (sum(popVec((2002 - startYear) * stepsPerYear , ageNegInds))) * 100;
end

% McDonald 2014
cinPosAct(: , 1) = [0.125
0.054
0.128
0.154
0.081
0.054
0.079
0.071
0.077
0.077]; % mean

cinPosAct(: , 2) = [0.03
0.02
0.09
0.10
0.05
0.02
0.02
0.00
0.00
0.00
]; % lb

cinPosAct(: , 3) = [0.22
0.08
0.17
0.21
0.11
0.09
0.14
0.17
0.22
0.22]; % ub

%figure()
hold all;
subplot(2 , 1 , 1)
cinPosAct = cinPosAct .* 100; % convert to %
yPosError = abs(cinPosAct(: , 3) - cinPosAct(: , 1));
yNegError = abs(cinPosAct(: , 2) - cinPosAct(: , 1));
plot(1 : length(cinPos2002) , cinPos2002 ,'ko-')
hold on
errorbar(1 : length(cinPosAct) , cinPosAct(: , 1) , yNegError , yPosError , 'rs')
%legend('9v' , 'non-9v' , 'McDonald 2014')
legend('HIV-Positive (year 2002)' , 'Observed HIV-Positive: McDonald 2014' , 'HIV-Positive, 9v (year 2002)' , 'HIV-Positive, non-9v (year 2002)')
set(gca , 'xtick' , 1 : length(cinPosAct) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN2/3 Prevalence (%)')
title('HIV-Positive')
ylim([0 25])

cinNegAct(: , 1) = [0.016
0.027
0.021
0.036
0.029
0.031
0.031
0.021
0.014
0.014]; % mean

cinNegAct(: , 2) = [0.00
0.02
0.01
0.02
0.02
0.02
0.02
0.01
0.00
0.00]; % lb

cinNegAct(: , 3) = [0.03
0.04
0.03
0.05
0.04
0.04
0.04
0.03
0.03
0.03]; % ub

subplot(2 , 1 , 2)
cinNegAct = cinNegAct .* 100; % convert to %
plot(1 : length(cinNeg2002) , cinNeg2002 , 'ko-')
hold on
yPosError = abs(cinNegAct(: , 3) - cinNegAct(: , 1));
yNegError = abs(cinNegAct(: , 2) - cinNegAct(: , 1));
errorbar(1 : length(cinNegAct) , cinNegAct(: , 1) , yNegError , yPosError , 'rs')
%legend('9v' , 'non-9v' , 'McDonald 2014')
legend('HIV-Negative (year 2002)' , 'Observed HIV-Negative: McDonald 2014' , 'HIV-Negative, 9v (year 2002)' , 'HIV-Negative, non-9v (year 2002)')
set(gca , 'xtick' , 1 : length(cinNegAct) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN2/3 Prevalence (%)')
title('HIV-Negative')
ylim([0 25])

% CIN Prevalence among HIV+ Women
% cinHiv = [0.46	0.18 0.09] .* 100; % Observed, Johannesburg, Cynthia Firnhaber
% cinHiv2017 = zeros(3 , 1);
% aVec = {16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
% for aInd = 1 : 13
%     a = aVec{aInd};
%     for cin = 3 : 5
%         % HIV+
%         cinInds = toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvVaxStates , cin - 1, ...
%             1 : endpoints , 2 , a , 1 : risk));
%         ageInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
%             2 , a , 1 : risk));
%         cinArtInds = toInd(allcomb(10 , 6 , 2 : hpvVaxStates , cin - 1, ...
%             1 : endpoints , 2 , a , 1 : risk));
%         ageArtInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , cin - 1 , 1 : endpoints , ...
%             2 , a , 1 : risk));
%         cinHiv2017(cin - 2) = (sum(popVec((2017 - startYear) * stepsPerYear , cinInds)) + sum(popVec((2017 - startYear) * stepsPerYear , cinArtInds)))...
%             ./ (sum(popVec((2017 - startYear) * stepsPerYear , ageInds)) + sum(popVec((2017 - startYear) * stepsPerYear , ageArtInds)))* 100;
%     end
% end
% 
% figure()
% % cinGroup = {'CIN 1' , 'CIN 2' , 'CIN 3'};
% plot(1 : length(cinHiv2017) , cinHiv2017 , 'o-' , 1 : length(cinHiv) , cinHiv , 'rs');
% ylabel('Prevalence (%)')
% % set(gca , 'xtickLabel' , cinGroup);
% xlabel('CIN Stage')
% legend('Model' , 'Firnhaber')
% title('CIN Prevalence in HIV Positive Women in 2017')
 
%% CIN2/3 prevalence for All HR HPV types combined by HIV status and age in 2002 vs. McDonald 2014 data
cinPos2002 = zeros(10 , 1);
cinNeg2002 = cinPos2002;
ageGroup = {'17-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};
%aVec = {18:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
for a = 4 : 13  %note, age group 4 is 17-19 in the data
    %a = aVec{aInd};
    % HIV-positive (on and not on ART)
    cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cinPos2002(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , cinInds)))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
    % HIV-negative
    cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cinNeg2002(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , cinNegInds)))...
        ./ (sum(popVec((2002 - startYear) * stepsPerYear , ageNegInds))) * 100;
end

% McDonald 2014
cinPosAct(: , 1) = [0.125
0.054
0.128
0.154
0.081
0.054
0.079
0.071
0.077
0.077]; % mean

cinPosAct(: , 2) = [0.03
0.02
0.09
0.10
0.05
0.02
0.02
0.00
0.00
0.00
]; % lb

cinPosAct(: , 3) = [0.22
0.08
0.17
0.21
0.11
0.09
0.14
0.17
0.22
0.22]; % ub

figure();
subplot(2 , 1 , 1);
cinPosAct = cinPosAct .* 100; % convert to %
yPosError = abs(cinPosAct(: , 3) - cinPosAct(: , 1));
yNegError = abs(cinPosAct(: , 2) - cinPosAct(: , 1));
plot(1 : length(cinPos2002) , cinPos2002 ,'o-')
hold on
errorbar(1 : length(cinPosAct) , cinPosAct(: , 1) , yNegError , yPosError , 'rs')
legend('HR HPV CIN 2/3' , 'McDonald 2014')
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('Age Specific CIN 2/3 Prevalence Among HIV+ in 2002')
ylim([0 25])

cinNegAct(: , 1) = [0.016
0.027
0.021
0.036
0.029
0.031
0.031
0.021
0.014
0.014]; % mean

cinNegAct(: , 2) = [0.00
0.02
0.01
0.02
0.02
0.02
0.02
0.01
0.00
0.00]; % lb

cinNegAct(: , 3) = [0.03
0.04
0.03
0.05
0.04
0.04
0.04
0.03
0.03
0.03]; % ub

subplot(2 , 1 , 2)
cinNegAct = cinNegAct .* 100; % convert to %
plot(1 : length(cinNeg2002) , cinNeg2002 , 'o-')
hold on
yPosError = abs(cinNegAct(: , 3) - cinNegAct(: , 1));
yNegError = abs(cinNegAct(: , 2) - cinNegAct(: , 1));
errorbar(1 : length(cinNegAct) , cinNegAct(: , 1) , yNegError , yPosError , 'rs')
legend('HR HPV CIN 2/3' , 'McDonald 2014')
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('Age Specific CIN 2/3 Prevalence Among HIV- in 2002')
ylim([0 25])

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
plot(tVec((2011 - startYear) * stepsPerYear) , 58.2 , 'ko')
hold all;
plot(tVec , hpv_nonVax);
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 41.8 , 'o');
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('HPV');
ylim([0 100]);
xlim([2000 2015]);
%legend('9v-type HPV' , 'Guan 2012: 9v' , 'Non-9v-type HPV' , 'Guan 2012: non-9v');

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
plot(tVec((2011 - startYear) * stepsPerYear) , 67.1 , 'ko')
hold all;
plot(tVec , cin_nonVax);
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 32.9 , 'o')
ylim([0 100]);
xlim([2000 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('CIN1')
%legend('9v-type HPV' , 'Guan 2012: 9v' , 'Non-9v-type HPV' , 'Guan 2012: non-9v');

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
plot(tVec((2011 - startYear) * stepsPerYear) , 75.0 , 'ko')
hold all;
plot(tVec , cin_nonVax);
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 25.0 , 'o')
ylim([0 100]);
xlim([2000 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('CIN2')
%legend('9v-type HPV' , 'Guan 2012: 9v' , 'Non-9v-type HPV' , 'Guan 2012: non-9v');

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
plot(tVec((2011 - startYear) * stepsPerYear) , 82.9 , 'ko')
hold all;
plot(tVec , cin_nonVax);
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 17.1 , 'o')
ylim([0 100]);
xlim([2000 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('CIN3')
%legend('9v-type HPV' , 'Guan 2012: 9v' , 'Non-9v-type HPV' , 'Guan 2012: non-9v');

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
plot(tVec((2011 - startYear) * stepsPerYear) , 94.8 , 'ko')
hold all;
plot(tVec , cc_nonVax);
hold all;
plot(tVec((2011 - startYear) * stepsPerYear) , 5.2 , 'o')
ylim([0 100]);
xlim([2000 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('Cervical Cancer')
legend('9v-type HPV' , 'Guan 2012: 9v' , 'Non-9v-type HPV' , 'Guan 2012: non-9v');

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
% % ccIncYears = [2017 , 2003 , 1994 , 2012];
% % ccAgeRel = zeros(16 , length(ccIncYears));
% % ccAgeNegRel = ccAgeRel;
% % ccAgePosRel = zeros(16 , 4 , length(ccIncYears));
% % ccNegPosArt = zeros(16 , 3 , length(ccIncYears));
% % ccNegPosArtTot = ccNegPosArt;
% % ccArtRel = ccAgeRel;
% % fScale = 10^5;
% % ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
% %     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
% %     '60-64' , '65-69' , '70-74' , '75-79'};
% % annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 
% % ccYrs = ((ccIncYears - startYear) * stepsPerYear :...
% %     (ccIncYears + 1 - startYear) * stepsPerYear);
% % 
% % aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
% % for aInd = 1 : 16
% %     a = aVec{aInd};
% %     for y = 1 : length(ccIncYears)
% %         % Year
% %         yr_start = (ccIncYears(y) - 1 - startYear)  .* stepsPerYear;
% %         yr_end = (ccIncYears(y) - startYear) .* stepsPerYear - 1;
% %         % Total population
% %         ageInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %             2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 8 : 10 , 1 : endpoints , ...
% %             2 , a , 1 : risk))];
% %         ccAgeRel(aInd , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end , ...
% %             1 : disease , 1 : hpvVaxStates , a) , 2) , 3) , 4)) ...
% %             ./ (annlz(sum(popVec(yr_start : yr_end , ageInds) , 2)) ...
% %             ./ stepsPerYear) * fScale;
% % 
% %         % HIV Negative
% %         ageNegInds = [toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %             2 , a , 1 : risk)); toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 8 : 10 , 1 : endpoints , ...
% %             2 , a , 1 : risk))];
% %         ccAgeNegRel(aInd , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end...
% %             , 1 , 1 : hpvVaxStates , a) , 2) , 3) , 4)) ...
% %             ./ (annlz(sum(popVec(yr_start : yr_end , ageNegInds) , 2)) ...
% %             ./ stepsPerYear) * fScale;
% % 
% % 
% %         % Acute and CD4 > 500
% %         agePosInds = [toInd(allcomb(2 : 3 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %             2 , a , 1 : risk)); toInd(allcomb(2 : 3 , 1 : viral , 1 : hpvVaxStates , 8 : 10 , 1 : endpoints , ...
% %             2 , a , 1 : risk))];
% %         ccAgePosRel(aInd , 1 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end ...
% %             , 2 : 3 , 1 : hpvVaxStates , a), 2) , 3) , 4)) ...
% %             ./ (annlz(sum(popVec(yr_start : yr_end , agePosInds) , 2)) ...
% %             ./ stepsPerYear) * fScale;
% % 
% %         for d = 4 : 6
% %         % HIV Positive CD4 500-350 -> CD4 < 200
% %         agePosInds = [toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %             2 , a , 1 : risk)); toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 8 : 10 , 1 : endpoints , ...
% %             2 , a , 1 : risk))];
% %         ccAgePosRel(aInd , d - 2 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end...
% %             , d , 1 : hpvVaxStates , a), 2) , 3) , 4)) ...
% %             ./ (annlz(sum(popVec(yr_start : yr_end , agePosInds) , 2)) ...
% %             ./ stepsPerYear) * fScale;
% %         end
% % 
% %         % On ART
% %         ageArtInds = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %             2 , a , 1 : risk)); toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 8 : 10 , 1 : endpoints , ...
% %             2 , a , 1 : risk))];
% %         ccArtRel(aInd , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end ...
% %             , 10 , 1 : hpvVaxStates , a) , 2) , 3) , 4)) ...
% %             ./ (annlz(sum(popVec(yr_start : yr_end , ageArtInds) , 2)) ...
% %             ./ stepsPerYear) * fScale;
% % 
% %         % Proportion of cervical cancers by HIV/ART status and age
% %         % Total by age
% %         ageTotal = annlz(sum(sum(sum(popVec(yr_start : yr_end , ...
% %             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %             2 , a , 1 : risk))), 2 ) , 3), 4) + sum(sum(sum(popVec(yr_start : yr_end , ...
% %             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 8 : 10 , 1 : endpoints , ...
% %             2 , a , 1 : risk))), 2) , 3) , 4)) ./ stepsPerYear;
% % 
% %         % HIV-
% %         ccNegPosArt(aInd , 1 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end...
% %             , 1 , 1 : hpvVaxStates , a), 2) , 3) , 4)) ...
% %             ./ ageTotal;
% %         % HIV+
% %         ccNegPosArt(aInd , 2 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end ...
% %             , 2 : 6 , 1 : hpvVaxStates , a) , 2) , 3) , 4)) ...
% %             ./ ageTotal;
% %         % ART
% %         ccNegPosArt(aInd , 3 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end...
% %             , 10 , 1 : hpvVaxStates , a), 2) , 3) , 4)) ...
% %             ./ ageTotal;
% % 
% % 
% %          ageAllPosInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : endpoints , ...
% %             2 , a , 1 : risk)); toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 8 : 10 , 1 : endpoints , ...
% %             2 , a , 1 : risk))];    
% %         % HIV-
% %         ccNegPosArtTot(aInd , 1 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end...
% %             , 1 , 1 : hpvVaxStates , a), 2) , 3) , 4)) ...
% %             ./ (annlz(sum(popVec(yr_start : yr_end , ageNegInds) , 2)) ...
% %             ./ stepsPerYear) * fScale;
% %         % HIV+
% %         ccNegPosArtTot(aInd , 2 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end ...
% %             , 2 : 6 , 1 : hpvVaxStates , a) , 2) , 3) , 4)) ...
% %             ./ (annlz(sum(popVec(yr_start : yr_end , ageAllPosInds) , 2)) ...
% %             ./ stepsPerYear) * fScale;
% %         % ART
% %         ccNegPosArtTot(aInd , 3 , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end ...
% %             , 10 , 1 : hpvVaxStates , a), 2) , 3) , 4)) ...
% %             ./ (annlz(sum(popVec(yr_start : yr_end , ageArtInds) , 2)) ...
% %             ./ stepsPerYear) * fScale;
% %     end
% % end
% % 
% % globocan =[0
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
% % globocan_ub = [0.00
% % 2.65
% % 8.85
% % 45.19
% % 53.41
% % 67.05
% % 80.83
% % 78.97
% % 128.87
% % 105.27
% % 118.70
% % 111.81
% % 95.20];
% % 
% % globocan_lb = [0.00
% % 0.00
% % 0.41
% % 9.97
% % 12.61
% % 25.00
% % 45.69
% % 36.31
% % 50.55
% % 57.08
% % 62.69
% % 42.43
% % 52.01];
% % 
% % ccNegPosArt = ccNegPosArt .* fScale;
% % 
% % for y = 1 : length(ccIncYears)
% %     ccIncYear = ccIncYears(y);
% %     figure()
% %     plot(1 : size(ccAgeRel , 1) , ccAgeRel(: , y) , '-o' , 1 : length(ccAgeNegRel(: , y)) , ...
% %         ccAgeNegRel(: , y) , '-o' , 1 : size(ccAgePosRel , 1) , ccAgePosRel(: , : , y) , '-o' , ...
% %         1 : size(ccArtRel , 1) , ccArtRel(: , y) , '-o');
% %     xlabel('Age Group'); ylabel('Incidence per 100,000')
% %     set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
% %     title(['Cervical Cancer Incidence in ' num2str(ccIncYear)])
% % %     legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Location' , 'NorthWest');
% % %     legend('General' , 'HIV-' , ' Acute and CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , ...
% % %         'CD4 < 200' , 'ART' , 'Location' , 'NorthWest')
% %     hold on
% %     % globocan data
% %     plot(4 : 16 , globocan , '-' , 4 : 16 , globocan_ub , 'r--' , 4 : 16 , globocan_lb , 'r--')
% % %     legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Globocan' , 'Upper Bound' , 'Lower Bound' , ...
% % %         'Location' , 'NorthEastOutside')
% %     legend('General' , 'HIV-' , ' Acute and CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , ...
% %         'CD4 < 200' , 'ART' , 'Globocan' , 'Upper Bound' , 'Lower Bound' , ...
% %         'Location' , 'NorthEastOutside')
% % 
% % %     figure()
% % %     bar(1 : length(ccNegPosArt(: , : , y)) , ccNegPosArt(: , : , y), 'stacked')
% % %     xlabel('Age Group'); ylabel('Incidence per 100,000')
% % %     set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
% % %     title(['Cervical Cancer Incidence Distribution in ' , num2str(ccIncYear)])
% % %     legend('HIV-' , 'HIV+' , 'ART')
% %     
% %     figure()
% %     plot(1 : size(ccAgeRel , 1) , ccAgeRel(: , y) , '-ko' , 1 : size(ccNegPosArtTot, 1) , ccNegPosArtTot(: , 1 , y) , '-kp' , 1 : size(ccNegPosArtTot, 1) , ...
% %         ccNegPosArtTot(: , 2 , y) , '-k+' , 1 : size(ccNegPosArtTot, 1) , ccNegPosArtTot(: , 3 , y) , '-k^');
% %     hold on
% %     plot(4 : 16 , globocan , 'r-' , 4 : 16 , globocan_ub , 'r--' , 4 : 16 , globocan_lb , 'r--')
% %     title(['Cervical Cancer Incidence Distribution in ' , num2str(ccIncYear)])
% %     legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Globocan SA' , 'Upper Bound' , 'Lower Bound')
% %     xlabel('Age Group'); ylabel('Incidence per 100,000')
% %     set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
% %     title(['Cervical Cancer Incidence in ' num2str(ccIncYear)])
% % end
% 
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
% % popByAge = zeros(length(tVec) , age);
% % for i = 1 : length(tVec)
% %     for a = 1 : age
% %         ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
% %             1 : endpoints , 1 : gender , a , 1 : risk));
% %         popByAge(i , a) = sum(popVec(i , ageInds) , 2);
% %     end
% % end
% % %%
% % figure()
% % area(tVec , bsxfun(@rdivide , popByAge , sum(popVec , 2)))
% % xlabel('Year'); ylabel('Relative Population Size'); title('Population')
% % legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
% %     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
% %     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')
% % 
% % %%
% % figure()
% % ages = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
% %     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
% %     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
% % surf(1 : age , tVec , popByAge);
% % set(gca , 'xtickLabel' , ages);
% % set(gca , 'xLim' , [1 age]);
% % set(gca , 'yLim' , [tVec(1) tVec(end)]);
% % xlabel('Ages'); ylabel('Year'); zlabel('Population Size'); title('Population by Age')
% 
%% HIV incidence by age and gender
% % figure()
% % for a = 1 : age
% %     subplot(9,9,a)
% %     for r = 1 : risk
% % %     for g = 1 : 2
% %         hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , 2 , a , r)); ...
% %             toInd(allcomb(7 : 9 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %             1 : endpoints , 2 , a , r))];
% %         hivSus = annlz(sum(popVec(1:end-1 , hivSusInds) , 2)) ./ stepsPerYear;
% %         plot(tVec(1 : stepsPerYear : end-1) , ...
% %             annlz(sum(sum(newHiv(1:end-1 , 2 , a , r) ...
% %             , 3) , 4)) ./ hivSus * 100)
% %         %axis([startYear , endYear , 0 , 100])
% %         hold all
% % %     end
% %     end
% %     xlabel('Year'); ylabel('Rate Per 100'); title('HIV Incidence: Females')
% % end
% % %legend('Male' , 'Female')
% % legend('lr' , 'mr' , 'hr')
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
     1 : hpvVaxStates , 6 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     1 : hpvVaxStates , 6 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     1 : hpvVaxStates , 6 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
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

%% CIN1 prevalence by HIV group and HPV type
% Vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , 3 , 1 : hpvNonVaxStates , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , 3 , 1 : hpvNonVaxStates , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , 3 , 1 : hpvNonVaxStates , ...
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
     1 : hpvVaxStates , 3 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     1 : hpvVaxStates , 3 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     1 : hpvVaxStates , 3 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
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

%% CIN2 prevalence by HIV group and HPV type
% Vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , 4 , 1 : hpvNonVaxStates , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , 4 , 1 : hpvNonVaxStates , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , 4 , 1 : hpvNonVaxStates , ...
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
     1 : hpvVaxStates , 4 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     1 : hpvVaxStates , 4 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     1 : hpvVaxStates , 4 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
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

%% HPV prevalence by HIV group and HPV type
% Vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , 2 , 1 : hpvNonVaxStates , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , 2 , 1 : hpvNonVaxStates , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , 2 , 1 : hpvNonVaxStates , ...
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
     1 : hpvVaxStates , 2 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(popVec(: , ccHivInds) , 2);
popHivTot = popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     1 : hpvVaxStates , 2 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(popVec(: , ccArtInds) , 2);
popArtTot = popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     1 : hpvVaxStates , 2 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
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
