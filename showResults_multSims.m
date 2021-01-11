function[] = showResults_multSims()

%% Load parameters and results
paramDir = [pwd , '\Params\'];

[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd , ...
    annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4 , ...
    mue , mue2 , mue3 , mue4 , epsA_vec , epsR_vec , ...
    yr , ...
    hivOn , betaHIV_mod , muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , c2c3Mults , c1c2Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear , circNatStartYear , vaxStartYear , ...
    baseline , cisnet , who , whob , circProtect , condProtect , MTCTRate , ...
    hyst , OMEGA , ...
    ccInc2012_dObs , ccInc2018_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    cinPos2015_dObs , cinNeg2015_dObs , hpv_hiv_dObs , hpv_hivNeg_dObs , ...
    hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
    popAgeDist_dObs , totPopSize_dObs , ...
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
    infhpvNonVaxInds , fromVaxNoScrnInds , fromVaxScrnInds , toNonVaxNoScrnInds , ...
    toNonVaxScrnInds , ageInd , riskInd , deathInds , ...
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3 , dFertMat3 , deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2(1 , 0 , [] , [] , []);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1)

% Indices of calib runs to plot
fileInds = {'6_1' , '6_2' , '6_3' , '6_6' , '6_8' , '6_9' , '6_11' , ...
    '6_12' , '6_13' , '6_15' , '6_20' , '6_21' , '6_22' , '6_26' , ...
    '6_27' , '6_32' , '6_34' , '6_35' , '6_38' , '6_39' , '6_40' , ...
    '6_41' , '6_42' , '6_45' , '6_47'};    % 22Apr20Ph2V11
% '6_26' , '6_27' , ...
%     '6_28' , '6_29' , '6_30' , '6_31' , '6_32' , '6_33' , '6_34' , ...
%     '6_35' , '6_36' , '6_37' , '6_38' , '6_39' , '6_40' , '6_41' , ...
%     '6_42' , '6_43' , '6_44' , '6_45' , '6_46' , '6_47' , '6_48' , ...
%     '6_49' , '6_50'};    % 22Apr20Ph2V11
% '6_1' , '6_2' , '6_3' , '6_4' , '6_5' , '6_6' , ...
%     '6_7' , '6_8' , '6_9' , '6_10' , '6_11' , '6_12' , '6_13' , ...
%     '6_14' , '6_15' , '6_16' , '6_17' , '6_18' , '6_19' , '6_20' , ...
%     '6_21' , '6_22' , '6_23' , '6_24' , '6_25'};    % 22Apr20Ph2V11
% fileInds = {'4_1966' , '6_669' , '2_1058' , '5_853' , '2_903' , '0_5236' , '1_246' , ...
%     '2_1437' , '5_1416' , '2_1269' , '0_1947' , '3_2323' , '2_1959' , '4_297' , ...
%     '0_1430' , '5_2108' , '6_2204' , '6_2699' , '4_3256' , '3_147' , '6_3332' , '5_2419' , '3_3141' , '5_932' , ...
%     '0_3742'}; % , '5_1344' , '4_702' , '0_3204' , '5_799'};    % 22Apr20Ph2V11
% fileInds = {'1_1439' , '0_1611' , '0_4873' , '0_8201' , '1_5425' , ...
%     '0_8745' , '0_7117' , '1_1474' , '0_5506' , '0_408' , '1_5948' , ...
%     '2_1116' , '2_1631' , '2_6083' , '2_1512' , '2_4652' , '0_7729' , ...
%     '2_2947' , '2_2752' , '0_5935' , '0_2471' , '2_3927' , ...
%     '0_6494' , '1_2102' , '2_281' , '0_8976'};    % 22Apr20Ph2V7
% fileInds = {'0_3292' , '0_7871' , '0_4325' , '1_3259' , '0_1474' , ...
%     '1_4214' , '0_2709' , '1_1367' , ...
%     '2_3899' , '2_4116' , '2_413' , '2_4577'};    % 22Apr20Ph2V3
% fileInds = {'11_2946' , '11_4738' , '7_1476' , '2_1779' , '6_1657' , '11_2200' , ...
%     '10_3755' , '5_1279' , '0_6202' , '6_452' , '4_4034' , '5_2720' , '10_1939' , ...
%     '0_2605' , '11_2511' , '9_3353' , '4_3559' , '6_2571' , '10_2890' , '7_180' , ...
%     '5_5825' , '1_4859' , '4_4147' , '2_2314' , '7_741'}; % , ... %  % 22Apr20Ph2V2
%     '10_3270' , '10_3664' , '8_2827' , '11_2175' , '7_1426' , '2_1833' , '3_2610' , ...
%     '7_240' , '6_3972' , '10_4629' , '7_5509' , '6_1931' , '8_4473' , '0_2709' , ...
%     '7_3753' , '0_845' , '2_3468' , '8_2502' , '9_3299' , '11_669' , '9_4007' , ...
%     '8_5698' , '7_2004' , '5_33' , '8_1372'};
% fileInds = {'12_3346' , '12_2618' , '11_932' , '16_3038' , '8_597' , '12_2550' , ...  % 22Apr20Ph2
%     '12_3895' , '22_487' , '8_2705' , '22_3250' , '15_2550' , '14_563' , ...
%     '4_1887' , '10_688' , '18_3391' , '14_2659' , '19_2814' , '18_903' , ...
%     '22_2697' , '4_1676' , '4_2471' , '15_2517' , '16_1709' , '12_2481' , '16_3992'};
% fileInds = {'2_846' , '14_947' , '24_2579' , '25_3388' , '16_3127' , '12_689' , ...  % 22Apr20Ph1, top 50 sets
%     '10_2727' , '17_3986' , '16_2194' , '15_3850' , '9_334' , '0_6657' , ...
%     '16_2364' , '25_921' , '4_711' , '19_1017' , '24_1595' , '4_2361' , ...
%     '15_2155' , '21_1211' , '25_1470' , '17_594' , '19_1779' , '11_1541' , ...
%     '21_751' , '12_3055' , '6_746' , '20_944' , '13_3012' , '18_387' , ...
%     '17_1649' , '23_1150' , '25_982' , '17_3242' , '14_2649' , '21_2744' , ...
%     '16_2701' , '23_3179' , '20_1864' , '19_3788' , '24_2009' , '10_425' , ...
%     '11_3176' , '25_438' , '18_45' , '15_532' , '20_3201' , ...
%     '5_1822' , '7_2249' , '15_2965'};
nRuns = length(fileInds);

resultsDir = [pwd , '\HHCoM_Results\'];
for j = 1 : nRuns
    % Load results
    pathModifier = ['toNow_22Apr20Ph2V11_noBaseVax_baseScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_' , fileInds{j}];
    load([resultsDir , pathModifier])
   

    %% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

    %% Population size over time vs. UN/SSA data
    % Ages 0-79
    % All HIV-negatives
    hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
    % HIV-positives not on ART
    hivNoART = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
    % HIV-positives on ART
    art = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
    genArray = {hivNeg , hivNoART , art};
    totalPop0_79 = sum(popVec(:,genArray{1}),2) + sum(popVec(:,genArray{2}),2) + sum(popVec(:,genArray{3}),2);

    % Ages 0-69 (future projections only up to age 69)
    % All HIV-negatives
    hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 1 : 14 , 1 : risk));
    % HIV-positives not on ART
    hivNoART = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 1 : 14 , 1 : risk));
    % HIV-positives on ART
    art = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 1 : 14 , 1 : risk));
    genArray = {hivNeg , hivNoART , art};
    totalPop0_69 = sum(popVec(:,genArray{1}),2) + sum(popVec(:,genArray{2}),2) + sum(popVec(:,genArray{3}),2);

    % Load calibration data from Excel (years, values)
    file = [pwd , '/Config/Population_validation_targets.xlsx'];
    historicalPop0_69 = zeros(5,2);
    futurePop0_69 = zeros(16,2);
    historicalPop0_69(:,1) = xlsread(file , 'Demographics' , 'B91:F91'); % years
    historicalPop0_69(:,2) = xlsread(file , 'Demographics' , 'B130:F130') .* 1000; % estimates
    futurePop0_69(:,1) = xlsread(file , 'Demographics' , 'C144:R144'); % years
    futurePop0_69(:,2) = xlsread(file , 'Demographics' , 'C146:R146') .* 1000; % projections

    % Calibration error bars
    mean = totPopSize_dObs(: , 2);
    sdev = (totPopSize_dObs(: , 3).^(1/2));

    if j ==1
        fig1 = figure;
    else
        figure(fig1);
    end
    
    hold on;
    p1 = plot(tVec , totalPop0_79 , 'b-');
    p1.Color(4) = 0.6;
    hold all;
    errorbar(totPopSize_dObs(: , 1) , mean , sdev , 'ks')
    %hold all;
    %plot(tVec , totalPop0_69 , 'b-');
    hold all;
    plot(historicalPop0_69(:,1) , historicalPop0_69(:,2) , 'ro');
    hold all;
    plot(futurePop0_69(:,1) , futurePop0_69(:,2) , 'mo');
    title('KZN Population Size Ages 0-69')
    xlabel('Year'); ylabel('Individuals')
    xlim([1950 2120]);
    legend('Model prediction, ages 0-79' , 'Calibration SD, ages 0-79' , ...
        'KZN historical estimates (SSA), ages 0-69' , ... %'Model prediction, ages 0-69' , 
        'KZN future projections (UN & SSA), ages 0-69');
    hold off;

    %% Population size by broad age groups over time vs. SSA data

    % Load calibration data from Excel
    file = [pwd , '/Config/Population_validation_targets.xlsx'];
    years = xlsread(file , 'Demographics' , 'B91:F91');    % years
    kzn_popByage_yrs(: , :) = xlsread(file , 'Demographics' , 'M92:Q107').*1000;    % males and females by age in 1996-2019
    
    ageGroup = {'10-14' , '15-24' , '25-34' , '35-49' , '50-74'};
    popPropYrs = zeros(length(tVec) , 6);
    popPropYrs_obs = zeros(5 , 5);
    ageVec = {3 , [4:5] , [6:7] , [8:10] , [11:15]};
    for y = 1 : length(years)
        yearCurr = years(y);
        for aInd = 1 : length(ageVec)
            a = ageVec{aInd};
            popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , 3 : 15 , 1 : risk));
            popPropYrs(:,aInd) = sum(popVec(: , popAge),2) ./ sum(popVec(: , popTot),2);

            popPropYrs_obs(y,aInd) = sum(kzn_popByage_yrs(a , y)) / sumall(kzn_popByage_yrs(3:15 , y));  
        end
    end

    if j ==1
        fig2 = figure;
    else
        figure(fig2);
    end
    
    hold on;
    set(gca,'ColorOrderIndex',1)
    plot(tVec , popPropYrs);
    hold on;
    set(gca,'ColorOrderIndex',1)
    %set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
    plot(years , popPropYrs_obs , 'o');
    ylim([0.1 0.3]);
    xlim([1995 2020]);
    ylabel('Population proportion'); xlabel('Year'); title('KZN age distribution in broad groups');
    legend('10-14, Model' , '15-24' , '25-34' , '35-49' , '50-74' , ...
        '10-14, Observed' , '15-24' , '25-34' , '35-49' , '50-74' , 'Location' , 'EastOutside');
    %legend('Model 2019' , 'SSA KZN observed data 2019');
    %legend('Model 1919' , 'SSA KZN observed data 2019' , 'Model 1960' , ...
    %    'SSA KZN observed data 2019' ,'Model 1990' , 'SSA KZN observed data 2019' ,'Model 2019' , 'SSA KZN observed data 2019');

    %% Population size by 5-year age groups over time vs. SSA data

    % Load calibration data from Excel
    file = [pwd , '/Config/Population_validation_targets.xlsx'];
    years = xlsread(file , 'Demographics' , 'B91:F91');    % years
    kzn_popByage_yrs(: , :) = xlsread(file , 'Demographics' , 'M92:Q107').*1000;    % males and females by age in 1996-2019

    % Calibration error bars
    mean = [popAgeDist_dObs(1:16 , 2) , popAgeDist_dObs(17:32 , 2) , popAgeDist_dObs(33:48 , 2)]';
    sdev = ([popAgeDist_dObs(1:16 , 3) , popAgeDist_dObs(17:32 , 3) , popAgeDist_dObs(33:48 , 3)]'.^(1/2));

    popPropYrs = zeros(5,age);
    popPropYrs_obs = zeros(5,age);
    for y = 1 : length(years)
        yearCurr = years(y);
        for a = 1 : age
            popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
            popPropYrs(y,a) = sum(popVec(((yearCurr - startYear) * stepsPerYear +1) , popAge),2) ./ sum(popVec(((yearCurr - startYear) * stepsPerYear +1) , popTot),2);

            popPropYrs_obs(y,a) = sum(kzn_popByage_yrs(a , y)) / sumall(kzn_popByage_yrs(1 : end , y));
        end
    end

    if j ==1
        fig3 = figure;
    else
        figure(fig3);
    end
    
    subplot(1,3,1);
    plot(years , popPropYrs(: , 1:7));
    set(gca,'ColorOrderIndex',1)
    %set(gca, 'ColorOrder', circshift(get(gca, 'ColorOrder'), numel(h)))
    hold on;
    plot(years , popPropYrs_obs(: , 1:7) , 'o');
    hold on;
    calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , ...
        unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) ,];
    errorbar(calibYrs , mean(: , 1:7) , sdev(: , 1:7) , 'ks')
    ylim([0.05 0.15]);
    ylabel('Population proportion by age'); xlabel('Year');
    legend('0-4, Model' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' , '30-34' , ...
        '0-4, Observed' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' , '30-34' , ...
        'Calibration SD' , 'Location' , 'EastOutside');

    subplot(1,3,2);
    plot(years , popPropYrs(: , 8:14));
    set(gca,'ColorOrderIndex',1)
    hold on;
    plot(years , popPropYrs_obs(: , 8:14) , 'o');
    hold on;
    calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , ...
        unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) ,];
    errorbar(calibYrs , mean(: , 8:14) , sdev(: , 8:14) , 'ks')
    ylim([0.0 0.1]);
    ylabel('Population proportion by age'); xlabel('Year');
    legend('35-39, Model' , '40-44' , '45-49' , '50-54' , '55-59' , '60-64' , '65-69' , ...
        '35-39, Observed' , '40-44' , '45-49' , '50-54' , '55-59' , '60-64' , '65-69' , ...
        'Calibration SD' , 'Location' , 'EastOutside');

    subplot(1,3,3);
    set(gca,'ColorOrderIndex',1)
    plot(years , popPropYrs(: , 15:16));
    set(gca,'ColorOrderIndex',1)
    hold on;
    plot(years , popPropYrs_obs(: , 15:16) , 'o');
    hold on;
    calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1))];
    errorbar(calibYrs , mean(: , 15:16) , sdev(: , 15:16) , 'ks')
    ylim([0.0 0.02]);
    ylabel('Population proportion by age'); xlabel('Year'); %title('KZN age distribution in 5-year groups');
    legend('70-74, Model' , '75-79' , '70-74, Observed' , '75-79' , ...
        'Calibration SD' , 'Location' , 'EastOutside');




    %% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

    %% HIV prevalence by age over time vs. Africa Centre data
    hivAge = zeros(length(tVec) , age);
    ageGroup = {'15 - 19' , '20 - 24' , '25 - 29' ,...
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

    hivM(: , 1) = hivPrevM_dObs(: , 2) .* 100; % mean
    hivM(: , 2) = (hivPrevM_dObs(: , 3).^(1/2)) .* 100; % calibration SD
    hivF(: , 1) = hivPrevF_dObs(: , 2) .* 100; % mean
    hivF(: , 2) = (hivPrevF_dObs(: , 3).^(1/2)) .* 100; % calibration SD
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
        
        if j == 1
            if g == 1 
                fig4 = figure;
            elseif g == 2
                fig5 = figure;
            end
        else
            if g == 1 
                figure(fig4);
            elseif g == 2
                figure(fig5);
            end
        end
        
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
                p1 = plot(tVec' , hivAgeRel , 'b-'); 
                p1.Color(4) = 0.6;
                hold on;
                plot(prevYears , hivPrevs(((a-3) : 7 : end) , 1) , 'ro' , ...
                    prevYears2 , hivPrevs2((a-3) : 7 : end) , 'bo');
                hold on;
                errorbar(prevYears , hivPrevs(((a-3) : 7 : end) , 1) , hivPrevs(((a-3) : 7 : end) , 2) , 'ks')
            else
                p1 = plot(tVec' , hivAgeRel , 'b-');
                p1.Color(4) = 0.6;
            end
            xlabel('Year'); ylabel('Prevalence (%)'); title([gen{g} , 's ages ' , ageGroup{a-3}]) % , ' HIV Prevalence'])
            xlim([1980 2019])
        end
        legend('Model' , 'Observed KZN: AHRI DHHS' , 'Observed KZN: AHRI DHHS, VALIDATION' , 'Observed KZN: mean, SD')
    end

    %% HIV prevalence by gender over time vs. Africa Centre data
    if j ==1
        fig6 = figure;
    else
        figure(fig6);
    end
    
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
        hold on;
        p1 = plot(tVec' , hivPopPrev , 'b');
        p1.Color(4) = 0.6;
        hold all;
        plot(prevYears , hivData(:,:,g) , 'ro');
        xlabel('Year'); ylabel('Prevalence (%)'); title(gen{g});
        xlim([1980 2020])
        legend('Model' , 'Africa Center Data (Calibration)')
    end

%% HIV incidence by gender --> ***VALIDATION***
if j ==1
    fig13 = figure;
else
    figure(fig13);
end

hivInc_obs(: , : , 1) = [2005 2.14 1.57 2.93; % AHRI KZN: (Vandormael, 2019)
                         2006 2.24 1.69 2.96;
                         2007 2.30 1.74 3.05;
                         2008 2.35 1.78 3.09;
                         2009 2.45 1.85 3.24;
                         2010 2.45 1.85 3.25;
                         2011 2.30 1.70 3.11;
                         2012 2.49 1.83 3.37;
                         2013 2.22 1.64 3.01;
                         2014 1.83 1.29 2.59;
                         2015 1.39 0.94 2.07;
                         2016 1.24 0.79 1.95;
                         2017 1.01 0.58 1.76];
hivInc_obs(: , : , 2) = [2005 4.08 3.40 4.90;
                         2006 4.45 3.77 5.27;
                         2007 4.56 3.86 5.39;
                         2008 4.58 3.89 5.40;
                         2009 4.58 3.85 5.44;
                         2010 4.72 3.98 5.61;
                         2011 4.59 3.85 5.47;
                         2012 4.95 4.14 5.92;
                         2013 4.85 4.05 5.81;
                         2014 4.89 4.09 5.84;
                         2015 4.31 3.58 5.20;
                         2016 3.74 3.04 4.61;
                         2017 3.06 2.38 3.94];

gen = {'Male' , 'Female'};
for g = 1 : gender
    hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , g , 4 : 10 , 1 : risk));
    hivSus = annlz(sum(popVec(1:end-1 , hivSusInds) , 2)) ./ stepsPerYear;
    hivInc = annlz(sum(sum(sum(sum(sum(newHiv(1:end-1 , : , : , : , g , 4 : 10 , ...
        1 : risk), 2), 3), 4), 6), 7)) ./ hivSus * 100;

    subplot(1,2,g)
    plot(tVec(1:stepsPerYear:end-1) , hivInc , 'b-');
    hold all;
    errorbar(hivInc_obs(: , 1 , g) , hivInc_obs(: , 2 , g) , ...
        hivInc_obs(: , 2 , g) - hivInc_obs(: , 3 , g) , hivInc_obs(: , 4 , g) - hivInc_obs(: , 2 , g) , ...
        'rs' , 'LineWidth' , 1.5);
    xlabel('Year'); ylabel('Incidence per 100'); title(gen{g});
    xlim([1980 2020])
    if g == 1 
        legend('Model: ages 15-49' , '(Vandormael, 2019) Observed KZN, ages 15-54: 95% CI');
    elseif g == 2
        legend('Model: ages 15-49' , '(Vandormael, 2019) Observed KZN, ages 15-49: 95% CI');
    end
end

%% HIV incidence by age and gender vs. Barnighausen data --> ***VALIDATION***
    hivAge = zeros(length(tVec) , age);
    ageGroup = {'15 - 19' , '20 - 24' , '25 - 29' ,...
        '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
        '60 - 64' , '65 - 69' , '70 - 74'};

    hivMinc = [2.2 3.8 9.1 10.7 8.1 7.6 6.3 6.9]';
    hivFinc = [4.7 8.8 12.5 10.4 9.2 7.5 6.7]';
    
    gen = {'Male' , 'Female'};
    for g = 1 : gender
        hivPrevs = hivMinc;
        if g == 2
            hivPrevs = hivFinc;
        end
        
        if j == 1
            if g == 1 
                fig14 = figure;
            elseif g == 2
                fig15 = figure;
            end
        else
            if g == 1 
                figure(fig14);
            elseif g == 2
                figure(fig15);
            end
        end
        
        aVec = {16:20,21:25,26:30,31:35,36:40,41:45,46:50}; %,51:55,56:60,61:65,66:70,71:75};
        for a = 4 : 10
            hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , g , a , 1 : risk));
            hivSus = annlz(sum(popVec(1:end-1 , hivSusInds) , 2)) ./ stepsPerYear;
            hivInc = annlz(sum(sum(sum(sum(newHiv(1:end-1 , : , : , : , g , a , ...
                1 : risk), 2), 3), 4), 7)) ./ hivSus * 100;
            subplot(3 , 3 , a-3)
            if ((g == 1) && (a <= 11)) || ((g == 2) && (a <= 10))
                p1 = plot(tVec(1:stepsPerYear:end-1) , hivInc , 'b-'); 
                p1.Color(4) = 0.6;
                hold on;
                plot(2004 , hivPrevs((a-3) , 1) , 'ro')
            else
                p1 = plot(tVec(1:stepsPerYear:end-1) , hivInc , 'b-');
                p1.Color(4) = 0.6;
            end
            xlabel('Year'); ylabel('Incidence per 100'); title([gen{g} , 's ages ' , ageGroup{a-3}])
            xlim([1980 2020])
        end
        legend('Model' , 'Observed KZN: Barnighausen et al. - VALIDATION')
    end





    %% ********************************** HPV FIGURES **********************************************************************************************

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

    % McDonald 2014 (not sure where these bounds came from)
    hpvHivObs(: , 1) = hpv_hiv_dObs(: , 2); % mean
    hpvHivObs(: , 2) = [0.63 0.54 0.54 0.47 0.42 0.34 0.32 0.35 0.16]'; % lb
    hpvHivObs(: , 3) = [0.87 0.67 0.66 0.62 0.51 0.50 0.55 0.72 0.53]'; % ub
    hpvNegObs(: , 1) = hpv_hivNeg_dObs(: , 2); % mean
    hpvNegObs(: , 2) = [0.53 0.34 0.21 0.17 0.18 0.16 0.11 0.14 0.12]'; % lb
    hpvNegObs(: , 3) = [0.67 0.41 0.27 0.23 0.21 0.20 0.15 0.19 0.18]'; % ub
    hpvHivObs = hpvHivObs * 100;
    hpvNegObs = hpvNegObs * 100;

    % Calibration error bars
    mean = hpv_hiv_dObs(: , 2) .* 100;
    sdev = (hpv_hiv_dObs(: , 3).^(1/2)) .* 100;
    meanNeg = hpv_hivNeg_dObs(: , 2) .* 100;
    sdevNeg = (hpv_hivNeg_dObs(: , 3).^(1/2)) .* 100;

    if j == 1
        fig7 = figure;
    else
        figure(fig7);
    end
    
    subplot(2,1,1);
    hold on;
    % plot(1 : length(hpv2002) , hpv2002 , 'o-')
    % hold all;
    p1 = plot(1 : length(hpvHIV2002) , hpvHIV2002 , 'b-');
    p1.Color(4) = 0.6;
    hold on;
    % general
    % yPosError = abs(hrHpvObs(: , 3) - hrHpvObs(: , 1));
    % yNegError = abs(hrHpvObs(: , 2) - hrHpvObs(: , 1));
    % errorbar(1 : length(hrHpvObs) , hrHpvObs(: , 1) , yNegError , yPosError , 'rs')
    % HIV+
    yPosError = abs(hpvHivObs(: , 3) - hpvHivObs(: , 1));
    yNegError = abs(hpvHivObs(: , 2) - hpvHivObs(: , 1));
    %errorbar(1 : length(hpvHivObs) , hpvHivObs(: , 1) , yNegError , yPosError , 'bs')
    if j == 1
        errorbar(1 : length(mean) , mean , sdev , 'ks');
    end
    set(gca , 'xtickLabel' , ageGroup);
    set(gca , 'xtick' , 1 : length(hpvNegObs) , 'xtickLabel' , ageGroup);
    xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
    %legend('HIV-Positive (year 2002)' , 'Observed HIV-Positive: Calibration SD');    %'Observed HIV-Positive: McDonald 2014 bounds ?' , 
    legend(fileInds{1}, 'Observed HIV-Positive: Calibration SD' , fileInds{2:end});
    
    subplot(2,1,2);
    hold on;
    p1 = plot(1 : length(hpvNeg2002) , hpvNeg2002 , 'r-');
    p1.Color(4) = 0.6;
    hold on;
    %HIV-
    yPosError = abs(hpvNegObs(: , 3) - hpvNegObs(: , 1));
    yNegError = abs(hpvNegObs(: , 2) - hpvNegObs(: , 1));
    %errorbar(1 : length(hpvNegObs) , hpvNegObs(: , 1) , yNegError , yPosError , 'rs')
    if j == 1
        errorbar(1 : length(mean) , meanNeg , sdevNeg , 'ks');
    end
    set(gca , 'xtickLabel' , ageGroup);
    set(gca , 'xtick' , 1 : length(hpvNegObs) , 'xtickLabel' , ageGroup);
    xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
    %legend('HIV-Negative (year 2002)' , 'Observed HIV-Negative: Calibration SD');    %'Observed HIV-Negative: McDonald 2014 bounds ?' , 
    legend(fileInds{1} , 'Observed HIV-Negative: Calibration SD' , fileInds{2:end});
    %title('Age Specific hrHPV Prevalence in 2002')

    %% HPV prevalence by age and HIV status in 2008 vs. Mbulawa data
    yearPrev = 2008;
    ageGroup = {'15-24' , '25-34' , '35-44' , '45-64'};
    ageVec = {[4:5],[6:7],[8:9],[10:13]};
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

    % Calibration error bars
    mean = hpv_hivM2008_dObs(: , 2) .* 100;
    sdev = (hpv_hivM2008_dObs(: , 3).^(1/2)) .* 100;
    meanNeg = hpv_hivMNeg2008_dObs(: , 2) .* 100;
    sdevNeg = (hpv_hivMNeg2008_dObs(: , 3).^(1/2)) .* 100;

    if j ==1
        fig8 = figure;
    else
        figure(fig8);
    end
    
    subplot(2,1,1)
    hold on;
    p1 = plot([1 : length(ageVec)] , hpv_hivM2008(: , 1)' , 'b-');
    p1.Color(4) = 0.6;
    hold all;
    plot([1 : length(ageVec)] , hpv_hivM2008_dObs(: , 2)' .* 100 , 'bo');
    hold all;
    errorbar(1 : length(mean) , mean , sdev , 'ks')
    set(gca , 'xtick' , [1 : length(ageVec)] , 'xtickLabel' , ageGroup);
    legend('HIV-Positive Males (year 2008)' , 'Observed HIV-Positive Males: Mbulawa 2008' , ...
        'Observed HIV-Positive Males: Calibration SD');
    xlabel('Age Group'); ylabel('hrHPV Prevalence (%)'); ylim([0 100]);

    subplot(2,1,2)
    hold on;
    p1 = plot([1 : length(ageVec)] , hpv_hivM2008(: , 2)' , 'r-')
    p1.Color(4) = 0.6;
    hold all;
    plot([1 : length(ageVec)] , hpv_hivMNeg2008_dObs(: , 2)' .* 100 , 'ro');
    hold all;
    errorbar(1 : length(mean) , meanNeg , sdevNeg , 'ks')
    set(gca , 'xtick' , [1 : length(ageVec)] , 'xtickLabel' , ageGroup);
    legend('HIV-Negative Males (year 2008)' , 'Observed HIV-Negative Males: Mbulawa 2008' , ...
        'Observed HIV-Negative Males: Calibration SD');
    xlabel('Age Group'); ylabel('hrHPV Prevalence (%)'); ylim([0 100]);

    %% ********************************** CIN FIGURES *********************************************************************************************

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
    
    % McDonald 2014, HIV-positive (not sure where these bounds came from)
    cinPosAct(: , 1) = cinPos2002_dObs(: , 2); % mean
    cinPosAct(: , 2) = [0.03 0.02 0.09 0.10 0.05 0.02 0.02 0.00 0.00 0.00]'; % lb
    cinPosAct(: , 3) = [0.22 0.08 0.17 0.21 0.11 0.09 0.14 0.17 0.22 0.22]'; % ub
    cinPosAct = cinPosAct .* 100; % convert to %
    yPosError = abs(cinPosAct(: , 3) - cinPosAct(: , 1));
    yNegError = abs(cinPosAct(: , 2) - cinPosAct(: , 1));
    
    % Calibration error bars
    mean = cinPos2002_dObs(: , 2) .* 100;
    sdev = (cinPos2002_dObs(: , 3).^(1/2)) .* 100;
    
    if j ==1
        fig9 = figure;
    else
        figure(fig9);
    end
    
    subplot(2 , 1 , 1);
    p1 = plot(1 : length(cinPos2002) , cinPos2002 ,'r-');
    p1.Color(4) = 0.6;
    %hold on;
    %errorbar(1 : length(cinPosAct) , cinPosAct(: , 1) , yNegError , yPosError , 'rs')
    hold on;
    if j == 1
        errorbar(1 : length(mean) , mean , sdev , 'ks')
    end
    %legend('HR HPV CIN 2/3' , 'Calibration SD'); %'McDonald 2014 bounds ?' ,
    legend(fileInds{1} , 'Calibration SD' , fileInds{2:end});
    set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
    xlabel('Age Group'); ylabel('Prevalence (%)')
    title('Age Specific CIN 2/3 Prevalence Among HIV+ in 2002')
    ylim([0 25])
    
    % McDonald 2014, HIV-negative (not sure where these bounds came from)
    cinNegAct(: , 1) = cinNeg2002_dObs(: , 2); % mean
    cinNegAct(: , 2) = [0.00 0.02 0.01 0.02 0.02 0.02 0.02 0.01 0.00 0.00]'; % lb
    cinNegAct(: , 3) = [0.03 0.04 0.03 0.05 0.04 0.04 0.04 0.03 0.03 0.03]'; % ub
    cinNegAct = cinNegAct .* 100; % convert to %
    yPosError = abs(cinNegAct(: , 3) - cinNegAct(: , 1));
    yNegError = abs(cinNegAct(: , 2) - cinNegAct(: , 1));
    
    % Calibration error bars
    mean = cinNeg2002_dObs(: , 2) .* 100;
    sdev = (cinNeg2002_dObs(: , 3).^(1/2)) .* 100;
    
    subplot(2 , 1 , 2)
    p1 = plot(1 : length(cinNeg2002) , cinNeg2002 , 'b-');
    p1.Color(4) = 0.6;
    %hold on;
    %errorbar(1 : length(cinNegAct) , cinNegAct(: , 1) , yNegError , yPosError , 'rs')
    hold on;
    if j == 1
        errorbar(1 : length(mean) , mean , sdev , 'ks')
    end
    %legend('HR HPV CIN 2/3' , 'Calibration SD'); %'McDonald 2014 bounds ?' ,
    legend(fileInds{1} , 'Calibration SD' , fileInds{2:end});
    set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
    xlabel('Age Group'); ylabel('Prevalence (%)')
    title('Age Specific CIN 2/3 Prevalence Among HIV- in 2002')
    ylim([0 25])
    
    %% CIN1, CIN2, CIN3 prevalence for All HR HPV types combined by HIV status in 2015 vs. Kuhn 2015 data (calibration)
    if j ==1
        fig16 = figure;
    else
        figure(fig16);
    end
    
    % HIV-positive (on and not on ART)
    cin1Inds = unique([toInd(allcomb(3 : 8 , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 7 : 13 , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 3 , 7] , 3 , 1 , 1 : intervens , 2 , 7 : 13 , 1 : risk))]);
    cin2Inds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 , [1 : 4 , 7] , ...
        1 , 1 : intervens , 2 , 7 : 13 , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 4 , 7] , 4 , 1 , 1 : intervens , 2 , 7 : 13 , 1 : risk))]);
    cin3Inds = unique([toInd(allcomb(3 : 8 , 1 : viral , 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 7 : 13 , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 5 , 1 , 1 : intervens , 2 , 7 : 13 , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 7 : 13 , 1 : risk));
    cin1Pos2015 = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin1Inds)))...
        ./ sum(popVec((2015 - startYear) * stepsPerYear +1 , ageInds));
    cin2Pos2015 = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin2Inds)))...
        ./ sum(popVec((2015 - startYear) * stepsPerYear +1 , ageInds));
    cin3Pos2015 = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin3Inds)))...
        ./ sum(popVec((2015 - startYear) * stepsPerYear +1 , ageInds));
    % HIV-negative
    cin1NegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 7 : 13 , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 3 , 7] , 3 , 1 , 1 : intervens , 2 , 7 : 13 , 1 : risk))]);
    cin2NegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 , [1 : 4 , 7] , ...
        1 , 1 : intervens , 2 , 7 : 13 , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 4 , 7] , 4 , 1 , 1 : intervens , 2 , 7 : 13 , 1 : risk))]);
    cin3NegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 7 : 13 , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 5 , 1 , 1 : intervens , 2 , 7 : 13 , 1 : risk))]);
    ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 7 : 13 , 1 : risk));
    cin1Neg2015 = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin1NegInds)))...
        ./ (sum(popVec((2015 - startYear) * stepsPerYear +1 , ageNegInds)));
    cin2Neg2015 = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin2NegInds)))...
        ./ (sum(popVec((2015- startYear) * stepsPerYear +1 , ageNegInds)));
    cin3Neg2015 = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin3NegInds)))...
        ./ (sum(popVec((2015 - startYear) * stepsPerYear +1 , ageNegInds)));
    
    % Calibration error bars
    cinPos_mean = cinPos2015_dObs(: , 2) .* 100;
    cinPos_sdev = (cinPos2015_dObs(: , 3).^(1/2)).*2 .* 100;
    cinNeg_mean = cinNeg2015_dObs(: , 2) .* 100;
    cinNeg_sdev = (cinNeg2015_dObs(: , 3).^(1/2)).*2 .* 100;

    subplot(1 , 2 , 1);
    if j == 1
        errorbar(1 : length(cinPos_mean) , cinPos_mean , cinPos_sdev , ...
            'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
        hold all;
    end
    hold all;
    plot([1 2 3] , [(cin1Pos2015 * 100) (cin2Pos2015 * 100) (cin3Pos2015 * 100)] , 'k-' );
%         legend('(Kuhn, 2020) Observed CIN1,2,3 Cape Town: mean, 2SD' , 'Model, CIN1: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
%         'Model, CIN2: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
%         'Model, CIN3: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
    legend('(Kuhn, 2020) Observed CIN1,2,3 Cape Town: mean, 2SD' , ...
        fileInds{:});
    set(gca , 'xtick' , 1 : 3 , 'xtickLabel' , {'CIN1' , 'CIN2' , 'CIN3'});
    ylabel('Prevalence (%)')
    title('Prevalence in 2015 - Females aged 30-65, HIV+')
    ylim([0 25])
    grid on;

    subplot(1 , 2 , 2)
    if j == 1
        errorbar(1 : length(cinNeg_mean) , cinNeg_mean , cinNeg_sdev , ...
            'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
        hold all;
    end
    hold all;
    plot([1 2 3] , [(cin1Neg2015 * 100) (cin2Neg2015 * 100) (cin3Neg2015 * 100)] , 'k-' );
%         legend('(Kuhn, 2020) Observed CIN1,2,3 Cape Town: mean, 2SD' , 'Model, CIN1: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
%         'Model, CIN2: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
%         'Model, CIN3: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
    legend('(Kuhn, 2020) Observed CIN1,2,3 Cape Town: mean, 2SD' , ...
        fileInds{:});
    set(gca , 'xtick' , 1 : 3 , 'xtickLabel' , {'CIN1' , 'CIN2' , 'CIN3'});
    ylabel('Prevalence (%)')
    title('Prevalence in 2015 - Females aged 30-65, HIV-')
    ylim([0 10])
    grid on;

    %% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

    %% Cervical cancer incidence in 2012 by age vs. Globocan 2012 data and other sources
    %ccIncYears = [2017 , 2003 , 1994 , 2012];
    ccIncYears = [2012];
    ccCrude = zeros(age , length(ccIncYears));
    ccAgeNegRel = ccCrude;
    ccAgePosRel = zeros(age , 5 , length(ccIncYears));
    ccArtRel = ccCrude;
    ccNegPosArt = zeros(age , 3 , length(ccIncYears));
    fScale = 10^5;
    ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
        '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
        '60-64' , '65-69' , '70-74' , '75-79'};
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
            ccCrude(a , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end , ...
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
    
    globocan = [0.00
    2.646467154
    8.848389036
    45.1937379
    53.40682334
    63.4
    68.3
    70.7
    73
    77.4
    82.7
    88.6
    95.2];
    
    combined_ub = [0.00
    2.65
    8.85
    45.19
    53.41
    67.05
    80.83
    78.97
    128.87
    105.27
    118.70
    111.81
    95.20];
    
    combined_lb = [0.00
    0.00
    0.41
    9.97
    12.61
    25.00
    45.69
    36.31
    50.55
    57.08
    62.69
    42.43
    52.01];
    
    % Calibration error bars
    mean = ccInc2012_dObs(: , 2);
    sdev = (ccInc2012_dObs(: , 3).^(1/2));
    
    for y = 1 : length(ccIncYears)
        ccIncYear = ccIncYears(y);
        
        if j ==1
            fig10 = figure;
        else
            figure(fig10);
        end     
    
        % Plot model outputs
        p1 = plot(1 : size(ccCrude , 1) , ccCrude(: , y) , 'k-'); % , 1 : size(ccAgeNegRel(: , y) , 1) , ...
            %ccAgeNegRel(: , y) , '-kp' , 1 : size(ccAgePosRel , 1) , ccAgePosRel(: , 5 , y) , '-k+' , ...
            %1 : size(ccArtRel , 1) , ccArtRel(: , y) , '-k^');
        p1.Color(4) = 0.6;
        hold on;
        % Plot observed data
        if j == 1
            plot(4 : age , combined_ub , 'r--' , 4 : age , combined_lb , 'r--' , 4 : age , ccInc2012_dObs(: , 2) , 'r-');
            hold on;
            errorbar(4 : age , mean , sdev , 'rs')
        end
        xlabel('Age Group'); ylabel('Incidence per 100,000');
        set(gca , 'xtick' , 1 : length(ccCrude) , 'xtickLabel' , ageGroup);
        title(['Cervical Cancer Incidence in ' num2str(ccIncYear)]);
        %legend('General' , 'Combined SA: upper bound' , 'Combined SA: lower bound' , ...
        %    'Globocan SA (2012)' , 'Calibration SD'); % 'HIV-' , 'HIV+' , 'ART' , 
        legend(fileInds{1} , 'Combined SA: upper bound' , 'Combined SA: lower bound' , ...
            'Globocan SA (2012)' , 'Calibration SD' , fileInds{2:end});
        % legend('General' , 'HIV-' , ' Acute and CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , ...
        %     'CD4 < 200' , 'ART' , 'Globocan' , 'Upper Bound' , 'Lower Bound' , ...
        %     'Location' , 'NorthEastOutside')
    
    %     figure()
    %     bar(1 : length(ccNegPosArt(: , : , y)) , ccNegPosArt(: , : , y), 'stacked')
    %     xlabel('Age Group'); ylabel('Incidence per 100,000')
    %     set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
    %     title(['Cervical Cancer Incidence Distribution in ' , num2str(ccIncYear)])
    %     legend('HIV-' , 'HIV+' , 'ART')
    end
    
    %% Cervical cancer incidence in 2018 by age vs. Globocan 2018 data
    ccIncYears = [2018];
    ccCrude = zeros(age , length(ccIncYears));
    ccAgeNegRel = ccCrude;
    ccAgePosRel = zeros(age , 5 , length(ccIncYears));
    ccArtRel = ccCrude;
    ccNegPosArt = zeros(age , 3 , length(ccIncYears));
    fScale = 10^5;
    ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
        '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
        '60-64' , '65-69' , '70-74' , '75-79'};
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
            ccCrude(a , y) = annlz(sum(sum(sum(newCC(yr_start : yr_end , ...
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
    
    % Calibration error bars
    mean = ccInc2018_dObs(: , 2);
    sdev = (ccInc2018_dObs(: , 3).^(1/2));
    
    for y = 1 : length(ccIncYears)
        ccIncYear = ccIncYears(y);
        
        if j ==1
            fig17 = figure;
        else
            figure(fig17);
        end     
    
        % Plot model outputs
        p1 = plot(1 : size(ccCrude , 1) , ccCrude(: , y) , 'k-'); % , 1 : size(ccAgeNegRel(: , y) , 1) , ...
            %ccAgeNegRel(: , y) , '-kp' , 1 : size(ccAgePosRel , 1) , ccAgePosRel(: , 5 , y) , '-k+' , ...
            %1 : size(ccArtRel , 1) , ccArtRel(: , y) , '-k^');
        p1.Color(4) = 0.6;
        hold on;
        % Plot observed data
        if j == 1
            plot(4 : age-1 , ccInc2018_dObs(: , 2) , 'r-');
            hold on;
            errorbar(4 : age-1 , mean , sdev , 'rs')
        end
        xlabel('Age Group'); ylabel('Incidence per 100,000');
        set(gca , 'xtick' , 1 : length(ccCrude) , 'xtickLabel' , ageGroup);
        title(['Cervical Cancer Incidence in ' num2str(ccIncYear)]);
        %legend('General' , 'Combined SA: upper bound' , 'Combined SA: lower bound' , ...
        %    'Globocan SA (2012)' , 'Calibration SD'); % 'HIV-' , 'HIV+' , 'ART' , 
        legend(fileInds{1} , ...
            'Globocan SA (2018)' , 'Calibration SD' , fileInds{2:end});
        % legend('General' , 'HIV-' , ' Acute and CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , ...
        %     'CD4 < 200' , 'ART' , 'Globocan' , 'Upper Bound' , 'Lower Bound' , ...
        %     'Location' , 'NorthEastOutside')
    
    %     figure()
    %     bar(1 : length(ccNegPosArt(: , : , y)) , ccNegPosArt(: , : , y), 'stacked')
    %     xlabel('Age Group'); ylabel('Incidence per 100,000')
    %     set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
    %     title(['Cervical Cancer Incidence Distribution in ' , num2str(ccIncYear)])
    %     legend('HIV-' , 'HIV+' , 'ART')
    end
    
    %% Cervical cancer incidence over time
    ccIncYears = [2011 , 2012 , 2013 , 2014 , 2015 , 2016 , 2017 , 2018 , 2019];
    ccCrude = zeros(1,length(ccIncYears));
    ccCrudeObs = zeros(1,length(ccIncYears));
    ccCrudeObs_2018 = zeros(1,length(ccIncYears));
    fScale = 10^5;
    ccYrs = ((ccIncYears - startYear) * stepsPerYear : ...
        (ccIncYears + 1 - startYear) * stepsPerYear);
    
    globocan = [0.00
    0.00
    2.646467154
    8.848389036
    45.1937379
    53.40682334
    63.4
    68.3
    70.7
    73
    77.4
    82.7
    88.6
    95.2];

    globocan2018 = [0.00
    3.91
    19.72
    35.82
    53.75
    71.18
    85.09
    95.44
    96.10
    94.85
    97.49
    101.99
    110.43
    110.43]; %set to be the same as 70-74 age group
    
    
    for y = 1 : length(ccIncYears)
        % Year
        yr_start = (ccIncYears(y) - 1 - startYear)  .* stepsPerYear;
        yr_end = (ccIncYears(y) - startYear) .* stepsPerYear - 1;

        % Total population
        ageIndsAll = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
            1 : intervens , 2 , 3 : age , 1 : risk));
        ccCrude(1,y) = annlz(sum(sum(sum(newCC(yr_start : yr_end , ...
            1 : disease , 3 : age , :) , 2) , 3) , 4)) ...
            ./ (annlz(sum(popVec(yr_start : yr_end , ageIndsAll) , 2)) ...
            ./ stepsPerYear) * fScale;
        
        ccCases = 0;
        for a = 3 : age
            ageInds = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
            1 : intervens , 2 , a , 1 : risk));
            ccCases = ccCases + ...
                (annlz(sum(popVec(yr_start : yr_end , ageInds) , 2)) ./ stepsPerYear) ...
                .* globocan(a-2 , 1);
        end
        ccCrudeObs(1,y) = (ccCases ./ ...
            (annlz(sum(popVec(yr_start : yr_end , ageIndsAll) , 2)) ...
            ./ stepsPerYear));
        
        ccCases2018 = 0;
        for a = 3 : age
            ageInds = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
            1 : intervens , 2 , a , 1 : risk));
            ccCases2018 = ccCases2018 + ...
                (annlz(sum(popVec(yr_start : yr_end , ageInds) , 2)) ./ stepsPerYear) ...
                .* globocan2018(a-2 , 1);
        end
        ccCrudeObs_2018(1,y) = (ccCases2018 ./ ...
            (annlz(sum(popVec(yr_start : yr_end , ageIndsAll) , 2)) ...
            ./ stepsPerYear));
    end
            
    if j ==1
        fig11 = figure;
    else
        figure(fig11);
    end     

    % Plot model outputs
    hold all;
    p1 = plot(ccIncYears , ccCrude(1 , :) , 'b-');
    p1.Color(4) = 0.6;
    hold all;
    plot(2012 , ccCrudeObs(1,2) , 'ok')
    hold all;
    plot(2018 , ccCrudeObs_2018(1,8) , 'ok')
    xlabel('Year'); ylabel('Incidence per 100,000');
    xlim([2011 2019]); ylim([0 200]);
    legend('Model' , 'Globocan 2012, crude' , 'Globocan 2018, crude');
    title(['Cervical cancer incidence over time']);

        
    
    

    %% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************

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

    % Calibration error bars
    mean = hpv_dist_dObs(: , 2) .* 100;
    sdev = (hpv_dist_dObs(: , 3).^(1/2)) .* 100;

    if j ==1
        fig12 = figure;
    else
        figure(fig12);
    end
    
    subplot(2,3,1)
    p1 = plot(tVec , hpv_vax , 'k');
    p1.Color(4) = 0.6;
    hold all;
    plot(tVec((2011 - startYear) * stepsPerYear) , 46.82 , 'ko');
    hold all;
    p1 = plot(tVec , hpv_nonVax , 'b');
    p1.Color(4) = 0.6;
    hold all;
    plot(tVec((2011 - startYear) * stepsPerYear) , 53.18 , 'bo');
    hold all;
    errorbar([2011, 2011] , mean , sdev , 'ks');
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
    
    % Calibration error bars
    mean = cin1_dist_dObs(: , 2) .* 100;
    sdev = (cin1_dist_dObs(: , 3).^(1/2)) .* 100;
    
    subplot(2,3,2)
    p1 = plot(tVec , cin_vax , 'k');
    p1.Color(4) = 0.6;
    hold all;
    plot(tVec((2011 - startYear) * stepsPerYear) , 51.92 , 'ko');
    hold all;
    p1 = plot(tVec , cin_nonVax ,'b');
    p1.Color(4) = 0.6;
    hold all;
    plot(tVec((2011 - startYear) * stepsPerYear) , 48.08 , 'bo');
    hold all;
    errorbar([2011, 2011] , mean , sdev , 'ks');
    ylim([0 100]);
    xlim([2000 2015]);
    xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
    title('CIN1');
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
    
    subplot(2,3,3);
    p1 = plot(tVec , cin_vax , 'k');
    p1.Color(4) = 0.6;
    hold all;
    plot(tVec((2011 - startYear) * stepsPerYear) , 62.81 , 'ko');
    hold all;
    p1 = plot(tVec , cin_nonVax ,'b');
    p1.Color(4) = 0.6;
    hold all;
    plot(tVec((2011 - startYear) * stepsPerYear) , 37.19 , 'bo');
    ylim([0 100]);
    xlim([2000 2015]);
    xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
    title('CIN2');
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
    
    % Calibration error bars
    mean = cin3_dist_dObs(: , 2) .* 100;
    sdev = (cin3_dist_dObs(: , 3).^(1/2)) .* 100;
    
    subplot(2,3,4)
    p1 = plot(tVec , cin_vax , 'k');
    p1.Color(4) = 0.6;
    hold all;
    plot(tVec((2011 - startYear) * stepsPerYear) , 73.71 , 'ko');
    hold all;
    p1 = plot(tVec , cin_nonVax ,'b');
    p1.Color(4) = 0.6;
    hold all;
    plot(tVec((2011 - startYear) * stepsPerYear) , 26.29 , 'bo');
    hold all;
    errorbar([2011, 2011] , mean , sdev , 'ks');
    ylim([0 100]);
    xlim([2000 2015]);
    xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
    title('CIN3');
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
    
    % Calibration error bars
    mean = cc_dist_dObs(: , 2) .* 100;
    sdev = (cc_dist_dObs(: , 3).^(1/2)) .* 100;
    
    subplot(2,3,5)
    p1 = plot(tVec , cc_vax , 'k');
    p1.Color(4) = 0.6;
    hold all;
    plot(tVec((2011 - startYear) * stepsPerYear) , 85.78 , 'ko');
    hold all;
    p1 = plot(tVec , cc_nonVax ,'b');
    p1.Color(4) = 0.6;
    hold all;
    plot(tVec((2011 - startYear) * stepsPerYear) , 14.22 , 'bo');
    hold all;
    errorbar([2011, 2011] , mean , sdev , 'ks');
    ylim([0 100]);
    xlim([2000 2015]);
    xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
    title('Cervical Cancer')
    legend('9v-type HPV' , 'Observed 2011- 9v: mean, SD' , 'Non-9v-type HPV' , 'Observed 2011- non-9v: mean, SD');

end
