% function[] = showResults_multSims_CIs()

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
    c3c2Mults , c2c1Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , waning , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear , circNatStartYear , vaxStartYear , ...
    baseline , cisnet , who , whob , circProtect , condProtect , MTCTRate , ...
    hyst , OMEGA , ...
    ccInc2012_dObs , cc_dist_dObs , cin3_dist_dObs , ...
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
    infhpvNonVaxInds , ageInd , riskInd , ...
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
set(0 , 'defaultlinelinewidth' , 1.5)

% Indices of calib runs to plot
fileInds = {'3_2187' , '4_2584' , '1_3643' , '2_5021' , '1_1769' , ...
    '0_5382' , '3_2003' , '2_1876' , '3_4037' , '5_3628' , '3_3747' , ...
    '3_824' , '0_8037' , '5_293' , '2_1404' , '4_4681' , '5_2716' , ...
    '4_5427' , '1_4768' , '0_1909' , '5_2437' , '2_216' , '3_4120' , ...
    '5_2021' , '4_4099'}; % , '2_731' , '4_1001' , '3_5484' , ...
    %'4_58' , '4_1715' , '3_550' , '1_5534' , '3_14' , '2_3242' , ...
    %'4_2903' , '4_2416' , '5_1444' , '3_2813' , '5_5068'};    % 22Apr20Ph2V4
% fileInds = {'0_3292' , '0_7871' , '0_4325' , '1_3259' , '0_1474' , ...
%     '1_4214' , '0_2709' , '1_1367' , ...
%     '2_3899' , '2_4116' , '2_413' , '2_4577'};    % 22Apr20Ph2V3
% fileInds = {'11_1' , '11_2' , '11_3' , '11_4' , '11_5' , '11_6' , ...
%     '11_7' , '11_8' , '11_9' , '11_10' , '11_11' , '11_12' , '11_13' , ...
%     '11_14' , '11_15' , '11_16' , '11_17' , '11_18' , '11_19' , '11_20' , ...
%     '11_21' , '11_22' , '11_23' , '11_24' , '11_25'};
% fileInds = {'11_2946' , '6_1657' , ...
%     '10_3755' , '5_1279' , '0_6202' , '6_452' , '5_2720' , '10_1939' , ...
%     '11_2511' , '9_3353' , '7_180' , ...
%     '4_4147' , '2_2314' , '7_741' , ...
%     '10_3664' , '11_2175' , '3_2610' , ...
%     '7_240' , '10_4629' , '7_5509' , '8_4473' , '0_2709' , ...
%     '2_3468' , '9_3299' , '11_669'};  % 22Apr20Ph2V2, t11, 25 best-fitting CC inc sets
% fileInds = {'11_2946' , '11_4738' , '7_1476' , '2_1779' , '6_1657' , '11_2200' , ...
%     '10_3755' , '5_1279' , '0_6202' , '6_452' , '4_4034' , '5_2720' , '10_1939' , ...
%     '0_2605' , '11_2511' , '9_3353' , '4_3559' , '6_2571' , '10_2890' , '7_180' , ...
%     '5_5825' , '1_4859' , '4_4147' , '2_2314' , '7_741'};  % 22Apr20Ph2V2
%     '10_3270' , '10_3664' , '8_2827' , '11_2175' , '7_1426' , '2_1833' , '3_2610' , ...
%     '7_240' , '6_3972' , '10_4629' , '7_5509' , '6_1931' , '8_4473' , '0_2709' , ...
%     '7_3753' , '0_845' , '2_3468' , '8_2502' , '9_3299' , '11_669' , '9_4007' , ...
%     '8_5698' , '7_2004' , '5_33' , '8_1372'};
% fileInds = {'12_3346' , '12_2618' , '11_932' , '16_3038' , '8_597' , '12_2550' , ...  % 22Apr20Ph2, top 25 sets
%     '12_3895' , '22_487' , '8_2705' , '22_3250' , '15_2550' , '14_563' , ...
%     '4_1887' , '10_688' , '18_3391' , '14_2659' , '19_2814' , '18_903' , ...
%     '22_2697' , '4_1676' , '4_2471' , '15_2517' , '16_1709' , '12_2481' , '16_3992'};
% {'2_846' , '14_947' , '16_3127' , '12_689' , ...  % 22Apr20Ph1, top 50 sets
%     '10_2727' , '17_3986' , '16_2194' , '15_3850' , '9_334' , '0_6657' , ...
%     '16_2364' , '4_711' , '19_1017' , '4_2361' , '15_2155' , '17_594' , ...
%     '19_1779' , '11_1541' , '12_3055' , '6_746' , '20_944' , '13_3012' , ...
%     '18_387' , '17_1649' , '17_3242' , '14_2649' , '16_2701' , '20_1864' , ...
%     '19_3788' , '10_425' , '11_3176' , '18_45' , '15_532' , '20_3201' , ...
%     '5_1822' , '7_2249' , '15_2965' , '19_2304' , '10_2162' , '11_301' , ...
%     '14_1533' , '18_1043' , '0_8616' , '5_4018' , '11_3613' , '18_2578' , ...
%     '7_1869' , '15_1004' , '12_2230' , '8_3465'};
nRuns = length(fileInds);

% Initialize model output plots
% Total population size
popYearVec = unique(totPopSize_dObs(: ,1));
popSize = zeros(nRuns , length(popYearVec));
% Population age distribution
popYearVec = unique(popAgeDist_dObs(: ,1));
popProp = zeros(nRuns , length(popYearVec) , age);
% HIV prevalence
hivYearVec = [unique(hivPrevM_dObs(: ,1)) ; [2010 : 2016]'];
hivAgeM = zeros(nRuns , 7 , length(hivYearVec));
hivAgeF = hivAgeM;
hivPrev = zeros(nRuns , gender , length(hivYearVec));
hivYearVec2 = [2002 , 2005 , 2008 , 2012 , 2017];
hivPrevTot = zeros(nRuns , length(hivYearVec2));
hivPrevTot2 = zeros(nRuns , length(hivYearVec2));
hivPrevTot3 = zeros(nRuns , length(hivYearVec2));
% HIV incidence
hivIncYearVec = [2005 : 2017];
hivInc = zeros(nRuns , gender , length(hivIncYearVec));
hivIncYoung = hivInc;
hivIncOlder = hivInc;
hivIncYearVec2 = [2000 : 2017];
hivIncM = zeros(nRuns , 7 , length(hivIncYearVec2));
hivIncF = zeros(nRuns , 7 , length(hivIncYearVec2));
% Female HPV prevalence
hpv_hiv = zeros(nRuns , 9);
hpv_hivNeg = hpv_hiv;
% Male HPV prevalence
hpv_hivM2008 = zeros(nRuns , 4);
hpv_hivMNeg2008 = hpv_hivM2008;
% CIN2/3 prevalence
cinPos2002 = zeros(nRuns , 10);
cinNeg2002 = cinPos2002;
% CIN1,2,3, prevalence by age
cin1Pos2002 = cinPos2002;
cin1Neg2002 = cinPos2002;
cin2Pos2002 = cinPos2002;
cin2Neg2002 = cinPos2002;
cin3Pos2002 = cinPos2002;
cin3Neg2002 = cinPos2002;
% CIN1,2,3 prevalence
cin1Pos2015 = zeros(nRuns , 1);
cin1Neg2015 = cin1Pos2015;
cin2Pos2015 = cin1Pos2015;
cin2Neg2015 = cin1Pos2015;
cin3Pos2015 = cin1Pos2015;
cin3Neg2015 = cin1Pos2015;
% CC incidence
ccInc2012 = zeros(nRuns , age);
ccInc2012neg = ccInc2012;
ccInc2012pos = ccInc2012;
ccInc2012art = ccInc2012;
% HPV/CIN/CC type distribution
typeDistYearVec = [2010 , 2011 , 2012 , 2013 , 2014 , 2015];
cc_vax = zeros(nRuns , length(typeDistYearVec));
cc_nonVax = cc_vax;
cin3_vax = cc_vax;
cin3_nonVax = cc_vax;
cin1_vax = cc_vax;
cin1_nonVax = cc_vax;
hpv_vax = cc_vax;
hpv_nonVax = cc_vax;


resultsDir = [pwd , '\HHCoM_Results\'];
for j = 1 : nRuns
    % Load results
    pathModifier = ['toNow_22Apr20Ph2V4_noBaseVax_baseScreen_hpvHIVcalib_' , fileInds{j}];
    load([resultsDir , pathModifier])
   
    %% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

    %% Population size over time vs. Statistics South Africa data (calibration)
    popYearVec = unique(totPopSize_dObs(: ,1));
    for t = 1 : length(popYearVec)
        popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
        popSize(j , t) = sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTot),2);
    end

    %% Population size by 5-year age groups over time vs. Statistics South Africa data (calibration)
    popYearVec = unique(popAgeDist_dObs(: ,1));
    for t = 1 : length(popYearVec)
        for a = 1 : age
            popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
            popProp(j , t , a) = sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popAge),2) ./ ...
                sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTot),2);
        end
    end

    %% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

    %% HIV prevalence by age over time vs. AHRI (calibration, validation) and IHME model data (validation)
    hivYearVec = [unique(hivPrevM_dObs(: ,1)) ; [2010 : 2016]'];
    
    for t = 1 : length(hivYearVec)
        for a = 4 : 10
            hivMInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            artMInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            totMInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            hivFInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            artFInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            totFInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            hivAgeM(j , a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivMInds)) ...
                + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artMInds))) ...
                / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totMInds));
            hivAgeF(j , a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivFInds)) ...
                + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artFInds))) ...
                / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totFInds));
        end
    end
    
    %% HIV prevalence by gender over time vs. AHRI (validation) and (Vandormael, 2019) AHRI data (validation)
    hivYearVec = [unique(hivPrevM_dObs(: ,1)) ; [2010 : 2016]'];
    
    for g = 1 : gender
        for t = 1 : length(hivYearVec)
            prevYear = (hivYearVec(t) - startYear) * stepsPerYear +1;
            hivInds = [toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , g , 4 : 10 , 1 : risk)); toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , g , 4 : 10 , 1 : risk))];
            totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , g , 4 : 10 , 1 : risk));
            hivPrev(j , g , t) = (sum(popVec(prevYear , hivInds) , 2) / sum(popVec(prevYear , totInds) , 2)) * 100;
        end
    end
    
    %% HIV prevalence over time for ages 15-49 vs. AHRI (validation) and HSRC data (validation)
    hivYearVec2 = [2002 , 2005 , 2008 , 2012 , 2017];
    
    for t = 1 : length(hivYearVec2)
        prevYear = (hivYearVec2(t) - startYear) * stepsPerYear +1;
        hivInds = [toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 4 : 10 , 1 : risk)); toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 4 : 10 , 1 : risk))];
        totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 4 : 10 , 1 : risk));
        hivPrevTot(j , t) = (sum(popVec(prevYear , hivInds) , 2) / sum(popVec(prevYear , totInds) , 2)) * 100;
    end
    
    %% HIV prevalence over time for ages 25+ vs. HSRC data (validation)
    hivYearVec2 = [2002 , 2005 , 2008 , 2012 , 2017];
    
    for t = 1 : length(hivYearVec2)
        prevYear = (hivYearVec2(t) - startYear) * stepsPerYear +1;
        hivInds = [toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 6 : age , 1 : risk)); toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 6 : age , 1 : risk))];
        totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 6 : age , 1 : risk));
        hivPrevTot2(j , t) = (sum(popVec(prevYear , hivInds) , 2) / sum(popVec(prevYear , totInds) , 2)) * 100;
    end
    
    %% HIV prevalence over time for ages 50+ vs. HSRC data (validation)
    hivYearVec2 = [2002 , 2005 , 2008 , 2012 , 2017];
    
    for t = 1 : length(hivYearVec2)
        prevYear = (hivYearVec2(t) - startYear) * stepsPerYear +1;
        hivInds = [toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 11 : age , 1 : risk)); toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 11 : age , 1 : risk))];
        totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 11 : age , 1 : risk));
        hivPrevTot3(j , t) = (sum(popVec(prevYear , hivInds) , 2) / sum(popVec(prevYear , totInds) , 2)) * 100;
    end
    
    %% HIV incidence by gender over time for ages 15-49/54 vs. (Vandormael, 2019) AHRI data (validation)
    hivIncYearVec = [2005 : 2017];
    
    hivIncAgeInds = {4:11 , 4:10};
    for g = 1 : gender
        for t = 1 : length(hivIncYearVec)
            incYear = hivIncYearVec(t);
            incTimeSpan = [((incYear - startYear) * stepsPerYear +1) : ((incYear - startYear) * stepsPerYear +6)];
            hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , g , hivIncAgeInds{g} , 1 : risk));
            hivSus = annlz(sum(popVec(incTimeSpan , hivSusInds) , 2)) ./ stepsPerYear;
            hivInc(j , g , t) = annlz(sum(sum(sum(sum(sum(newHiv(incTimeSpan , ...
                : , : , : , g , hivIncAgeInds{g} , 1 : risk), 2), 3), 4), 6), 7)) ./ hivSus * 100;
        end
    end
    
    %% HIV incidence by gender over time for ages 15-29 vs. (Vandormael, 2019) AHRI data (validation)
    hivIncYearVec = [2005 : 2017];
    
    hivIncAgeInds = {4:6 , 4:6};
    for g = 1 : gender
        for t = 1 : length(hivIncYearVec)
            incYear = hivIncYearVec(t);
            incTimeSpan = [((incYear - startYear) * stepsPerYear +1) : ((incYear - startYear) * stepsPerYear +6)];
            hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , g , hivIncAgeInds{g} , 1 : risk));
            hivSus = annlz(sum(popVec(incTimeSpan , hivSusInds) , 2)) ./ stepsPerYear;
            hivIncYoung(j , g , t) = annlz(sum(sum(sum(sum(sum(newHiv(incTimeSpan , ...
                : , : , : , g , hivIncAgeInds{g} , 1 : risk), 2), 3), 4), 6), 7)) ./ hivSus * 100;
        end
    end
        
    %% HIV incidence by gender over time for ages 30-49/54 vs. (Vandormael, 2019) AHRI data (validation)
    hivIncYearVec = [2005 : 2017];
    
    hivIncAgeInds = {7:11 , 7:10};
    for g = 1 : gender
        for t = 1 : length(hivIncYearVec)
            incYear = hivIncYearVec(t);
            incTimeSpan = [((incYear - startYear) * stepsPerYear +1) : ((incYear - startYear) * stepsPerYear +6)];
            hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , g , hivIncAgeInds{g} , 1 : risk));
            hivSus = annlz(sum(popVec(incTimeSpan , hivSusInds) , 2)) ./ stepsPerYear;
            hivIncOlder(j , g , t) = annlz(sum(sum(sum(sum(sum(newHiv(incTimeSpan , ...
                : , : , : , g , hivIncAgeInds{g} , 1 : risk), 2), 3), 4), 6), 7)) ./ hivSus * 100;
        end
    end
    
    %% HIV incidence by gender and 5-year age groups over time vs. IHME model data (validation)
    hivIncYearVec2 = [2000 : 2017];
 
    for t = 1 : length(hivIncYearVec2)
        for a = 4 : 10
            incYear = hivIncYearVec2(t);
            incTimeSpan = [((incYear - startYear) * stepsPerYear +1) : ((incYear - startYear) * stepsPerYear +6)];
            hivSusIndsM = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , 1 , a , 1 : risk));
            hivSusM = annlz(sum(popVec(incTimeSpan , hivSusIndsM) , 2)) ./ stepsPerYear;
            hivIncM(j , a - 3 , t) = annlz(sum(sum(sum(sum(newHiv(incTimeSpan , ...
                : , : , : , 1 , a , 1 : risk), 2), 3), 4), 7)) ./ hivSusM * 100;
            hivSusIndsF = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , 2 , a , 1 : risk));
            hivSusF = annlz(sum(popVec(incTimeSpan , hivSusIndsF) , 2)) ./ stepsPerYear;
            hivIncF(j , a - 3 , t) = annlz(sum(sum(sum(sum(newHiv(incTimeSpan , ...
                : , : , : , 2 , a , 1 : risk), 2), 3), 4), 7)) ./ hivSusF * 100;
        end
    end
    
    %% ********************************** HPV FIGURES **********************************************************************************************

    %% HPV Prevalence by age in 2002 vs. McDonald 2014 data (calibration)
    yr = 2002;
    for a = 4 : 12 % 15-19 -> 55-65
        hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpv_hiv(j , a - 3) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
            ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
        
        hpvInds_hivNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds_hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpv_hivNeg(j , a - 3) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivNeg))...
            ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivNeg));
    end

    %% HPV prevalence by age and HIV status in 2008 vs. Mbulawa data (calibration)
    ageVec = {[4:5],[6:7],[8:9],[10:13]};
    for aV = 1 : length(ageVec)
        a = ageVec{aV};
        hpvInds_hivM = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 , [1 : 2 , 7] , ...
            1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
        ageInds_hivM = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        hpv_hivM2008(j , aV) = sum(popVec((2008 - startYear) * stepsPerYear +1 , hpvInds_hivM))...
            ./ sum(popVec((2008 - startYear) * stepsPerYear+1 , ageInds_hivM));
        
        hpvInds_hivMNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 , [1 : 2 , 7] , ...
            1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
        ageInds_hivMNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        hpv_hivMNeg2008(j , aV) = sum(popVec((2008 - startYear) * stepsPerYear +1 , hpvInds_hivMNeg))...
        ./ sum(popVec((2008 - startYear) * stepsPerYear +1 , ageInds_hivMNeg));
    end

    %% ********************************** CIN FIGURES *********************************************************************************************

    %% CIN2/3 prevalence for All HR HPV types combined by HIV status and age in 2002 vs. McDonald 2014 data (calibration)
    for a = 4 : 13 % 15-19 -> 60-64
        % HIV-positive (on and not on ART)
        cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        cinPos2002(j , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cinInds)))...
            ./ sum(popVec((2002 - startYear) * stepsPerYear +1 , ageInds));
        % HIV-negative
        cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        cinNeg2002(j , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cinNegInds)))...
            ./ (sum(popVec((2002 - startYear) * stepsPerYear +1 , ageNegInds)));
    end
    
    %% CIN1, CIN2, CIN3 prevalence for All HR HPV types combined by HIV status and age in 2002 vs. McDonald 2014 data (calibration)
    for a = 4 : 13 % 15-19 -> 60-64
        % HIV-positive (on and not on ART)
        cin1Inds = unique([toInd(allcomb(3 : 8 , 1 : viral , 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 3 , 7] , 3 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        cin2Inds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 , [1 : 4 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 4 , 7] , 4 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        cin3Inds = unique([toInd(allcomb(3 : 8 , 1 : viral , 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        cin1Pos2002(j , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cin1Inds)))...
            ./ sum(popVec((2002 - startYear) * stepsPerYear +1 , ageInds));
        cin2Pos2002(j , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cin2Inds)))...
            ./ sum(popVec((2002 - startYear) * stepsPerYear +1 , ageInds));
        cin3Pos2002(j , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cin3Inds)))...
            ./ sum(popVec((2002 - startYear) * stepsPerYear +1 , ageInds));
        % HIV-negative
        cin1NegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 3 , 7] , 3 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        cin2NegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 , [1 : 4 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 4 , 7] , 4 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        cin3NegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 5 , 7] , 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        cin1Neg2002(j , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cin1NegInds)))...
            ./ (sum(popVec((2002 - startYear) * stepsPerYear +1 , ageNegInds)));
        cin2Neg2002(j , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cin2NegInds)))...
            ./ (sum(popVec((2002 - startYear) * stepsPerYear +1 , ageNegInds)));
        cin3Neg2002(j , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cin3NegInds)))...
            ./ (sum(popVec((2002 - startYear) * stepsPerYear +1 , ageNegInds)));
    end
    
    %% CIN1, CIN2, CIN3 prevalence for All HR HPV types combined by HIV status in 2015 vs. Kuhn 2015 data (calibration)
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
    cin1Pos2015(j , 1) = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin1Inds)))...
        ./ sum(popVec((2015 - startYear) * stepsPerYear +1 , ageInds));
    cin2Pos2015(j , 1) = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin2Inds)))...
        ./ sum(popVec((2015 - startYear) * stepsPerYear +1 , ageInds));
    cin3Pos2015(j , 1) = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin3Inds)))...
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
    cin1Neg2015(j , 1) = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin1NegInds)))...
        ./ (sum(popVec((2015 - startYear) * stepsPerYear +1 , ageNegInds)));
    cin2Neg2015(j , 1) = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin2NegInds)))...
        ./ (sum(popVec((2015- startYear) * stepsPerYear +1 , ageNegInds)));
    cin3Neg2015(j , 1) = (sum(popVec((2015 - startYear) * stepsPerYear +1 , cin3NegInds)))...
        ./ (sum(popVec((2015 - startYear) * stepsPerYear +1 , ageNegInds)));

    %% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

    %% Cervical cancer incidence in 2012 by age vs. Globocan 2012 data and other sources (calibration)
    incTimeSpan = [((2012 - startYear) * stepsPerYear +1) : ((2012 - startYear) * stepsPerYear +6)];
    fac = 10 ^ 5;

    for a = 1 : age
        % General population
        allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012(j , a) = ...
            (annlz(sum(sum(sum(newCC(incTimeSpan , : , a , :),2),3),4)) ./ ...
            (annlz(sum(popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);
        
        % HIV-negative
        allFneg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012neg(j , a) = ...
            (annlz(sum(sum(sum(newCC(incTimeSpan , 1 : 2 , a , :),2),3),4)) ./ ...
            (annlz(sum(popVec(incTimeSpan , allFneg) , 2) ./ stepsPerYear)) * fac);
        
        % HIV-positive untreated
        allFpos = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012pos(j , a) = ...
            (annlz(sum(sum(sum(newCC(incTimeSpan , 3 : 7 , a , :),2),3),4)) ./ ...
            (annlz(sum(popVec(incTimeSpan , allFpos) , 2) ./ stepsPerYear)) * fac);
        
        % HIV-positive on ART
        allFart = toInd(allcomb(8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012art(j , a) = ...
            (annlz(sum(sum(sum(newCC(incTimeSpan , 8 , a , :),2),3),4)) ./ ...
            (annlz(sum(popVec(incTimeSpan , allFart) , 2) ./ stepsPerYear)) * fac);
    end
   
    %% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************
    
    %% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)
    typeDistYearVec = [2010 , 2011 , 2012 , 2013 , 2014 , 2015];
    ccInds_vax = toInd(allcomb(1 : disease , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
        1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
    ccInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , 6 , ...
        1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
    ccInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
            1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
        
    cin3Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin3Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 4 , 7] , 5 , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin3Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
        
    cin1Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin1Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 2 , 7] , 3 , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin1Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
   
    hpvInds_vax = toInd(allcomb(1 : disease , 1 : viral , 2 , [1 : 2 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    hpvInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 , 7] , 2 , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    hpvInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 2 , [1 : 2 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 , 7] , 2 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

    for i = 1 : length(typeDistYearVec)
        cc_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);
        cc_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);

        cin3_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);
        cin3_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);

        cin1_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_tot) , 2);
        cin1_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_tot) , 2);

        hpv_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_tot) , 2);
        hpv_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_tot) , 2);
    end
    
end

%% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

%% Population size over time vs. Statistics South Africa data (calibration)
popYearVec = unique(totPopSize_dObs(: ,1));
% Load calibration data from Excel (years, values)
file = [pwd , '/Config/Population_validation_targets.xlsx'];
historicalPop0_69 = zeros(5,2);
futurePop0_69 = zeros(16,2);
historicalPop0_69(:,1) = xlsread(file , 'Demographics' , 'B91:F91'); % years
historicalPop0_69(:,2) = xlsread(file , 'Demographics' , 'B130:F130') .* 1000; % estimates
futurePop0_69(:,1) = xlsread(file , 'Demographics' , 'C144:R144'); % years
futurePop0_69(:,2) = xlsread(file , 'Demographics' , 'C146:R146') .* 1000; % projections

% Calibration error bars
meanObs = totPopSize_dObs(: , 2);
sdevObs = (totPopSize_dObs(: , 3).^(1/2)).*2;

figure;
errorbar(totPopSize_dObs(: , 1) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
%boxplot(popSize , 'Positions' , popYearVec , 'Labels' , popYearVec , 'Color' , 'k' , 'Whisker' , 5);
plot(totPopSize_dObs(: , 1) , mean(popSize,1)' , 'k-' , ...
    totPopSize_dObs(: , 1) , min(popSize,[],1)' , 'k--' , ...
    totPopSize_dObs(: , 1) , max(popSize,[],1)' , 'k--' , 'LineWidth' , 1.5);
title('KZN Population Size Ages 0-79')
xlabel('Year'); ylabel('Individuals')
xlim([2000 2020]); ylim([0 (14*10^6)]);
legend('(Statistics SA) Observed KZN, ages 0-79: mean, 2SD' , 'Model, ages 0-79: 25-sets mean' , 'Model, ages 0-79: 25-sets minimum' , 'Model, ages 0-79: 25-sets maximum');
grid on;

%% Population size by 5-year age groups over time vs. Statistics South Africa data (calibration)
% Load calibration data from Excel
file = [pwd , '/Config/Population_validation_targets.xlsx'];
years = xlsread(file , 'Demographics' , 'B91:F91');    % years
kzn_popByage_yrs(: , :) = xlsread(file , 'Demographics' , 'M92:Q107').*1000;    % males and females by age in 1996-2019
popProp_obs = zeros(5,age);
for y = 1 : length(years)
    yearCurr = years(y);
    for a = 1 : age
        popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
        popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
        popProp_obs(y,a) = sum(kzn_popByage_yrs(a , y)) / sumall(kzn_popByage_yrs(1 : end , y));
    end
end

popYearVec = unique(popAgeDist_dObs(: ,1));

% Calibration error bars
meanObs = [popAgeDist_dObs(1:16 , 2) , popAgeDist_dObs(17:32 , 2) , popAgeDist_dObs(33:48 , 2)]';
sdevObs = ([popAgeDist_dObs(1:16 , 3) , popAgeDist_dObs(17:32 , 3) , popAgeDist_dObs(33:48 , 3)]'.^(1/2)).*2;

figure;
subplot(1,3,1);
set(gca,'ColorOrderIndex',1)
calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , ...
    unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) ,];
errorbar(calibYrs , meanObs(: , 1:7) , sdevObs(: , 1:7) , ...
    's' , 'LineWidth' , 1.5); %'Color' , [0.9290, 0.6940, 0.1250])
hold on;
% set(gca,'ColorOrderIndex',1)
% plot(years , popProp_obs(: , 1:7) , 'o');
% hold on;
for a = 1 : 7
    set(gca,'ColorOrderIndex',a)
    plot(popYearVec , mean(popProp(: , : , a),1)' , '-' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a)
    plot(popYearVec , min(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a)
    plot(popYearVec , max(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
end
ylim([0.05 0.18]);
ylabel('Population proportion by age'); xlabel('Year');
legend('(Statistics SA) Observed KZN, ages 0-4: mean, 2SD' , 'ages 5-9: mean, 2SD' , 'ages 10-14: mean, 2SD' , 'ages 15-19: mean, 2SD' , 'ages 20-24: mean, 2SD' , 'ages 25-29: mean, 2SD' , 'ages 30-34: mean, 2SD' , ...
    'Model, ages 0-4: 25-sets mean' , 'ages 0-4: 25-sets minimum' , 'ages 0-4: 25-sets maximum' , '...' , 'Location' , 'north');

subplot(1,3,2);
set(gca,'ColorOrderIndex',1)
calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , ...
    unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) ,];
errorbar(calibYrs , meanObs(: , 8:14) , sdevObs(: , 8:14) , ...
    's' , 'LineWidth' , 1.5); %'Color' , [0.9290, 0.6940, 0.1250])
hold on;
% set(gca,'ColorOrderIndex',1)
% plot(years , popProp_obs(: , 8:14) , 'o');
% hold on;
for a = 8 : 14
    set(gca,'ColorOrderIndex',a-7)
    plot(popYearVec , mean(popProp(: , : , a),1)' , '-' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-7)
    plot(popYearVec , min(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-7)
    plot(popYearVec , max(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
end
ylim([0.0 0.13]);
ylabel('Population proportion by age'); xlabel('Year');
legend('(Statistics SA) Observed KZN, ages 35-39: mean, 2SD' , 'ages 40-44: mean, 2SD' , 'ages 45-49: mean, 2SD' , 'ages 50-54: mean, 2SD' , 'ages 55-59: mean, 2SD' , 'ages 60-64: mean, 2SD' , 'ages 65-69: mean, 2SD' , ...
    'Model, ages 35-39: 25-sets mean' , 'ages 35-39: 25-sets minimum' , 'ages 35-39: 25-sets maximum' , '...' , 'Location' , 'north');


subplot(1,3,3);
set(gca,'ColorOrderIndex',1)
calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1))];
errorbar(calibYrs , meanObs(: , 15:16) , sdevObs(: , 15:16) , ...
    's' , 'LineWidth' , 1.5); %'Color' , [0.9290, 0.6940, 0.1250])
hold on;
% set(gca,'ColorOrderIndex',1)
% plot(years , popProp_obs(: , 15:16) , 'o');
% hold on;
for a = 15 : 16
    set(gca,'ColorOrderIndex',a-14)
    plot(popYearVec , mean(popProp(: , : , a),1)' , '-' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-14)
    plot(popYearVec , min(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-14)
    plot(popYearVec , max(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
end
ylim([0.0 0.04]);
ylabel('Population proportion by age'); xlabel('Year'); %title('KZN age distribution in 5-year groups');
legend('(Statistics SA) Observed KZN, ages 70-74: mean, 2SD' , 'ages 75-79: mean, 2SD' , ...
    'Model, ages 70-74: 25-sets mean' , 'ages 70-74: 25-sets minimum' , 'ages 70-74: 25-sets maximum' , '...' , 'Location' , 'north');

%% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

%% HIV prevalence by age over time vs. AHRI (calibration, validation) and IHME model data (validation)
hivYearVec = [unique(hivPrevM_dObs(: ,1)) ; [2010 : 2016]'];
prevYears2 = [2010 : 2016];

% 2010-2016 AC data
hivPrevF_val = [9.29	9.02	10.45	9.33	10.37	11.00	9.35
    31.41	31.68	30.64	33.95	34.56	34.12	33.42
    53.27	51.72	50.80	51.33	51.94	53.98	52.41
    59.18	61.35	58.66	64.90	62.57	64.71	63.09
    53.97	54.08	58.77	65.12	65.28	64.66	66.95
    42.69	43.27	45.29	49.16	54.25	56.37	61.28
    32.34	34.30	39.18	41.47	48.21	49.57	50.23];
hivPrevM_val = [1.60	1.85	2.75	3.46	2.87	3.95	4.50
    9.56	8.02	9.87	9.65	11.86	7.19	8.02
    28.99	21.92	24.88	29.84	35.40	27.65	27.31
    46.47	44.51	39.49	47.22	46.35	41.64	42.08
    52.03	44.30	49.61	63.33	51.41	52.05	51.35
    41.73	41.53	51.55	51.64	59.40	52.69	51.18
    36.64	37.12	33.01	40.00	40.54	44.52	52.17];

% 2012 & 2017 HSRC SABSSM data for SA, ages 15-59 in 5-year age groups
hivPrevF_val2 = [5.6 17.4 28.4 36.0 31.6 28.0 19.7 14.8 9.7 ... % 2012
                 5.8 15.6 27.5 34.7 39.4 35.9 30.3 22.2 17.6]; % 2017
hivPrevM_val2 = [0.7 5.1 17.3 25.6 28.8 15.8 13.4 15.5 5.5 ... % 2012
                 4.7 4.8 12.4 18.4 23.7 22.4 24.8 20.2 14.8]; % 2017
             
% IHME model prevalence estimates for KZN, ages 15-49 in 5-year age groups
file = [pwd , '/Config/ihme_hiv_kzn_validation.xlsx'];
hivPrevF_val3 = xlsread(file , 'Prevalence' , 'L128:N253').*100; % 2000-2017 (val, upper, lower)
hivPrevM_val3 = xlsread(file , 'Prevalence' , 'L2:N127').*100; % 2000-2017 (val, upper, lower)

% Calibration error bars
hivM(: , 1) = hivPrevM_dObs(: , 2) .* 100; % mean
hivM(: , 2) = (hivPrevM_dObs(: , 3).^(1/2)).*2 .* 100; % calibration SD
hivF(: , 1) = hivPrevF_dObs(: , 2) .* 100; % mean
hivF(: , 2) = (hivPrevF_dObs(: , 3).^(1/2)).*2 .* 100; % calibration SD

ageGroup = {'15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74'};

gen = {'Male' , 'Female'};
for g = 1 : gender
    hivPrevs = hivM;
    hivPrevs2 = hivPrevM_val;
    hivPrevs3 = hivPrevM_val2;
    hivPrevs4 = hivPrevM_val3;
    hivModel = hivAgeM;
    if g == 2
        hivPrevs = hivF;
        hivPrevs2 = hivPrevF_val;
        hivPrevs3 = hivPrevF_val2;
        hivPrevs4 = hivPrevF_val3;
        hivModel = hivAgeF;
    end

    figure; 
    for a = 4 : 10
        subplot(3 , 3 , a-3)
        hold all;
        if a <= 11            
            errorbar(unique(hivPrevM_dObs(: ,1)) , hivPrevs(((a-3) : 7 : end) , 1) , hivPrevs(((a-3) : 7 : end) , 2) , ...
                'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
            hold on;
            plot(prevYears2 , hivPrevs2((a-3) : 7 : end) , 'bo');
            hold on;
            %plot([2012 2017] , hivPrevs3((a-3) : 9 : end) , 'co');
            errorbar([2000 : 2017] , hivPrevs4(((a-3) : 7 : end) , 1)' , ...
                hivPrevs4(((a-3) : 7 : end) , 1)' - hivPrevs4(((a-3) : 7 : end) , 3)' , ...
                hivPrevs4(((a-3) : 7 : end) , 2)' - hivPrevs4(((a-3) : 7 : end) , 1)' , ...
                'cs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
            hold on;
            %boxplot((squeeze(hivModel(: , a-3 , :)) .* 100) , ...
            %    'Positions' , hivYearVec , 'Labels' , hivYearVec , 'Color' , 'k' , 'Whisker' , 5);
            plot(hivYearVec , mean((squeeze(hivModel(: , a-3 , :)) .* 100),1)' , 'k-' , ...
                hivYearVec , min((squeeze(hivModel(: , a-3 , :)) .* 100),[],1)' , 'k--' , ...
                hivYearVec , max((squeeze(hivModel(: , a-3 , :)) .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
        else
            %plot([2012 2017] , hivPrevs3((a-3) : 9 : end) , 'co');
            errorbar([2000 : 2017] , hivPrevs4(((a-3) : 7 : end) , 1)' , ...
                hivPrevs4(((a-3) : 7 : end) , 1)' - hivPrevs4(((a-3) : 7 : end) , 3)' , ...
                hivPrevs4(((a-3) : 7 : end) , 2)' - hivPrevs4(((a-3) : 7 : end) , 1)' , ...
                'cs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
            hold on;
            %boxplot((squeeze(hivModel(: , a-3 , :)) .* 100) , ...
            %    'Positions' , hivYearVec , 'Labels' , hivYearVec , 'Color' , 'k' , 'Whisker' , 5);
            plot(hivYearVec , mean((squeeze(hivModel(: , a-3 , :)) .* 100),1)' , 'k-' , ...
                hivYearVec , min((squeeze(hivModel(: , a-3 , :)) .* 100),[],1)' , 'k--' , ...
                hivYearVec , max((squeeze(hivModel(: , a-3 , :)) .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
        end
        xlabel('Year'); ylabel('HIV Prevalence (%)'); title([gen{g} , 's ages ' , ageGroup{a-3}]) % , ' HIV Prevalence'])
        xlim([2000 2020]); 
        if g == 1 
            ylim([0 60]);
        elseif g == 2
            ylim([0 80]);
        end
        grid on;
    end
    legend('(AHRI data request) Observed KZN: mean, 2SD' , '(AHRI data request) Observed KZN' , '(IHME model) Observed KZN: val, lower lb, upper lb' , ...
        'Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum')
end

%% HIV prevalence by gender over time vs. AHRI (validation) and (Vandormael, 2019) AHRI data (validation)
hivYearVec = [unique(hivPrevM_dObs(: ,1)) ; [2010 : 2016]'];

hivData(: , : , 1) = zeros(length(unique(hivPrevM_dObs(: ,1))) , 1);
hivData(: , : , 2) = zeros(length(unique(hivPrevM_dObs(: ,1))) , 1);

hivRaw(:,:,1) = hivPrevM_dObs(: , 4:5);
hivRaw(:,:,2) = hivPrevF_dObs(: , 4:5);

for i = 1 : length(unique(hivPrevM_dObs(: ,1)))
    for g = 1 : gender
        hivData(i,1,g) = (sumall(hivRaw(((i-1)*7+1):(i*7) , 1 , g)) ./ sumall(hivRaw(((i-1)*7+1):(i*7) , 2 , g))) .* 100;
    end
end

hivYearVecAlt = [2005 : 2017]';
hivDataAlt(: , 1 , 1) = [13.83; 13.30; 14.78; 14.52; 17.30; 16.61; 16.87; 17.48; 19.13; 20.87; 18.96; 19.91; 20.17]; % AHRI KZN, males, ages 15-54: (Vandormael, 2019)
hivDataAlt(: , 1 , 2) = [25.49; 25.21; 27.41; 28.90; 31.72; 34.71; 35.14; 36.19; 38.12; 40.76; 41.81; 42.06; 40.45]; % AHRI KZN, females , ages 15-49: (Vandormael, 2019)

figure;
gen = {'Male' , 'Female'};
for g = 1 : gender
    subplot(1,2,g)
    plot(unique(hivPrevM_dObs(: ,1)) , hivData(:,:,g) , 'bo');
    hold on;
    plot(hivYearVecAlt , hivDataAlt(:,:,g) , 'co');
    hold on;
    plot(hivYearVec , mean(squeeze(hivPrev(: , g , :)),1)' , 'k-' , ...
        hivYearVec , min(squeeze(hivPrev(: , g , :)),[],1)' , 'k--' , ...
        hivYearVec , max(squeeze(hivPrev(: , g , :)),[],1)' , 'k--' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
    xlabel('Year'); ylabel('HIV Prevalence (%)'); title(gen{g});
    xlim([2000 2020]); ylim([0 50]);
    if g == 1 
        legend('(AHRI data request) Observed KZN, ages 15-49' , '(Vandormael, 2019) Observed KZN, ages 15-54' , ...
            'Model, ages 15-49: 25-sets mean' , 'Model, ages 15-49: 25-sets minimum' , 'Model, ages 15-49: 25-sets maximum');
    elseif g == 2
        legend('(AHRI data request) Observed KZN, ages 15-49' , '(Vandormael, 2019) Observed KZN, ages 15-49' , ...
            'Model, ages 15-49: 25-sets mean' , 'Model, ages 15-49: 25-sets minimum' , 'Model, ages 15-49: 25-sets maximum');
    end
end

%% HIV prevalence over time for ages 15-49 vs. AHRI (validation) and HSRC data (validation)
hivYearVec2 = [2002 , 2005 , 2008 , 2012 , 2017];
hivPrev_obs = [2002 15.7 11.6 21.1
               2005 21.9 18.3 25.9
               2008 25.8 22.1 29.8
               2012 27.9 25.2 30.8
               2017 27.0 23.9 30.4];

prevYears2 = [2010 : 2016];
hivPrev_val = [29.04 29.05 29.32 31.19 33.81 33.82 34.45];

figure;
errorbar(hivPrev_obs(: , 1) , hivPrev_obs(: , 2) , ...
        hivPrev_obs(: , 2) - hivPrev_obs(: , 3) , hivPrev_obs(: , 4) - hivPrev_obs(: , 2) , ...
        'bs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
plot(prevYears2 , hivPrev_val , 'co');
hold on;
plot(hivYearVec2 , mean(hivPrevTot(: , :),1)' , 'k-' , ...
    hivYearVec2 , min(hivPrevTot(: , :),[],1)' , 'k--' , ...
    hivYearVec2 , max(hivPrevTot(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
xlabel('Year'); ylabel('HIV Prevalence (%)');
xlim([2000 2020]); ylim([0 50]);
legend('(HSRC SABSSM) Observed KZN, ages 15-49: 95% CI' , '(AHRI data request) Observed KZN, ages 15-49' , ...
    'Model, ages 15-49: 25-sets mean' , 'Model, ages 15-49: 25-sets minimum' , 'Model, ages 15-49: 25-sets maximum');

%% HIV prevalence over time for ages 25+ vs. HSRC data (validation)
hivYearVec2 = [2002 , 2005 , 2008 , 2012 , 2017];

hivPrev_obs = [2002 14.9 10.1 21.5
               2005 20.5 16.8 24.6
               2008 23.5 19.7 27.8
               2012 30.1 26.9 33.6
               2017 31.2 27.6 35.0];

figure;
errorbar(hivPrev_obs(: , 1) , hivPrev_obs(: , 2) , ...
        hivPrev_obs(: , 2) - hivPrev_obs(: , 3) , hivPrev_obs(: , 4) - hivPrev_obs(: , 2) , ...
        'bs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
plot(hivYearVec2 , mean(hivPrevTot2(: , :),1)' , 'k-' , ...
    hivYearVec2 , min(hivPrevTot2(: , :),[],1)' , 'k--' , ...
    hivYearVec2 , max(hivPrevTot2(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
xlabel('Year'); ylabel('HIV Prevalence (%)');
xlim([2000 2020]); ylim([0 50]);
legend('(HSRC SABSSM) Observed KZN, ages 25+: 95% CI' , 'Model, ages 25-79: 25-sets mean' , ...
    'Model, ages 25-79: 25-sets minimum' , 'Model, ages 25-79: 25-sets maximum');

%% HIV prevalence over time for ages 50+ vs. HSRC data (validation)
hivYearVec2 = [2002 , 2005 , 2008 , 2012 , 2017];

hivPrev_obs = [2002 11.0 4.5 24.3
               2005 9.5 5.9 14.8
               2008 6.1 3.7 10.1
               2012 9.8 7.4 12.8
               2017 17.9 13.9 22.8];

figure;
errorbar(hivPrev_obs(: , 1) , hivPrev_obs(: , 2) , ...
        hivPrev_obs(: , 2) - hivPrev_obs(: , 3) , hivPrev_obs(: , 4) - hivPrev_obs(: , 2) , ...
        'bs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
plot(hivYearVec2 , mean(hivPrevTot3(: , :),1)' , 'k-' , ...
    hivYearVec2 , min(hivPrevTot3(: , :),[],1)' , 'k--' , ...
    hivYearVec2 , max(hivPrevTot3(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
xlabel('Year'); ylabel('HIV Prevalence (%)');
xlim([2000 2020]); ylim([0 50]);
legend('(HSRC SABSSM) Observed KZN, ages 50+: 95% CI' , 'Model, ages 50-79: 25-sets mean' , ...
    'Model, ages 50-79: 25-sets minimum' , 'Model, ages 50-79: 25-sets maximum');

%% HIV incidence by gender and broad age groups over time vs. (Vandormael, 2019) AHRI data (validation)
hivIncYearVec = [2005 : 2017];

hivInc_obs(: , : , 1) = [2005 2.14 1.57 2.93; % AHRI KZN, males, ages 15-54: (Vandormael, 2019)
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
hivInc_obs(: , : , 2) = [2005 4.08 3.40 4.90; % AHRI KZN, females , ages 15-49: (Vandormael, 2019)
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
                     
% Import from Excel incidence validation targets by gender and broad age groups over time
file = [pwd , '/Config/AHRI_KZN_HIV_incidence-byGenderAge-validation.xlsx'];
hivIncYoung_obs(: , : , 1) = xlsread(file , 'HIV_incidence' , 'A21:D33'); % AHRI KZN, males, ages 15-29: (Vandormael, 2019)
hivIncYoung_obs(: , : , 2) = xlsread(file , 'HIV_incidence' , 'A5:D17'); % AHRI KZN, females, ages 15-29: (Vandormael, 2019)           
hivIncOlder_obs(: , : , 1) = xlsread(file , 'HIV_incidence' , 'F21:I33'); % AHRI KZN, males, ages 30-54: (Vandormael, 2019)   
hivIncOlder_obs(: , : , 2) = xlsread(file , 'HIV_incidence' , 'F5:I17'); % AHRI KZN, females, ages 30-49: (Vandormael, 2019)

figure;
titlesHIVinc = {'Male: ages 15-54' , 'Male: ages 15-29' , 'Male: ages 30-54' , 'Female: ages 15-49' , 'Female: ages 15-29' , 'Female: ages 30-49'};
for g = 1 : gender
    subplot(2,3,(1 + ((g-1)*3)))
    errorbar(hivInc_obs(: , 1 , g) , hivInc_obs(: , 2 , g) , ...
        hivInc_obs(: , 2 , g) - hivInc_obs(: , 3 , g) , hivInc_obs(: , 4 , g) - hivInc_obs(: , 2 , g) , ...
        'bs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
    hold on;
    plot(hivIncYearVec , mean(squeeze(hivInc(: , g , :)),1)' , 'k-' , ...
        hivIncYearVec , min(squeeze(hivInc(: , g , :)),[],1)' , 'k--' , ...
        hivIncYearVec , max(squeeze(hivInc(: , g , :)),[],1)' , 'k--' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
    xlabel('Year'); ylabel('HIV incidence per 100'); 
    xlim([2000 2020]); ylim([0 10]); title(titlesHIVinc{(1 + ((g-1)*3))});
    legend('(Vandormael, 2019) Observed KZN: 95% CI' , 'Model: 25-sets mean' , ...
        'Model: 25-sets minimum' , 'Model: 25-sets maximum');
    
    subplot(2,3,(2 + ((g-1)*3)))
    errorbar(hivIncYoung_obs(: , 1 , g) , hivIncYoung_obs(: , 2 , g) , ...
        hivIncYoung_obs(: , 2 , g) - hivIncYoung_obs(: , 3 , g) , hivIncYoung_obs(: , 4 , g) - hivIncYoung_obs(: , 2 , g) , ...
        'bs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
    hold on;
    plot(hivIncYearVec , mean(squeeze(hivIncYoung(: , g , :)),1)' , 'k-' , ...
        hivIncYearVec , min(squeeze(hivIncYoung(: , g , :)),[],1)' , 'k--' , ...
        hivIncYearVec , max(squeeze(hivIncYoung(: , g , :)),[],1)' , 'k--' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
    xlabel('Year'); ylabel('HIV incidence per 100'); 
    xlim([2000 2020]); ylim([0 10]); title(titlesHIVinc{(2 + ((g-1)*3))});
    legend('(Vandormael, 2019) Observed KZN: 95% CI' , 'Model: 25-sets mean' , ...
        'Model: 25-sets minimum' , 'Model: 25-sets maximum');
    
    subplot(2,3,(3 + ((g-1)*3)))
    errorbar(hivIncOlder_obs(: , 1 , g) , hivIncOlder_obs(: , 2 , g) , ...
        hivIncOlder_obs(: , 2 , g) - hivIncOlder_obs(: , 3 , g) , hivIncOlder_obs(: , 4 , g) - hivIncOlder_obs(: , 2 , g) , ...
        'bs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
    hold on;
    plot(hivIncYearVec , mean(squeeze(hivIncOlder(: , g , :)),1)' , 'k-' , ...
        hivIncYearVec , min(squeeze(hivIncOlder(: , g , :)),[],1)' , 'k--' , ...
        hivIncYearVec , max(squeeze(hivIncOlder(: , g , :)),[],1)' , 'k--' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
    xlabel('Year'); ylabel('HIV incidence per 100'); 
    xlim([2000 2020]); ylim([0 10]); title(titlesHIVinc{(3 + ((g-1)*3))});
    legend('(Vandormael, 2019) Observed KZN: 95% CI' , 'Model: 25-sets mean' , ...
        'Model: 25-sets minimum' , 'Model: 25-sets maximum');
end

%% HIV incidence by gender and 5-year age groups over time vs. IHME model data (validation)
hivIncYearVec2 = [2000 : 2017];

% IHME model prevalence estimates for KZN, ages 15-49 in 5-year age groups
file = [pwd , '/Config/ihme_hiv_kzn_validation.xlsx'];
hivIncF_val = xlsread(file , 'Incidence' , 'L128:N253').*100; % 2000-2017 (val, upper, lower)
hivIncM_val = xlsread(file , 'Incidence' , 'L2:N127').*100; % 2000-2017 (val, upper, lower)

ageGroup = {'15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74'};

gen = {'Male' , 'Female'};
for g = 1 : gender
    hivIncR = hivIncM_val;
    hivModel = hivIncM;
    if g == 2
        hivIncR = hivIncF_val;
        hivModel = hivIncF;
    end

    figure; 
    for a = 4 : 10
        subplot(3 , 3 , a-3)
        hold all;
        if a <= 11            
            errorbar(hivIncYearVec2 , hivIncR(((a-3) : 7 : end) , 1)' , ...
                hivIncR(((a-3) : 7 : end) , 1)' - hivIncR(((a-3) : 7 : end) , 3)' , ...
                hivIncR(((a-3) : 7 : end) , 2)' - hivIncR(((a-3) : 7 : end) , 1)' , ...
                'cs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
            hold on;
            plot(hivIncYearVec2 , mean(squeeze(hivModel(: , a-3 , :)),1)' , 'k-' , ...
                hivIncYearVec2 , min(squeeze(hivModel(: , a-3 , :)),[],1)' , 'k--' , ...
                hivIncYearVec2 , max(squeeze(hivModel(: , a-3 , :)),[],1)' , 'k--' , 'LineWidth' , 1.5);
        else
            errorbar(hivIncYearVec2 , hivIncR(((a-3) : 7 : end) , 1)' , ...
                hivIncR(((a-3) : 7 : end) , 1)' - hivIncR(((a-3) : 7 : end) , 3)' , ...
                hivIncR(((a-3) : 7 : end) , 2)' - hivIncR(((a-3) : 7 : end) , 1)' , ...
                'cs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
            hold on;
            plot(hivIncYearVec2 , mean(squeeze(hivModel(: , a-3 , :)),1)' , 'k-' , ...
                hivIncYearVec2 , min(squeeze(hivModel(: , a-3 , :)),[],1)' , 'k--' , ...
                hivIncYearVec2 , max(squeeze(hivModel(: , a-3 , :)),[],1)' , 'k--' , 'LineWidth' , 1.5);
        end
        xlabel('Year'); ylabel('HIV incidence per 100'); title([gen{g} , 's ages ' , ageGroup{a-3}])
        xlim([2000 2020]); ylim([0 10]);
        grid on;
    end
    legend('(IHME model) Observed KZN: val, lower lb, upper lb' , ...
        'Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum')
end

%% ********************************** HPV FIGURES **********************************************************************************************

%% HPV Prevalence by age in 2002 vs. McDonald 2014 data (calibration)
ageGroup = {'17 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};

% Calibration error bars
meanObs = hpv_hiv_dObs(: , 2) .* 100;
sdevObs = (hpv_hiv_dObs(: , 3).^(1/2)).*2 .* 100;
meanNeg = hpv_hivNeg_dObs(: , 2) .* 100;
sdevNeg = (hpv_hivNeg_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
subplot(2,1,1);
errorbar(1 : length(meanObs) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
%boxplot((hpv_hiv .* 100) , 'Color' , 'k' , 'Whisker' , 5);
plot(1 : length(meanObs) , mean((hpv_hiv .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_hiv .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_hiv .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('(McDonald, 2014) Observed Cape Town: mean, 2SD' , 'Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
title('hrHPV Prevalence in 2002 - Females, HIV+');

subplot(2,1,2);
errorbar(1 : length(meanObs) , meanNeg , sdevNeg , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
%boxplot((hpv_hivNeg .* 100) , 'Color' , 'k' , 'Whisker' , 5);
plot(1 : length(meanObs) , mean((hpv_hivNeg .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_hivNeg .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_hivNeg .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('(McDonald, 2014) Observed Cape Town: mean, 2SD' , 'Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
title('hrHPV Prevalence in 2002 - Females, HIV-');
grid on;

%% HPV prevalence by age and HIV status in 2008 vs. Mbulawa data (calibration)
ageGroup = {'15-24' , '25-34' , '35-44' , '45-64'};

% Calibration error bars
meanObs = hpv_hivM2008_dObs(: , 2) .* 100;
sdevObs = (hpv_hivM2008_dObs(: , 3).^(1/2)).*2 .* 100;
meanNeg = hpv_hivMNeg2008_dObs(: , 2) .* 100;
sdevNeg = (hpv_hivMNeg2008_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
subplot(2,1,1)
errorbar(1 : length(meanObs) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
%boxplot((hpv_hivM2008 .* 100) , 'Color' , 'k' , 'Whisker' , 5)
plot(1 : length(meanObs) , mean((hpv_hivM2008 .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_hivM2008 .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_hivM2008 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtick' , [1 : length(ageGroup)] , 'xtickLabel' , ageGroup);
legend('(Mbulawa, 2015) Observed SA: mean, 2SD' , 'Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)'); ylim([0 100]);
title('hrHPV Prevalence in 2008 - Males, HIV+');
grid on;

subplot(2,1,2)
errorbar(1 : length(meanObs) , meanNeg , sdevNeg , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
%boxplot((hpv_hivMNeg2008 .* 100) , 'Color' , 'k' , 'Whisker' , 5)
plot(1 : length(meanObs) , mean((hpv_hivMNeg2008 .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_hivMNeg2008 .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_hivMNeg2008 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtick' , [1 : length(ageGroup)] , 'xtickLabel' , ageGroup);
legend('(Mbulawa, 2015) Observed SA: mean, 2SD' , 'Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)'); ylim([0 100]);
title('hrHPV Prevalence in 2008 - Males, HIV-');
grid on;

%% ********************************** CIN FIGURES *********************************************************************************************

%% CIN2/3 prevalence for All HR HPV types combined by HIV status and age in 2002 vs. McDonald 2014 data (calibration)
ageGroup = {'17-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

% Calibration error bars
cinPos_mean = cinPos2002_dObs(: , 2) .* 100;
cinPos_sdev = (cinPos2002_dObs(: , 3).^(1/2)).*2 .* 100;
cinNeg_mean = cinNeg2002_dObs(: , 2) .* 100;
cinNeg_sdev = (cinNeg2002_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
subplot(2 , 1 , 1);
errorbar(1 : length(cinPos_mean) , cinPos_mean , cinPos_sdev , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
%boxplot((cinPos2002 .* 100) , 'Color' , 'k' , 'Whisker' , 5)
plot(1 : length(cinPos_mean) , mean((cinPos2002 .* 100),1)' , 'k-' , ...
    1 : length(cinPos_mean) , min((cinPos2002 .* 100),[],1)' , 'k--' , ...
    1 : length(cinPos_mean) , max((cinPos2002 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
legend('(McDonald, 2014) Observed Cape Town: mean, 2SD' , 'Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence (%)')
title('CIN 2/3 Prevalence in 2002 - Females, HIV+')
ylim([0 25])
grid on;

subplot(2 , 1 , 2)
errorbar(1 : length(cinNeg_mean) , cinNeg_mean , cinNeg_sdev , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
%boxplot((cinNeg2002 .* 100) , 'Color' , 'k' , 'Whisker' , 5);
plot(1 : length(cinNeg_mean) , mean((cinNeg2002 .* 100),1)' , 'k-' , ...
    1 : length(cinNeg_mean) , min((cinNeg2002 .* 100),[],1)' , 'k--' , ...
    1 : length(cinNeg_mean) , max((cinNeg2002 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
legend('(McDonald, 2014) Observed Cape Town: mean, 2SD' , 'Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence (%)')
title('CIN 2/3 Prevalence in 2002 - Females, HIV-')
ylim([0 10])
grid on;

%% CIN1, CIN2, CIN3 prevalence for All HR HPV types combined by HIV status and age in 2002 vs. McDonald 2014 data (calibration)
ageGroup = {'17-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

% Calibration error bars
cinPos_mean = cinPos2002_dObs(: , 2) .* 100;
cinPos_sdev = (cinPos2002_dObs(: , 3).^(1/2)).*2 .* 100;
cinNeg_mean = cinNeg2002_dObs(: , 2) .* 100;
cinNeg_sdev = (cinNeg2002_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
subplot(2 , 1 , 1);
errorbar(1 : length(cinPos_mean) , cinPos_mean , cinPos_sdev , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
%boxplot((cinPos2002 .* 100) , 'Color' , 'k' , 'Whisker' , 5)
plot(1 : length(cinPos_mean) , mean((cin1Pos2002 .* 100),1)' , 'g-' , ...
    1 : length(cinPos_mean) , min((cin1Pos2002 .* 100),[],1)' , 'g--' , ...
    1 : length(cinPos_mean) , max((cin1Pos2002 .* 100),[],1)' , 'g--' , 'LineWidth' , 1.5);
plot(1 : length(cinPos_mean) , mean((cin2Pos2002 .* 100),1)' , 'c-' , ...
    1 : length(cinPos_mean) , min((cin2Pos2002 .* 100),[],1)' , 'c--' , ...
    1 : length(cinPos_mean) , max((cin2Pos2002 .* 100),[],1)' , 'c--' , 'LineWidth' , 1.5);
plot(1 : length(cinPos_mean) , mean((cin3Pos2002 .* 100),1)' , 'm-' , ...
    1 : length(cinPos_mean) , min((cin3Pos2002 .* 100),[],1)' , 'm--' , ...
    1 : length(cinPos_mean) , max((cin3Pos2002 .* 100),[],1)' , 'm--' , 'LineWidth' , 1.5);
hold all;
formatTemp = ones(1 , length(1 : length(cinPos_mean)));
plot(1 : length(cinPos_mean) , formatTemp.*14 , ':g'); hold all;
plot(1 : length(cinPos_mean) , formatTemp.*8 , ':b'); hold all;
plot(1 : length(cinPos_mean) , formatTemp.*7 , ':y');
legend('(McDonald, 2014) Observed CIN2/3 Cape Town: mean, 2SD' , 'Model, CIN1: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, CIN2: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, CIN3: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    '(Kuhn, 2020) CIN1 Cape Town' , '(Kuhn, 2020) CIN2 Cape Town' , '(Kuhn 2020) CIN3 Cape Town');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence (%)')
title('CIN 2/3 Prevalence in 2002 - Females, HIV+')
ylim([0 25])
grid on;

subplot(2 , 1 , 2)
errorbar(1 : length(cinNeg_mean) , cinNeg_mean , cinNeg_sdev , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
%boxplot((cinNeg2002 .* 100) , 'Color' , 'k' , 'Whisker' , 5);
plot(1 : length(cinNeg_mean) , mean((cin1Neg2002 .* 100),1)' , 'g-' , ...
    1 : length(cinNeg_mean) , min((cin1Neg2002 .* 100),[],1)' , 'g--' , ...
    1 : length(cinNeg_mean) , max((cin1Neg2002 .* 100),[],1)' , 'g--' , 'LineWidth' , 1.5);
plot(1 : length(cinNeg_mean) , mean((cin2Neg2002 .* 100),1)' , 'c-' , ...
    1 : length(cinNeg_mean) , min((cin2Neg2002 .* 100),[],1)' , 'c--' , ...
    1 : length(cinNeg_mean) , max((cin2Neg2002 .* 100),[],1)' , 'c--' , 'LineWidth' , 1.5);
plot(1 : length(cinNeg_mean) , mean((cin3Neg2002 .* 100),1)' , 'm-' , ...
    1 : length(cinNeg_mean) , min((cin3Neg2002 .* 100),[],1)' , 'm--' , ...
    1 : length(cinNeg_mean) , max((cin3Neg2002 .* 100),[],1)' , 'm--' , 'LineWidth' , 1.5);
hold all;
formatTemp = ones(1 , length(1 : length(cinNeg_mean)));
plot(1 : length(cinNeg_mean) , formatTemp.*7 , ':g'); hold all;
plot(1 : length(cinNeg_mean) , formatTemp.*2 , ':b'); hold all;
plot(1 : length(cinNeg_mean) , formatTemp.*3 , ':y');
legend('(McDonald, 2014) Observed CIN2/3 Cape Town: mean, 2SD' , 'Model, CIN1: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, CIN2: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, CIN3: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    '(Kuhn, 2020) CIN1 Cape Town' , '(Kuhn, 2020) CIN2 Cape Town' , '(Kuhn 2020) CIN3 Cape Town');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence (%)')
title('CIN 2/3 Prevalence in 2002 - Females, HIV-')
ylim([0 10])
grid on;

%% CIN1, CIN2, CIN3 prevalence for All HR HPV types combined by HIV status in 2015 vs. Kuhn 2020 data (calibration)

% Calibration error bars
cinPos_mean = cinPos2015_dObs(: , 2) .* 100;
cinPos_sdev = (cinPos2015_dObs(: , 3).^(1/2)).*2 .* 100;
cinNeg_mean = cinNeg2015_dObs(: , 2) .* 100;
cinNeg_sdev = (cinNeg2015_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
subplot(1 , 2 , 1);
errorbar(1 : length(cinPos_mean) , cinPos_mean , cinPos_sdev , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
plot(1 , mean((cin1Pos2015 .* 100),1)' , 'ko' , ...
   1 , min((cin1Pos2015 .* 100),[],1)' , 'k*' , ...
    1 , max((cin1Pos2015 .* 100),[],1)' , 'k*');
plot(2 , mean((cin2Pos2015 .* 100),1)' , 'ko' , ...
    2 , min((cin2Pos2015 .* 100),[],1)' , 'k*' , ...
    2 , max((cin2Pos2015 .* 100),[],1)' , 'k*');
plot(3 , mean((cin3Pos2015 .* 100),1)' , 'ko' , ...
    3 , min((cin3Pos2015 .* 100),[],1)' , 'k*' , ...
    3 , max((cin3Pos2015 .* 100),[],1)' , 'k*');
legend('(Kuhn, 2020) Observed CIN1,2,3 Cape Town: mean, 2SD' , 'Model, CIN1: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, CIN2: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, CIN3: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
set(gca , 'xtick' , 1 : 3 , 'xtickLabel' , {'CIN1' , 'CIN2' , 'CIN3'});
ylabel('Prevalence (%)')
title('Prevalence in 2015 - Females aged 30-65, HIV+')
ylim([0 25])
grid on;

subplot(1 , 2 , 2)
errorbar(1 : length(cinNeg_mean) , cinNeg_mean , cinNeg_sdev , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
plot(1 , mean((cin1Neg2015 .* 100),1)' , 'ko' , ...
    1 , min((cin1Neg2015 .* 100),[],1)' , 'k*' , ...
    1 , max((cin1Neg2015 .* 100),[],1)' , 'k*');
plot(2 , mean((cin2Neg2015 .* 100),1)' , 'ko' , ...
    2 , min((cin2Neg2015 .* 100),[],1)' , 'k*' , ...
    2 , max((cin2Neg2015 .* 100),[],1)' , 'k*');
plot(3 , mean((cin3Neg2015 .* 100),1)' , 'ko' , ...
    3 , min((cin3Neg2015 .* 100),[],1)' , 'k*' , ...
    3 , max((cin3Neg2015 .* 100),[],1)' , 'k*');
legend('(Kuhn, 2020) Observed CIN1,2,3 Cape Town: mean, 2SD' , 'Model, CIN1: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, CIN2: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, CIN3: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
set(gca , 'xtick' , 1 : 3 , 'xtickLabel' , {'CIN1' , 'CIN2' , 'CIN3'});
ylabel('Prevalence (%)')
title('Prevalence in 2015 - Females aged 30-65, HIV-')
ylim([0 10])
grid on;


%% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

%% Cervical cancer incidence in 2012 by age vs. Globocan 2012 data and other sources (calibration)
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

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
meanObs = ccInc2012_dObs(: , 2);
sdevObs = (ccInc2012_dObs(: , 3).^(1/2)).*2;

figure;    
% Plot observed data
plot(4 : age , combined_ub , 'r-' , 4 : age , combined_lb , 'r-'); % , ...
    %'Color' , [0.9290, 0.6940, 0.1250]);
hold on; 
errorbar(8 : age , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
%boxplot(ccInc2012 , 'Color' , 'k' , 'Whisker' , 5);
% General
plot(1 : age , mean(ccInc2012,1)' , 'k-' , ...
    1 : age , min(ccInc2012,[],1)' , 'k--' , ...
    1 : age , max(ccInc2012,[],1)' , 'k--' , 'LineWidth' , 1.5);
hold all;
% HIV-negative
plot(1 : age , mean(ccInc2012neg,1)' , 'g-' , ...
    1 : age , min(ccInc2012neg,[],1)' , 'g--' , ...
    1 : age , max(ccInc2012neg,[],1)' , 'g--' , 'LineWidth' , 1.5);
hold all;
% HIV-positive untreated
plot(1 : age , mean(ccInc2012pos,1)' , 'm-' , ...
    1 : age , min(ccInc2012pos,[],1)' , 'm--' , ...
    1 : age , max(ccInc2012pos,[],1)' , 'm--' , 'LineWidth' , 1.5);
hold all;
% HIV-positive on ART
plot(1 : age , mean(ccInc2012art,1)' , 'c-' , ...
    1 : age , min(ccInc2012art,[],1)' , 'c--' , ...
    1 : age , max(ccInc2012art,[],1)' , 'c--' , 'LineWidth' , 1.5);
hold all;
xlabel('Age Group'); ylabel('Cervical cancer incidence per 100K');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
%ylim([0 160]);
title(['Cervical Cancer Incidence in 2012']);
legend('Combined SA: upper bound' , 'Combined SA: lower bound' , ...
    '(Globocan, 2012) Observed SA: mean, 2SD (ages 15-39 comb)' , 'Model, general: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, HIV-negative: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, WLWHIV untreated: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, WLWHIV on ART: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
grid on;

%% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************

%% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)

% HPV infected
% Calibration error bars
meanObs = hpv_dist_dObs(: , 2) .* 100;
sdevObs = (hpv_dist_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
subplot(2,3,1)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
plot(2012 , 46.82 , 'k*');
hold on;
plot(2012 , 53.18 , 'b*');
hold on;
%boxplot(hpv_vax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
%    'Colors' , 'k' , 'Whisker' , 5);
plot(typeDistYearVec , mean((hpv_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((hpv_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((hpv_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
%boxplot(hpv_nonVax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
%    'Colors' , 'b' , 'Whisker' , 5);
plot(typeDistYearVec , mean((hpv_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((hpv_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((hpv_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('HPV');
ylim([0 100]);
xlim([2010 2015]);
grid on;

% CIN1
% Calibration error bars
meanObs = cin1_dist_dObs(: , 2) .* 100;
sdevObs = (cin1_dist_dObs(: , 3).^(1/2)).*2 .* 100;

subplot(2,3,2)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
plot(2012 , 51.92 , 'k*');
hold on;
plot(2012 , 48.08 , 'b*');
hold on;
%boxplot(cin1_vax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
%    'Colors' , 'k' , 'Whisker' , 5);
plot(typeDistYearVec , mean((cin1_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cin1_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cin1_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
%boxplot(cin1_nonVax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
%    'Colors' , 'b' , 'Whisker' , 5);
plot(typeDistYearVec , mean((cin1_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((cin1_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((cin1_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('CIN1');
grid on;

% CIN2
subplot(2,3,3);
plot(2012 , 62.81 , 'k*');
hold on;
plot(2012 , 37.19 , 'b*');
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('CIN2');
grid on;

% CIN3
% Calibration error bars
meanObs = cin3_dist_dObs(: , 2) .* 100;
sdevObs = (cin3_dist_dObs(: , 3).^(1/2)).*2 .* 100;

subplot(2,3,4)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
plot(2012 , 73.71 , 'k*');
hold on;
plot(2012 , 26.29 , 'b*');
hold on;
%boxplot(cin3_vax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
%    'Colors' , 'k' , 'Whisker' , 5);
plot(typeDistYearVec , mean((cin3_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cin3_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cin3_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
%boxplot(cin3_nonVax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
%    'Colors' , 'b' , 'Whisker' , 5);
plot(typeDistYearVec , mean((cin3_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((cin3_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((cin3_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('CIN3');
grid on;

% CC
% Calibration error bars
meanObs = cc_dist_dObs(: , 2) .* 100;
sdevObs = (cc_dist_dObs(: , 3).^(1/2)).*2 .* 100;

subplot(2,3,5)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
plot(2012 , 85.78 , 'k*');
hold on;
plot(2012 , 14.22 , 'b*');
hold on;
%boxplot(cc_vax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
%    'Colors' , 'k' , 'Whisker' , 5);
plot(typeDistYearVec , mean((cc_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cc_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cc_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
%boxplot(cc_nonVax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
%    'Colors' , 'b' , 'Whisker' , 5);
plot(typeDistYearVec , mean((cc_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((cc_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((cc_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('Cervical Cancer')
legend('Observed 2012: mean, 2SD' , 'Observed 2012- 9v' , 'Observed 2012- non-9v' , ...
    'Model- 9v: 25-sets mean' , 'Model- 9v: 25-sets minimum' , 'Model- 9v: 25-sets maximum' , ...
    'Model- non-9v: 25-sets mean' , 'Model- non-9v: 25-sets minimum' , 'Model- non-9v: 25-sets maximum');
grid on;
