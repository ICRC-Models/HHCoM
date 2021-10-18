function[] = showResults_multSims_CIs(fileList , nBestFits , tstep_abc , date)

%% Load parameters and results
paramDir = [pwd , '\Params\'];

[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd, annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , partnersMmult, maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4,...
    mue , mue2 , mue3 , mue4 , epsA_vec , epsR_vec , yr , ...
    hivOn , betaHIV_mod ,  muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , c2c3Mults , c1c2Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , ...
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
    infhpvNonVaxInds , fromVaxNoScrnInds , fromVaxScrnInds , toNonVaxNoScrnInds , ...
    toNonVaxScrnInds , ageInd , riskInd , ...
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3  , dFertMat3, d_partnersMmult, riskAdj, d_riskAdj, ...
    deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2(1 , 0 , [] , [] , []);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)

% Indices of calib runs to plot
files = fileList;
nRuns = length(fileList);

% Initialize model output plots
% Total population size
popYearVec = unique(totPopSize_dObs(: ,1));
popYearVec_length = length(popYearVec);
popSize = zeros(nRuns , length(popYearVec));
% Population age distribution
popYearVecAge = unique(popAgeDist_dObs(: ,1));
popYearVecAge_length = length(popYearVecAge);
popProp = zeros(nRuns , length(popYearVecAge) , age);
% HIV prevalence
hivPrevAFvec = {4 , 5 , 6 , 7 , 8 , 9 , 10};
hivPrevAFvec_length = length(hivPrevAFvec);
hivPrevAMvec = {4 , 5 , 6 , 7 , 8 , 9 , 10 , 11};
hivPrevAMvec_length = length(hivPrevAMvec);
hivPrevACvec = {4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13};
hivPrevACvec_length = length(hivPrevACvec);
hivYearVec = unique(hivPrevM_dObs(: ,1));
hivYearVec_length = length(hivYearVec);
hivAgeM = zeros(nRuns , 8 , length(hivYearVec));
hivAgeF = zeros(nRuns , 7 , length(hivYearVec));
hivAgeAll = zeros(nRuns , 10 , 1);
hiv_prev = zeros(nRuns,length([startYear:timeStep:currYear]));
% HPV Prevalence in high risk women
ageVecHpvHR = {[4:5],[6],[7:8],[9:10]};
ageVecHpvHR_length = length(ageVecHpvHR);
hpv_hiv = zeros(nRuns , 4);
hpv_hivNeg = zeros(nRuns , 4);
genPopHPV = zeros(nRuns,length([startYear:timeStep:currYear]));
hivPosHPV = genPopHPV; 
artHPV = genPopHPV; 
hivNegHPV = genPopHPV; 
% HPV prevalence in all women (not including CIN2+)
hivPrevFnoCinvec = {6 , 7 , 8 , 9 , 10 , 11};
hivPrevFnoCinvec_length = length(hivPrevFnoCinvec);
hpv_all = zeros(nRuns , 6);
% HPV prevalence in HIV+ women (including CIN)
ageVecHpvWCin = {[4:6], 7, 8, 9, [10:11]};
ageVecHpvWCin_length = length(ageVecHpvWCin);
hpv_hiv2009 = zeros(nRuns , 5);
% HPV prevalence by HIV status
hpvPrevAvec = {4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12};
hpvPrevAvec_length = length(hpvPrevAvec);
hpv2005 = zeros(nRuns,9,1);
hpvHIV2005 = hpv2005;
hpvNeg2005 = hpv2005;
% CIN1 prevalence among women aged 20-50 in 2010 by HIV status 
dVec = {[1:2],[3:8]};
dVec_length = length(dVec);
cin1_2010 = zeros(nRuns , 2);
% CIN2/3 prevalence among women aged 20-50 in 2010 by HIV status 
cin2_2010 = zeros(nRuns , 2);
% CIN2/3 prevalence among HIV+ women aged 20-50 in 2007
ageVecCin23Hiv = {[5:6],[7:8],[9:10]};
ageVecCin23Hiv_length = length(ageVecCin23Hiv);
cinPos2007 = zeros(nRuns , 3);
% CC incidence
ccInc2012 = zeros(nRuns , age);
ccInc2012neg = ccInc2012;
ccInc2012pos = ccInc2012;
ccInc2012art = ccInc2012;
% HPV/CIN/CC type distribution
typeDistYearVec = [2010 , 2011 , 2012 , 2013 , 2014 , 2015];
typeDistYearVec_length = length(typeDistYearVec);
cc_vax = zeros(nRuns , length(typeDistYearVec));
cc_nonVax = cc_vax;
cin3_vax = cc_vax;
cin3_nonVax = cc_vax;
cin1_vax = cc_vax;
cin1_nonVax = cc_vax;
hpv_vax = cc_vax;
hpv_nonVax = cc_vax;


resultsDir = [pwd , '/HHCoM_Results/'];
figuresDir = [pwd, '/HHCoM_Figures/'];
t = datetime();
groupDir = strcat(figuresDir , 'Figures_' , date , '_' , num2str(tstep_abc) , '_' , num2str(nBestFits) , 'fits' , '-' , ...
    int2str(t.Month) , '_' , int2str(t.Day) , '_' , int2str(t.Hour) , '_' , int2str(t.Minute));
mkdir(groupDir);


for j = 1 : nRuns
    % Load results
    pathModifier = files(j);
    histResult = load(strcat(resultsDir , pathModifier));

    %% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

    %% Population size over time vs. Kenya Population Census data (calibration)
    for tInd = 1 : popYearVec_length
        t = popYearVec(tInd);

        popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
        popSize(j , tInd) = sum(histResult.popVec(((t - startYear) * stepsPerYear +1) , popTot),2);
    end

    %% Population size by 5-year age groups over time vs. Kenya Population Census data (calibration)
    for tInd = 1 : popYearVecAge_length
        t = popYearVecAge(tInd);
        for aV = 1 : age
            a = aV;
            popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
            popProp(j , tInd , aV) = sum(histResult.popVec(((t - startYear) * stepsPerYear +1) , popAge),2) ./ ...
                sum(histResult.popVec(((t - startYear) * stepsPerYear +1) , popTot),2);
        end
    end

    %% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

    %% HIV prevalence by age and gender over time vs. KAIS and DHS (calibration)
    for tInd = 1 : hivYearVec_length
        t = hivYearVec(tInd);
        for aInd = 1 : hivPrevAFvec_length
            a = hivPrevAFvec{aInd};
            hivFInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            artFInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            totFInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            hivAgeF(j , aInd , tInd) =  (sum(histResult.popVec((t - startYear) * stepsPerYear +1 , hivFInds)) ...
                + sum(histResult.popVec((t - startYear) * stepsPerYear +1 , artFInds))) ...
                / sum(histResult.popVec((t - startYear) * stepsPerYear +1 , totFInds));
        end
    end

    for tInd = 1 : hivYearVec_length
        t = hivYearVec(tInd);
        for aInd = 1 : hivPrevAMvec_length
            a = hivPrevAMvec{aInd};
            hivMInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            artMInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            totMInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            hivAgeM(j , aInd , tInd) =  (sum(histResult.popVec((t - startYear) * stepsPerYear +1 , hivMInds)) ...
                + sum(histResult.popVec((t - startYear) * stepsPerYear +1 , artMInds))) ...
                / sum(histResult.popVec((t - startYear) * stepsPerYear +1 , totMInds));
        end
    end

    for tInd = 1 : hivYearVec_length
        t = hivYearVec(tInd);
        for aInd = 1 : hivPrevACvec_length
            a = hivPrevACvec{aInd};
            hivAllInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            artAllInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            totAllInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            hivAgeAll(j , aInd , 1) =  (sum(histResult.popVec((t - startYear) * stepsPerYear +1 , hivAllInds)) ...
                + sum(histResult.popVec((t - startYear) * stepsPerYear +1 , artAllInds))) ...
                / sum(histResult.popVec((t - startYear) * stepsPerYear +1 , totAllInds));
        end
    end

    %% Total HIV Prevalence
    % Total HIV positive
    hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , 1 : 2 , 4 : 11, 1 : risk));
    hivPop = sum(histResult.popVec(: , hivInds) , 2);
    artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , 1 : 2 , 4 : 11 , 1 : risk));
    art = sum(histResult.popVec(: , artInds) , 2);
    popTot = histResult.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 1 : 2 , 4 : 11 , 1 : risk)));
    hiv_prev(j,:) = ((hivPop + art) ./ sum(popTot , 2) * 100)';


    %% ********************************** HPV FIGURES **********************************************************************************************
    %% HPV Prevalence in high risk women
    yr = 2006;
    for aV = 1 : ageVecHpvHR_length
        a = ageVecHpvHR{aV};
        %HIV-positive
        hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 3)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 3))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 3));
        hpv_hiv(j , aV) = sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
            ./ sum(histResult.popVec((yr - startYear) * stepsPerYear+1 , ageInds));
        %HIV-negative
        hpvInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 3)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5, 1 , 1 : intervens , 2 , a , 3))]);
        ageInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 3));
        hpv_hivNeg(j , aV) = sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
            ./ sum(histResult.popVec((yr - startYear) * stepsPerYear+1 , ageInds));
    end

    %% HPV prevalence in all women (not including CIN2+) by age
    yr = 2000;
    for aInd = 1 : hivPrevFnoCinvec_length
        a = hivPrevFnoCinvec{aInd};
        hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 3 , 7] , 2 : 3 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(1 : disease, 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpv_all(j , aInd) = sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
            ./ sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    end

    %% HPV prevalence in HIV+ women (including CIN)
    yr = 2009;
    for aV = 1 : ageVecHpvWCin_length
        a = ageVecHpvWCin{aV}; % 15-19 ->  55-65
        hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpv_hiv2009(j , aV) = sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
            ./ sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    end

    %% HPV Prevalence by age in 2005 vs. Yamada data
    ageGroup = {'17 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
    yr = 2005;
    aVec = {18:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
    for aInd = 1 : hpvPrevAvec_length
        a = hpvPrevAvec{aInd};
        %a = aVec{aInd};

        hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpv2005(j , aInd , 1) = sum(histResult.popVec((yr - startYear) * stepsPerYear , hpvInds))...
            ./ sum(histResult.popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;

        % HIV+
        hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 :5  , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpvHIV2005(j , aInd , 1) = sum(histResult.popVec((yr - startYear) * stepsPerYear , hpvInds))...
            ./ sum(histResult.popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;

        % HIV-
        hpvInds = unique([toInd(allcomb(1 : 2 , 1  , 2  : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 , ...
            [1 : 5 , 7] , 2 : 5, 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(1 : 2 , 1  , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpvNeg2005(j , aInd , 1) = sum(histResult.popVec((yr - startYear) * stepsPerYear , hpvInds))...
            ./ sum(histResult.popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;
    end

    %% HPV prevalence over time by HIV status among women
    % General pop
    hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 4 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 4 : age , 1 : risk))]);
    hpvPop = sum(histResult.popVec(: , hpvInds) , 2);
    popTot = histResult.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 4 : age , 1 : risk)));
    genPopHPV(j, :) = 100 * hpvPop ./ sum(popTot , 2);
    
    % HIV+, not on ART 
    hpvHivInds = unique([toInd(allcomb(3 : 7 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 4 : age , 1 : risk)); toInd(allcomb(3 : 7 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 4 : age , 1 : risk))]);
    hpvHivPop = sum(histResult.popVec(: , hpvHivInds) , 2);
    popHivTot = histResult.popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 4 : age , 1 : risk)));
    hivPosHPV(j, :) = 100 * hpvHivPop ./ sum(popHivTot , 2);
    
    %ART
    hpvArtInds = unique([toInd(allcomb(8 , 6 , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 4 : age , 1 : risk)); toInd(allcomb(8 , 6 , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 4 : age , 1 : risk))]);
    hpvArtPop = sum(histResult.popVec(: , hpvArtInds) , 2);
    popArtTot = histResult.popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 4 : age , 1 : risk)));
    artHPV(j, :) =  100 * hpvArtPop ./ sum(popArtTot , 2);
    
    %HIV-
    hpvHivNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 4 : age , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 4 : age , 1 : risk))]);
    hpvHivNegPop = sum(histResult.popVec(: , hpvHivNegInds) , 2);
    popHivNegTot = histResult.popVec(: , toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 4 : age , 1 : risk)));
    hivNegHPV(j, :) = 100 * hpvHivNegPop ./ sum(popHivNegTot , 2);
    
    %% ********************************** CIN FIGURES *********************************************************************************************
    %% CIN2/3 prevalence among HIV+ women aged 20-50 in 2007
    yr = 2007;
    for aV = 1 : ageVecCin23Hiv_length
        a = ageVecCin23Hiv{aV};
        % HIV-positive (on and not on ART)
        cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        cinPos2007(j , aV) = (sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
            ./ sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    end

    %% CIN1 prevalence among women aged 20-50 in 2010 by HIV status 
    yr = 2010;
    for dV = 1 : dVec_length
        d = dVec{dV};
        cinInds = unique([toInd(allcomb(d , 1 : viral , 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , 5 : 10 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
            [1 : 3 , 7] , 3 , 1 , 1 : intervens , 2 , 5 : 10 , 1 : risk))]);
        ageInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 5 : 10 , 1 : risk));
        cin1_2010(j , dV) = (sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
            ./ sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    end

    %% CIN2/3 prevalence among women aged 20-50 in 2010 by HIV status 
    yr = 2010;
    for dV = 1 : dVec_length
        d = dVec{dV};
        cinInds = unique([toInd(allcomb(d , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 5 : 10 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 5 : 10 , 1 : risk))]);
        ageInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 5 : 10 , 1 : risk));
        cin2_2010(j , dV) = (sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
            ./ sum(histResult.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    end


    %% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

    %% Cervical cancer incidence in 2012 by age vs. Globocan 2012 data and other sources (calibration)
    incTimeSpan = [((2012 - startYear) * stepsPerYear +1) : ((2012 - startYear) * stepsPerYear +6)];
    fac = 10 ^ 5;

    for aV = 1 : age
        a = aV;
        % General population
        allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012(j , aV) = ...
            (annlz(sum(sum(sum(histResult.newCC(incTimeSpan , : , a , :),2),3),4)) ./ ...
            (annlz(sum(histResult.popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);

        % HIV-negative
        allFneg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012neg(j , aV) = ...
            (annlz(sum(sum(sum(histResult.newCC(incTimeSpan , 1 : 2 , a , :),2),3),4)) ./ ...
            (annlz(sum(histResult.popVec(incTimeSpan , allFneg) , 2) ./ stepsPerYear)) * fac);

        % HIV-positive untreated
        allFpos = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012pos(j , aV) = ...
            (annlz(sum(sum(sum(histResult.newCC(incTimeSpan , 3 : 7 , a , :),2),3),4)) ./ ...
            (annlz(sum(histResult.popVec(incTimeSpan , allFpos) , 2) ./ stepsPerYear)) * fac);

        % HIV-positive on ART
        allFart = toInd(allcomb(8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012art(j , aV) = ...
            (annlz(sum(sum(sum(histResult.newCC(incTimeSpan , 8 , a , :),2),3),4)) ./ ...
            (annlz(sum(histResult.popVec(incTimeSpan , allFart) , 2) ./ stepsPerYear)) * fac);
    end


    %% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************

    %% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)
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

    cin2Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 4 , [1 : 4 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin2Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 3 , 7] , 4 , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin2Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 4 , [1 : 4 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 3 , 7] , 4 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

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

    for i = 1 : typeDistYearVec_length
        cc_vax(j , i) = sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_vax) , 2)...
            ./ sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);
        cc_nonVax(j , i) = sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_nonVax) , 2)...
            ./ sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);

        cin3_vax(j , i) = sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_vax) , 2)...
            ./ sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);
        cin3_nonVax(j , i) = sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_nonVax) , 2)...
            ./ sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);

        cin2_vax(j , i) = sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin2Inds_vax) , 2)...
            ./ sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin2Inds_tot) , 2);
        cin2_nonVax(j , i) = sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin2Inds_nonVax) , 2)...
            ./ sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin2Inds_tot) , 2);

        cin1_vax(j , i) = sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_vax) , 2)...
            ./ sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_tot) , 2);
        cin1_nonVax(j , i) = sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_nonVax) , 2)...
            ./ sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_tot) , 2);

        hpv_vax(j , i) = sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_vax) , 2)...
            ./ sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_tot) , 2);
        hpv_nonVax(j , i) = sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_nonVax) , 2)...
            ./ sum(histResult.popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_tot) , 2);
    end
end

%% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

%% Population size over time vs. Kenya Population Census data (calibration)
% Calibration error bars
meanObs = totPopSize_dObs(: , 2);
sdevObs = (totPopSize_dObs(: , 3).^(1/2)).*2;

figure;
errorbar(totPopSize_dObs(: , 1) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold all;
plot(totPopSize_dObs(: , 1) , mean(popSize,1)' , 'k-' , ...
    totPopSize_dObs(: , 1) , min(popSize,[],1)' , 'k--' , ...
    totPopSize_dObs(: , 1) , max(popSize,[],1)' , 'k--' , 'LineWidth' , 1.5);
title('Kenya Population Size Ages 0-79')
xlabel('Year'); ylabel('Individuals')
xlim([1970 2020]); ylim([0 (4*10^7)]);
legend('(Kenya Population Census) Observed Kenya, ages 0-79: mean, 2SD' , ...
    ['Model, ages 0-79: ' , num2str(nRuns) , '-sets mean'] , ['Model, ages 0-79: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model, ages 0-79: ' , num2str(nRuns) , '-sets maximum'] , 'Location' , 'NorthWest');
grid on;
saveas(gcf,strcat(groupDir,'/populationByYear.fig'));
close(gcf);

%% Population size by 5-year age groups over time vs. Kenya Population Census data (calibration)
% Calibration error bars
meanObs = [popAgeDist_dObs(1:16 , 2) , popAgeDist_dObs(17:32 , 2) , popAgeDist_dObs(33:48 , 2) , popAgeDist_dObs(49:64 , 2)]';
sdevObs = ([popAgeDist_dObs(1:16 , 3) , popAgeDist_dObs(17:32 , 3) , popAgeDist_dObs(33:48 , 3) , popAgeDist_dObs(49:64 , 3)]'.^(1/2)).*2;

figure;
subplot(1,3,1);
set(gca,'ColorOrderIndex',1)
calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , ...
    unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) ,];
errorbar(calibYrs , meanObs(: , 1:7) , sdevObs(: , 1:7) , ...
    's' , 'LineWidth' , 1.5);
hold on;
for a = 1 : 7
    set(gca,'ColorOrderIndex',a)
    plot(popYearVecAge , mean(popProp(: , : , a),1)' , '-' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a)
    plot(popYearVecAge , min(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a)
    plot(popYearVecAge , max(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
end
ylim([0.05 0.25]);
ylabel('Population proportion by age'); xlabel('Year');
legend('(Kenya Population Census) Observed Kenya, ages 0-4: mean, 2SD' , ...
    'ages 5-9: mean, 2SD' , 'ages 10-14: mean, 2SD' , 'ages 15-19: mean, 2SD' , ...
    'ages 20-24: mean, 2SD' , 'ages 25-29: mean, 2SD' , 'ages 30-34: mean, 2SD' , ...
    ['Model, ages 0-4: ' , num2str(nRuns) , '-sets mean'] , ...
    ['ages 0-4: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['ages 0-4: ' , num2str(nRuns) , '-sets maximum'] , ...
    '...' , 'Location' , 'north');

subplot(1,3,2);
set(gca,'ColorOrderIndex',1)
calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , ...
    unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) ,];
errorbar(calibYrs , meanObs(: , 8:14) , sdevObs(: , 8:14) , ...
    's' , 'LineWidth' , 1.5);
hold on;
for a = 8 : 14
    set(gca,'ColorOrderIndex',a-7)
    plot(popYearVecAge , mean(popProp(: , : , a),1)' , '-' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-7)
    plot(popYearVecAge , min(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-7)
    plot(popYearVecAge , max(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
end
ylim([0.0 0.12]);
ylabel('Population proportion by age'); xlabel('Year');
legend('(Kenya Population Census) Observed Kenya, ages 35-39: mean, 2SD' , ...
    'ages 40-44: mean, 2SD' , 'ages 45-49: mean, 2SD' , 'ages 50-54: mean, 2SD' , ...
    'ages 55-59: mean, 2SD' , 'ages 60-64: mean, 2SD' , 'ages 65-69: mean, 2SD' , ...
    ['Model, ages 35-39: ' , num2str(nRuns) , '-sets mean'] , ...
    ['ages 35-39: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['ages 35-39: ' , num2str(nRuns) , '-sets maximum'] , '...' , 'Location' , 'north');


subplot(1,3,3);
set(gca,'ColorOrderIndex',1)
calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1))];
errorbar(calibYrs , meanObs(: , 15:16) , sdevObs(: , 15:16) , ...
    's' , 'LineWidth' , 1.5);
hold on;
for a = 15 : 16
    set(gca,'ColorOrderIndex',a-14)
    plot(popYearVecAge , mean(popProp(: , : , a),1)' , '-' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-14)
    plot(popYearVecAge , min(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-14)
    plot(popYearVecAge , max(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
end
ylim([0.0 0.04]);
ylabel('Population proportion by age'); xlabel('Year'); %title('KZN age distribution in 5-year groups');
legend('(Kenya Population Census) Observed Kenya, ages 70-74: mean, 2SD' , ...
    'ages 75-79: mean, 2SD' , ['Model, ages 70-74: ' , num2str(nRuns) , '-sets mean'] , ...
    ['ages 70-74: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['ages 70-74: ' , num2str(nRuns) , '-sets maximum'] , '...' , 'Location' , 'north');

saveas(gcf,strcat(groupDir,'/ageGroupSizeByTime.fig'));
close(gcf);

%% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

%% HIV prevalence by age and gender over time vs. DHS and KAIS (calibration)
% Calibration error bars
hivM(: , 1) = hivPrevM_dObs(1:24 , 2) .* 100; % mean
hivM(: , 2) = (hivPrevM_dObs(1:24 , 3).^(1/2)).*2 .* 100; % calibration SD
hivF(: , 1) = hivPrevF_dObs(1:21 , 2) .* 100; % mean
hivF(: , 2) = (hivPrevF_dObs(1:21 , 3).^(1/2)).*2 .* 100; % calibration SD

hivAll(: , 1) = hivPrevAll_dObs(: , 2) .* 100; % mean
hivAll(: , 2) = (hivPrevAll_dObs(: , 3).^(1/2)).*2 .* 100; % calibration SD

ageGroup = {'15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74'};
ageEndVec = {11-3 , 10-3 , 13-3};
genVecLabels = {'Male' , 'Female' , 'Male+Female Combined'};
genVec = {1 , 2 , 1:2};
yrColorVec = {'k' , 'b' , 'r'};

figure;
for gInd = 1 : length(genVec)
    g = genVec{gInd};
    if gInd == 1
        hivPrevs = hivM;
        hivModel = hivAgeM(: , : , (1:(end-1)));
    elseif gInd == 2
        hivPrevs = hivF;
        hivModel = hivAgeF(: , : , (1:(end-1)));
    elseif gInd == 3
        hivPrevs = hivAll;
        hivModel = hivAgeAll;
    end

    subplot(3 , 1 , gInd);
    for i = 1 : size(hivModel , 3)
        hold all;      
        errorbar(1 : ageEndVec{gInd} , hivPrevs(((i-1)*ageEndVec{gInd})+1 : ageEndVec{gInd}*i , 1) , ...
            hivPrevs(((i-1)*ageEndVec{gInd})+1 : ageEndVec{gInd}*i , 2) , ...
            's' , 'LineWidth' , 1.5 , 'Color' , yrColorVec{i});
        hold all;
        plot(1 : ageEndVec{gInd} , mean((squeeze(hivModel(: , : , i)) .* 100),1)' , '-' , 'LineWidth' , 1.5 , 'Color' , yrColorVec{i});
        hold all;
        plot(1 : ageEndVec{gInd} , min((squeeze(hivModel(: , : , i)) .* 100),[],1)' , '--' , 'LineWidth' , 1.5 , 'Color' , yrColorVec{i});
        hold all;
        plot(1 : ageEndVec{gInd} , max((squeeze(hivModel(: , : , i)) .* 100),[],1)' , '--' , 'LineWidth' , 1.5 , 'Color' , yrColorVec{i});
        set(gca , 'xtickLabel' , ageGroup); set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
    end
   xlabel('Age Group'); ylabel('HIV Prevalence (%)'); title(genVecLabels{gInd});
   ylim([0 50])
   grid on;
if (gInd == 1) || (gInd == 2)
      legend('(DHS/KAIS) Observed Kenya, 2003: mean, 2SD' , ...
      ['Model, 2003: ' , num2str(nRuns) , '-sets mean'] , ...
      ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
      ['Model: ' , num2str(nRuns) , '-sets maximum'] , ...
       '(DHS/KAIS) Observed Kenya, 2007: mean, 2SD' , ...
      ['Model, 2007: ' , num2str(nRuns) , '-sets mean'] , ...
      ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
      ['Model: ' , num2str(nRuns) , '-sets maximum'] , ...
       '(DHS/KAIS) Observed Kenya, 2008: mean, 2SD' , ...
       ['Model, 2009: ' , num2str(nRuns) , '-sets mean'] , ...
       ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
       ['Model: ' , num2str(nRuns) , '-sets maximum']);
elseif gInd == 3
        legend('(DHS/KAIS) Observed Kenya: mean, 2SD' , ...
        ['Model, 2012: ' , num2str(nRuns) , '-sets mean'] , ...
        ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
        ['Model: ' , num2str(nRuns) , '-sets maximum']);
end
end

saveas(gcf,strcat(groupDir,'/hivPrevByAgeGender.fig'));
close(gcf);

%% HIV Prevalence
% Compared to observed HIV data

file = [pwd , '/Config/Kenya_parameters_Feb20.xlsx'];
HIV_Ken_spectrum = xlsread(file , ['HIV ' 'prevalence'] , 'B184:E213');
DHS_KAIS = [2003 6.7 5.8 7.6;
    2007 7.1 6.6 7.9;
    2009 6.4 5.4 7.3;
    2012 5.6 4.9 6.3];

figure()

plot([startYear : timeStep : currYear] , mean(hiv_prev) , HIV_Ken_spectrum(: , 1)' , HIV_Ken_spectrum(: , 2)' , '+', ...
DHS_KAIS(:, 1)',  DHS_KAIS(:, 2)', 'o')

hold on
plot([startYear : timeStep : currYear], max(hiv_prev,[],1),'--')
hold on
plot([startYear : timeStep : currYear], min(hiv_prev,[],1),'--')
hold on 
xlabel('Year'); ylabel('Proportion of Population (%)'); 
title({'HIV Prevalence (Ages 15-54)'})
legend('Model Mean' , 'Kenya (Spectrum)', 'Kenya (DHS/KAIS)','Model Max','Model Min')
xlim([1980 2020])
ylim([0, 15])

saveas(gcf,strcat(groupDir,'/hivPrevByTime.fig'));
close(gcf);


%% ********************************** HPV FIGURES **********************************************************************************************

%% HPV prevalence over time among women 
figure;
hpvTimeVec = [startYear:timeStep:currYear]';
plot(hpvTimeVec , mean(genPopHPV,1)' , 'k-' , ...
    hpvTimeVec , min(genPopHPV,[],1)' , 'k--' , ...
    hpvTimeVec , max(genPopHPV,[],1)' , 'k--' , ...
    hpvTimeVec , mean(hivPosHPV, 1)', 'r-' , ...
    hpvTimeVec , min(hivPosHPV,[],1)' , 'r--' , ...
    hpvTimeVec , max(hivPosHPV,[],1)' , 'r--' , ...
    hpvTimeVec , mean(artHPV, 1)', 'g-' , ...
    hpvTimeVec , min(artHPV,[],1)' , 'g--' , ...
    hpvTimeVec , max(artHPV,[],1)' , 'g--' , ...
    hpvTimeVec , mean(hivNegHPV, 1)', 'b-' , ...
    hpvTimeVec , min(hivNegHPV,[],1)' , 'b--' , ...
    hpvTimeVec , max(hivNegHPV,[],1)' , 'b--' , 'LineWidth' , 1.5);
title('HPV prevalence over time in women by HIV status')
xlabel('Year'); ylabel('Prevalence')
xlim([1970 2020]); 
legend(['General pop: ' , num2str(nRuns) , '-sets mean'] , ['General pop: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['General pop: ' , num2str(nRuns) , '-sets maximum'] , ...
    ['HIV-positive: ' , num2str(nRuns) , '-sets mean'] , ['HIV-positive: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['HIV-positive: ' , num2str(nRuns) , '-sets maximum'] , ...
    ['On ART: ' , num2str(nRuns) , '-sets mean'] , ['On ART: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['On ART: ' , num2str(nRuns) , '-sets maximum'] , ...
    ['HIV-negative: ' , num2str(nRuns) , '-sets mean'] , ['HIV-negative: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['HIV-negative: ' , num2str(nRuns) , '-sets maximum'] , 'Location' , 'NorthWest');
grid on;
saveas(gcf,strcat(groupDir,'/hpvPrevByTime.fig'));
close(gcf);


%% HPV Prevalence in high risk women
ageGroup = {'15 -24' , '25 - 29' , '30 - 39' , '40 - 49'};

% Calibration error bars
meanObs = hpv_hiv_dObs(: , 2) .* 100;
sdevObs = (hpv_hiv_dObs(: , 3).^(1/2)).*2 .* 100;
meanNeg = hpv_hivNeg_dObs(: , 2) .* 100;
sdevNeg = (hpv_hivNeg_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
subplot(2,1,1);
errorbar(1 : length(meanObs) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(meanObs) , mean((hpv_hiv .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_hiv .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_hiv .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('(Luchter, 2010) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
title('hrHPV Prevalence in 2006 - High-risk females, HIV+');
grid on;

subplot(2,1,2);
errorbar(1 : length(meanObs) , meanNeg , sdevNeg , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(meanObs) , mean((hpv_hivNeg .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_hivNeg .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_hivNeg .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('(Luchter, 2010) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
title('hrHPV Prevalence in 2006 - High-risk females, HIV-');
grid on;

saveas(gcf,strcat(groupDir,'/hpvPrevHighRiskWomen.fig'));
close(gcf);

%% HPV prevalence in all women (not including CIN2+) by age
ageGroup = {'25 - 29' , '30 - 34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54'};

% Calibration error bars
meanObs = hpv_all_dObs(: , 2) .* 100;
sdevObs = (hpv_all_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
errorbar(1 : length(meanObs) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(meanObs) , mean((hpv_all .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_all .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_all .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('(DeVuyst, 2003) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
title('hrHPV Prevalence in 2000 - All females, not including CIN2+');
grid on;

saveas(gcf,strcat(groupDir,'/hpvPrevAllWomen.fig'));
close(gcf);
%% HPV prevalence in HIV+ women (including CIN)
ageGroup = {'15-19' , '20-24' , '25-29' , '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' }; 

% Calibration error bars
meanObs = hpv_hiv2009_dObs(: , 2) .* 100;
sdevObs = (hpv_hiv2009_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
errorbar(1 : length(meanObs) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(meanObs) , mean((hpv_hiv2009 .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_hiv2009 .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_hiv2009 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('(DeVuyst, 2012; BJC) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
title('hrHPV Prevalence in 2009 - Females, HIV+');
grid on;

saveas(gcf,strcat(groupDir,'/hpvPrevForHivPos.fig'));
close(gcf);

%Yamada,2009: HIV+ family planning/gyn health center, any HPV prevalence, but >50% types tested for

% were HR: 16,18,31,33,35,39,45,51,52,56,58,59,67,68,82,(6,11,13,44),(26, 30, 53, 66,70)

hpvHivObs = [NaN NaN NaN
0.66 0.58 0.75
0.67 0.56 0.79
0.57 0.45 0.71
0.57 0.45 0.71
0.56 0.43 0.68
NaN NaN NaN
NaN NaN NaN
NaN NaN NaN];

%Yamada 2009 : HIV- family planning/gyn health center

hpvNegObs = [NaN NaN NaN
0.39 0.25 0.52
0.30 0.12 0.47
0.26 0.13 0.39
0.26 0.13 0.39
0.24 0.10 0.39
NaN NaN NaN
NaN NaN NaN
NaN NaN NaN];

hpvHivObs = hpvHivObs * 100;
hpvNegObs = hpvNegObs * 100;

% hpvNegObs = hpvNegObs * 100;

figure()

% plot(1 : length(hpv2002) , hpv2002 , 'o-')

% hold all;

plot(1 : size(hpvHIV2005,2) , mean(hpvHIV2005,'omitnan') , 'o-');
hold all;
plot(1 : size(hpvNeg2005,2) , mean(hpvNeg2005) , 'o-')
hold all;
set(gca , 'xtickLabel' , ageGroup); 
hold on;
errorbar([1 : length(hpvHivObs)]' , hpvHivObs(:, 1) , ...
    hpvHivObs(:,1)-hpvHivObs(:, 2) , ...
    hpvHivObs(:,3)-hpvHivObs(:,1) , ...
    'x--' , 'LineWidth' , 1.5);
hold on;
errorbar([1 : length(hpvNegObs)]' , hpvNegObs(:, 1) , ...
    hpvNegObs(:,1)-hpvNegObs(:, 2) , ...
    hpvNegObs(:,3)-hpvNegObs(:,1) , ...
    'x--' , 'LineWidth' , 1.5);

set(gca , 'xtick' , 1 : length(hpvHivObs) , 'xtickLabel' , ageGroup);

legend('Model HIV-pos' , 'Model HIV-neg' ,...
'Obs HIV-pos: Yamada/Luchters','Obs HIV-neg: Yamada/Luchters')

xlabel('Age Group'); ylabel('hrHPV Prevalence (%)')

title(['Kenya HPV prevalence in women in 2005'])

ylim([0 100])

saveas(gcf,strcat(groupDir,'/hpvPrevByHiv.fig'));
close(gcf);



%% ********************************** CIN FIGURES *********************************************************************************************

%% CIN2/3 prevalence among HIV+ women aged 20-50 in 2007
ageGroup = {'20 - 29' , '30 - 39' , '40 - 49'};

% Calibration error bars
cinPos_mean = cinPos2007_dObs(: , 2) .* 100;
cinPos_sdev = (cinPos2007_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
errorbar(1 : length(cinPos_mean) , cinPos_mean , cinPos_sdev , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(cinPos_mean) , mean((cinPos2007 .* 100),1)' , 'k-' , ...
    1 : length(cinPos_mean) , min((cinPos2007 .* 100),[],1)' , 'k--' , ...
    1 : length(cinPos_mean) , max((cinPos2007 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
legend('(Hutchko) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence (%)')
title('CIN 2/3 Prevalence in 2007 - Females, HIV+')
ylim([0 25])
grid on;

saveas(gcf,strcat(groupDir,'/cin23Prev2007.fig'));
close(gcf);

%% CIN1 prevalence among women aged 20-50 in 2010 by HIV status 
diseaseGroup = {'HIV-negative' , 'HIV-positive'};

% Calibration error bars
cin1_mean = cin1_2010_dObs(: , 2) .* 100;
cin1_sdev = (cin1_2010_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
errorbar(1 : length(cin1_mean) , cin1_mean , cin1_sdev , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(cin1_mean) , mean((cin1_2010 .* 100),1)' , 'k-' , ...
    1 : length(cin1_mean) , min((cin1_2010 .* 100),[],1)' , 'k--' , ...
    1 : length(cin1_mean) , max((cin1_2010 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
legend('(Sweet, 2020) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
set(gca , 'xtick' , 1 : length(diseaseGroup) , 'xtickLabel' , diseaseGroup);
xlabel('HIV Status'); ylabel('CIN1 Prevalence (%)')
title('CIN1 Prevalence in 2010 - Females aged 20-50')
ylim([0 25])
grid on;

saveas(gcf,strcat(groupDir,'/cin1ByHiv2007.fig'));
close(gcf);

%% CIN2/3 prevalence among women aged 20-50 in 2010 by HIV status 
diseaseGroup = {'HIV-negative' , 'HIV-positive'};

% Calibration error bars
cin2_mean = cin2_2010_dObs(: , 2) .* 100;
cin2_sdev = (cin2_2010_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
errorbar(1 : length(cin2_mean) , cin2_mean , cin2_sdev , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(cin2_mean) , mean((cin2_2010 .* 100),1)' , 'k-' , ...
    1 : length(cin2_mean) , min((cin2_2010 .* 100),[],1)' , 'k--' , ...
    1 : length(cin2_mean) , max((cin2_2010 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
legend('(Sweet, 2020) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
set(gca , 'xtick' , 1 : length(diseaseGroup) , 'xtickLabel' , diseaseGroup);
xlabel('HIV Status'); ylabel('CIN2/3 Prevalence (%)')
title('CIN2/3 Prevalence in 2010 - Females aged 20-50')
ylim([0 25])
grid on;

saveas(gcf,strcat(groupDir,'/cin23ByHiv2007.fig'));
close(gcf);


%% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

%% Cervical cancer incidence in 2012 by age vs. Globocan 2012 data and other sources (calibration)
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

% Calibration error bars
meanObs = ccInc2012_dObs(: , 2);
sdevObs = (ccInc2012_dObs(: , 3).^(1/2)).*2;

figure;    
% Plot observed data
errorbar(4 : age , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
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
ylim([0 450]);
title(['Cervical Cancer Incidence in 2012']);
legend('(Globocan, 2012) Observed SA: mean, 2SD' , ...
    ['Model, general: ' , num2str(nRuns) , '-sets mean'] , ['Model: ' , num2str(nRuns) , '-sets minimum'] , ['Model: ' , num2str(nRuns) , '-sets maximum'] , ...
    ['Model, HIV-negative: ' , num2str(nRuns) , '-sets mean'] , ['Model: ' , num2str(nRuns) , '-sets minimum'] , ['Model: ' , num2str(nRuns) , '-sets maximum'] , ...
    ['Model, WLWHIV untreated: ' , num2str(nRuns) , '-sets mean'] , ['Model: ' , num2str(nRuns) , '-sets minimum'] , ['Model: ' , num2str(nRuns) , '-sets maximum'] , ...
    ['Model, WLWHIV on ART: ' , num2str(nRuns) , '-sets mean'] , ['Model: ' , num2str(nRuns) , '-sets minimum'] , ['Model: ' , num2str(nRuns) , '-sets maximum'] , ...
    'Location' , 'Northwest');
grid on;

saveas(gcf,strcat(groupDir,'/ccByAge.fig'));
close(gcf);


%% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************

%% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)

% HPV infected
% Calibration error bars
meanObs = hpv_dist_dObs(: , 2) .* 100;
sdevObs = (hpv_dist_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
subplot(2,3,1)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 46.82 , 'k*');
hold on;
plot(2012 , 53.18 , 'b*');
hold on;
plot(typeDistYearVec , mean((hpv_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((hpv_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((hpv_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
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
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 51.92 , 'k*');
hold on;
plot(2012 , 48.08 , 'b*');
hold on;
plot(typeDistYearVec , mean((cin1_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cin1_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cin1_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
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
hold on;
plot(typeDistYearVec , mean((cin2_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cin2_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cin2_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(typeDistYearVec , mean((cin2_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((cin2_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((cin2_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
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
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 73.71 , 'k*');
hold on;
plot(2012 , 26.29 , 'b*');
hold on;
plot(typeDistYearVec , mean((cin3_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cin3_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cin3_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
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
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 85.78 , 'k*');
hold on;
plot(2012 , 14.22 , 'b*');
hold on;
plot(typeDistYearVec , mean((cc_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cc_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cc_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(typeDistYearVec , mean((cc_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((cc_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((cc_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('Cervical Cancer')
legend('Observed 2012: mean, 2SD' , 'Observed 2012- 9v' , 'Observed 2012- non-9v' , ...
    ['Model- 9v: ' , num2str(nRuns) , '-sets mean'] , ['Model- 9v: ' , num2str(nRuns) , '-sets minimum'] , ['Model- 9v: ' , num2str(nRuns) , '-sets maximum'] , ...
    ['Model- non-9v: ' , num2str(nRuns) , '-sets mean'] , ['Model- non-9v: ' , num2str(nRuns) , '-sets minimum'] , ['Model- non-9v: ' , num2str(nRuns) , '-sets maximum']);
grid on;
saveas(gcf,strcat(groupDir,'/hpvTypeDistribution.fig'));
close(gcf);
