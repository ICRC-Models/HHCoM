% function [] = vaxCEA_multSims_CI()

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
    condUse , screenYrs , hpvScreenStartYear , waning , ...
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

lastYear = 2101;

% Indices of calib runs to plot
fileInds = {'11_1' , '11_2'}; % , '11_3' , '11_4' , '11_5' , '11_6' , ...
%     '11_7' , '11_8' , '11_9' , '11_10' , '11_11' , '11_12' , '11_13' , ...
%     '11_14' , '11_15' , '11_16' , '11_17' , '11_18' , '11_19' , '11_20' , ...
%     '11_21' , '11_22' , '11_23' , '11_24' , '11_25'};  % DO ART, 22Apr20Ph2V2, t=11
% fileInds = {'7_1' , '7_2' , '7_3' , '7_4' , '7_5' , '7_6' , '7_7' , '7_8' , ...
%     '7_9' , '7_10' , '7_11' , '7_12' , '7_13' , '7_14' , '7_15' , '7_16' , ...
%     '7_17' , '7_18' , '7_19' , '7_20' , '7_21' , '7_22' , '7_23' , '7_24' , '7_25'};  % DO ART, 22Apr20Ph2V2, t=7
nRuns = length(fileInds);

% Initialize model output plots
% Annual timespan
monthlyTimespan = [startYear : (1/6) : lastYear];
monthlyTimespan = monthlyTimespan(1 : end-1);
annualTimespan = [startYear : lastYear-1];
% Total population size
popSize = zeros(nRuns , length(monthlyTimespan));
% Population age distribution
popYearVec = [2018 2100];
popPropF = zeros(nRuns , length(popYearVec) , age);
% HIV prevalence
hivAgeM = zeros(nRuns , 7 , length(monthlyTimespan));
hivAgeF = hivAgeM;
hivPrev = zeros(nRuns , gender , length(monthlyTimespan));
% HIV deaths
hivDeathsM = zeros(nRuns , length(annualTimespan));
hivDeathsF = hivDeathsM;
% ART coverage
artCovM = zeros(nRuns , length(monthlyTimespan));
artCovF = artCovM;
artCovAge = zeros(nRuns , age , length(monthlyTimespan));
% Female HPV prevalence
hpv_hiv = zeros(nRuns , age , 2);
hpv_hivNeg = hpv_hiv;
hpv_hivTot = zeros(nRuns , age , 1);
% Male HPV prevalence
hpv_hivM = zeros(nRuns , 4 , 2);
hpv_hivMNeg = hpv_hivM;
hpv_hivMtot = zeros(nRuns , age , 1);
% HPV prevalence over time
hpv_hivTimeF = zeros(nRuns , length(monthlyTimespan));
hpv_hivNegTimeF = hpv_hivTimeF;
hpv_time = zeros(nRuns , 2 , length(monthlyTimespan));
% CIN2/3 prevalence
cinPosAge = zeros(nRuns , 2 , age);
cinNegAge = cinPosAge;
cinGenAge = cinPosAge;
cinGenTime = zeros(nRuns , length(monthlyTimespan));
cinPosTime = cinGenTime;
cinNegTime = cinGenTime;
% CC incidence
ccIncAge = zeros(nRuns , 3 , age);
ccIncTime = zeros(nRuns , length(annualTimespan));
ccIncTimeNeg = ccIncTime;
ccIncTimePos = ccIncTime;
ccIncTimeArt = ccIncTime;
% HPV/CIN/CC type distribution
cc_vax = zeros(nRuns , 4 , length(monthlyTimespan));
cc_nonVax = cc_vax;
cin23_vax = cc_vax;
cin23_nonVax = cc_vax;
cin3_vax = cc_vax;
cin3_nonVax = cc_vax;
cin2_vax = cc_vax;
cin2_nonVax = cc_vax;
cin1_vax = cc_vax;
cin1_nonVax = cc_vax;
hpv_vax = cc_vax;
hpv_nonVax = cc_vax;

resultsDir = [pwd , '\HHCoM_Results\'];
for j = 1 : nRuns
    % Load results
    pathModifier = ['22Apr20Ph2V2_noBaseVax_baseScreen_hpvHIVcalib_fertDec042-076-052_' , fileInds{j}]; % ***SET ME***: name for simulation output file
    nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
    curr = load([pwd , '/HHCoM_Results/toNow_22Apr20Ph2V2_noBaseVax_baseScreen_hpvHIVcalib_fertDec042-076_' , fileInds{j}]); % ***SET ME***: name for historical run output file 
    
    vaxResult = cell(nSims , 1);
    resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxSimResult'];
    if waning
        resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxWaneSimResult'];
    end
    for n = nSims
        % load results from vaccine run into cell array
        vaxResult{n} = load([resultFileName , num2str(3), '.mat']);
        % concatenate vectors/matrices of population up to current year to population
        % matrices for years past current year
        vaxResult{n}.popVec = [curr.popVec(1 : end  , :); vaxResult{n}.popVec(2 : end , :)];
        vaxResult{n}.ccDeath = [curr.ccDeath(1 : end , : , : , :) ; vaxResult{n}.ccDeath(2 : end , : , : , :)];
        vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :); vaxResult{n}.newCC(2 : end , : , : ,:)];
        vaxResult{n}.newScreen = [curr.newScreen(1 : end , : , : , : , : , : , : , : , :); vaxResult{n}.newScreen(2 : end , : , : , : , : , : , : , : , :)];
        vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , : , : , : , :); vaxResult{n}.newHiv(2 : end , : , : , : , : , : , :)];
        vaxResult{n}.hivDeaths = [curr.hivDeaths(1 : end , : , : , :); vaxResult{n}.hivDeaths(2 : end , : , : , :)];
        vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :); vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , :)];
        vaxResult{n}.tVec = [curr.tVec(1 : end), vaxResult{n}.tVec(2 : end)];
    end
    noVaxInd = nSims;
    noV = vaxResult{noVaxInd};
    tVec = noV.tVec;
    tVecYr = tVec(1 : stepsPerYear : end);

    %% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

    %% Population size over time vs. Statistics South Africa data (calibration)
    popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
    popSize(j , :) = sum(noV.popVec(: , popTot),2);
    
    %% Female population size by 5-year age groups over time vs. Statistics South Africa data (internal validation)
    for t = 1 : length(popYearVec)
        for a = 1 : age
            popAgeF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            popTotF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
            popPropF(j , t , a) = sum(noV.popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popAgeF),2) ./ ...
                sum(noV.popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTotF),2);
        end
    end

    %% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

    %% HIV prevalence by age over time vs. AHRI (calibration)  
    for a = 1 : age
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
        hivAgeM(j , a , :) =  (sum(noV.popVec(: , hivMInds) , 2) + sum(noV.popVec(: , artMInds) , 2)) ...
            ./ sum(noV.popVec(: , totMInds) , 2);
        hivAgeF(j , a , :) =  (sum(noV.popVec(: , hivFInds) , 2) + sum(noV.popVec(: , artFInds) , 2)) ...
            ./ sum(noV.popVec(: , totFInds) , 2);
    end

    %% HIV prevalence by gender over time vs. AHRI (validation) and (Vandormael, 2019) AHRI data (validation)
    for g = 1 : gender
            hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , g , 4 : 10 , 1 : risk));
            totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , g , 4 : 10 , 1 : risk));
            hivPrev(j , g , :) = (sum(noV.popVec(: , hivInds) , 2) ./ sum(noV.popVec(: , totInds) , 2));
    end
    
    %% HIV-associated deaths by gender over time
    hivDeathsM(j , :) = annlz(sum(sum(noV.hivDeaths(: , : , 1 , :), 2), 4));
    hivDeathsF(j , :) = annlz(sum(sum(noV.hivDeaths(: , : , 2 , :), 2), 4));
    
    %% Proportion of total HIV+ population on ART and VS (denominator: CD4-eligible and ineligible)
    artIndsF = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 4 : age , 1 : risk));
    hivAllIndsF = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , 2 , 4 : age , 1 : risk));
    artIndsM = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 , 4 : age , 1 : risk));
    hivAllIndsM = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
        1 : endpoints , 1 : intervens , 1 , 4 : age , 1 : risk));
    
    artCovF(j , :) = sum(noV.popVec(: , artIndsF) , 2) ./ (sum(noV.popVec(: , hivAllIndsF) , 2) + sum(noV.popVec(: , artIndsF) , 2));
    artCovM(j , :) = sum(noV.popVec(: , artIndsM) , 2) ./ (sum(noV.popVec(: , hivAllIndsM) , 2) + sum(noV.popVec(: , artIndsM) , 2));
    
    %% Proportion of total HIV+ population on ART and VS by age
    for a = 1 : age
        artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
        artPop = sum(noV.popVec(: , artInds) , 2);
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
            1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
        hivPop = sum(noV.popVec(: , hivInds) , 2);
        artCovAge(j , a , :) = artPop ./ hivPop;
    end

    
    %% ********************************** HPV FIGURES **********************************************************************************************

    %% Female HPV Prevalence by age and HIV status in 2002 and 2018 vs. McDonald 2014 data (calibration)
    hpvYearVec = [2002 2018];
    for i = 1 : length(hpvYearVec)
        yr = hpvYearVec(i);
        for a = 1 : age % 15-19 -> 55-65
            hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
                [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
            ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            hpv_hiv(j , a , i) = sum(noV.popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
                ./ sum(noV.popVec((yr - startYear) * stepsPerYear +1 , ageInds));

            hpvInds_hivNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
                [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
            ageInds_hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            hpv_hivNeg(j , a , i) = sum(noV.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivNeg))...
                ./ sum(noV.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivNeg));
        end
    end
    
    %% Female HPV Prevalence by age 2018
    hpvYearVec = [2018];
    for i = 1 : length(hpvYearVec)
        yr = hpvYearVec(i);
        for a = 1 : age
            hpvInds_hivTot = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
                [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
            ageInds_hivTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            hpv_hivTot(j , a , i) = sum(noV.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivTot))...
                ./ sum(noV.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivTot));
        end
    end

    %% Male HPV prevalence by age and HIV status in 2008 vs. Mbulawa data (calibration)
    hpvYearVec = [2008 2018];
    ageVec = {[4:5],[6:7],[8:9],[10:13]};
    for i = 1 : length(hpvYearVec)
        yr = hpvYearVec(i);
        for aV = 1 : length(ageVec)
            a = ageVec{aV};
            hpvInds_hivM = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 , [1 : 2 , 7] , ...
                1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
                [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
            ageInds_hivM = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            hpv_hivM(j , aV , i) = sum(noV.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivM))...
                / sum(noV.popVec((yr - startYear) * stepsPerYear+1 , ageInds_hivM));

            hpvInds_hivMNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 , [1 : 2 , 7] , ...
                1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
                [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
            ageInds_hivMNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            hpv_hivMNeg(j , aV , i) = sum(noV.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivMNeg))...
                / sum(noV.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivMNeg));
        end
    end
    
    %% Female HPV Prevalence over time by HIV status
    hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
    hpv_hivTimeF(j , :) = sum(noV.popVec(: , hpvInds) , 2) ./ sum(noV.popVec(: , ageInds) , 2);

    hpvInds_hivNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    ageInds_hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
    hpv_hivNegTimeF(j , :) = sum(noV.popVec(: , hpvInds_hivNeg) , 2) ./ sum(noV.popVec(: , ageInds_hivNeg) , 2);

    %% HPV Prevalence over time by sex       
    for g = 1 : gender
        hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , g , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , g , 1 : age , 1 : risk))]);
        popInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , 1 : age , 1 : risk));
        hpv_time(j , g , :) = sum(noV.popVec(: , hpvInds) , 2)...
            ./ sum(noV.popVec(: , popInds) , 2);
    end


    %% ********************************** CIN FIGURES *********************************************************************************************

    %% CIN2/3 prevalence for All HR HPV types combined by HIV status and age in 2002 vs. McDonald 2014 data (calibration)
    cinYearVec = [2002 2018];
    for i = 1 : length(cinYearVec)
        yr = cinYearVec(i);
        for a = 1 : age % date is 15-19 -> 60-64
            % HIV-positive (on and not on ART)
            cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
                [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
            ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            cinPosAge(j , i , a) = (sum(noV.popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
                ./ sum(noV.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
            % HIV-negative
            cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
                [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
            ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            cinNegAge(j , i , a) = (sum(noV.popVec((yr - startYear) * stepsPerYear +1 , cinNegInds)))...
                ./ (sum(noV.popVec((yr - startYear) * stepsPerYear +1 , ageNegInds)));
            % General
            cinGenInds = unique([toInd(allcomb(1 : disease , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
                [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
            ageGenInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            cinGenAge(j , i , a) = (sum(noV.popVec((yr - startYear) * stepsPerYear +1 , cinGenInds)))...
                ./ (sum(noV.popVec((yr - startYear) * stepsPerYear +1 , ageGenInds)));
        end
    end
    
    %% CIN2/3 prevalence for All HR HPV types combined over time
    % General
    cinGenInds = unique([toInd(allcomb(1 : disease , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    ageGenInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
    cinGenTime(j , :) = sum(noV.popVec(: , cinGenInds) , 2)...
        ./ sum(noV.popVec(: , ageGenInds) , 2);
    % HIV-positive (on and not on ART)
    cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
    cinPosTime(j , :) = sum(noV.popVec(: , cinInds) , 2)...
        ./ sum(noV.popVec(: , ageInds) , 2);
    % HIV-negative
    cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
    cinNegTime(j , :) = sum(noV.popVec(: , cinNegInds) , 2)...
        ./ sum(noV.popVec(: , ageNegInds) , 2);

    
    %% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

    %% Cervical cancer incidence in 2012 by age vs. Globocan 2012 data and other sources (calibration)
    fac = 10 ^ 5;
    ccYearVec = [2005 2012 2018];
    for i = 1 : length(ccYearVec)
        yr = ccYearVec(i);
        incTimeSpan = [((yr - startYear) * stepsPerYear +1) : ((yr - startYear) * stepsPerYear +6)];
        for a = 1 : age
            allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            % Calculate incidence
            ccIncAge(j , i , a) = ...
                (annlz(sum(sum(sum(noV.newCC(incTimeSpan , : , a , :),2),3),4)) ./ ...
                (annlz(sum(noV.popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);
        end
    end
    
    %% Cervical cancer incidence over time
    fac = 10 ^ 5;
    % General population
    allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
    % Calculate incidence
    ccIncTime(j , :) = ...
        (annlz(sum(sum(sum(noV.newCC(: , : , 1 : age , :),2),3),4)) ./ ...
        (annlz(sum(noV.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);

    % HIV-negative
    allFneg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
    % Calculate incidence
    ccIncTimeNeg(j , :) = ...
        (annlz(sum(sum(sum(noV.newCC(: , 1 : 2 , 1 : age , :),2),3),4)) ./ ...
        (annlz(sum(noV.popVec(: , allFneg) , 2) ./ stepsPerYear)) * fac);

    % HIV-positive untreated
    allFpos = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
    % Calculate incidence
    ccIncTimePos(j , :) = ...
        (annlz(sum(sum(sum(noV.newCC(: , 3 : 7 , 1 : age , :),2),3),4)) ./ ...
        (annlz(sum(noV.popVec(: , allFpos) , 2) ./ stepsPerYear)) * fac);

    % HIV-positive on ART
    allFart = toInd(allcomb(8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
    % Calculate incidence
    ccIncTimeArt(j , :) = ...
        (annlz(sum(sum(sum(noV.newCC(: , 8 , 1 : age , :),2),3),4)) ./ ...
        (annlz(sum(noV.popVec(: , allFart) , 2) ./ stepsPerYear)) * fac);
    
   
    %% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************
    
    %% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)
    diseaseInds = {[1 : disease] , [1 : 2] , [3 : 7] , 8};
    for dInd = 1 : length(diseaseInds)
        d = diseaseInds{dInd};
        ccInds_vax = toInd(allcomb(d , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
            1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
        ccInds_nonVax = toInd(allcomb(d , 1 : viral , [1 : 5 , 7] , 6 , ...
            1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
        ccInds_tot = unique([toInd(allcomb(d , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
                1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

       cin23Inds_vax = [toInd(allcomb(d , 1 : viral , 4 , [1 : 4 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(d , 1 : viral , 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk))];
        cin23Inds_nonVax = [toInd(allcomb(d , 1 : viral , [1 : 3 , 7] , 4 , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(d , 1 : viral , [1 : 4 , 7] , 5 , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk))];
        cin23Inds_tot = unique([toInd(allcomb(d , 1 : viral , 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(d , 1 : viral , 4 , [1 : 4 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); 
            toInd(allcomb(d , 1 : viral , ...
            [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(d , 1 : viral , ...
            [1 3 , 7] , 4 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);     

        cin3Inds_vax = toInd(allcomb(d , 1 : viral , 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk));
        cin3Inds_nonVax = toInd(allcomb(d , 1 : viral , [1 : 4 , 7] , 5 , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk));
        cin3Inds_tot = unique([toInd(allcomb(d , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

        cin2Inds_vax = toInd(allcomb(d , 1 : viral , 4 , [1 : 4 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk));
        cin2Inds_nonVax = toInd(allcomb(d , 1 : viral , [1 : 3 , 7] , 4 , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk));
        cin2Inds_tot = unique([toInd(allcomb(d , 1 : viral , 4 , [1 : 4 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(d, 1 : viral , ...
                [1 : 3 , 7] , 4 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

        cin1Inds_vax = toInd(allcomb(d , 1 : viral , 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk));
        cin1Inds_nonVax = toInd(allcomb(d , 1 : viral , [1 : 2 , 7] , 3 , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk));
        cin1Inds_tot = unique([toInd(allcomb(d , 1 : viral , 3 , [1 : 3 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

        hpvInds_vax = toInd(allcomb(d , 1 : viral , 2 , [1 : 2 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk));
        hpvInds_nonVax = toInd(allcomb(d , 1 : viral , [1 , 7] , 2 , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk));
        hpvInds_tot = unique([toInd(allcomb(d , 1 : viral , 2 , [1 : 2 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                [1 , 7] , 2 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

        cc_vax(j , dInd , :) = sum(noV.popVec(: , ccInds_vax) , 2)...
            ./ sum(noV.popVec(: , ccInds_tot) , 2);
        cc_nonVax(j , dInd , :) = sum(noV.popVec(: , ccInds_nonVax) , 2)...
            ./ sum(noV.popVec(: , ccInds_tot) , 2);

        cin23_vax(j , dInd , :) = sum(noV.popVec(: , cin23Inds_vax) , 2)...
            ./ sum(noV.popVec(: , cin23Inds_tot) , 2);
        cin23_nonVax(j , dInd , :) = sum(noV.popVec(: , cin23Inds_nonVax) , 2)...
            ./ sum(noV.popVec(: , cin23Inds_tot) , 2);

        cin3_vax(j , dInd , :) = sum(noV.popVec(: , cin3Inds_vax) , 2)...
            ./ sum(noV.popVec(: , cin3Inds_tot) , 2);
        cin3_nonVax(j , dInd , :) = sum(noV.popVec(: , cin3Inds_nonVax) , 2)...
            ./ sum(noV.popVec(: , cin3Inds_tot) , 2);

        cin2_vax(j , dInd , :) = sum(noV.popVec(: , cin2Inds_vax) , 2)...
            ./ sum(noV.popVec(: , cin2Inds_tot) , 2);
        cin2_nonVax(j , dInd , :) = sum(noV.popVec(: , cin2Inds_nonVax) , 2)...
            ./ sum(noV.popVec(: , cin2Inds_tot) , 2);

        cin1_vax(j , dInd , :) = sum(noV.popVec(: , cin1Inds_vax) , 2)...
            ./ sum(noV.popVec(: , cin1Inds_tot) , 2);
        cin1_nonVax(j , dInd , :) = sum(noV.popVec(: , cin1Inds_nonVax) , 2)...
            ./ sum(noV.popVec(: , cin1Inds_tot) , 2);

        hpv_vax(j , dInd , :) = sum(noV.popVec(: , hpvInds_vax) , 2)...
            ./ sum(noV.popVec(: , hpvInds_tot) , 2);
        hpv_nonVax(j , dInd , :) = sum(noV.popVec(: , hpvInds_nonVax) , 2)...
            ./ sum(noV.popVec(: , hpvInds_tot) , 2);
    end

end

%% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

%% Population size over time vs. Statistics South Africa data (calibration/validation)
% Load calibration/validation data from Excel (years, values)
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
    'rs' , 'LineWidth' , 1.5);
hold all;
plot(futurePop0_69(:,1) , futurePop0_69(:,2) , 'bo');
plot(monthlyTimespan , mean(popSize,1)' , 'k-' , ...
    monthlyTimespan , min(popSize,[],1)' , 'k--' , ...
    monthlyTimespan , max(popSize,[],1)' , 'k--' , 'LineWidth' , 1.5);
title('KZN Population Size, Ages 0-79')
xlabel('Year'); ylabel('Individuals')
xlim([1980 2120]); ylim([0 (20*10^6)]);
legend('(Statistics SA) Observed KZN, ages 0-79: mean, 2SD' , ...
    '(UN World Population Prospects) Future SA adjusted for KZN, ages 0-69' , ...
    'Model, ages 0-79: 25-sets mean' , 'Model, ages 0-79: 25-sets minimum' , 'Model, ages 0-79: 25-sets maximum');
grid on;

%% Female population size by 5-year age groups over time vs. Statistics South Africa data (internal validation)
% Load calibration data from Excel
file = [pwd , '/Config/Population_validation_targets.xlsx'];
kzn_popByage_yrs(: , 1) = xlsread(file , 'Demographics' , 'F112:F127').*1000;    % females by age in 2019
popPropF_obs = zeros(age , 1);
for a = 1 : age
    popPropF_obs(a , 1) = sum(kzn_popByage_yrs(a , 1)) / sumall(kzn_popByage_yrs(1 : end , 1));
end

ageGroup = {'0-4' , '5-9' , '10-14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
figure;
plot(1 : length(popPropF_obs) , popPropF_obs , 'ro')
hold all;
plot(1 : length(popPropF) , squeeze(mean(popPropF(: , 1 , :),1)) , 'k-' , 'LineWidth' , 1.5);
hold all;
plot(1 : length(popPropF) , squeeze(min(popPropF(: , 1 , :),[],1)) , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(1 : length(popPropF) , squeeze(max(popPropF(: , 1 , :),[],1)) , 'k--' , 'LineWidth' , 1.5);
hold all;

plot(1 : length(popPropF) , squeeze(mean(popPropF(: , 2 , :),1)) , 'b-' , 'LineWidth' , 1.5);
hold all;
plot(1 : length(popPropF) , squeeze(min(popPropF(: , 2 , :),[],1)) , 'b--' , 'LineWidth' , 1.5);
hold all;
plot(1 : length(popPropF) , squeeze(max(popPropF(: , 2 , :),[],1)) , 'b--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Population proportion'); title('Female Age Distribution');
ylim([0 0.2]); grid on;
legend('(Statistics SA) Observed KZN, 2019' , 'Model, 2019: 25-sets mean' , ...
    'Model, 2019: 25-sets minimum' , 'Model, 2019: 25-sets maximum' , ...
    'Model, 2100: 25-sets mean' , 'Model, 2100: 25-sets minimum' , ...
    'Model, 2100: 25-sets maximum' , 'Location' , 'northeast');

%% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

%% HIV prevalence by age over time vs. AHRI (calibration, validation) and IHME model data (validation)
% Calibration error bars
hivM(: , 1) = hivPrevM_dObs(: , 2); % mean
hivM(: , 2) = (hivPrevM_dObs(: , 3).^(1/2)).*2; % calibration SD
hivF(: , 1) = hivPrevF_dObs(: , 2); % mean
hivF(: , 2) = (hivPrevF_dObs(: , 3).^(1/2)).*2; % calibration SD

ageGroup = {'15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74'};
gen = {'Male' , 'Female'};
for g = 1 : gender
    hivPrevs = hivM;
    hivModel = hivAgeM;
    if g == 2
        hivPrevs = hivF;
        hivModel = hivAgeF;
    end
    figure; 
    for a = 4 : 10
        subplot(3 , 3 , a-3)
        hold all;
        if a <= 11            
            errorbar(unique(hivPrevM_dObs(: ,1)) , hivPrevs(((a-3) : 7 : end) , 1) , hivPrevs(((a-3) : 7 : end) , 2) , ...
                'rs' , 'LineWidth' , 1.5);
            hold on;
            plot(monthlyTimespan , mean((squeeze(hivModel(: , a , :))),1)' , 'k-' , ...
                monthlyTimespan , min((squeeze(hivModel(: , a , :))),[],1)' , 'k--' , ...
                monthlyTimespan , max((squeeze(hivModel(: , a , :))),[],1)' , 'k--' , 'LineWidth' , 1.5);
        else
            plot(monthlyTimespan , mean((squeeze(hivModel(: , a , :))),1)' , 'k-' , ...
                monthlyTimespan , min((squeeze(hivModel(: , a , :))),[],1)' , 'k--' , ...
                monthlyTimespan , max((squeeze(hivModel(: , a , :))),[],1)' , 'k--' , 'LineWidth' , 1.5);
        end
        xlabel('Year'); ylabel('HIV Prevalence'); title([gen{g} , 's ages ' , ageGroup{a-3}])
        xlim([1985 2030]); 
        if g == 1 
            ylim([0 0.6]);
        elseif g == 2
            ylim([0 0.8]);
        end
        grid on;
    end
    legend('(AHRI data request) Observed KZN: mean, 2SD' , ...
        'Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum')
end

%% HIV prevalence by gender over time vs. AHRI (calibration)
hivData(: , : , 1) = zeros(length(unique(hivPrevM_dObs(: ,1))) , 1);
hivData(: , : , 2) = zeros(length(unique(hivPrevM_dObs(: ,1))) , 1);
hivRaw(:,:,1) = hivPrevM_dObs(: , 4:5);
hivRaw(:,:,2) = hivPrevF_dObs(: , 4:5);
for i = 1 : length(unique(hivPrevM_dObs(: ,1)))
    for g = 1 : gender
        hivData(i,1,g) = (sumall(hivRaw(((i-1)*7+1):(i*7) , 1 , g)) ./ sumall(hivRaw(((i-1)*7+1):(i*7) , 2 , g)));
    end
end

figure;
gen = {'Males aged 15-49' , 'Female aged 15-49'};
for g = 1 : gender
    subplot(1,2,g)
    plot(unique(hivPrevM_dObs(: ,1)) , hivData(:,:,g) , 'ro');
    hold on;
    plot(monthlyTimespan , mean(squeeze(hivPrev(: , g , :)),1)' , 'k-' , ...
        monthlyTimespan , min(squeeze(hivPrev(: , g , :)),[],1)' , 'k--' , ...
        monthlyTimespan , max(squeeze(hivPrev(: , g , :)),[],1)' , 'k--' , 'LineWidth' , 1.5);
    xlabel('Year'); ylabel('HIV Prevalence'); title(gen{g});
    xlim([1985 2030]); ylim([0 0.5]);
    if g == 1 
        legend('(AHRI data request) Observed KZN, ages 15-49' , ...
            'Model, ages 15-49: 25-sets mean' , 'Model, ages 15-49: 25-sets minimum' , 'Model, ages 15-49: 25-sets maximum');
    elseif g == 2
        legend('(AHRI data request) Observed KZN, ages 15-49' , ...
            'Model, ages 15-49: 25-sets mean' , 'Model, ages 15-49: 25-sets minimum' , 'Model, ages 15-49: 25-sets maximum');
    end
end

%% HIV Prevalence by age in 2009, 2018 vs. AHRI (calibration)
ageGroup = {'0-4' , '5-9' , '10-14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
% Calibration error bars
hivM_2009(: , 1) = hivPrevM_dObs(end-6:end , 2); % mean
hivM_2009(: , 2) = (hivPrevM_dObs(end-6:end , 3).^(1/2)).*2; % calibration SD
hivF_2009(: , 1) = hivPrevF_dObs(end-6:end , 2); % mean
hivF_2009(: , 2) = (hivPrevF_dObs(end-6:end , 3).^(1/2)).*2; % calibration SD

figure;
gen = {'Males' , 'Females'};
for g = 1 : gender
    hivPrevs = hivM_2009;
    hivModel = hivAgeM;
    if g == 2
        hivPrevs = hivF_2009;
        hivModel = hivAgeF;
    end
    subplot(2 , 1 , g)
    hold all;            
    errorbar([4 : 10] , hivPrevs(: , 1)' , hivPrevs(: , 2)' , 'rs' , 'LineWidth' , 1.5);
    hold all;
    plot(1 : length(ageGroup) , mean((squeeze(hivModel(: , : , ((2009 - startYear) * stepsPerYear +1)))),1) , 'k-' , ...
        1 : length(ageGroup) , min((squeeze(hivModel(: , : , ((2009 - startYear) * stepsPerYear +1)))),[],1) , 'k--' , ...
        1 : length(ageGroup) , max((squeeze(hivModel(: , : , ((2009 - startYear) * stepsPerYear +1)))),[],1) , 'k--' , 'LineWidth' , 1.5);
    hold all;
    plot(1 : length(ageGroup) , mean((squeeze(hivModel(: , : , ((2018 - startYear) * stepsPerYear +1)))),1) , 'b-' , ...
        1 : length(ageGroup) , min((squeeze(hivModel(: , : , ((2018 - startYear) * stepsPerYear +1)))),[],1) , 'b--' , ...
        1 : length(ageGroup) , max((squeeze(hivModel(: , : , ((2018 - startYear) * stepsPerYear +1)))),[],1) , 'b--' , 'LineWidth' , 1.5);
    hold all;
    set(gca , 'xtickLabel' , ageGroup);
    set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
    xlabel('Age Group'); ylabel('HIV Prevalence'); title(gen{g});
    if g == 1 
        ylim([0 0.8]);
    elseif g == 2
        ylim([0 0.8]);
    end
    grid on;
    legend('(AHRI data request) Observed KZN, 2009: mean, 2SD' , ...
        'Model, 2009: 25-sets mean' , 'Model, 2009: 25-sets minimum' , 'Model, 2009: 25-sets maximum' , ...
        'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum')
end

%% HIV-associated deaths by gender over time
figure;
plot(annualTimespan , mean(hivDeathsM(: , :),1)' , 'k-' , ...
    annualTimespan , min(hivDeathsM(: , :),[],1)' , 'k--' , ...
    annualTimespan , max(hivDeathsM(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(annualTimespan , mean(hivDeathsF(: , :),1)' , 'b-' , ...
    annualTimespan , min(hivDeathsF(: , :),[],1)' , 'b--' , ...
    annualTimespan , max(hivDeathsF(: , :),[],1)' , 'b--' , 'LineWidth' , 1.5);
xlabel('Year'); ylabel('HIV-associated deaths'); title('HIV-associated deaths over time'); grid on;
xlim([1980 2030]);
legend('Model, males aged 0-79: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, females aged 0-79: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');

%% Proportion of total HIV+ population on ART and VS (denominator: CD4-eligible and ineligible)
figure;
subplot(1,2,1);
plot((artYr + 1) , maxRateM , 'ro');
hold all;
plot(monthlyTimespan , mean(artCovM,1)' , 'k-' , ...
    monthlyTimespan , min(artCovM,[],1)' , 'k--' , ...
    monthlyTimespan , max(artCovM,[],1)' , 'k--' , 'LineWidth' , 1.5);
hold all;
xlim([1985 2030]); ylim([0 1]);
xlabel('Year'); ylabel('Proportion MLWHIV on ART + VS')
title('Males aged 15-79'); grid on;
legend('Observed KZN' , 'Model, ages 15-79: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');

subplot(1,2,2);
plot((artYr + 1) , maxRateF , 'ro');
hold all;
plot(monthlyTimespan , mean(artCovF,1)' , 'k-' , ...
    monthlyTimespan , min(artCovF,[],1)' , 'k--' , ...
    monthlyTimespan , max(artCovF,[],1)' , 'k--' , 'LineWidth' , 1.5);
xlim([1985 2030]); ylim([0 1]);
xlabel('Year'); ylabel('Proportion WLWHIV on ART + VS')
title('Females aged 15-79'); grid on;
legend('Model, ages 15-79: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');

%% Proportion of total HIV+ population on ART and VS by age
figure;
ageGroup = {'0-4' , '5-9' ,'10-14' , '15-19' , '20-24' , '25-29' ,...
     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
     '60-64' , '65-69' , '70-74' , '75-79'};
plot([1:age] , mean(artCovAge(:,:,((2018 - startYear) * stepsPerYear +1)),1) , 'k-' , ...
    [1:age] , min(artCovAge(:,:,((2018 - startYear) * stepsPerYear +1)),[],1) , 'k--' , ...
    [1:age] , max(artCovAge(:,:,((2018 - startYear) * stepsPerYear +1)),[],1) , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age'); ylabel('Proportion PLWHIV on ART + VS');
ylim([0 1]); grid on;
legend('Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');

%% ********************************** HPV FIGURES **********************************************************************************************

%% HPV Prevalence by age in 2002 and 2018 vs. McDonald 2014 data (calibration)
ageGroup = {'0 - 4' , '5 - 9' , '10 - 14' , '*17* - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};

% Calibration error bars
meanObs = hpv_hiv_dObs(: , 2);
sdevObs = (hpv_hiv_dObs(: , 3).^(1/2)).*2;
meanNeg = hpv_hivNeg_dObs(: , 2);
sdevNeg = (hpv_hivNeg_dObs(: , 3).^(1/2)).*2;

figure;
subplot(3,1,1);
errorbar(4 : length(meanObs)+4-1 , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , mean(hpv_hiv(: , : , 1),1)' , 'k-' , ...
    1 : age , min(hpv_hiv(: , : , 1),[],1)' , 'k--' , ...
    1 : age , max(hpv_hiv(: , : , 1),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , mean(hpv_hiv(: , : , 2),1)' , 'b-' , ...
    1 : age , min(hpv_hiv(: , : , 2),[],1)' , 'b--' , ...
    1 : age , max(hpv_hiv(: , : , 2),[],1)' , 'b--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence');
ylim([0 1]);
legend('(McDonald, 2014) Observed Cape Town: mean, 2SD' , 'Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');
title('hrHPV Prevalence - Females, HIV+ (includes CIN)'); grid on;

subplot(3,1,2);
errorbar(4 : length(meanObs)+4-1 , meanNeg , sdevNeg , ...
    'rs' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , mean(hpv_hivNeg(: , : , 1),1)' , 'k-' , ...
    1 : age , min(hpv_hivNeg(: , : , 1),[],1)' , 'k--' , ...
    1 : age , max(hpv_hivNeg(: , : , 1),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , mean(hpv_hivNeg(: , : , 2),1)' , 'b-' , ...
    1 : age , min(hpv_hivNeg(: , : , 2),[],1)' , 'b--' , ...
    1 : age , max(hpv_hivNeg(: , : , 2),[],1)' , 'b--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence');
ylim([0 1]);
legend('(McDonald, 2014) Observed Cape Town: mean, 2SD' , ...
    'Model, 2002: 25-sets mean' , 'Model, 2002: 25-sets minimum' , 'Model, 2002: 25-sets maximum' , ...
    'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');
title('hrHPV Prevalence - Females, HIV- (includes CIN)');
grid on;

subplot(3,1,3);
ageGroup = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 - 34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
plot(1 : age , mean(hpv_hivTot(: , : , 1),1)' , 'b-' , ...
    1 : age , min(hpv_hivTot(: , : , 1),[],1)' , 'b--' , ...
    1 : age , max(hpv_hivTot(: , : , 1),[],1)' , 'b--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence');
ylim([0 1]);
legend('Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');
title('hrHPV Prevalence - Females, All');
grid on;

%% HPV prevalence by age and HIV status in 2008 vs. Mbulawa data (calibration)
ageGroup = {'15-24' , '25-34' , '35-44' , '45-64'};

% Calibration error bars
meanObs = hpv_hivM2008_dObs(: , 2);
sdevObs = (hpv_hivM2008_dObs(: , 3).^(1/2)).*2;
meanNeg = hpv_hivMNeg2008_dObs(: , 2);
sdevNeg = (hpv_hivMNeg2008_dObs(: , 3).^(1/2)).*2;

figure;
subplot(3,1,1)
errorbar(1 : length(meanObs) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold all;
plot(1 : length(meanObs) , mean(hpv_hivM(: , : , 1),1) , 'k-' , ...
    1 : length(meanObs) , min(hpv_hivM(: , : , 1),[],1) , 'k--' , ...
    1 : length(meanObs) , max(hpv_hivM(: , : , 1),[],1) , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(1 : length(meanObs) , mean(hpv_hivM(: , : , 2),1) , 'b-' , ...
    1 : length(meanObs) , min(hpv_hivM(: , : , 2),[],1) , 'b--' , ...
    1 : length(meanObs) , max(hpv_hivM(: , : , 2),[],1) , 'b--' , 'LineWidth' , 1.5);
set(gca , 'xtick' , [1 : length(ageGroup)] , 'xtickLabel' , ageGroup);
legend('(Mbulawa, 2015) Observed SA: mean, 2SD' , ...
    'Model, 2008: 25-sets mean' , 'Model, 2008: 25-sets minimum' , 'Model, 2008: 25-sets maximum' , ...
    'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');
xlabel('Age Group'); ylabel('hrHPV Prevalence'); ylim([0 1]);
title('hrHPV Prevalence - Males, HIV+');
grid on;

subplot(3,1,2)
errorbar(1 : length(meanObs) , meanNeg , sdevNeg , ...
    'rs' , 'LineWidth' , 1.5);
hold all;
plot(1 : length(meanObs) , mean(hpv_hivMNeg(: , : , 1),1) , 'k-' , ...
    1 : length(meanObs) , min(hpv_hivMNeg(: , : , 1),[],1) , 'k--' , ...
    1 : length(meanObs) , max(hpv_hivMNeg(: , : , 1),[],1) , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(1 : length(meanObs) , mean(hpv_hivMNeg(: , : , 2),1) , 'b-' , ...
    1 : length(meanObs) , min(hpv_hivMNeg(: , : , 2),[],1) , 'b--' , ...
    1 : length(meanObs) , max(hpv_hivMNeg(: , : , 2),[],1) , 'b--' , 'LineWidth' , 1.5);
set(gca , 'xtick' , [1 : length(ageGroup)] , 'xtickLabel' , ageGroup);
legend('(Mbulawa, 2015) Observed SA: mean, 2SD' , ...
    'Model, 2008: 25-sets mean' , 'Model, 2008: 25-sets minimum' , 'Model, 2008: 25-sets maximum' , ...
    'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');
xlabel('Age Group'); ylabel('hrHPV Prevalence'); ylim([0 1]);
title('hrHPV Prevalence - Males, HIV-');
grid on;

subplot(3,1,3);
plot(1 : length(meanObs) , mean(hpv_hivMtot(: , : , 1),1) , 'b-' , ...
    1 : length(meanObs) , min(hpv_hivMtot(: , : , 1),[],1) , 'b--' , ...
    1 : length(meanObs) , max(hpv_hivMtot(: , : , 1),[],1) , 'b--' , 'LineWidth' , 1.5);
set(gca , 'xtick' , [1 : length(ageGroup)] , 'xtickLabel' , ageGroup);
legend('Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');
xlabel('Age Group'); ylabel('hrHPV Prevalence'); ylim([0 1]);
title('hrHPV Prevalence - Males, HIV-');
grid on;

%% Female HPV Prevalence over time by HIV status
figure;
subplot(1,2,2);
plot(monthlyTimespan , mean(hpv_hivTimeF(: , :),1)' , 'k-' , ...
    monthlyTimespan , min(hpv_hivTimeF(: , :),[],1)' , 'k--' , ...
    monthlyTimespan , max(hpv_hivTimeF(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
    xlabel('Year'); ylabel('HPV Prevalence'); title('WLWHIV ages 0-79  (includes CIN)');
    xlim([1985 2030]); ylim([0 1]); grid on;
    legend('Model, ages 0-79: 25-sets mean' , 'Model, ages 0-79: 25-sets minimum' , ...
        'Model, ages 0-79: 25-sets maximum');
    
subplot(1,2,1);
plot(monthlyTimespan , mean(hpv_hivNegTimeF(: , :),1)' , 'k-' , ...
    monthlyTimespan , min(hpv_hivNegTimeF(: , :),[],1)' , 'k--' , ...
    monthlyTimespan , max(hpv_hivNegTimeF(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
xlabel('Year'); ylabel('HPV Prevalence'); title('HIV-negative women ages 0-79 (includes CIN)');
xlim([1985 2030]); ylim([0 1]); grid on;
legend('Model, ages 0-79: 25-sets mean' , 'Model, ages 0-79: 25-sets minimum' , ...
    'Model, ages 0-79: 25-sets maximum');

%% HPV Prevalence over time by sex    
figure;
gen = {'Males aged 0-79' , 'Females aged 0-79 (includes CIN)'};
for g = 1 : gender
    subplot(1,2,g)
    plot(monthlyTimespan , mean(squeeze(hpv_time(: , g , :)),1)' , 'k-' , ...
        monthlyTimespan , min(squeeze(hpv_time(: , g , :)),[],1)' , 'k--' , ...
        monthlyTimespan , max(squeeze(hpv_time(: , g , :)),[],1)' , 'k--' , 'LineWidth' , 1.5);
    xlabel('Year'); ylabel('HPV Prevalence'); title(gen{g});
    xlim([1985 2030]); ylim([0 1]); grid on;
    legend('Model, ages 0-79: 25-sets mean' , 'Model, ages 0-79: 25-sets minimum' , 'Model, ages 0-79: 25-sets maximum');
end

%% ********************************** CIN FIGURES *********************************************************************************************

%% CIN2/3 prevalence for All HR HPV types combined by HIV status and age in 2002 vs. McDonald 2014 data (calibration)
ageGroup = {'0-4' , '5-9' , '10-14' , '*17*-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

% Calibration error bars
cinPos_mean = cinPos2002_dObs(: , 2);
cinPos_sdev = (cinPos2002_dObs(: , 3).^(1/2)).*2;
cinNeg_mean = cinNeg2002_dObs(: , 2);
cinNeg_sdev = (cinNeg2002_dObs(: , 3).^(1/2)).*2;

figure;
subplot(3 , 1 , 1);
errorbar(4 : length(cinPos_mean)+4-1 , cinPos_mean' , cinPos_sdev' , ...
    'rs' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , mean(squeeze(cinPosAge(:,1,:)),1)' , 'k-' , ...
    1 : age , min(squeeze(cinPosAge(:,1,:)),[],1)' , 'k--' , ...
    1 : age , max(squeeze(cinPosAge(:,1,:)),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , mean(squeeze(cinPosAge(:,2,:)),1)' , 'b-' , ...
    1 : age , min(squeeze(cinPosAge(:,2,:)),[],1)' , 'b--' , ...
    1 : age , max(squeeze(cinPosAge(:,2,:)),[],1)' , 'b--' , 'LineWidth' , 1.5);
legend('(McDonald, 2014) Observed Cape Town: mean, 2SD' , ...
    'Model, 2002: 25-sets mean' , 'Model, 2002: 25-sets minimum' , 'Model, 2002: 25-sets maximum' , ...
    'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence')
title('CIN 2/3 Prevalence - WLWHIV')
ylim([0 0.25])
grid on;

subplot(3 , 1 , 2)
errorbar(4 : length(cinNeg_mean)+4-1 , cinNeg_mean' , cinNeg_sdev' , ...
    'rs' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , mean(squeeze(cinNegAge(:,1,:)),1)' , 'k-' , ...
    1 : age , min(squeeze(cinNegAge(:,1,:)),[],1)' , 'k--' , ...
    1 : age , max(squeeze(cinNegAge(:,1,:)),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , mean(squeeze(cinNegAge(:,2,:)),1)' , 'b-' , ...
    1 : age , min(squeeze(cinNegAge(:,2,:)),[],1)' , 'b--' , ...
    1 : age , max(squeeze(cinNegAge(:,2,:)),[],1)' , 'b--' , 'LineWidth' , 1.5);
hold all;
legend('(McDonald, 2014) Observed Cape Town: mean, 2SD' , ...
    'Model, 2002: 25-sets mean' , 'Model, 2002: 25-sets minimum' , 'Model, 2002: 25-sets maximum' , ...
    'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence')
title('CIN 2/3 Prevalence - HIV-negative women')
ylim([0 0.1])
grid on;

subplot(3 , 1 , 3);
plot(1 : age , mean(squeeze(cinGenAge(:,1,:)),1)' , 'k-' , ...
    1 : age , min(squeeze(cinGenAge(:,1,:)),[],1)' , 'k--' , ...
    1 : age , max(squeeze(cinGenAge(:,1,:)),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , mean(squeeze(cinGenAge(:,2,:)),1)' , 'b-' , ...
    1 : age , min(squeeze(cinGenAge(:,2,:)),[],1)' , 'b--' , ...
    1 : age , max(squeeze(cinGenAge(:,2,:)),[],1)' , 'b--' , 'LineWidth' , 1.5);
legend('Model, 2002: 25-sets mean' , 'Model, 2002: 25-sets minimum' , 'Model, 2002: 25-sets maximum' , ...
    'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence')
title('CIN 2/3 Prevalence - General')
ylim([0 0.1])
grid on;

%% CIN2/3 prevalence for All HR HPV types combined by HIV status over time
figure;
subplot(1,2,2);
plot(monthlyTimespan , mean(cinPosTime(: , :),1)' , 'k-' , ...
    monthlyTimespan , min(cinPosTime(: , :),[],1)' , 'k--' , ...
    monthlyTimespan , max(cinPosTime(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
    xlabel('Year'); ylabel('CIN 2/3 Prevalence'); title('WLWHIV ages 0-79');
    xlim([1985 2030]); ylim([0 0.1]); grid on;
    legend('Model, ages 0-79: 25-sets mean' , 'Model, ages 0-79: 25-sets minimum' , ...
        'Model, ages 0-79: 25-sets maximum');
    
subplot(1,2,1);
plot(monthlyTimespan , mean(cinNegTime(: , :),1)' , 'k-' , ...
    monthlyTimespan , min(cinNegTime(: , :),[],1)' , 'k--' , ...
    monthlyTimespan , max(cinNegTime(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
xlabel('Year'); ylabel('CIN 2/3 Prevalence'); title('HIV-negative women ages 0-79');
xlim([1985 2030]); ylim([0 0.1]); grid on;
legend('Model, ages 0-79: 25-sets mean' , 'Model, ages 0-79: 25-sets minimum' , ...
    'Model, ages 0-79: 25-sets maximum');

%% CIN2/3 prevalence for All HR HPV types combined over time
figure;
plot(monthlyTimespan , mean(cinGenTime(: , :),1)' , 'k-' , ...
    monthlyTimespan , min(cinGenTime(: , :),[],1)' , 'k--' , ...
    monthlyTimespan , max(cinGenTime(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
    xlabel('Year'); ylabel('CIN 2/3 Prevalence'); title('All females aged 0-79');
    xlim([1985 2030]); ylim([0 0.1]); grid on;
    legend('Model, ages 0-79: 25-sets mean' , 'Model, ages 0-79: 25-sets minimum' , ...
        'Model, ages 0-79: 25-sets maximum');
    

%% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

%% Cervical cancer incidence in 2005 by age
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

figure;    
% General
plot(1 : age , mean(squeeze(ccIncAge(: , 1 , :)),1)' , 'k-' , ...
    1 : age , min(squeeze(ccIncAge(: , 1 , :)),[],1)' , 'k--' , ...
    1 : age , max(squeeze(ccIncAge(: , 1 , :)),[],1)' , 'k--' , 'LineWidth' , 1.5);
xlabel('Age Group'); ylabel('Cervical cancer incidence per 100K');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
ylim([0 50]); grid on;
title(['Cervical Cancer Incidence in 2005']);
legend('Model, general: 25-sets mean' , ...
    'Model: 25-sets minimum' , 'Model: 25-sets maximum');

%% Cervical cancer incidence in 2012 by age vs. Globocan 2012 data
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

% Calibration error bars
meanObs = ccInc2012_dObs(: , 2);
sdevObs = (ccInc2012_dObs(: , 3).^(1/2)).*2;

figure;    
% Plot observed data
errorbar(4 : age , meanObs , sdevObs , 'rs' , 'LineWidth' , 1.5);
hold all;
% General
plot(1 : age , mean(squeeze(ccIncAge(: , 2 , :)),1) , 'k-' , ...
    1 : age , min(squeeze(ccIncAge(: , 2 , :)),[],1) , 'k--' , ...
    1 : age , max(squeeze(ccIncAge(: , 2 , :)),[],1) , 'k--' , 'LineWidth' , 1.5);
xlabel('Age Group'); ylabel('Cervical cancer incidence per 100K');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
ylim([0 150]); grid on;
title(['Cervical Cancer Incidence in 2012']);
legend('(Globocan, 2012) Observed SA: mean, 2SD' , 'Model, general: 25-sets mean' , ...
    'Model: 25-sets minimum' , 'Model: 25-sets maximum');
   
%% Cervical cancer incidence in 2018 by age vs. Globocan 2018 data (calibration)
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

% Calibration error bars
meanObs = ccInc2018_dObs(: , 2);
sdevObs = (ccInc2018_dObs(: , 3).^(1/2)).*2;

figure;    
% Plot observed data
errorbar(4 : age-1 , meanObs , sdevObs , 'rs' , 'LineWidth' , 1.5);
hold all;
% General
plot(1 : age , mean(squeeze(ccIncAge(: , 3 , :)),1) , 'k-' , ...
    1 : age , min(squeeze(ccIncAge(: , 3 , :)),[],1) , 'k--' , ...
    1 : age , max(squeeze(ccIncAge(: , 3 , :)),[],1) , 'k--' , 'LineWidth' , 1.5);
xlabel('Age Group'); ylabel('Cervical cancer incidence per 100K');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
ylim([0 150]); grid on;
title(['Cervical Cancer Incidence in 2018']);
legend('(Globocan, 2018) Observed SA: mean, 2SD' , 'Model, general: 25-sets mean' , ...
    'Model: 25-sets minimum' , 'Model: 25-sets maximum');

%% Cervical cancer incidence over time
figure;   
% Plot observed data
% Find crude Globocan values
% General
plot(annualTimespan , mean(ccIncTime,1) , 'k-' , ...
    annualTimespan , min(ccIncTime,[],1) , 'k--' , ...
    annualTimespan , max(ccIncTime,[],1) , 'k--' , 'LineWidth' , 1.5);
hold all;
% HIV-negative
plot(annualTimespan , mean(ccIncTimeNeg,1) , 'g-' , ...
    annualTimespan , min(ccIncTimeNeg,[],1) , 'g--' , ...
    annualTimespan , max(ccIncTimeNeg,[],1) , 'g--' , 'LineWidth' , 1.5);
hold all;
% HIV-positive untreated
plot(annualTimespan , mean(ccIncTimePos,1) , 'r-' , ...
    annualTimespan , min(ccIncTimePos,[],1) , 'r--' , ...
    annualTimespan , max(ccIncTimePos,[],1) , 'r--' , 'LineWidth' , 1.5);
hold all;
% HIV-positive on ART
plot(annualTimespan , mean(ccIncTimeArt,1) , 'b-' , ...
    annualTimespan , min(ccIncTimeArt,[],1) , 'b--' , ...
    annualTimespan , max(ccIncTimeArt,[],1) , 'b--' , 'LineWidth' , 1.5);
xlabel('Time'); ylabel('Cervical cancer incidence per 100K');
xlim([1990 2020]); ylim([0 60]); grid on;
title(['Cervical Cancer Incidence']);
legend(... %'(Globocan, 2005) Observed SA crude' , '(Globocan, 2012) Observed SA crude' , ...
    ... %'(Globocan, 2018) Observed SA crude' , ...
    'Model, general: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, HIV-negative: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, WLWHIV untreated: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, WLWHIV on ART: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Location' , 'northwest');


%% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************

%% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)

% HPV infected
% Calibration error bars
meanObs = hpv_dist_dObs(: , 2);
sdevObs = (hpv_dist_dObs(: , 3).^(1/2)).*2;
figure;
subplot(2,3,1)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 0.4682 , 'k*');
hold on;
plot(2012 , 0.5318 , 'b*');
hold on;
plot(monthlyTimespan , mean(squeeze(hpv_vax(:,1,:)),1) , 'k-' , ...
    monthlyTimespan , min(squeeze(hpv_vax(:,1,:)),[],1) , 'k--' , ...
    monthlyTimespan , max(squeeze(hpv_vax(:,1,:)),[],1) , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(monthlyTimespan , mean(squeeze(hpv_nonVax(:,1,:)),1) , 'b-' , ...
    monthlyTimespan , min(squeeze(hpv_nonVax(:,1,:)),[],1) , 'b--' , ...
    monthlyTimespan , max(squeeze(hpv_nonVax(:,1,:)),[],1) , 'b--' , 'LineWidth' , 1.5);
xlabel('Year'); ylabel('Prevalence Proportion by Type');
title('HPV');
ylim([0 1]);
xlim([2010 2015]);
grid on;

% CIN1
% Calibration error bars
meanObs = cin1_dist_dObs(: , 2);
sdevObs = (cin1_dist_dObs(: , 3).^(1/2)).*2;
subplot(2,3,2)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 0.5192 , 'k*');
hold on;
plot(2012 , 0.4808 , 'b*');
hold on;
plot(monthlyTimespan , mean(squeeze(cin1_vax(:,1,:)),1) , 'k-' , ...
    monthlyTimespan , min(squeeze(cin1_vax(:,1,:)),[],1) , 'k--' , ...
    monthlyTimespan , max(squeeze(cin1_vax(:,1,:)),[],1) , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(monthlyTimespan , mean(squeeze(cin1_nonVax(:,1,:)),1) , 'b-' , ...
    monthlyTimespan , min(squeeze(cin1_nonVax(:,1,:)),[],1) , 'b--' , ...
    monthlyTimespan , max(squeeze(cin1_nonVax(:,1,:)),[],1) , 'b--' , 'LineWidth' , 1.5);
ylim([0 1]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type');
title('CIN1');
grid on;

% CIN2
subplot(2,3,3);
plot(2012 , 0.6281 , 'k*');
hold on;
plot(2012 , 0.3719 , 'b*');
hold on;
plot(monthlyTimespan , mean(squeeze(cin2_vax(:,1,:)),1) , 'k-' , ...
    monthlyTimespan , min(squeeze(cin2_vax(:,1,:)),[],1) , 'k--' , ...
    monthlyTimespan , max(squeeze(cin2_vax(:,1,:)),[],1) , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(monthlyTimespan , mean(squeeze(cin2_nonVax(:,1,:)),1) , 'b-' , ...
    monthlyTimespan , min(squeeze(cin2_nonVax(:,1,:)),[],1) , 'b--' , ...
    monthlyTimespan , max(squeeze(cin2_nonVax(:,1,:)),[],1) , 'b--' , 'LineWidth' , 1.5);
ylim([0 1]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type');
title('CIN2');
grid on;

% CIN3
% Calibration error bars
meanObs = cin3_dist_dObs(: , 2);
sdevObs = (cin3_dist_dObs(: , 3).^(1/2)).*2;
subplot(2,3,4)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 0.7371 , 'k*');
hold on;
plot(2012 , 0.2629 , 'b*');
hold on;
plot(monthlyTimespan , mean(squeeze(cin3_vax(:,1,:)),1) , 'k-' , ...
    monthlyTimespan , min(squeeze(cin3_vax(:,1,:)),[],1) , 'k--' , ...
    monthlyTimespan , max(squeeze(cin3_vax(:,1,:)),[],1) , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(monthlyTimespan , mean(squeeze(cin3_nonVax(:,1,:)),1) , 'b-' , ...
    monthlyTimespan , min(squeeze(cin3_nonVax(:,1,:)),[],1) , 'b--' , ...
    monthlyTimespan , max(squeeze(cin3_nonVax(:,1,:)),[],1) , 'b--' , 'LineWidth' , 1.5);
ylim([0 1]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type');
title('CIN3');
grid on;

% CC
% Calibration error bars
meanObs = cc_dist_dObs(: , 2);
sdevObs = (cc_dist_dObs(: , 3).^(1/2)).*2;
subplot(2,3,5)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
plot(2012 , 0.8578 , 'k*');
hold on;
plot(2012 , 0.1422 , 'b*');
hold on;
plot(monthlyTimespan , mean(squeeze(cc_vax(:,1,:)),1) , 'k-' , ...
    monthlyTimespan , min(squeeze(cc_vax(:,1,:)),[],1) , 'k--' , ...
    monthlyTimespan , max(squeeze(cc_vax(:,1,:)),[],1) , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(monthlyTimespan , mean(squeeze(cc_nonVax(:,1,:)),1) , 'b-' , ...
    monthlyTimespan , min(squeeze(cc_nonVax(:,1,:)),[],1) , 'b--' , ...
    monthlyTimespan , max(squeeze(cc_nonVax(:,1,:)),[],1) , 'b--' , 'LineWidth' , 1.5);
ylim([0 1]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type')
title('Cervical Cancer')
legend('Observed 2012: mean, 2SD' , 'Observed 2012- 9v' , 'Observed 2012- non-9v' , ...
    'Model- 9v: 25-sets mean' , 'Model- 9v: 25-sets minimum' , 'Model- 9v: 25-sets maximum' , ...
    'Model- non-9v: 25-sets mean' , 'Model- non-9v: 25-sets minimum' , 'Model- non-9v: 25-sets maximum');
grid on;

%% HPV type distribution by state in 2000 and 2018 (coinfections grouped as 9v-type HPV)
% HPV infected
% Calibration error bars
meanObs = [hpv_dist_dObs(: , 2) , [((0.6281+0.7371)/2); ((0.3719+0.2629)/2)] , cc_dist_dObs(: , 2)];
sdevObs = [((hpv_dist_dObs(: , 3).^(1/2)).*2) , [0;0] , ...
    ((cc_dist_dObs(: , 3).^(1/2)).*2)];

figure;
errorbar([1 : 3] , meanObs(1 , :) , sdevObs(1 , :) , 'ks' , 'LineWidth' , 1.5);
hold all;
errorbar([1 : 3] , meanObs(2 , :) , sdevObs(2 , :) , 'bs' , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [mean(hpv_vax(:,1,((2000 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cin23_vax(:,1,((2000 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_vax(:,1,((2000 - startYear) * stepsPerYear +1)),1)'] , 'k-' , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [min(hpv_vax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cin23_vax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_vax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)'] , 'k--' , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [max(hpv_vax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cin23_vax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_vax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)'] , 'k--' , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [mean(hpv_vax(:,1,((2018 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cin23_vax(:,1,((2018 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_vax(:,1,((2018 - startYear) * stepsPerYear +1)),1)'] , '-' , 'Color' , [0.7 0.7 0.7] , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [min(hpv_vax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cin23_vax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_vax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)'] , '--' , 'Color' , [0.7 0.7 0.7] , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [max(hpv_vax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cin23_vax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_vax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)'] , '--' , 'Color' , [0.7 0.7 0.7] , 'LineWidth' , 1.5);
hold all;

plot([1 : 3] , [mean(hpv_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cin23_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),1)'] , 'b-' , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [min(hpv_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cin23_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)'] , 'b--' , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [max(hpv_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cin23_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)'] , 'b--' , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [mean(hpv_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cin23_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),1)'] , 'c-' , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [min(hpv_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cin23_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)'] , 'c--' , 'LineWidth' , 1.5);
hold all;
plot([1 : 3] , [max(hpv_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cin23_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)'] , 'c--' , 'LineWidth' , 1.5);
xlabel('Stage');
set(gca , 'xtick' , 1 : 3 , 'xtickLabel' , {'HPV' , 'CIN2/3' , 'CC'});
ylabel('Prevalence Proportion by Type');
title('Type Distribution by stage (coinfections grouped as 9v)');
ylim([0 1]); grid on;
legend('Observed- 9v: mean, 2SD' , 'Observed- non-9v: mean, 2SD' , ...
    'Model- 9v, 2000: 25-sets mean' , 'Model- 9v, 2000: 25-sets minimum' , 'Model- 9v, 2000: 25-sets maximum' , ...
    'Model- 9v, 2018: 25-sets mean' , 'Model- 9v, 2018: 25-sets minimum' , 'Model- 9v, 2018: 25-sets maximum' , ...
    'Model- non-9v, 2000: 25-sets mean' , 'Model- non-9v, 2000: 25-sets minimum' , 'Model- non-9v, 2000: 25-sets maximum' , ...
    'Model- non-9v, 2018: 25-sets mean' , 'Model- non-9v, 2018: 25-sets minimum' , 'Model- non-9v, 2018: 25-sets maximum');

%% HPV type distribution of cervical cancer in 2000 and 2018 by HIV status(coinfections grouped as 9v-type HPV)
figure;
plot([1 : 4] , [mean(cc_vax(:,1,((2000 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_vax(:,2,((2000 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_vax(:,3,((2000 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_vax(:,4,((2000 - startYear) * stepsPerYear +1)),1)'] , 'k-' , 'LineWidth' , 1.5);
hold all;
plot([1 : 4] , [min(cc_vax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_vax(:,2,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_vax(:,3,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_vax(:,4,((2000 - startYear) * stepsPerYear +1)),[],1)'] , 'k--' , 'LineWidth' , 1.5);
hold all;
plot([1 : 4] , [max(cc_vax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_vax(:,2,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_vax(:,3,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_vax(:,4,((2000 - startYear) * stepsPerYear +1)),[],1)'] , 'k--' , 'LineWidth' , 1.5);
hold all;
plot([1 : 4] , [mean(cc_vax(:,1,((2018 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_vax(:,2,((2018 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_vax(:,3,((2018 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_vax(:,4,((2018 - startYear) * stepsPerYear +1)),1)'] , '-' , 'Color' , [0.7 0.7 0.7] , 'LineWidth' , 1.5);
hold all;
plot([1 : 4] , [min(cc_vax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_vax(:,2,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_vax(:,3,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_vax(:,4,((2018 - startYear) * stepsPerYear +1)),[],1)'] , '--' , 'Color' , [0.7 0.7 0.7] , 'LineWidth' , 1.5);
hold all;
plot([1 : 4] , [max(cc_vax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_vax(:,2,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_vax(:,3,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_vax(:,4,((2018 - startYear) * stepsPerYear +1)),[],1)'] , '--' , 'Color' , [0.7 0.7 0.7] , 'LineWidth' , 1.5);
hold all;

plot([1 : 4] , [mean(cc_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_nonVax(:,2,((2000 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_nonVax(:,3,((2000 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_nonVax(:,4,((2000 - startYear) * stepsPerYear +1)),1)'] , 'b-' , 'LineWidth' , 1.5);
hold all;
plot([1 : 4] , [min(cc_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_nonVax(:,2,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_nonVax(:,3,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_nonVax(:,4,((2000 - startYear) * stepsPerYear +1)),[],1)'] , 'b--' , 'LineWidth' , 1.5);
hold all;
plot([1 : 4] , [max(cc_nonVax(:,1,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_nonVax(:,2,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_nonVax(:,3,((2000 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_nonVax(:,4,((2000 - startYear) * stepsPerYear +1)),[],1)'] , 'b--' , 'LineWidth' , 1.5);
hold all;
plot([1 : 4] , [mean(cc_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_nonVax(:,2,((2018 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_nonVax(:,3,((2018 - startYear) * stepsPerYear +1)),1)' , ...
    mean(cc_nonVax(:,4,((2018 - startYear) * stepsPerYear +1)),1)'] , 'c-' , 'LineWidth' , 1.5);
hold all;
plot([1 : 4] , [min(cc_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_nonVax(:,2,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_nonVax(:,3,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    min(cc_nonVax(:,4,((2018 - startYear) * stepsPerYear +1)),[],1)'] , 'c--' , 'LineWidth' , 1.5);
hold all;
plot([1 : 4] , [max(cc_nonVax(:,1,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_nonVax(:,2,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_nonVax(:,3,((2018 - startYear) * stepsPerYear +1)),[],1)' , ...
    max(cc_nonVax(:,4,((2018 - startYear) * stepsPerYear +1)),[],1)'] , 'c--' , 'LineWidth' , 1.5);
xlabel('Stage');
set(gca , 'xtick' , 1 : 4 , 'xtickLabel' , {'General' , 'HIV-negative' , 'HIV+, untreated' , 'HIV+ on ART'});
ylabel('Prevalence Proportion by Type');
title('Cervical cancer type distribution by HIV status (coinfections grouped as 9v)');
ylim([0 1]); grid on;
legend('Model- 9v, 2000: 25-sets mean' , 'Model- 9v, 2000: 25-sets minimum' , 'Model- 9v, 2000: 25-sets maximum' , ...
    'Model- 9v, 2018: 25-sets mean' , 'Model- 9v, 2018: 25-sets minimum' , 'Model- 9v, 2018: 25-sets maximum' , ...
    'Model- non-9v, 2000: 25-sets mean' , 'Model- non-9v, 2000: 25-sets minimum' , 'Model- non-9v, 2000: 25-sets maximum' , ...
    'Model- non-9v, 2018: 25-sets mean' , 'Model- non-9v, 2018: 25-sets minimum' , 'Model- non-9v, 2018: 25-sets maximum');

