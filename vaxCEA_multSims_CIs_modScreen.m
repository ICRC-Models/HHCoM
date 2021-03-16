function [] = vaxCEA_multSims_CIs_modScreen(vaxResultInd , sceNum , fileNameNums)
% example: vaxCEA_multSims_CIs_modScreen(1 , '0' , {'0' , '__'})

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
    baseline , who , spCyto , spHpvDna , spGentyp , spAve , spHpvAve , ...
    circProtect , condProtect , MTCTRate , hyst , ...
    OMEGA , ...
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
    toNonVaxScrnInds , ageInd , riskInd , ...
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

lastYear = 2121;

% Indices of calib runs to plot
fileInds = {'6_1' , '6_2' , '6_3' , '6_6' , '6_8' , '6_9' , '6_11' , ...
     '6_12' , '6_13' , '6_15' , '6_20' , '6_21' , '6_22' , '6_26' , ...
    '6_27' , '6_32' , '6_34' , '6_35' , '6_38' , '6_39' , '6_40' , ...
    '6_41' , '6_42' , '6_45' , '6_47'};    % 22Apr20Ph2V11
nRuns = length(fileInds);

% Initialize model output plots
% Timespans
monthlyTimespan = [startYear : (1/6) : lastYear];
monthlyTimespan = monthlyTimespan(1 : end-1);
annualTimespan = [startYear : lastYear-1];
futAnnualTimespan = [2020 : lastYear-1];
midAnnualTimespan = [(startYear+(3/stepsPerYear)) : ((lastYear-1)+(3/stepsPerYear))];
% Total population size
popSize = zeros(nRuns , length(monthlyTimespan));
% Population age distribution
popYearVec = [2018 2100];
popYearVecLength = length(popYearVec);
popPropF = zeros(nRuns , length(popYearVec) , age);
% HIV prevalence
hivAgeM = zeros(nRuns , age , length(monthlyTimespan));
hivAgeF = hivAgeM;
hivPrev = zeros(nRuns , gender , length(monthlyTimespan));
hivPrevW = zeros(nRuns , length(monthlyTimespan));
% ART coverage
artCovM = zeros(nRuns , length(monthlyTimespan));
artCovF = artCovM;
artCovW = artCovM;
artCovAge = zeros(nRuns , age , length(monthlyTimespan));
% VMMC coverage
ageVecCirc = {4 , 5 , [6:10] , [11:age]}; % Ages: (15-19), (20-24), (25-49), (50+)
ageVecCirc_length = length(ageVecCirc);
circProp = zeros(nRuns , ageVecCirc_length , length(monthlyTimespan));
% Female HPV prevalence
hpvYearVec_orig = [2002 2018];
hpvYearVec2018_orig = [2018];
hpvYearVecLength = length(hpvYearVec_orig);
hpvYearVec2018Length = length(hpvYearVec2018_orig);
hpv_hiv = zeros(nRuns , age , 2);
hpv_hivNeg = hpv_hiv;
hpv_hivTot = zeros(nRuns , age , 1);
% Male HPV prevalence
hpvYearVecMale = [2008 2018];
hpvYearVecMaleLength = length(hpvYearVecMale);
hpv_hivM = zeros(nRuns , 4 , 2);
hpv_hivMNeg = hpv_hivM;
hpv_hivMtot = hpv_hivM;
ageVec_mHPV = {[4:5],[6:7],[8:9],[10:13]};
ageVecLength_mHPV = length(ageVec_mHPV);
% Combined HPV prevalence
hpv_hivW = zeros(nRuns , 5 , length(monthlyTimespan));
hpv_hivAgeW = zeros(nRuns , 5 , age , length(monthlyTimespan));
% CIN2/3 prevalence
cinPosAge = zeros(nRuns , 2 , age);
cinNegAge = cinPosAge;
cinGenAge = cinPosAge;
% CC incidence
ccYearVec = [2005 2012 2018];
ccYearVecLength = length(ccYearVec);
ccIncAge = zeros(nRuns , 3 , age);
ccIncTime = zeros(nRuns , length(annualTimespan));
ccIncTimeNeg = ccIncTime;
ccIncTimePos = ccIncTime;
ccIncTimeArt = ccIncTime;
ccIncTimePosAll = ccIncTime;
ccIncHivAgeTime = zeros(nRuns , 5 , age , length(annualTimespan));
ccIncHivTime = zeros(nRuns , 5 , length(annualTimespan));
ccCumHivTime = zeros(nRuns , 5 , length(futAnnualTimespan));
ccAnlHivTime = zeros(nRuns , 5 , length(annualTimespan));
ccAnlHivAgeTime = zeros(nRuns , 5 , age , length(annualTimespan));
diseaseVec_ccInc = {[1 : disease] , [1 : 2] , [3 : 8] , [3 : 7] , 8};
diseaseVecLength_ccInc = length(diseaseVec_ccInc);
% HPV/CIN/CC type distribution
diseaseInds_typeDist = {[1 : disease] , [1 : 2] , [3 : 7] , 8 , [3 : 8]};
diseaseIndsLength_typeDist = length(diseaseInds_typeDist);
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
% HPV vaccination and screening
overScreenTotAnnualCum = zeros(nRuns , 5 , 5 , length(futAnnualTimespan));
vaxCoverage = zeros(nRuns , length(monthlyTimespan));
vaxCoverageAge = zeros(nRuns , age , length(monthlyTimespan));

resultsDir = [pwd , '\HHCoM_Results\'];
fileKey = {'sim1' , 'sim0'};
fileKeyNums = fileNameNums;
n = vaxResultInd;
baseFileName = ['22Apr20Ph2V11_2v57BaseVax_spCytoScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_CISNET-S' , sceNum , '_'];
loopSegments = {0 , round(nRuns/2) , nRuns};
loopSegmentsLength = length(loopSegments);
for k = 1 : loopSegmentsLength-1
    parfor j = loopSegments{k}+1 : loopSegments{k+1}
        % Load results
        pathModifier = [baseFileName , fileInds{j}]; % ***SET ME***: name for simulation output file
        nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
        curr = load([pwd , '/HHCoM_Results/toNow_22Apr20Ph2V11_2v57BaseVax_spCytoScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_obsHist_' , fileInds{j}]); % ***SET ME***: name for historical run output file 

        vaxResult = cell(nSims , 1);
        resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxSimResult'];
        % load results from vaccine run into cell array
        vaxResult{n} = load([resultFileName , num2str(n), '.mat']);
        % concatenate vectors/matrices of population up to current year to population
        % matrices for years past current year
        vaxResult{n}.popVec = [curr.popVec(1 : end  , :); vaxResult{n}.popVec(2 : end , :)];
        vaxResult{n}.ccDeath = [curr.ccDeath(1 : end , : , : , :) ; vaxResult{n}.ccDeath(2 : end , : , : , :)];
        vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :); vaxResult{n}.newCC(2 : end , : , : , :)];
        vaxResult{n}.newHpvVax = [curr.newHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvVax(2 : end , : , : , : , : , :)];
        vaxResult{n}.newImmHpvVax = [curr.newImmHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvVax(2 : end , : , : , : , : , :)];
        vaxResult{n}.newHpvNonVax = [curr.newHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvNonVax(2 : end , : , : , : , : , :)];
        vaxResult{n}.newImmHpvNonVax = [curr.newImmHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvNonVax(2 : end , : , : , : , : , :)];
        vaxResult{n}.newScreen = [vaxResult{n}.newScreen(1 : end , : , : , : , : , : , :)]; %[curr.newScreen(1 : end , : , : , : , : , : , : ); vaxResult{n}.newScreen(2 : end , : , : , : , : , : , :)];
        vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , : , : , : , :); vaxResult{n}.newHiv(2 : end , : , : , : , : , : , :)];
        vaxResult{n}.hivDeaths = [curr.hivDeaths(1 : end , : , : , :); vaxResult{n}.hivDeaths(2 : end , : , : , :)];
        vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :); vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , :)];
        vaxResult{n}.tVec = [curr.tVec(1 : end), vaxResult{n}.tVec(2 : end)];

    %     noVaxInd = nSims;
    %     noV = vaxResult{noVaxInd};
        tVec = vaxResult{n}.tVec;
        tVecYr = tVec(1 : stepsPerYear : end);
        
        % Initiatlize variables
        hpvYearVec = hpvYearVec_orig;
        hpvYearVec2018 = hpvYearVec2018_orig;

        %% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

        %% Population size over time vs. Statistics South Africa data (calibration)
        popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
        popSize(j , :) = sum(vaxResult{n}.popVec(: , popTot),2);
        
        %% Female population proportion by 5-year age groups over time vs. Statistics South Africa data (internal validation)
        for t = 1 : popYearVecLength
            for a = 1 : age
                popAgeF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                popTotF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
                popPropF(j , t , a) = sum(vaxResult{n}.popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popAgeF),2) ./ ...
                    sum(vaxResult{n}.popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTotF),2);
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
            hivAgeM(j , a , :) =  (sum(vaxResult{n}.popVec(: , hivMInds) , 2) + sum(vaxResult{n}.popVec(: , artMInds) , 2)) ...
                ./ sum(vaxResult{n}.popVec(: , totMInds) , 2);
            hivAgeF(j , a , :) =  (sum(vaxResult{n}.popVec(: , hivFInds) , 2) + sum(vaxResult{n}.popVec(: , artFInds) , 2)) ...
                ./ sum(vaxResult{n}.popVec(: , totFInds) , 2);
        end
    
        %% HIV prevalence by gender over time vs. AHRI (validation) and (Vandormael, 2019) AHRI data (validation)
        for g = 1 : gender
                hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                    1 : intervens , g , 4 : 10 , 1 : risk));
                totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                    1 : intervens , g , 4 : 10 , 1 : risk));
                hivPrev(j , g , :) = (sum(vaxResult{n}.popVec(: , hivInds) , 2) ./ sum(vaxResult{n}.popVec(: , totInds) , 2)) .* 100;
        end
        
        %% HIV prevalence in women over time for all ages
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 1 : age , 1 : risk));
        totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 1 : age , 1 : risk));
        hivPrevW(j , :) = (sum(vaxResult{n}.popVec(: , hivInds) , 2) ./ sum(vaxResult{n}.popVec(: , totInds) , 2)) .* 100;
               
        %% Proportion of total age-eligible HIV+ population on ART and VS (denominator: CD4-eligible and ineligible)
        artIndsF = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 3 : age , 1 : risk));
        hivAllIndsF = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
            1 : endpoints , 1 : intervens , 2 , 3 : age , 1 : risk));
        artIndsM = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , 3 : age , 1 : risk));
        hivAllIndsM = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
            1 : endpoints , 1 : intervens , 1 , 3 : age , 1 : risk));
        
        artCovF(j , :) = (sum(vaxResult{n}.popVec(: , artIndsF) , 2) ./ (sum(vaxResult{n}.popVec(: , hivAllIndsF) , 2) + sum(vaxResult{n}.popVec(: , artIndsF) , 2))) .* 100;
        artCovM(j , :) = (sum(vaxResult{n}.popVec(: , artIndsM) , 2) ./ (sum(vaxResult{n}.popVec(: , hivAllIndsM) , 2) + sum(vaxResult{n}.popVec(: , artIndsM) , 2))) .* 100;
        
        %% Proportion of total HIV+ population on ART and VS (denominator: CD4-eligible and ineligible); all ages
        artIndsW = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        hivAllIndsW = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        artCovW(j , :) = (sum(vaxResult{n}.popVec(: , artIndsW) , 2) ./ (sum(vaxResult{n}.popVec(: , hivAllIndsW) , 2) + sum(vaxResult{n}.popVec(: , artIndsW) , 2))) .* 100;
        
        %% Proportion of total HIV+ population on ART and VS by age
        for a = 1 : age
            artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            artPop = sum(vaxResult{n}.popVec(: , artInds) , 2);
            hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            hivPop = sum(vaxResult{n}.popVec(: , hivInds) , 2);
            artCovAge(j , a , :) = artPop ./ hivPop;
        end
        
        %% Proportion HIV-negative males circumcised by broad age groups over time
        for aInd = 1 : ageVecCirc_length
            aGroup = ageVecCirc{aInd};
            circInds = toInd(allcomb(2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , aGroup , 1 : risk));
            circPop = sum(vaxResult{n}.popVec(: , circInds) , 2);
            hivNegInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
                1 : endpoints , 1 : intervens , 1 , aGroup , 1 : risk));
            hivNegPop = sum(vaxResult{n}.popVec(: , hivNegInds) , 2);
            circProp(j , aInd , :) = 100 .* (circPop ./ hivNegPop);
        end
    
        
        %% ********************************** HPV FIGURES **********************************************************************************************
    
        %% Female HPV Prevalence by age and HIV status in 2002 and 2018 vs. McDonald 2014 data (calibration)
        for i = 1 : hpvYearVecLength
            yr = hpvYearVec(i);
            for a = 1 : age % 15-19 -> 55-65
                hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
                    [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
                ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                hpv_hiv(j , a , i) = sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
                    ./ sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    
                hpvInds_hivNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
                    [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
                ageInds_hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                hpv_hivNeg(j , a , i) = sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivNeg))...
                    ./ sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivNeg));
            end
        end
        
        %% Female HPV Prevalence by age 2018
        for i = 1 : hpvYearVec2018Length
            yr = hpvYearVec2018(i);
            for a = 1 : age
                hpvInds_hivTot = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
                    [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
                ageInds_hivTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                hpv_hivTot(j , a , i) = sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivTot))...
                    ./ sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivTot));
            end
        end
       
        %% Male HPV prevalence by age and HIV status in 2008 vs. Mbulawa data (calibration)
        for i = 1 : hpvYearVecMaleLength
            yr = hpvYearVecMale(i);
            for aV = 1 : ageVecLength_mHPV
                aGroup = ageVec_mHPV{aV};
                hpvInds_hivM = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 , [1 : 2 , 7] , ...
                    1 , 1 : intervens , 1 , aGroup , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
                    [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , aGroup , 1 : risk))]);
                ageInds_hivM = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 1 , aGroup , 1 : risk));
                hpv_hivM(j , aV , i) = sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivM))...
                    / sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear+1 , ageInds_hivM));
    
                hpvInds_hivMNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 , [1 : 2 , 7] , ...
                    1 , 1 : intervens , 1 , aGroup , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
                    [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , aGroup , 1 : risk))]);
                ageInds_hivMNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 1 , aGroup , 1 : risk));
                hpv_hivMNeg(j , aV , i) = sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivMNeg))...
                    / sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivMNeg));
                
                hpvInds_hivMtot = unique([toInd(allcomb(1 : disease , 1 : viral , 2 , [1 : 2 , 7] , ...
                    1 , 1 : intervens , 1 , aGroup , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
                    [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , aGroup , 1 : risk))]);
                ageInds_hivMtot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 1 , aGroup , 1 : risk));
                hpv_hivMtot(j , aV , i) = sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivMNeg))...
                    / sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivMNeg));
            end
        end
        
        %% HPV prevalence by HIV status over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
                hpvInds = unique([toInd(allcomb(d , 1 : viral , 2 : 6 , 1 : hpvNonVaxStates , ...
                    1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                    1 : hpvVaxStates , 2 : 5 , 1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
                ageInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
                hpv_hivW(j , dInd , :) = sum(vaxResult{n}.popVec(: , hpvInds) , 2)...
                    ./ sum(vaxResult{n}.popVec(: , ageInds) , 2);
        end
        
        %% HPV prevalence by HIV status and age over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                hpvInds = unique([toInd(allcomb(d , 1 : viral , 2 : 6 , 1 : hpvNonVaxStates , ...
                    1 : 3 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                    1 : hpvVaxStates , 2 : 5 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk))]);
                ageInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                hpv_hivAgeW(j , dInd , a , :) = sum(vaxResult{n}.popVec(: , hpvInds) , 2)...
                    ./ sum(vaxResult{n}.popVec(: , ageInds) , 2);
            end
        end
        
        
        %% ********************************** CIN FIGURES *********************************************************************************************
    
        %% CIN2/3 prevalence for All HR HPV types combined by HIV status and age in 2002 vs. McDonald 2014 data (calibration)
        for i = 1 : hpvYearVecLength
            yr = hpvYearVec(i);
            for a = 1 : age % date is 15-19 -> 60-64
                % HIV-positive (on and not on ART)
                cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
                    [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
                ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                cinPosAge(j , i , a) = (sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
                    ./ sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
                % HIV-negative
                cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
                    [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
                ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                cinNegAge(j , i , a) = (sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , cinNegInds)))...
                    ./ (sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , ageNegInds)));
                % General
                cinGenInds = unique([toInd(allcomb(1 : disease , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
                    [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
                ageGenInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                cinGenAge(j , i , a) = (sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , cinGenInds)))...
                    ./ (sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , ageGenInds)));
            end
        end
    
        
        %% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************
    
        %% Cervical cancer incidence in 2005, 2012, 2018 by age vs. Globocan data and other sources (calibration)
        fac = 10 ^ 5;
        for i = 1 : ccYearVecLength
            yr = ccYearVec(i);
            incTimeSpan = [((yr - startYear) * stepsPerYear +1) : ((yr - startYear) * stepsPerYear +6)];
            for a = 1 : age
                allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                % Calculate incidence
                ccIncAge(j , i , a) = ...
                    (annlz(sum(sum(sum(vaxResult{n}.newCC(incTimeSpan , : , a , :),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{n}.popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);
            end
        end
        
        %% Cervical cancer incidence over time
        fac = 10 ^ 5;
        % General population
        allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 15 , 1 : risk));
        % Calculate incidence
        ccIncTime(j , :) = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , : , 4 : 15 , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-negative
        allFneg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 15 , 1 : risk));
        % Calculate incidence
        ccIncTimeNeg(j , :) = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , 1 : 2 , 4 : 15 , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , allFneg) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-positive untreated
        allFpos = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 15 , 1 : risk));
        % Calculate incidence
        ccIncTimePos(j , :) = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , 3 : 7 , 4 : 15 , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , allFpos) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-positive on ART
        allFart = toInd(allcomb(8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 15 , 1 : risk));
        % Calculate incidence
        ccIncTimeArt(j , :) = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , 8 , 4 : 15 , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , allFart) , 2) ./ stepsPerYear)) * fac);
        
        % HIV-positive all
        allFposAll = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 15 , 1 : risk));
        % Calculate incidence
        ccIncTimePosAll(j , :) = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , 3 : 8 , 4 : 15 , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , allFposAll) , 2) ./ stepsPerYear)) * fac);
        
        %% Cervical cancer incidence by HIV status and age over time
        fac = 10 ^ 5;
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                allF = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                % Calculate incidence
                ccIncHivAgeTime(j , dInd , a , :) = ...
                    (annlz(sum(sum(sum(vaxResult{n}.newCC(: , d , a , :),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{n}.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);
            end
        end
        
        %% Cervical cancer incidence by HIV status over time
        fac = 10 ^ 5;
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            allF = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
            % Calculate incidence
            ccIncHivTime(j , dInd , :) = ...
                (annlz(sum(sum(sum(vaxResult{n}.newCC(: , d , 3 : age , :),2),3),4)) ./ ...
                (annlz(sum(vaxResult{n}.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);
        end
        
        %% Cumulative cervical cancer cases by HIV status over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
                ccCumHivTime(j , dInd , :) = ...
                    cumsum(squeeze(annlz(sum(sum(sum(vaxResult{n}.newCC(((2020 - startYear) * stepsPerYear +1):end , d , : , :),2),3),4))),2);
        end
        
        %% Annual cervical cancer cases by HIV status over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            ccAnlHivTime(j , dInd , :) = ...
                squeeze(annlz(sum(sum(sum(vaxResult{n}.newCC(: , d , : , :),2),3),4)));
        end
        
        %% Annual cervical cancer cases by HIV status and age over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                ccAnlHivAgeTime(j , dInd , a , :) = ...
                    squeeze(annlz(sum(sum(sum(vaxResult{n}.newCC(: , d , a , :),2),3),4)));
            end
        end
        
       
        %% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************
        
        %% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)
        for dInd = 1 : diseaseIndsLength_typeDist
            d = diseaseInds_typeDist{dInd};
            ccInds_vax = toInd(allcomb(d , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
                1 : 3 , 1 : intervens , 2 , 4 : 15 , 1 : risk));
            ccInds_nonVax = toInd(allcomb(d , 1 : viral , [1 : 5 , 7] , 6 , ...
                1 : 3 , 1 : intervens , 2 , 4 : 15 , 1 : risk));
            ccInds_tot = unique([toInd(allcomb(d , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
                    1 : 3 , 1 : intervens , 2 , 4 : 15 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                    [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 4 : 15 , 1 : risk))]);
    
           cin23Inds_vax = [toInd(allcomb(d , 1 : viral , 4 , [1 : 4 , 7] , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk)); ...
                toInd(allcomb(d , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk))];
            cin23Inds_nonVax = [toInd(allcomb(d , 1 : viral , [1 : 3 , 7] , 4 , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk)); ...
                toInd(allcomb(d , 1 : viral , [1 : 4 , 7] , 5 , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk))];
            cin23Inds_tot = unique([toInd(allcomb(d , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk)); ...
                toInd(allcomb(d , 1 : viral , 4 , [1 : 4 , 7] , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk)); 
                toInd(allcomb(d , 1 : viral , ...
                [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 4 : 15 , 1 : risk)); ...
                toInd(allcomb(d , 1 : viral , ...
                [1 3 , 7] , 4 , 1 , 1 : intervens , 2 , 4 : 15 , 1 : risk))]);     
    
            cin3Inds_vax = toInd(allcomb(d , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk));
            cin3Inds_nonVax = toInd(allcomb(d , 1 : viral , [1 : 4 , 7] , 5 , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk));
            cin3Inds_tot = unique([toInd(allcomb(d , 1 : viral , 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , 4 : 15 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                    [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 4 : 15 , 1 : risk))]);
    
            cin2Inds_vax = toInd(allcomb(d , 1 : viral , 4 , [1 : 4 , 7] , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk));
            cin2Inds_nonVax = toInd(allcomb(d , 1 : viral , [1 : 3 , 7] , 4 , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk));
            cin2Inds_tot = unique([toInd(allcomb(d , 1 : viral , 4 , [1 : 4 , 7] , ...
                    1 , 1 : intervens , 2 , 4 : 15 , 1 : risk)); toInd(allcomb(d, 1 : viral , ...
                    [1 : 3 , 7] , 4 , 1 , 1 : intervens , 2 , 4 : 15 , 1 : risk))]);
    
            cin1Inds_vax = toInd(allcomb(d , 1 : viral , 3 , [1 : 3 , 7] , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk));
            cin1Inds_nonVax = toInd(allcomb(d , 1 : viral , [1 : 2 , 7] , 3 , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk));
            cin1Inds_tot = unique([toInd(allcomb(d , 1 : viral , 3 , [1 : 3 , 7] , ...
                    1 , 1 : intervens , 2 , 4 : 15 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                    [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 4 : 15 , 1 : risk))]);
    
            hpvInds_vax = toInd(allcomb(d , 1 : viral , 2 , [1 : 2 , 7] , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk));
            hpvInds_nonVax = toInd(allcomb(d , 1 : viral , [1 , 7] , 2 , ...
                1 , 1 : intervens , 2 , 4 : 15 , 1 : risk));
            hpvInds_tot = unique([toInd(allcomb(d , 1 : viral , 2 , [1 : 2 , 7] , ...
                    1 , 1 : intervens , 2 , 4 : 15 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                    [1 , 7] , 2 , 1 , 1 : intervens , 2 , 4 : 15 , 1 : risk))]);
    
            cc_vax(j , dInd , :) = sum(vaxResult{n}.popVec(: , ccInds_vax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , ccInds_tot) , 2);
            cc_nonVax(j , dInd , :) = sum(vaxResult{n}.popVec(: , ccInds_nonVax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , ccInds_tot) , 2);
   
            cin23_vax(j , dInd , :) = sum(vaxResult{n}.popVec(: , cin23Inds_vax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , cin23Inds_tot) , 2);
            cin23_nonVax(j , dInd , :) = sum(vaxResult{n}.popVec(: , cin23Inds_nonVax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , cin23Inds_tot) , 2);
    
            cin3_vax(j , dInd , :) = sum(vaxResult{n}.popVec(: , cin3Inds_vax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , cin3Inds_tot) , 2);
            cin3_nonVax(j , dInd , :) = sum(vaxResult{n}.popVec(: , cin3Inds_nonVax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , cin3Inds_tot) , 2);
    
            cin2_vax(j , dInd , :) = sum(vaxResult{n}.popVec(: , cin2Inds_vax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , cin2Inds_tot) , 2);
            cin2_nonVax(j , dInd , :) = sum(vaxResult{n}.popVec(: , cin2Inds_nonVax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , cin2Inds_tot) , 2);
    
            cin1_vax(j , dInd , :) = sum(vaxResult{n}.popVec(: , cin1Inds_vax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , cin1Inds_tot) , 2);
            cin1_nonVax(j , dInd , :) = sum(vaxResult{n}.popVec(: , cin1Inds_nonVax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , cin1Inds_tot) , 2);
    
            hpv_vax(j , dInd , :) = sum(vaxResult{n}.popVec(: , hpvInds_vax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , hpvInds_tot) , 2);
            hpv_nonVax(j , dInd , :) = sum(vaxResult{n}.popVec(: , hpvInds_nonVax) , 2)...
                ./ sum(vaxResult{n}.popVec(: , hpvInds_tot) , 2);
        end
        
        
        %% ************************** SCREENING & VACCINATION FIGURES *******************************************************************************
        
        %% Total number of women overscreened annually (SUS, HPV:9v+ , HPV:9v-, CIN1:9v+ , CIN1:9v-) by HIV disease status
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            overScreenTotAnnualCum(j , dInd , : , :) = [annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , [1,7] , [1,7] , 1 , : , :),2),3),4),5),6),7)) ; ...
                annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , 2 , [1:2 , 7] , 1 , : , :),2),3),4),5),6),7)) ; ...
                annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , [1 , 7] , 2 , 1 , : , :),2),3),4),5),6),7)) ; ...
                (annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , 3 , [1:3 , 7] , 1 , : , :),2),3),4),5),6),7)) + annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , 2 , 3 , 1 , : , :),2),3),4),5),6),7))); ...
                annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , [1 , 7] , 3 , 1 , : , :),2),3),4),5),6),7))];
        end     
            
        %% Vaccine coverage overall
        % Overall
        vaxInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , [2 , 4] , 2 , 3 : 16 , 1 : risk));
        popInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 3 : 16 , 1 : risk));
        vaxCoverage(j , :) = sum(vaxResult{n}.popVec(: , vaxInds) , 2) ./ sum(vaxResult{n}.popVec(: , popInds) , 2);
        
        %% Vaccine coverage by age
        for a = 1 : age
            vaxInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , [2 , 4] , 2 , a , 1 : risk));
            popInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            vaxCoverageAge(j , a , :) = sum(vaxResult{n}.popVec(: , vaxInds) , 2) ./ sum(vaxResult{n}.popVec(: , popInds) , 2);
        end
        
        
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

%% Female population proportion by 5-year age groups over time vs. Statistics South Africa data (internal validation)
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
plot(1 : age , mean(squeeze(popPropF(: , 1 , :)),1) , 'k-' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , min(squeeze(popPropF(: , 1 , :)),[],1) , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , max(squeeze(popPropF(: , 1 , :)),[],1) , 'k--' , 'LineWidth' , 1.5);
hold all;

plot(1 : age , mean(squeeze(popPropF(: , 2 , :)),1) , 'b-' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , min(squeeze(popPropF(: , 2 , :)),[],1) , 'b--' , 'LineWidth' , 1.5);
hold all;
plot(1 : age , max(squeeze(popPropF(: , 2 , :)),[],1) , 'b--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Population proportion'); title('Female Age Distribution');
ylim([0 0.2]); grid on;
legend('(Statistics SA) Observed KZN, 2019' , 'Model, 2018: 25-sets mean' , ...
    'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum' , ...
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
    sgtitle('HIV prevalence by age over time');
end

%% HIV prevalence by gender over time vs. AHRI (calibration)
hivData(: , : , 1) = zeros(length(unique(hivPrevM_dObs(: ,1))) , 2);
hivData(: , : , 2) = zeros(length(unique(hivPrevM_dObs(: ,1))) , 2);
hivRaw(:,:,1) = hivPrevM_dObs(: , 4:5);
hivRaw(:,:,2) = hivPrevF_dObs(: , 4:5);
for i = 1 : length(unique(hivPrevM_dObs(: ,1)))
    for g = 1 : gender
        prev = (sumall(hivRaw(((i-1)*7+1):(i*7) , 1 , g)) ./ sumall(hivRaw(((i-1)*7+1):(i*7) , 2 , g)));
        hivData(i,1,g) = prev * 100;
        var = (prev*(1-prev)) / sumall(hivRaw(((i-1)*7+1):(i*7) , 2 , g));
        hivData(i,2,g) = (var ^ (1/2)) * 2 * 100;
    end
end

figure('DefaultAxesFontSize' , 18);
gen = {'Males' , 'Females'};
genFlipInd = {2 , 1};
for gInd = 1 : gender
    g = genFlipInd{gInd};
    subplot(1,2,g)
    errorbar(unique(hivPrevM_dObs(: ,1)) , hivData(:,1,g) , hivData(:,2,g) , ...
        'rs' , 'LineWidth' , 1.5);
    hold on;
    plot(monthlyTimespan , median(squeeze(hivPrev(: , g , :)),1)' , 'k-' , 'LineWidth' , 1.5);
    x2 = [monthlyTimespan , fliplr(monthlyTimespan)];
    inBetween = [max(squeeze(hivPrev(: , g , :)),[],1) , fliplr(min(squeeze(hivPrev(: , g , :)),[],1))];
    h = fill(x2 , inBetween , 'k');
    h.FaceAlpha = 0.3;
    h.LineStyle = '--';
    xlabel('Year'); ylabel('HIV Prevalence (%)'); title(gen{g});
    xlim([1990 2020]); ylim([0 100]);
    if g == 1 
        legend('(AHRI data request) Observed KZN, ages 15-49: mean, 2SD' , ...
            'Model, ages 15-49: 25-sets median' , 'Model, ages 15-49: 25-sets range' , ...
            'Location' , 'SouthEast');
    elseif g == 2
        legend('(AHRI data request) Observed KZN, ages 15-49: mean, 2SD' , ...
            'Model, ages 15-49: 25-sets median' , 'Model, ages 15-49: 25-sets range' , ...
            'Location' , 'SouthEast');
    end
end
sgtitle('HIV prevalence over time');

%% HIV Prevalence by age in all years with data vs. AHRI (calibration)
ageGroup = {'15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49'};
hivPrevYrVec = [2003 , 2005:2009];

gen = {'Male' , 'Female'};
for g = 1 : gender
    figure('DefaultAxesFontSize' , 13);
    for yrInd = 1 : length(hivPrevYrVec)
        yr = hivPrevYrVec(yrInd);
        
        % Calibration error bars
        hivM2(: , 1) = hivPrevM_dObs(((yrInd-1)*7+1):(yrInd*7) , 2) .* 100; % mean
        hivM2(: , 2) = (hivPrevM_dObs(((yrInd-1)*7+1):(yrInd*7) , 3).^(1/2)).*2 .* 100; % calibration SD
        hivF2(: , 1) = hivPrevF_dObs(((yrInd-1)*7+1):(yrInd*7) , 2) .* 100; % mean
        hivF2(: , 2) = (hivPrevF_dObs(((yrInd-1)*7+1):(yrInd*7) , 3).^(1/2)).*2 .* 100; % calibration SD
        
        subplot(2 , 3 , yrInd);
        hivPrevs = hivM2;
        hivModel = hivAgeM;
        if g == 2
            hivPrevs = hivF2;
            hivModel = hivAgeF;
        end
        hold all;            
        errorbar(1 : length(ageGroup) , hivPrevs(: , 1)' , hivPrevs(: , 2)' , 'rs' , 'LineWidth' , 1.5);
        hold all;
        plot(1 : length(ageGroup) , median((squeeze(hivModel(: , [4:10] , ((yr - startYear) * stepsPerYear +1)).*100)),1) , 'k-' , 'LineWidth' , 1.5);
        hold all;
        x2 = [[1:length(ageGroup)] , fliplr([1:length(ageGroup)])];
        inBetween = [max(squeeze(hivModel(: , [4:10] , ((yr - startYear) * stepsPerYear +1)).* 100),[],1) , ...
            fliplr(min(squeeze(hivModel(: , [4:10] , ((yr - startYear) * stepsPerYear +1)).* 100),[],1))];
        h = fill(x2 , inBetween , 'k');
        h.FaceAlpha = 0.3;
        h.LineStyle = '--';
        hold all;
        set(gca , 'xtickLabel' , ageGroup);
        set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
        xlabel('Age Group'); ylabel('HIV Prevalence (%)'); title(num2str(yr));
        if g == 1 
            ylim([0 100]);
        elseif g == 2
            ylim([0 100]);
        end
        grid on;
        legend('(AHRI data request) Observed KZN: mean, 2SD' , ...
            'Model: 25-sets median' , 'Model: 25-sets range')
    end
    sgtitle([gen{g} ,' HIV prevalence']);
end

%% Write crude HIV prevalence in women over time
firstYrInd = ((1982 - startYear)*stepsPerYear +1);
t82on = (1982:(lastYear-1))';
outputVec = [t82on , mean(squeeze(hivPrevW(: , (firstYrInd:stepsPerYear:end))) , 1)'];   
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'ART_comparative_modeling_outcome_templates_030421.xlsx'];
writematrix(outputVec , fname , 'Sheet' , 'HIVprev-Crude');

%% Proportion of total HIV+ population on ART and VS (denominator: CD4-eligible and ineligible)
figure('DefaultAxesFontSize' , 18);
subplot(1,2,1);
plot((artYr + 1) , maxRateM.*100 , 'ro');
hold all;
plot(monthlyTimespan , median(artCovM,1)' , 'k-' , 'LineWidth' , 1.5);
hold all;
x2 = [monthlyTimespan , fliplr(monthlyTimespan)];
inBetween = [max(squeeze(artCovM),[],1) , fliplr(min(squeeze(artCovM),[],1))];
h = fill(x2 , inBetween , 'k');
h.FaceAlpha = 0.3;
h.LineStyle = '--';
xlim([1990 2020]); ylim([0 100]);
xlabel('Year'); ylabel('MLWHIV on ART + VS (%)')
title('Males'); grid on;
legend('Observed KZN' , 'Model, ages 15-79: 25-sets median' , 'Model: 25-sets range');

subplot(1,2,2);
plot((artYr + 1) , maxRateF.*100 , 'ro');
hold all;
plot(monthlyTimespan , median(artCovF,1)' , 'k-' , 'LineWidth' , 1.5);
hold all;
x2 = [monthlyTimespan , fliplr(monthlyTimespan)];
inBetween = [max(squeeze(artCovF),[],1) , fliplr(min(squeeze(artCovF),[],1))];
h = fill(x2 , inBetween , 'k');
h.FaceAlpha = 0.3;
h.LineStyle = '--';
xlim([1990 2020]); ylim([0 100]);
xlabel('Year'); ylabel('WLWHIV on ART + VS (%)')
grid on; title('Females');
legend('Observed KZN' , 'Model, ages 10-79: 25-sets median' , 'Model: 25-sets range');
sgtitle('ART + VS Coverage');

%% Write crude ART coverage in women over time
firstYrInd = ((1982 - startYear)*stepsPerYear +1);
t82on = (1982:(lastYear-1))';
outputVec = [t82on , mean(squeeze(artCovW(: , (firstYrInd:stepsPerYear:end))) , 1)'];   
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'ART_comparative_modeling_outcome_templates_030421.xlsx'];
writematrix(outputVec , fname , 'Sheet' , 'ARTprev-Crude');

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
title('Proportion on ART + VS by age');
ylim([0 1]); grid on;
legend('Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');

%% Proportion HIV-negative males circumcised by broad age groups over time
circPropYr_obs = vmmcYr;
circProp_obs = vmmcRate' .* 100;

figure('DefaultAxesFontSize' , 18);
plot(circPropYr_obs , circProp_obs , 'o');
hold on;
set(gca,'ColorOrderIndex',1)
p = plot(monthlyTimespan , squeeze(median(circProp,1)) , '-' , 'LineWidth' , 1.5);
hold all;
x2 = [monthlyTimespan , fliplr(monthlyTimespan)];
inBetween = [squeeze(max(squeeze(circProp),[],1)) , fliplr(squeeze(min(squeeze(circProp),[],1)))];
colorVecP = get(p,'Color');
h1 = fill(x2 , inBetween(1,:) , colorVecP{1});
h1.FaceAlpha = 0.3;
h1.LineStyle = '--';
h2 = fill(x2 , inBetween(2,:) , colorVecP{2});
h2.FaceAlpha = 0.3;
h2.LineStyle = '--';
h3 = fill(x2 , inBetween(3,:) , colorVecP{3});
h3.FaceAlpha = 0.3;
h3.LineStyle = '--';
h4 = fill(x2 , inBetween(4,:) , colorVecP{4});
h4.FaceAlpha = 0.3;
h4.LineStyle = '--';
xlim([1990 2020]); ylim([0 100]);
xlabel('Year'); ylabel('HIV-Negative Males Circumcised(%)')
grid on;
legend('Observed KZN, ages 15-19' , ...
    'Observed KZN, ages 20-24' , 'Observed KZN, ages 25-49' , ...
    'Observed KZN, ages 50+' , ...
    'Model, ages 15-19: 25-sets median' , ...
    'Model, ages 20-24: 25-sets median' , ...
    'Model, ages 25-49: 25-sets median'  , ...
    'Model, ages 50+: 25-sets median' , ...
    'Model: 25-sets range' , ...
    'Model: 25-sets range' , ...
    'Model: 25-sets range' , ...
    'Model: 25-sets range');


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
legend('(McDonald, 2014) Observed Cape Town: mean, 2SD' , ...
    'Model, 2002: 25-sets mean' , 'Model, 2002: 25-sets minimum' , 'Model, 2002: 25-sets maximum' , ...
    'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');title('HIV+'); grid on;

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
title('HIV-');
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
title('All');
grid on;
sgtitle('Female hrHPV Prevalence (includes CIN) by HIV status');

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
title('HIV+');
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
title('HIV-');
grid on;

subplot(3,1,3);
plot(1 : length(meanObs) , mean(hpv_hivMtot(: , : , 1),1) , 'b-' , ...
    1 : length(meanObs) , min(hpv_hivMtot(: , : , 1),[],1) , 'b--' , ...
    1 : length(meanObs) , max(hpv_hivMtot(: , : , 1),[],1) , 'b--' , 'LineWidth' , 1.5);
set(gca , 'xtick' , [1 : length(ageGroup)] , 'xtickLabel' , ageGroup);
legend('Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum');
xlabel('Age Group'); ylabel('hrHPV Prevalence'); ylim([0 1]);
title('General');
grid on;
sgtitle('Male hrHPV Prevalence by HIV status');

%% Write crude HPV prevalence by HIV status over time
firstYrInd = ((1982 - startYear)*stepsPerYear +1);
t82on = (1982:(lastYear-1))';
outputVec = [];
for dInd = 1 : diseaseVecLength_ccInc    
    outputVecA = [];
    for a = 1 : age
        outputVecA = [outputVecA , mean(squeeze(hpv_hivAgeW(: , dInd , a , (firstYrInd:stepsPerYear:end))) , 1)'];
    end
    outputVec = [outputVec ; ...
        [t82on , ones(length(t82on),1).*dInd , ...
        mean(squeeze(hpv_hivW(: , dInd , (firstYrInd:stepsPerYear:end))) , 1)' , ...
        outputVecA]];   
end
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'ART_comparative_modeling_outcome_templates_030421.xlsx'];
writematrix(outputVec , fname , 'Sheet' , 'HPVprev-Crude');

%% Write age-standardized HPV prevalence by HIV status over time
% % Note: the age-standardization process shifts the prevalence of the
% % last modelled age group to the next age group in the following year.
% % However, prevalence is NaN prior to HIV introduction in the
% % HIV-positive no ART group, and NaN prior to ART introduction in the
% % HIV-positive ART group. Since we have four age groups past the 16 we
% % model, a NaN value is present for four years past the introduction of
% % HIV/ART, leading to a NaN value for summed HPV infected during these 
% % years. We therefore lack data in this four-year interval in the
% % saved/plotted results.
% worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
%     247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
%     23477 9261 2155];
% 
% diseaseLabels = {'Tot (ICC)' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'};
% firstYrRange2 = (lastYear-1) - 1982;
% t82on = (1982:(lastYear-1))';
% outputVec = [];
% for dInd = 1 : length(diseaseLabels)
%     hpv_hivAgeW_dis = squeeze(hpv_hivAgeW(: , dInd , : , (1:stepsPerYear:end)));
% 
%     numHpvTot = zeros(size(hpv_hivAgeW_dis,1) , 1 , size(hpv_hivAgeW_dis,3));       
%     for aInd = 1:age+4
%         a = aInd;
%         if aInd >= age
%             a = age;
%         end
% 
%         if aInd <= age    
%             numHpv = hpv_hivAgeW_dis(: , a , :) .* worldStandard_WP2015(aInd);
%             if (a < 3)
%                 numHpv = zeros(size(hpv_hivAgeW_dis,1) , 1 , size(hpv_hivAgeW_dis,3));
%             end
%         elseif aInd > age
%             numHpv = hpv_hivAgeW_dis(: , a , :);
%             numHpv = cat(3 , (ones(size(numHpv,1),1,aInd-a).*numHpv(:,1,1)) , numHpv(: , 1 ,1:end-(aInd-a)));
%             numHpv = numHpv .* worldStandard_WP2015(aInd);
%         end
%         numHpvTot = numHpvTot + numHpv;
%     end
%     hpvPrevTot = numHpvTot ./ (sum(worldStandard_WP2015(1:age+4)));
%   
%     outputVec = [outputVec ; ...
%         [t82on , ones(length(t82on),1).*dInd , ...
%         squeeze(mean(squeeze(hpvPrevTot(: , : , (end-firstYrRange2):end)) , 1))' , ...
%         squeeze(min(squeeze(hpvPrevTot(: , : , (end-firstYrRange2):end)) , [] , 1))' , ...
%         squeeze(max(squeeze(hpvPrevTot(: , : , (end-firstYrRange2):end)) , [] , 1))' , ...
%         squeeze(hpvPrevTot(: , : , (end-firstYrRange2):end))']];    
% end
% 
% fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
%     'ART_comparative_modeling_outcome_templates_020221-fullData.xlsx'];  
% writematrix(outputVec , fname , 'Sheet' , 'HPVprev-AS')


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
    'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum' , ...
    'Location' , 'northwest');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence')
title('WLWHIV')
ylim([0 0.3])
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
    'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum' , ...
    'Location' , 'northwest');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence')
title('HIV-negative women')
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
    'Model, 2018: 25-sets mean' , 'Model, 2018: 25-sets minimum' , 'Model, 2018: 25-sets maximum' , ...
    'Location' , 'northwest');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence')
title('General')
ylim([0 0.2])
grid on;
sgtitle('CIN2/3 Prevalence by HIV status');


%% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

%% Cervical cancer incidence in 2018 by age vs. Globocan 2018 data (calibration)
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

% Load adjusted Globocan 2018 rates for KZN
file = [pwd , '/Config/Reweighted_GlobocanCC_rates.xlsx'];
ccInc2018adjKZN(:,1) = xlsread(file , 'CC rates' , 'AB4:AB15');

% Calibration error bars
meanObs = ccInc2018_dObs(: , 2);
sdevObs = (ccInc2018_dObs(: , 3).^(1/2)).*2;

figure('DefaultAxesFontSize' , 18); 
% Plot observed data
errorbar(4 : age-1 , meanObs , sdevObs , 'rs' , 'LineWidth' , 1.5);
hold all;
plot(4 : age-1 , ccInc2018adjKZN , 'r*');
hold all;
% General
plot(1 : age , median(squeeze(ccIncAge(: , 3 , :)),1) , 'k-' , 'LineWidth' , 1.5);
hold all;
x2 = [[1 : age] , fliplr([1 : age])];
inBetween = [max(squeeze(ccIncAge(: , 3 , :)),[],1) , fliplr(min(squeeze(ccIncAge(: , 3 , :)),[],1))];
h = fill(x2 , inBetween , 'k');
h.FaceAlpha = 0.3;
h.LineStyle = '--';
xlabel('Age Group'); ylabel('Cervical cancer incidence per 100K');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
ylim([0 300]); grid on;
title(['Cervical Cancer Incidence in 2018 by age']);
legend('(Globocan, 2018) Observed SA: mean, 2SD' , 'Estimated KZN, adjusted Globocan 2018' , ...
    'Model, general: 25-sets median' , ...
    'Model: 25-sets range' , 'Location' , 'NorthWest');

%% Cervical cancer incidence over time by HIV status
figure;   
% Plot observed data
plot([2005 , 2012 , 2018] , [26.6 , 31.7 , 44.4] , 'ro');
hold all;
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
xlim([1980 2120]); ylim([0 250]); grid on;
title(['Cervical Cancer Incidence over time, ages 1-79']);
legend('(Globocan, 2005, 2012, 2018) Observed SA crude ??????' , ...
    'Model, general: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, HIV-negative: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, WLWHIV untreated: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, WLWHIV on ART: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Location' , 'northwest');

%% Cervical cancer incidence over time
figure('DefaultAxesFontSize' , 18);  
% Plot observed data
plot([2005 , 2012 , 2018] , [26.6 , 31.7 , 44.4] , 'ro');
hold all;
% General
plot(annualTimespan , median(ccIncTime,1) , 'k-' , 'LineWidth' , 1.5);
hold all;
x2 = [annualTimespan , fliplr(annualTimespan)];
inBetween = [max(squeeze(ccIncTime),[],1) , fliplr(min(squeeze(ccIncTime),[],1))];
h = fill(x2 , inBetween , 'k');
h.FaceAlpha = 0.3;
h.LineStyle = '--';
xlabel('Time'); ylabel('Cervical cancer incidence per 100K');
xlim([1990 2020]); ylim([0 100]); grid on;
title(['Cervical Cancer Incidence over Time']);
legend('(Globocan, 2005, 2012, 2018) Observed SA crude ??????' , ...
    'Model, ages 1-79: 25-sets median' , 'Model, ages 1-79: 25-sets range' , ...
    'Location' , 'northwest');

%% Write age-standardized cervical cancer incidence rates by HIV status over time (2020-2120) into existing template
% Note: the age-standardization process shifts the incidence rate of the
% last modelled age group to the next age group in the following year.
% However, CC incidence is NaN prior to HIV introduction in the
% HIV-positive no ART group, and NaN prior to ART introduction in the
% HIV-positive ART group. Since we have four age groups past the 16 we
% model, a NaN value is present for four years past the introduction of
% HIV/ART, leading to a NaN value for summed incidence during these 
% years. We therefore lack data in this four-year interval in the
% saved/plotted results.
fac = 10 ^ 5;
worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
    247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
    23477 9261 2155];

diseaseLabels = {'Tot (ICC)' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'};
firstYrRange = (lastYear-1) - currYear;
firstYrRange2 = (lastYear-1) - 1982;
t82on = (1982:(lastYear-1))';
outputVec = [];
for dInd = 1 : length(diseaseLabels)
    ccIncHivAgeTime_dis = squeeze(ccIncHivAgeTime(: , dInd , : , :));

    ccIncRefTot = zeros(size(ccIncHivAgeTime_dis,1) , 1 , size(ccIncHivAgeTime_dis,3));       
    for aInd = 1:age+4
        a = aInd;
        if aInd >= age
            a = age;
        end

        if aInd <= age    
            ccIncRef = ccIncHivAgeTime_dis(: , a , :) .* worldStandard_WP2015(aInd);
            if (a < 3)
                ccIncRef = zeros(size(ccIncHivAgeTime_dis,1) , 1 , size(ccIncHivAgeTime_dis,3));
            end
        elseif aInd > age
            ccIncRef = ccIncHivAgeTime_dis(: , a , :);
            ccIncRef = cat(3 , (ones(size(ccIncRef,1),1,aInd-a).*ccIncRef(:,1,1)) , ccIncRef(: , 1 ,1:end-(aInd-a)));
            ccIncRef = ccIncRef .* worldStandard_WP2015(aInd);
        end
        ccIncRefTot = ccIncRefTot + ccIncRef;
    end
    ccInc = ccIncRefTot ./ (sum(worldStandard_WP2015(1:age+4)));

%     fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
%         'SA_screening_S' , fileKeyNums{n} , '.xlsx'];  
%     writematrix([squeeze(median(squeeze(ccInc(: , : , (end-firstYrRange):end)) , 1))' , ...
%         squeeze(min(squeeze(ccInc(: , : , (end-firstYrRange):end)) , [] , 1))' , ...
%         squeeze(max(squeeze(ccInc(: , : , (end-firstYrRange):end)) , [] , 1))' , ...
%         squeeze(ccInc(: , : , (end-firstYrRange):end))'] , ...
%         fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B3') 
    
    outputVec = [outputVec; ...
        [t82on , ones(length(t82on),1).*dInd , ...
        squeeze(mean(squeeze(ccInc(: , : , (end-firstYrRange2):end)) , 1))' , ...
        squeeze(min(squeeze(ccInc(: , : , (end-firstYrRange2):end)) , [] , 1))' , ...
        squeeze(max(squeeze(ccInc(: , : , (end-firstYrRange2):end)) , [] , 1))' , ...
        squeeze(ccInc(: , : , (end-firstYrRange2):end))']];
end
% fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
%     'ART_comparative_modeling_outcome_templates_020221-fullData.xlsx'];  
% writematrix(outputVec , fname , 'Sheet' , 'ICC-AS')

%% Write crude cervical cancer incidence by HIV status and age over time
firstYrInd = ((1982 - startYear) +1);
t82on = (1982:(lastYear-1))';
outputVec = [];
for dInd = 1 : diseaseVecLength_ccInc
    outputVecA = [];
    for a = 1 : age
        outputVecA = [outputVecA , mean(squeeze(ccIncHivAgeTime(: , dInd , a , (firstYrInd:end))) , 1)'];
    end
    outputVec = [outputVec ; ...
        [t82on , ones(length(t82on),1).*dInd , ...
        mean(squeeze(ccIncHivTime(: , dInd , (firstYrInd:end))) , 1)' , ...
        outputVecA]];   
end
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'ART_comparative_modeling_outcome_templates_030421.xlsx'];
writematrix(outputVec , fname , 'Sheet' , 'ICC-Crude');

%% Write cumulative cervical cancer cases by by HIV status over time (2020-2120) into existing template
% diseaseLabels = {'Tot (CCC)' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)'};
% for dInd = 1 : length(diseaseLabels)
%     fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
%         'SA_screening_S' , fileKeyNums{n} , '.xlsx'];
%     writematrix([squeeze(median(squeeze(ccCumHivTime(: , dInd , :)) , 1))' , ...
%         squeeze(min(squeeze(ccCumHivTime(: , dInd , :)) , [] , 1))' , ...
%         squeeze(max(squeeze(ccCumHivTime(: , dInd , :)) , [] , 1))' , ...
%         squeeze(ccCumHivTime(: , dInd , :))'] , ...
%         fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B3')   
% end

%% Write annual cervical cancer cases by HIV status and age over time
firstYrInd = ((1982 - startYear) +1);
t82on = (1982:(lastYear-1))';
outputVec = [];
for dInd = 1 : diseaseVecLength_ccInc
    outputVecA = [];
    for a = 1 : age
        outputVecA = [outputVecA , mean(squeeze(ccAnlHivAgeTime(: , dInd , a , (firstYrInd:end))) , 1)'];
    end
    outputVec = [outputVec ; ...
        [t82on , ones(length(t82on),1).*dInd , ...
        mean(squeeze(ccAnlHivTime(: , dInd , (firstYrInd:end))) , 1)' , ...
        outputVecA]];   
end
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'ART_comparative_modeling_outcome_templates_030421.xlsx'];
writematrix(outputVec , fname , 'Sheet' , 'CCC-Crude');


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
sgtitle('Type distribution by state (coinfections grouped as 9v-type HPV), ages 0-79');


%% ************************** SCREENING & VACCINATION FIGURES *******************************************************************************

%% Total number of women overscreened annually (SUS, HPV:9v+ , HPV:9v-, CIN1:9v+ , CIN1:9v-) by HIV disease status
% diseaseLabels = {'Tot (OS)' , 'HIV- (OS)' , 'HIV+ (OS)' , 'HIV+ no ART (OS)' , 'HIV+ ART (OS)'};
% for dInd = 1 : length(diseaseLabels)
%     fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
%         'SA_screening_S' , fileKeyNums{n} , '.xlsx'];
%     writematrix([[squeeze(median(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 1 , :)),2) , 1))' , ...
%             squeeze(min(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 1 , :)),2) , [] , 1))' , ...
%             squeeze(max(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 1 , :)),2) , [] , 1))' , ...
%             cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 1 , :)),2)'] ; ...
%         [squeeze(median(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 2 , :)),2) , 1))' , ...
%             squeeze(min(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 2 , :)),2) , [] , 1))' , ...
%             squeeze(max(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 2 , :)),2) , [] , 1))' , ...
%             cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 2 , :)),2)'] ; ...
%         [squeeze(median(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 3 , :)),2) , 1))' , ...
%             squeeze(min(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 3 , :)),2) , [] , 1))' , ...
%             squeeze(max(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 3 , :)),2) , [] , 1))' , ...
%             cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 3 , :)),2)'] ; ...
%         [squeeze(median(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 4 , :)),2) , 1))' , ...
%             squeeze(min(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 4 , :)),2) , [] , 1))' , ...
%             squeeze(max(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 4 , :)),2) , [] , 1))' , ...
%             cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 4 , :)),2)'] ; ...
%         [squeeze(median(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 5 , :)),2) , 1))' , ...
%             squeeze(min(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 5 , :)),2) , [] , 1))' , ...
%             squeeze(max(cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 5 , :)),2) , [] , 1))' , ...
%             cumsum(squeeze(overScreenTotAnnualCum(: , dInd , 5 , :)),2)']] , ...
%         fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'C3')   
% end

%% Vaccine coverage overall
figure;   
plot(midAnnualTimespan , mean(vaxCoverage(: , (4 : stepsPerYear : end)),1) , 'k-' , ...
    midAnnualTimespan , min(vaxCoverage(: , (4 : stepsPerYear : end)),[],1) , 'k--' , ...
    midAnnualTimespan , max(vaxCoverage(: , (4 : stepsPerYear : end)),[],1) , 'k--' , 'LineWidth' , 1.5);
xlabel('Time'); ylabel('Vaccine coverage');
xlim([2020 2100]); ylim([0 1]); grid on;
title(['Vaccine coverage']);
legend('Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Location' , 'northwest');

%% Vaccine coverage by age
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

figure;    
for a = 1 : age
    plot(midAnnualTimespan , mean(squeeze(vaxCoverageAge(: , a , (4 : stepsPerYear : end))),1)' , '-' , 'LineWidth' , 1.5);
    hold all;
end
xlabel('Year'); ylabel('Vaccine coverage');
xlim([2020 2100]); ylim([0 1]); grid on;
title(['Vaccine coverage by age']);
legend('Model, 0-4: 25-sets mean' , ...
    'Model, 5-9: 25-sets mean' , ...
    'Model, 10-14: 25-sets mean' , ...
    'Model, 15-19: 25-sets mean' , ...
    'Model, 20-24: 25-sets mean' , ...
    'Model, 25-29: 25-sets mean' , ...
    'Model, 30-34: 25-sets mean' , ...
    'Model, 35-39: 25-sets mean' , ...
    'Model, 40-44: 25-sets mean' , ...
    'Model, 45-49: 25-sets mean' , ...
    'Model, 50-54: 25-sets mean' , ...
    'Model, 55-59: 25-sets mean' , ...
    'Model, 60-64: 25-sets mean' , ...
    'Model, 65-69: 25-sets mean' , ...
    'Model, 70-74: 25-sets mean' , ...
    'Model, 75-59: 25-sets mean');

