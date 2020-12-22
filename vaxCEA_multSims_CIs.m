function [] = vaxCEA_multSims_CIs(vaxResultInd , sceNum)

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
% fileInds = {'6_2' , '6_3' , '6_4' , '6_6' , ...
%     '6_7' , '6_8' , '6_9' , '6_10' , '6_11' , '6_12' , '6_13' , ...
%     '6_14' , '6_15' , '6_16' , '6_17' , '6_18' , '6_19' , '6_20' , ...
%     '6_21' , '6_22' , '6_23' , '6_24' , '6_25'};    % 22Apr20Ph2V11, t=6
% fileInds = {'5_1' , '5_3' , '5_4' , '5_19' , '5_21' , '5_22' , '5_24'};    % 22Apr20Ph2V11, t=5
% fileInds = {'11_1' , '11_2' , '11_3' , '11_4' , '11_5' , '11_6' , ...
%     '11_7' , '11_8' , '11_9' , '11_10' , '11_11' , '11_12' , '11_13' , ...
%     '11_14' , '11_15' , '11_16' , '11_17' , '11_18' , '11_19' , '11_20' , ...
%     '11_21' , '11_22' , '11_23' , '11_24' , '11_25'};   % DO ART, 22Apr20Ph2V2, t=11
% fileInds = {'7_1' , '7_2' , '7_3' , '7_4' , '7_5' , '7_6' , '7_7' , '7_8' , ...
%     '7_9' , '7_10' , '7_11' , '7_12' , '7_13' , '7_14' , '7_15' , '7_16' , ...
%     '7_17' , '7_18' , '7_19' , '7_20' , '7_21' , '7_22' , '7_23' , '7_24' , '7_25'};   % DO ART, 22Apr20Ph2V2, t=7
nRuns = length(fileInds);

% Initialize model output plots
% Timespans
monthlyTimespan = [startYear : (1/6) : lastYear];
monthlyTimespan = monthlyTimespan(1 : end-1);
annualTimespan = [startYear : lastYear-1];
futAnnualTimespan = [2019 : lastYear-1];
midAnnualTimespan = [(startYear+(3/stepsPerYear)) : ((lastYear-1)+(3/stepsPerYear))];
screenAnnualTimespan = [(2020+(3/stepsPerYear)) : ((lastYear-1)+(3/stepsPerYear))];
screenMonthlyTimespan = [2020 : (1/6) : lastYear];
screenMonthlyTimespan = screenMonthlyTimespan(1 : end-1);
% Total population size
popSize = zeros(nRuns , length(monthlyTimespan));
popSizeAgeF = zeros(nRuns , 5 , age , length(monthlyTimespan));
% Population age distribution
popYearVec = [2018 2100];
popYearVecLength = length(popYearVec);
popPropF = zeros(nRuns , length(popYearVec) , age);
popYearVec2 = [2019 2100];
popYearVecLength2 = length(popYearVec2);
ageVec_cPopDist = {3 , [4:5] , [6:7] , [8:9] , [10:11] , [12:13] , [14:15]};
ageVecLength_cPopDist = length(ageVec_cPopDist);
popPropBroadC = zeros(nRuns , length(popYearVec2) , ageVecLength_cPopDist);
% Female risk distribution
popRiskDistF = zeros(nRuns , risk , length(monthlyTimespan));
popRiskDistHivF = popRiskDistF;
popRiskDistCinF = popRiskDistF;
% Risk distribution by HIV status and gender
popRiskDistHivProp = zeros(nRuns , risk , gender , length(monthlyTimespan));
% HIV prevalence
hivAgeM = zeros(nRuns , age , length(monthlyTimespan));
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
hpvYearVec_orig = [2002 2018];
hpvYearVec2018_orig = [2018];
hpvYearVecLength = length(hpvYearVec_orig);
hpvYearVec2018Length = length(hpvYearVec2018_orig);
hpv_hiv = zeros(nRuns , age , 2);
hpv_hivNeg = hpv_hiv;
hpv_hivTot = zeros(nRuns , age , 1);
diseaseVec_fHpv = {[1 : disease] , [1 : 2] , [3 : 8] , [3 : 7] , 8};
diseaseVecLength_fHpv = length(diseaseVec_fHpv);
ageVec_fHPV = {3 , [4:5] , [6:7] , [8:9] , [10:11] , [12:13] , [14:15] , 16 };
ageVecLength_fHPV = length(ageVec_fHPV);
hpv_hivTot2020 = zeros(nRuns , diseaseVecLength_fHpv , ageVecLength_fHPV);
% Male HPV prevalence
hpvYearVecMale = [2008 2018];
hpvYearVecMaleLength = length(hpvYearVecMale);
hpv_hivM = zeros(nRuns , 4 , 2);
hpv_hivMNeg = hpv_hivM;
hpv_hivMtot = hpv_hivM;
ageVec_mHPV = {[4:5],[6:7],[8:9],[10:13]};
ageVecLength_mHPV = length(ageVec_mHPV);
% HPV prevalence over time
hpv_hivTimeF = zeros(nRuns , length(monthlyTimespan));
hpv_hivNegTimeF = hpv_hivTimeF;
hpv_time = zeros(nRuns , 2 , length(monthlyTimespan));
% HPV incidence
hpvIncTime = zeros(nRuns , length(annualTimespan));
hpvIncTimeNeg = hpvIncTime;
hpvIncTimePos = hpvIncTime;
hpvIncTimeArt = hpvIncTime;
hpvIncTimePosAll = hpvIncTime;
hpvIncHivRiskAgeTime = zeros(nRuns , 5 , risk , age , length(annualTimespan));
% CIN2/3 prevalence
cinPosAge = zeros(nRuns , 2 , age);
cinNegAge = cinPosAge;
cinGenAge = cinPosAge;
cinGenTime = zeros(nRuns , length(monthlyTimespan));
cinPosTime = cinGenTime;
cinNegTime = cinGenTime;
cin_hivTot2020 = zeros(nRuns , diseaseVecLength_fHpv , ageVecLength_fHPV);
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
ccCumHivAgeTime = zeros(nRuns , 5 , age , length(futAnnualTimespan));
diseaseVec_ccInc = {[1 : disease] , [1 : 2] , [3 : 8] , [3 : 7] , 8};
diseaseVecLength_ccInc = length(diseaseVec_ccInc);
% Prevalence ratios
diseaseVec_hpv = {[1 : 2] , [3 : 8] , [3 : 7] , 8};
hpvPrevYearVec = [2005 2018 2019];
diseaseVecLength_hpv = length(diseaseVec_hpv);
hpvPrevYearVecLength = length(hpvPrevYearVec);
hpv_prev_ratios = zeros(nRuns , gender , 4 , 2);
cin_prev_ratios = zeros(nRuns , 4 , 2);
cc_prev_ratios = zeros(nRuns , 4 , 2);
cc_inc_ratios = zeros(nRuns , 4 , 2);
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
newScreenTime = zeros(nRuns , length(screenAnnualTimespan));
screenCovTime = zeros(nRuns , length(screenMonthlyTimespan));
screenTotAnnual = zeros(nRuns , 5 , length([(1925+(3/stepsPerYear)) : ((lastYear-1)+(3/stepsPerYear))]));
vaxCoverage = zeros(nRuns , length(monthlyTimespan));
vaxCoverageAge = zeros(nRuns , age , length(monthlyTimespan));
vaxTotAge = zeros(nRuns , age , 5 , length(monthlyTimespan));

resultsDir = [pwd , '\HHCoM_Results\'];
fileKey = {'sim1' , 'sim2' , 'sim0'};
n = vaxResultInd;
baseFileName = ['22Apr20Ph2V11_noBaseVax_baseScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_WHO-SCES' , sceNum , '_'];
loopSegments = {0 , round(nRuns/2) , nRuns};
loopSegmentsLength = length(loopSegments);
for k = 1 : loopSegmentsLength-1
    parfor j = loopSegments{k}+1 : loopSegments{k+1}
        % Load results
        pathModifier = [baseFileName , fileInds{j}]; % ***SET ME***: name for simulation output file
        nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
        curr = load([pwd , '/HHCoM_Results/toNow_22Apr20Ph2V11_noBaseVax_baseScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_' , fileInds{j}]); % ***SET ME***: name for historical run output file 

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
        vaxResult{n}.newScreen = [vaxResult{n}.newScreen(1 : end , : , : , : , : , : , : , : , :)]; %curr.newScreen(1 : end , : , : , : , : , : , : , : , :); vaxResult{n}.newScreen(2 : end , : , : , : , : , : , : , : , :)]; 
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
        
        %% Population proportion by broad age groups over time
        for t = 1 : popYearVecLength2
            for aV = 1 : ageVecLength_cPopDist
                aGroup = ageVec_cPopDist{aV};
                popAgeF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 1 : gender , aGroup , 1 : risk));
                popTotF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 1 : gender , 3 : 15 , 1 : risk));
                popPropBroadC(j , t , aV) = (sum(vaxResult{n}.popVec(((popYearVec2(t) - startYear) * stepsPerYear +1) , popAgeF),2) ./ ...
                    sum(vaxResult{n}.popVec(((popYearVec2(t) - startYear) * stepsPerYear +1) , popTotF),2)) .* 100;
            end
        end
        
        %% Female total population size by 5-year age groups over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                popAgeF = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                popSizeAgeF(j , dInd , a , :) = sum(vaxResult{n}.popVec(: , popAgeF),2);
            end
        end
        
        %% Risk distribution of HIV-negative women with CIN3 over time
        for rInd = 1 : risk
            r = rInd;
            popRiskF = unique([toInd(allcomb(1 : 2 , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 4 : age , r)); toInd(allcomb(1 : 2 , 1 : viral , ...
                [1 : 5 , 7] , 5 , 1 , 1 : intervens , 2 , 4 : age , r))]);
            popAllF = unique([toInd(allcomb(1 : 2 , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 4 : age , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
                [1 : 5 , 7] , 5 , 1 , 1 : intervens , 2 , 4 : age , 1 : risk))]);
            popRiskDistCinF(j , rInd , :) = sum(vaxResult{n}.popVec(: , popRiskF),2) ./ sum(vaxResult{n}.popVec(: , popAllF),2);
        end
        
        %% Risk distribution of HIV-negative women over time
        for rInd = 1 : risk
            r = rInd;
            popRiskF = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , 5 : 7 , r));
            popAllF = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , 5 : 7 , 1 : risk)); 
            popRiskDistF(j , rInd , :) = sum(vaxResult{n}.popVec(: , popRiskF),2) ./ sum(vaxResult{n}.popVec(: , popAllF),2);
        end
        
        %% Risk distribution of HIV-positive women over time
        for rInd = 1 : risk
            r = rInd;
            popRiskF = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , 4 : age , r));
            popAllF = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , 4 : age , 1 : risk)); 
            popRiskDistHivF(j , rInd , :) = sum(vaxResult{n}.popVec(: , popRiskF),2) ./ sum(vaxResult{n}.popVec(: , popAllF),2);
        end
        
        %% Proportion of each risk group that is HIV-positive over time
        for g = 1 : gender
            for rInd = 1 : risk
                r = rInd;
                popRiskF = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , g , 4 : age , r));
                popAllF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , g , 4 : age , r)); 
                popRiskDistHivProp(j , rInd , g , :) = sum(vaxResult{n}.popVec(: , popRiskF),2) ./ sum(vaxResult{n}.popVec(: , popAllF),2);
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
        
        %% HIV-associated deaths by gender over time
        hivDeathsM(j , :) = annlz(sum(sum(vaxResult{n}.hivDeaths(: , : , 1 , :), 2), 4));
        hivDeathsF(j , :) = annlz(sum(sum(vaxResult{n}.hivDeaths(: , : , 2 , :), 2), 4));
        
        %% Proportion of total HIV+ population on ART and VS (denominator: CD4-eligible and ineligible)
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
        
        %% Female HPV Prevalence by broad age groups 2019
        for dInd = 1 : diseaseVecLength_fHpv
            d = diseaseVec_fHpv{dInd};                 
            for aV = 1 : ageVecLength_fHPV
                aGroup = ageVec_fHPV{aV};
                hpvInds_hivTot = unique([toInd(allcomb(d , 1 : viral , 2 : 6 , 1 : hpvNonVaxStates , ...
                    1 : 3 , 1 : intervens , 2 , aGroup , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                    1 : hpvVaxStates , 2 : 6 , 1 : 3 , 1 : intervens , 2 , aGroup , 1 : risk))]);
                ageInds_hivTot = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , aGroup , 1 : risk));
                hpv_hivTot2020(j , dInd , aV) = sum(vaxResult{n}.popVec((2019 - startYear) * stepsPerYear +1 , hpvInds_hivTot))...
                    ./ sum(vaxResult{n}.popVec((2019 - startYear) * stepsPerYear +1 , ageInds_hivTot));
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
        
        %% Female HPV Prevalence over time by HIV status
        hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 4 : 13 , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk));
        hpv_hivTimeF(j , :) = sum(vaxResult{n}.popVec(: , hpvInds) , 2) ./ sum(vaxResult{n}.popVec(: , ageInds) , 2);
    
        hpvInds_hivNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 4 : 13 , 1 : risk))]);
        ageInds_hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk));
        hpv_hivNegTimeF(j , :) = sum(vaxResult{n}.popVec(: , hpvInds_hivNeg) , 2) ./ sum(vaxResult{n}.popVec(: , ageInds_hivNeg) , 2);
    
        %% HPV Prevalence over time by sex       
        for g = 1 : gender
            hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , g , 4 : 13 , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
                [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , g , 4 : 13 , 1 : risk))]);
            popInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , g , 4 : 13 , 1 : risk));
            hpv_time(j , g , :) = sum(vaxResult{n}.popVec(: , hpvInds) , 2)...
                ./ sum(vaxResult{n}.popVec(: , popInds) , 2);
        end
        
        %% HPV prevalence ratios in 2005, 2018, 2019
        for i = 1 : hpvPrevYearVecLength
            yr = hpvPrevYearVec(i); 
            for g = 1 : gender
                for dInd = 1 : diseaseVecLength_hpv
                    d = diseaseVec_hpv{dInd};
                    hpvInds = unique([toInd(allcomb(d , 1 : viral , 2 : 6 , 1 : hpvNonVaxStates , ...
                        1 : 3 , 1 : intervens , g , 4 : 13 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                        1 : hpvVaxStates , 2 : 6 , 1 : 3 , 1 : intervens , g , 4 : 13 , 1 : risk))]);
                    popInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                        1 : endpoints , 1 : intervens , g , 4 : 13 , 1 : risk));
                    hpv_prev_ratios(j , g , dInd , i) = sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds) , 2)...
                        ./ sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , popInds) , 2);
                end
            end
        end
    
        %% HPV incidence over time
        fac = 100;
        % General population
        allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        hpvIncTime(j , :) = ...
            ((annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvVax(: , 2 , 1 : disease , : , : , :),2),3),4),5),6)) + ...
            annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvVax(: , 2 , 1 : disease , : , : , :),2),3),4),5),6)) + ...
            annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvNonVax(: , 2 , 1 : disease , : , : , :),2),3),4),5),6)) ...
            ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvNonVax(: , 2 , 1 : disease , : , : , :),2),3),4),5),6))
            (annlz(sum(vaxResult{n}.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-negative
        allFneg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        hpvIncTimeNeg(j , :) = ...
            ((annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvVax(: , 2 , 1 : 2 , : , : , :),2),3),4),5),6)) + ...
            annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvVax(: , 2 , 1 : 2 , : , : , :),2),3),4),5),6)) + ...
            annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvNonVax(: , 2 , 1 : 2 , : , : , :),2),3),4),5),6)) ... 
            ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvNonVax(: , 2 , 1 : 2 , : , : , :),2),3),4),5),6))
            (annlz(sum(vaxResult{n}.popVec(: , allFneg) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-positive untreated
        allFpos = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        hpvIncTimePos(j , :) = ...
            ((annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvVax(: , 2 , 3 : 7 , : , : , :),2),3),4),5),6)) + ...
            annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvVax(: , 2 , 3 : 7 , : , : , :),2),3),4),5),6)) + ...
            annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvNonVax(: , 2 , 3 : 7 , : , : , :),2),3),4),5),6)) ...
            ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvNonVax(: , 2 , 3 : 7 , : , : , :),2),3),4),5),6))
            (annlz(sum(vaxResult{n}.popVec(: , allFpos) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-positive on ART
        allFart = toInd(allcomb(8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        hpvIncTimeArt(j , :) = ...
            ((annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvVax(: , 2 , 8 , : , : , :),2),3),4),5),6)) + ...
            annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvVax(: , 2 , 8 , : , : , :),2),3),4),5),6)) + ...
            annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvNonVax(: , 2 , 8 , : , : , :),2),3),4),5),6)) ...
            ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvNonVax(: , 2 , 8 , : , : , :),2),3),4),5),6))
            (annlz(sum(vaxResult{n}.popVec(: , allFart) , 2) ./ stepsPerYear)) * fac);
        
        % HIV-positive all
        allFposAll = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        hpvIncTimePosAll(j , :) = ...
            ((annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvVax(: , 2 , 3 : 8 , : , : , :),2),3),4),5),6)) + ...
            annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvVax(: , 2 , 3 : 8 , : , : , :),2),3),4),5),6)) + ...
            annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvNonVax(: , 2 , 3 : 8 , : , : , :),2),3),4),5),6)) ...
            ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvNonVax(: , 2 , 3 : 8 , : , : , :),2),3),4),5),6))
            (annlz(sum(vaxResult{n}.popVec(: , allFposAll) , 2) ./ stepsPerYear)) * fac);
        
        %% HPV incidence by HIV status and age over time
        fac = 100;
        for dInd = 1 : diseaseVecLength_ccInc;
            d = diseaseVec_ccInc{dInd};
            for rInd = 1 : risk
                r = rInd;
                for a = 1 : age
                    allF = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                        1 : endpoints , 1 : intervens , 2 , a , r));
                    % Calculate incidence
                    hpvIncHivRiskAgeTime(j , dInd , rInd , a , :) = ...
                        ((annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvVax(: , 2 , d , a , r , :),2),3),4),5),6)) + ...
                        annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvVax(: , 2 , d , a , r , :),2),3),4),5),6)) + ...
                        annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHpvNonVax(: , 2 , d , a , r , :),2),3),4),5),6)) ...
                        ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{n}.newImmHpvNonVax(: , 2 , d , a , r , :),2),3),4),5),6))
                        (annlz(sum(vaxResult{n}.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);
                end
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
        
        %% CIN2/3 prevalence for All HR HPV types combined over time
        % General
        cinGenInds = unique([toInd(allcomb(1 : disease , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 4 : 13 , 1 : risk))]);
        ageGenInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk));
        cinGenTime(j , :) = sum(vaxResult{n}.popVec(: , cinGenInds) , 2)...
            ./ sum(vaxResult{n}.popVec(: , ageGenInds) , 2);
        % HIV-positive (on and not on ART)
        cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 4 : 13 , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk));
        cinPosTime(j , :) = sum(vaxResult{n}.popVec(: , cinInds) , 2)...
            ./ sum(vaxResult{n}.popVec(: , ageInds) , 2);
        % HIV-negative
        cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 4 : 13 , 1 : risk))]);
        ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk));
        cinNegTime(j , :) = sum(vaxResult{n}.popVec(: , cinNegInds) , 2)...
            ./ sum(vaxResult{n}.popVec(: , ageNegInds) , 2);
        
        %% CIN2/3 prevalence by broad age groups 2019
        for dInd = 1 : diseaseVecLength_fHpv
            d = diseaseVec_fHpv{dInd};                 
            for aV = 1 : ageVecLength_fHPV
                aGroup = ageVec_fHPV{aV};
                cinInds = unique([toInd(allcomb(d , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , aGroup , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                    [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , aGroup , 1 : risk))]);
                ageInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , aGroup , 1 : risk));
                cin_hivTot2020(j , dInd , aV) = (sum(vaxResult{n}.popVec((2019 - startYear) * stepsPerYear +1 , cinInds)))...
                    ./ (sum(vaxResult{n}.popVec((2019 - startYear) * stepsPerYear +1 , ageInds)));
            end
        end
               
        %% CIN2/3 prevalence ratios in 2005, 2018, 2019
        for i = 1 : hpvPrevYearVecLength
            yr = hpvPrevYearVec(i);  
            for dInd = 1 : diseaseVecLength_hpv
                d = diseaseVec_hpv{dInd};
                cinInds = unique([toInd(allcomb(d , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                    [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 4 : 13 , 1 : risk))]);
                popInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk));
                cin_prev_ratios(j , dInd , i) = sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , cinInds) , 2)...
                    ./ sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , popInds) , 2);
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
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        ccIncTime(j , :) = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , : , 1 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-negative
        allFneg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        ccIncTimeNeg(j , :) = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , 1 : 2 , 1 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , allFneg) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-positive untreated
        allFpos = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        ccIncTimePos(j , :) = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , 3 : 7 , 1 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , allFpos) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-positive on ART
        allFart = toInd(allcomb(8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        ccIncTimeArt(j , :) = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , 8 , 1 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , allFart) , 2) ./ stepsPerYear)) * fac);
        
        % HIV-positive all
        allFposAll = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        ccIncTimePosAll(j , :) = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , 3 : 8 , 1 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , allFposAll) , 2) ./ stepsPerYear)) * fac);
        
        %% Cervical cancer incidence by HIV status and age over time
        fac = 10 ^ 5;
        for dInd = 1 : diseaseVecLength_ccInc;
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                allF = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                % Calculate incidence
                ccIncHivAgeTime(j , dInd , a , :) = ...
                    (annlz(sum(sum(vaxResult{n}.newCC(: , d , a , :),2),4)) ./ ...
                    (annlz(sum(vaxResult{n}.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);
            end
        end
        
        %% Cumulative cervical cancer cases by HIV status and age over time
        for dInd = 1 : diseaseVecLength_ccInc;
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                % Calculate incidence
                ccCumHivAgeTime(j , dInd , a , :) = ...
                    cumsum(squeeze(annlz(sum(sum(vaxResult{n}.newCC(((2019 - startYear) * stepsPerYear +1):end , d , a , :),2),4))),2);
            end
        end

        %% Cervical cancer prevalence ratios in 2005, 2018, 2019
        for i = 1 : hpvPrevYearVecLength
            yr = hpvPrevYearVec(i);  
            for dInd = 1 : diseaseVecLength_hpv
                d = diseaseVec_hpv{dInd};
                ccInds = unique([toInd(allcomb(d , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                    1 : hpvVaxStates , 6 , 1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk))]);
                popInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk));
                cc_prev_ratios(j , dInd , i) = sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , ccInds) , 2)...
                    ./ sum(vaxResult{n}.popVec((yr - startYear) * stepsPerYear +1 , popInds) , 2);
            end
        end
        
       
        %% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************
        
        %% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)
        for dInd = 1 : diseaseIndsLength_typeDist
            d = diseaseInds_typeDist{dInd};
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
        
%         %% Screening "coverage"
%         allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : endpoints , 1 : intervens , 2 , 7 : age , 1 : risk));
%         % Calculate incidence
%         newScreenTime(j , :) = ...
%             annlz(sum(sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , : , : , : , : , : , : , : , :),2),3),4),5),6),7),8),9)) ./ ...
%             (annlz(sum(vaxResult{n}.popVec(((2020 - startYear) * stepsPerYear +1):end , allF) , 2) ./ stepsPerYear) * 0.1);
%         
        %% Screening coverage ages 35-39
%         allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : endpoints , 1 : intervens , 2 , 8 , 1 : risk));
%         screenF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : endpoints , [3 , 4] , 2 , 8 , 1 : risk));
%         
%         screenCovTime(j , :) = ...
%             sum(vaxResult{n}.popVec(((2020 - startYear) * stepsPerYear +1):end , screenF) , 2) ./ ...
%             sum(vaxResult{n}.popVec(((2020 - startYear) * stepsPerYear +1):end , allF) , 2);
        
        %% Total number of women screened annually by age and disease status
%         for dInd = 1 : diseaseVecLength_ccInc
%             d = diseaseVec_ccInc{dInd};
%             screenTotAnnual(j , dInd , :) = annlz(sum(sum(sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , : , : , : , : , : , : , :),2),3),4),5),6),7),8),9));
%         end
        
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
        
        %% Total number of women vaccinated by age and disease status
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                vaxInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , [2 , 4] , 2 , a , 1 : risk));
                vaxTotAge(j , a , dInd , :) = sum(vaxResult{n}.popVec(: , vaxInds) , 2);
            end
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

%% Population proportion by broad age groups over time
% Load calibration data from Excel
file = [pwd , '/Config/Population_validation_targets.xlsx'];
kzn_popByageC_yrs(: , 1) = xlsread(file , 'Demographics' , 'Q92:Q107').*1000;    % females by age in 2019
popPropBroadC_obs = zeros(ageVecLength_cPopDist , 1);
for aV = 1 : ageVecLength_cPopDist
    a = ageVec_cPopDist{aV};
    popPropBroadC_obs(aV , 1) = (sum(kzn_popByageC_yrs(a , 1)) / sumall(kzn_popByageC_yrs(3 : 15 , 1))) * 100;
end

ageGroup = {'10-14' , '15 - 24' , '25 - 34' , '35 - 44' , '45 - 54' , '55 - 64' ,...
    '65 - 74'};
figure('DefaultAxesFontSize' , 18);
subplot(2 , 1 , 1)
plot(1 : ageVecLength_cPopDist , popPropBroadC_obs , 'ro')
hold all;
plot(1 : ageVecLength_cPopDist , median(squeeze(popPropBroadC(: , 1 , :)),1) , 'k-' , 'LineWidth' , 1.5);
hold all;
x2 = [[1 : ageVecLength_cPopDist] , fliplr([1 : ageVecLength_cPopDist])];
inBetween = [max(squeeze(popPropBroadC(: , 1 , :)),[],1) , fliplr(min(squeeze(popPropBroadC(: , 1 , :)),[],1))];
h = fill(x2 , inBetween , 'k');
h.FaceAlpha = 0.3;
h.LineStyle = '--';
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
ylabel('Population Proportion (%)');
ylim([0 50]); grid on;
legend('(Statistics SA) Observed KZN, 2019' , 'Model, 2019: 25-sets median' , ...
    'Model, 2019: 25-sets range' , 'Location' , 'northeast');

subplot(2 , 1 , 2)
plot(1 : ageVecLength_cPopDist , median(squeeze(popPropBroadC(: , 2 , :)),1) , 'k-' , 'LineWidth' , 1.5);
hold all;
x2 = [[1 : ageVecLength_cPopDist] , fliplr([1 : ageVecLength_cPopDist])];
inBetween = [max(squeeze(popPropBroadC(: , 2 , :)),[],1) , fliplr(min(squeeze(popPropBroadC(: , 2 , :)),[],1))];
h = fill(x2 , inBetween , 'k');
h.FaceAlpha = 0.3;
h.LineStyle = '--';
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Population Proportion (%)');
ylim([0 50]); grid on;
legend('Model, 2100: 25-sets median' , ...
    'Model, 2100: 25-sets range' , 'Location' , 'northeast');
sgtitle('Age Distribution');

%% Female total population size by 5-year age groups over time
diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'WmnYrs_' , diseaseLabels{dInd} , '_' , fileKey{n} , '.csv'];
    writematrix([[0 ; (1:age)' ; (1:age)' ; (1:age)'] , ...
        [midAnnualTimespan ;
        [squeeze(median(squeeze(popSizeAgeF(: , dInd , : , (3 : stepsPerYear : end))) , 1)) ; ...
        squeeze(min(squeeze(popSizeAgeF(: , dInd , : , (3 : stepsPerYear : end))) , [] , 1)) ; ...
        squeeze(max(squeeze(popSizeAgeF(: , dInd , : , (3 : stepsPerYear : end))) , [] , 1))]]] , fname)   
end

%% Risk distribution of HIV-negative women with CIN over time
figure;
for r = 1 : risk
    plot(monthlyTimespan , mean(squeeze(popRiskDistCinF(:,r,:)),1)' , 'k-' , ...
        monthlyTimespan , min(squeeze(popRiskDistCinF(:,r,:)),[],1)' , 'k--' , ...
        monthlyTimespan , max(squeeze(popRiskDistCinF(:,r,:)),[],1)' , 'k--' , 'LineWidth' , 1.5);
    hold all;
end
title('Risk distribution of HIV-negative females with CIN3 over time')
xlabel('Year'); ylabel('Proportion')
xlim([1980 2120]); ylim([0 1]);
legend('Model, ages 15-79, LR: 25-sets mean' , 'Model, ages 15-79: 25-sets minimum' , 'Model, ages 15-79: 25-sets maximum' , ...
    'Model, ages 15-79, MR: 25-sets mean' , 'Model, ages 15-79: 25-sets minimum' , 'Model, ages 15-79: 25-sets maximum' , ...
    'Model, ages 15-79, hr: 25-sets mean' , 'Model, ages 15-79: 25-sets minimum' , 'Model, ages 15-79: 25-sets maximum');
grid on;

%% Risk distribution of HIV-negative women over time
figure;
for r = 1 : risk
    plot(monthlyTimespan , mean(squeeze(popRiskDistF(:,r,:)),1)' , 'k-' , ...
        monthlyTimespan , min(squeeze(popRiskDistF(:,r,:)),[],1)' , 'k--' , ...
        monthlyTimespan , max(squeeze(popRiskDistF(:,r,:)),[],1)' , 'k--' , 'LineWidth' , 1.5);
    hold all;
end
title('Risk distribution of HIV-negative females over time')
xlabel('Year'); ylabel('Proportion')
xlim([1980 2120]); ylim([0 1]);
legend('Model, ages 20-34, LR: 25-sets mean' , 'Model, ages 20-34: 25-sets minimum' , 'Model, ages 20-34: 25-sets maximum' , ...
    'Model, ages 20-34, MR: 25-sets mean' , 'Model, ages 20-34: 25-sets minimum' , 'Model, ages 20-34: 25-sets maximum' , ...
    'Model, ages 20-34, hr: 25-sets mean' , 'Model, ages 20-34: 25-sets minimum' , 'Model, ages 20-34: 25-sets maximum');
grid on;

%% Risk distribution of HIV-positive women over time
figure;
for r = 1 : risk
    plot(monthlyTimespan , mean(squeeze(popRiskDistHivF(:,r,:)),1)' , 'k-' , ...
        monthlyTimespan , min(squeeze(popRiskDistHivF(:,r,:)),[],1)' , 'k--' , ...
        monthlyTimespan , max(squeeze(popRiskDistHivF(:,r,:)),[],1)' , 'k--' , 'LineWidth' , 1.5);
    hold all;
end
title('Risk distribution of HIV-positive females over time')
xlabel('Year'); ylabel('Proportion')
xlim([1980 2120]); ylim([0 1]);
legend('Model, ages 15-79, LR: 25-sets mean' , 'Model, ages 15-79: 25-sets minimum' , 'Model, ages 15-79: 25-sets maximum' , ...
    'Model, ages 15-79, MR: 25-sets mean' , 'Model, ages 15-79: 25-sets minimum' , 'Model, ages 15-79: 25-sets maximum' , ...
    'Model, ages 15-79, hr: 25-sets mean' , 'Model, ages 15-79: 25-sets minimum' , 'Model, ages 15-79: 25-sets maximum');
grid on;

%% Risk distribution by HIV status and gender over time
figure;
gen = {'Male' , 'Female'};
for g = 1 : gender
    subplot(1 , 2 , g);
    for r = 1 : risk
        plot(monthlyTimespan , mean(squeeze(popRiskDistHivProp(:,r,g,:)),1)' , 'k-' , ...
            monthlyTimespan , min(squeeze(popRiskDistHivProp(:,r,g,:)),[],1)' , 'k--' , ...
            monthlyTimespan , max(squeeze(popRiskDistHivProp(:,r,g,:)),[],1)' , 'k--' , 'LineWidth' , 1.5);
        hold all;
    end
    title(gen{g});
    xlabel('Year'); ylabel('Proportion')
    xlim([1980 2120]); ylim([0 1]); grid on;
    
end
%suptitle('Proportion of risk group HIV-positive over time')
legend('Model, ages 15-79, LR: 25-sets mean' , 'Model, ages 15-79: 25-sets minimum' , 'Model, ages 15-79: 25-sets maximum' , ...
    'Model, ages 15-79, MR: 25-sets mean' , 'Model, ages 15-79: 25-sets minimum' , 'Model, ages 15-79: 25-sets maximum' , ...
    'Model, ages 15-79, hr: 25-sets mean' , 'Model, ages 15-79: 25-sets minimum' , 'Model, ages 15-79: 25-sets maximum');


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
hivData(: , : , 1) = zeros(length(unique(hivPrevM_dObs(: ,1))) , 1);
hivData(: , : , 2) = zeros(length(unique(hivPrevM_dObs(: ,1))) , 1);
hivRaw(:,:,1) = hivPrevM_dObs(: , 4:5);
hivRaw(:,:,2) = hivPrevF_dObs(: , 4:5);
for i = 1 : length(unique(hivPrevM_dObs(: ,1)))
    for g = 1 : gender
        hivData(i,1,g) = (sumall(hivRaw(((i-1)*7+1):(i*7) , 1 , g)) ./ sumall(hivRaw(((i-1)*7+1):(i*7) , 2 , g))) .* 100;
    end
end

figure('DefaultAxesFontSize' , 18);
gen = {'Males' , 'Females'};
genFlipInd = {2 , 1};
for gInd = 1 : gender
    g = genFlipInd{gInd};
    subplot(1,2,g)
    plot(unique(hivPrevM_dObs(: ,1)) , hivData(:,:,g) , 'ro');
    hold on;
    plot(monthlyTimespan , median(squeeze(hivPrev(: , g , :)),1)' , 'k-' , 'LineWidth' , 1.5);
    x2 = [monthlyTimespan , fliplr(monthlyTimespan)];
    inBetween = [max(squeeze(hivPrev(: , g , :)),[],1) , fliplr(min(squeeze(hivPrev(: , g , :)),[],1))];
    h = fill(x2 , inBetween , 'k');
    h.FaceAlpha = 0.3;
    h.LineStyle = '--';
    xlabel('Year'); ylabel('HIV Prevalence (%)'); title(gen{g});
    xlim([1990 2020]); ylim([0 50]);
    if g == 1 
        legend('(AHRI data request) Observed KZN, ages 15-49' , ...
            'Model, ages 15-49: 25-sets median' , 'Model, ages 15-49: 25-sets range' , ...
            'Location' , 'SouthEast');
    elseif g == 2
        legend('(AHRI data request) Observed KZN, ages 15-49' , ...
            'Model, ages 15-49: 25-sets median' , 'Model, ages 15-49: 25-sets range' , ...
            'Location' , 'SouthEast');
    end
end
sgtitle('HIV prevalence over time');

%% HIV Prevalence by age in 2009, 2018 vs. AHRI (calibration)
ageGroup = {'0-4' , '5-9' , '10-14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
% Calibration error bars
hivM_2009(: , 1) = hivPrevM_dObs(end-6:end , 2); % mean
hivM_2009(: , 2) = (hivPrevM_dObs(end-6:end , 3).^(1/2)) .* 2; % calibration SD
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
    errorbar([4 : 10], hivPrevs(: , 1)' , hivPrevs(: , 2)' , 'rs' , 'LineWidth' , 1.5);
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
sgtitle('HIV prevalence by gender and age over time');

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


%% HIV-associated deaths by gender over time
figure;
plot(annualTimespan , mean(hivDeathsM(: , :),1)' , 'k-' , ...
    annualTimespan , min(hivDeathsM(: , :),[],1)' , 'k--' , ...
    annualTimespan , max(hivDeathsM(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold all;
plot(annualTimespan , mean(hivDeathsF(: , :),1)' , 'b-' , ...
    annualTimespan , min(hivDeathsF(: , :),[],1)' , 'b--' , ...
    annualTimespan , max(hivDeathsF(: , :),[],1)' , 'b--' , 'LineWidth' , 1.5);
xlabel('Year'); ylabel('HIV-associated deaths'); title('HIV-associated deaths by gender over time'); grid on;
xlim([1980 2030]);
legend('Model, males aged 0-79: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, females aged 0-79: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum');

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

%% Female HPV Prevalence by broad age groups and HIV status in 2019
diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'HPVprevF2019_' , diseaseLabels{dInd} , '_' , fileKey{n} , '.csv'];
    writematrix([[0 ; (1:8)' ; (1:8)' ; (1:8)' ; (1:8)' ; (1:8)'] , ...
        [2019 ;
        [squeeze(median(squeeze(hpv_hivTot2020(: , dInd , :)) , 1))' ; ...
        squeeze(prctile(squeeze(hpv_hivTot2020(: , dInd , :)) , 10 , 1))' ; ...
        squeeze(prctile(squeeze(hpv_hivTot2020(: , dInd , :)) , 90 , 1))' ; ...
        squeeze(min(squeeze(hpv_hivTot2020(: , dInd , :)) , [] , 1))' ; ...
        squeeze(max(squeeze(hpv_hivTot2020(: , dInd , :)) , [] , 1))']]] , fname)
end

%% Female HPV Prevalence over time by HIV status
figure;
subplot(1,2,2);
plot(monthlyTimespan , mean(hpv_hivTimeF(: , :),1)' , 'k-' , ...
    monthlyTimespan , min(hpv_hivTimeF(: , :),[],1)' , 'k--' , ...
    monthlyTimespan , max(hpv_hivTimeF(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
    xlabel('Year'); ylabel('HPV Prevalence'); title('WLWHIV ages 15-64  (includes CIN)');
    xlim([1985 2120]); ylim([0 1]); grid on;
    legend('Model, ages 15-64: 25-sets mean' , 'Model, ages 15-64: 25-sets minimum' , ...
        'Model, ages 15-64: 25-sets maximum');
    
subplot(1,2,1);
plot(monthlyTimespan , mean(hpv_hivNegTimeF(: , :),1)' , 'k-' , ...
    monthlyTimespan , min(hpv_hivNegTimeF(: , :),[],1)' , 'k--' , ...
    monthlyTimespan , max(hpv_hivNegTimeF(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
xlabel('Year'); ylabel('HPV Prevalence'); title('HIV-negative women ages 15-64 (includes CIN)');
xlim([1985 2120]); ylim([0 1]); grid on;
legend('Model, ages 15-64: 25-sets mean' , 'Model, ages 15-64: 25-sets minimum' , ...
    'Model, ages 15-64: 25-sets maximum');
sgtitle('Female hrHPV Prevalence (includes CIN) by HIV status');

%% HPV Prevalence over time by sex    
figure;
gen = {'Males aged 15-64' , 'Females aged 15-64 (includes CIN)'};
for g = 1 : gender
    subplot(1,2,g)
    plot(monthlyTimespan , mean(squeeze(hpv_time(: , g , :)),1)' , 'k-' , ...
        monthlyTimespan , min(squeeze(hpv_time(: , g , :)),[],1)' , 'k--' , ...
        monthlyTimespan , max(squeeze(hpv_time(: , g , :)),[],1)' , 'k--' , 'LineWidth' , 1.5);
    xlabel('Year'); ylabel('HPV Prevalence'); title(gen{g});
    xlim([1985 2120]); ylim([0 1]); grid on;
    legend('Model, ages 15-64: 25-sets mean' , 'Model, ages 15-64: 25-sets minimum' , 'Model, ages 15-64: 25-sets maximum');
end
sgtitle('hrHPV Prevalence (includes CIN) by gender');

%% HPV prevalence ratios in 2005 and 2018
hpvRatioHivStatusF = squeeze(hpv_prev_ratios(: , 2 , 2 , :) ./ hpv_prev_ratios(: , 2 , 1 , :));
hpvRatioHivStatusM = squeeze(hpv_prev_ratios(: , 1 , 2 , :) ./ hpv_prev_ratios(: , 1 , 1 , :));
hpvRatioArtStatusF = squeeze(hpv_prev_ratios(: , 2 , 4 , :) ./ hpv_prev_ratios(: , 2 , 3 , :));
hpvRatioArtStatusM = squeeze(hpv_prev_ratios(: , 1 , 4 , :) ./ hpv_prev_ratios(: , 1 , 3 , :));

% HPV prevalence women (HIV+/HIV-)
disp(['HPV prevalence women ages 15-64 in 2005: Ratio (HIV+/HIV-) = ' , num2str(median(hpvRatioHivStatusF(: , 1) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioHivStatusF(: , 1) , 1)) , ' ( ' , num2str(min(hpvRatioHivStatusF(: , 1),[],1)) , ' ' , ...
    num2str(max(hpvRatioHivStatusF(: , 1),[],1)) , ' )'])
disp(['HPV prevalence women ages 15-64 in 2018: Ratio (HIV+/HIV-) = ' , num2str(median(hpvRatioHivStatusF(: , 2) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioHivStatusF(: , 2) , 1)) , ' ( ' , num2str(min(hpvRatioHivStatusF(: , 2),[],1)) , ' ' , ...
    num2str(max(hpvRatioHivStatusF(: , 2),[],1)) , ' )'])

% HPV prevalence women (HIV+ on ART/HIV+ no ART)
disp(['HPV prevalence women ages 15-64 in 2005: Ratio (ART/noART) = ' , num2str(median(hpvRatioArtStatusF(: , 1) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioArtStatusF(: , 1) , 1)) , ' ( ' , num2str(min(hpvRatioArtStatusF(: , 1),[],1)) , ' ' , ...
    num2str(max(hpvRatioArtStatusF(: , 1),[],1)) , ' )'])
disp(['HPV prevalence women ages 15-64 in 2018: Ratio (ART/noART) = ' , num2str(median(hpvRatioArtStatusF(: , 2) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioArtStatusF(: , 2) , 1)) , ' ( ' , num2str(min(hpvRatioArtStatusF(: , 2),[],1)) , ' ' , ...
    num2str(max(hpvRatioArtStatusF(: , 2),[],1)) , ' )'])
    
% HPV prevalence men (HIV+/HIV-)
disp(['HPV prevalence men ages 15-64 in 2005: Ratio (HIV+/HIV-) = ' , num2str(median(hpvRatioHivStatusM(: , 1) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioHivStatusM(: , 1) , 1)) , ' ( ' , num2str(min(hpvRatioHivStatusM(: , 1),[],1)) , ' ' , ...
    num2str(max(hpvRatioHivStatusM(: , 1),[],1)) , ' )'])
disp(['HPV prevalence men ages 15-64 in 2018: Ratio (HIV+/HIV-) = ' , num2str(median(hpvRatioHivStatusM(: , 2) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioHivStatusM(: , 2) , 1)) , ' ( ' , num2str(min(hpvRatioHivStatusM(: , 2),[],1)) , ' ' , ...
    num2str(max(hpvRatioHivStatusM(: , 2),[],1)) , ' )'])

% HPV prevalence men (HIV+ on ART/HIV+ no ART)
disp(['HPV prevalence men ages 15-64 in 2005: Ratio (ART/noART) = ' , num2str(median(hpvRatioArtStatusM(: , 1) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioArtStatusM(: , 1) , 1)) , ' ( ' , num2str(min(hpvRatioArtStatusM(: , 1),[],1)) , ' ' , ...
    num2str(max(hpvRatioArtStatusM(: , 1),[],1)) , ' )'])
disp(['HPV prevalence men ages 15-64 in 2018: Ratio (ART/noART) = ' , num2str(median(hpvRatioArtStatusM(: , 2) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioArtStatusM(: , 2) , 1)) , ' ( ' , num2str(min(hpvRatioArtStatusM(: , 2),[],1)) , ' ' , ...
    num2str(max(hpvRatioArtStatusM(: , 2),[],1)) , ' )'])

%% HPV prevalence ratios in 2019
hpvRatioHivStatusF = squeeze(hpv_prev_ratios(: , 2 , 2 , :) ./ hpv_prev_ratios(: , 2 , 1 , :));
hpvRatioHivStatusM = squeeze(hpv_prev_ratios(: , 1 , 2 , :) ./ hpv_prev_ratios(: , 1 , 1 , :));
hpvRatioArtStatusF = squeeze(hpv_prev_ratios(: , 2 , 4 , :) ./ hpv_prev_ratios(: , 2 , 3 , :));
hpvRatioArtStatusM = squeeze(hpv_prev_ratios(: , 1 , 4 , :) ./ hpv_prev_ratios(: , 1 , 3 , :));

% HPV prevalence women (HIV+/HIV-)
disp(['HPV prevalence women ages 15-64 in 2019: Ratio (HIV+/HIV-) = ' , num2str(median(hpvRatioHivStatusF(: , 3) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioHivStatusF(: , 3) , 1)) , ...
    ' ( ' , num2str(prctile(hpvRatioHivStatusF(: , 3),10,1)) , ' ' , num2str(prctile(hpvRatioHivStatusF(: , 3),90,1)) , ' ) ' , ...
    ' ( ' , num2str(min(hpvRatioHivStatusF(: , 3),[],1)) , ' ' , num2str(max(hpvRatioHivStatusF(: , 3),[],1)) , ' )'])

% HPV prevalence women (HIV+ on ART/HIV+ no ART)
disp(['HPV prevalence women ages 15-64 in 2019: Ratio (ART/noART) = ' , num2str(median(hpvRatioArtStatusF(: , 3) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioArtStatusF(: , 3) , 1)) , ...
    ' ( ' , num2str(prctile(hpvRatioArtStatusF(: , 3),10,1)) , ' ' , num2str(prctile(hpvRatioArtStatusF(: , 3),90,1)) , ' ) ' , ...
    ' ( ' , num2str(min(hpvRatioArtStatusF(: , 3),[],1)) , ' ' , num2str(max(hpvRatioArtStatusF(: , 3),[],1)) , ' )'])
    
% HPV prevalence men (HIV+/HIV-)
disp(['HPV prevalence men ages 15-64 in 2019: Ratio (HIV+/HIV-) = ' , num2str(median(hpvRatioHivStatusM(: , 3) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioHivStatusM(: , 3) , 1)) , ...
    ' ( ' , num2str(prctile(hpvRatioHivStatusM(: , 3),10,1)) , ' ' , num2str(prctile(hpvRatioHivStatusM(: , 3),90,1)) , ' ) ' , ...
    ' ( ' , num2str(min(hpvRatioHivStatusM(: , 3),[],1)) , ' ' , num2str(max(hpvRatioHivStatusM(: , 3),[],1)) , ' )'])

% HPV prevalence men (HIV+ on ART/HIV+ no ART)
disp(['HPV prevalence men ages 15-64 in 2019: Ratio (ART/noART) = ' , num2str(median(hpvRatioArtStatusM(: , 3) , 1)) , ' / ' , ...
    num2str(mean(hpvRatioArtStatusM(: , 3) , 1)) , ...
    ' ( ' , num2str(prctile(hpvRatioArtStatusM(: , 3),10,1)) , ' ' , num2str(prctile(hpvRatioArtStatusM(: , 3),90,1)) , ' ) ' , ...
    ' ( ' , num2str(min(hpvRatioArtStatusM(: , 3),[],1)) , ' ' , num2str(max(hpvRatioArtStatusM(: , 3),[],1)) , ' )'])

%% HPV incidence over time
figure;   
% General
plot(annualTimespan , mean(hpvIncTime,1) , 'k-' , ...
    annualTimespan , min(hpvIncTime,[],1) , 'k--' , ...
    annualTimespan , max(hpvIncTime,[],1) , 'k--' , 'LineWidth' , 1.5);
hold all;
% HIV-negative
plot(annualTimespan , mean(hpvIncTimeNeg,1) , 'g-' , ...
    annualTimespan , min(hpvIncTimeNeg,[],1) , 'g--' , ...
    annualTimespan , max(hpvIncTimeNeg,[],1) , 'g--' , 'LineWidth' , 1.5);
hold all;
% HIV-positive untreated
plot(annualTimespan , mean(hpvIncTimePos,1) , 'r-' , ...
    annualTimespan , min(hpvIncTimePos,[],1) , 'r--' , ...
    annualTimespan , max(hpvIncTimePos,[],1) , 'r--' , 'LineWidth' , 1.5);
hold all;
% HIV-positive on ART
plot(annualTimespan , mean(hpvIncTimeArt,1) , 'b-' , ...
    annualTimespan , min(hpvIncTimeArt,[],1) , 'b--' , ...
    annualTimespan , max(hpvIncTimeArt,[],1) , 'b--' , 'LineWidth' , 1.5);
xlabel('Time'); ylabel('HPV incidence per 100');
xlim([1980 2120]); ylim([0 100]); grid on;
title(['HPV incidence over time, ages 1-79']);
legend('Model, general: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, HIV-negative: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, WLWHIV untreated: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Model, WLWHIV on ART: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
    'Location' , 'northwest');

%% HPV incidence by HIV status and age over time
riskLabels = {'lr' , 'mr' , 'hr'};
diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};
for r = 1 : risk    
    for dInd = 1 : length(diseaseLabels);
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
            'HPVinc_' , diseaseLabels{dInd} , '_' , riskLabels{r} , '_' , fileKey{n} , '.csv'];
        writematrix([[0 ; (1:age)' ; (1:age)' ; (1:age)'] , ...
            [annualTimespan ;
            [squeeze(median(squeeze(hpvIncHivRiskAgeTime(: , dInd , r , : , :)) , 1)) ; ...
            squeeze(min(squeeze(hpvIncHivRiskAgeTime(: , dInd , r , : , :)) , [] , 1)) ; ...
            squeeze(max(squeeze(hpvIncHivRiskAgeTime(: , dInd , r , : , :)) , [] , 1))]]] , fname)
    end
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

%% CIN2/3 prevalence for All HR HPV types combined by HIV status over time
figure;
subplot(1,2,2);
plot(monthlyTimespan , mean(cinPosTime(: , :),1)' , 'k-' , ...
    monthlyTimespan , min(cinPosTime(: , :),[],1)' , 'k--' , ...
    monthlyTimespan , max(cinPosTime(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
    xlabel('Year'); ylabel('CIN 2/3 Prevalence'); title('WLWHIV ages 15-64');
    xlim([1985 2030]); ylim([0 0.3]); grid on;
    legend('Model, ages 15-64: 25-sets mean' , 'Model, ages 15-64: 25-sets minimum' , ...
        'Model, ages 15-64: 25-sets maximum');
    
subplot(1,2,1);
plot(monthlyTimespan , mean(cinNegTime(: , :),1)' , 'k-' , ...
    monthlyTimespan , min(cinNegTime(: , :),[],1)' , 'k--' , ...
    monthlyTimespan , max(cinNegTime(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
xlabel('Year'); ylabel('CIN 2/3 Prevalence'); title('HIV-negative women ages 15-64');
xlim([1985 2030]); ylim([0 0.3]); grid on;
legend('Model, ages 15-64: 25-sets mean' , 'Model, ages 15-64: 25-sets minimum' , ...
    'Model, ages 15-64: 25-sets maximum');
sgtitle('CIN2/3 Prevalence by HIV status');

%% CIN2/3 prevalence for All HR HPV types combined over time
figure;
plot(monthlyTimespan , mean(cinGenTime(: , :),1)' , 'k-' , ...
    monthlyTimespan , min(cinGenTime(: , :),[],1)' , 'k--' , ...
    monthlyTimespan , max(cinGenTime(: , :),[],1)' , 'k--' , 'LineWidth' , 1.5);
    xlabel('Year'); ylabel('CIN 2/3 Prevalence'); title('CIN2/3 prevalence for women aged 15-64');
    xlim([1985 2030]); ylim([0 0.1]); grid on;
    legend('Model, ages 15-64: 25-sets mean' , 'Model, ages 15-64: 25-sets minimum' , ...
        'Model, ages 15-64: 25-sets maximum');
    
%% CIN2/3 Prevalence by broad age groups and HIV status in 2019
diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'CINprevF2019_' , diseaseLabels{dInd} , '_' , fileKey{n} , '.csv'];
    writematrix([[0 ; (1:8)' ; (1:8)' ; (1:8)' ; (1:8)' ; (1:8)'] , ...
        [2019 ;
        [squeeze(median(squeeze(cin_hivTot2020(: , dInd , :)) , 1))' ; ...
        squeeze(prctile(squeeze(cin_hivTot2020(: , dInd , :)) , 10 , 1))' ; ...
        squeeze(prctile(squeeze(cin_hivTot2020(: , dInd , :)) , 90 , 1))' ; ...
        squeeze(min(squeeze(cin_hivTot2020(: , dInd , :)) , [] , 1))' ; ...
        squeeze(max(squeeze(cin_hivTot2020(: , dInd , :)) , [] , 1))']]] , fname)
end
    
%% CIN2/3 prevalence ratios in 2005 and 2018
cinRatioHivStatusF = squeeze(cin_prev_ratios(: , 2 , :) ./ cin_prev_ratios(: , 1 , :));
cinRatioArtStatusF = squeeze(cin_prev_ratios(: , 4 , :) ./ cin_prev_ratios(: , 3 , :));

% CIN prevalence women (HIV+/HIV-)
disp(['CIN2/3 prevalence women aged 15-64 in 2005: Ratio (HIV+/HIV-) = ' , num2str(median(cinRatioHivStatusF(: , 1) , 1)) , ' / ' , ...
    num2str(mean(cinRatioHivStatusF(: , 1) , 1)) , ' ( ' , num2str(min(cinRatioHivStatusF(: , 1),[],1)) , ' ' , ...
    num2str(max(cinRatioHivStatusF(: , 1),[],1)) , ' )'])
disp(['CIN2/3 prevalence women aged 15-64 in 2018: Ratio (HIV+/HIV-) = ' , num2str(median(cinRatioHivStatusF(: , 2) , 1)) , ' / ' , ...
    num2str(mean(cinRatioHivStatusF(: , 2) , 1)) , ' ( ' , num2str(min(cinRatioHivStatusF(: , 2),[],1)) , ' ' , ...
    num2str(max(cinRatioHivStatusF(: , 2),[],1)) , ' )'])  

% CIN prevalence women (HIV+ on ART/HIV+ no ART)
disp(['CIN2/3 prevalence women aged 15-64 in 2005: Ratio (ART/noART) = ' , num2str(median(cinRatioArtStatusF(: , 1) , 1)) , ' / ' , ...
    num2str(mean(cinRatioArtStatusF(: , 1) , 1)) , ' ( ' , num2str(min(cinRatioArtStatusF(: , 1),[],1)) , ' ' , ...
    num2str(max(cinRatioArtStatusF(: , 1),[],1)) , ' )'])
disp(['CIN2/3 prevalence women aged 15-64 in 2018: Ratio (ART/noART) = ' , num2str(median(cinRatioArtStatusF(: , 2) , 1)) , ' / ' , ...
    num2str(mean(cinRatioArtStatusF(: , 2) , 1)) , ' ( ' , num2str(min(cinRatioArtStatusF(: , 2),[],1)) , ' ' , ...
    num2str(max(cinRatioArtStatusF(: , 2),[],1)) , ' )'])  

%% CIN2/3 prevalence ratios in 2019
cinRatioHivStatusF = squeeze(cin_prev_ratios(: , 2 , :) ./ cin_prev_ratios(: , 1 , :));
cinRatioArtStatusF = squeeze(cin_prev_ratios(: , 4 , :) ./ cin_prev_ratios(: , 3 , :));

% CIN prevalence women (HIV+/HIV-)
disp(['CIN2/3 prevalence women aged 15-64 in 2019: Ratio (HIV+/HIV-) = ' , num2str(median(cinRatioHivStatusF(: , 3) , 1)) , ' / ' , ...
    num2str(mean(cinRatioHivStatusF(: , 3) , 1)) , ...
    ' ( ' , num2str(prctile(cinRatioHivStatusF(: , 3),10,1)) , ' ' , num2str(prctile(cinRatioHivStatusF(: , 3),90,1)) , ' ) ' , ...
    ' ( ' , num2str(min(cinRatioHivStatusF(: , 3),[],1)) , ' ' , num2str(max(cinRatioHivStatusF(: , 3),[],1)) , ' )'])

% CIN prevalence women (HIV+ on ART/HIV+ no ART)
disp(['CIN2/3 prevalence women aged 15-64 in 2019: Ratio (ART/noART) = ' , num2str(median(cinRatioArtStatusF(: , 3) , 1)) , ' / ' , ...
    num2str(mean(cinRatioArtStatusF(: , 3) , 1)) , ...
    ' ( ' , num2str(prctile(cinRatioArtStatusF(: , 3),10,1)) , ' ' , num2str(prctile(cinRatioArtStatusF(: , 3),90,1)) , ' ) ' , ...
    ' ( ' , num2str(min(cinRatioArtStatusF(: , 3),[],1)) , ' ' , num2str(max(cinRatioArtStatusF(: , 3),[],1)) , ' )'])


%% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

%% Cervical cancer incidence in 2005 by age
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

figure;    
plot(1 : age , mean(squeeze(ccIncAge(: , 1 , :)),1) , 'k-' , ...
    1 : age , min(squeeze(ccIncAge(: , 1 , :)),[],1) , 'k--' , ...
    1 : age , max(squeeze(ccIncAge(: , 1 , :)),[],1) , 'k--' , 'LineWidth' , 1.5);
xlabel('Age Group'); ylabel('Cervical cancer incidence per 100K');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
ylim([0 175]); grid on;
title(['Cervical Cancer Incidence in 2005 by age']);
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
ylim([0 175]); grid on;
title(['Cervical Cancer Incidence in 2012 by age']);
legend('(Globocan, 2012) Observed SA: mean, 2SD' , 'Model, general: 25-sets mean' , ...
    'Model: 25-sets minimum' , 'Model: 25-sets maximum');
   
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

%% Cervical cancer incidence by HIV status and age over time
diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};
for dInd = 1 : length(diseaseLabels);
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'ICC_' , diseaseLabels{dInd} , '_' , fileKey{n} , '.csv'];
    writematrix([[0 ; (1:age)' ; (1:age)' ; (1:age)'] , ...
        [annualTimespan ;
        [squeeze(median(squeeze(ccIncHivAgeTime(: , dInd , : , :)) , 1)) ; ...
        squeeze(min(squeeze(ccIncHivAgeTime(: , dInd , : , :)) , [] , 1)) ; ...
        squeeze(max(squeeze(ccIncHivAgeTime(: , dInd , : , :)) , [] , 1))]]] , fname)
end

%% CC incidence - age standardized
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

diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};
for dInd = 1 : length(diseaseLabels);
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'ICC-medAS_' , diseaseLabels{dInd} , '_' , fileKey{n} , '.csv'];
    ccIncHivAgeTime_med = squeeze(median(squeeze(ccIncHivAgeTime(: , dInd , : , :)) , 1));

    ccIncRefTot = zeros(1 , size(ccIncHivAgeTime_med,2));       
    for aInd = 1:age+4
        a = aInd;
        if aInd >= age
            a = age;
        end

        if aInd <= age    
            ccIncRef = ccIncHivAgeTime_med(a , :) .* worldStandard_WP2015(aInd);
            if (a < 3)
                ccIncRef = zeros(1 , size(ccIncHivAgeTime_med,2));
            end
        elseif aInd > age
            ccIncRef = ccIncHivAgeTime_med(a , :);
            ccIncRef = [(ones(1,aInd-a).*ccIncRef(1,1)) , ccIncRef(1,1:end-(aInd-a))];
            ccIncRef = ccIncRef .* worldStandard_WP2015(aInd);
        end
        ccIncRefTot = ccIncRefTot + ccIncRef;
    end
    ccInc = ccIncRefTot ./ (sum(worldStandard_WP2015(1:age+4)));

    writematrix([annualTimespan ; ccInc] , fname)
end

%% Cumulative cervical cancer cases by HIV status and age over time
diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};
for dInd = 1 : length(diseaseLabels);
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'CumCC_' , diseaseLabels{dInd} , '_' , fileKey{n} , '.csv'];
    writematrix([[0 ; (1:age)' ; (1:age)' ; (1:age)'] , ...
        [futAnnualTimespan ;
        [squeeze(median(squeeze(ccCumHivAgeTime(: , dInd , : , :)) , 1)) ; ...
        squeeze(min(squeeze(ccCumHivAgeTime(: , dInd , : , :)) , [] , 1)) ; ...
        squeeze(max(squeeze(ccCumHivAgeTime(: , dInd , : , :)) , [] , 1))]]] , fname)
end

%% Cervical cancer prevalence ratios in 2005 and 2018
ccRatioHivStatusF = squeeze(cc_prev_ratios(: , 2 , :) ./ cc_prev_ratios(: , 1 , :));
ccRatioArtStatusF = squeeze(cc_prev_ratios(: , 4 , :) ./ cc_prev_ratios(: , 3 , :));

% CC prevalence women (HIV+/HIV-)
disp(['Cervical cancer prevalence women ages 15-64 in 2005: Ratio (HIV+/HIV-) = ' , num2str(median(ccRatioHivStatusF(: , 1) , 1)) , ' / ' , ...
    num2str(mean(ccRatioHivStatusF(: , 1) , 1)) , ' ( ' , num2str(min(ccRatioHivStatusF(: , 1),[],1)) , ' ' , ...
    num2str(max(ccRatioHivStatusF(: , 1),[],1)) , ' )'])
disp(['Cervical cancer prevalence women 15-64 in 2018: Ratio (HIV+/HIV-) = ' , num2str(median(ccRatioHivStatusF(: , 2) , 1)) , ' / ' , ...
    num2str(mean(ccRatioHivStatusF(: , 2) , 1)) , ' ( ' , num2str(min(ccRatioHivStatusF(: , 2),[],1)) , ' ' , ...
    num2str(max(ccRatioHivStatusF(: , 2),[],1)) , ' )'])  

% CC prevalence women (HIV+ on ART/HIV+ no ART)
disp(['Cervical cancer prevalence women 15-64 in 2005: Ratio (ART/noART) = ' , num2str(median(ccRatioArtStatusF(: , 1) , 1)) , ' / ' , ...
    num2str(mean(ccRatioArtStatusF(: , 1) , 1)) , ' ( ' , num2str(min(ccRatioArtStatusF(: , 1),[],1)) , ' ' , ...
    num2str(max(ccRatioArtStatusF(: , 1),[],1)) , ' )'])
disp(['Cervical cancer prevalence women 15-64 in 2018: Ratio (ART/noART) = ' , num2str(median(ccRatioArtStatusF(: , 2) , 1)) , ' / ' , ...
    num2str(mean(ccRatioArtStatusF(: , 2) , 1)) , ' ( ' , num2str(min(ccRatioArtStatusF(: , 2),[],1)) , ' ' , ...
    num2str(max(ccRatioArtStatusF(: , 2),[],1)) , ' )'])  

%% Cervical cancer incidence ratios in 2005 and 2018
ccRatioHivStatusF_2005 = ccIncTimePosAll(: , ((2005 - startYear) +1)) ./ ccIncTimeNeg(: , ((2005 - startYear) +1));
ccRatioArtStatusF_2005 = ccIncTimeArt(: , ((2005 - startYear) +1)) ./ ccIncTimePos(: , ((2005 - startYear) +1));
ccRatioHivStatusF_2018 = ccIncTimePosAll(: , ((2018 - startYear) +1)) ./ ccIncTimeNeg(: , ((2018 - startYear) +1));
ccRatioArtStatusF_2018 = ccIncTimeArt(: , ((2018 - startYear) +1)) ./ ccIncTimePos(: , ((2018 - startYear) +1));

% CC prevalence women (HIV+/HIV-)
disp(['Cervical cancer incidence women ages 15-79 in 2005: Ratio (HIV+/HIV-) = ' , num2str(median(ccRatioHivStatusF_2005(: , 1) , 1)) , ' / ' , ...
    num2str(mean(ccRatioHivStatusF_2005(: , 1) , 1)) , ' ( ' , num2str(min(ccRatioHivStatusF_2005(: , 1),[],1)) , ' ' , ...
    num2str(max(ccRatioHivStatusF_2005(: , 1),[],1)) , ' )'])
disp(['Cervical cancer incidence women ages 15-79 in 2018: Ratio (HIV+/HIV-) = ' , num2str(median(ccRatioHivStatusF_2018(: , 1) , 1)) , ' / ' , ...
    num2str(mean(ccRatioHivStatusF_2018(: , 1) , 1)) , ' ( ' , num2str(min(ccRatioHivStatusF_2018(: , 1),[],1)) , ' ' , ...
    num2str(max(ccRatioHivStatusF_2018(: , 1),[],1)) , ' )'])  

% CC prevalence women (HIV+ on ART/HIV+ no ART)
disp(['Cervical cancer incidence women ages 15-79 in 2005: Ratio (ART/noART) = ' , num2str(median(ccRatioArtStatusF_2005(: , 1) , 1)) , ' / ' , ...
    num2str(mean(ccRatioArtStatusF_2005(: , 1) , 1)) , ' ( ' , num2str(min(ccRatioArtStatusF_2005(: , 1),[],1)) , ' ' , ...
    num2str(max(ccRatioArtStatusF_2005(: , 1),[],1)) , ' )'])
disp(['Cervical cancer incidence women ages 15-79 in 2018: Ratio (ART/noART) = ' , num2str(median(ccRatioArtStatusF_2018(: , 1) , 1)) , ' / ' , ...
    num2str(mean(ccRatioArtStatusF_2018(: , 1) , 1)) , ' ( ' , num2str(min(ccRatioArtStatusF_2018(: , 1),[],1)) , ' ' , ...
    num2str(max(ccRatioArtStatusF_2018(: , 1),[],1)) , ' )'])  

%% Cervical cancer prevalence ratios in 2019
ccRatioHivStatusF = squeeze(cc_prev_ratios(: , 2 , :) ./ cc_prev_ratios(: , 1 , :));
ccRatioArtStatusF = squeeze(cc_prev_ratios(: , 4 , :) ./ cc_prev_ratios(: , 3 , :));

% CC prevalence women (HIV+/HIV-)
disp(['Cervical cancer prevalence women ages 15-64 in 2019: Ratio (HIV+/HIV-) = ' , num2str(median(ccRatioHivStatusF(: , 3) , 1)) , ' / ' , ...
    num2str(mean(ccRatioHivStatusF(: , 3) , 1)) , ...
    ' ( ' , num2str(prctile(ccRatioHivStatusF(: , 3),10,1)) , ' ' , num2str(prctile(ccRatioHivStatusF(: , 3),90,1)) , ' ) ' , ...
    ' ( ' , num2str(min(ccRatioHivStatusF(: , 3),[],1)) , ' ' , num2str(max(ccRatioHivStatusF(: , 3),[],1)) , ' )'])

% CC prevalence women (HIV+ on ART/HIV+ no ART)
disp(['Cervical cancer prevalence women 15-64 in 2019: Ratio (ART/noART) = ' , num2str(median(ccRatioArtStatusF(: , 3) , 1)) , ' / ' , ...
    num2str(mean(ccRatioArtStatusF(: , 3) , 1)) , ...
    ' ( ' , num2str(prctile(ccRatioArtStatusF(: , 3),10,1)) , ' ' , num2str(prctile(ccRatioArtStatusF(: , 3),90,1)) , ' ) ' , ...
    ' ( ' , num2str(min(ccRatioArtStatusF(: , 3),[],1)) , ' ' , num2str(max(ccRatioArtStatusF(: , 3),[],1)) , ' )'])

%% Cervical cancer incidence ratios in 2019
% ccRatioHivStatusF_2020 = ccIncTimePosAll(: , ((2019 - startYear) +1)) ./ ccIncTimeNeg(: , ((2019 - startYear) +1));
% ccRatioArtStatusF_2020 = ccIncTimeArt(: , ((2019 - startYear) +1)) ./ ccIncTimePos(: , ((2019 - startYear) +1));
% 
% % CC prevalence women (HIV+/HIV-)
% disp(['Cervical cancer incidence women ages 15-74 in 2019: Ratio (HIV+/HIV-) = ' , num2str(median(ccRatioHivStatusF_2020(: , 1) , 1)) , ' / ' , ...
%     num2str(mean(ccRatioHivStatusF_2020(: , 1) , 1)) , ...
%     ' ( ' , num2str(prctile(ccRatioHivStatusF_2020(: , 1),10,1)) , ' ' , num2str(prctile(ccRatioHivStatusF_2020(: , 1),90,1)) , ' ) ' , ...
%     ' ( ' , num2str(min(ccRatioHivStatusF_2020(: , 1),[],1)) , ' ' , num2str(max(ccRatioHivStatusF_2020(: , 1),[],1)) , ' )'])
% 
% % CC prevalence women (HIV+ on ART/HIV+ no ART)
% disp(['Cervical cancer incidence women ages 15-74 in 2019: Ratio (ART/noART) = ' , num2str(median(ccRatioArtStatusF_2020(: , 1) , 1)) , ' / ' , ...
%     num2str(mean(ccRatioArtStatusF_2020(: , 1) , 1)) , ...
%     ' ( ' , num2str(prctile(ccRatioArtStatusF_2020(: , 1),10,1)) , ' ' , num2str(prctile(ccRatioArtStatusF_2020(: , 1),90,1)) , ' ) ' , ...
%     ' ( ' , num2str(min(ccRatioArtStatusF_2020(: , 1),[],1)) , ' ' , num2str(max(ccRatioArtStatusF_2020(: , 1),[],1)) , ' )'])


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
title('Type Distribution by stage (coinfections grouped as 9v), ages 0-79');
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
title('Cervical cancer type distribution by HIV status (coinfections grouped as 9v), ages 0-79');
ylim([0 1]); grid on;
legend('Model- 9v, 2000: 25-sets mean' , 'Model- 9v, 2000: 25-sets minimum' , 'Model- 9v, 2000: 25-sets maximum' , ...
    'Model- 9v, 2018: 25-sets mean' , 'Model- 9v, 2018: 25-sets minimum' , 'Model- 9v, 2018: 25-sets maximum' , ...
    'Model- non-9v, 2000: 25-sets mean' , 'Model- non-9v, 2000: 25-sets minimum' , 'Model- non-9v, 2000: 25-sets maximum' , ...
    'Model- non-9v, 2018: 25-sets mean' , 'Model- non-9v, 2018: 25-sets minimum' , 'Model- non-9v, 2018: 25-sets maximum');

%% HPV type distribution in 2019 by state and HIV status(coinfections grouped as 9v-type HPV)
% diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll'};
% diseaseLabelsInds = {1 , 2 , 5};
% for dInd = 1 : length(diseaseLabels)
%     d = diseaseLabelsInds{dInd};
%     fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
%         'HPVtypeDist2019_' , diseaseLabels{dInd} , '_' , fileKey{n} , '.csv'];
%     writematrix([[0 ; (1:8)' ; (1:8)' ; (1:8)' ; (1:8)' ; (1:8)'] , ...
%         [2019 ;
%         [squeeze(median(squeeze(hpv_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 1)) ; ...
%         squeeze(median(squeeze(hpv_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 1)) ; ...
%         squeeze(median(squeeze(cin1_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 1)) ; ...
%         squeeze(median(squeeze(cin1_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 1)) ; ...
%         squeeze(median(squeeze(cin23_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 1)) ; ...
%         squeeze(median(squeeze(cin23_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 1)) ; ...
%         squeeze(median(squeeze(cc_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 1)) ; ...
%         squeeze(median(squeeze(cc_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 1)) ; ...
%         
%         squeeze(prctile(squeeze(hpv_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 10 , 1)) ; ...
%         squeeze(prctile(squeeze(hpv_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 10 , 1)) ; ...
%         squeeze(prctile(squeeze(cin1_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 10 , 1)) ; ...
%         squeeze(prctile(squeeze(cin1_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 10 , 1)) ; ...
%         squeeze(prctile(squeeze(cin23_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 10 , 1)) ; ...
%         squeeze(prctile(squeeze(cin23_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 10 , 1)) ; ...
%         squeeze(prctile(squeeze(cc_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 10 , 1)) ; ...
%         squeeze(prctile(squeeze(cc_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 10 , 1)) ; ...
%         
%         squeeze(prctile(squeeze(hpv_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 90 , 1)) ; ...
%         squeeze(prctile(squeeze(hpv_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 90 , 1)) ; ...
%         squeeze(prctile(squeeze(cin1_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 90 , 1)) ; ...
%         squeeze(prctile(squeeze(cin1_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 90 , 1)) ; ...
%         squeeze(prctile(squeeze(cin23_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 90 , 1)) ; ...
%         squeeze(prctile(squeeze(cin23_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 90 , 1)) ; ...
%         squeeze(prctile(squeeze(cc_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , 90 , 1)) ; ...
%         squeeze(prctile(squeeze(cc_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , 90 , 1)) ; ...
%         
%         squeeze(min(squeeze(hpv_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(min(squeeze(hpv_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(min(squeeze(cin1_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(min(squeeze(cin1_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(min(squeeze(cin23_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(min(squeeze(cin23_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(min(squeeze(cc_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(min(squeeze(cc_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         
%         squeeze(max(squeeze(hpv_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(max(squeeze(hpv_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(max(squeeze(cin1_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(max(squeeze(cin1_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(max(squeeze(cin23_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(max(squeeze(cin23_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(max(squeeze(cc_vax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1)) ; ...
%         squeeze(max(squeeze(cc_nonVax(:,d,((2019 - startYear) * stepsPerYear +1))) , [] , 1))]]] , fname)
% end


%% ************************** SCREENING & VACCINATION FIGURES *******************************************************************************

%% Screening "coverage"
% figure;   
% plot(screenAnnualTimespan , mean(newScreenTime,1) , 'k-' , ...
%     screenAnnualTimespan , min(newScreenTime,[],1) , 'k--' , ...
%     screenAnnualTimespan , max(newScreenTime,[],1) , 'k--' , 'LineWidth' , 1.5);
% xlabel('Time'); ylabel('Screening coverage');
% xlim([2020 2100]); ylim([0 0.3]); grid on;
% title(['Screening coverage']);
% legend('Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
%     'Location' , 'northwest');

%% Screening coverage ages 35-39
% figure;   
% plot(screenMonthlyTimespan , mean(screenCovTime,1) , 'k-' , ...
%     screenMonthlyTimespan , min(screenCovTime,[],1) , 'k--' , ...
%     screenMonthlyTimespan , max(screenCovTime,[],1) , 'k--' , 'LineWidth' , 1.5);
% xlabel('Time'); ylabel('Screening coverage');
% xlim([2020 2100]); ylim([0 0.3]); grid on;
% title(['Screening coverage ages 35-39']);
% legend('Model: 25-sets mean' , 'Model: 25-sets minimum' , 'Model: 25-sets maximum' , ...
%     'Location' , 'northwest');

%% Total number of women screened annually by age and disease status
% diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};
% for dInd = 1 : length(diseaseLabels)
%     fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
%         'ScreenTot_' , diseaseLabels{dInd} , '_' , fileKey{n} , '.csv'];
%     writematrix([[0 ; 35 ; 35 ; 35] , ...
%         [[(1925+(3/stepsPerYear)) : ((lastYear-1)+(3/stepsPerYear))] ;
%         [squeeze(median(squeeze(screenTotAnnual(: , dInd , :)) , 1)) ; ...
%         squeeze(min(squeeze(screenTotAnnual(: , dInd , :)) , [] , 1)) ; ...
%         squeeze(max(squeeze(screenTotAnnual(: , dInd , :)) , [] , 1))]]] , fname)
% end

%% Vaccine coverage overall
figure;   
plot(midAnnualTimespan , mean(vaxCoverage(: , (3 : stepsPerYear : end)),1) , 'k-' , ...
    midAnnualTimespan , min(vaxCoverage(: , (3 : stepsPerYear : end)),[],1) , 'k--' , ...
    midAnnualTimespan , max(vaxCoverage(: , (3 : stepsPerYear : end)),[],1) , 'k--' , 'LineWidth' , 1.5);
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
    plot(midAnnualTimespan , mean(squeeze(vaxCoverageAge(: , a , (3 : stepsPerYear : end))),1)' , '-' , 'LineWidth' , 1.5);
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

%% Vaccine coverage by age and disease status
diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'VaxTot_' , diseaseLabels{dInd} , '_' , fileKey{n} , '.csv'];
    writematrix([[0 ; (1:age)' ; (1:age)' ; (1:age)'] , ...
        [midAnnualTimespan ;
        [squeeze(median(squeeze(vaxTotAge(: , : , dInd , (3 : stepsPerYear : end))) , 1)) ; ...
        squeeze(min(squeeze(vaxTotAge(: , : , dInd , (3 : stepsPerYear : end))) , [] , 1)) ; ...
        squeeze(max(squeeze(vaxTotAge(: , : , dInd , (3 : stepsPerYear : end))) , [] , 1))]]] , fname)
end 
