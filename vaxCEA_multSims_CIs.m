function [] = vaxCEA_multSims_CIs(vaxResultInd , sceNum)

%% Load parameters and results
paramDir = [pwd , '\Params\'];

[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd, annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , partnersMmult, maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4,...
    mue , mue2 , mue3 , mue4 , mue5, epsA_vec , epsR_vec , yr , ...
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
    prepStartYear , prepCoverage , prepProtect , ...
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
    deathMat , deathMat2 , deathMat3 , deathMat4 , deathMat5,...
    dDeathMat , dDeathMat2 , dDeathMat3 , dDeathMat4, dMue] = loadUp2(1 , 0 , [] , [] , []);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)

lastYear = 2071;

% Indices of calib runs to plot
fileKey = {'noVax' , 'SAC', 'Routine' , 'CU50', 'CU80', 'FASTER70', 'WHO'};
fileInds = num2cell(1:50);
nRuns = length(fileInds);
% 
% Initialize model output plots
popYearVec = [1990:2020 2070];
popYearVecLength = length(popYearVec);
adultsAgeVec = {4, 5, 6, 7, 8, 9,10, 11, 12, 13};
adultsAgeVec_length = length(adultsAgeVec);
ageVec_cPopDist = {[1:2], [3:4] , [5:6] , [7:8] , [9:10] , [11:12], [13:16]};
ageVecLength_cPopDist = length(ageVec_cPopDist);
hpvYearVec = [2005 2009 2020 2070];
hpvYearVec2018_orig = [2018];
hpvYearVecLength = length(hpvYearVec);
hpvYearVec2018Length = length(hpvYearVec2018_orig);
diseaseVec_fHpv = {[1 : disease] , [1 : 2] , [3 : 8] , [3 : 7] , 8};
diseaseVecLength_fHpv = length(diseaseVec_fHpv);
ageVec_fHPV = {[5],[6],[7:8],[9:10]};
ageVecLength_fHPV = length(ageVec_fHPV);
ageVec_fHPV2 = {[3],[4:5],[6:7],[8:9], [10:11], [12:13], [14:15], 16};
ageVecLength_fHPV2 = length(ageVec_fHPV2);
hpvYearVecMale = [2005 2009 2020];
hpvYearVecMaleLength = length(hpvYearVecMale);
ageVec_mHPV = {[4:5],[6:7],[8:9],[10:13]};
ageVecLength_mHPV = length(ageVec_mHPV);
ccYearVec = [2009 2012 2018 2020];
ccYearVecLength = length(ccYearVec);
diseaseVec_ccInc = {[1 : disease] , [1 : 2] , [3 : 8] , [3 : 7] , 8};
diseaseVecLength_ccInc = length(diseaseVec_ccInc);
diseaseVec_hpv = {[1 : 2] , [3 : 8] , [3 : 7] , 8};
hpvPrevYearVec = [2005 2018 2019];
diseaseVecLength_hpv = length(diseaseVec_hpv);
hpvPrevYearVecLength = length(hpvPrevYearVec);
diseaseInds_typeDist = {[1 : disease] , [1 : 2] , [3 : 7] , 8 , [3 : 8]};
diseaseIndsLength_typeDist = length(diseaseInds_typeDist);

% Timespans
monthlyTimespan = [startYear : (1/6) : lastYear];
monthlyTimespan = monthlyTimespan(1 : end-1);
annualTimespan = [startYear : lastYear-1];
futAnnualTimespan = [2020 : lastYear-1];
midAnnualTimespan = [(startYear+(3/stepsPerYear)) : ((lastYear-1)+(3/stepsPerYear))];
screenAnnualTimespan = [(2020+(3/stepsPerYear)) : ((lastYear-1)+(3/stepsPerYear))];
screenMonthlyTimespan = [2020 : (1/6) : lastYear];
screenMonthlyTimespan = screenMonthlyTimespan(1 : end-1);

% Total population size
% popSize = zeros(length(fileKey), nRuns , gender, length(monthlyTimespan));
% popSizeAgeF = zeros(length(fileKey),nRuns , 5 , age , length(monthlyTimespan));
% % Population age distribution
% popPropF = zeros(length(fileKey),nRuns , length(popYearVec) , age);
% popPropBroadC = zeros(length(fileKey),nRuns , length(popYearVec) , ageVecLength_cPopDist);
% % Female risk distribution
% popRiskDistF = zeros(length(fileKey),nRuns , risk , length(monthlyTimespan));
% popRiskDistHivF = popRiskDistF;
% popRiskDistCinF = popRiskDistF;
% % Risk distribution by HIV status and gender
% popRiskDistHivProp = zeros(length(fileKey),nRuns , risk , gender , length(monthlyTimespan));
% 
% % HIV prevalence
% hivAgeM = zeros(length(fileKey),nRuns , age , length(monthlyTimespan));
% hivAgeF = hivAgeM;
% hivPrevGen = zeros(length(fileKey),nRuns , gender , length(monthlyTimespan));
% hivPrev = zeros(length(fileKey), nRuns, length(monthlyTimespan));
% hivPrevRiskF = zeros(length(fileKey),nRuns , risk , length(monthlyTimespan));
% % HIV deaths
% hivDeathsM = zeros(length(fileKey),nRuns , length(annualTimespan));
% hivDeathsF = hivDeathsM;
% hivInc = zeros(length(fileKey), nRuns, gender, length(annualTimespan));
% hivIncAgeM = zeros(length(fileKey), nRuns, age, length(annualTimespan));
% hivIncAgeF = hivIncAgeM;
% hivIncRiskF = zeros(length(fileKey), nRuns, risk, length(annualTimespan));
% hivCumCasesF = zeros(length(fileKey), nRuns, age , risk , length(annualTimespan));
% hivCumCasesM = hivCumCasesF;
% % ART coverage
% artCovM = zeros(length(fileKey),nRuns , length(monthlyTimespan));
% artCovF = artCovM;
% artCovAge = zeros(length(fileKey),nRuns , age , length(monthlyTimespan));

% % % Female HPV prevalence
% hpv_hiv = zeros(length(fileKey),nRuns , age , hpvYearVecLength);
% hpv_hivNeg = hpv_hiv;
% hpv_hivTot = zeros(length(fileKey),nRuns , age , 1);
% hpv_hiv_risk = zeros(length(fileKey),nRuns , length(ageVec_fHPV) );
% hpv_hivNeg_risk = hpv_hiv_risk;
% hpv_hivTot2020 = zeros(length(fileKey),nRuns , diseaseVecLength_fHpv , length(ageVec_fHPV2));
% % Male HPV prevalence
% hpv_hivM = zeros(length(fileKey),nRuns , 4 , 2);
% hpv_hivMNeg = hpv_hivM;
% hpv_hivMtot = hpv_hivM;

% % HPV prevalence over time
% hpv_hivTimeF = zeros(length(fileKey),nRuns , length(monthlyTimespan));
% hpv_hivNegTimeF = hpv_hivTimeF;
% hpv_time = zeros(length(fileKey),nRuns , 2 , length(monthlyTimespan));
% % HPV incidence
% hpvIncTime = zeros(length(fileKey),nRuns , length(annualTimespan));
% hpvIncTimeNeg = hpvIncTime;
% hpvIncTimePos = hpvIncTime;
% hpvIncTimeArt = hpvIncTime;
% hpvIncTimePosAll = hpvIncTime;
% hpvIncHivRiskAgeTime = zeros(length(fileKey),nRuns , 5 , risk , age , length(annualTimespan));
% % CIN2/3 prevalence
% cinPosAge = zeros(length(fileKey),nRuns , 3);
% cinNegAge = cinPosAge;
% % cinGenAge = cinPosAge;
% cin1Age = zeros(length(fileKey),nRuns , length(4:13));
% cin3Age = cin1Age;
% cin_hivTot2020 = zeros(length(fileKey),nRuns , diseaseVecLength_fHpv , ageVecLength_fHPV);
% % CC incidence
% ccIncAge_all = zeros(length(fileKey),nRuns , 3 , age);
% ccIncAge_pos = ccIncAge_all;
% ccIncAge_neg = ccIncAge_all;
% ccIncAge_art = ccIncAge_all;
% ccIncTime = zeros(length(fileKey),nRuns , length(annualTimespan));
% ccIncTimeNeg = ccIncTime;
% ccIncTimePos = ccIncTime;
% ccIncTimeArt = ccIncTime;
% ccIncTimePosAll = ccIncTime;
% ccIncHivAgeTime = zeros(length(fileKey),nRuns , 5 , age , length(annualTimespan));
% ccCumHivAgeTime = zeros(length(fileKey),nRuns , 5 , age , length(annualTimespan));

segi = [0.1212 0.1010 0.0909 0.0909 0.0808 0.0808 0.0606 ...
    0.0606 0.0606 0.0606 0.0505 0.0404 0.0404 0.0303 0.0202 0.0101];
who = [0.0886 0.0869 0.086 0.0847  0.0822 0.0793 0.0761 0.0715 0.0659 ...
    0.0604 0.0537 0.0455 0.0372 0.0296 0.0221 0.0152];
% ccIncStd = zeros(length(fileKey),nRuns , length(diseaseVec_ccInc), length(annualTimespan));
% ccIncStd_who = ccIncStd;
% % Prevalence ratios
% hpv_prev_ratios = zeros(length(fileKey),nRuns , gender , 4 , 2);
% cin_prev_ratios = zeros(length(fileKey),nRuns , 4 , 2);
% cc_prev_ratios = zeros(length(fileKey),nRuns , 4 , 2);
% cc_inc_ratios = zeros(length(fileKey),nRuns , 4 , 2);
% % HPV/CIN/CC type distribution
% 
% cc_vax = zeros(length(fileKey),nRuns , 4 , length(monthlyTimespan));
% cc_nonVax = cc_vax;
% cin23_vax = cc_vax;
% cin23_nonVax = cc_vax;
% cin3_vax = cc_vax;
% cin3_nonVax = cc_vax;
% cin2_vax = cc_vax;
% cin2_nonVax = cc_vax;
% cin1_vax = cc_vax;
% cin1_nonVax = cc_vax;
% hpv_vax = cc_vax;
% hpv_nonVax = cc_vax;
% % HPV vaccination and screening
% % % newScreenTime = zeros(length(fileKey),nRuns , length(screenAnnualTimespan));
% % % screenCovTime = zeros(length(fileKey),nRuns , length(screenMonthlyTimespan));
% % % screenTotAnnual = zeros(length(fileKey),nRuns , 5 , length([(1925+(3/stepsPerYear)) : ((lastYear-1)+(3/stepsPerYear))]));
% vaxCoverage = zeros(length(fileKey),nRuns , length(monthlyTimespan));
% vaxCoverageAge = zeros(length(fileKey),nRuns , age , length(monthlyTimespan));
% vaxTotAge = zeros(length(fileKey),nRuns , age , 5 , length(monthlyTimespan));

resultsDir = ['C:\HPV-HIVacq_model_9v'];

%
% baseFileName = ['14DEC20_stochMod_'];
loopSegments = (0 : round(nRuns/5) : nRuns);
loopSegmentsLength = length(loopSegments);
for z = 6:7
    pathModifier = ['\Vaccine15Jan_stochMod9v_' , fileKey{z}]; % ***SET ME***: name of the folder containing futureSim results
for k = 1 : loopSegmentsLength-1
    parfor j = loopSegments(k)+1 : loopSegments(k+1)
%     for j = 1 : nRuns
        % Load results
        
        nSims = size(dir([resultsDir ,  pathModifier, '\' , '*.mat']) , 1);
        curr = load(['C:\HPV-HIVacq_model_4v\toNow_15Jan_stochMod_' , num2str(j)]); % ***SET ME***: name for historical run output file 

        vaxResult = cell(nSims , 1);
        resultFileName = [resultsDir , pathModifier, '\vaxSimResult']; %, fileInds{j}];
        % load results from vaccine run into cell array
        vaxResult{j} = load([resultFileName , num2str(j), '.mat']);
        % concatenate vectors/matrices of population up to current year to population
        % matrices for years past current year
        vaxResult{j}.popVec = [curr.popVec(1 : end  , :); vaxResult{j}.popVec(2 : end , :)];
        vaxResult{j}.ccDeath = [curr.ccDeath(1 : end , : , : , :) ; vaxResult{j}.ccDeath(2 : end , : , : , :)];
        vaxResult{j}.newCC = [curr.newCC(1 : end , : , : , :); vaxResult{j}.newCC(2 : end , : , : , :)];
        vaxResult{j}.newHpvVax = [curr.newHpvVax(1 : end , : , : , : , : , :); vaxResult{j}.newHpvVax(2 : end , : , : , : , : , :)];
        vaxResult{j}.newImmHpvVax = [curr.newImmHpvVax(1 : end , : , : , : , : , :); vaxResult{j}.newImmHpvVax(2 : end , : , : , : , : , :)];
        vaxResult{j}.newHpvNonVax = [curr.newHpvNonVax(1 : end , : , : , : , : , :); vaxResult{j}.newHpvNonVax(2 : end , : , : , : , : , :)];
        vaxResult{j}.newImmHpvNonVax = [curr.newImmHpvNonVax(1 : end , : , : , : , : , :); vaxResult{j}.newImmHpvNonVax(2 : end , : , : , : , : , :)];
        vaxResult{j}.newScreen = [vaxResult{j}.newScreen(1 : end , : , : , : , : , : , : , : , :)]; %curr.newScreen(1 : end , : , : , : , : , : , : , : , :); vaxResult{j}.newScreen(2 : end , : , : , : , : , : , : , : , :)]; 
        vaxResult{j}.newHiv = [curr.newHiv(1 : end , : , : , : , : , : , :); vaxResult{j}.newHiv(2 : end , : , : , : , : , : , :)];
        vaxResult{j}.hivDeaths = [curr.hivDeaths(1 : end , : , : , :); vaxResult{j}.hivDeaths(2 : end , : , : , :)];
%         vaxResult{j}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :); vaxResult{j}.artTreatTracker(2 : end , : , : , : , : , :)];
        vaxResult{j}.tVec = [curr.tVec(1 : end), vaxResult{j}.tVec(2 : end)];

    %     noVaxInd = nSims;
    %     noV = vaxResult{noVaxInd};
        tVec = vaxResult{j}.tVec;
        tVecYr = tVec(1 : stepsPerYear : end);
        

        %% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

        %% total population size over time in by gender 
        for g = 1 : gender 
            popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , g , 1 : age , 1 : risk));
            popSize(z, j , g, :) = sum(vaxResult{j}.popVec(: , popTot),2);
        end
        
        %% Female population proportion by 5-year age groups in 2020 and 2070 
        for t = 1 : popYearVecLength
            for a = 1 : age
                popAgeF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                popTotF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
                popPropF(z, j , t , a) = sum(vaxResult{j}.popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popAgeF),2) ./ ...
                    sum(vaxResult{j}.popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTotF),2);
            end
        end
%          
        % Population proportion by 10-year age groups in 2020 and 2070
        for t = 1 : popYearVecLength
            for aV = 1 : ageVecLength_cPopDist
                aGroup = ageVec_cPopDist{aV};
                popAgeF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 1 : gender , aGroup , 1 : risk));
                popTotF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
                popPropBroadC(z, j , t , aV) = (sum(vaxResult{j}.popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popAgeF),2) ./ ...
                    sum(vaxResult{j}.popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTotF),2)) .* 100;
            end
        end
        
        % Female total population size by 5-year age groups over time
        for d = 1 : diseaseVecLength_ccInc
            dInds = diseaseVec_ccInc{d};
            for a = 1 : age
                popAgeF = toInd(allcomb(dInds , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                popSizeAgeF(z, j , d , a , :) = sum(vaxResult{j}.popVec(: , popAgeF),2);
            end
        end
        
        %% Risk distribution of HIV-negative women with CIN3 over time
        for rInd = 1 : risk
%             r = rInd;
            popRiskF = unique([toInd(allcomb(1 : 2 , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 4 : age , rInd)); toInd(allcomb(1 : 2 , 1 : viral , ...
                [1 : 5 , 7] , 5 , 1 , 1 : intervens , 2 , 4 : age , rInd))]);
            popAllF = unique([toInd(allcomb(1 : 2 , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 4 : age , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
                [1 : 5 , 7] , 5 , 1 , 1 : intervens , 2 , 4 : age , 1 : risk))]);
            popRiskDistCinF(z, j , rInd , :) = sum(vaxResult{j}.popVec(: , popRiskF),2) ./ sum(vaxResult{j}.popVec(: , popAllF),2);
        end
        
        % Risk distribution of HIV-negative women over time
        for rInd = 1 : risk
%             r = rInd;
            popRiskF = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , 5 : 7 , rInd));
            popAllF = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , 5 : 7 , 1 : risk)); 
            popRiskDistF(z, j , rInd , :) = sum(vaxResult{j}.popVec(: , popRiskF),2) ./ sum(vaxResult{j}.popVec(: , popAllF),2);
        end
        
        % Risk distribution of HIV-positive women over time
        for rInd = 1 : risk
%             r = rInd;
            popRiskF = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , 4 : age , rInd));
            popAllF = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , 4 : age , 1 : risk)); 
            popRiskDistHivF(z,j , rInd , :) = sum(vaxResult{j}.popVec(: , popRiskF),2) ./ sum(vaxResult{j}.popVec(: , popAllF),2);
        end
        
        % Proportion of each risk group that is HIV-positive over time
        for g = 1 : gender
            for rInd = 1 : risk
%                 r = rInd;
                popRiskF = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , g , 4 : age , rInd));
                popAllF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , g , 4 : age , rInd)); 
                popRiskDistHivProp(z, j , rInd , g , :) = sum(vaxResult{j}.popVec(: , popRiskF),2) ./ sum(vaxResult{j}.popVec(: , popAllF),2);
            end
        end
        
            
        %% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************
    
        %% HIV prevalence by age over time  
        for a = 1 : age
            hivMInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ... %untreated HIV+ males 
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            artMInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ... %treated HIV+ males
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            totMInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ... % all males
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            hivFInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ... %untreated HIV+ females
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            artFInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...  %treated HIV+ females
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            totFInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...   %all females
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            hivAgeM(z, j , a , :) =  (sum(vaxResult{j}.popVec(: , hivMInds) , 2) + sum(vaxResult{j}.popVec(: , artMInds) , 2)) ...
                ./ sum(vaxResult{j}.popVec(: , totMInds) , 2);
            hivAgeF(z, j , a , :) =  (sum(vaxResult{j}.popVec(: , hivFInds) , 2) + sum(vaxResult{j}.popVec(: , artFInds) , 2)) ...
                ./ sum(vaxResult{j}.popVec(: , totFInds) , 2);
        end
        
         %% HIV prevalence overall over time 
       
        hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 4 : 10 , 1 : risk));
        totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : gender , 4 : 10 , 1 : risk));
        hivPrev(z, j , :) = (sum(vaxResult{j}.popVec(: , hivInds) , 2) ./ sum(vaxResult{j}.popVec(: , totInds) , 2)) .* 100;

    
        %% HIV prevalence by gender over time 
        for g = 1 : gender
                hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                    1 : intervens , g , 4 : 10 , 1 : risk));
                totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                    1 : intervens , g , 4 : 10 , 1 : risk));
                hivPrevGen(z, j , g , :) = (sum(vaxResult{j}.popVec(: , hivInds) , 2) ./ sum(vaxResult{j}.popVec(: , totInds) , 2)) .* 100;
        end
        %% HIV prevalence by risk over time 
        for r = 1 : risk 
            hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                    1 : intervens , 2 , 4 : 10 , r));
            totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                    1 : intervens , 2 , 4 : 10 , r));  
            hivPrevRiskF(z, j, r, :) = (sum(vaxResult{j}.popVec(: , hivInds) , 2) ./ sum(vaxResult{j}.popVec(: , totInds) , 2)) .* 100;
        end
        %% HIV-associated deaths by gender over time
        hivDeathsM(z, j , :) = annlz(sum(sum(vaxResult{j}.hivDeaths(: , : , 1 , :), 2), 4));
        hivDeathsF(z, j , :) = annlz(sum(sum(vaxResult{j}.hivDeaths(: , : , 2 , :), 2), 4));
        
        %% HIV incidence over time 
        for g = 1 : gender
            hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens, g , 1 : age , 1: risk));
            hivInc(z, j, g, :) = (annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHiv(: , : , : , :, g , :, :), 2), 3), 4), 6), 7))) ./ ...
            (annlz(sum(vaxResult{j}.popVec(:, hivSusInds), 2) ./ stepsPerYear)) * 100;
        end
                
        %% HIV incidence by age over time 
        for a = 1 : age
           hivSusIndsM = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens, 1 , a , 1: risk));
           hivIncAgeM(z, j, a, :) = (annlz(sum(sum(sum(sum(vaxResult{j}.newHiv(: , : , : , :, 1 , a, :), 2), 3), 4), 7))) ./ ...
               (annlz(sum(vaxResult{j}.popVec(:, hivSusIndsM), 2) ./ stepsPerYear)) * 100;
           hivSusIndsF = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens, 2 , a , 1: risk));
           hivIncAgeF(z, j, a, :) = (annlz(sum(sum(sum(sum(vaxResult{j}.newHiv(: , : , : , :, 2 , a, :), 2), 3), 4), 7))) ./ ...
               (annlz(sum(vaxResult{j}.popVec(:, hivSusIndsF), 2) ./ stepsPerYear)) * 100;
        end
        
        %% HIV incidence by risk among women over time 
        for r = 1 : risk 
           hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens, 2 , 1 : age , r)); 
        hivIncRiskF(z, j, r, :) = (annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHiv(: , : , : , :, 2 , :, r), 2), 3), 4), 5), 6))) ./ ...
            (annlz(sum(vaxResult{j}.popVec(:, hivSusInds), 2) ./ stepsPerYear)) * 100;
        end
        %% cumulative number of HIV cases over time among women 
        for a = 1 : age
            for r = 1 : risk
                % Calculate incidence
                hivCumCasesF(z, j , a , r ,:) = ...
                    cumsum(squeeze(annlz(sum(sum(sum(sum(vaxResult{j}.newHiv(: , : , : , :, 2 , a, r), 2), 3), 4), 5))),2);
            end
        end
        
        for a = 1 : age
            for r = 1 : risk
                % Calculate incidence
                hivCumCasesM(z, j , a , r ,:) = ...
                    cumsum(squeeze(annlz(sum(sum(sum(sum(vaxResult{j}.newHiv(: , : , : , :, 1 , a, r), 2), 3), 4), 5))),2);
            end
        end
        
        %% Proportion of total HIV+ population on ART and VS (denominator: CD4-eligible and ineligible)
        artIndsF = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 3 : age , 1 : risk));
        hivAllIndsF = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
            1 : endpoints , 1 : intervens , 2 , 3 : age , 1 : risk));
        artIndsM = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , 3 : age , 1 : risk));
        hivAllIndsM = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
            1 : endpoints , 1 : intervens , 1 , 3 : age , 1 : risk));
        
        artCovF(z, j , :) = (sum(vaxResult{j}.popVec(: , artIndsF) , 2) ./ (sum(vaxResult{j}.popVec(: , hivAllIndsF) , 2) + sum(vaxResult{j}.popVec(: , artIndsF) , 2))) .* 100;
        artCovM(z, j , :) = (sum(vaxResult{j}.popVec(: , artIndsM) , 2) ./ (sum(vaxResult{j}.popVec(: , hivAllIndsM) , 2) + sum(vaxResult{j}.popVec(: , artIndsM) , 2))) .* 100;
        
        %% Proportion of total HIV+ population on ART and VS by age
        for a = 1 : age
            artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            artPop = sum(vaxResult{j}.popVec(: , artInds) , 2);
            hivInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            hivPop = sum(vaxResult{j}.popVec(: , hivInds) , 2);
            artCovAge(z, j , a , :) = artPop ./ hivPop;
        end
    
        
        %% ********************************** HPV FIGURES **********************************************************************************************
    
        %% Female HPV Prevalence by age and HIV status in 2005 and 2009 
        for i = 1 : hpvYearVecLength
            yr = hpvYearVec(i);          
            for a = 1 : ageVecLength_fHPV % 15-19 -> 55-65
                aGroup = ageVec_fHPV{a}
                hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , aGroup , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
                    [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , aGroup , 1 : risk))]);
                ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , aGroup , 1 : risk));
                hpv_hiv(z, j , a , i) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
                    ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    
                hpvInds_hivNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , aGroup , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
                    [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , aGroup , 1 : risk))]);
                ageInds_hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , aGroup , 1 : risk));
                hpv_hivNeg(z, j , a , i) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivNeg))...
                    ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivNeg));
            end
        end
        
%         %% Female HPV Prevalence by age 2018
%         for i = 1 : hpvYearVec2018Length
%             yr = hpvYearVec2018(i);
%             for a = 1 : age
%                 hpvInds_hivTot = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
%                     1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
%                     [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
%                 ageInds_hivTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
%                 hpv_hivTot(z, j , a , i) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivTot))...
%                     ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivTot));
%             end
%         end
        
        %% Female HPV Prevalence (not including cancer) by 5-year age groups 2020
%         for d = 1 : diseaseVecLength_fHpv
%             dInds = diseaseVec_fHpv{d};                 
%             for a = 4 : age          
%                 hpvInds_hivTot = unique([toInd(allcomb(dInds , 1 : viral ,  2 : 5 , [1 : 5 , 7], ...
%                     1 : 3 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(dInds , 1 : viral , ...
%                     [1 : 5 , 7] , 2 : 5 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk))]);
%                 ageInds_hivTot = toInd(allcomb(dInds , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
%                 hpv_hivTot2020(z, j , d , a) = sum(vaxResult{j}.popVec((2020 - startYear) * stepsPerYear +1 , hpvInds_hivTot))...
%                     ./ sum(vaxResult{j}.popVec((2020 - startYear) * stepsPerYear +1 , ageInds_hivTot));
%             end
%         end
%        
        %% Male HPV prevalence by age and HIV status in 2008 vs. Mbulawa data (calibration)
%         for i = 1 : hpvYearVecMaleLength
%             yr = hpvYearVecMale(i);
%             for a = 4 : age
%                 hpvInds_hivM = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 , [1 : 2 , 7] , ...
%                     1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
%                     [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
%                 ageInds_hivM = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
%                 hpv_hivM(z, j , a , i) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivM))...
%                     / sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear+1 , ageInds_hivM));
%     
%                 hpvInds_hivMNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 , [1 : 2 , 7] , ...
%                     1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
%                     [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
%                 ageInds_hivMNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
%                 hpv_hivMNeg(z, j , a , i) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivMNeg))...
%                     / sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivMNeg));
%                 
%                 hpvInds_hivMtot = unique([toInd(allcomb(1 : disease , 1 : viral , 2 , [1 : 2 , 7] , ...
%                     1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
%                     [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
%                 ageInds_hivMtot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
%                 hpv_hivMtot(z, j , a , i) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivMtot))...
%                     / sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivMtot));
%             end
%         end
        
        %% Female HPV Prevalence over time by HIV status
        hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 4 : 16 , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 4 : 16 , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 16 , 1 : risk));
        hpv_hivTimeF(z, j , :) = sum(vaxResult{j}.popVec(: , hpvInds) , 2) ./ sum(vaxResult{j}.popVec(: , ageInds) , 2);
    
        hpvInds_hivNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 4 : 16 , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , 4 : 16 , 1 : risk))]);
        ageInds_hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : 16 , 1 : risk));
        hpv_hivNegTimeF(z, j , :) = sum(vaxResult{j}.popVec(: , hpvInds_hivNeg) , 2) ./ sum(vaxResult{j}.popVec(: , ageInds_hivNeg) , 2);
    
        %% HPV Prevalence over time by sex       
        for g = 1 : gender
            hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , g , 4 : 16 , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
                [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , g , 4 : 16 , 1 : risk))]);
            popInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , g , 4 : 16 , 1 : risk));
            hpv_time(z, j , g , :) = sum(vaxResult{j}.popVec(: , hpvInds) , 2)...
                ./ sum(vaxResult{j}.popVec(: , popInds) , 2);
        end
        
        %% HPV prevalence in high risk group and by HIV status in 2006    
        yr = 2010;
            for a = 1 : ageVecLength_fHPV % 15-19 -> 55-65
                aGroup = ageVec_fHPV{a};
                hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , aGroup , 3)); toInd(allcomb(3 : 8 , 1 : viral , ...
                    [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , aGroup , 3))]);
                ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , aGroup , 3));
                hpv_hiv_risk(z, j , a) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
                    ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    
                hpvInds_hivNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , aGroup , 3)); toInd(allcomb(1 : 2 , 1 : viral , ...
                    [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , aGroup , 3))]);
                ageInds_hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , aGroup , 3));
                hpv_hivNeg_risk(z, j , a) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivNeg))...
                    ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivNeg));
            end
        %% HPV prevalence ratios in 2005, 2018, 2019
%         for i = 1 : hpvPrevYearVecLength
%             yr = hpvPrevYearVec(i); 
%             for g = 1 : gender
%                 for dInd = 1 : diseaseVecLength_hpv
%                     d = diseaseVec_hpv{dInd};
%                     hpvInds = unique([toInd(allcomb(d , 1 : viral , 2 : 6 , 1 : hpvNonVaxStates , ...
%                         1 : 3 , 1 : intervens , g , 4 : 13 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
%                         1 : hpvVaxStates , 2 : 6 , 1 : 3 , 1 : intervens , g , 4 : 13 , 1 : risk))]);
%                     popInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                         1 : endpoints , 1 : intervens , g , 4 : 13 , 1 : risk));
%                     hpv_prev_ratios(z, j , g , dInd , i) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , hpvInds) , 2)...
%                         ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , popInds) , 2);
%                 end
%             end
%         end
    
%         %% HPV incidence over time
%         fac = 100;
%         % General population
%         allSusF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 , ...
%         1 : endpoints , 1 : intervens, 2 , 1 : age , 1 : risk)); ...
%         toInd(allcomb(1 : disease , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
%         1 : endpoints , 1 : intervens, 2 , 1 : age, 1 : risk))];
%         % Calculate incidence
%         hpvIncTime(z, j , :) = ...
%             ((annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvVax(: , 2 , 1 : disease , : , : , :),2),3),4),5),6)) + ...
%             annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvVax(: , 2 , 1 : disease , : , : , :),2),3),4),5),6)) + ...
%             annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvNonVax(: , 2 , 1 : disease , : , : , :),2),3),4),5),6)) ...
%             ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvNonVax(: , 2 , 1 : disease , : , : , :),2),3),4),5),6))
%             (annlz(sum(vaxResult{j}.popVec(: , allSusF) , 2) ./ stepsPerYear)) * fac);
%     
%         % HIV-negative
%         allFneg = [toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1  , ...
%             1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk)); ... 
%             toInd(allcomb(1 : 2 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
%             1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk))]
%         % Calculate incidence
%         hpvIncTimeNeg(z, j , :) = ...
%             ((annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvVax(: , 2 , 1 : 2 , : , : , :),2),3),4),5),6)) + ...
%             annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvVax(: , 2 , 1 : 2 , : , : , :),2),3),4),5),6)) + ...
%             annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvNonVax(: , 2 , 1 : 2 , : , : , :),2),3),4),5),6)) ... 
%             ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvNonVax(: , 2 , 1 : 2 , : , : , :),2),3),4),5),6))
%             (annlz(sum(vaxResult{j}.popVec(: , allFneg) , 2) ./ stepsPerYear)) * fac);
%     
%         % HIV-positive untreated
%         allFpos = [toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1  , ...
%             1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
%             toInd(allcomb(3 : 7 , 1 : viral , 1  , 1 : hpvNonVaxStates , ...
%             1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk))];
%         % Calculate incidence
%         hpvIncTimePos(z, j , :) = ...
%             ((annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvVax(: , 2 , 3 : 7 , : , : , :),2),3),4),5),6)) + ...
%             annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvVax(: , 2 , 3 : 7 , : , : , :),2),3),4),5),6)) + ...
%             annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvNonVax(: , 2 , 3 : 7 , : , : , :),2),3),4),5),6)) ...
%             ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvNonVax(: , 2 , 3 : 7 , : , : , :),2),3),4),5),6))
%             (annlz(sum(vaxResult{j}.popVec(: , allFpos) , 2) ./ stepsPerYear)) * fac);
%     
%         % HIV-positive on ART
%         allFart = [toInd(allcomb(8 , 1 : viral , 1 : hpvVaxStates , 1  , ...
%             1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
%             toInd(allcomb(8 , 1 : viral , 1  , 1 : hpvNonVaxStates , ...
%             1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk))] ;
%         % Calculate incidence
%         hpvIncTimeArt(z, j , :) = ...
%             ((annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvVax(: , 2 , 8 , : , : , :),2),3),4),5),6)) + ...
%             annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvVax(: , 2 , 8 , : , : , :),2),3),4),5),6)) + ...
%             annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvNonVax(: , 2 , 8 , : , : , :),2),3),4),5),6)) ...
%             ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvNonVax(: , 2 , 8 , : , : , :),2),3),4),5),6))
%             (annlz(sum(vaxResult{j}.popVec(: , allFart) , 2) ./ stepsPerYear)) * fac);
%         
%         % HIV-positive all
%         allFposAll = [toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1  , ...
%             1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
%             toInd(allcomb(3 : 8 , 1 : viral , 1 , 1 : hpvNonVaxStates , ...
%             1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk))];
%         % Calculate incidence
%         hpvIncTimePosAll(z, j , :) = ...
%             ((annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvVax(: , 2 , 3 : 8 , : , : , :),2),3),4),5),6)) + ...
%             annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvVax(: , 2 , 3 : 8 , : , : , :),2),3),4),5),6)) + ...
%             annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvNonVax(: , 2 , 3 : 8 , : , : , :),2),3),4),5),6)) ...
%             ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvNonVax(: , 2 , 3 : 8 , : , : , :),2),3),4),5),6))
%             (annlz(sum(vaxResult{j}.popVec(: , allFposAll) , 2) ./ stepsPerYear)) * fac);
%         
%         %% HPV incidence by HIV status and age over time
%         fac = 100;
%         for d = 1 : diseaseVecLength_ccInc
%             dInds = diseaseVec_ccInc{d};
%             for rInd = 1 : risk
%                 riskInds = rInd;
%                 for a = 1 : age
%                     allF = toInd(allcomb(dInds , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                         1 : endpoints , 1 : intervens , 2 , a , riskInds));
%                     % Calculate incidence
%                     hpvIncHivRiskAgeTime(z, j , d , rInd , a , :) = ...
%                         ((annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvVax(: , 2 , dInds , a , riskInds , :),2),3),4),5),6)) + ...
%                         annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvVax(: , 2 , dInds , a , riskInds , :),2),3),4),5),6)) + ...
%                         annlz(sum(sum(sum(sum(sum(vaxResult{j}.newHpvNonVax(: , 2 , dInds , a , riskInds , :),2),3),4),5),6)) ...
%                         ) ./ ... %annlz(sum(sum(sum(sum(sum(vaxResult{j}.newImmHpvNonVax(: , 2 , dInds , a , r , :),2),3),4),5),6))
%                         (annlz(sum(vaxResult{j}.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);
%                 end
%             end
%         end
% 
%         
        %% ********************************** CIN FIGURES *********************************************************************************************
    
        %% CIN 1,2,3 prevalence by HIV status in 2010 
       
            yr = 2010;
            cInd = [3, 4, 5]; 
            for cin = 1 : 3
                cinlvl = cInd(cin); 
                % date is 15-19 -> 60-64
                % HIV-positive (on and not on ART)
                cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , cinlvl , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , 5 : 12 , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
                    [1 : 5 , 7] , cinlvl , 1 , 1 : intervens , 2 , 5 : 12 , 1 : risk))]);
                ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , 5 : 12 , 1 : risk));
                cinPosAge(z, j , cin) = (sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
                    ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageInds));
                % HIV-negative
                cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , cinlvl , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , 5 : 12 , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
                    [1 : 5 , 7] , cinlvl , 1 , 1 : intervens , 2 , 5 : 12 , 1 : risk))]);
                ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , 5 : 12 , 1 : risk));
                cinNegAge(z, j , cin) = (sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , cinNegInds)))...
                    ./ (sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageNegInds)));
                % General
%                 cinGenInds = unique([toInd(allcomb(1 : disease , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
%                     1 , 1 : intervens , 2 , 5 : 12 , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
%                     [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 5 : 12 , 1 : risk))]);
%                 ageGenInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
%                 cinGenAge(z, j , i , a) = (sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , cinGenInds)))...
%                     ./ (sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageGenInds)));
            end
          
        %% CIN1 and 2/3 prevalence by age for all HIV status combined in 2000
        yr = 2000;
        for a = 1 : adultsAgeVec_length
            aGroup = adultsAgeVec{a};
        % CIN 1 by age 
        cinGenInds = unique([toInd(allcomb(1 : disease , 1 : viral , 1 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , aGroup , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 5 , 7] , 1 , 1 , 1 : intervens , 2 , aGroup , 1 : risk))]);
        ageGenInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , aGroup , 1 : risk));
        cin1Age(z, j , a) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , cinGenInds) , 2)...
            ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageGenInds) , 2);
        % CIN 2/3 by age
        cinInds = unique([toInd(allcomb(1 : disease , 1 : viral , 4:5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , aGroup , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , aGroup , 1 : risk))]);
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , aGroup, 1 : risk));
        cin3Age(z, j , a) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , cinInds) , 2)...
            ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ageInds) , 2);
        end
        
        %% CIN2/3 prevalence by broad age groups 2020
        for d = 1 : diseaseVecLength_fHpv
            dInds = diseaseVec_fHpv{d};                 
            for aV = 1 : ageVecLength_fHPV2
                aGroup = ageVec_fHPV2{aV};
                cinInds = unique([toInd(allcomb(dInds , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , aGroup , 1 : risk)); toInd(allcomb(dInds , 1 : viral , ...
                    [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , aGroup , 1 : risk))]);
                ageInds = toInd(allcomb(dInds , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , aGroup , 1 : risk));
                cin_hivTot2020(z, j , d , aV) = (sum(vaxResult{j}.popVec((2020 - startYear) * stepsPerYear +1 , cinInds)))...
                    ./ (sum(vaxResult{j}.popVec((2020 - startYear) * stepsPerYear +1 , ageInds)));
            end
        end
               
        %% CIN2/3 prevalence ratios in 2005, 2018, 2019
%         for i = 1 : hpvPrevYearVecLength
%             yr = hpvPrevYearVec(i);  
%             for dInd = 1 : diseaseVecLength_hpv
%                 d = diseaseVec_hpv{dInd};
%                 cinInds = unique([toInd(allcomb(d , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
%                     1 , 1 : intervens , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
%                     [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 4 : 13 , 1 : risk))]);
%                 popInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk));
%                 cin_prev_ratios(z, j , dInd , i) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , cinInds) , 2)...
%                     ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , popInds) , 2);
%             end
%         end
    
        
        %% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************
    
        %% Cervical cancer incidence in 2005, 2012, 2018 by age vs. Globocan data 
        fac = 10 ^ 5;
        allF = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        for i = 1 : ccYearVecLength
            yr = ccYearVec(i);
            incTimeSpan = [((yr - startYear) * stepsPerYear +1) : ((yr - startYear) * stepsPerYear +6)];
            for a = 1 : age
                % All  
                allF = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7], [1 : 5 , 7] , ...
                    1  , 1 : intervens , 2 , a , 1 : risk));
                ccIncAge_all(z, j , i , a) = ...
                    (annlz(sum(sum(sum(vaxResult{j}.newCC(incTimeSpan , : , a , :),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{j}.popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);
               
                % HIV negative
                ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
                    1 : intervens , 2 , a , 1 : risk));
                ccIncAge_neg(z, j , i , a) = ...
                    (annlz(sum(sum(sum(vaxResult{j}.newCC(incTimeSpan , 1:2 , a , :),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{j}.popVec(incTimeSpan , ageNegInds) , 2) ./ stepsPerYear)) * fac);
                
                % HIV pos not on ART
                agePosInds = toInd(allcomb(3 : 7 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
                    1 : intervens , 2 , a , 1 : risk));
                ccIncAge_pos(z, j , i , a) = ...
                    (annlz(sum(sum(sum(vaxResult{j}.newCC(incTimeSpan , 3:7 , a , :),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{j}.popVec(incTimeSpan , agePosInds) , 2) ./ stepsPerYear)) * fac);
                
                %HIV pos on ART 
                ageArtInds = toInd(allcomb(8 , 6 , [1 : 5 , 7] , [1 : 5 , 7] , 1 , ...
                    1 : intervens , 2 , a , 1 : risk));
                ccIncAge_art(z, j , i , a) = ...
                    (annlz(sum(sum(sum(vaxResult{j}.newCC(incTimeSpan , 8 , a , :),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{j}.popVec(incTimeSpan , ageArtInds) , 2) ./ stepsPerYear)) * fac);
            end
        end
        
        %% Cervical cancer incidence over time - all ages, not standardized
        fac = 10 ^ 5;
        % General population
        allF = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        ccIncTime(z, j , :) = ...
            (annlz(sum(sum(sum(vaxResult{j}.newCC(: , : , 1 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{j}.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-negative
        allFneg = toInd(allcomb(1 : 2 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        ccIncTimeNeg(z, j , :) = ...
            (annlz(sum(sum(sum(vaxResult{j}.newCC(: , 1 : 2 , 1 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{j}.popVec(: , allFneg) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-positive untreated
        allFpos = toInd(allcomb(3 : 7 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        ccIncTimePos(z, j , :) = ...
            (annlz(sum(sum(sum(vaxResult{j}.newCC(: , 3 : 7 , 1 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{j}.popVec(: , allFpos) , 2) ./ stepsPerYear)) * fac);
    
        % HIV-positive on ART
        allFart = toInd(allcomb(8 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        ccIncTimeArt(z, j , :) = ...
            (annlz(sum(sum(sum(vaxResult{j}.newCC(: , 8 , 1 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{j}.popVec(: , allFart) , 2) ./ stepsPerYear)) * fac);
        
        % HIV-positive all
        allFposAll = toInd(allcomb(3 : 8 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 : endpoints , 1 : intervens , 2 , 1 : age , 1 : risk));
        % Calculate incidence
        ccIncTimePosAll(z, j , :) = ...
            (annlz(sum(sum(sum(vaxResult{j}.newCC(: , 3 : 8 , 1 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{j}.popVec(: , allFposAll) , 2) ./ stepsPerYear)) * fac);
        
        %% Cervical cancer incidence by HIV status over time, all ages, Segi pop standardized
        fac = 10 ^ 5;
        for d = 1 : diseaseVecLength_ccInc
            dInds = diseaseVec_ccInc{d};
            ccIncRefVec = zeros( length(annualTimespan),1)';
        for a = 1 : age
            % General
            allF = toInd(allcomb(1 : disease, 1 : viral, [1 : 5, 7], [1 : 5, 7], ...
                1 , 1 : intervens , 2 , a , 1 : risk));
            % All HIV-negative women
            hivNeg = toInd(allcomb(1 : 2, 1 : viral, [1 : 5 , 7], [1 : 5 , 7], 1, 1 : intervens, 2, a, 1 : risk));
            % All HIV-positive women
            hivAllF = toInd(allcomb(3 : 8 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , a , 1 : risk));
            % HIV-positive women not on ART
            hivNoARTF = toInd(allcomb(3 : 7, 1 : viral, [1 : 5, 7], [1 : 5, 7], ...
                1, 1 : intervens, 2, a, 1 : risk));
            % Women on ART
            artF = toInd(allcomb(8 , 6 , [1 : 5 , 7] , [1 : 5 , 7] , 1 , 1 : intervens , 2 , a , 1 : risk));
            
            genArray = {allF , hivNeg , hivAllF, hivNoARTF , artF};

            % Calculate incidence
                ccIncRef = ...
                    (annlz(sum(sum(vaxResult{j}.newCC(: , dInds , a , :),2),4)) ./ ...
                    (annlz(sum(vaxResult{j}.popVec(: , genArray{d}) , 2) ./ stepsPerYear))* fac) .* segi(a);
%                 if (d == 5) && (a < 3) && (max(annlz(sum(sum(vaxResult{j}.newCC(: , dInds , a , :),2),4))) == 0.0)
%                     ccIncRef = zeros(length(vaxResult{j}.popVec(:),1))';
%                 end
            ccIncRefVec = ccIncRefVec + ccIncRef;
             ccIncStd(z, j, d, :) = ccIncRefVec;
        end      
        end
        
         %% Cervical cancer incidence by HIV status over time, all ages, WHO standardized
        fac = 10 ^ 5;
        for d = 1 : diseaseVecLength_ccInc
            dInds = diseaseVec_ccInc{d};
            ccIncRefVec = zeros( length(annualTimespan),1)';
        for a = 1 : age
            % General
            allF = toInd(allcomb(1 : disease, 1 : viral, [1 : 5, 7], [1 : 5, 7], ...
                1 , 1 : intervens , 2 , a , 1 : risk));
            % All HIV-negative women
            hivNeg = toInd(allcomb(1 : 2, 1 : viral, [1 : 5 , 7], [1 : 5 , 7], 1, 1 : intervens, 2, a, 1 : risk));
            % All HIV-positive women
            hivAllF = toInd(allcomb(3 : 8 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , a , 1 : risk));
            % HIV-positive women not on ART
            hivNoARTF = toInd(allcomb(3 : 7, 1 : viral, [1 : 5, 7], [1 : 5, 7], ...
                1, 1 : intervens, 2, a, 1 : risk));
            % Women on ART
            artF = toInd(allcomb(8 , 6 , [1 : 5 , 7] , [1 : 5 , 7] , 1 , 1 : intervens , 2 , a , 1 : risk));
            
            genArray = {allF , hivNeg , hivAllF, hivNoARTF , artF};

            % Calculate incidence
                ccIncRef = ...
                    (annlz(sum(sum(vaxResult{j}.newCC(: , dInds , a , :),2),4)) ./ ...
                    (annlz(sum(vaxResult{j}.popVec(: , genArray{d}) , 2) ./ stepsPerYear))* fac) .* who(a);
%                 if (d == 5) && (a < 3) && (max(annlz(sum(sum(vaxResult{j}.newCC(: , dInds , a , :),2),4))) == 0.0)
%                     ccIncRef = zeros(length(vaxResult{j}.popVec(:),1))';
%                 end
            ccIncRefVec = ccIncRefVec + ccIncRef;
            ccIncStd_who(z, j, d, :) = ccIncRefVec;
        end      
        end
        
        %% Cervical cancer incidence by HIV status and age over time
        fac = 10 ^ 5;
        for d = 1 : diseaseVecLength_ccInc;
            dInds = diseaseVec_ccInc{d};
            for a = 1 : age
                allF = toInd(allcomb(dInds , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                % Calculate incidence
                ccIncHivAgeTime(z, j , d , a , :) = ...
                    (annlz(sum(sum(vaxResult{j}.newCC(: , dInds , a , :),2),4)) ./ ...
                    (annlz(sum(vaxResult{j}.popVec(: , allF) , 2) ./ stepsPerYear)) * fac);
            end
        end
        
        %% Cumulative cervical cancer cases by HIV status and age over time
        for d = 1 : diseaseVecLength_ccInc
            dInds = diseaseVec_ccInc{d};
            for a = 1 : age
                % Calculate incidence
                ccCumHivAgeTime(z, j , d , a , :) = ...
                    cumsum(squeeze(annlz(sum(sum(vaxResult{j}.newCC(: , dInds , a , :),2),4))),2);
            end
        end
        
        %% Yearly new CC cases by HIV 
        for d = 1 : diseaseVecLength_ccInc
            dInds = diseaseVec_ccInc{d};         
                ccHivTime(z, j , d , :) = ...
                    squeeze(sum(sum(sum(vaxResult{j}.newCC(: , dInds , : , :),2),3),4));   
        end 

        %% Cervical cancer prevalence ratios in 2005, 2018, 2019
%         for i = 1 : hpvPrevYearVecLength
%             yr = hpvPrevYearVec(i);  
%             for dInd = 1 : diseaseVecLength_hpv
%                 d = diseaseVec_hpv{dInd};
%                 ccInds = unique([toInd(allcomb(d , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
%                     1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
%                     1 : hpvVaxStates , 6 , 1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk))]);
%                 popInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : endpoints , 1 : intervens , 2 , 4 : 13 , 1 : risk));
%                 cc_prev_ratios(z, j , dInd , i) = sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , ccInds) , 2)...
%                     ./ sum(vaxResult{j}.popVec((yr - startYear) * stepsPerYear +1 , popInds) , 2);
%             end
%         end
%         
       
        %% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************
        
        %% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)
        for dInd = 1 : diseaseIndsLength_typeDist
            dInds = diseaseInds_typeDist{dInd};
            ccInds_vax = toInd(allcomb(dInds , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
                1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
            ccInds_nonVax = toInd(allcomb(dInds , 1 : viral , [1 : 5 , 7] , 6 , ...
                1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
            ccInds_tot = unique([toInd(allcomb(dInds , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
                    1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(dInds , 1 : viral , ...
                    [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    
           cin23Inds_vax = [toInd(allcomb(dInds , 1 : viral , 4 , [1 : 4 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
                toInd(allcomb(dInds , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk))];
            cin23Inds_nonVax = [toInd(allcomb(dInds , 1 : viral , [1 : 3 , 7] , 4 , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
                toInd(allcomb(dInds , 1 : viral , [1 : 4 , 7] , 5 , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk))];
            cin23Inds_tot = unique([toInd(allcomb(dInds , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
                toInd(allcomb(dInds , 1 : viral , 4 , [1 : 4 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk)); 
                toInd(allcomb(dInds , 1 : viral , ...
                [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk)); ...
                toInd(allcomb(dInds , 1 : viral , ...
                [1 3 , 7] , 4 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);     
    
            cin3Inds_vax = toInd(allcomb(dInds , 1 : viral , 5 , [1 : 5 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk));
            cin3Inds_nonVax = toInd(allcomb(dInds , 1 : viral , [1 : 4 , 7] , 5 , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk));
            cin3Inds_tot = unique([toInd(allcomb(dInds , 1 : viral , 5 , [1 : 5 , 7] , ...
                    1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(dInds , 1 : viral , ...
                    [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    
            cin2Inds_vax = toInd(allcomb(dInds , 1 : viral , 4 , [1 : 4 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk));
            cin2Inds_nonVax = toInd(allcomb(dInds , 1 : viral , [1 : 3 , 7] , 4 , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk));
            cin2Inds_tot = unique([toInd(allcomb(dInds , 1 : viral , 4 , [1 : 4 , 7] , ...
                    1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(dInds, 1 : viral , ...
                    [1 : 3 , 7] , 4 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    
            cin1Inds_vax = toInd(allcomb(dInds , 1 : viral , 3 , [1 : 3 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk));
            cin1Inds_nonVax = toInd(allcomb(dInds , 1 : viral , [1 : 2 , 7] , 3 , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk));
            cin1Inds_tot = unique([toInd(allcomb(dInds , 1 : viral , 3 , [1 : 3 , 7] , ...
                    1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(dInds , 1 : viral , ...
                    [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    
            hpvInds_vax = toInd(allcomb(dInds , 1 : viral , 2 , [1 : 2 , 7] , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk));
            hpvInds_nonVax = toInd(allcomb(dInds, 1 : viral , [1 , 7] , 2 , ...
                1 , 1 : intervens , 2 , 1 : age , 1 : risk));
            hpvInds_tot = unique([toInd(allcomb(dInds , 1 : viral , 2 , [1 : 2 , 7] , ...
                    1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(dInds , 1 : viral , ...
                    [1 , 7] , 2 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    
            cc_vax(z, j , dInd , :) = sum(vaxResult{j}.popVec(: , ccInds_vax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , ccInds_tot) , 2);
            cc_nonVax(z, j , dInd , :) = sum(vaxResult{j}.popVec(: , ccInds_nonVax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , ccInds_tot) , 2);
   
            cin23_vax(z,j , dInd , :) = sum(vaxResult{j}.popVec(: , cin23Inds_vax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , cin23Inds_tot) , 2);
            cin23_nonVax(z,j , dInd , :) = sum(vaxResult{j}.popVec(: , cin23Inds_nonVax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , cin23Inds_tot) , 2);
    
            cin3_vax(z,j , dInd , :) = sum(vaxResult{j}.popVec(: , cin3Inds_vax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , cin3Inds_tot) , 2);
            cin3_nonVax(z,j , dInd , :) = sum(vaxResult{j}.popVec(: , cin3Inds_nonVax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , cin3Inds_tot) , 2);
    
            cin2_vax(z, j , dInd , :) = sum(vaxResult{j}.popVec(: , cin2Inds_vax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , cin2Inds_tot) , 2);
            cin2_nonVax(z, j , dInd , :) = sum(vaxResult{j}.popVec(: , cin2Inds_nonVax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , cin2Inds_tot) , 2);
    
            cin1_vax(z, j , dInd , :) = sum(vaxResult{j}.popVec(: , cin1Inds_vax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , cin1Inds_tot) , 2);
            cin1_nonVax(z, j , dInd , :) = sum(vaxResult{j}.popVec(: , cin1Inds_nonVax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , cin1Inds_tot) , 2);
    
            hpv_vax(z, j , dInd , :) = sum(vaxResult{j}.popVec(: , hpvInds_vax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , hpvInds_tot) , 2);
            hpv_nonVax(z, j , dInd , :) = sum(vaxResult{j}.popVec(: , hpvInds_nonVax) , 2)...
                ./ sum(vaxResult{j}.popVec(: , hpvInds_tot) , 2);
        end
        
        
        %% ************************** SCREENING & VACCINATION FIGURES *******************************************************************************
        
%         %% Screening "coverage"
%         allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : endpoints , 1 : intervens , 2 , 7 : age , 1 : risk));
%         % Calculate incidence
%         newScreenTime(z,j , :) = ...
%             annlz(sum(sum(sum(sum(sum(sum(sum(sum(vaxResult{j}.newScreen(: , : , : , : , : , : , : , : , :),2),3),4),5),6),7),8),9)) ./ ...
%             (annlz(sum(vaxResult{j}.popVec(((2020 - startYear) * stepsPerYear +1):end , allF) , 2) ./ stepsPerYear) * 0.1);
%         
        %% Screening coverage ages 35-39
%         allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : endpoints , 1 : intervens , 2 , 8 , 1 : risk));
%         screenF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : endpoints , [3 , 4] , 2 , 8 , 1 : risk));
%         
%         screenCovTime(z,j , :) = ...
%             sum(vaxResult{j}.popVec(((2020 - startYear) * stepsPerYear +1):end , screenF) , 2) ./ ...
%             sum(vaxResult{j}.popVec(((2020 - startYear) * stepsPerYear +1):end , allF) , 2);
        
        %% Total number of women screened annually by age and disease status
%         for dInd = 1 : diseaseVecLength_ccInc
%             d = diseaseVec_ccInc{dInd};
%             screenTotAnnual(z,j , dInd , :) = annlz(sum(sum(sum(sum(sum(sum(sum(sum(vaxResult{j}.newScreen(: , d , : , : , : , : , : , : , :),2),3),4),5),6),7),8),9));
%         end
        
        %% Vaccine coverage overall
        % Overall
        vaxInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , [2 , 4] , 2 , 3 : 16 , 1 : risk));
        popInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 3 : 16 , 1 : risk));
        vaxCoverage(z, j , :) = sum(vaxResult{j}.popVec(: , vaxInds) , 2) ./ sum(vaxResult{j}.popVec(: , popInds) , 2);
        
        %% Vaccine coverage by age
        for a = 1 : age
            vaxInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , [2 , 4] , 2 , a , 1 : risk));
            popInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            vaxCoverageAge(z, j , a , :) = sum(vaxResult{j}.popVec(: , vaxInds) , 2) ./ sum(vaxResult{j}.popVec(: , popInds) , 2);
        end
        
        %% Total number of women vaccinated by age and disease status
        for dInd = 1 : diseaseVecLength_ccInc
            dInds = diseaseVec_ccInc{dInd};
            for a = 1 : age
                vaxInds = toInd(allcomb(dInds , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , [2 , 4] , 2 , a , 1 : risk));
                vaxTotAge(z, j , a , dInd , :) = sum(vaxResult{j}.popVec(: , vaxInds) , 2);
            end
        end
        
        
    end
end
end

savDir = [pwd , '/HHCoM_Results/'];
filename = ['stochmMod_outputs_9v'];
save(fullfile(savDir, filename), 'monthlyTimespan', 'annualTimespan', 'futAnnualTimespan', ...
    'popSize', 'popSizeAgeF', 'popYearVec', 'popYearVecLength', 'adultsAgeVec', 'adultsAgeVec_length', ...
    'popPropF', 'ageVec_cPopDist', 'ageVecLength_cPopDist', 'popPropBroadC', ...
    'hivAgeM', 'hivAgeF', 'hivPrevGen', 'hivPrev', 'hivPrevRiskF', ...
    'hivDeathsM', 'hivDeathsF', 'hivInc', 'hivIncAgeM', 'hivIncAgeF', 'hivIncRiskF', 'hivCumCasesF', 'hivCumCasesM', 'artCovM', ...
    'artCovF', 'artCovM', 'hpvYearVec', 'hpvYearVec2018_orig', 'hpvYearVecLength', 'hpvYearVec2018Length', ... 
    'hpv_hiv' , 'hpv_hivNeg', 'diseaseVec_fHpv', 'diseaseVecLength_fHpv', 'ageVec_fHPV', 'ageVecLength_fHPV', ...
    'ageVec_fHPV2', 'ageVecLength_fHPV2', 'hpvYearVecMale', ...
    'hpvYearVecMaleLength', 'ageVec_mHPV', 'ageVecLength_mHPV', 'hpv_hivTimeF', ...
    'hpv_hivNegTimeF', 'hpv_time', 'hpvIncTime', 'hpvIncTimeNeg', 'hpvIncTimePos', 'hpvIncTimeArt', 'hpvIncTimePosAll', ...
    'cinPosAge', 'cinNegAge','ccYearVec', 'ccYearVecLength', ...
    'ccIncAge_all', 'ccIncAge_pos', 'ccIncAge_neg', 'ccIncAge_art', 'ccIncTime', 'ccIncTimeNeg', 'ccIncTimePos', 'ccIncTimeArt', ...
    'ccIncTimePosAll', 'ccCumHivAgeTime', 'diseaseVec_ccInc', 'diseaseVecLength_ccInc', 'segi', ...
    'who', 'ccIncStd', 'ccIncStd_who', 'diseaseVec_hpv', 'hpvPrevYearVec', 'diseaseVecLength_hpv', 'hpvPrevYearVecLength', ...
    'hpv_prev_ratios', 'cin_prev_ratios', 'cc_prev_ratios', 'cc_inc_ratios', 'diseaseInds_typeDist', ...
    'diseaseIndsLength_typeDist', 'cc_vax', 'cc_nonVax', 'cin23_vax', 'cin23_nonVax', 'cin3_vax', 'cin3_nonVax', 'cin2_vax', ...
    'cin2_nonVax', 'cin1_vax', 'cin1_nonVax', 'hpv_vax', 'hpv_nonVax',  ...
    'vaxCoverage', 'vaxCoverageAge', 'vaxTotAge')

%% *************************** Write output to excel  **********************************************************************************************

%% HIV and CC reduction 
 hivPrev_21 = zeros(nRuns, gender, length(monthlyTimespan));
 hivPrev_31 = hivPrev_21; 
 hivPrev_41 = hivPrev_21;
 hivPrev_51 = hivPrev_21;
 hivPrev_61 = hivPrev_21;
 hivPrev_71 = hivPrev_21;
 
 hivPrevRiskF_21 = zeros(nRuns, length(monthlyTimespan));
 hivPrevRiskF_31 = hivPrevRiskF_21;
 hivPrevRiskF_41 = hivPrevRiskF_21;
 hivPrevRiskF_51 = hivPrevRiskF_21;
 hivPrevRiskF_61 = hivPrevRiskF_21;
 hivPrevRiskF_71 = hivPrevRiskF_21;
 
 hivInc_21 = zeros(nRuns, gender, length(annualTimespan));
 hivInc_31 = hivInc_21;
 hivInc_41 = hivInc_21;
 hivInc_51 = hivInc_21;
 hivInc_61 = hivInc_21;
 hivInc_71 = hivInc_21;
 
 hivIncRisk_21 = zeros(nRuns, length(annualTimespan));
 hivIncRisk_31 = hivIncRisk_21;
 hivIncRisk_41 = hivIncRisk_21;
 hivIncRisk_51 = hivIncRisk_21;
 hivIncRisk_61 = hivIncRisk_21;
 hivIncRisk_71 = hivIncRisk_21;
 
 ccInc_21 = hivIncRisk_21;
 ccInc_31 = hivIncRisk_21;
 ccInc_41 = hivIncRisk_21;
 ccInc_51 = hivIncRisk_21;
 ccInc_61 = hivIncRisk_21;
 ccInc_71 = hivIncRisk_21;
 
 ccIncStd_21 = hivIncRisk_21;
 ccIncStd_31 = hivIncRisk_21;
 ccIncStd_41 = hivIncRisk_21;
 ccIncStd_51 = hivIncRisk_21;
 ccIncStd_61 = hivIncRisk_21;
 ccIncStd_71 = hivIncRisk_21;
 
 ccCases_21 = hivIncRisk_21;
 ccCases_31 = hivIncRisk_21;
 ccCases_41 = hivIncRisk_21;
 ccCases_51 = hivIncRisk_21;
 ccCases_61 = hivIncRisk_21;
 ccCases_71 = hivIncRisk_21;
 
 ccCasesPct_21 = hivIncRisk_21;
 ccCasesPct_31 = hivIncRisk_21;
 ccCasesPct_41 = hivIncRisk_21;
 ccCasesPct_51 = hivIncRisk_21;
 ccCasesPct_61 = hivIncRisk_21;
 ccCasesPct_71 = hivIncRisk_21;
 
 hivCases_21 = zeros(nRuns, length(annualTimespan));
 hivCases_31 = hivCases_21;
 hivCases_41 = hivCases_21;
 hivCases_51 = hivCases_21;
 hivCases_61 = hivCases_21;
 hivCases_71 = hivCases_21;
 
 hivCasesM_21 = hivCases_21;
 hivCasesM_31 = hivCases_21;
 hivCasesM_41 = hivCases_21;
 hivCasesM_51 = hivCases_21;
 hivCasesM_61 = hivCases_21;
 hivCasesM_71 = hivCases_21;
 
 hivCasesRisk_21 = zeros(nRuns, length(annualTimespan));
 hivCasesRisk_31 = hivCasesRisk_21;
 hivCasesRisk_41 = hivCasesRisk_21;
 hivCasesRisk_51 = hivCasesRisk_21;
 hivCasesRisk_61 = hivCasesRisk_21; 
 hivCasesRisk_71 = hivCasesRisk_21;
 

for j = 1 : nRuns
    for g = 1 :gender
        hivPrev_21(j, g, :) = (squeeze(hivPrevGen(1, j, g, :)) - squeeze(hivPrevGen(2, j, g, :))) ./ squeeze(hivPrevGen(1, j, g, :)) * 100;
        hivPrev_31(j, g, :) = (squeeze(hivPrevGen(1, j, g, :)) - squeeze(hivPrevGen(3, j, g, :))) ./ squeeze(hivPrevGen(1, j, g, :)) * 100;
        hivPrev_41(j, g, :) = (squeeze(hivPrevGen(1, j, g, :)) - squeeze(hivPrevGen(4, j, g, :))) ./ squeeze(hivPrevGen(1, j, g, :)) * 100;
        hivPrev_51(j, g, :) = (squeeze(hivPrevGen(1, j, g, :)) - squeeze(hivPrevGen(5, j, g, :))) ./ squeeze(hivPrevGen(1, j, g, :)) * 100;
        hivPrev_61(j, g, :) = (squeeze(hivPrevGen(1, j, g, :)) - squeeze(hivPrevGen(6, j, g, :))) ./ squeeze(hivPrevGen(1, j, g, :)) * 100;
        hivPrev_71(j, g, :) = (squeeze(hivPrevGen(1, j, g, :)) - squeeze(hivPrevGen(7, j, g, :))) ./ squeeze(hivPrevGen(2, j, g, :)) * 100;
    end
    hivPrevRiskF_21(j, :) = (squeeze(hivPrevRiskF(1, j, 3, :)) - squeeze(hivPrevRiskF(2, j, 3, :))) ./ squeeze(hivPrevRiskF(1, j, 3, :)) * 100;
    hivPrevRiskF_31(j, :) = (squeeze(hivPrevRiskF(1, j, 3, :)) - squeeze(hivPrevRiskF(3, j, 3, :))) ./ squeeze(hivPrevRiskF(1, j, 3, :)) * 100;
    hivPrevRiskF_41(j, :) = (squeeze(hivPrevRiskF(1, j, 3, :)) - squeeze(hivPrevRiskF(4, j, 3, :))) ./ squeeze(hivPrevRiskF(1, j, 3, :)) * 100;
    hivPrevRiskF_51(j, :) = (squeeze(hivPrevRiskF(1, j, 3, :)) - squeeze(hivPrevRiskF(5, j, 3, :))) ./ squeeze(hivPrevRiskF(1, j, 3, :)) * 100;
    hivPrevRiskF_61(j, :) = (squeeze(hivPrevRiskF(1, j, 3, :)) - squeeze(hivPrevRiskF(6, j, 3, :))) ./ squeeze(hivPrevRiskF(1, j, 3, :)) * 100;
    hivPrevRiskF_71(j, :) = (squeeze(hivPrevRiskF(1, j, 3, :)) - squeeze(hivPrevRiskF(7, j, 3, :))) ./ squeeze(hivPrevRiskF(1, j, 3, :)) * 100;
    
    for g = 1:gender
        hivInc_21(j, g, :) = squeeze(hivInc(2, j, g, :)) ./ squeeze(hivInc(1, j, g, :));
        hivInc_31(j, g, :) = squeeze(hivInc(3, j, g, :)) ./ squeeze(hivInc(1, j, g, :));
        hivInc_41(j, g, :) = squeeze(hivInc(4, j, g, :)) ./ squeeze(hivInc(1, j, g, :));
        hivInc_51(j, g, :) = squeeze(hivInc(5, j, g, :)) ./ squeeze(hivInc(1, j, g, :));
        hivInc_61(j, g, :) = squeeze(hivInc(6, j, g, :)) ./ squeeze(hivInc(1, j, g, :));
        hivInc_71(j, g, :) = squeeze(hivInc(7, j, g, :)) ./ squeeze(hivInc(1, j, g, :));
    end
    
    hivIncRisk_21(j, :) = squeeze(hivIncRiskF(2, j, 3, :)) ./ squeeze(hivIncRiskF(1, j, 3, :));
    hivIncRisk_31(j, :) = squeeze(hivIncRiskF(3, j, 3, :)) ./ squeeze(hivIncRiskF(1, j, 3, :));
    hivIncRisk_41(j, :) = squeeze(hivIncRiskF(4, j, 3, :)) ./ squeeze(hivIncRiskF(1, j, 3, :));
    hivIncRisk_51(j, :) = squeeze(hivIncRiskF(5, j, 3, :)) ./ squeeze(hivIncRiskF(1, j, 3, :));
    hivIncRisk_61(j, :) = squeeze(hivIncRiskF(6, j, 3, :)) ./ squeeze(hivIncRiskF(1, j, 3, :));
    hivIncRisk_71(j, :) = squeeze(hivIncRiskF(7, j, 3, :)) ./ squeeze(hivIncRiskF(1, j, 3, :));
    
    ccInc_21(j, :) = squeeze(ccIncTime(2, j, :)) ./ squeeze(ccIncTime(1, j, :));
    ccInc_31(j, :) = squeeze(ccIncTime(3, j, :)) ./ squeeze(ccIncTime(1, j, :));
    ccInc_41(j, :) = squeeze(ccIncTime(4, j, :)) ./ squeeze(ccIncTime(1, j, :));
    ccInc_51(j, :) = squeeze(ccIncTime(5, j, :)) ./ squeeze(ccIncTime(1, j, :));
    ccInc_61(j, :) = squeeze(ccIncTime(6, j, :)) ./ squeeze(ccIncTime(1, j, :));
    ccInc_71(j, :) = squeeze(ccIncTime(7, j, :)) ./ squeeze(ccIncTime(1, j, :));
    
    ccIncStd_21(j, :) = squeeze(ccIncStd(2, j, 1, :)) ./ squeeze(ccIncStd(1, j, 1, :));
    ccIncStd_31(j, :) = squeeze(ccIncStd(3, j, 1, :)) ./ squeeze(ccIncStd(1, j, 1, :));
    ccIncStd_41(j, :) = squeeze(ccIncStd(4, j, 1, :)) ./ squeeze(ccIncStd(1, j, 1, :));
    ccIncStd_51(j, :) = squeeze(ccIncStd(5, j, 1, :)) ./ squeeze(ccIncStd(1, j, 1, :));
    ccIncStd_61(j, :) = squeeze(ccIncStd(6, j, 1, :)) ./ squeeze(ccIncStd(1, j, 1, :));
    ccIncStd_71(j, :) = squeeze(ccIncStd(7, j, 1, :)) ./ squeeze(ccIncStd(1, j, 1, :));
    
    ccIncStdWho_21(j, :) = squeeze(ccIncStd_who(2, j, 1, :)) ./ squeeze(ccIncStd_who(1, j, 1, :));
    ccIncStdWho_31(j, :) = squeeze(ccIncStd_who(3, j, 1, :)) ./ squeeze(ccIncStd_who(1, j, 1, :));
    ccIncStdWho_41(j, :) = squeeze(ccIncStd_who(4, j, 1, :)) ./ squeeze(ccIncStd_who(1, j, 1, :));
    ccIncStdWho_51(j, :) = squeeze(ccIncStd_who(5, j, 1, :)) ./ squeeze(ccIncStd_who(1, j, 1, :));
    ccIncStdWho_61(j, :) = squeeze(ccIncStd_who(6, j, 1, :)) ./ squeeze(ccIncStd_who(1, j, 1, :));
    ccIncStdWho_71(j, :) = squeeze(ccIncStd_who(7, j, 1, :)) ./ squeeze(ccIncStd_who(1, j, 1, :));
    
    ccCases_21(j, :) = squeeze(sum(ccCumHivAgeTime(2, j, 1, :, :), 4)) - squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4));
    ccCases_31(j, :) = squeeze(sum(ccCumHivAgeTime(3, j, 1, :, :), 4)) - squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4));
    ccCases_41(j, :) = squeeze(sum(ccCumHivAgeTime(4, j, 1, :, :), 4)) - squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4));
    ccCases_51(j, :) = squeeze(sum(ccCumHivAgeTime(5, j, 1, :, :), 4)) - squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4));
    ccCases_61(j, :) = squeeze(sum(ccCumHivAgeTime(6, j, 1, :, :), 4)) - squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4));
    ccCases_71(j, :) = squeeze(sum(ccCumHivAgeTime(7, j, 1, :, :), 4)) - squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4));
    
    ccCasesPct_21(j, :) = ccCases_21(j, :) ./ squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4))' * 100;
    ccCasesPct_31(j, :) = ccCases_31(j, :) ./ squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4))' * 100;
    ccCasesPct_41(j, :) = ccCases_41(j, :) ./ squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4))' * 100;
    ccCasesPct_51(j, :) = ccCases_51(j, :) ./ squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4))' * 100;
    ccCasesPct_61(j, :) = ccCases_61(j, :) ./ squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4))' * 100;
    ccCasesPct_71(j, :) = ccCases_71(j, :) ./ squeeze(sum(ccCumHivAgeTime(1, j, 1, :, :), 4))' * 100;
    
    hivCases_21(j, :) = squeeze(sum(sum(hivCumCasesF(2, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, :, :), 3), 4));
    hivCases_31(j, :) = squeeze(sum(sum(hivCumCasesF(3, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, :, :), 3), 4));
    hivCases_41(j, :) = squeeze(sum(sum(hivCumCasesF(4, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, :, :), 3), 4));
    hivCases_51(j, :) = squeeze(sum(sum(hivCumCasesF(5, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, :, :), 3), 4));
    hivCases_61(j, :) = squeeze(sum(sum(hivCumCasesF(6, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, :, :), 3), 4));
    hivCases_71(j, :) = squeeze(sum(sum(hivCumCasesF(7, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, :, :), 3), 4));
    
    hivCasesM_21(j, :) = squeeze(sum(sum(hivCumCasesM(2, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesM(1, j, :, :, :), 3), 4));
    hivCasesM_31(j, :) = squeeze(sum(sum(hivCumCasesM(3, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesM(1, j, :, :, :), 3), 4));
    hivCasesM_41(j, :) = squeeze(sum(sum(hivCumCasesM(4, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesM(1, j, :, :, :), 3), 4));
    hivCasesM_51(j, :) = squeeze(sum(sum(hivCumCasesM(5, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesM(1, j, :, :, :), 3), 4));
    hivCasesM_61(j, :) = squeeze(sum(sum(hivCumCasesM(6, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesM(1, j, :, :, :), 3), 4));
    hivCasesM_71(j, :) = squeeze(sum(sum(hivCumCasesM(7, j, :, :, :), 3), 4)) - squeeze(sum(sum(hivCumCasesM(1, j, :, :, :), 3), 4));
    
    hivCasesRisk_21(j, :) = squeeze(sum(sum(hivCumCasesF(2, j, :, 3, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, 3, :), 3), 4));
    hivCasesRisk_31(j, :) = squeeze(sum(sum(hivCumCasesF(3, j, :, 3, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, 3, :), 3), 4));
    hivCasesRisk_41(j, :) = squeeze(sum(sum(hivCumCasesF(4, j, :, 3, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, 3, :), 3), 4));
    hivCasesRisk_51(j, :) = squeeze(sum(sum(hivCumCasesF(5, j, :, 3, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, 3, :), 3), 4));
    hivCasesRisk_61(j, :) = squeeze(sum(sum(hivCumCasesF(6, j, :, 3, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, 3, :), 3), 4));
    hivCasesRisk_71(j, :) = squeeze(sum(sum(hivCumCasesF(7, j, :, 3, :), 3), 4)) - squeeze(sum(sum(hivCumCasesF(1, j, :, 3, :), 3), 4));
    
end
%%
for g = 1 : gender
 pct_hivPrev_21(:, g, :) = prctile(hivPrev_21(:, g, 571:stepsPerYear:end), [25 50 75], 1);
 pct_hivPrev_31(:, g, :) = prctile(hivPrev_31(:, g, 571:stepsPerYear:end), [25 50 75], 1); 
 pct_hivPrev_41(:, g, :) = prctile(hivPrev_41(:, g, 571:stepsPerYear:end), [25 50 75], 1);
 pct_hivPrev_51(:, g, :) = prctile(hivPrev_51(:, g, 571:stepsPerYear:end), [25 50 75], 1);
 pct_hivPrev_61(:, g, :) = prctile(hivPrev_61(:, g, 571:stepsPerYear:end), [25 50 75], 1);
 pct_hivPrev_71(:, g, :) = prctile(hivPrev_71(:, g, 571:stepsPerYear:end), [25 50 75], 1);
end

 pct_hivPrevRisk_21 = prctile(hivPrevRiskF_21(:, 571:stepsPerYear:end), [25 50 75], 1);
 pct_hivPrevRisk_31 = prctile(hivPrevRiskF_31(:, 571:stepsPerYear:end), [25 50 75], 1);
 pct_hivPrevRisk_41 = prctile(hivPrevRiskF_41(:, 571:stepsPerYear:end), [25 50 75], 1);
 pct_hivPrevRisk_51 = prctile(hivPrevRiskF_51(:, 571:stepsPerYear:end), [25 50 75], 1);
 pct_hivPrevRisk_61 = prctile(hivPrevRiskF_61(:, 571:stepsPerYear:end), [25 50 75], 1);
 pct_hivPrevRisk_71 = prctile(hivPrevRiskF_71(:, 571:stepsPerYear:end), [25 50 75], 1);
 
 for g = 1 : gender
     pct_hivInc_21(:, g, :) = prctile(hivInc_21(:, g, 96:end), [25 50 75], 1);
     pct_hivInc_31(:, g, :) = prctile(hivInc_31(:, g, 96:end), [25 50 75], 1);
     pct_hivInc_41(:, g, :) = prctile(hivInc_41(:, g, 96:end), [25 50 75], 1);
     pct_hivInc_51(:, g, :) = prctile(hivInc_51(:, g, 96:end), [25 50 75], 1);
     pct_hivInc_61(:, g, :) = prctile(hivInc_61(:, g, 96:end), [25 50 75], 1);
     pct_hivInc_71(:, g, :) = prctile(hivInc_71(:, g, 96:end), [25 50 75], 1);
 end
 pct_hivIncRisk_21 = prctile(hivIncRisk_21(:, 96:end), [25 50 75], 1);
 pct_hivIncRisk_31 = prctile(hivIncRisk_31(:, 96:end), [25 50 75], 1);
 pct_hivIncRisk_41 = prctile(hivIncRisk_41(:, 96:end), [25 50 75], 1);
 pct_hivIncRisk_51 = prctile(hivIncRisk_51(:, 96:end), [25 50 75], 1);
 pct_hivIncRisk_61 = prctile(hivIncRisk_61(:, 96:end), [25 50 75], 1);
 pct_hivIncRisk_71 = prctile(hivIncRisk_71(:, 96:end), [25 50 75], 1);
 
 pct_ccInc_21 = prctile(ccInc_21(:, 96:end), [25 50 75], 1);
 pct_ccInc_31 = prctile(ccInc_31(:, 96:end), [25 50 75], 1);
 pct_ccInc_41 = prctile(ccInc_41(:, 96:end), [25 50 75], 1);
 pct_ccInc_51 = prctile(ccInc_51(:, 96:end), [25 50 75], 1);
 pct_ccInc_61 = prctile(ccInc_61(:, 96:end), [25 50 75], 1);
 pct_ccInc_71 = prctile(ccInc_71(:, 96:end), [25 50 75], 1);
 
 pct_ccIncStd_21 = prctile(ccIncStd_21(:, 96:end), [25 50 75], 1);
 pct_ccIncStd_31 = prctile(ccIncStd_31(:, 96:end), [25 50 75], 1);
 pct_ccIncStd_41 = prctile(ccIncStd_41(:, 96:end), [25 50 75], 1);
 pct_ccIncStd_51 = prctile(ccIncStd_51(:, 96:end), [25 50 75], 1);
 pct_ccIncStd_61 = prctile(ccIncStd_61(:, 96:end), [25 50 75], 1);
 pct_ccIncStd_71 = prctile(ccIncStd_71(:, 96:end), [25 50 75], 1);
 
 pct_ccIncStdWho_21 = prctile(ccIncStdWho_21(:, 96:end), [25 50 75], 1);
 pct_ccIncStdWho_31 = prctile(ccIncStdWho_31(:, 96:end), [25 50 75], 1);
 pct_ccIncStdWho_41 = prctile(ccIncStdWho_41(:, 96:end), [25 50 75], 1);
 pct_ccIncStdWho_51 = prctile(ccIncStdWho_51(:, 96:end), [25 50 75], 1);
 pct_ccIncStdWho_61 = prctile(ccIncStdWho_61(:, 96:end), [25 50 75], 1);
 pct_ccIncStdWho_71 = prctile(ccIncStdWho_71(:, 96:end), [25 50 75], 1);
 
 pct_ccCases_21 = prctile(ccCases_21(:, 96:end), [25 50 75], 1);
 pct_ccCases_31 = prctile(ccCases_31(:, 96:end), [25 50 75], 1);
 pct_ccCases_41 = prctile(ccCases_41(:, 96:end), [25 50 75], 1);
 pct_ccCases_51 = prctile(ccCases_51(:, 96:end), [25 50 75], 1);
 pct_ccCases_61 = prctile(ccCases_61(:, 96:end), [25 50 75], 1);
 pct_ccCases_71 = prctile(ccCases_71(:, 96:end), [25 50 75], 1);
 
 pct_ccPercent_21 = prctile(ccCasesPct_21(:, 96:end), [25 50 75], 1);
 pct_ccPercent_31 = prctile(ccCasesPct_31(:, 96:end), [25 50 75], 1);
 pct_ccPercent_41 = prctile(ccCasesPct_41(:, 96:end), [25 50 75], 1);
 pct_ccPercent_51 = prctile(ccCasesPct_51(:, 96:end), [25 50 75], 1);
 pct_ccPercent_61 = prctile(ccCasesPct_61(:, 96:end), [25 50 75], 1);
 pct_ccPercent_71 = prctile(ccCasesPct_71(:, 96:end), [25 50 75], 1);
 
 pct_hivCases_21 = prctile(hivCases_21(:, 96:end), [25 50 75], 1);
 pct_hivCases_31 = prctile(hivCases_31(:, 96:end), [25 50 75], 1);
 pct_hivCases_41 = prctile(hivCases_41(:, 96:end), [25 50 75], 1);
 pct_hivCases_51 = prctile(hivCases_51(:, 96:end), [25 50 75], 1);
 pct_hivCases_61 = prctile(hivCases_61(:, 96:end), [25 50 75], 1);
 pct_hivCases_71 = prctile(hivCases_71(:, 96:end), [25 50 75], 1);
 
 pct_hivCasesM_21 = prctile(hivCasesM_21(:, 96:end), [25 50 75], 1);
 pct_hivCasesM_31 = prctile(hivCasesM_31(:, 96:end), [25 50 75], 1);
 pct_hivCasesM_41 = prctile(hivCasesM_41(:, 96:end), [25 50 75], 1);
 pct_hivCasesM_51 = prctile(hivCasesM_51(:, 96:end), [25 50 75], 1);
 pct_hivCasesM_61 = prctile(hivCasesM_61(:, 96:end), [25 50 75], 1);
 pct_hivCasesM_71 = prctile(hivCasesM_71(:, 96:end), [25 50 75], 1);
 
 pct_hivCasesRisk_21 = prctile(hivCasesRisk_21(:, 96:end), [25 50 75], 1);
 pct_hivCasesRisk_31 = prctile(hivCasesRisk_31(:, 96:end), [25 50 75], 1);
 pct_hivCasesRisk_41 = prctile(hivCasesRisk_41(:, 96:end), [25 50 75], 1);
 pct_hivCasesRisk_51 = prctile(hivCasesRisk_51(:, 96:end), [25 50 75], 1);
 pct_hivCasesRisk_61 = prctile(hivCasesRisk_61(:, 96:end), [25 50 75], 1);
 pct_hivCasesRisk_71 = prctile(hivCasesRisk_71(:, 96:end), [25 50 75], 1);

 fname = [pwd, '\HHCoM_Results\Results_9v\Redux.xlsx'];
 Reduxtable1 = array2table([annualTimespan(96:end)' squeeze(pct_hivPrev_21(:, 1, :))' squeeze(pct_hivPrev_31(:, 1, :))' squeeze(pct_hivPrev_41(:, 1, :))' ...
     squeeze(pct_hivPrev_51(:, 1, :))' squeeze(pct_hivPrev_61(:, 1, :))' squeeze(pct_hivPrev_71(:, 1, :))' ]); 
 Reduxtable1.Properties.VariableNames(1:size(Reduxtable1, 2)) = {'Year', 'P25_hivPrevM_21', 'Med_hivPrevM_21', 'P75_hivPrevM_21', ...
     'P25_hivPrevM_31', 'Med_hivPrevM_31', 'P75_hivPrevM_31', 'P25_hivPrevM_41', 'Med_hivPrevM_41', 'P75_hivPrevM_41', ... 
     'P25_hivPrevM_51', 'Med_hivPrevM_51', 'P75_hivPrevM_51', 'P25_hivPrevM_61', 'Med_hivPrevM_61', 'P75_hivPrevM_61', ...
     'P25_hivPrevM_71', 'Med_hivPrevM_71', 'P75_hivPrevM_71'};
 sheet = 'hivPrevM'; 
 writetable(Reduxtable1, fname, 'Sheet', sheet);
 
 Reduxtable2 = array2table([annualTimespan(96:end)' squeeze(pct_hivPrev_21(:, 2, :))' squeeze(pct_hivPrev_31(:, 2, :))' squeeze(pct_hivPrev_41(:, 2, :))' ...
     squeeze(pct_hivPrev_51(:, 2, :))' squeeze(pct_hivPrev_61(:, 2, :))' squeeze(pct_hivPrev_71(:, 2, :))' ]); 
 Reduxtable2.Properties.VariableNames(1:size(Reduxtable2, 2)) = {'Year', 'P25_hivPrevF_21', 'Med_hivPrevF_21', 'P75_hivPrevF_21', ...
     'P25_hivPrevF_31', 'Med_hivPrevF_31', 'P75_hivPrevF_31', 'P25_hivPrevF_41', 'Med_hivPrevF_41', 'P75_hivPrevF_41', ... 
     'P25_hivPrevF_51', 'Med_hivPrevF_51', 'P75_hivPrevF_51', 'P25_hivPrevF_61', 'Med_hivPrevF_61', 'P75_hivPrevF_61', ...
     'P25_hivPrevF_71', 'Med_hivPrevF_71', 'P75_hivPrevF_71'};
 sheet = 'hivPrevF'; 
 writetable(Reduxtable2, fname, 'Sheet', sheet);
 
 Reduxtable3 = array2table([annualTimespan(96:end)' pct_hivPrevRisk_21' pct_hivPrevRisk_31' pct_hivPrevRisk_41' ...
     pct_hivPrevRisk_51' pct_hivPrevRisk_61' pct_hivPrevRisk_71']);
 Reduxtable3.Properties.VariableNames(1:size(Reduxtable3, 2)) = {'Year', 'P25_hivPrevRisk_21', 'Med_hivPrevRisk_21', 'P75_hivPrevRisk_21', ...
     'P25_hivPrevRisk_31', 'Med_hivPrevRisk_31', 'P75_hivPrevRisk_31', 'P25_hivPrevRisk_41', 'Med_hivPrevRisk_41', 'P75_hivPrevRisk_41', ... 
     'P25_hivPrevRisk_51', 'Med_hivPrevRisk_51', 'P75_hivPrevRisk_51', 'P25_hivPrevRisk_61', 'Med_hivPrevRisk_61', 'P75_hivPrevRisk_61', ...
     'P25_hivPrevRisk_71', 'Med_hivPrevRisk_71', 'P75_hivPrevRisk_71'};
 sheet = 'hivPrevRisk'; 
 writetable(Reduxtable3, fname, 'Sheet', sheet);
 
 Reduxtable4 = array2table([annualTimespan(96:end)' squeeze(pct_hivInc_21(:, 1, :))' squeeze(pct_hivInc_31(:, 1, :))' squeeze(pct_hivInc_41(:, 1, :))' ...
     squeeze(pct_hivInc_51(:, 1, :))' squeeze(pct_hivInc_61(:, 1, :))' squeeze(pct_hivInc_71(:, 1, :))' ]);
 Reduxtable4.Properties.VariableNames(1:size(Reduxtable4, 2)) = {'Year', 'P25_hivIncM_21', 'Med_hivIncM_21', 'P75_hivIncM_21', ...
     'P25_hivIncM_31', 'Med_hivIncM_31', 'P75_hivIncM_31', 'P25_hivIncM_41', 'Med_hivIncM_41', 'P75_hivIncM_41', ... 
     'P25_hivIncM_51', 'Med_hivIncM_51', 'P75_hivIncM_51', 'P25_hivIncM_61', 'Med_hivIncM_61', 'P75_hivIncM_61', ...
     'P25_hivIncM_71', 'Med_hivIncM_71', 'P75_hivIncM_71'};
 sheet = 'hivIncM'; 
 writetable(Reduxtable4, fname, 'Sheet', sheet);
 
 Reduxtable5 = array2table([annualTimespan(96:end)' squeeze(pct_hivInc_21(:, 2, :))' squeeze(pct_hivInc_31(:, 2, :))' squeeze(pct_hivInc_41(:, 2, :))' ...
     squeeze(pct_hivInc_51(:, 2, :))' squeeze(pct_hivInc_61(:, 2, :))'  squeeze(pct_hivInc_71(:, 2, :))'  ]);
 Reduxtable5.Properties.VariableNames(1:size(Reduxtable5, 2)) = {'Year', 'P25_hivIncF_21', 'Med_hivIncF_21', 'P75_hivIncF_21', ...
     'P25_hivIncF_31', 'Med_hivIncF_31', 'P75_hivIncF_31', 'P25_hivIncF_41', 'Med_hivIncF_41', 'P75_hivIncF_41', ... 
     'P25_hivIncF_51', 'Med_hivIncF_51', 'P75_hivIncF_51', 'P25_hivIncF_61', 'Med_hivIncF_61', 'P75_hivIncF_61', ...
     'P25_hivIncF_71', 'Med_hivIncF_71', 'P75_hivIncF_71'};
 sheet = 'hivIncF'; 
 writetable(Reduxtable5, fname, 'Sheet', sheet);
 
 Reduxtable6 = array2table([annualTimespan(96:end)' pct_hivIncRisk_21' pct_hivIncRisk_31' pct_hivIncRisk_41' ...
     pct_hivIncRisk_51' pct_hivIncRisk_61'  pct_hivIncRisk_71']);
 Reduxtable6.Properties.VariableNames(1:size(Reduxtable6, 2)) = {'Year', 'P25_hivIncRisk_21', 'Med_hivIncRisk_21', 'P75_hivIncRisk_21', ...
     'P25_hivIncRisk_31', 'Med_hivIncRisk_31', 'P75_hivIncRisk_31', 'P25_hivIncRisk_41', 'Med_hivIncRisk_41', 'P75_hivIncRisk_41', ... 
     'P25_hivIncRisk_51', 'Med_hivIncRisk_51', 'P75_hivIncRisk_51', 'P25_hivIncRisk_61', 'Med_hivIncRisk_61', 'P75_hivIncRisk_61', ...
     'P25_hivIncRisk_71', 'Med_hivIncRisk_71', 'P75_hivIncRisk_71'};
 sheet = 'hivIncRisk'; 
 writetable(Reduxtable6, fname, 'Sheet', sheet);
 
 Reduxtable7 = array2table([annualTimespan(96:end)' pct_ccInc_21' pct_ccInc_31' pct_ccInc_41' ...
     pct_ccInc_51' pct_ccInc_61' pct_ccInc_71' ]);
 Reduxtable7.Properties.VariableNames(1:size(Reduxtable7, 2)) = {'Year', 'P25_ccInc_21', 'Med_ccInc_21', 'P75_ccInc_21', ...
     'P25_ccInc_31', 'Med_ccInc_31', 'P75_ccInc_31', 'P25_ccInc_41', 'Med_ccInc_41', 'P75_ccInc_41', ... 
     'P25_ccInc_51', 'Med_ccInc_51', 'P75_ccInc_51', 'P25_ccInc_61', 'Med_ccInc_61', 'P75_ccInc_61', ...
     'P25_ccInc_71', 'Med_ccInc_71', 'P75_ccInc_71'};
 sheet = 'ccInc'; 
 writetable(Reduxtable7, fname, 'Sheet', sheet);

 Reduxtable8 = array2table([annualTimespan(96:end)' pct_ccIncStd_21' pct_ccIncStd_31' pct_ccIncStd_41' ...
     pct_ccIncStd_51' pct_ccIncStd_61' pct_ccIncStd_71']);
 Reduxtable8.Properties.VariableNames(1:size(Reduxtable8, 2)) = {'Year', 'P25_ccIncStd_21', 'Med_ccIncStd_21', 'P75_ccIncStd_21', ...
     'P25_ccIncStd_31', 'Med_ccIncStd_31', 'P75_ccIncStd_31', 'P25_ccIncStd_41', 'Med_ccIncStd_41', 'P75_ccIncStd_41', ... 
     'P25_ccIncStd_51', 'Med_ccIncStd_51', 'P75_ccIncStd_51', 'P25_ccIncStd_61', 'Med_ccIncStd_61', 'P75_ccIncStd_61', ...
     'P25_ccIncStd_71', 'Med_ccIncStd_71', 'P75_ccIncStd_71'};
 sheet = 'ccIncStd'; 
 writetable(Reduxtable8, fname, 'Sheet', sheet);
 
 Reduxtable9 = array2table([annualTimespan(96:end)' pct_ccCases_21' pct_ccCases_31' pct_ccCases_41' ...
     pct_ccCases_51' pct_ccCases_61' pct_ccCases_71']);
 Reduxtable9.Properties.VariableNames(1:size(Reduxtable9, 2)) = {'Year', 'P25_ccCases_21', 'Med_ccCases_21', 'P75_ccCases_21', ...
     'P25_ccCases_31', 'Med_ccCases_31', 'P75_ccCases_31', 'P25_ccCases_41', 'Med_ccCases_41', 'P75_ccCases_41', ... 
     'P25_ccCases_51', 'Med_ccCases_51', 'P75_ccCases_51', 'P25_ccCases_61', 'Med_ccCases_61', 'P75_ccCases_61', ...
     'P25_ccCases_71', 'Med_ccCases_71', 'P75_ccCases_71'};
 sheet = 'ccCases'; 
 writetable(Reduxtable9, fname, 'Sheet', sheet);
 
 Reduxtable10 = array2table([annualTimespan(96:end)' pct_ccPercent_21' pct_ccPercent_31' pct_ccPercent_41' ...
     pct_ccPercent_51' pct_ccPercent_61' pct_ccPercent_71']);
 Reduxtable10.Properties.VariableNames(1:size(Reduxtable10, 2)) = {'Year', 'P25_ccPercent_21', 'Med_ccPercent_21', 'P75_ccPercent_21', ...
     'P25_ccPercent_31', 'Med_ccPercent_31', 'P75_ccPercent_31', 'P25_ccPercent_41', 'Med_ccPercent_41', 'P75_ccPercent_41', ... 
     'P25_ccPercent_51', 'Med_ccPercent_51', 'P75_ccPercent_51', 'P25_ccPercent_61', 'Med_ccPercent_61', 'P75_ccPercent_61', ...
     'P25_ccPercent_71', 'Med_ccPercent_71', 'P75_ccPercent_71'};
 sheet = 'ccPercent'; 
 writetable(Reduxtable10, fname, 'Sheet', sheet);

 Reduxtable11 = array2table([annualTimespan(96:end)' pct_hivCases_21' pct_hivCases_31' pct_hivCases_41' ...
     pct_hivCases_51' pct_hivCases_61' pct_hivCases_71']);
 Reduxtable11.Properties.VariableNames(1:size(Reduxtable11, 2)) = {'Year', 'P25_hivCases_21', 'Med_hivCases_21', 'P75_hivCases_21', ...
     'P25_hivCases_31', 'Med_hivCases_31', 'P75_hivCases_31', 'P25_hivCases_41', 'Med_hivCases_41', 'P75_hivCases_41', ... 
     'P25_hivCases_51', 'Med_hivCases_51', 'P75_hivCases_51', 'P25_hivCases_61', 'Med_hivCases_61', 'P75_hivCases_61', ...
     'P25_hivCases_71', 'Med_hivCases_71', 'P75_hivCases_71'}; 
 sheet = 'hivCases'; 
 writetable(Reduxtable11, fname, 'Sheet', sheet);
  
 Reduxtable12 = array2table([annualTimespan(96:end)' pct_hivCasesRisk_21' pct_hivCasesRisk_31' pct_hivCasesRisk_41' ...
     pct_hivCasesRisk_51' pct_hivCasesRisk_61' pct_hivCasesRisk_71']);
 Reduxtable12.Properties.VariableNames(1:size(Reduxtable12, 2)) = {'Year', 'P25_hivCasesRisk_21', 'Med_hivCasesRisk_21', 'P75_hivCasesRisk_21', ...
     'P25_hivCasesRisk_31', 'Med_hivCasesRisk_31', 'P75_hivCasesRisk_31', 'P25_hivCasesRisk_41', 'Med_hivCasesRisk_41', 'P75_hivCasesRisk_41', ... 
     'P25_hivCasesRisk_51', 'Med_hivCasesRisk_51', 'P75_hivCasesRisk_51', 'P25_hivCasesRisk_61', 'Med_hivCasesRisk_61', 'P75_hivCasesRisk_61', ...
     'P25_hivCasesRisk_71', 'Med_hivCasesRisk_71', 'P75_hivCasesRisk_71'};
 sheet = 'hivCasesRisk'; 
 writetable(Reduxtable12, fname, 'Sheet', sheet);
 
  Reduxtable13 = array2table([annualTimespan(96:end)' pct_ccIncStdWho_21' pct_ccIncStdWho_31' pct_ccIncStdWho_41' ...
     pct_ccIncStdWho_51' pct_ccIncStdWho_61'  pct_ccIncStdWho_71']);
 Reduxtable13.Properties.VariableNames(1:size(Reduxtable13, 2)) = {'Year', 'P25_ccIncStdWho_21', 'Med_ccIncStdWho_21', 'P75_ccIncStdWho_21', ...
     'P25_ccIncStdWho_31', 'Med_ccIncStdWho_31', 'P75_ccIncStdWho_31', 'P25_ccIncStdWho_41', 'Med_ccIncStdWho_41', 'P75_ccIncStdWho_41', ... 
     'P25_ccIncStdWho_51', 'Med_ccIncStdWho_51', 'P75_ccIncStdWho_51', 'P25_ccIncStdWho_61', 'Med_ccIncStdWho_61', 'P75_ccIncStdWho_61', ...
     'P25_ccIncStdWho_71', 'Med_ccIncStdWho_71', 'P75_ccIncStdWho_71'};
 sheet = 'ccIncStdWho'; 
 writetable(Reduxtable13, fname, 'Sheet', sheet);
 
  Reduxtable14 = array2table([annualTimespan(96:end)' pct_hivCasesM_21' pct_hivCasesM_31' pct_hivCasesM_41' ...
     pct_hivCasesM_51' pct_hivCasesM_61' pct_hivCasesM_71']);
 Reduxtable14.Properties.VariableNames(1:size(Reduxtable14, 2)) = {'Year', 'P25_hivCasesM_21', 'Med_hivCasesM_21', 'P75_hivCasesM_21', ...
     'P25_hivCasesM_31', 'Med_hivCasesM_31', 'P75_hivCasesM_31', 'P25_hivCasesM_41', 'Med_hivCasesM_41', 'P75_hivCasesM_41', ... 
     'P25_hivCasesM_51', 'Med_hivCasesM_51', 'P75_hivCasesM_51', 'P25_hivCasesM_61', 'Med_hivCasesM_61', 'P75_hivCasesM_61', ...
     'P25_hivCasesM_71', 'Med_hivCasesM_71', 'P75_hivCasesM_71'}; 
 sheet = 'hivCasesM'; 
 writetable(Reduxtable14, fname, 'Sheet', sheet);

%% Overall cancer incidence over time 
for z = 1 : length(fileKey)
%mean_ccIncTime = mean(squeeze(ccIncTime(z, :, :)), 1, 'omitnan');
pct_ccIncTime = prctile(squeeze(ccIncTime(z, :, :)), [25 50 75], 1); %10th and 90th percentiles

%mean_ccIncTimeNeg =  mean(squeeze(ccIncTimeNeg(z, :, :)), 1, 'omitnan'); 
pct_ccIncTimeNeg = prctile(squeeze(ccIncTimeNeg(z, :, :)), [25 50 75], 1); 

%mean_ccIncTimePos = mean(squeeze(ccIncTimePos(z, :, :)), 1, 'omitnan') ;
pct_ccIncTimePos = prctile(squeeze(ccIncTimePos(z, :, :)), [25 50 75], 1);

%mean_ccIncTimeArt = mean(squeeze(ccIncTimeArt(z, :, :)), 1, 'omitnan');
pct_ccIncTimeArt = prctile(squeeze(ccIncTimeArt(z, :, :)), [25 50 75], 1);
%
CCtable = array2table([annualTimespan'  pct_ccIncTime'  pct_ccIncTimeNeg' ...
     pct_ccIncTimePos'  pct_ccIncTimeArt']);
CCtable.Properties.VariableNames(1:13) = {'Year', [fileKey{z},'_P25_ccInc'], [fileKey{z},'_Med_ccInc'], [fileKey{z},'_P75_ccInc'], ...
    [fileKey{z},'_P25_ccIncNeg'], [fileKey{z},'_Med_ccIncNeg'], [fileKey{z},'_P75_ccIncNeg'], ...
    [fileKey{z},'_P25_ccIncPos'], [fileKey{z},'_Med_ccIncPos'], [fileKey{z},'_P75_ccIncPos'], ...
    [fileKey{z},'_P25_ccIncArt'], [fileKey{z},'_Med_ccIncArt'], [fileKey{z},'_P75_ccIncArt']};
fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_ccInc.xlsx'];
sheet = [fileKey{z}, '_incByYr'];
writetable(CCtable, fname, 'Sheet', sheet);
end

% CC incidence by age in 2009, 2012, 2018, and 2020
mean_ccIncAge = zeros(3, 16);
pct_ccIncAge = zeros(3, 16, 3);
mean_ccIncAgePos = mean_ccIncAge;
pct_ccIncAgePos = pct_ccIncAge;
mean_ccIncAgeNeg = mean_ccIncAge;
pct_ccIncAgeNeg = pct_ccIncAge;
mean_ccIncAgeArt = mean_ccIncAge;
pct_ccIncAgeArt = pct_ccIncAge;
colname = cell(3, 12);
ageGroup = ["0-4" ; "5-9"; "10-14"; "15-19" ; "20-24"; "25-29"; "30-34"; 
    "35-39"; "40-44"; "45-49"; "50-54"; "55-59" ; "60-64"; "65-69"; "70-74"; "75-79"];
yr = {'2009', '2012', '2018', '2020'};
for z = 1 : length(fileKey)
    for y = 1 : length(ccYearVec)
       %mean_ccIncAge(y, :) = mean(squeeze(ccIncAge_all(z, :, y, :)), 1, 'omitnan')';
       pct_ccIncAge(y, :, :) = prctile(squeeze(ccIncAge_all(z, :, y, :)), [25 50 75], 1)';
       %mean_ccIncAgePos(y, :) = mean(squeeze(ccIncAge_pos(z, :, y, :)), 1, 'omitnan')';
       pct_ccIncAgePos(y, :, :) = prctile(squeeze(ccIncAge_pos(z, :, y, :)), [25 50 75], 1)';
       %mean_ccIncAgeNeg(y, :) = mean(squeeze(ccIncAge_neg(z, :, y, :)), 1, 'omitnan')';
       pct_ccIncAgeNeg( y, :, :) = prctile(squeeze(ccIncAge_neg(z, :, y, :)), [25 50 75], 1)';
       %mean_ccIncAgeArt(y, :) = mean(squeeze(ccIncAge_art(z, :, y, :)), 1, 'omitnan')';
       pct_ccIncAgeArt(y, :, :) = prctile(squeeze(ccIncAge_art(z, :, y, :)), [25 50 75], 1)';
    end
    table = [];
    for x = 1 : length(ccYearVec)
	table = [table squeeze(pct_ccIncAge(x, :, :)) squeeze(pct_ccIncAgePos(x, :, :)) ...
        squeeze(pct_ccIncAgeNeg(x, :, :)) squeeze(pct_ccIncAgeArt(x, :, :)) ];
    end  
    for b =  1 : length(ccYearVec) 
           colname(b, 1:3) = {[fileKey{z}, '_', yr{b}, '_P25_ccIncAge']; [fileKey{z}, '_', yr{b}, '_Med_ccIncAge']; [fileKey{z}, '_', yr{b}, '_P75_ccIncAge']};
           colname(b, 4:6) = {[fileKey{z}, '_',yr{b}, '_P25_ccIncAgePos']; [fileKey{z}, '_',yr{b}, '_Med_ccIncAgePos']; [fileKey{z}, '_',yr{b}, '_P75_ccIncAgePos']};
           colname(b, 7:9) = {[fileKey{z}, '_',yr{b}, '_P25_ccIncAgeNeg']; [fileKey{z}, '_',yr{b}, '_Med_ccIncAgeNeg']; [fileKey{z}, '_',yr{b}, '_P75_ccIncAgeNeg']};
           colname(b, 10:12) = {[fileKey{z}, '_',yr{b}, '_P25_ccIncAgeArt']; [fileKey{z}, '_',yr{b}, '_Med_ccIncAgeArt']; [fileKey{z}, '_',yr{b}, '_P75_ccIncAgeArt']};
    end
    CCtable2 = array2table([ageGroup table]);
    CCtable2.Properties.VariableNames(1:size(CCtable2, 2)) = {'Age', colname{1, :}, colname{2, :}, colname{3, :}, colname{4, :}};
    fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_ccInc.xlsx'];
    sheet = [fileKey{z}, '_incByAge'];
    writetable(CCtable2, fname, 'Sheet', sheet);
end


% Age standardized CC over time 
for z = 1 : length(fileKey)
    %mean_IncStdAll = mean(squeeze(ccIncStd(z, :, 1, :)), 1, 'omitnan');
    pct_IncStdAll = prctile(squeeze(ccIncStd(z, :, 1, :)), [25 50 75], 1);
    %mean_IncStdNeg = mean(squeeze(ccIncStd(z, :, 2, :)), 1, 'omitnan');
    pct_IncStdNeg = prctile(squeeze(ccIncStd(z, :, 2, :)), [25 50 75], 1);
   % mean_IncStdPos = mean(squeeze(ccIncStd(z, :, 3, :)), 1, 'omitnan');
    pct_IncStdPos = prctile(squeeze(ccIncStd(z, :, 3, :)), [25 50 75], 1);
   % mean_IncStdnoArt = mean(squeeze(ccIncStd(z, :, 4, :)), 1, 'omitnan');
    pct_IncStdnoArt = prctile(squeeze(ccIncStd(z, :, 4, :)), [25 50 75], 1);
    %mean_IncStdArt = mean(squeeze(ccIncStd(z, :, 5, :)), 1, 'omitnan');
    pct_IncStdArt = prctile(squeeze(ccIncStd(z, :, 5, :)), [25 50 75], 1);
    
    CCtable3 = array2table([annualTimespan' pct_IncStdAll'  pct_IncStdNeg' ...
        pct_IncStdPos' pct_IncStdnoArt'  pct_IncStdArt']); 
    CCtable3.Properties.VariableNames(1:size(CCtable3, 2)) = {'Year', [fileKey{z}, '_P25_IncStdAll'], [fileKey{z}, '_Med_IncStdAll'], ...
        [fileKey{z}, '_P75_IncStdAll'], ...
        [fileKey{z}, '_P25_IncStdNeg'], [fileKey{z}, '_Med_IncStdNeg'], [fileKey{z}, '_P75_IncStdNeg'], ...
        [fileKey{z}, '_P25_IncStdPos'], [fileKey{z}, '_Med_IncStdPos'], [fileKey{z}, '_P75_IncStdPos'], ...
        [fileKey{z}, '_P25_IncStdnoArt'], [fileKey{z}, '_Med_IncStdnoArt'], [fileKey{z}, '_P75_IncStdnoArt'], ...
        [fileKey{z}, '_P25_IncStdArt'], [fileKey{z}, '_Med_IncStdArt'], [fileKey{z}, '_P75_IncStdArt'] };
    fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_ccInc.xlsx'];
    sheet = [fileKey{z}, '_IncStd'];
    writetable(CCtable3, fname, 'Sheet', sheet);
end
% WHO standardized 
for z = 1 : length(fileKey)
    pct_IncStdWhoAll = prctile(squeeze(ccIncStd_who(z, :, 1, :)), [25 50 75], 1);
    pct_IncStdWhoNeg = prctile(squeeze(ccIncStd_who(z, :, 2, :)), [25 50 75], 1);
    pct_IncStdWhoPos = prctile(squeeze(ccIncStd_who(z, :, 3, :)), [25 50 75], 1);
    pct_IncStdWhonoArt = prctile(squeeze(ccIncStd_who(z, :, 4, :)), [25 50 75], 1);
    pct_IncStdWhoArt = prctile(squeeze(ccIncStd_who(z, :, 5, :)), [25 50 75], 1);
    
    CCtable4 = array2table([annualTimespan' pct_IncStdWhoAll'  pct_IncStdWhoNeg' ...
        pct_IncStdWhoPos' pct_IncStdWhonoArt'  pct_IncStdWhoArt']); 
    CCtable4.Properties.VariableNames(1:size(CCtable4, 2)) = {'Year', [fileKey{z}, '_P25_IncStdWhoAll'], [fileKey{z}, '_Med_IncStdWhoAll'], ...
        [fileKey{z}, '_P75_IncStdWhoAll'], ...
        [fileKey{z}, '_P25_IncStdWhoNeg'], [fileKey{z}, '_Med_IncStdWhoNeg'], [fileKey{z}, '_P75_IncStdWhoNeg'], ...
        [fileKey{z}, '_P25_IncStdWhoPos'], [fileKey{z}, '_Med_IncStdWhoPos'], [fileKey{z}, '_P75_IncStdWhoPos'], ...
        [fileKey{z}, '_P25_IncStdWhonoArt'], [fileKey{z}, '_Med_IncStdWhonoArt'], [fileKey{z}, '_P75_IncStdWhonoArt'], ...
        [fileKey{z}, '_P25_IncStdWhoArt'], [fileKey{z}, '_Med_IncStdWhoArt'], [fileKey{z}, '_P75_IncStdWhoArt'] };
    fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_ccInc.xlsx'];
    sheet = [fileKey{z}, '_IncStdWho'];
    writetable(CCtable4, fname, 'Sheet', sheet);
end

% Annual number of new CC cases
pct_ccCases_all = zeros(3, 146);
for z = 1 : length(fileKey)
   pct_ccCases_all = prctile(squeeze(sum(ccCumHivAgeTime(z, :, 2, :, :), 4)), [25 50 75], 1); 
   CCtable5 = array2table([annualTimespan' pct_ccCases_all']); 
   CCtable5.Properties.VariableNames(1:size(CCtable5, 2)) = {'Year', [fileKey{z}, '_P25_ccCases_all'], ...
       [fileKey{z}, '_Med_ccCases_all'], [fileKey{z}, '_P75_ccCases_all'] };
   fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_ccInc.xlsx'];
    sheet = [fileKey{z}, '_ccCasesYr'];
    writetable(CCtable5, fname, 'Sheet', sheet);
end

 %% Overall population over time 
 popSizeTot = zeros(length(fileKey), length(fileInds), length(monthlyTimespan));
 pct_popSizeTot = zeros(3, length(monthlyTimespan));
 for z = 1 : length(fileKey)
    popSizeTot(z, :, :) = squeeze(sum( popSize(z, : , :, :) , 3));
    mean_popSizeTot = mean(squeeze(popSizeTot(z, :, :)), 1, 'omitnan');
    pct_popSizeTot = prctile(squeeze(popSizeTot(z, :, :)), [25 50 75], 1); 
 
 Poptable = array2table([annualTimespan' pct_popSizeTot(:, 1 : stepsPerYear : end)' ]);
 Poptable.Properties.VariableNames(1:4) = {'Year', [fileKey{z},'_P25_pop'], ...
     [fileKey{z},'_Med_pop'], [fileKey{z},'_P75_pop']};
 fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_Pop.xlsx'];
 sheet = [fileKey{z}, '_popByYr'];
 writetable(Poptable, fname, 'Sheet', sheet);
 end
 
 % Female population by 10-year age groups in 1990-2020 and 2070
 ageGroup = ["0-9" ; "10-19" ; "20-29"; "30-39"; "40-49"; "50-59" ; "60+"];
 for z = 1 : length(fileKey) 
    pct_popPropF2010 = prctile(squeeze(popPropBroadC(z, :, 21, :)), [25 50 75], 1);
    
    mean_popPropF2020 = mean(squeeze(popPropBroadC(z, :, 31, :)), 1, 'omitnan');
    pct_popPropF2020 = prctile(squeeze(popPropBroadC(z, :, 31, :)), [25 50 75], 1);
    
    mean_popPropF2070 = mean(squeeze(popPropBroadC(z, :, 32, :)), 1, 'omitnan');
    pct_popPropF2070 = prctile(squeeze(popPropBroadC(z, :, 32, :)), [25 50 75], 1);
    
    Poptable2 = array2table([ageGroup  pct_popPropF2010' pct_popPropF2020' pct_popPropF2070']);
 
    Poptable2.Properties.VariableNames(1:size(Poptable2, 2)) = {'Age', [fileKey{z},'_P25_popPropF2010'], [fileKey{z},'_Med_popProp2010'], ...
     [fileKey{z},'_P75_popProp2010'], [fileKey{z},'_P25_popProp2020'], [fileKey{z},'_Med_popProp2020'], ...
     [fileKey{z},'_P75_popProp2020'], [fileKey{z},'_P25_popProp2070'], [fileKey{z},'_Med_popProp2070'], ...
     [fileKey{z},'_P75_popProp2070']};
 fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_Pop.xlsx'];
 sheet = [fileKey{z}, '_popProp'];
 writetable(Poptable2, fname, 'Sheet', sheet);
 end
 
 %% Vaccine coverage 
  % by 5- year age groups 
  ageGroup = ["10-14"; "15-19" ; "20-24"];
  pct_vaxCovAge = zeros(length(age), length(monthlyTimespan));
  for z = 1 : length(fileKey)
      for a = 1 : 3
        pct_vaxCovAge(a, :) = prctile(squeeze(vaxCoverageAge(z, :, a+2, :)), [50], 1); 
      end
      Vaxtable1 = array2table([annualTimespan(96:end)' squeeze(pct_vaxCovAge(:, 576 : stepsPerYear : end))']);
      
      Vaxtable1.Properties.VariableNames(1:size(Vaxtable1, 2)) = {'Year', [fileKey{z},'_Med_vaxCov_10_14'], ...
          [fileKey{z},'_Med_vaxCov_15_19'], [fileKey{z},'_Med_vaxCov_20_24']};
      fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_Vax.xlsx'];
      sheet = [fileKey{z}, '_vaxCovAge'];
      writetable(Vaxtable1, fname, 'Sheet', sheet);
  end
  
  % by total eligible pop size 
  vaxCovEligRout = zeros(length(fileKey), length(fileInds), length(monthlyTimespan));
  vaxCovEligCU = vaxCovEligRout;
  
  for z = 1 : length(fileKey)
      for j = 1 : length(fileInds)
         vaxCovEligRout(z, j, :) = squeeze(vaxTotAge(z, j, 3, 1, :) ./ popSizeAgeF(z, j, 1, 3, :) * 100);
         vaxCovEligCU(z, j, :) = squeeze(sum(vaxTotAge(z, j, 3:5, 1, :), 3) ./ sum(popSizeAgeF(z, j, 1, 3:5, :), 4) * 100); 
      end
  end
 
  for z = 1 : length(fileKey)
      pct_vaxCovEligRout = prctile(squeeze(vaxCovEligRout(z, :, :)), [50], 1);
      pct_vaxCovEligCU = prctile(squeeze(vaxCovEligCU(z, :, :)), [50], 1);
      
      Vaxtable2 = array2table([annualTimespan(96:end)' squeeze(pct_vaxCovEligRout(:, 576:stepsPerYear :end))' ... 
          squeeze(pct_vaxCovEligCU(:, 576:stepsPerYear:end))']);
      Vaxtable2.Properties.VariableNames(1:size(Vaxtable2, 2)) = {'Year', [fileKey{z}, '_Med_vaxCovEligRout'], ...
          [fileKey{z}, '_Med_vaxCovEligCU']};
      
      fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_Vax.xlsx'];
      sheet = [fileKey{z}, '_vaxCovElig'];
      writetable(Vaxtable2, fname, 'Sheet', sheet);
  end
  
  % vaccine coverage in the total population 
  for z = 1 : length(fileKey)
      pct_vaxCov = prctile(squeeze(vaxCoverage(z,:,:)), [50], 1);
      
      Vaxtable3 = array2table([futAnnualTimespan' squeeze(pct_vaxCov(:, 576:stepsPerYear :end))' ]);
      Vaxtable3.Properties.VariableNames(1:size(Vaxtable3, 2)) = {'Year', [fileKey{z}, '_Med_vaxCovTotPop']};
      fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_Vax.xlsx'];
      sheet = [fileKey{z}, '_vaxCovTot'];
      writetable(Vaxtable3, fname, 'Sheet', sheet);
   end
  
 %% Overall HIV and HIV by sex over year
 for z = 1 : length(fileKey)
     mean_hivPrevM = mean(squeeze(hivPrevGen(z, :, 1, :)), 1, 'omitnan');
     pct_hivPrevM = prctile(squeeze(hivPrevGen(z, :, 1, :)), [25 50 75], 1); 
     mean_hivPrevF = mean(squeeze(hivPrevGen(z, :, 2, :)), 1, 'omitnan');
     pct_hivPrevF = prctile(squeeze(hivPrevGen(z, :, 2, :)), [25 50 75], 1); 
     mean_hivPrev = mean(squeeze(hivPrev(z, :, :)), 1, 'omitnan');
     pct_hivPrev = prctile(squeeze(hivPrev(z, :, :)), [25 50 75]);
     
     HIVtable = array2table([annualTimespan(56:end)' ...
     pct_hivPrev(:, 331 : stepsPerYear : end)' pct_hivPrevM(:, 331 : stepsPerYear : end)'  ...
     pct_hivPrevF(:, 331 : stepsPerYear : end)' ]);
     HIVtable.Properties.VariableNames(1:10) = {'Year', [fileKey{z},'_P25_hiv'], [fileKey{z},'_Med_hiv'],  ...
     [fileKey{z},'_P75_hiv'], [fileKey{z},'_P25_hivM'], [fileKey{z},'_Med_hivM'], ...
     [fileKey{z},'_P75_hivM'], [fileKey{z},'_P25_hivF'], [fileKey{z},'_Med_hivF'], [fileKey{z},'_P75_hivF']};
     
     fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HIV.xlsx'];
     sheet = [fileKey{z}, '_HIVbySexYr'];
     writetable(HIVtable, fname, 'Sheet', sheet);
 end
 
 % HIV by age
 ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' ,...
    '60-64' , '65-69' , '70-74' , '75-79'};

 mean_hivAgeM = zeros(16, length(monthlyTimespan));
 pct_hivAgeM = zeros(3, 16, length(monthlyTimespan));
 mean_hivAgeF = mean_hivAgeM;
 pct_hivAgeF = pct_hivAgeM;
 colname = cell(6, length(age));
 for z = 1 : length(fileKey)
    for a = 1 : age
        mean_hivAgeM(a, :) = mean(squeeze(hivAgeM(z, :, a, :)), 1, 'omitnan');
        pct_hivAgeM(:, a, :) = prctile(squeeze(hivAgeM(z, :, a, :)), [25 50 75]);
        mean_hivAgeF(a, :) = mean(squeeze(hivAgeF(z, :, a, :)), 1, 'omitnan');
        pct_hivAgeF(:, a, :) = prctile(squeeze(hivAgeF(z, :, a, :)), [25 50 75]);
    end
        HIVtable2 = array2table([annualTimespan(56:end)'  squeeze(pct_hivAgeM(1, :, 331 : stepsPerYear : end))' ...
            squeeze(pct_hivAgeM(2, :, 331 : stepsPerYear : end))' squeeze(pct_hivAgeM(3, :, 331 : stepsPerYear : end))' squeeze(pct_hivAgeF(1, :, 331 : stepsPerYear : end))' ...
            squeeze(pct_hivAgeF(2, :, 331 : stepsPerYear : end))' squeeze(pct_hivAgeF(3, :, 331 : stepsPerYear : end))']);
        
        for b =  1 : age 
           colname(1, b)  = {[fileKey{z}, '_', ageGroup{b}, '_P25_hivAgeM']};
           colname(2, b) = {[fileKey{z}, '_', ageGroup{b}, '_Med_hivAgeM']};
           colname(3, b) = {[fileKey{z}, '_', ageGroup{b}, '_P75_hivAgeM']};
           colname(4, b)  = {[fileKey{z}, '_', ageGroup{b}, '_P25_hivAgeF']};
           colname(5, b) = {[fileKey{z}, '_', ageGroup{b}, '_Med_hivAgeF']};
           colname(6, b) = {[fileKey{z}, '_', ageGroup{b}, '_P75_hivAgeF']};
        end
        cellstr(colname);
        HIVtable2.Properties.VariableNames(1:size(HIVtable2, 2)) = {'Year', colname{1,:}, colname{2,:}, colname{3,:}, ...
           colname{4,:}, colname{5,:}, colname{6,:}};
        fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HIV.xlsx'];
        sheet = [fileKey{z}, '_HIVbyAge'];
        writetable(HIVtable2, fname, 'Sheet', sheet);
 end
 
 % HIV prevalence by risk among women over time 
 riskGroup = {'Low', 'Medium', 'High'};
 mean_hivRiskF = zeros(3, length(monthlyTimespan));
 pct_hivRiskF = zeros(3, 3, length(monthlyTimespan));
 colname = cell(3, length(risk));
 
 for z = 1 : length(fileKey)
    for r = 1 : risk
        mean_hivRiskF(r, :) = mean(squeeze(hivPrevRiskF(z, :, r, :)), 1, 'omitnan');
        pct_hivRiskF(:, r, :) = prctile(squeeze(hivPrevRiskF(z, :, r, :)), [25 50 75]);
    end
        HIVtable7 = array2table([annualTimespan(56:end)'  squeeze(pct_hivRiskF(1, :, 331 : stepsPerYear : end))' ...
            squeeze(pct_hivRiskF(2, :, 331 : stepsPerYear : end))' squeeze(pct_hivRiskF(3, :, 331 : stepsPerYear : end))']);
        
        for b =  1 : risk 
           colname(1, b)  = {[fileKey{z}, '_', riskGroup{b}, '_P25_hivRiskF']};
           colname(2, b) = {[fileKey{z}, '_', riskGroup{b}, '_Med_hivRiskF']};
           colname(3, b) = {[fileKey{z}, '_', riskGroup{b}, '_P75_hivRiskF']};
        end
        cellstr(colname);
        HIVtable7.Properties.VariableNames(1:size(HIVtable7, 2)) = {'Year', colname{1,:}, colname{2,:}, colname{3,:}};
        fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HIV.xlsx'];
        sheet = [fileKey{z}, '_HIVbyRisk'];
        writetable(HIVtable7, fname, 'Sheet', sheet);
 end
 
 % HIV incidence over time 
 for z = 1 : length(fileKey)
    mean_hivIncM = mean(squeeze(hivInc(z, :, 1, :)), 1, 'omitnan');
    pct_hivIncM = prctile(squeeze(hivInc(z, :, 1, :)), [25 50 75], 1);
    mean_hivIncF = mean(squeeze(hivInc(z, :, 2, :)), 1, 'omitnan');
    pct_hivIncF = prctile(squeeze(hivInc(z, :, 2, :)), [25 50 75], 1);
    
    HIVtable3 = array2table([annualTimespan(56:end)' pct_hivIncM(:, 56:end)' pct_hivIncF(:, 56:end)' ]);
    HIVtable3.Properties.VariableNames(1:size(HIVtable3, 2)) = {'Year', [fileKey{z},'_P25_hivIncM'], [fileKey{z}, '_Med_hivIncM'], ...
     [fileKey{z},'_P75_hivIncM'], [fileKey{z},'_P25_hivIncF'],  [fileKey{z},'_Med_hivIncF'], [fileKey{z},'_P75_hivIncF']}; 
     
    fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HIV.xlsx'];
    sheet = [fileKey{z}, '_HIVInc'];
    writetable(HIVtable3, fname, 'Sheet', sheet);
 end
 
 % HIV incidence by age over time 
 ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' ,...
    '60-64' , '65-69' , '70-74' , '75-79'};

 mean_hivIncAgeM = zeros(16, length(annualTimespan));
 pct_hivIncAgeM = zeros(3, 16, length(annualTimespan));
 mean_hivIncAgeF = mean_hivIncAgeM;
 pct_hivIncAgeF = pct_hivIncAgeM;
 colname = cell(6, length(age));
 for z = 1 : length(fileKey)
     for a = 1 : age
        mean_hivIncAgeM(a, :) = mean(squeeze(hivIncAgeM(z, :, a, :)), 1, 'omitnan');
        pct_hivIncAgeM(:, a, :) = prctile(squeeze(hivIncAgeM(z, :, a, :)), [25 50 75], 1);
        mean_hivIncAgeF(a, :) = mean(squeeze(hivIncAgeF(z, :, a, :)), 1, 'omitnan');
        pct_hivIncAgeF(:, a, :) = prctile(squeeze(hivIncAgeF(z, :, a, :)), [25 50 75], 1);
     end
     
     HIVtable4 = array2table([annualTimespan(56:end)'  squeeze(pct_hivIncAgeM(1, :, 56 : end))' squeeze(pct_hivIncAgeM(2, :, 56 : end))'  ...
            squeeze(pct_hivIncAgeM(3, :, 56 : end))' squeeze(pct_hivIncAgeF(1, :, 56 : end))' squeeze(pct_hivIncAgeF(2, :, 56 : end))' ...
            squeeze(pct_hivIncAgeF(3, :, 56 : end))']);
     for b =  1 : age 
           colname(1, b)  = {[fileKey{z}, '_', ageGroup{b}, '_P25_hivIncAgeM']};
           colname(2, b) = {[fileKey{z}, '_', ageGroup{b}, '_Med_hivIncAgeM']};
           colname(3, b) = {[fileKey{z}, '_', ageGroup{b}, '_P75_hivIncAgeM']};
           colname(4, b)  = {[fileKey{z}, '_', ageGroup{b}, '_P25_hivIncAgeF']};
           colname(5, b) = {[fileKey{z}, '_', ageGroup{b}, '_Med_hivIncAgeF']};
           colname(6, b) = {[fileKey{z}, '_', ageGroup{b}, '_P75_hivIncAgeF']};
     end
        cellstr(colname);
        HIVtable4.Properties.VariableNames(1:size(HIVtable4, 2)) = {'Year', colname{1,:}, colname{2,:}, colname{3,:}, ...
           colname{4,:}, colname{5,:}, colname{6,:}};
        fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HIV.xlsx'];
        sheet = [fileKey{z}, '_HIVIncbyAge'];
        writetable(HIVtable4, fname, 'Sheet', sheet);
 end
 
 % HIV incidence by risk among women over time 
 mean_hivIncRiskF = zeros(3, length(annualTimespan));
 pct_hivIncRiskF = zeros(3, 3, length(annualTimespan));
 colname = cell(3, length(risk));
 riskGroup = {'Low', 'Medium', 'High'};
 for z = 1 : length(fileKey)
     for r = 1 : risk
        mean_hivIncRiskF(r, :) = mean(squeeze(hivIncRiskF(z, :, r, :)), 1, 'omitnan');
        pct_hivIncRiskF(:, r, :) = prctile(squeeze(hivIncRiskF(z, :, r, :)), [25 50 75], 1);        
     end
     HIVtable6 = array2table([annualTimespan(56:end)' squeeze(pct_hivIncRiskF(1, :, 56 : end))' squeeze(pct_hivIncRiskF(2, :, 56 : end))' ...
            squeeze(pct_hivIncRiskF(3, :, 56 : end))']);
      for b =  1 : risk 
           colname(1, b)  = {[fileKey{z}, '_', riskGroup{b}, '_P25_hivIncRiskF']};
           colname(2, b) = {[fileKey{z}, '_', riskGroup{b}, '_Med_hivIncRiskF']};
           colname(3, b) = {[fileKey{z}, '_', riskGroup{b}, '_P75_hivIncRiskF']};
     end
        cellstr(colname);
        HIVtable6.Properties.VariableNames(1:size(HIVtable6, 2)) = {'Year', colname{1,:}, colname{2,:}, colname{3,:}};
        fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HIV.xlsx'];
        sheet = [fileKey{z}, '_HIVIncbyRisk'];
        writetable(HIVtable6, fname, 'Sheet', sheet);
 end
 % HIV deaths over time 
 for z = 1 : length(fileKey)
     mean_hivDeathsM = mean(squeeze(hivDeathsM(z, :, :)), 1, 'omitnan');
     pct_hivDeathsM = prctile(squeeze(hivDeathsM(z, :, :)), [25 50 75], 1); 
     mean_hivDeathsF = mean(squeeze(hivDeathsF(z, :, :)), 1, 'omitnan');
     pct_hivDeathsF = prctile(squeeze(hivDeathsF(z, :, :)), [25 50 75], 1); 
     
     HIVtable5 = array2table([annualTimespan(56:end)' pct_hivDeathsM(:, 56:end)' pct_hivDeathsF(:, 56:end)' ]);
     HIVtable5.Properties.VariableNames(1:size(HIVtable5,2)) = {'Year', [fileKey{z},'_P25_hivDeathsM'], [fileKey{z},'_Med_hivDeathsM'], ...
     [fileKey{z},'_P75_hivDeathsM'], [fileKey{z},'_P25_hivDeathsF'], [fileKey{z},'_Med_hivDeathsF'], [fileKey{z},'_P75_hivDeathsF']};
     
     fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HIV.xlsx'];
     sheet = [fileKey{z}, '_HIVdeaths'];
     writetable(HIVtable5, fname, 'Sheet', sheet);
 end
 
 %% HPV prevalence by age and HIV status in 2005, 2009, 2020, and 2070
 mean_hpv_hiv = zeros(length(age), length(hpvYearVec));
 pct_hpv_hiv = zeros(3, length(age), length(hpvYearVec));
 mean_hpv_hivNeg = mean_hpv_hiv;
 pct_hpv_hivNeg = pct_hpv_hiv;
 ageGroup = {'20-24' , '25-29' , '30-39' , '40-49' };
 colname=cell(6, length(ageGroup));

 for z = 1 : length(fileKey)
     for a = 1 : ageVecLength_fHPV
    %mean_hpv_hiv(a, :) = mean(squeeze(hpv_hiv(z, :, a, :)), 1, 'omitnan');
    pct_hpv_hiv(:, a, :) = prctile(squeeze(hpv_hiv(z, :, a, :)), [25 50 75], 1);
    %mean_hpv_hivNeg(a, :) = mean(squeeze(hpv_hivNeg(z, :, a, :)), 1, 'omitnan');
    pct_hpv_hivNeg(:, a, :) = prctile(squeeze(hpv_hivNeg(z, :, a, :)), [25 50 75], 1);      
     end
     
     HPVtable1 = array2table([hpvYearVec' squeeze(pct_hpv_hiv(1, :, :))' squeeze(pct_hpv_hiv(2, :, :))' squeeze(pct_hpv_hiv(3, :, :))' ...
         squeeze(pct_hpv_hivNeg(1, :, :))' squeeze(pct_hpv_hivNeg(2, :, :))' squeeze(pct_hpv_hivNeg(3, :, :))']);
     
        for b =  1 : ageVecLength_fHPV
           colname(1, b)  = {[fileKey{z}, '_', ageGroup{b}, '_P25_hpvPrevPos']};
           colname(2, b) = {[fileKey{z}, '_', ageGroup{b}, '_Med_hpvPrevPos']};
           colname(3, b) = {[fileKey{z}, '_', ageGroup{b}, '_P75_hpvPrevPos']};
           colname(4, b)  = {[fileKey{z}, '_', ageGroup{b}, '_P25_hpvPrevNeg']};
           colname(5, b) = {[fileKey{z}, '_', ageGroup{b}, '_Med_hpvPrevNeg']};
           colname(6, b) = {[fileKey{z}, '_', ageGroup{b}, '_P75_hpvPrevNeg']};
        end
        cellstr(colname);
        HPVtable1.Properties.VariableNames(1:size(HPVtable1,2)) = {'Year', colname{1,:}, colname{2,:}, colname{3,:}, ...
           colname{4,:}, colname{5,:}, colname{6,:}};
        fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HPV.xlsx'];
        sheet = [fileKey{z}, '_HPVbyAge'];
        writetable(HPVtable1, fname, 'Sheet', sheet);
 end 
 
 
% HPV prevalence over time 
pct_hpvHivTimeF = zeros(3, length(monthlyTimespan));
pct_hpvHivNegTimeF = pct_hpvHivTimeF;
pct_hpvAllTimeF = pct_hpvHivTimeF;
pct_hpvAllTimeM = pct_hpvHivTimeF;
for z = 1 : length(fileKey)
   mean_hpvHivTimeF = mean(squeeze(hpv_hivTimeF(z, :, :)), 1, 'omitnan');
   pct_hpvHivTimeF = prctile(squeeze(hpv_hivTimeF(z, :, :)), [25 50 75], 1);
   mean_hpvHivNegTimeF = mean(squeeze(hpv_hivNegTimeF(z, :, :)), 1, 'omitnan');
   pct_hpvHivNegTimeF = prctile(squeeze(hpv_hivNegTimeF(z, :, :)), [25 50 75], 1);
   mean_hpvAllTimeF = mean(squeeze(hpv_time(z, :, 2, :)), 1, 'omitnan');
   pct_hpvAllTimeF = prctile(squeeze(hpv_time(z, :, 2, :)), [25 50 75], 1);
   mean_hpvAllTimeM = mean(squeeze(hpv_time(z, :, 1, :)), 1, 'omitnan');
   pct_hpvAllTimeM = prctile(squeeze(hpv_time(z, :, 1, :)), [25 50 75], 1);
   
   HPVtable2 = array2table([annualTimespan(56:end)' squeeze(pct_hpvHivTimeF(:, 331:stepsPerYear:end))' ...
       squeeze(pct_hpvHivNegTimeF(:, 331:stepsPerYear:end))' squeeze(pct_hpvAllTimeF(:, 331:stepsPerYear:end))' ...
       squeeze(pct_hpvAllTimeM(:, 331:stepsPerYear:end))' ]);
   HPVtable2.Properties.VariableNames(1:size(HPVtable2, 2)) = {'Year', [fileKey{z}, '_P25_hpvHivTimeF'], [fileKey{z}, '_Med_hpvHivTimeF'], ...
       [fileKey{z}, '_P75_hpvHivTimeF'], [fileKey{z}, '_P25_hpvHivNegTimeF'], [fileKey{z}, '_Med_hpvHivNegTimeF'], [fileKey{z}, '_P75_hpvHivNegTimeF'],...
       [fileKey{z}, '_P25_hpvAllTimeF'], [fileKey{z}, '_Med_hpvAllTimeF'], [fileKey{z}, '_P75_hpvAllTimeF'], ...
       [fileKey{z}, '_P25_hpvAllTimeM'], [fileKey{z}, '_Med_hpvAllTimeM'], [fileKey{z}, '_P75_hpvAllTimeM']};
   fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HPV.xlsx'];
        sheet = [fileKey{z}, '_HPVPrevbyHIV'];
        writetable(HPVtable2, fname, 'Sheet', sheet);
end

% HPV incidence by HIV status among women  
% for z = 1 : length(fileKey)
%    mean_hpvIncAll = mean(squeeze(hpvIncTime(z, :, :)), 1, 'omitnan');
%    pct_hpvIncAll = prctile(squeeze(hpvIncTime(z, :, :)), [25 50 75], 1);
%    mean_hpvIncNeg = mean(squeeze(hpvIncTimeNeg(z, :, :)), 1, 'omitnan');
%    pct_hpvIncNeg = prctile(squeeze(hpvIncTimeNeg(z, :, :)), [25 50 75], 1);
%    mean_hpvIncNoArt = mean(squeeze(hpvIncTimePos(z, :, :)), 1, 'omitnan');
%    pct_hpvIncNoArt = prctile(squeeze(hpvIncTimePos(z, :, :)), [25 50 75], 1);
%    mean_hpvIncPosAll = mean(squeeze(hpvIncTimePosAll(z, :, :)), 1, 'omitnan');
%    pct_hpvIncPosAll = prctile(squeeze(hpvIncTimePosAll(z, :, :)), [25 50 75], 1);
%    mean_hpvIncArt = mean(squeeze(hpvIncTimeArt(z, :, :)), 1, 'omitnan');
%    pct_hpvIncArt = prctile(squeeze(hpvIncTimeArt(z, :, :)), [25 50 75], 1);
%    
%    HPVtable3 = array2table([annualTimespan(56:end)' squeeze(pct_hpvIncAll(:,56:end))' squeeze(pct_hpvIncNeg(:,56:end))' ...
%        squeeze(pct_hpvIncNoArt(:,56:end))' squeeze(pct_hpvIncPosAll(:,56:end))' squeeze(pct_hpvIncArt(:,56:end))' ]);
%    
%    HPVtable3.Properties.VariableNames(1:size(HPVtable3, 2)) = {'Year', ...
%        [fileKey{z}, '_P25_hpvIncAll'] , [fileKey{z}, '_Med_hpvIncAll'], [fileKey{z}, '_P75_hpvIncAll'], ...
%        [fileKey{z}, '_P25_hpvIncNeg'] , [fileKey{z}, '_Med_hpvIncNeg'], [fileKey{z}, '_P75_hpvIncNeg'], ... 
%        [fileKey{z}, '_P25_hpvIncNoArt'] , [fileKey{z}, '_Med_hpvIncNoArt'], [fileKey{z}, '_P75_hpvIncNoArt'], ... 
%        [fileKey{z}, '_P25_hpvIncPosAll'] , [fileKey{z}, '_Med_hpvIncPosAll'], [fileKey{z}, '_P75_hpvIncPosAll'], ... 
%        [fileKey{z}, '_P25_hpvIncArt'] , [fileKey{z}, '_Med_hpvIncArt'], [fileKey{z}, '_P75_hpvIncArt'] };
%     
%     fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HPV.xlsx'];
%         sheet = [fileKey{z}, '_HPVInc'];
%         writetable(HPVtable3, fname, 'Sheet', sheet);
%    
% end

%% CIN prevalence by HIV status 
cinlvl = ["CIN1"; "CIN2"; "CIN3"];
for z = 1 : length(fileKey)
   mean_cinPos = mean(squeeze(cinPosAge(z, :, :)), 1, 'omitnan' ); 
   pct_cinPos = prctile(squeeze(cinPosAge(z, :, :)), [25 50 75], 1);
   mean_cinNeg = mean(squeeze(cinNegAge(z, :, :)), 1, 'omitnan');
   pct_cinNeg = prctile(squeeze(cinNegAge(z, :, :)), [25 50 75], 1); 
  
   CINtable = array2table([cinlvl squeeze(pct_cinPos)' squeeze(pct_cinNeg)']);
   CINtable.Properties.VariableNames(1:size(CINtable, 2)) = {'CIN_level', ...
       [fileKey{z}, '_P25_cinPos'] , [fileKey{z}, '_Med_cinPos'], [fileKey{z}, '_P75_cinPos'], ...
       [fileKey{z}, '_P25_cinNeg'] ,  [fileKey{z}, '_Med_cinNeg'], [fileKey{z}, '_P75_cinNeg'] };
   fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HPV.xlsx'];
        sheet = [fileKey{z}, '_CINprev'];
        writetable(CINtable, fname, 'Sheet', sheet);
end

% HPV prevalence by HIV status among high risk women 
ageGroup = ["20-24"; "25-29"; "30-39"; "40-49"];
for z = 1 : length(fileKey)
        mean_hpv_hiv_risk = mean(squeeze(hpv_hiv_risk(z, :, :)), 1, 'omitnan' );
        pct_hpv_hiv_risk = prctile(squeeze(hpv_hiv_risk(z, :, :)), [25 50 75], 1);
        mean_hpv_hivNeg_risk = mean(squeeze(hpv_hivNeg_risk(z, :, :)), 1, 'omitnan' );
        pct_hpv_hivNeg_risk = prctile(squeeze(hpv_hivNeg_risk(z, :, :)), [25 50 75], 1);
    
    HPVtable4 = array2table([ageGroup squeeze(pct_hpv_hiv_risk)' squeeze(pct_hpv_hivNeg_risk)']);
    HPVtable4.Properties.VariableNames(1:size(HPVtable4, 2))= {'Age', [fileKey{z}, '_P25_hpv_hiv_risk'] , [fileKey{z}, '_Med_hpv_hiv_risk'], ...
        [fileKey{z}, '_P75_hpv_hiv_risk'], [fileKey{z}, '_P25_hpv_hivNeg_risk'], [fileKey{z}, '_Med_hpv_hivNeg_risk'], [fileKey{z}, '_P75_hpv_hivNeg_risk']};
    fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HPV.xlsx'];
        sheet = [fileKey{z}, '_HPVprevRisk'];
        writetable(HPVtable4, fname, 'Sheet', sheet);
end
% Type distribution 
yrs = ["2005"; "2020"; "2070"];
stat = ["P25"; "Med"; "P75"];
for z = 1 : length(fileKey)
        pct_hpv_vax = prctile(squeeze(hpv_vax(z, :, 1, :)), [25 50 75], 1 );
        pct_hpv_nonVax = prctile(squeeze(hpv_nonVax(z, :, 1, :)), [25 50 75], 1);
        pct_cin1_vax = prctile(squeeze(cin1_vax(z, :, 1, :)), [25 50 75], 1 );
        pct_cin1_nonVax = prctile(squeeze(cin1_nonVax(z, :, 1, :)), [25 50 75], 1);
        pct_cin2_vax = prctile(squeeze(cin2_vax(z, :, 1, :)), [25 50 75], 1 );
        pct_cin2_nonVax = prctile(squeeze(cin2_nonVax(z, :, 1, :)), [25 50 75], 1);
        pct_cin3_vax = prctile(squeeze(cin3_vax(z, :, 1, :)), [25 50 75], 1 );
        pct_cin3_nonVax = prctile(squeeze(cin3_nonVax(z, :, 1, :)), [25 50 75], 1);
        pct_cc_vax = prctile(squeeze(cc_vax(z, :, 1, :)), [25 50 75], 1 );
        pct_cc_nonVax = prctile(squeeze(cc_nonVax(z, :, 1, :)), [25 50 75], 1);
    
    HPVtable5 = array2table([repmat(yrs(1), 30,1) repmat(stat, 10, 1) [squeeze(pct_hpv_vax(: , 481)); ...
        squeeze(pct_hpv_nonVax(: , 481)); squeeze(pct_cin1_vax(: , 481)); ...
        squeeze(pct_cin1_nonVax(: , 481)); squeeze(pct_cin2_vax(: , 481)); squeeze(pct_cin2_nonVax(: , 481)); ...
        squeeze(pct_cin3_vax(: , 481)); squeeze(pct_cin3_nonVax(: , 481)); squeeze(pct_cc_vax(: , 481)); ...
        squeeze(pct_cc_nonVax(: , 481)) ]]);
    HPVtable5.Properties.VariableNames(1:size(HPVtable5, 2))= {'Year', 'Stat', [fileKey{z}, '_Estimate']};
    fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_HPV.xlsx'];
        sheet = [fileKey{z}, '_typeDist'];
        writetable(HPVtable5, fname, 'Sheet', sheet);
end

%% ART coverage among men and women by year 
pct_artcovM = zeros(length(fileKey), length(monthlyTimespan));
pct_artcovF = pct_artcovM;
for z = 1 : length(fileKey)
   pct_artcovM = prctile(squeeze(artCovM(z, :, :)), 50);
   pct_artcovF = prctile(squeeze(artCovF(z, :, :)), 50);
   
   ARTtable = array2table([annualTimespan(56:end)' pct_artcovM(:, 331:stepsPerYear:end)' pct_artcovF(:, 331:stepsPerYear:end)']); 
   ARTtable.Properties.VariableNames(1:size(ARTtable, 2)) = {'Year', [fileKey{z}, '_artCov_Male'], [fileKey{z}, '_artCov_Female']};
   fname = [pwd, '\HHCoM_Results\Results_9v\StochMod_output_ART.xlsx'];
   sheet = [fileKey{z}, '_artCov'];
   writetable(ARTtable, fname, 'Sheet', sheet);
end



