function [] = vaxCEA_multSims_CIs_IPVC(vaxResultInd , sceNum , fileNameNums)
% example: vaxCEA_multSims_CIs_IPVC(1 , '34' , {'3' , '4' , '0'})

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

lastYear = 2122;

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
monthlyTimespan = [startYear : timeStep : lastYear];
monthlyTimespan = monthlyTimespan(1 : end-1);
annualTimespan = [startYear : lastYear-1];
futAnnualTimespan = [2019 : lastYear-1];
% Total population size
popSizeAgeF = zeros(nRuns , 5 , age , length(monthlyTimespan));
popSizeAgeC = popSizeAgeF;
% Female HPV prevalence
hpvYearVec_orig = [2002 2018];
hpvYearVec2018_orig = [2018];
hpvHivAgeF = zeros(nRuns , 5 , age , length(monthlyTimespan));
hpv9VHivAgeF = hpvHivAgeF;
hpv9VHivAgeC = hpvHivAgeF;
hpvAllHivAgeC = hpvHivAgeF;
% CC incidence
ccIncHivAgeTime = zeros(nRuns , 5 , age , length(annualTimespan));
ccCumHivAgeTime = zeros(nRuns , 5 , age , length(futAnnualTimespan));
ccCumHivTime = zeros(nRuns , 5 , length([2021 : 2121]));
diseaseVec_ccInc = {[1 : disease] , [1 : 2] , [3 : 8] , [3 : 7] , 8};
diseaseVecLength_ccInc = length(diseaseVec_ccInc);
% HPV vaccination and screening
cumVaxSchool = zeros(nRuns , length([2021 : 2121]));
cumVaxCU = zeros(nRuns , length([2021 : 2121]));

resultsDir = [pwd , '\HHCoM_Results\'];
fileKey = {'sim1' , 'sim2' , 'sim0'};
fileKeyNums = fileNameNums;
n = vaxResultInd;
baseFileName = ['Vaccine22Apr20Ph2V11_baseVax057_baseScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_8ts-2021_WHO-SCES' , sceNum , '_057_'];
loopSegments = {0 , round(nRuns/2) , nRuns};
loopSegmentsLength = length(loopSegments);
for k = 1 : loopSegmentsLength-1
    parfor j = loopSegments{k}+1 : loopSegments{k+1}
        % Load results
        pathModifier = [baseFileName , fileInds{j}]; % ***SET ME***: name for simulation output file
        nSims = size(dir([pwd , '\HHCoM_Results\' , pathModifier, '\' , '*.mat']) , 1);
        curr = load([pwd , '/HHCoM_Results/toNow_22Apr20Ph2V11_baseVax057_baseScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_8ts-2021_' , fileInds{j}]); % ***SET ME***: name for historical run output file 

        vaxResult = cell(nSims , 1);
        resultFileName = [pwd , '\HHCoM_Results\' , pathModifier, '\' , 'vaxSimResult'];
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
        vaxResult{n}.newScreen = [vaxResult{n}.newScreen(1 : end , : , : , : , : , : , :)]; %[curr.newScreen(1 : end , : , : , : , : , : , : , : , :); vaxResult{n}.newScreen(2 : end , : , : , : , : , : , : , : , :)];
        %vaxResult{n}.newTreatImm = [vaxResult{n}.newTreatImm(1 : end , : , : , : , : , : , : , : , :)];
        %vaxResult{n}.newTreatHpv = [vaxResult{n}.newTreatHpv(1 : end , : , : , : , : , : , : , : , :)];
        vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , : , : , : , :); vaxResult{n}.newHiv(2 : end , : , : , : , : , : , :)];
        vaxResult{n}.hivDeaths = [curr.hivDeaths(1 : end , : , : , :); vaxResult{n}.hivDeaths(2 : end , : , : , :)];
        vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :); vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , :)];
        vaxResult{n}.vaxdSchool = vaxResult{n}.vaxdSchool(1 : end , :);
        vaxResult{n}.vaxdCU = vaxResult{n}.vaxdCU(1 : end , :);
        vaxResult{n}.tVec = [curr.tVec(1 : end), vaxResult{n}.tVec(2 : end)];

    %     noVaxInd = nSims;
    %     noV = vaxResult{noVaxInd};
        tVec = vaxResult{n}.tVec;
        tVecYr = tVec(1 : stepsPerYear : end);
        
        % Initialize variables
        hpvYearVec = hpvYearVec_orig;
        hpvYearVec2018 = hpvYearVec2018_orig;

        %% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************
        
        %% Female total population size by 5-year age groups over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                popAgeF = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
                popSizeAgeF(j , dInd , a , :) = sum(vaxResult{n}.popVec(: , popAgeF),2);
            end
        end
        
        %% Male+Female total population size by 5-year age groups over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                popAgeF = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                    1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
                popSizeAgeC(j , dInd , a , :) = sum(vaxResult{n}.popVec(: , popAgeF),2);
            end
        end
        
        
        %% ********************************** HPV FIGURES **********************************************************************************************
    
        %% Female hrHPV prevalence by 5-year age groups over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                hpvInds = unique([toInd(allcomb(d , 1 : viral , 2 : 6 , 1 : hpvNonVaxStates , ...
                1 : 3 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
                1 : hpvVaxStates , 2 : 6 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk))]);
                hpvHivAgeF(j , dInd , a , :) = sum(vaxResult{n}.popVec(: , hpvInds),2);
            end
        end
        
        %% Female VT-hrHPV prevalence by 5-year age groups over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                hpvInds = toInd(allcomb(d , 1 : viral , 2 : 6 , 1 : hpvNonVaxStates , ...
                1 : 3 , 1 : intervens , 2 , a , 1 : risk));
                hpv9VHivAgeF(j , dInd , a , :) = sum(vaxResult{n}.popVec(: , hpvInds),2);
            end
        end
        
        %% Male+Female Combined VT-hrHPV prevalence by 5-year age groups over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                hpvInds = toInd(allcomb(d , 1 : viral , 2 : 6 , 1 : hpvNonVaxStates , ...
                1 : 3 , 1 : intervens , 1 : gender , a , 1 : risk));
                hpv9VHivAgeC(j , dInd , a , :) = sum(vaxResult{n}.popVec(: , hpvInds),2);
            end
        end
        
        %% Male+Female Combined all-hrHPV prevalence by 5-year age groups over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                hpvInds = unique([toInd(allcomb(d , 1 : viral , 2 : 6 , 1 : hpvNonVaxStates , ...
                    1 : 3 , 1 : intervens , 1 : gender , a , 1 : risk));
                    toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 2 : 6 , ...
                    1 : 3 , 1 : intervens , 1 : gender , a , 1 : risk))]);
                hpvAllHivAgeC(j , dInd , a , :) = sum(vaxResult{n}.popVec(: , hpvInds),2);
            end
        end
        
        
        %% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************
        
        %% Cervical cancer incidence by HIV status and age over time
        fac = 10 ^ 5;
        for dInd = 1 : diseaseVecLength_ccInc
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
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            for a = 1 : age
                ccCumHivAgeTime(j , dInd , a , :) = ...
                    cumsum(squeeze(annlz(sum(sum(vaxResult{n}.newCC(((2019 - startYear) * stepsPerYear +1):end , d , a , :),2),4))),2);
            end
        end
        

        %% Cumulative cervical cancer cases by HIV status over time
        for dInd = 1 : diseaseVecLength_ccInc
            d = diseaseVec_ccInc{dInd};
            ccCumHivTime(j , dInd , :) = ...
                cumsum(squeeze(annlz(sum(sum(sum(vaxResult{n}.newCC(((2021 - startYear) * stepsPerYear +1):end , d , : , :),2),3),4))),2);
        end
        
       
        %% ************************** SCREENING & VACCINATION FIGURES *******************************************************************************
        
        %% Total number school-based vaccines administered
        cumVaxSchool(j , :) = annlz(vaxResult{n}.vaxdSchool(:,1)) * 2;

        %% Total number CU vaccines administered
        cumVaxCU(j , :) = annlz(vaxResult{n}.vaxdCU(:,1)) * 3;
        
    end
end

%% ***************************** DEMOGRAPHY FIGURES *********************************************************************************************

%% Write female total population size by 5-year age groups over time (2019-2120) into existing template
diseaseLabels = {'All - Pop (P-Y)' , 'HIV- (P-Y)' , 'HIV+ (P-Y)' , 'HIV+ no ART (P-Y)' , 'HIV+ ART (P-Y)'};
firstYrInd = ((2019 - startYear) * stepsPerYear +((stepsPerYear/2)+1));
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_CC_IncidenceRates-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(popSizeAgeF(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end

diseaseLabels = {'All - Pop (P-Y)' , 'HIV- (P-Y)' , 'HIV+ (P-Y)' , 'HIV+ no ART (P-Y)' , 'HIV+ ART (P-Y)'};
firstYrInd = ((2020 - startYear) * stepsPerYear +((stepsPerYear/2)+1));
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_CIN2+_Prevalence-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(popSizeAgeF(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(popSizeAgeF(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(popSizeAgeF(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_VThrHPV-MWcombined_Prevalence-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(popSizeAgeC(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(prctile(squeeze(popSizeAgeF(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 5 , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(prctile(squeeze(popSizeAgeF(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 95 , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end

diseaseLabels = {'All - Pop (N)' , 'HIV- (N)' , 'HIV+ (N)' , 'HIV+ no ART (N)' , 'HIV+ ART (N)'};
firstYrInd = ((2019 - startYear) * stepsPerYear +4);
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_CumulativeImpact_CC-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(popSizeAgeF(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end


%% ********************************** HPV FIGURES **********************************************************************************************

%% Female hrHPV prevalence by 5-year age groups over time               
diseaseLabels = {'Pop(All) (hrHPV)' , 'HIV-  (hrHPV)' , 'HIV+   (hrHPV)' , 'HIV+ no ART  (hrHPV)' , 'HIV+ ART  (hrHPV)'};
firstYrInd = ((2020 - startYear) * stepsPerYear +((stepsPerYear/2)+1));
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(hpvHivAgeF(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end

%% Write age-standardized hrHPV cases by HIV status over time
% Note: the age-standardization process shifts the prevalence of the
% last modelled age group to the next age group in the following year.
% However, prevalence is NaN prior to HIV introduction in the
% HIV-positive no ART group, and NaN prior to ART introduction in the
% HIV-positive ART group. Since we have four age groups past the 16 we
% model, a NaN value is present for four years past the introduction of
% HIV/ART, leading to a NaN value for summed HPV infected during these 
% years. We therefore lack data in this four-year interval in the
% saved/plotted results.
worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
    247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
    23477 9261 2155];

diseaseLabels = {'Pop(All) (hrHPV)' , 'HIV-  (hrHPV)' , 'HIV+   (hrHPV)' , 'HIV+ no ART  (hrHPV)' , 'HIV+ ART  (hrHPV)'};
firstYrRange2 = (lastYear-1) - 2020;
for dInd = 1 : length(diseaseLabels)
    hpvHivAgeF_dis = squeeze(hpvHivAgeF(: , dInd , : , (((stepsPerYear/2)+1):stepsPerYear:end))) ./ squeeze(popSizeAgeF(: , dInd , : , (((stepsPerYear/2)+1):stepsPerYear:end)));

    numHpvTot = zeros(size(hpvHivAgeF_dis,1) , 1 , size(hpvHivAgeF_dis,3));       
    for aInd = 1:age+4
        a = aInd;
        if aInd >= age
            a = age;
        end

        if aInd <= age    
            numHpv =hpvHivAgeF_dis(: , a , :) .* worldStandard_WP2015(aInd);
            if (a < 3)
                numHpv = zeros(size(hpvHivAgeF_dis,1) , 1 , size(hpvHivAgeF_dis,3));
            end
        elseif aInd > age
            numHpv = hpvHivAgeF_dis(: , a , :);
            numHpv = cat(3 , (ones(size(numHpv,1),1,aInd-a).*numHpv(:,1,1)) , numHpv(: , 1 ,1:end-(aInd-a)));
            numHpv = numHpv .* worldStandard_WP2015(aInd);
        end
        numHpvTot = numHpvTot + numHpv;
    end
    hpvCumTot = numHpvTot;
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_hrHPV_Prevalence-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(hpvCumTot(: , : , (end-firstYrRange2):end)) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B4')   
end
   
%% Female VT-hrHPV prevalence by 5-year age groups over time
diseaseLabels = {'Pop(All) (vtHPV)' , 'HIV-  (vtHPV)' , 'HIV+   (vtHPV)' , 'HIV+ no ART  (vtHPV)' , 'HIV+ ART  (vtHPV)'};
firstYrInd = ((2020 - startYear) * stepsPerYear +((stepsPerYear/2)+1));   
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(hpv9VHivAgeF(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end

%% Male+Female Combined VT-hrHPV prevalence by 5-year age groups over time
diseaseLabels = {'Pop(All) (vtHPV)' , 'HIV-  (vtHPV)' , 'HIV+   (vtHPV)' , 'HIV+ no ART  (vtHPV)' , 'HIV+ ART  (vtHPV)'};
firstYrInd = ((2020 - startYear) * stepsPerYear +((stepsPerYear/2)+1));   
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_VThrHPV-MWcombined_Prevalence-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(hpv9VHivAgeC(: , dInd , : , (firstYrInd : stepsPerYear : end))) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B6')   
end

%% Write age-standardized VT-hrHPV cases by HIV status over time
% Note: the age-standardization process shifts the prevalence of the
% last modelled age group to the next age group in the following year.
% However, prevalence is NaN prior to HIV introduction in the
% HIV-positive no ART group, and NaN prior to ART introduction in the
% HIV-positive ART group. Since we have four age groups past the 16 we
% model, a NaN value is present for four years past the introduction of
% HIV/ART, leading to a NaN value for summed HPV infected during these 
% years. We therefore lack data in this four-year interval in the
% saved/plotted results.
worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
    247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
    23477 9261 2155];

diseaseLabels = {'Pop(All) (vtHPV)' , 'HIV-  (vtHPV)' , 'HIV+   (vtHPV)' , 'HIV+ no ART  (vtHPV)' , 'HIV+ ART  (vtHPV)'};
firstYrRange2 = (lastYear-1) - 2020;
for dInd = 1 : length(diseaseLabels)
    hpvHivAgeF_dis = squeeze(hpv9VHivAgeF(: , dInd , : , (((stepsPerYear/2)+1):stepsPerYear:end))) ./ squeeze(popSizeAgeF(: , dInd , : , (((stepsPerYear/2)+1):stepsPerYear:end)));

    numHpvTot = zeros(size(hpvHivAgeF_dis,1) , 1 , size(hpvHivAgeF_dis,3));       
    for aInd = 1:age+4
        a = aInd;
        if aInd >= age
            a = age;
        end

        if aInd <= age    
            numHpv =hpvHivAgeF_dis(: , a , :) .* worldStandard_WP2015(aInd);
            if (a < 3)
                numHpv = zeros(size(hpvHivAgeF_dis,1) , 1 , size(hpvHivAgeF_dis,3));
            end
        elseif aInd > age
            numHpv = hpvHivAgeF_dis(: , a , :);
            numHpv = cat(3 , (ones(size(numHpv,1),1,aInd-a).*numHpv(:,1,1)) , numHpv(: , 1 ,1:end-(aInd-a)));
            numHpv = numHpv .* worldStandard_WP2015(aInd);
        end
        numHpvTot = numHpvTot + numHpv;
    end
    hpvCumTot = numHpvTot;
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_VThrHPV_Prevalence-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(hpvCumTot(: , : , (end-firstYrRange2):end)) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B4')   
end

%% Write age-standardized male+female combined VT-hrHPV cases by HIV status over time
% Note: the age-standardization process shifts the prevalence of the
% last modelled age group to the next age group in the following year.
% However, prevalence is NaN prior to HIV introduction in the
% HIV-positive no ART group, and NaN prior to ART introduction in the
% HIV-positive ART group. Since we have four age groups past the 16 we
% model, a NaN value is present for four years past the introduction of
% HIV/ART, leading to a NaN value for summed HPV infected during these 
% years. We therefore lack data in this four-year interval in the
% saved/plotted results.
worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
    247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
    23477 9261 2155];

diseaseLabels = {'Pop(All) (vtHPV)' , 'HIV-  (vtHPV)' , 'HIV+   (vtHPV)' , 'HIV+ no ART  (vtHPV)' , 'HIV+ ART  (vtHPV)'};
firstYrRange2 = (lastYear-1) - 2020;
for dInd = 1 : length(diseaseLabels)
    hpvHivAgeF_dis = squeeze(hpv9VHivAgeC(: , dInd , : , (((stepsPerYear/2)+1):stepsPerYear:end))) ./ squeeze(popSizeAgeC(: , dInd , : , (((stepsPerYear/2)+1):stepsPerYear:end)));

    numHpvTot = zeros(size(hpvHivAgeF_dis,1) , 1 , size(hpvHivAgeF_dis,3));       
    for aInd = 1:age+4
        a = aInd;
        if aInd >= age
            a = age;
        end

        if aInd <= age    
            numHpv =hpvHivAgeF_dis(: , a , :) .* worldStandard_WP2015(aInd);
            if (a < 3)
                numHpv = zeros(size(hpvHivAgeF_dis,1) , 1 , size(hpvHivAgeF_dis,3));
            end
        elseif aInd > age
            numHpv = hpvHivAgeF_dis(: , a , :);
            numHpv = cat(3 , (ones(size(numHpv,1),1,aInd-a).*numHpv(:,1,1)) , numHpv(: , 1 ,1:end-(aInd-a)));
            numHpv = numHpv .* worldStandard_WP2015(aInd);
        end
        numHpvTot = numHpvTot + numHpv;
    end
    hpvCumTot = numHpvTot;
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_VThrHPV-MWcombined_Prevalence-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(hpvCumTot(: , : , (end-firstYrRange2):end)) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B4')   
end

%% Write age-standardized male+female combined VT-hrHPV prevalence for ages 15+ by HIV status over time
% Note: the age-standardization process shifts the prevalence of the
% last modelled age group to the next age group in the following year.
% However, prevalence is NaN prior to HIV introduction in the
% HIV-positive no ART group, and NaN prior to ART introduction in the
% HIV-positive ART group. Since we have four age groups past the 16 we
% model, a NaN value is present for four years past the introduction of
% HIV/ART, leading to a NaN value for summed HPV infected during these 
% years. We therefore lack data in this four-year interval in the
% saved/plotted results.
worldStandard_WP2015 = [287187 291738 299655 272348 ... %325428 311262 295693 
    247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
    23477 9261 2155];

diseaseLabels = {'Pop(All) (vtHPV)' , 'HIV-  (vtHPV)' , 'HIV+   (vtHPV)' , 'HIV+ no ART  (vtHPV)' , 'HIV+ ART  (vtHPV)'};
firstYrRange2 = (lastYear-1) - 2020;
for dInd = 1 : length(diseaseLabels)
    hpvHivAgeF_dis = squeeze(hpv9VHivAgeC(: , dInd , 4 : age , (((stepsPerYear/2)+1):stepsPerYear:end))) ./ squeeze(popSizeAgeC(: , dInd , 4 : age , (((stepsPerYear/2)+1):stepsPerYear:end)));

    numHpvTot = zeros(size(hpvHivAgeF_dis,1) , 1 , size(hpvHivAgeF_dis,3));       
    for aInd = 1:age-3+4
        a = aInd;
        if aInd >= age-3
            a = age-3;
        end

        if aInd <= age-3    
            numHpv =hpvHivAgeF_dis(: , a , :) .* worldStandard_WP2015(aInd);
            if (a < 3)
                numHpv = zeros(size(hpvHivAgeF_dis,1) , 1 , size(hpvHivAgeF_dis,3));
            end
        elseif aInd > age-3
            numHpv = hpvHivAgeF_dis(: , a , :);
            numHpv = cat(3 , (ones(size(numHpv,1),1,aInd-a).*numHpv(:,1,1)) , numHpv(: , 1 ,1:end-(aInd-a)));
            numHpv = numHpv .* worldStandard_WP2015(aInd);
        end
        numHpvTot = numHpvTot + numHpv;
    end
    hpvCumTot = numHpvTot ./ sum(worldStandard_WP2015);
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'VThrHPV-MWcombined_Prevalence_wUncert_ages15plus-standardised-(2020-2121)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix([[2020 : lastYear-1] ; ...
        squeeze(median(squeeze(hpvCumTot(: , : , (end-firstYrRange2):end)) , 1)) ; ...
        squeeze(min(squeeze(hpvCumTot(: , : , (end-firstYrRange2):end)) , [] , 1)) ; ...
        squeeze(max(squeeze(hpvCumTot(: , : , (end-firstYrRange2):end)) , [] , 1))] , fname , 'Sheet' , diseaseLabels{dInd})   
end

%% Write age-standardized male+female combined all-hrHPV prevalence for ages 15+ by HIV status over time
% Note: the age-standardization process shifts the prevalence of the
% last modelled age group to the next age group in the following year.
% However, prevalence is NaN prior to HIV introduction in the
% HIV-positive no ART group, and NaN prior to ART introduction in the
% HIV-positive ART group. Since we have four age groups past the 16 we
% model, a NaN value is present for four years past the introduction of
% HIV/ART, leading to a NaN value for summed HPV infected during these 
% years. We therefore lack data in this four-year interval in the
% saved/plotted results.
worldStandard_WP2015 = [287187 291738 299655 272348 ... %325428 311262 295693 
    247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
    23477 9261 2155];

diseaseLabels = {'Pop(All) (vtHPV)' , 'HIV-  (vtHPV)' , 'HIV+   (vtHPV)' , 'HIV+ no ART  (vtHPV)' , 'HIV+ ART  (vtHPV)'};
firstYrRange2 = (lastYear-1) - 2020;
for dInd = 1 : length(diseaseLabels)
    hpvHivAgeF_dis = squeeze(hpvAllHivAgeC(: , dInd , 4 : age , (((stepsPerYear/2)+1):stepsPerYear:end))) ./ squeeze(popSizeAgeC(: , dInd , 4 : age , (((stepsPerYear/2)+1):stepsPerYear:end)));

    numHpvTot = zeros(size(hpvHivAgeF_dis,1) , 1 , size(hpvHivAgeF_dis,3));       
    for aInd = 1:age-3+4
        a = aInd;
        if aInd >= age-3
            a = age-3;
        end

        if aInd <= age-3    
            numHpv =hpvHivAgeF_dis(: , a , :) .* worldStandard_WP2015(aInd);
            if (a < 3)
                numHpv = zeros(size(hpvHivAgeF_dis,1) , 1 , size(hpvHivAgeF_dis,3));
            end
        elseif aInd > age-3
            numHpv = hpvHivAgeF_dis(: , a , :);
            numHpv = cat(3 , (ones(size(numHpv,1),1,aInd-a).*numHpv(:,1,1)) , numHpv(: , 1 ,1:end-(aInd-a)));
            numHpv = numHpv .* worldStandard_WP2015(aInd);
        end
        numHpvTot = numHpvTot + numHpv;
    end
    hpvCumTot = numHpvTot ./ sum(worldStandard_WP2015);
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'AllHrHPV-MWcombined_Prevalence_wUncert_ages15plus-standardised-(2020-2121)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix([[2020 : lastYear-1] ; ...
        squeeze(median(squeeze(hpvCumTot(: , : , (end-firstYrRange2):end)) , 1)) ; ...
        squeeze(min(squeeze(hpvCumTot(: , : , (end-firstYrRange2):end)) , [] , 1)) ; ...
        squeeze(max(squeeze(hpvCumTot(: , : , (end-firstYrRange2):end)) , [] , 1))] , fname , 'Sheet' , diseaseLabels{dInd})   
end


%% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

%% Write cervical cancer incidence rates by 5-year age groups over time (2019-2120) into existing template
diseaseLabels = {'Pop(All) ICC' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'};
firstYrInd = ((2019 - startYear) +1);
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_CC_IncidenceRates-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(ccIncHivAgeTime(: , dInd , 3:age , firstYrInd:end)) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B8')   
end

diseaseLabels = {'Pop(All) ICC' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'};
firstYrInd = ((2020 - startYear) +1);
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(prctile(squeeze(ccIncHivAgeTime(: , dInd , 3:age , firstYrInd:end)) , 5 , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B8')   
end
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(prctile(squeeze(ccIncHivAgeTime(: , dInd , 3:age , firstYrInd:end)) , 95 , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B8')   
end

if (vaxResultInd == 3) && contains(baseFileName , 'noBaseVax_baseScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_WHO')
    diseaseLabels = {'Pop(All) (ICC)' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'};
    firstYrInd = ((1990 - startYear) +1);
    lastYrInd = ((2020 - startYear) +1);
    for dInd = 1 : length(diseaseLabels)
        fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
            'UofW_Pre-Impact_CC_IncidenceRates-standardised-(Before_2020)_S' , fileKeyNums{n} , 'f.xlsx'];
        writematrix(squeeze(median(squeeze(ccIncHivAgeTime(: , dInd , 3:age , firstYrInd:lastYrInd)) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B8')   
    end  
end

%% Write cervical cancer age-standardized incidence rates over time (2019-2120) into existing template
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

diseaseLabels = {'Pop(All) ICC' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'};
diseaseLabels2 = {'Pop(All) (ICC)' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'};
firstYrInd = ((2019 - startYear) +1);
firstYrInd2 = ((1990 - startYear) +1);
firstYrInd3 = ((2020 - startYear) +1);
lastYrInd = ((2020 - startYear) +1);
for dInd = 1 : length(diseaseLabels)
    ccIncHivAgeTime_med = squeeze(median(squeeze(ccIncHivAgeTime(: , dInd , : , :)) , 1));
    ccIncHivAgeTime_lb = squeeze(prctile(squeeze(ccIncHivAgeTime(: , dInd , : , :)) , 5 , 1));
    ccIncHivAgeTime_ub = squeeze(prctile(squeeze(ccIncHivAgeTime(: , dInd , : , :)) , 95 , 1));

    ccIncRefTot = zeros(1 , size(ccIncHivAgeTime_med,2)); 
    ccIncRefTot_lb = zeros(1 , size(ccIncHivAgeTime_lb,2)); 
    ccIncRefTot_ub = zeros(1 , size(ccIncHivAgeTime_ub,2)); 
    for aInd = 1:age+4
        a = aInd;
        if aInd >= age
            a = age;
        end

        if aInd <= age    
            ccIncRef = ccIncHivAgeTime_med(a , :) .* worldStandard_WP2015(aInd);
            ccIncRef_lb = ccIncHivAgeTime_lb(a , :) .* worldStandard_WP2015(aInd);
            ccIncRef_ub = ccIncHivAgeTime_ub(a , :) .* worldStandard_WP2015(aInd);
            if (a < 3)
                ccIncRef = zeros(1 , size(ccIncHivAgeTime_med,2));
                ccIncRef_lb = zeros(1 , size(ccIncHivAgeTime_lb,2));
                ccIncRef_ub = zeros(1 , size(ccIncHivAgeTime_ub,2));
            end
        elseif aInd > age
            ccIncRef = ccIncHivAgeTime_med(a , :);
            ccIncRef = [(ones(1,aInd-a).*ccIncRef(1,1)) , ccIncRef(1,1:end-(aInd-a))];
            ccIncRef = ccIncRef .* worldStandard_WP2015(aInd);
            ccIncRef_lb = ccIncHivAgeTime_lb(a , :);
            ccIncRef_lb = [(ones(1,aInd-a).*ccIncRef_lb(1,1)) , ccIncRef_lb(1,1:end-(aInd-a))];
            ccIncRef_lb = ccIncRef_lb .* worldStandard_WP2015(aInd);
            ccIncRef_ub = ccIncHivAgeTime_ub(a , :);
            ccIncRef_ub = [(ones(1,aInd-a).*ccIncRef_ub(1,1)) , ccIncRef_ub(1,1:end-(aInd-a))];
            ccIncRef_ub = ccIncRef_ub .* worldStandard_WP2015(aInd);
        end
        ccIncRefTot = ccIncRefTot + ccIncRef;
        ccIncRefTot_lb = ccIncRefTot_lb + ccIncRef_lb;
        ccIncRefTot_ub = ccIncRefTot_ub + ccIncRef_ub;
    end
    ccInc = ccIncRefTot ./ (sum(worldStandard_WP2015(1:age+4)));
    ccInc_lb = ccIncRefTot_lb ./ (sum(worldStandard_WP2015(1:age+4)));
    ccInc_ub = ccIncRefTot_ub ./ (sum(worldStandard_WP2015(1:age+4)));

    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_CC_IncidenceRates-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(ccInc(firstYrInd:end) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B4')
    
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_CC_IncidenceRates-standardised-(2020-2120)LUB_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(ccInc_lb(firstYrInd3:end) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B4')
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_Impact_CC_IncidenceRates-standardised-(2020-2120)HUB_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(ccInc_ub(firstYrInd3:end) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B4')
end

%% Write age-standardized cervical cancer incidence rates by HIV status over time (2020-2121) into existing template
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
firstYrRange = (lastYear-1) - (currYear-1);
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

    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'CCinc_wUncert_standardized_2020-2121_S' , fileKeyNums{n} , '.xlsx'];  
    writematrix([[2020 : lastYear-1] ; ...
        squeeze(median(squeeze(ccInc(: , : , (end-firstYrRange):end)) , 1)) ; ...
        squeeze(min(squeeze(ccInc(: , : , (end-firstYrRange):end)) , [] , 1)) ; ...
        squeeze(max(squeeze(ccInc(: , : , (end-firstYrRange):end)) , [] , 1))] , ...
        fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B3') 
end

%% Cumulative cervical cancer cases by HIV status and age over time
diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'CumCC_' , diseaseLabels{dInd} , '_' , fileKey{n} , '.csv'];
    writematrix([[0 ; (1:age)' ; (1:age)' ; (1:age)'] , ...
        [futAnnualTimespan ;
        [squeeze(median(squeeze(ccCumHivAgeTime(: , dInd , : , :)) , 1)) ; ...
        squeeze(min(squeeze(ccCumHivAgeTime(: , dInd , : , :)) , [] , 1)) ; ...
        squeeze(max(squeeze(ccCumHivAgeTime(: , dInd , : , :)) , [] , 1))]]] , fname)
end

%% Write cumulative cervical cancer cases by 5-year age groups over time (2019-2120) into existing template
diseaseLabels = {'Pop(All) CCC' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)'};
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'UofW_CumulativeImpact_CC-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix(squeeze(median(squeeze(ccCumHivAgeTime(: , dInd , 3:age , :)) , 1)) , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B8')   
end

%% Write cumulative cervical cancer cases over time (2021-2121) into existing template
diseaseLabels = {'Pop(All) CCC' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)'};
for dInd = 1 : length(diseaseLabels)
    fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
        'CumulativeImpact_CC-standardised_wUncert-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
    writematrix([[currYear : lastYear-1] ; ...
        squeeze(median(squeeze(ccCumHivTime(: , dInd , :)) , 1)) ; ...
        squeeze(min(squeeze(ccCumHivTime(: , dInd , :)) , [] , 1)) ; ...
        squeeze(max(squeeze(ccCumHivTime(: , dInd , :)) , [] , 1)) ; ...
        squeeze(ccCumHivTime(: , dInd , :))] , fname , 'Sheet' , diseaseLabels{dInd})   
end

%% Write age-standardized cumulative cervical cancer cases by HIV status over time (2019-2120) into existing template
% Note: the age-standardization process shifts the incidence rate of the
% last modelled age group to the next age group in the following year.
% However, CC incidence is NaN prior to HIV introduction in the
% HIV-positive no ART group, and NaN prior to ART introduction in the
% HIV-positive ART group. Since we have four age groups past the 16 we
% model, a NaN value is present for four years past the introduction of
% HIV/ART, leading to a NaN value for summed incidence during these 
% years. We therefore lack data in this four-year interval in the
% saved/plotted results.
figure;

fac = 10 ^ 5;
worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
    247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
    23477 9261 2155];

firstYrInd = ((2019 - startYear) +1);
fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
    'UofW_CumulativeImpact_CC-standardised-(2020-2120)_S' , fileKeyNums{n} , '.xlsx'];
diseaseLabels = {'Pop(All) CCC' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)'};
for dInd = 1 : length(diseaseLabels)
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
        ccIncRefTot = ccIncRefTot + (ccIncRef ./ fac);
    end
    ccCum = cumsum(ccIncRefTot(firstYrInd:end));
    writematrix(ccCum , fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'B4')
    
%     hold all;
%     plot(futAnnualTimespan , ccCum);
end     
% xlim([2000 2120]);
% legend('Pop(All) CCC' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)');


%% ************************** SCREENING & VACCINATION FIGURES ******************************************************************************

%% Total number school-based vaccines administered
fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
    'RoutineSchoolVax_wUncert-(byYear)_S' , fileKeyNums{n} , '.xlsx'];
writematrix([squeeze(median(cumVaxSchool(: , :) , 1)) ; ...
    squeeze(min(cumVaxSchool(: , :) , [] , 1)) ; ...
    squeeze(max(cumVaxSchool(: , :) , [] , 1)) ; ...
    cumVaxSchool(: , :)] , fname , 'Sheet' , diseaseLabels{dInd})

%% Total number CU vaccines administered
fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
    'CUVax_wUncert-(byYear)_S' , fileKeyNums{n} , '.xlsx'];
writematrix([squeeze(median(cumVaxCU(: , :) , 1)) ; ...
    squeeze(min(cumVaxCU(: , :) , [] , 1)) ; ...
    squeeze(max(cumVaxCU(: , :) , [] , 1)) ; ...
    cumVaxCU(: , :)] , fname , 'Sheet' , diseaseLabels{dInd})
