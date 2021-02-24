function vaxCEA_plotSaveDoArt_071720(fileInd)

%% SET RUN-TIME VARIABLES
lastYear = 2061;
baseFileNameShort = '22Apr20Ph2V11_baseVax057_baseScreen_baseVMMC_fertDec042-076-052_2020ARTfxd_trackCD4-Discont_discontFxd';    % **** SET ME ****
baseFileName = [baseFileNameShort , '_DoART_S' , num2str(fileInd) , '_']; %'_diagHiv075_DoART_S' ,
%  Note: if not saving outputs for a simulation with explicit HIV diagnosis,
%  need to comment out sections below referencing the following variables:
%  nTestedHivNegC_multSims , nTestedNeg
%  nTestedHivUndiagC_multSims , nTestedUndiag
%  propHivDiag_multSims , propHivDiag
n = 2; % use 2nd, or vaxResult scenario with baseline HPV vaccination (does not affect this analysis)
baseHistFileNameShort = 'toNow_22Apr20Ph2V11_baseVax057_baseScreen_baseVMMC_fertDec042-076_2020ARTfxd_trackCD4-Discont_discontFxd_DoART_S1_';    % **** SET ME ****
% Indices of calib runs to plot
fileInds = {'6_1' , '6_2' , '6_3' , '6_6' , '6_8' , '6_9' , '6_11' , ...
     '6_12' , '6_13' , '6_15' , '6_20' , '6_21' , '6_22' , '6_26' , ...
    '6_27' , '6_32' , '6_34' , '6_35' , '6_38' , '6_39' , '6_40' , ...
    '6_41' , '6_42' , '6_45' , '6_47'};    % 22Apr20Ph2V11    % **** SET ME ****
% fileInds = {'11_1' , '11_2' , '11_3' , '11_4' , '11_5' , '11_6' , '11_7' , '11_8' , ...
%     '11_9' , '11_10' , '11_11' , '11_12' , '11_13' , '11_14' , '11_15' , '11_16' , ...
%     '11_17' , '11_18' , '11_19' , '11_20' , '11_21' , '11_22' , '11_23' , '11_24' , '11_25'};  % DO ART, 22Apr20Ph2V2, t=11
% fileInds = {'7_1'}; % , '7_2' , '7_3' , '7_4' , '7_5'}; % , '7_6' , '7_7' , '7_8' , ...
%     '7_9' , '7_10' , '7_11' , '7_12' , '7_13' , '7_14' , '7_15' , '7_16' , ...
%     '7_17' , '7_18' , '7_19' , '7_20' , '7_21' , '7_22' , '7_23' , '7_24' , '7_25'};  % DO ART, 22Apr20Ph2V2, t=11
nRuns = length(fileInds);
% Plot specifications and scenarios
colorList = {'k' , 'b'}; % , 'r'};
iIndList = {1 , 2}; % , 3};
sceFileNameList = {'_DoART' , '_diagHiv075_DoART'}; % , '_DoART'};

%% LOAD PARAMETERS
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

%% INITIALIZE OUTPUT VECTORS
resultsDir = [pwd , '\HHCoM_Results\'];
% GENERAL LOOP VARIABLES
agesEligVec = {4 , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14 , 15 , 16};
agesEligVecLength = length(agesEligVec);
% TIME
tVec = [startYear : timeStep : lastYear-timeStep];
tVecYr = tVec(1 : stepsPerYear : end);
% HIV INCIDENCE
hivIncF_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
hivIncM_multSims = hivIncF_multSims;
hivIncC_multSims = hivIncF_multSims;
hivIncF1559_multSims = hivIncF_multSims;
hivIncM1559_multSims = hivIncF_multSims;
hivIncAgeC_multSims = zeros(length([startYear : lastYear-1]) , nRuns , age-3);
hivIncAllAgeC_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
hivIncAllAgeF_multSims = hivIncAllAgeC_multSims;
hivIncAllAgeM_multSims = hivIncAllAgeC_multSims;
% HIV PREVALENCE
cIndsHivPrev = {3 : 8 , 3 : 5 , 6 , 7 , 8};
cIndsHivPrevLength = length(cIndsHivPrev);
hivPrevF_multSims = zeros(length([startYear : lastYear-1]) , nRuns , 5);
hivPrevM_multSims = hivPrevF_multSims;
hivPrevC_multSims = hivPrevF_multSims;
hivPrevF1559_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
hivPrevM1559_multSims = hivPrevF1559_multSims;
artPrevF_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
artPrevM_multSims = artPrevF_multSims;
artPrevC_multSims = artPrevF_multSims;
% POPULATION CD4 DISTRIBUTION
cd4DistVec = {7 , [3 : 6] , 8}; % CD4<200, CD4>=200, on ART
cd4DistVecLength = length(cd4DistVec);
cd4PropMF1559 = zeros(length(tVecYr) , nRuns , cd4DistVecLength , gender);
% ART INIT CD4 DISTRIBUTION/TRANSITIONS
cIndsCD4Dist = {3 : 5 , 6 , 7};
cIndsCD4DistLength = length(cIndsCD4Dist);
cd4DistF_multSims = zeros(length([startYear : lastYear-1]) , nRuns , 3);
cd4DistM_multSims = cd4DistF_multSims;
cd4DistC_multSims = cd4DistF_multSims;
cIndsArtOnOff = {3 : 7 , 3 : 5 , 6 , 7};
cIndsArtOnOffLength = length(cIndsArtOnOff);
cd4DistAgeC_multSims = zeros(length([startYear : lastYear-1]) , nRuns , age-3 , 4);
cd4DiscontAgeC_multSims = zeros(length([startYear : lastYear-1]) , nRuns , age-3 , 4);
cIndsCD4Trans = {6 , 7};
cIndsCD4TransLength = length(cIndsCD4Trans);
cd4TransAgeC_multSims = zeros(length([startYear : lastYear-1]) , nRuns , age-3 , 2);
% HIV-ASSOCIATED & ALL-CAUSE MORTALITY
cIndsHivMort = {3 : 8 , 3 : 5 , 6 , 7 , 8};
cIndsHivMortLength = length(cIndsHivMort);
hivMortF_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
hivMortM_multSims = hivMortF_multSims;
hivMortC_multSims = hivMortF_multSims;
hivMortAgeC_multSims = zeros(length([startYear : lastYear-1]) , nRuns , age-3 , 5);
cIndsAllMort = {1 : 2 , 3 : 8 , 3 : 5 , 6 , 7 , 8};
cIndsAllMortLength = length(cIndsAllMort);
allCauseMortAgeC_multSims = zeros(length([startYear : lastYear-1]) , nRuns , age-3 , 5);
popSizeGAD_disInds = {1 : 2 , 3 : 5 , 6 , 7 , 8};
popSizeGAD_disIndsLength = length(popSizeGAD_disInds);
allCauseDeathGAD_multSims = zeros(length([startYear : lastYear-1]) , nRuns , popSizeGAD_disIndsLength , gender , agesEligVecLength);
hivMortAllAgeC_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
hivMortAllAgeF_multSims = hivMortAllAgeC_multSims;
hivMortAllAgeM_multSims = hivMortAllAgeC_multSims;
% POPULATION SIZE
popSizeF_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
popSizeM_multSims = popSizeF_multSims;
popSizeC_multSims = popSizeF_multSims;
popSizeGAD_multSims = zeros(length([startYear : lastYear-1]) , nRuns , popSizeGAD_disIndsLength , gender , agesEligVecLength);
% VMMC
vmmcM_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
nTestedHivNegC_multSims = zeros(length([currYear : lastYear-1]) , nRuns);
nTestedHivUndiagC_multSims = nTestedHivNegC_multSims;
propHivDiag_multSims = zeros(length([currYear : timeStep : lastYear-timeStep]) , gender , nRuns);

%% LOOP THROUGH NRUNS
loopSegments = {0 , round(nRuns/2) , nRuns};
loopSegmentsLength = length(loopSegments);
for k = 1 : loopSegmentsLength-1
    parfor j = loopSegments{k}+1 : loopSegments{k+1}
        % Load results
        pathModifier = [baseFileName , fileInds{j}]; % ***SET ME***: name for simulation output file
        nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
        curr = load([pwd , '/HHCoM_Results/' , baseHistFileNameShort , fileInds{j}]); % ***SET ME***: name for historical run output file 
        vaxResult = cell(nSims , 1);
        resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxSimResult'];
        % load results from vaccine run into cell array
        vaxResult{n} = load([resultFileName , num2str(n), '.mat']);
        % concatenate vectors/matrices of population up to current year to population
        % matrices for years past current year
        vaxResult{n}.popVec = [curr.popVec(1 : end  , :); vaxResult{n}.popVec(2 : end , :)];
        vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , : , : , : , :); vaxResult{n}.newHiv(2 : end , : , : , : , : , : , :)];
        vaxResult{n}.transCD4 = [curr.transCD4(1 : end , : , : , :); vaxResult{n}.transCD4(2 : end , : , : , :)];
        vaxResult{n}.hivDeaths = [curr.hivDeaths(1 : end , : , : , :); vaxResult{n}.hivDeaths(2 : end , : , : , :)];
        vaxResult{n}.deaths = [curr.deaths(1 : end , : , : , :); vaxResult{n}.deaths(2 : end , : , : , :)];
        vaxResult{n}.ccDeath = [curr.ccDeath(1 : end , : , : , :); vaxResult{n}.ccDeath(2 : end , : , : , :)];
        vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , : , : , : , : , :); vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , :)];
        vaxResult{n}.artDiscont = [curr.artDiscont(1 : end , : , : , : , : , :); vaxResult{n}.artDiscont(2 : end , : , : , : , : , :)];
        vaxResult{n}.menCirc = [curr.menCirc(1 : end  , :); vaxResult{n}.menCirc(2 : end , :)];
%         vaxResult{n}.nTestedNeg = vaxResult{n}.nTestedNeg (1 : end , :);
%         vaxResult{n}.nTestedUndiag = vaxResult{n}.nTestedUndiag (1 : end , :);
%         vaxResult{n}.propHivDiag = vaxResult{n}.propHivDiag(1 : end , :);
        vaxResult{n}.tVec = [curr.tVec(1 : end), vaxResult{n}.tVec(2 : end)];

        %% HIV INCIDENCE
        % Validation data
%         hivInc_obs(: , : , 1) = [2005 2.14 1.57 2.93; % AHRI KZN: (Vandormael, 2019)
%                              2006 2.24 1.69 2.96;
%                              2007 2.30 1.74 3.05;
%                              2008 2.35 1.78 3.09;
%                              2009 2.45 1.85 3.24;
%                              2010 2.45 1.85 3.25;
%                              2011 2.30 1.70 3.11;
%                              2012 2.49 1.83 3.37;
%                              2013 2.22 1.64 3.01;
%                              2014 1.83 1.29 2.59;
%                              2015 1.39 0.94 2.07;
%                              2016 1.24 0.79 1.95;
%                              2017 1.01 0.58 1.76];
%         hivInc_obs(: , : , 2) = [2005 4.08 3.40 4.90;
%                              2006 4.45 3.77 5.27;
%                              2007 4.56 3.86 5.39;
%                              2008 4.58 3.89 5.40;
%                              2009 4.58 3.85 5.44;
%                              2010 4.72 3.98 5.61;
%                              2011 4.59 3.85 5.47;
%                              2012 4.95 4.14 5.92;
%                              2013 4.85 4.05 5.81;
%                              2014 4.89 4.09 5.84;
%                              2015 4.31 3.58 5.20;
%                              2016 3.74 3.04 4.61;
%                              2017 3.06 2.38 3.94];

        % Calculate female HIV incidence
        hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 4 : age , 1 : risk));
        hivSus = annlz(sum(vaxResult{n}.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
        hivIncF = annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHiv(: , : , : , : , 2 , 4 : age , ...
            1 : risk), 2), 3), 4), 6), 7)) ./ hivSus * 100;
        hivIncF_multSims(: , j) = hivIncF(1 : end)';

        % Plot female HIV incidence  
%         if (j == 1)
%             fig1 = figure;
%             %errorbar(hivInc_obs(: , 1 , 2) , hivInc_obs(: , 2 , 2) , ...
%             %    hivInc_obs(: , 2 , 2) - hivInc_obs(: , 3 , 2) , hivInc_obs(: , 4 , 2) - hivInc_obs(: , 2 , 2) , ...
%             %    'rs' , 'LineWidth' , 1.5);
%         else
%             figure(fig1);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , hivIncF(1 : end)' , 'b-')
%         xlabel('Year'); ylabel('Incidence per 100'); title('Female HIV incidence, ages 15-79');
%         xlim([1980 2060]); ylim([0 10]);
%         legend('Model: ages 15-79'); %'(Vandormael, 2019) Observed KZN, ages 15-49: 95% CI' , 
%         grid on;

        % Calculate male HIV incidence
        hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 , 4 : age , 1 : risk));
        hivSus = annlz(sum(vaxResult{n}.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
        hivIncM = annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHiv(: , : , : , : , 1 , 4 : age , ...
            1 : risk), 2), 3), 4), 6), 7)) ./ hivSus * 100;
        hivIncM_multSims(: , j) = hivIncM(1 : end)';

        % Plot male HIV incidence            
%         if (j == 1)
%             fig2 = figure;
%             %errorbar(hivInc_obs(: , 1 , 1) , hivInc_obs(: , 2 , 1) , ...
%             %    hivInc_obs(: , 2 , 1) - hivInc_obs(: , 3 , 1) , hivInc_obs(: , 4 , 1) - hivInc_obs(: , 2 , 1) , ...
%             %    'rs' , 'LineWidth' , 1.5);
%         else
%             figure(fig2);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , hivIncM(1 : end)' , 'b-')
%         xlabel('Year'); ylabel('Incidence per 100'); title('Male HIV incidence, ages 15-79');
%         xlim([1980 2060]); ylim([0 10]);
%         legend('Model: ages 15-79'); %'(Vandormael, 2019) Observed KZN, ages 15-49: 95% CI' , 
%         grid on;

        % Calculate combined HIV incidence
        hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : 2 , 4 : age , 1 : risk));
        hivSus = annlz(sum(vaxResult{n}.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
        hivIncC = annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newHiv(: , : , : , : , 1 : 2 , 4 : age , ...
            1 : risk), 2), 3), 4), 5), 6), 7)) ./ hivSus * 100;
        hivIncC_multSims(: , j) = hivIncC(1 : end)';

        % Plot combined HIV incidence            
%         if (j == 1)
%             fig13 = figure;
%         %    errorbar(hivInc_obs(: , 1 , 1) , hivInc_obs(: , 2 , 1) , ...
%         %        hivInc_obs(: , 2 , 1) - hivInc_obs(: , 3 , 1) , hivInc_obs(: , 4 , 1) - hivInc_obs(: , 2 , 1) , ...
%         %        'rs' , 'LineWidth' , 1.5);
%         else
%             figure(fig13);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , hivIncC(1 : end)' , 'b-')
%         xlabel('Year'); ylabel('Incidence per 100'); title('Combined HIV incidence, ages 15-79');
%         xlim([1980 2060]); ylim([0 10]);
%         legend('Model: ages 15-79'); %'(Vandormael, 2019) Observed KZN, ages 15-49: 95% CI' , 
%         grid on;
            
        %% HIV INCIDENCE BY GENDER, AGES 15-59, FOR HIV-TB MODEL
        % Calculate female HIV incidence
        hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 4 : 12 , 1 : risk));
        hivSus = annlz(sum(vaxResult{n}.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
        hivIncF = annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHiv(: , : , : , : , 2 , 4 : 12 , ...
            1 : risk), 2), 3), 4), 6), 7)) ./ hivSus * 100;
        hivIncF1559_multSims(: , j) = hivIncF(1 : end)';
        
        % Calculate male HIV incidence
        hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 , 4 : 12 , 1 : risk));
        hivSus = annlz(sum(vaxResult{n}.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
        hivIncM = annlz(sum(sum(sum(sum(sum(vaxResult{n}.newHiv(: , : , : , : , 1 , 4 : 12 , ...
            1 : risk), 2), 3), 4), 6), 7)) ./ hivSus * 100;
        hivIncM1559_multSims(: , j) = hivIncM(1 : end)';
        
        %% CD4 DISTRIBUTION BY GENDER, AGES 15-59 FOR HIV-TB MODEL        
        for gInd = 1 : gender
            gGroup = gInd;
            for cInd = 1 : cd4DistVecLength
                c = cd4DistVec{cInd};
                cd4Inds = toInd(allcomb(c , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                    1 : intervens , gGroup , 4 : 12 , 1 : risk));
                hivInds = toInd(allcomb(3 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                    1 : intervens , gGroup , 4 : 12 , 1 : risk));
                cd4PropMF1559(: , j , cInd , gInd) = sum(vaxResult{n}.popVec([1 : stepsPerYear : end] , cd4Inds) , 2) ./ ...
                    sum(vaxResult{n}.popVec([1 : stepsPerYear : end] , hivInds) , 2);
            end
        end

        %% HIV INCIDENCE CASES BY AGE
        % Calculate combined HIV incidence
        for aInd = 1 : agesEligVecLength
            a = agesEligVec{aInd};
            hivIncC = annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newHiv(: , : , : , : , 1 : 2 , a , 1 : risk), 2), 3), 4), 5), 6), 7));
            hivIncAgeC_multSims(: , j , aInd) = hivIncC(1 : end)';
        end

        %% HIV TOTAL CASES ACROSS AGES
        % Calculate female HIV cases
        hivIncF = annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newHiv(: , : , : , : , 2 , 4 : age , 1 : risk), 2), 3), 4), 5), 6), 7));
        hivIncAllAgeF_multSims(: , j) = hivIncF(1 : end)';

        % Calculate male HIV cases
        hivIncM = annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newHiv(: , : , : , : , 1 , 4 : age , 1 : risk), 2), 3), 4), 5), 6), 7));
        hivIncAllAgeM_multSims(: , j) = hivIncM(1 : end)';

        % Calculate combined HIV cases
        hivIncC = annlz(sum(sum(sum(sum(sum(sum(vaxResult{n}.newHiv(: , : , : , : , 1 : 2 , 4 : age , 1 : risk), 2), 3), 4), 5), 6), 7));
        hivIncAllAgeC_multSims(: , j) = hivIncC(1 : end)';

        %% NUMBER OF CIRCUMCISIONS ACROSS AGES
        circNewM = annlz(vaxResult{n}.menCirc(: , :));
        vmmcM_multSims(: , j) = circNewM(1 : end)';

        %% HIV PREVALENCE
        cTits = {'all HIV' , 'CD4 350plus noART' , 'CD4 200-350 noART' , 'CD4 below200 noART' , 'onART'};

        % Plot female HIV prevalence      
%         if (j == 1)
%             fig3 = figure;
%             %insert observed data
%         else
%             figure(fig3);
%         end
        for cInd = 1 : cIndsHivPrevLength
            cGroup = cIndsHivPrev{cInd};
            % Calculate female HIV prevalence
            hivIndsF = toInd(allcomb(cGroup , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , 2 , 4 : age , 1 : risk));
            totIndsF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , 2 , 4 : age , 1 : risk));
            hivPopF = sum(vaxResult{n}.popVec(: , hivIndsF) , 2);
            totPopF = sum(vaxResult{n}.popVec(: , totIndsF) , 2);
            hivPrevF = bsxfun(@rdivide , hivPopF , totPopF);
            hivPrevF_multSims(: , j , cInd) = hivPrevF(1 : stepsPerYear : end);

%             subplot(2 , 3 , cInd);
%             hold all;
%             plot(tVec(1 : stepsPerYear : end)' , hivPrevF(1 : stepsPerYear : end) , 'b-')
%             xlabel('Year'); ylabel('Female HIV prevalence, ages 15-79'); title(cTits{cInd});
%             xlim([1980 2060]); ylim([0.0 0.5]);
%             grid on;
        end

        % Plot male HIV prevalence      
%         if (j == 1)
%             fig4 = figure;
%             %insert observed data
%         else
%             figure(fig4);
%         end
        for cInd = 1 : cIndsHivPrevLength
            cGroup = cIndsHivPrev{cInd};
            % Calculate male HIV prevalence
            hivIndsM = toInd(allcomb(cGroup , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , 1 , 4 : age , 1 : risk));
            totIndsM = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , 1 , 4 : age , 1 : risk));
            hivPopM = sum(vaxResult{n}.popVec(: , hivIndsM) , 2);
            totPopM = sum(vaxResult{n}.popVec(: , totIndsM) , 2);
            hivPrevM = bsxfun(@rdivide , hivPopM , totPopM);
            hivPrevM_multSims(: , j , cInd) = hivPrevM(1 : stepsPerYear : end);

%             subplot(2 , 3 , cInd);
%             hold all;
%             plot(tVec(1 : stepsPerYear : end)' , hivPrevM(1 : stepsPerYear : end) , 'b-')
%             xlabel('Year'); ylabel('Male HIV prevalence, ages 15-79'); title(cTits{cInd});
%             xlim([1980 2060]); ylim([0.0 0.5]);
%             grid on;
        end

        % Plot combined HIV prevalence      
%         if (j == 1)
%             fig14 = figure;
%             %insert observed data
%         else
%             figure(fig14);
%         end
        for cInd = 1 : cIndsHivPrevLength
            cGroup = cIndsHivPrev{cInd};
            % Calculate combined HIV prevalence
            hivIndsC = toInd(allcomb(cGroup , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , 1 : 2 , 4 : age , 1 : risk));
            totIndsC = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
                1 : intervens , 1 : 2 , 4 : age , 1 : risk));
            hivPopC = sum(vaxResult{n}.popVec(: , hivIndsC) , 2);
            totPopC = sum(vaxResult{n}.popVec(: , totIndsC) , 2);
            hivPrevC = bsxfun(@rdivide , hivPopC , totPopC);
            hivPrevC_multSims(: , j , cInd) = hivPrevC(1 : stepsPerYear : end);

%             subplot(2 , 3 , cInd);
%             hold all;
%             plot(tVec(1 : stepsPerYear : end)' , hivPrevC(1 : stepsPerYear : end) , 'b-')
%             xlabel('Year'); ylabel('Combined HIV prevalence, ages 15-79'); title(cTits{cInd});
%             xlim([1980 2060]); ylim([0.0 0.5]);
%             grid on;
        end
        
        
        
        %% HIV PREVALENCE BY GENDER, AGES 15-59, FOR HIV-TB MODEL
        % Calculate female HIV prevalence
        hivIndsF = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 4 : 12 , 1 : risk));
        totIndsF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 4 : 12 , 1 : risk));
        hivPopF = sum(vaxResult{n}.popVec(: , hivIndsF) , 2);
        totPopF = sum(vaxResult{n}.popVec(: , totIndsF) , 2);
        hivPrevF = bsxfun(@rdivide , hivPopF , totPopF);
        hivPrevF1559_multSims(: , j) = hivPrevF(1 : stepsPerYear : end).*100;
        
        % Calculate male HIV prevalence
        hivIndsM = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 , 4 : 12 , 1 : risk));
        totIndsM = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 , 4 : 12 , 1 : risk));
        hivPopM = sum(vaxResult{n}.popVec(: , hivIndsM) , 2);
        totPopM = sum(vaxResult{n}.popVec(: , totIndsM) , 2);
        hivPrevM = bsxfun(@rdivide , hivPopM , totPopM);
        hivPrevM1559_multSims(: , j) = hivPrevM(1 : stepsPerYear : end).*100;

        %% ART PREVALENCE 
        % Calculate female ART prevalence
        artIndsF = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 4 : age , 1 : risk));
        hivIndsF = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 4 : age , 1 : risk));
        artPopF = sum(vaxResult{n}.popVec(: , artIndsF) , 2);
        hivPopF = sum(vaxResult{n}.popVec(: , hivIndsF) , 2);
        artPrevF = bsxfun(@rdivide , artPopF , hivPopF);
        artPrevF_multSims(: , j) = artPrevF(1 : stepsPerYear : end);

        % Plot female ART prevalence      
%         if (j == 1)
%             fig5 = figure;
%             %insert observed data
%         else
%             figure(fig5);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , artPrevF(1 : stepsPerYear : end) , 'b-')
%         xlabel('Year'); ylabel('Proportion WLWHIV on ART, ages 15-79'); title('Proportion WLWHIV on ART, ages 15-79');
%         xlim([1980 2060]); ylim([0.0 1.0]);
%         legend('Model: ages 15-79'); % Insert observed data description
%         grid on;

        % Calculate male ART prevalence
        artIndsM = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 , 4 : age , 1 : risk));
        hivIndsM = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 , 4 : age , 1 : risk));
        artPopM = sum(vaxResult{n}.popVec(: , artIndsM) , 2);
        hivPopM = sum(vaxResult{n}.popVec(: , hivIndsM) , 2);
        artPrevM = bsxfun(@rdivide , artPopM , hivPopM);
        artPrevM_multSims(: , j) = artPrevM(1 : stepsPerYear : end);

        % Plot male ART prevalence      
%         if (j == 1)
%             fig6 = figure;
%             %insert observed data
%         else
%             figure(fig6);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , artPrevM(1 : stepsPerYear : end) , 'b-')
%         xlabel('Year'); ylabel('Proportion MLWHIV on ART, ages 15-79'); title('Proportion MLWHIV on ART, ages 15-79');
%         xlim([1980 2060]); ylim([0.0 1.0]);
%         legend('Model: ages 15-79'); % Insert observed data description
%         grid on;

        % Calculate combined ART prevalence
        artIndsC = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : 2 , 4 : age , 1 : risk));
        hivIndsC = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : 2 , 4 : age , 1 : risk));
        artPopC = sum(vaxResult{n}.popVec(: , artIndsC) , 2);
        hivPopC = sum(vaxResult{n}.popVec(: , hivIndsC) , 2);
        artPrevC = bsxfun(@rdivide , artPopC , hivPopC);
        artPrevC_multSims(: , j) = artPrevC(1 : stepsPerYear : end);

        % Plot combined ART prevalence      
%         if (j == 1)
%             fig15 = figure;
%             %insert observed data
%         else
%             figure(fig15);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , artPrevC(1 : stepsPerYear : end) , 'b-')
%         xlabel('Year'); ylabel('Proportion PLWHIV on ART, ages 15-79'); title('Proportion PLWHIV on ART, ages 15-79');
%         xlim([1980 2060]); ylim([0.0 1.0]);
%         legend('Model: ages 15-79'); % Insert observed data description
%         grid on;

        %% CD4 DISTRIBUTION AT ART INITIATION
        cTits = {'CD4 350plus' , 'CD4 200-350' , 'CD4 below200'};

        % Plot female distribution
%         if (j == 1)
%             fig7 = figure;
%             %insert observed data
%         else
%             figure(fig7);
%         end
        for cInd = 1 : cIndsCD4DistLength
            cGroup = cIndsCD4Dist{cInd};
            % Calculate female distribution
            init_subF = sum(sum(sum(sum(vaxResult{n}.artTreatTracker(: , cGroup , : , 2 , 4 : age , 1 : risk), 2), 3), 5), 6);
            init_allCD4F = sum(sum(sum(sum(vaxResult{n}.artTreatTracker(: , 3 : 7 , : , 2 , 4 : age , 1 : risk), 2), 3), 5), 6);
            cd4DistF = init_subF ./ init_allCD4F;
            cd4DistF_multSims(: , j , cInd) = cd4DistF(1 : stepsPerYear : end)';

%             subplot(1 ,3 , cInd);
%             hold all;
%             plot(tVec(1 : stepsPerYear : end)' , cd4DistF(1 : stepsPerYear : end)' , 'b-')
%             xlabel('Year'); ylabel('Proportion of females initiating ART, ages 15-79'); title(['Proportion CD4: ' , cTits{cInd}]);
%             xlim([1980 2060]); ylim([0.0 1.0]);
%             grid on;
        end

        % Plot male distribution
%         if (j == 1)
%             fig8 = figure;
%             %insert observed data
%         else
%             figure(fig8);
%         end
        for cInd = 1 : cIndsCD4DistLength
            cGroup = cIndsCD4Dist{cInd};
            % Calculate male distribution
            init_subM = sum(sum(sum(sum(vaxResult{n}.artTreatTracker(: , cGroup , : , 1 , 4 : age , 1 : risk), 2), 3), 5), 6);
            init_allCD4M = sum(sum(sum(sum(vaxResult{n}.artTreatTracker(: , 3 : 7 , : , 1 , 4 : age , 1 : risk), 2), 3), 5), 6);
            cd4DistM = init_subM ./ init_allCD4M;
            cd4DistM_multSims(: , j , cInd) = cd4DistM(1 : stepsPerYear : end)';

%             subplot(1 , 3 , cInd);
%             hold all;
%             plot(tVec(1 : stepsPerYear : end)' , cd4DistM(1 : stepsPerYear : end)' , 'b-')
%             xlabel('Year'); ylabel('Proportion of males initiating ART, ages 15-79'); title(['Proportion CD4: ' , cTits{cInd}]);
%             xlim([1980 2060]); ylim([0.0 1.0]);
%             grid on;
        end

        % Plot combined distribution
%         if (j == 1)
%             fig16 = figure;
%             %insert observed data
%         else
%             figure(fig16);
%         end
        for cInd = 1 : cIndsCD4DistLength
            cGroup = cIndsCD4Dist{cInd};
            % Calculate combined distribution
            init_subC = sum(sum(sum(sum(sum(vaxResult{n}.artTreatTracker(: , cGroup , : , 1 : 2 , 4 : age , 1 : risk), 2), 3), 4), 5), 6);
            init_allCD4C = sum(sum(sum(sum(sum(vaxResult{n}.artTreatTracker(: , 3 : 7 , : , 1 : 2 , 4 : age , 1 : risk), 2), 3), 4), 5), 6);
            cd4DistC = init_subC ./ init_allCD4C;
            cd4DistC_multSims(: , j , cInd) = cd4DistC(1 : stepsPerYear : end)';

%             subplot(1 , 3 , cInd);
%             hold all;
%             plot(tVec(1 : stepsPerYear : end)' , cd4DistC(1 : stepsPerYear : end)' , 'b-')
%             xlabel('Year'); ylabel('Proportion of persons initiating ART, ages 15-79'); title(['Proportion CD4: ' , cTits{cInd}]);
%             xlim([1980 2060]); ylim([0.0 1.0]);
%             grid on;
        end

        %% ART INITIATION BY CD4 AND AGE
        % Calculate combined numbers initiating ART
        for cInd = 1 : cIndsArtOnOffLength
            cGroup = cIndsArtOnOff{cInd};
            for aInd = 1 : agesEligVecLength
                a = agesEligVec{aInd};
                init_subC = annlz(sum(sum(sum(sum(sum(vaxResult{n}.artTreatTracker(: , cGroup , : , 1 : 2 , a , 1 : risk), 2), 3), 4), 5), 6));
                cd4DistAgeC_multSims(: , j , aInd , cInd) = init_subC(1 : end)';
            end
        end

        %% ART DISCONTINUATION BY CD4 AND AGE
        % Calculate combined numbers discontinuing ART
        for cInd = 1 : cIndsArtOnOffLength
            cGroup = cIndsArtOnOff{cInd};
            for aInd = 1 : agesEligVecLength
                a = agesEligVec{aInd};
                discont_subC = annlz(sum(sum(sum(sum(sum(vaxResult{n}.artDiscont(: , cGroup , : , 1 : 2 , a , 1 : risk), 2), 3), 4), 5), 6));
                cd4DiscontAgeC_multSims(: , j , aInd , cInd) = discont_subC(1 : end)';
            end
        end

        %% CD4 TRANSITIONS BY AGE
        % Calculate combined
        for cInd = 1 : cIndsCD4TransLength
            cGroup = cIndsCD4Trans{cInd};
            for aInd = 1 : agesEligVecLength
                a = agesEligVec{aInd};
                init_subC = annlz(sum(sum(sum(vaxResult{n}.transCD4(: , cGroup , 1 : 2 , a), 2), 3), 4));
                cd4TransAgeC_multSims(: , j , aInd , cInd) = init_subC(1 : end)';
            end
        end

        %% HIV-ASSOCIATED MORTALITY

        % Calculate female HIV-associated mortality
        popIndsF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 4 : age , 1 : risk));
        popTotF = annlz(sum(vaxResult{n}.popVec(: , popIndsF) , 2)) ./ stepsPerYear;
        hivMortF = annlz(sum(sum(vaxResult{n}.hivDeaths(: , : , 2 , 4 : age), 2), 4)) ./ popTotF * 100000;
        hivMortF_multSims(: , j) = hivMortF(1 : end)';

        % Plot female HIV-associated mortality
%         if (j == 1)
%             fig9 = figure;
%             %insert observed data
%         else
%             figure(fig9);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , hivMortF(1 : end)' , 'b-')
%         xlabel('Year'); ylabel('Mortality per 100K'); title('Female HIV-associated mortality, ages 15-79');
%         xlim([1980 2060]); ylim([0 5000]);
%         legend('Model: ages 15-79');
%         grid on;

        % Calculate male HIV-associated mortality
        popIndsM = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 , 4 : age , 1 : risk));
        popTotM = annlz(sum(vaxResult{n}.popVec(: , popIndsM) , 2)) ./ stepsPerYear;
        hivMortM = annlz(sum(sum(vaxResult{n}.hivDeaths(: , : , 1 , 4 : age), 2), 4)) ./ popTotM * 100000;
        hivMortM_multSims(: , j) = hivMortM(1 : end)';

        % Plot male HIV-associated mortality
%         if (j == 1)
%             fig10 = figure;
%             %insert observed data
%         else
%             figure(fig10);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , hivMortM(1 : end)' , 'b-')
%         xlabel('Year'); ylabel('Mortality per 100K'); title('Male HIV-associated mortality, ages 15-79');
%         xlim([1980 2060]); ylim([0 5000]);
%         legend('Model: ages 15-79');
%         grid on;

        % Calculate combined HIV-associated mortality
        popIndsC = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : 2 , 4 : age , 1 : risk));
        popTotC = annlz(sum(vaxResult{n}.popVec(: , popIndsC) , 2)) ./ stepsPerYear;
        hivMortC = annlz(sum(sum(sum(vaxResult{n}.hivDeaths(: , : , 1 : 2 , 4 : age), 2), 3), 4)) ./ popTotC * 100000;
        hivMortC_multSims(: , j) = hivMortC(1 : end)';

        % Plot combined HIV-associated mortality
%         if (j == 1)
%             fig17 = figure;
%             %insert observed data
%         else
%             figure(fig17);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , hivMortC(1 : end)' , 'b-')
%         xlabel('Year'); ylabel('Mortality per 100K'); title('Combined HIV-associated mortality, ages 15-79');
%         xlim([1980 2060]); ylim([0 5000]);
%         legend('Model: ages 15-79');
%         grid on;

        %% HIV-ASSOCIATED MORTALITY BY AGE
        % Calculate combined HIV-associated mortality
        for cInd = 1 : cIndsHivMortLength
            cGroup = cIndsHivMort{cInd};
            for aInd = 1 : agesEligVecLength
                a = agesEligVec{aInd};
                hivMortC = annlz(sum(sum(sum(vaxResult{n}.hivDeaths(: , cGroup , 1 : 2 , a), 2), 3), 4));
                hivMortAgeC_multSims(: , j , aInd , cInd) = hivMortC(1 : end)';
            end
        end

        %% ALL-CAUSE DEATHS BY AGE
        % Calculate combined HIV-associated mortality
        for cInd = 1 : cIndsAllMortLength
            cGroup = cIndsAllMort{cInd};
            for aInd = 1 : agesEligVecLength
                a = agesEligVec{aInd};
                allCauseMortC = annlz(sum(sum(sum(vaxResult{n}.deaths(: , cGroup , 1 : 2 , a) ,2) ,3) ,4)) + annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , cGroup , a , 1 : hpvTypeGroups), 2), 3), 4)) + ...
                    annlz(sum(sum(sum(vaxResult{n}.hivDeaths(: , cGroup , 1 : 2 , a), 2), 3), 4));
                allCauseMortAgeC_multSims(: , j , aInd , cInd) = allCauseMortC(1 : end)';
            end
        end
        
        %% ALL-CAUSE DEATHS BY GENDER, AGE, HIV/ART STATUS   
        for gInd = 1 : 2
            g = gInd;
            for aInd = 1 : agesEligVecLength
                a = agesEligVec{aInd};
                for dInd = 1 : popSizeGAD_disIndsLength
                    d = popSizeGAD_disInds{dInd};
                    allCauseDeathGAD = annlz(sum(sum(sum(vaxResult{n}.deaths(: , d , g , a) ,2) ,3) ,4)) + ...
                        (g-1).*annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , d , a , 1 : hpvTypeGroups), 2), 3), 4)) + ...
                        annlz(sum(sum(sum(vaxResult{n}.hivDeaths(: , d , g , a), 2), 3), 4));
                    allCauseDeathGAD_multSims(: , j , dInd , gInd , aInd) = allCauseDeathGAD';
                end
            end
        end

        %% HIV-ASSOCIATED DEATHS ACROSS ALL AGES
        % Calculate female HIV-associated mortality 
        hivMortF = annlz(sum(sum(sum(vaxResult{n}.hivDeaths(: , 3 : 8 , 2 , 4 : age), 2), 3), 4));
        hivMortAllAgeF_multSims(: , j) = hivMortF(1 : end)';

        % Calculate male HIV-associated mortality 
        hivMortM = annlz(sum(sum(sum(vaxResult{n}.hivDeaths(: , 3 : 8 , 1 , 4 : age), 2), 3), 4));
        hivMortAllAgeM_multSims(: , j) = hivMortM(1 : end)';

        % Calculate combined HIV-associated mortality 
        hivMortC = annlz(sum(sum(sum(vaxResult{n}.hivDeaths(: , 3 : 8 , 1 : 2 , 4 : age), 2), 3), 4));
        hivMortAllAgeC_multSims(: , j) = hivMortC(1 : end)';
        
%         %% NUMBER ADDITIONAL TESTED DURING HOME TESTING CAMPAIGNS
%         
%         nTestedHivNegC_multSims(: , j) = annlz(sum(vaxResult{n}.nTestedNeg(: , 1 : gender) , 2));
%         nTestedHivUndiagC_multSims(: , j) = annlz(sum(vaxResult{n}.nTestedUndiag(: , 1 : gender) , 2));
%         
%         %% PROPORTION TESTED DURING HOME TESTING CAMPAIGNS
%         
%         for g = 1 : gender
%             propHivDiag_multSims(: , g , j) = vaxResult{n}.propHivDiag(: , g);
%         end
        
        %% TOTAL POPULATION SIZE

        % Calculate female population size
        popTotF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 4 : age , 1 : risk));
        popSizeF = sum(vaxResult{n}.popVec(: , popTotF) , 2);
        popSizeF_multSims(: , j) = popSizeF(1 : stepsPerYear : end);

        % Plot female population size
%         if (j == 1)
%             fig11 = figure;
%             %insert observed data
%         else
%             figure(fig11);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , popSizeF(1 : stepsPerYear : end) , 'b-')
%         xlabel('Year'); ylabel('Individuals'); title('Female population size, ages 15-79');
%         xlim([1980 2060])
%         legend('Model: ages 15-79');
%         grid on;

        % Calculate male population size
        popTotM = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , 4 : age , 1 : risk));
        popSizeM = sum(vaxResult{n}.popVec(: , popTotM) , 2);
        popSizeM_multSims(: , j) = popSizeM(1 : stepsPerYear : end);

        % Plot female population size
%         if (j == 1)
%             fig12 = figure;
%             %insert observed data
%         else
%             figure(fig12);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , popSizeM(1 : stepsPerYear : end) , 'b-')
%         xlabel('Year'); ylabel('Individuals'); title('Male population size, ages 15-79');
%         xlim([1980 2060])
%         legend('Model: ages 15-79');
%         grid on;

        % Calculate combined population size
        popTotC = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : 2 , 4 : age , 1 : risk));
        popSizeC = sum(vaxResult{n}.popVec(: , popTotC) , 2);
        popSizeC_multSims(: , j) = popSizeC(1 : stepsPerYear : end);

        % Plot combined population size
%         if (j == 1)
%             fig18 = figure;
%             %insert observed data
%         else
%             figure(fig18);
%         end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , popSizeC(1 : stepsPerYear : end) , 'b-')
%         xlabel('Year'); ylabel('Individuals'); title('Combined population size, ages 15-79');
%         xlim([1980 2060])
%         legend('Model: ages 15-79');
%         grid on;

        %% TOTAL POPULATION SIZE BY GENDER, AGE, HIV/ART STATUS   
        for gInd = 1 : 2
            g = gInd;
            for aInd = 1 : agesEligVecLength
                a = agesEligVec{aInd};
                for dInd = 1 : popSizeGAD_disIndsLength
                    d = popSizeGAD_disInds{dInd};
                    popSizeGADinds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                        1 : endpoints , 1 : intervens , g , a , 1 : risk));
                    popSizeGAD_multSims(: , j , dInd , gInd , aInd) = sum(vaxResult{n}.popVec((1 : stepsPerYear:end) , popSizeGADinds) , 2);
                end
            end
        end

    end
end
    
%% Save HIV incidence
% female
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'HIV_incidence_females_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivIncF_multSims , 2) , ...
    min(hivIncF_multSims , [] , 2) , max(hivIncF_multSims , [] , 2) , ...
    hivIncF_multSims] , fname)
% male
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'HIV_incidence_males_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivIncM_multSims , 2) , ...
    min(hivIncM_multSims , [] , 2) , max(hivIncM_multSims , [] , 2) , ...
    hivIncM_multSims] , fname)
% combined
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'HIV_incidence_combined_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivIncC_multSims , 2) , ...
    min(hivIncC_multSims , [] , 2) , max(hivIncC_multSims , [] , 2) , ...
    hivIncC_multSims] , fname)

%% Save HIV incidence by gender, ages 15-59, for HIV-TB model
% female
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'HIV_incidence_females_aged15-59_forHivTB' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivIncF1559_multSims , 2) , ...
    min(hivIncF1559_multSims , [] , 2) , max(hivIncF1559_multSims , [] , 2) , ...
    hivIncF1559_multSims] , fname)
% male
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'HIV_incidence_males_aged15-59_forHivTB' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivIncM1559_multSims , 2) , ...
    min(hivIncM1559_multSims , [] , 2) , max(hivIncM1559_multSims , [] , 2) , ...
    hivIncM1559_multSims] , fname)  

%% Save CD4 distribution by gender, ages 15-59 for HIV-TB model
outputVec = [];
for g = 1 : gender
    for cInd = 1 : cd4DistVecLength
       outputVec = [outputVec; ...
           [tVecYr' , ones(length(tVecYr),1).*g , ...
           ones(length(tVecYr),1).*cInd , ...
           mean(squeeze(cd4PropMF1559(: , : , cInd , g)) , 2) , ...
           min(squeeze(cd4PropMF1559(: , : , cInd , g)) , [] , 2) , ...
           max(squeeze(cd4PropMF1559(: , : , cInd , g)) , [] , 2) , ...
           squeeze(cd4PropMF1559(: , : , cInd , g))]];
    end
end
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'CD4distribution_MF_aged15-59_forHivTB' , '.csv'];
writematrix(outputVec , fname) 

%% Save HIV incidence by age
ageTits = {'15-19' , '20-24' , '25-29' , '30-34' , '35-39' , '40-44' , '45-49' , ...
    '50-54' , '55-59' , '60-64' , '65-69' , '70-74' , '75-79'};

% combined
for a = 4 : age
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Raw_HIV_incidence_combined_' , ageTits{a-3} , '.csv'];
%     fullMatrix = [0 ; tVec(1 : stepsPerYear : end)'];
%     fullMatrix = [[[ageTits{a-3} , ': mean'] , 'min' , 'max' , '1' , '2' , '3' , '4' , '5' ,'6' , ...
%         '7' , '8' , '9' , '10' , '11' , '12' , '13' , '14' , '15' , '16' , '17' , ...
%         '18' , '19' , '20' , '21' , '22' , '23' , '24' , '25'] ; [fullMatrix , mean(squeeze(hivIncAgeC_multSims(: , : , a-3)) , 2) , ...
%     min(squeeze(hivIncAgeC_multSims(: , : , a-3)) , [] , 2) , max(squeeze(hivIncAgeC_multSims(: , : , a-3)) , [] , 2) , ...
%     hivIncAgeC_multSims(: , : , a-3)]]
    writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(hivIncAgeC_multSims(: , : , a-3)) , 2) , ...
        min(squeeze(hivIncAgeC_multSims(: , : , a-3)) , [] , 2) , max(squeeze(hivIncAgeC_multSims(: , : , a-3)) , [] , 2) , ...
        squeeze(hivIncAgeC_multSims(: , : , a-3))] , fname)  
end

%% Save total HIV cases across ages
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Raw_HIV_incidence_females_ages15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(hivIncAllAgeF_multSims(: , :)) , 2) , ...
    min(squeeze(hivIncAllAgeF_multSims(: , :)) , [] , 2) , max(squeeze(hivIncAllAgeF_multSims(: , :)) , [] , 2) , ...
    squeeze(hivIncAllAgeF_multSims(: , :))] , fname)  

fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Raw_HIV_incidence_males_ages15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(hivIncAllAgeM_multSims(: , :)) , 2) , ...
    min(squeeze(hivIncAllAgeM_multSims(: , :)) , [] , 2) , max(squeeze(hivIncAllAgeM_multSims(: , :)) , [] , 2) , ...
    squeeze(hivIncAllAgeM_multSims(: , :))] , fname)  

fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Raw_HIV_incidence_combined_ages15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(hivIncAllAgeC_multSims(: , :)) , 2) , ...
    min(squeeze(hivIncAllAgeC_multSims(: , :)) , [] , 2) , max(squeeze(hivIncAllAgeC_multSims(: , :)) , [] , 2) , ...
    squeeze(hivIncAllAgeC_multSims(: , :))] , fname)  

%% Save number of circumcision across ages
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Raw_VMMC_male_ages15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(vmmcM_multSims(: , :)) , 2) , ...
    min(squeeze(vmmcM_multSims(: , :)) , [] , 2) , max(squeeze(vmmcM_multSims(: , :)) , [] , 2) , ...
    squeeze(vmmcM_multSims(: , :))] , fname)

%% Save HIV prevalence
cInds = {3 : 8 , 3 : 5 , 6 , 7 , 8};
cTits = {'all_HIV' , 'CD4_350plus_noART' , 'CD4_200-350_noART' , 'CD4_below200_noART' , 'onART'};
for c = 1 : length(cInds)
    % female
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'HIV_prevalence_females_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(hivPrevF_multSims(: , : , c) , 2) , ...
        min(hivPrevF_multSims(: , : , c) , [] , 2) , max(hivPrevF_multSims(: , : , c) , [] , 2) , ...
        hivPrevF_multSims(: , : , c)] , fname)
    % male
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'HIV_prevalence_males_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(hivPrevM_multSims(: , : , c) , 2) , ...
        min(hivPrevM_multSims(: , : , c) , [] , 2) , max(hivPrevM_multSims(: , : , c) , [] , 2) , ...
        hivPrevM_multSims(: , : , c)] , fname)
    % combined
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'HIV_prevalence_combined_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(hivPrevC_multSims(: , : , c) , 2) , ...
        min(hivPrevC_multSims(: , : , c) , [] , 2) , max(hivPrevC_multSims(: , : , c) , [] , 2) , ...
        hivPrevC_multSims(: , : , c)] , fname)
end

%% Save HIV prevalence by gender, ages 15-59, for HIV-TB model
% female
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'HIV_prevalence_females_aged15-59_forHivTB' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivPrevF1559_multSims(: , :) , 2) , ...
    min(hivPrevF1559_multSims(: , :) , [] , 2) , max(hivPrevF1559_multSims(: , :) , [] , 2) , ...
    hivPrevF1559_multSims(: , :)] , fname)
% male
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'HIV_prevalence_males_aged15-59_forHivTB' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivPrevM1559_multSims(: , :) , 2) , ...
    min(hivPrevM1559_multSims(: , :) , [] , 2) , max(hivPrevM1559_multSims(: , :) , [] , 2) , ...
    hivPrevM1559_multSims(: , :)] , fname)

%% Save ART prevalence
% female
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Proportion_WLWHIV_onART_aged15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(artPrevF_multSims , 2) , ...
    min(artPrevF_multSims , [] , 2) , max(artPrevF_multSims , [] , 2) , ...
    artPrevF_multSims] , fname)
% male
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Proportion_MLWHIV_onART_aged15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(artPrevM_multSims , 2) , ...
    min(artPrevM_multSims , [] , 2) , max(artPrevM_multSims , [] , 2) , ...
    artPrevM_multSims] , fname)
% combined
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Proportion_PLWHIV_onART_aged15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(artPrevC_multSims , 2) , ...
    min(artPrevC_multSims , [] , 2) , max(artPrevC_multSims , [] , 2) , ...
    artPrevC_multSims] , fname)

%% Save CD4 distribution at ART initiation
cInds = {3 : 5 , 6 , 7};
cTits = {'CD4_350plus' , 'CD4_200-350' , 'CD4_below200'};
for c = 1 : length(cInds)
    % female
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'CD4_dist_initART_females_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(cd4DistF_multSims(: , : , c) , 2) , ...
        min(cd4DistF_multSims(: , : , c) , [] , 2) , max(cd4DistF_multSims(: , : , c) , [] , 2) , ...
        cd4DistF_multSims(: , : , c)] , fname)
    % male 
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'CD4_dist_initART_males_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(cd4DistM_multSims(: , : , c) , 2) , ...
        min(cd4DistM_multSims(: , : , c) , [] , 2) , max(cd4DistM_multSims(: , : , c) , [] , 2) , ...
        cd4DistM_multSims(: , : , c)] , fname)
    % combined 
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
        'CD4_dist_initART_combined_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(cd4DistC_multSims(: , : , c) , 2) , ...
        min(cd4DistC_multSims(: , : , c) , [] , 2) , max(cd4DistC_multSims(: , : , c) , [] , 2) , ...
        cd4DistC_multSims(: , : , c)] , fname)
end

%% Save ART initiation by CD4 and age
cInds = {3 : 7 , 3 : 5 , 6 , 7};
cTits = {'Any_CD4' , 'CD4_350plus' , 'CD4_200-350' , 'CD4_below200'};
ageTits = {'15-19' , '20-24' , '25-29' , '30-34' , '35-39' , '40-44' , '45-49' , ...
    '50-54' , '55-59' , '60-64' , '65-69' , '70-74' , '75-79'};

for c = 1 : length(cInds)
    for a = 4 : age 
        % combined 
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
            'Raw_ART_incidence_combined_' , ageTits{a-3} , '_' , cTits{c} , '.csv'];
        writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(cd4DistAgeC_multSims(: , : , a-3 , c)) , 2) , ...
            min(squeeze(cd4DistAgeC_multSims(: , : , a-3 , c)) , [] , 2) , max(squeeze(cd4DistAgeC_multSims(: , : , a-3 , c)) , [] , 2) , ...
            cd4DistAgeC_multSims(: , : , a-3 , c)] , fname)
    end
end

%% Save ART discontinuation by CD4 and age
cInds = {3 : 7 , 3 : 5 , 6 , 7};
cTits = {'Any_CD4' , 'CD4_350plus' , 'CD4_200-350' , 'CD4_below200'};
ageTits = {'15-19' , '20-24' , '25-29' , '30-34' , '35-39' , '40-44' , '45-49' , ...
    '50-54' , '55-59' , '60-64' , '65-69' , '70-74' , '75-79'};

for c = 1 : length(cInds)
    for a = 4 : age 
        % combined 
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
            'Raw_ART_discont_combined_' , ageTits{a-3} , '_' , cTits{c} , '.csv'];
        writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(cd4DiscontAgeC_multSims(: , : , a-3 , c)) , 2) , ...
            min(squeeze(cd4DiscontAgeC_multSims(: , : , a-3 , c)) , [] , 2) , max(squeeze(cd4DiscontAgeC_multSims(: , : , a-3 , c)) , [] , 2) , ...
            cd4DiscontAgeC_multSims(: , : , a-3 , c)] , fname)
    end
end

%% Save CD4 transitions by age
cInds = {6 , 7};
cTits = {'CD4_200-350' , 'CD4_below200'};
ageTits = {'15-19' , '20-24' , '25-29' , '30-34' , '35-39' , '40-44' , '45-49' , ...
    '50-54' , '55-59' , '60-64' , '65-69' , '70-74' , '75-79'};

for c = 1 : length(cInds)
    for a = 4 : age 
        % combined 
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
            'Raw_CD4_incidence_combined_' , ageTits{a-3} , '_' , cTits{c} , '.csv'];
        writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(cd4TransAgeC_multSims(: , : , a-3 , c)) , 2) , ...
            min(squeeze(cd4TransAgeC_multSims(: , : , a-3 , c)) , [] , 2) , max(squeeze(cd4TransAgeC_multSims(: , : , a-3 , c)) , [] , 2) , ...
            cd4TransAgeC_multSims(: , : , a-3 , c)] , fname)
    end
end

%% Save HIV mortality
% female
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'HIV_mortality_females_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivMortF_multSims , 2) , ...
    min(hivMortF_multSims , [] , 2) , max(hivMortF_multSims , [] , 2) , ...
    hivMortF_multSims] , fname)
% male
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'HIV_mortality_males_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivMortM_multSims , 2) , ...
    min(hivMortM_multSims , [] , 2) , max(hivMortM_multSims , [] , 2) , ...
    hivMortM_multSims] , fname)
% combined
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'HIV_mortality_combined_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivMortC_multSims , 2) , ...
    min(hivMortC_multSims , [] , 2) , max(hivMortC_multSims , [] , 2) , ...
    hivMortC_multSims] , fname)

%% Save HIV mortality by age
cInds = {3 : 8 , 3 : 5 , 6 , 7 , 8};
cTits = {'all_HIV' , 'CD4_350plus_noART' , 'CD4_200-350_noART' , 'CD4_below200_noART' , 'onART'};
ageTits = {'15-19' , '20-24' , '25-29' , '30-34' , '35-39' , '40-44' , '45-49' , ...
    '50-54' , '55-59' , '60-64' , '65-69' , '70-74' , '75-79'};

% combined
for c = 1 : length(cInds)
    for a = 4 : age
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
            'Raw_HIV_mortality_combined_' , ageTits{a-3} , '_' , cTits{c} , '.csv'];
        writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(hivMortAgeC_multSims(: , : , a-3 , c)) , 2) , ...
            min(squeeze(hivMortAgeC_multSims(: , : , a-3 , c)) , [] , 2) , max(squeeze(hivMortAgeC_multSims(: , : , a-3 , c)) , [] , 2) , ...
            squeeze(hivMortAgeC_multSims(: , : , a-3 , c))] , fname)
    end
end

%% All-cause mortality by age
cInds = {1 : 2 , 3 : 8 , 3 : 5 , 6 , 7 , 8};
cTits = {'HIV_neg' , 'all_HIV' , 'CD4_350plus_noART' , 'CD4_200-350_noART' , 'CD4_below200_noART' , 'onART'};
ageTits = {'15-19' , '20-24' , '25-29' , '30-34' , '35-39' , '40-44' , '45-49' , ...
    '50-54' , '55-59' , '60-64' , '65-69' , '70-74' , '75-79'};

% combined
for c = 1 : length(cInds)
    for a = 4 : age
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
            'Raw_allCause_mortality_combined_' , ageTits{a-3} , '_' , cTits{c} , '.csv'];
        writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(allCauseMortAgeC_multSims(: , : , a-3 , c)) , 2) , ...
            min(squeeze(allCauseMortAgeC_multSims(: , : , a-3 , c)) , [] , 2) , max(squeeze(allCauseMortAgeC_multSims(: , : , a-3 , c)) , [] , 2) , ...
            squeeze(allCauseMortAgeC_multSims(: , : , a-3 , c))] , fname)
    end
end

%% Save total HIV deaths across ages
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Raw_HIV_mortality_females_ages15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(hivMortAllAgeF_multSims(: , :)) , 2) , ...
    min(squeeze(hivMortAllAgeF_multSims(: , :)) , [] , 2) , max(squeeze(hivMortAllAgeF_multSims(: , :)) , [] , 2) , ...
    squeeze(hivMortAllAgeF_multSims(: , :))] , fname)

fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Raw_HIV_mortality_males_ages15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(hivMortAllAgeM_multSims(: , :)) , 2) , ...
    min(squeeze(hivMortAllAgeM_multSims(: , :)) , [] , 2) , max(squeeze(hivMortAllAgeM_multSims(: , :)) , [] , 2) , ...
    squeeze(hivMortAllAgeM_multSims(: , :))] , fname)

fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'Raw_HIV_mortality_combined_ages15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(squeeze(hivMortAllAgeC_multSims(: , :)) , 2) , ...
    min(squeeze(hivMortAllAgeC_multSims(: , :)) , [] , 2) , max(squeeze(hivMortAllAgeC_multSims(: , :)) , [] , 2) , ...
    squeeze(hivMortAllAgeC_multSims(: , :))] , fname)

% %% Save additional number tested with home testing campaigns
% fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
%     'Raw_HTC_combined_hivNeg.csv'];
% writematrix([[currYear : lastYear-1]' , mean(squeeze(nTestedHivNegC_multSims(: , :)) , 2) , ...
%     min(squeeze(nTestedHivNegC_multSims(: , :)) , [] , 2) , max(squeeze(nTestedHivNegC_multSims(: , :)) , [] , 2) , ...
%     squeeze(nTestedHivNegC_multSims(: , :))] , fname)
% 
% fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
%     'Raw_HTC_combined_hivUndiag.csv'];
% writematrix([[currYear : lastYear-1]' , mean(squeeze(nTestedHivUndiagC_multSims(: , :)) , 2) , ...
%     min(squeeze(nTestedHivUndiagC_multSims(: , :)) , [] , 2) , max(squeeze(nTestedHivUndiagC_multSims(: , :)) , [] , 2) , ...
%     squeeze(nTestedHivUndiagC_multSims(: , :))] , fname)

%% Save population size
% female
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'PopulationSize_females_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(popSizeF_multSims , 2) , ...
    min(popSizeF_multSims , [] , 2) , max(popSizeF_multSims , [] , 2) , ...
    popSizeF_multSims] , fname)
% male 
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'PopulationSize_males_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(popSizeM_multSims , 2) , ...
    min(popSizeM_multSims , [] , 2) , max(popSizeM_multSims , [] , 2) , ...
    popSizeM_multSims] , fname)
% combined
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
    'PopulationSize_combined_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(popSizeC_multSims , 2) , ...
    min(popSizeC_multSims , [] , 2) , max(popSizeC_multSims , [] , 2) , ...
    popSizeC_multSims] , fname)

%% Save population size and deaths by gender, age, and HIV/ART status to existing template
ageLabelVec = {17 , 22 , 27 , 32 , 37 , 42 , 47 , 52 , 57 , 62 , 67 , 72 , 77};
hivStatusVec = {0 , 1 , 1 , 1 , 1};
cd4CountVec = {0 , 3 , 2 , 1 , 0};
artStatusVec = {0 , 0 , 0 , 0 , 1};
firstYrInd = ((2020 - startYear)*stepsPerYear +1);
t20on = tVec(firstYrInd:stepsPerYear:end)';
t20onLen = length(t20on);
firstYrAnlInd = ((2020 - startYear) +1);
outputVec = [];

for gInd = 1 : gender
    for aInd = 1 : agesEligVecLength
        for dInd = 1 : popSizeGAD_disIndsLength
           outputVec = [outputVec; ...
               [t20on , ones(t20onLen,1).*(gInd-1) , ...
               ones(t20onLen,1).*ageLabelVec{aInd} , ones(t20onLen,1).*hivStatusVec{dInd} , ...
               ones(t20onLen,1).*cd4CountVec{dInd} , ones(t20onLen,1).*artStatusVec{dInd} , ...
               mean(squeeze(popSizeGAD_multSims((firstYrAnlInd:end) , : , dInd , gInd , aInd)),2) , ...
               mean(squeeze(allCauseDeathGAD_multSims((firstYrAnlInd:end) , : , dInd , gInd , aInd)),2) , ...
               zeros(t20onLen,1)]];
        end
    end
end

for n = 1 : nRuns
    for gInd = 1 : gender
        for aInd = 1 : agesEligVecLength
            for dInd = 1 : popSizeGAD_disIndsLength
               outputVec = [outputVec; ...
                   [t20on , ones(t20onLen,1).*(gInd-1) , ...
                   ones(t20onLen,1).*ageLabelVec{aInd} , ones(t20onLen,1).*hivStatusVec{dInd} , ...
                   ones(t20onLen,1).*cd4CountVec{dInd} , ones(t20onLen,1).*artStatusVec{dInd} , ...
                   squeeze(popSizeGAD_multSims((firstYrAnlInd:end) , n , dInd , gInd , aInd)) , ...
                   squeeze(allCauseDeathGAD_multSims((firstYrAnlInd:end) , n , dInd , gInd , aInd)) , ...
                   (ones(t20onLen,1).*n)]];
            end
        end
    end
end

fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , fileInds{1} , '\' , ...
'do_art_daly_template_4allen_S' , num2str(fileInd) , '.xlsx'];
writematrix(outputVec , fname , 'Range' , 'A2')


%% PLOT SCENARIOS COMPARED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot proportion tested during home testing campaigns
figure('DefaultAxesFontSize' , 18);

subplot(1 , 2 , 1);
plot([currYear : timeStep : lastYear-timeStep] , mean(squeeze(propHivDiag_multSims(: , 2 , :)) , 2)' , 'LineStyle' , '-' , 'Color' , colorList{2} , 'LineWidth' , 2);
hold all;
x2 = [[currYear : timeStep : lastYear-timeStep] , fliplr([currYear : timeStep : lastYear-timeStep])];
inBetween = [max(squeeze(propHivDiag_multSims(: , 2 , :)) , [] , 2)' , fliplr(min(squeeze(propHivDiag_multSims(: , 2 , :)) , [] , 2)')];
h = fill(x2 , inBetween , 'k');
h.FaceColor = colorList{2};
h.EdgeColor = colorList{2};
h.FaceAlpha = 0.3;
h.LineStyle = '--';
xlabel('Year'); ylabel('Proportion WLWHIV diagnosed - S2'); title('Females');
xlim([2020 2060]); ylim([0 1]);
grid on;

subplot(1 , 2 , 2);
plot([currYear : timeStep : lastYear-timeStep] , mean(squeeze(propHivDiag_multSims(: , 1 , :)) , 2)' , 'LineStyle' , '-' , 'Color' , colorList{2} , 'LineWidth' , 2);
hold all;
x2 = [[currYear : timeStep : lastYear-timeStep] , fliplr([currYear : timeStep : lastYear-timeStep])];
inBetween = [max(squeeze(propHivDiag_multSims(: , 1 , :)) , [] , 2)' , fliplr(min(squeeze(propHivDiag_multSims(: , 1 , :)) , [] , 2)')];
h = fill(x2 , inBetween , 'k');
h.FaceColor = colorList{2};
h.EdgeColor = colorList{2};
h.FaceAlpha = 0.3;
h.LineStyle = '--';
xlabel('Year'); ylabel('Proportion MLWHIV diagnosed - S2'); title('Males');
xlim([2020 2060]); ylim([0 1]);
grid on;

%% Plot HIV incidence
figure('DefaultAxesFontSize' , 18);

% female
subplot(1 , 4 , 1);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_incidence_females_aged15-79' , '.csv'];
    hivIncF = xlsread(fname);
    if i == 1
        hold all;
        plot(hivIncF(: , 1) , hivIncF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivIncF(: , 1)' , fliplr(hivIncF(: , 1)')];
        inBetween = [hivIncF(: , 4)' , fliplr(hivIncF(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.3;
        h.LineStyle = '--';
        %plot(hivIncF(: , 1) , hivIncF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivIncF(: , 1) , hivIncF(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    else
        hold all;
        plot(hivIncF((2019-startYear)+1 : end , 1) , hivIncF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivIncF((2019-startYear)+1 : end , 1)' , fliplr(hivIncF((2019-startYear)+1 : end , 1)')];
        inBetween = [hivIncF((2019-startYear)+1 : end , 4)' , fliplr(hivIncF((2019-startYear)+1 : end , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivIncF((2019-startYear)+1 : end , 1) , hivIncF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivIncF((2019-startYear)+1 : end , 1) , hivIncF((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    end
end
ylabel('HIV incidence per 100'); title('Women'); %xlabel('Year'); 
xlim([2020 2060]); ylim([0 4]);
% legend('Clinic ART, mean' , 'range' , ...
%     'HTC + Community ART, mean' , 'range'); % , 'Scenario 3, mean' , 'range');
grid on;

% male
subplot(1 , 4 , 2);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_incidence_males_aged15-79' , '.csv'];
    hivIncM = xlsread(fname);
    if i == 1
        hold all;
        plot(hivIncM(: , 1) , hivIncM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivIncM(: , 1)' , fliplr(hivIncM(: , 1)')];
        inBetween = [hivIncM(: , 4)' , fliplr(hivIncM(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivIncM(: , 1) , hivIncM(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivIncM(: , 1) , hivIncM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    else
        hold all;
        plot(hivIncM((2019-startYear)+1 : end , 1) , hivIncM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivIncM(: , 1)' , fliplr(hivIncM(: , 1)')];
        inBetween = [hivIncM(: , 4)' , fliplr(hivIncM(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivIncM((2019-startYear)+1 : end , 1) , hivIncM((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivIncM((2019-startYear)+1 : end , 1) , hivIncM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    end
end
title('Men'); %ylabel('HIV incidence per 100'); xlabel('Year'); 
xlim([2020 2060]); ylim([0 4]);
% legend('Clinic ART, mean' , 'range' , ...
%     'HTC + Community ART, mean' , 'range'); % , 'Scenario 3, mean' , 'range');
grid on;

% combined
subplot(1 , 4 , 3);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_incidence_combined_aged15-79' , '.csv'];
    hivIncC = xlsread(fname);
    if i == 1
        hold all;
        plot(hivIncC(: , 1) , hivIncC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivIncC(: , 1)' , fliplr(hivIncC(: , 1)')];
        inBetween = [hivIncC(: , 4)' , fliplr(hivIncC(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivIncC(: , 1) , hivIncC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivIncC(: , 1) , hivIncC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    else
        hold all;
        plot(hivIncC((2019-startYear)+1 : end , 1) , hivIncC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivIncC(: , 1)' , fliplr(hivIncC(: , 1)')];
        inBetween = [hivIncC(: , 4)' , fliplr(hivIncC(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivIncC((2019-startYear)+1 : end , 1) , hivIncC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivIncC((2019-startYear)+1 : end , 1) , hivIncC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    end
end
title('Total Population'); %ylabel('HIV incidence per 100'); xlabel('Year'); 
xlim([2020 2060]); ylim([0 4]);
legend('Clinic ART, mean' , 'range' , ...
    'HTC + Community ART, mean' , 'range'); % , 'Scenario 3, mean' , 'range');
grid on;
%sgtitle('HIV Incidence');

%% Plot HIV mortality
figure('DefaultAxesFontSize' , 18);

% female
subplot(1 , 4 , 1);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_mortality_females_aged15-79' , '.csv'];
    hivMortF = xlsread(fname);
    if i == 1
        hold all;
        plot(hivMortF(: , 1) , hivMortF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivMortF(: , 1)' , fliplr(hivMortF(: , 1)')];
        inBetween = [hivMortF(: , 4)' , fliplr(hivMortF(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivMortF(: , 1) , hivMortF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivMortF(: , 1) , hivMortF(: , 4) , 'LineStyle' ,  '--' , 'Color' , colorList{iInd});
    else
        hold all;
        plot(hivMortF((2019-startYear)+1 : end , 1) , hivMortF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivMortF((2019-startYear)+1 : end , 1)' , fliplr(hivMortF((2019-startYear)+1 : end , 1)')];
        inBetween = [hivMortF((2019-startYear)+1 : end , 4)' , fliplr(hivMortF((2019-startYear)+1 : end , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivMortF((2019-startYear)+1 : end , 1) , hivMortF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivMortF((2019-startYear)+1 : end , 1) , hivMortF((2019-startYear)+1 : end , 4) , 'LineStyle' ,  '--' , 'Color' , colorList{iInd});
    end
end
ylabel('HIV mortality per 100K'); %title('Women'); xlabel('Year'); 
xlim([2020 2060]); ylim([0 2000]);
% legend('Scenario 1, mean' , 'range' , ...
%     'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range');
grid on;

% male
subplot(1 , 4 , 2);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_mortality_males_aged15-79' , '.csv'];
    hivMortM = xlsread(fname);
    if i == 1
        hold all;
        plot(hivMortM(: , 1) , hivMortM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivMortM(: , 1)' , fliplr(hivMortM(: , 1)')];
        inBetween = [hivMortM(: , 4)' , fliplr(hivMortM(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivMortM(: , 1) , hivMortM(: , 3) , 'LineStyle' ,  '--' , 'Color' , colorList{iInd});
        %plot(hivMortM(: , 1) , hivMortM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    else
        hold all;
        plot(hivMortM((2019-startYear)+1 : end , 1) , hivMortM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivMortM((2019-startYear)+1 : end , 1)' , fliplr(hivMortM((2019-startYear)+1 : end , 1)')];
        inBetween = [hivMortM((2019-startYear)+1 : end , 4)' , fliplr(hivMortM((2019-startYear)+1 : end , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivMortM((2019-startYear)+1 : end , 1) , hivMortM((2019-startYear)+1 : end , 3) , 'LineStyle' ,  '--' , 'Color' , colorList{iInd});
        %plot(hivMortM((2019-startYear)+1 : end , 1) , hivMortM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    end
end
%title('Men'); ylabel('Mortality per 100K'); xlabel('Year'); 
xlim([2020 2060]); ylim([0 2000]);
% legend('Scenario 1, mean' , 'range' , ...
%     'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range');
grid on;

% combined
subplot(1 , 4 , 3);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_mortality_combined_aged15-79' , '.csv'];
    hivMortC = xlsread(fname);
    if i == 1
        hold all;
        plot(hivMortC(: , 1) , hivMortC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivMortC(: , 1)' , fliplr(hivMortC(: , 1)')];
        inBetween = [hivMortC(: , 4)' , fliplr(hivMortC(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivMortC(: , 1) , hivMortC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivMortC(: , 1) , hivMortC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    else
        hold all;
        plot(hivMortC((2019-startYear)+1 : end , 1) , hivMortC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivMortC((2019-startYear)+1 : end , 1)' , fliplr(hivMortC((2019-startYear)+1 : end , 1)')];
        inBetween = [hivMortC((2019-startYear)+1 : end , 4)' , fliplr(hivMortC((2019-startYear)+1 : end , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivMortC((2019-startYear)+1 : end , 1) , hivMortC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivMortC((2019-startYear)+1 : end , 1) , hivMortC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    end
end
%xlabel('Year'); ylabel('Mortality per 100K'); title('Total Population');
xlim([2020 2060]); ylim([0 2000]);
% legend('Scenario 1, mean' , 'range' , ...
%     'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range');
grid on;
% sgtitle('HIV-Associated Mortality');

%% Plot HIV prevalence overall
figure('DefaultAxesFontSize' , 18);

% female
subplot(1 , 4 , 1);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_prevalence_females_aged15-79_all_HIV' , '.csv'];
    hivPrevF = xlsread(fname);
    if i == 1
        hold all;
        plot(hivPrevF(: , 1) , hivPrevF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivPrevF(: , 1)' , fliplr(hivPrevF(: , 1)')];
        inBetween = [hivPrevF(: , 4)' , fliplr(hivPrevF(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivPrevF(: , 1) , hivPrevF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivPrevF(: , 1) , hivPrevF(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    else
        hold all;
        plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivPrevF((2019-startYear)+1 : end , 1)' , fliplr(hivPrevF((2019-startYear)+1 : end , 1)')];
        inBetween = [hivPrevF((2019-startYear)+1 : end , 4)' , fliplr(hivPrevF((2019-startYear)+1 : end , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    end
end
xlabel('Year'); ylabel('HIV prevalence'); %title('Women'); 
xlim([2020 2060]); ylim([0.0 0.5]);
% legend('Scenario 1, mean' , 'range' , ...
%     'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range');
grid on;

% male
subplot(1 , 4 , 2);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_prevalence_males_aged15-79_all_HIV' , '.csv'];
    hivPrevM = xlsread(fname);
    if i == 1
        hold all;
        plot(hivPrevM(: , 1) , hivPrevM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivPrevM(: , 1)' , fliplr(hivPrevM(: , 1)')];
        inBetween = [hivPrevM(: , 4)' , fliplr(hivPrevM(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivPrevM(: , 1) , hivPrevM(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivPrevM(: , 1) , hivPrevM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    else
        hold all;
        plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivPrevM(: , 1)' , fliplr(hivPrevM(: , 1)')];
        inBetween = [hivPrevM(: , 4)' , fliplr(hivPrevM(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    end
end
xlabel('Year'); %title('Men'); %ylabel('HIV prevalence'); %
xlim([2020 2060]); ylim([0.0 0.5]);
% legend('Scenario 1, mean' , 'range' , ...
%     'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range');
grid on;

% combined
subplot(1 , 4 , 3);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_prevalence_combined_aged15-79_all_HIV' , '.csv'];
    hivPrevC = xlsread(fname);
    if i == 1
        hold all;
        plot(hivPrevC(: , 1) , hivPrevC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivPrevC(: , 1)' , fliplr(hivPrevC(: , 1)')];
        inBetween = [hivPrevC(: , 4)' , fliplr(hivPrevC(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivPrevC(: , 1) , hivPrevC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivPrevC(: , 1) , hivPrevC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    else
        hold all;
        plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [hivPrevC(: , 1)' , fliplr(hivPrevC(: , 1)')];
        inBetween = [hivPrevC(: , 4)' , fliplr(hivPrevC(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
        %plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        %plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    end
end
xlabel('Year'); %ylabel('HIV prevalence'); %title('Total Population'); 
xlim([2020 2060]); ylim([0.0 0.5]);
% legend('Scenario 1, mean' , 'range' , ...
%     'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range');
grid on;
% sgtitle('HIV Prevalence');

%% Plot HIV cases
figure('DefaultAxesFontSize' , 18);

% female
subplot(1 , 4 , 1);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Raw_HIV_incidence_females_ages15-79' , '.csv'];
    hivIncF = xlsread(fname);
    hivIncF = hivIncF((2020-startYear)+1 : end , :);
    hivIncF(: , 2:end) = cumsum(hivIncF(: , 2:end));
    hold all;
    plot(hivIncF(: , 1) , hivIncF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
    hold all;
    x2 = [hivIncF(: , 1)' , fliplr(hivIncF(: , 1)')];
    inBetween = [hivIncF(: , 4)' , fliplr(hivIncF(: , 3)')];
    h = fill(x2 , inBetween , 'k');
    h.FaceColor = colorList{iInd};
    h.EdgeColor = colorList{iInd};
    h.FaceAlpha = 0.3;
    h.LineStyle = '--';
    %plot(hivIncF(: , 1) , hivIncF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    %plot(hivIncF(: , 1) , hivIncF(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
end
ylabel('Cumulative new HIV cases'); title('Women'); %xlabel('Year'); 
xlim([2020 2060]); ylim([0 4*10^6]);
% legend('Clinic ART, mean' , 'range' , ...
%     'HTC + Community ART, mean' , 'range'); % , 'Scenario 3, mean' , 'range');
grid on;

% male
subplot(1 , 4 , 2);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Raw_HIV_incidence_males_ages15-79' , '.csv'];
    hivIncM = xlsread(fname);
    hivIncM = hivIncM((2020-startYear)+1 : end , :);
    hivIncM(: , 2:end) = cumsum(hivIncM(: , 2:end));
    hold all;
    plot(hivIncM(: , 1) , hivIncM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
    hold all;
    x2 = [hivIncM(: , 1)' , fliplr(hivIncM(: , 1)')];
    inBetween = [hivIncM(: , 4)' , fliplr(hivIncM(: , 3)')];
    h = fill(x2 , inBetween , 'k');
    h.FaceColor = colorList{iInd};
    h.EdgeColor = colorList{iInd};
    h.FaceAlpha = 0.2;
    h.LineStyle = '--';
    %plot(hivIncM(: , 1) , hivIncM(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    %plot(hivIncM(: , 1) , hivIncM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
end
title('Men'); %ylabel('HIV incidence per 100'); xlabel('Year'); 
xlim([2020 2060]); ylim([0 4*10^6]);
% legend('Clinic ART, mean' , 'range' , ...
%     'HTC + Community ART, mean' , 'range'); % , 'Scenario 3, mean' , 'range');
grid on;

% combined
subplot(1 , 4 , 3);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Raw_HIV_incidence_combined_ages15-79' , '.csv'];
    hivIncC = xlsread(fname);
    hivIncC = hivIncC((2020-startYear)+1 : end , :);
    hivIncC(: , 2:end) = cumsum(hivIncC(: , 2:end));
    hold all;
    plot(hivIncC(: , 1) , hivIncC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
    hold all;
    x2 = [hivIncC(: , 1)' , fliplr(hivIncC(: , 1)')];
    inBetween = [hivIncC(: , 4)' , fliplr(hivIncC(: , 3)')];
    h = fill(x2 , inBetween , 'k');
    h.FaceColor = colorList{iInd};
    h.EdgeColor = colorList{iInd};
    h.FaceAlpha = 0.2;
    h.LineStyle = '--';
    %plot(hivIncC(: , 1) , hivIncC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    %plot(hivIncC(: , 1) , hivIncC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
end
title('Total Population'); %ylabel('HIV incidence per 100'); xlabel('Year'); 
xlim([2020 2060]); ylim([0 4*10^6]);
legend('Clinic ART, mean' , 'range' , ...
    'HTC + Community ART, mean' , 'range'); % , 'Scenario 3, mean' , 'range');
grid on;

%% Plot HIV mortality
figure('DefaultAxesFontSize' , 18);

% female
subplot(1 , 4 , 1);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Raw_HIV_mortality_females_ages15-79' , '.csv'];
    hivMortF = xlsread(fname);
    hivMortF = hivMortF((2020-startYear)+1 : end , :);
    hivMortF(: , 2:end) = cumsum(hivMortF(: , 2:end));
    hold all;
    plot(hivMortF(: , 1) , hivMortF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
    hold all;
    x2 = [hivMortF(: , 1)' , fliplr(hivMortF(: , 1)')];
    inBetween = [hivMortF(: , 4)' , fliplr(hivMortF(: , 3)')];
    h = fill(x2 , inBetween , 'k');
    h.FaceColor = colorList{iInd};
    h.EdgeColor = colorList{iInd};
    h.FaceAlpha = 0.2;
    h.LineStyle = '--';
    %plot(hivMortF(: , 1) , hivMortF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    %plot(hivMortF(: , 1) , hivMortF(: , 4) , 'LineStyle' ,  '--' , 'Color' , colorList{iInd});
end
ylabel('Cumulative new HIV deaths'); xlabel('Year'); %title('Women'); 
xlim([2020 2060]); ylim([0 4*10^6]);
% legend('Scenario 1, mean' , 'range' , ...
%     'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range');
grid on;

% male
subplot(1 , 4 , 2);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'raw_HIV_mortality_males_ages15-79' , '.csv'];
    hivMortM = xlsread(fname);
    hivMortM = hivMortM((2020-startYear)+1 : end , :);
    hivMortM(: , 2:end) = cumsum(hivMortM(: , 2:end));
    hold all;
    plot(hivMortM(: , 1) , hivMortM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
    hold all;
    x2 = [hivMortM(: , 1)' , fliplr(hivMortM(: , 1)')];
    inBetween = [hivMortM(: , 4)' , fliplr(hivMortM(: , 3)')];
    h = fill(x2 , inBetween , 'k');
    h.FaceColor = colorList{iInd};
    h.EdgeColor = colorList{iInd};
    h.FaceAlpha = 0.2;
    h.LineStyle = '--';
    %plot(hivMortM(: , 1) , hivMortM(: , 3) , 'LineStyle' ,  '--' , 'Color' , colorList{iInd});
    %plot(hivMortM(: , 1) , hivMortM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
end
xlabel('Year'); %title('Men'); ylabel('Mortality per 100K');  
xlim([2020 2060]); ylim([0 4*10^6]);
% legend('Scenario 1, mean' , 'range' , ...
%     'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range');
grid on;

% combined
subplot(1 , 4 , 3);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Raw_HIV_mortality_combined_ages15-79' , '.csv'];
    hivMortC = xlsread(fname);
    hivMortC = hivMortC((2020-startYear)+1 : end , :);
    hivMortC(: , 2:end) = cumsum(hivMortC(: , 2:end));
    hold all;
    plot(hivMortC(: , 1) , hivMortC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
    hold all;
    x2 = [hivMortC(: , 1)' , fliplr(hivMortC(: , 1)')];
    inBetween = [hivMortC(: , 4)' , fliplr(hivMortC(: , 3)')];
    h = fill(x2 , inBetween , 'k');
    h.FaceColor = colorList{iInd};
    h.EdgeColor = colorList{iInd};
    h.FaceAlpha = 0.2;
    h.LineStyle = '--';
    %plot(hivMortC(: , 1) , hivMortC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    %plot(hivMortC(: , 1) , hivMortC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
end
xlabel('Year'); %ylabel('Mortality per 100K'); title('Total Population');
xlim([2020 2060]); ylim([0 4*10^6]);
% legend('Scenario 1, mean' , 'range' , ...
%     'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range');
grid on;

%% Plot ART prevalence
figure('DefaultAxesFontSize' , 18);

% female
subplot(1 , 3 , 1);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Proportion_WLWHIV_onART_aged15-79.csv'];
    artPrevF = xlsread(fname);
    if i == 1
        hold all;
        plot(artPrevF(: , 1) , artPrevF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [artPrevF(: , 1)' , fliplr(artPrevF(: , 1)')];
        inBetween = [artPrevF(: , 4)' , fliplr(artPrevF(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
    else
        hold all;
        plot(artPrevF((2019-startYear)+1 : end , 1) , artPrevF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [artPrevF((2019-startYear)+1 : end , 1)' , fliplr(artPrevF((2019-startYear)+1 : end , 1)')];
        inBetween = [artPrevF((2019-startYear)+1 : end , 4)' , fliplr(artPrevF((2019-startYear)+1 : end , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
    end
end
xlabel('Year'); ylabel('Proportion WLWHIV VS on ART'); title('WLWHIV on ART, ages 15-79');
xlim([2000 2060]); ylim([0.0 1.0]);
legend('Scenario 1, mean' , 'range' , ...
    'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range' , 'Location' , 'SouthEast');
grid on;

% male
subplot(1 , 3 , 2);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Proportion_MLWHIV_onART_aged15-79.csv'];
    artPrevM = xlsread(fname);
    if i == 1
        hold all;
        plot(artPrevM(: , 1) , artPrevM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [artPrevM(: , 1)' , fliplr(artPrevM(: , 1)')];
        inBetween = [artPrevM(: , 4)' , fliplr(artPrevM(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
    else
        hold all;
        plot(artPrevM((2019-startYear)+1 : end , 1) , artPrevM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [artPrevM((2019-startYear)+1 : end , 1)' , fliplr(artPrevM((2019-startYear)+1 : end , 1)')];
        inBetween = [artPrevM((2019-startYear)+1 : end , 4)' , fliplr(artPrevM((2019-startYear)+1 : end , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
    end
end
xlabel('Year'); ylabel('Proportion MLWHIV VS on ART'); title('MLWHIV on ART, ages 15-79');
xlim([2000 2060]); ylim([0.0 1.0]);
legend('Scenario 1, mean' , 'range' , ...
    'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range' , 'Location' , 'SouthEast');
grid on;

% combined
subplot(1 , 3 , 3);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Proportion_PLWHIV_onART_aged15-79.csv'];
    artPrevC = xlsread(fname);
    if i == 1
        hold all;
        plot(artPrevC(: , 1) , artPrevC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [artPrevC(: , 1)' , fliplr(artPrevC(: , 1)')];
        inBetween = [artPrevC(: , 4)' , fliplr(artPrevC(: , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
    else
        hold all;
        plot(artPrevC((2019-startYear)+1 : end , 1) , artPrevC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        hold all;
        x2 = [artPrevC((2019-startYear)+1 : end , 1)' , fliplr(artPrevC((2019-startYear)+1 : end , 1)')];
        inBetween = [artPrevC((2019-startYear)+1 : end , 4)' , fliplr(artPrevC((2019-startYear)+1 : end , 3)')];
        h = fill(x2 , inBetween , 'k');
        h.FaceColor = colorList{iInd};
        h.EdgeColor = colorList{iInd};
        h.FaceAlpha = 0.2;
        h.LineStyle = '--';
    end
end
xlabel('Year'); ylabel('Proportion PLWHIV VS on ART'); title('PLWHIV on ART, ages 15-79');
xlim([2000 2060]); ylim([0.0 1.0]);
legend('Scenario 1, mean' , 'range' , ...
    'Scenario 2, mean' , 'range' , 'Scenario 3, mean' , 'range' , 'Location' , 'SouthEast');
grid on;
sgtitle('Coverage of PLWHIV on ART and virally suppressed');

%% Plot HIV prevalence by CD4
cInds = {3 : 8 , 3 : 5 , 6 , 7 , 8};
cTits = {'all_HIV' , 'CD4_350plus_noART' , 'CD4_200-350_noART' , 'CD4_below200_noART' , 'onART'};
cTitsPlot = {'all HIV' , 'CD4 350plus noART' , 'CD4 200-350 noART' , 'CD4 below200 noART' , 'onART'};

% female
figure;
for c = 1 : length(cInds)
    subplot(2 , 3 , c);
    for iInd = 1 : length(iIndList)
        i = iIndList{iInd};
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
            'HIV_prevalence_females_aged15-79_' , cTits{c} , '.csv'];
        hivPrevF = xlsread(fname);
        if i == 1
            hold all;
            plot(hivPrevF(: , 1) , hivPrevF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(hivPrevF(: , 1) , hivPrevF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(hivPrevF(: , 1) , hivPrevF(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        else
            hold all;
            plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        end
    end
    xlabel('Year'); ylabel('HIV prevalence'); title(cTitsPlot{c}); 
    xlim([1980 2060]); ylim([0.0 0.5]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2a, mean' , 'min' , 'max' , ...
        'Scenario 2, mean' , 'min' , 'max' , 'Scenario 3, mean' , 'min' , 'max');
    grid on;
end
sgtitle('HIV Prevalence: Females, aged 15-79');

% male
figure;
for c = 1 : length(cInds)
    subplot(2 , 3 , c);
    for iInd = 1 : length(iIndList)
        i = iIndList{iInd};
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
            'HIV_prevalence_males_aged15-79_' , cTits{c} , '.csv'];
        hivPrevM = xlsread(fname);
        if i == 1
            hold all;
            plot(hivPrevM(: , 1) , hivPrevM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(hivPrevM(: , 1) , hivPrevM(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(hivPrevM(: , 1) , hivPrevM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        else
            hold all;
            plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        end
    end
    xlabel('Year'); ylabel('HIV prevalence'); title(cTitsPlot{c}); 
    xlim([1980 2060]); ylim([0.0 0.3]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2a, mean' , 'min' , 'max' , ...
        'Scenario 2, mean' , 'min' , 'max' , 'Scenario 3, mean' , 'min' , 'max');
    grid on;
end   
sgtitle('HIV Prevalence: Males, aged 15-79');
    
% combined
figure;
for c = 1 : length(cInds)
    subplot(2 , 3 , c);
    for iInd = 1 : length(iIndList)
        i = iIndList{iInd};
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
            'HIV_prevalence_combined_aged15-79_' , cTits{c} , '.csv'];
        hivPrevC = xlsread(fname);
        if i == 1
            hold all;
            plot(hivPrevC(: , 1) , hivPrevC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(hivPrevC(: , 1) , hivPrevC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(hivPrevC(: , 1) , hivPrevC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        else
            hold all;
            plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        end
    end
    xlabel('Year'); ylabel('HIV prevalence'); title(cTitsPlot{c}); 
    xlim([1980 2060]); ylim([0.0 0.4]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2a, mean' , 'min' , 'max' , ...
        'Scenario 2, mean' , 'min' , 'max' , 'Scenario 3, mean' , 'min' , 'max');
    grid on;
end   
sgtitle('HIV Prevalence: Females + Males, aged 15-79');

%% Plot CD4 distribution at ART initiation
cInds = {3 : 5 , 6 , 7};
cTits = {'CD4_350plus' , 'CD4_200-350' , 'CD4_below200'};
cTitsPlot = {'CD4 350plus' , 'CD4 200-350' , 'CD4 below200'};

% female
figure;
for c = 1 : length(cInds)
    subplot(1 , 3 , c);
    for iInd = 1 : length(iIndList)
        i = iIndList{iInd};
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
            'CD4_dist_initART_females_aged15-79_' , cTits{c} , '.csv'];
        cd4DistF = xlsread(fname);
        if i == 1
            hold all;
            plot(cd4DistF(: , 1) , cd4DistF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(cd4DistF(: , 1) , cd4DistF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd})'
            plot(cd4DistF(: , 1) , cd4DistF(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        else
            hold all;
            plot(cd4DistF((2019-startYear)+1 : end , 1) , cd4DistF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(cd4DistF((2019-startYear)+1 : end , 1) , cd4DistF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(cd4DistF((2019-startYear)+1 : end , 1) , cd4DistF((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        end
    end
    xlabel('Year'); ylabel('Proportion of females initiating ART'); title(['Proportion CD4: ' , cTitsPlot{c}]);
    xlim([1980 2060]); ylim([0.0 1.0]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2a, mean' , 'min' , 'max' , ...
        'Scenario 2, mean' , 'Scenario 3, mean' , 'min' , 'max' , 'min' , 'max');
    grid on;
end
sgtitle('Proportion of persons initiating ART: Females, aged 15-79');
   
% male 
figure;
for c = 1 : length(cInds)
    subplot(1 , 3 , c);
    for iInd = 1 : length(iIndList)
        i = iIndList{iInd};
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
            'CD4_dist_initART_males_aged15-79_' , cTits{c} , '.csv'];
        cd4DistM = xlsread(fname);
        if i == 1 
            hold all;
            plot(cd4DistM(: , 1) , cd4DistM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(cd4DistM(: , 1) , cd4DistM(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(cd4DistM(: , 1) , cd4DistM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        else
            hold all;
            plot(cd4DistM((2019-startYear)+1 : end , 1) , cd4DistM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(cd4DistM((2019-startYear)+1 : end , 1) , cd4DistM((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(cd4DistM((2019-startYear)+1 : end , 1) , cd4DistM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        end
    end
    xlabel('Year'); ylabel('Proportion of males initiating ART'); title(['Proportion CD4: ' , cTitsPlot{c}]);
    xlim([1980 2060]); ylim([0.0 1.0]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2a, mean' , 'min' , 'max' , ...
        'Scenario 2, mean' , 'min' , 'max' , 'Scenario 3, mean' , 'min' , 'max');
    grid on;
end
sgtitle('Proportion of persons initiating ART: Males, aged 15-79');
   
% combined 
figure;
for c = 1 : length(cInds)
    subplot(1 , 3 , c);
    for iInd = 1 : length(iIndList)
        i = iIndList{iInd};
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
            'CD4_dist_initART_combined_aged15-79_' , cTits{c} , '.csv'];
        cd4DistC = xlsread(fname);
        if i == 1
            hold all;
            plot(cd4DistC(: , 1) , cd4DistC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(cd4DistC(: , 1) , cd4DistC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(cd4DistC(: , 1) , cd4DistC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        else
            hold all;
            plot(cd4DistC((2019-startYear)+1 : end , 1) , cd4DistC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
            plot(cd4DistC((2019-startYear)+1 : end , 1) , cd4DistC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
            plot(cd4DistC((2019-startYear)+1 : end , 1) , cd4DistC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        end
    end
    xlabel('Year'); ylabel('Proportion of persons initiating ART'); title(['Proportion CD4: ' , cTitsPlot{c}]);
    xlim([1980 2060]); ylim([0.0 1.0]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2a, mean' , 'min' , 'max' , ...
        'Scenario 2, mean' , 'min' , 'max' , 'Scenario 3, mean' , 'min' , 'max');
    grid on;
end
sgtitle('Proportion of persons initiating ART: Females + Males, aged 15-79');

%% Plot percent reduction in HIV incidence and mortality over time
fig = figure;
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0.5 0.5 0.5]]);
set(fig,'DefaultAxesFontSize' , 18);

% female
% subplot(1 , 3 , 1);
% H = bar([1] , [62.23 65.70-62.23] , 'stacked', 'FaceColor' ,'flat');
% H(1).CData = [0, 0.4470, 0.7410]; 
% H(2).CData = [0.8500, 0.3250, 0.0980];
% box on;
% set(gca,'XColor','none')
% xlim([0.5 1.5]); ylim([0 100]);
% ylabel('Percent of Women Living with HIV who are Virally Suppressed');
% legend('Clinic ART' , 'Clinic ART + Community-based ART');
% title('Impact of Community-based ART on Viral Suppression');

subplot(2 , 2 , [1,2]);
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_DoART' , '_S' , num2str(1) , '_' , fileInds{1} , '\' , ...
    'HIV_incidence_females_aged15-79' , '.csv'];
hivIncFbase = xlsread(fname);
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_DoART' , '_S' , num2str(1) , '_' , fileInds{1} , '\' , ...
    'HIV_mortality_females_aged15-79' , '.csv'];
hivMortFbase = xlsread(fname);
timeInd = (currYear-startYear)+1;
for iInd = 2
    i = 2;
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_diagHiv075_DoART' , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_incidence_females_aged15-79' , '.csv'];
    hivIncF = xlsread(fname);
    hold all;
    %yyaxis left
    percReduct = ((hivIncF(: , 5:29) - hivIncFbase(: , 5:29)) ./ hivIncFbase(: , 5:29)) .* 100;
    percReduct_CI = [mean(percReduct , 2) , min(percReduct , [] , 2) , max(percReduct , [] , 2)];
    disp('Incidence reduction female (2025, 2060):')
    percReduct_CI((2025-1925)+1 , :)
    percReduct_CI(end , :)
    plot(hivIncF(timeInd:end , 1) , percReduct_CI(timeInd:end , 1) , 'LineStyle' , '-' , 'Color' , 'k' , 'LineWidth' , 2);
    x2 = [hivIncF(timeInd:end , 1)' , fliplr(hivIncF(timeInd:end , 1)')];
    inBetween = [percReduct_CI(timeInd:end , 2)' , fliplr(percReduct_CI(timeInd:end , 3)')];
    h = fill(x2 , inBetween , 'k');
    h.FaceColor = 'k';
    h.EdgeColor = 'k';
    h.FaceAlpha = 0.3;
    h.LineStyle = '-';
    h.LineWidth = 0.5;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    ylabel('Percent reduction'); ylim([-55 0]); box on;
    
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_diagHiv075_DoART' , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_mortality_females_aged15-79' , '.csv'];
    hivMortF = xlsread(fname);
    hold all;
    %yyaxis right
    percReduct = ((hivMortF(: , 5:29) - hivMortFbase(: , 5:29)) ./ hivMortFbase(: , 5:29)) .* 100;
    percReduct_CI = [mean(percReduct , 2) , min(percReduct , [] , 2) , max(percReduct , [] , 2)];
    disp('Mortality reduction female (2025, 2060):')
    percReduct_CI((2025-1925)+1 , :)
    percReduct_CI(end , :)
    plot(hivMortF(: , 1) , percReduct_CI(: , 1) , 'LineStyle' , '-' , 'Color' , [0.5 0.5 0.5] , 'LineWidth' , 2);
    x2 = [hivMortF(timeInd:end , 1)' , fliplr(hivMortF(timeInd:end , 1)')];
    inBetween = [percReduct_CI(timeInd:end , 2)' , fliplr(percReduct_CI(timeInd:end , 3)')];
    h = fill(x2 , inBetween , 'k');
    h.FaceColor = [0.5 0.5 0.5];
    h.EdgeColor = [0.5 0.5 0.5];
    h.FaceAlpha = 0.2;
    h.LineStyle = '-';
    h.LineWidth = 0.5;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    ylim([-60 0]); %ylabel('Percent reduction'); 
end
ax.YAxis(1).Color = [0 0 0];
ax.YAxis(2).Color = [0.5 0.5 0.5];
title('Women'); %xlabel('Year');  
hlegend1 = legend('Incidence, mean [uncertainty range]' , 'Mortality, mean [uncertainty range]' , 'Location' , 'NorthEast'); % , 'Orientation' , 'horizontal');
%hlegend1.NumColumns=2;
xlim([2020 2060]); grid on;
%sgtitle('FEMALES');

% fig = figure;
% set(fig,'defaultAxesColorOrder',[[0 0 0]; [0.5 0.5 0.5]]);
% set(fig,'DefaultAxesFontSize' , 18);

% male
% subplot(1 , 3 , 1);
% H = bar([1] , [39.78 65.50-39.78] , 'stacked' , 'FaceColor' ,'flat');
% H(1).CData = [0, 0.4470, 0.7410]; 
% H(2).CData = [0.8500, 0.3250, 0.0980];
% set(gca,'XColor','none')
% xlim([0.5 1.5]); ylim([0 100]); box on;
% ylabel('Percent of Men Living with HIV who are Virally Suppressed');
% legend('Clinic ART' , 'Clinic ART + Community-based ART');
% title('Impact of Community-based ART on Viral Suppression');

subplot(2 , 2 , [3,4]);
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_DoART' , '_S' , num2str(1) , '_' , fileInds{1} , '\' , ...
    'HIV_incidence_males_aged15-79' , '.csv'];
hivIncMbase = xlsread(fname);
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_DoART' , '_S' , num2str(1) , '_' , fileInds{1} , '\' , ...
    'HIV_mortality_males_aged15-79' , '.csv'];
hivMortMbase = xlsread(fname);
timeInd = (currYear-startYear)+1;
for iInd = 2
    i = 2;
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_diagHiv075_DoART' , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_incidence_males_aged15-79' , '.csv'];
    hivIncM = xlsread(fname);
    hold all;
    %yyaxis left
    percReduct = ((hivIncM(: , 5:29) - hivIncMbase(: , 5:29)) ./ hivIncMbase(: , 5:29)) .* 100;
    percReduct_CI = [mean(percReduct , 2) , min(percReduct , [] , 2) , max(percReduct , [] , 2)];
    disp('Incidence reduction male (2025, 2060):')
    percReduct_CI((2025-1925)+1 , :)
    percReduct_CI(end , :)
    plot(hivIncM(: , 1) , percReduct_CI(: , 1) , 'LineStyle' , '-' , 'Color' , 'k' , 'LineWidth' , 2);
    x2 = [hivIncM(timeInd:end , 1)' , fliplr(hivIncM(timeInd:end , 1)')];
    inBetween = [percReduct_CI(timeInd:end , 2)' , fliplr(percReduct_CI(timeInd:end , 3)')];
    h = fill(x2 , inBetween , 'k');
    h.FaceColor = 'k';
    h.EdgeColor = 'k';
    h.FaceAlpha = 0.3;
    h.LineStyle = '-';
    h.LineWidth = 0.5;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    ylabel('Percent reduction'); ylim([-55 0]); box on;
    
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_diagHiv075_DoART' , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_mortality_males_aged15-79' , '.csv'];
    hivMortM = xlsread(fname);
    hold all;
    %yyaxis right
    percReduct = ((hivMortM(: , 5:29) - hivMortMbase(: , 5:29)) ./ hivMortMbase(: , 5:29)) .* 100;
    percReduct_CI = [mean(percReduct , 2) , min(percReduct , [] , 2) , max(percReduct , [] , 2)];
    disp('Mortality reduction male (2025, 2060):')
    percReduct_CI((2025-1925)+1 , :)
    percReduct_CI(end , :)
    plot(hivMortM(: , 1) , percReduct_CI(: , 1) , 'LineStyle' , '-' , 'Color' , [0.5 0.5 0.5] , 'LineWidth' , 2);
    x2 = [hivMortM(timeInd:end , 1)' , fliplr(hivMortM(timeInd:end , 1)')];
    inBetween = [percReduct_CI(timeInd:end , 2)' , fliplr(percReduct_CI(timeInd:end , 3)')];
    h = fill(x2 , inBetween , 'k');
    h.FaceColor = [0.5 0.5 0.5];
    h.EdgeColor = [0.5 0.5 0.5];
    h.FaceAlpha = 0.2;
    h.LineStyle = '-';
    h.LineWidth = 0.5;
    h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    ylim([-60 0]); %ylabel('HIV mortality reduction (%)'); 
end
xlabel('Year');  title('Men');
hlegend2 = legend('Incidence, mean [uncertainty range]' , 'Mortality, mean [uncertainty range]' , 'Location' , 'NorthEast'); % , 'Orientation' , 'horizontal');
%hlegend2.NumColumns=2;
xlim([2020 2060]); grid on;
%sgtitle('MALES');

%% Plot gender difference in HIV incidence over time
fig = figure;

iIndListGenDiff = {1 , 2};
sceFileNameListGenDiff = {'_DoART' , '_diagHiv075_DoART'};
for iInd = 1 : length(iIndListGenDiff)
    i = iIndListGenDiff{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameListGenDiff{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_incidence_females_aged15-79' , '.csv'];
    hivIncF = xlsread(fname);
       
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameListGenDiff{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_incidence_males_aged15-79' , '.csv'];
    hivIncM = xlsread(fname);
        
    hold all;
    genderGap = (hivIncF(: , 5:29) ./ hivIncM(: , 5:29));
    genderGap_CI = [mean(genderGap , 2) , min(genderGap , [] , 2) , max(genderGap , [] , 2)];
    disp(['Gender gap (2025, 2060) for S' , num2str(i) , ':'])
    genderGap_CI((2025-1925)+1 , :)
    genderGap_CI(end , :)
    plot(hivIncF(: , 1) , genderGap_CI(: , 1) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
    plot(hivIncF(: , 1) , genderGap_CI(: , 2) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    plot(hivIncF(: , 1) , genderGap_CI(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    ylabel('Gender gap');
end
xlabel('Year');  title('Gender Gap in HIV Incidence');
xlim([2019 2060]); ylim([0 5]); grid on;

%% Plot gender difference in HIV prevalence over time
fig = figure;

iIndListGenDiff = {1 , 2};
sceFileNameListGenDiff = {'_DoART' , '_diagHiv075_DoART'};
for iInd = 1 : length(iIndListGenDiff)
    i = iIndListGenDiff{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameListGenDiff{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_prevalence_females_aged15-79_all_HIV' , '.csv'];
    hivIncF = xlsread(fname);
       
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameListGenDiff{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'HIV_prevalence_males_aged15-79_all_HIV' , '.csv'];
    hivIncM = xlsread(fname);
        
    hold all;
    genderGap = (hivIncF(: , 5:29) ./ hivIncM(: , 5:29));
    genderGap_CI = [mean(genderGap , 2) , min(genderGap , [] , 2) , max(genderGap , [] , 2)];
    disp(['Gender gap (2020 , 2025, 2060) for S' , num2str(i) , ':'])
    genderGap_CI((2020-1925)+1 , :)
    genderGap_CI((2025-1925)+1 , :)
    genderGap_CI(end , :)
    plot(hivIncF(: , 1) , genderGap_CI(: , 1) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
    plot(hivIncF(: , 1) , genderGap_CI(: , 2) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    plot(hivIncF(: , 1) , genderGap_CI(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    ylabel('Gender gap');
end
xlabel('Year');  title('Gender Gap in HIV Prevalence');
xlim([2019 2060]); grid on; %ylim([0 0.5]);

%% Output cumulative HIV cases and deaths
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_DoART' , '_S' , num2str(1) , '_' , fileInds{1} , '\' , ...
    'Raw_HIV_incidence_combined_ages15-79' , '.csv'];
hivIncCbase = xlsread(fname);
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_DoART' , '_S' , num2str(1) , '_' , fileInds{1} , '\' , ...
    'Raw_HIV_mortality_combined_ages15-79' , '.csv'];
hivMortCbase = xlsread(fname);
for iInd = 2
    i = 2;
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_diagHiv075_DoART' , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Raw_HIV_incidence_combined_ages15-79' , '.csv'];
    hivIncC = xlsread(fname);
    casesAverted = sum(hivIncCbase(: , 5:29) , 1) - sum(hivIncC(: , 5:29) , 1);
    casesAvertedPer = (sum(hivIncCbase(96:end , 5:29) , 1) - sum(hivIncC(96:end , 5:29) , 1)) ./ sum(hivIncCbase(96:end , 5:29) , 1) .* 100; 
    casesAverted_CI = [mean(casesAverted , 2) , min(casesAverted , [] , 2) , max(casesAverted , [] , 2)];
    casesAvertedPer_CI = [mean(casesAvertedPer , 2) , min(casesAvertedPer , [] , 2) , max(casesAvertedPer , [] , 2)];
    disp('Cases averted by 2060:')
    num2str(casesAverted_CI(end , :))
    num2str(casesAvertedPer_CI(end,:))
    casesAverted = sum(hivIncCbase(1:((2025-1925)+1) , 5:29) , 1) - sum(hivIncC(1:((2025-1925)+1) , 5:29) , 1);
    casesAverted_CI = [mean(casesAverted , 2) , min(casesAverted , [] , 2) , max(casesAverted , [] , 2)];
    disp('Cases averted by 2025:')
    num2str(casesAverted_CI(end , :))
    
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_diagHiv075_DoART' , '_S' , num2str(i) , '_', fileInds{1} , '\' , ...
        'Raw_HIV_mortality_combined_ages15-79' , '.csv'];
    hivMortC = xlsread(fname);
    deathsAverted = sum(hivMortCbase(: , 5:29) , 1) - sum(hivMortC(: , 5:29) , 1);
    deathsAvertedPer = (sum(hivMortCbase(96:end , 5:29) , 1) - sum(hivMortC(96:end , 5:29) , 1)) ./ sum(hivMortCbase(96:end , 5:29) , 1) .* 100;
    deathsAverted_CI = [mean(deathsAverted , 2) , min(deathsAverted , [] , 2) , max(deathsAverted , [] , 2)];
    deathsAvertedPer_CI = [mean(deathsAvertedPer , 2) , min(deathsAvertedPer , [] , 2) , max(deathsAvertedPer , [] , 2)];
    disp('Deaths averted by 2060:')
    num2str(deathsAverted_CI(end , :))
    num2str(deathsAvertedPer_CI(end,:))
    deathsAverted = sum(hivMortCbase(1:((2025-1925)+1) , 5:29) , 1) - sum(hivMortC(1:((2025-1925)+1) , 5:29) , 1);
    deathsAverted_CI = [mean(deathsAverted , 2) , min(deathsAverted , [] , 2) , max(deathsAverted , [] , 2)];
    disp('Deaths averted by 2025:')
    num2str(deathsAverted_CI(end , :))
end

%% Output additional circumcisions
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , '22Apr20Ph2V11_baseVax057_baseScreen_baseVMMC_fertDec042-076-052_2020ARTfxd_trackCD4-Discont_discontFxd' , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Raw_VMMC_male_ages15-79.csv'];
    vmmcMbase = xlsread(fname);
    
    fname = [pwd , '\HHCoM_Results\Vaccine' , '22Apr20Ph2V11_baseVax057_baseScreen_scaleVMMC_fertDec042-076-052_2020ARTfxd_trackCD4-Discont_discontFxd' , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'Raw_VMMC_male_ages15-79.csv'];
    vmmcMscale = xlsread(fname);
    
    addCircs = sum(vmmcMscale(: , 5:29) , 1) - sum(vmmcMbase(: , 5:29) , 1);
    addCircs_CI = [mean(addCircs , 2) , min(addCircs , [] , 2) , max(addCircs , [] , 2)];
    disp(['Additional circumcisions by 2060 for scenario' , num2str(i) , ':'])
    num2str(addCircs_CI(end , :))
end

%% Plot population size
figure;

% combined
%subplot(1 ,3 , 3);
for iInd = 1 : length(iIndList)
    i = iIndList{iInd};
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , sceFileNameList{iInd} , '_S' , num2str(i) , '_' , fileInds{1} , '\' , ...
        'PopulationSize_combined_aged15-79' , '.csv'];
    popSizeC = xlsread(fname);
    if i == 1
        hold all;
        plot(popSizeC(: , 1) , popSizeC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        plot(popSizeC(: , 1) , popSizeC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        plot(popSizeC(: , 1) , popSizeC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    else
        hold all;
        plot(popSizeC((2019-startYear)+1 : end , 1) , popSizeC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{iInd} , 'LineWidth' , 2);
        plot(popSizeC((2019-startYear)+1 : end , 1) , popSizeC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
        plot(popSizeC((2019-startYear)+1 : end , 1) , popSizeC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{iInd});
    end
end
xlabel('Year'); ylabel('Individuals'); title('Population Size: Females + Males, aged 15-79');
xlim([1980 2060]); ylim([0 14*(10^6)]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2a, mean' , 'min' , 'max' , ...
    'Scenario 2, mean' , 'min' , 'max' , 'Scenario 3, mean' , 'min' , 'max' , ...
    'Location' , 'Northwest');
grid on;
%sgtitle('Population size');
