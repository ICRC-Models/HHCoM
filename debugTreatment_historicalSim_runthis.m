% DESCRIPTION: Looping vaxCEA_multSims_SACEA.m over all the scenario
% numbers, cleaning the output, turning into an array, and exporting to CSV for 
% processing in R. 

function debugTreatment_historicalSim_runthis(username)

clear;

% Initialize variables
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
    c3c2Mults , c2c1Mults , c2c3Mults , c1c2Mults , muCC , muCC_ud, muCC_d, kRL , kDR , artHpvMult , ...
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
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue, ...
    ccLochpvVaxIndsFrom_treat , ...
    ccReghpvVaxInds_treat , ccDisthpvVaxInds_treat] = loadUp2(1 , 0 , [] , [] , []);

% Indices of calib runs to plot
% Temporarily commenting out to only run one scenario first to test out
% code
% fileInds = {'6_1' , '6_2' , '6_3' , '6_6' , '6_8' , '6_9' , '6_11' , ...
%      '6_12' , '6_13' , '6_15' , '6_20' , '6_21' , '6_22' , '6_26' , ...
%     '6_27' , '6_32' , '6_34' , '6_35' , '6_38' , '6_39' , '6_40' , ...
%     '6_41' , '6_42' , '6_45' , '6_47'};    % 22Apr20Ph2V11 ***************SET ME****************
% fileInds = {'6_1', '6_2'}; % FORTESTING
% nRuns = length(fileInds);
nRuns = 1; 
fileInds = {'6_1'}; 

lastYear = 2021; % manually set in futureSim 
monthlyTimespan = [startYear : timeStep : lastYear]; % list all the timespans in a vector
% monthlyTimespan = monthlyTimespan(1 : end-1); % remove the very last date
monthlyTimespanFut = [endYear : timeStep : lastYear]; % screening time span starts at 2021
monthlyTimespanFut = monthlyTimespanFut(1 : end-1); 
nTimepoints = length(monthlyTimespan);
nTimepointsFut = length(monthlyTimespanFut); 

% parallelizing the for loop
% loopSegments = {0 , round(10/2) , 10}; % running 10 scenarios
% loopSegmentsLength = length(loopSegments);

% for k = 1 : loopSegmentsLength-1 
%     parfor j = loopSegments{k}+1 : loopSegments{k+1}
    
% for j = [1 4] % FORTESTING
j=1; 
    sceNum = j - 1; 
    sceString = num2str(sceNum); % turn sceNum into string sceString
    sce = sceNum + 1; % add one since indices start at 1 (so scenarios will be 1-10 in this case) 

    % Initialize result matrices 
    deaths = zeros(nTimepoints, age, 4, 7, nRuns); % time, age (1:17), 3 death data elements, number of parameters, 10 scenarios
    ccHealthState = zeros(nTimepoints, age, 10, nRuns); 
    newCC = zeros(nTimepoints, disease, age, 2, nRuns); 
    symptomatic = zeros(nTimepoints, age, 3, 3, nRuns);     
    treatment = zeros(nTimepoints, age, 3, 3, nRuns); 
    screening = zeros(nTimepoints, age, 3, 2, nRuns); 

    % Feeding in the zeroed result matrix, spitting out the same matrix but with all the counts added in for that scenario
    [deaths, newCC, ccHealthState, symptomatic, treatment, screening] = ...
        debugTreatment_historicalSim(1 , sceString , {'0'}, fileInds, deaths, newCC, ccHealthState, symptomatic, ...
            treatment, screening);  

% turn all the result matrices into 2D 
    for param = 1 : nRuns
        for a = 1 : age 
            % turning deaths into 2D 
            for var = 1 : 4 
                for x = 1 : 7
                    if (param == 1 && a ==1 && var == 1)
                        deathsReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
                                        x.*ones(nTimepoints,1), ...
                                        param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), deaths(:, a, var, x, param)];
                    else
                        deathsReshape = [deathsReshape; ...
                                 transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
                                 x.*ones(nTimepoints,1) ...
                                 param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), deaths(:, a, var, x, param)];
                    end
                end  
            end 

            for var = 1 : 10
                if (param == 1 && a ==1 && var == 1)
                    ccHealthStateReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
                                             param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ccHealthState(:, a, var, param)];
                else 
                    ccHealthStateReshape = [ccHealthStateReshape; 
                                            transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
                                             param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ccHealthState(:, a, var, param)];
                end 
            end 

            for var = 1 : 3
                for x = 1 : 3
                    if (param == 1 && a == 1 && var == 1 && x == 1)
                        sympReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), ...
                                             var.*ones(nTimepoints,1), ...
                                             param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), symptomatic(:, a, x, var, param)];
                    else
                        sympReshape = [sympReshape; 
                                            transpose(monthlyTimespan), a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), ...
                                            var.*ones(nTimepoints,1), ...
                                            param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), symptomatic(:, a, x, var, param)];
                    end
                end
            end

            for d = 1 : disease 
                for rand = 1 : 2

                    if (param == 1 && a == 1 && d == 1 && rand == 1)

                        newCCReshape = [transpose(monthlyTimespan), d.*ones(nTimepoints,1), a.*ones(nTimepoints,1), ...
                                            rand.*ones(nTimepoints,1), ...
                                            param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
                                                newCC(:, d, a, rand)]; 
                    else 
                        newCCReshape = [newCCReshape; 
                                                transpose(monthlyTimespan), d.*ones(nTimepoints,1), a.*ones(nTimepoints,1), ...
                                                rand.*ones(nTimepoints,1), ...
                                                param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
                                                newCC(:, d, a, rand)]; 
                    end 
                end 
            end 

            for var = 1 : 3
                for x = 1 : 3
                    if (param == 1 && a == 1 && var == 1 && x == 1)
                        treatReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), ...
                                             var.*ones(nTimepoints,1), ...
                                             param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), treatment(:, a, x, var, param)];
                    else
                        treatReshape = [treatReshape; 
                                            transpose(monthlyTimespan), a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), ...
                                            var.*ones(nTimepoints,1), ...
                                            param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), treatment(:, a, x, var, param)];
                    end
                end
            end     

            for var = 1 : 2
                for x = 1 : 3
                    if (param == 1 && a == 1 && var == 1 && x == 1)
                        screenReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), ...
                                             var.*ones(nTimepoints,1), ...
                                             param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), screening(:, a, x, var, param)];
                    else
                        screenReshape = [screenReshape; 
                                            transpose(monthlyTimespan), a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), ...
                                            var.*ones(nTimepoints,1), ...
                                            param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), screening(:, a, x, var, param)];
                    end
                end
            end

        end
        disp(['Complete Scenario ', num2str(sce), ', Parameter ', num2str(param)])
    end

% turn into arrays
deathsReshape1 = array2table(deathsReshape, 'VariableNames', {'year', 'age', 'categ', 'stage', 'paramNum', ...
        'sceNum', 'count'}); 
ccHealthStateReshape1 = array2table(ccHealthStateReshape, 'VariableNames', {'year', 'age', 'categ', 'paramNum', ...
        'sceNum', 'count'}); 
newCCReshape1 = array2table(newCCReshape, 'VariableNames', {'year', 'disease', 'age', 'vaxtype', 'paramNum', 'sceNum', 'count'}); 
sympReshape1 = array2table(sympReshape, 'VariableNames', {'year', 'age', 'stage', 'categ', 'paramNum', ...
        'sceNum', 'count'}); 
treatReshape1 = array2table(treatReshape, 'VariableNames', {'year', 'age', 'stage', 'categ', 'paramNum', ...
        'sceNum', 'count'}); 
screenReshape1 = array2table(screenReshape, 'VariableNames', {'year', 'age', 'stage', 'categ', 'paramNum', ...
        'sceNum', 'count'}); 

% spit out into CSV 
writetable(deathsReshape1,[pwd '/testing/deaths_S' sceString '.csv']);
writetable(ccHealthStateReshape1, [pwd '/testing/ccHealthState_S' sceString '.csv']);
writetable(newCCReshape1, [pwd '/testing/newCC_S' sceString '.csv']); 
writetable(sympReshape1, [pwd '/testing/symp_S' sceString '.csv']);
writetable(treatReshape1, [pwd '/testing/treat_S' sceString '.csv']);
writetable(screenReshape1, [pwd '/testing/screen_S' sceString '.csv']);