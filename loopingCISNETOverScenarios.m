% DESCRIPTION: Looping vaxCEA_multSims_CISNET_UWRuns.m over all the scenario
% numbers, cleaning the output, turning into an array, and exporting to CSV for 
% processing in R. 

function loopingCISNETOverScenarios(username)

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

% Indices of calib runs to plot
% Temporarily commenting out to only run one scenario first to test out
% code
fileInds = {'6_1' , '6_2' , '6_3' , '6_6' , '6_8' , '6_9' , '6_11' , ...
     '6_12' , '6_13' , '6_15' , '6_20' , '6_21' , '6_22' , '6_26' , ...
    '6_27' , '6_32' , '6_34' , '6_35' , '6_38' , '6_39' , '6_40' , ...
    '6_41' , '6_42' , '6_45' , '6_47'};    % 22Apr20Ph2V11 ***************SET ME****************
% fileInds = {'6_1', '6_2'}; % FORTESTING
nRuns = length(fileInds);

lastYear = 2122; % manually set in futureSim 
monthlyTimespan = [startYear : timeStep : lastYear]; % list all the timespans in a vector
monthlyTimespan = monthlyTimespan(1 : end-1); % remove the very last date
monthlyTimespanFut = [endYear : timeStep : lastYear]; % screening time span starts at 2021
monthlyTimespanFut = monthlyTimespanFut(1 : end-1); 
nTimepoints = length(monthlyTimespan);
nTimepointsFut = length(monthlyTimespanFut); 

% parallelizing the for loop
loopSegments = {0 , round(10/2) , 10}; % running 10 scenarios
loopSegmentsLength = length(loopSegments);

% for k = 1 : loopSegmentsLength-1 
%     parfor j = loopSegments{k}+1 : loopSegments{k+1}
    
% only run scenario 1
for j = [1] % scenario number 

    sceNum = j; 
    sceString = num2str(sceNum); % turn sceNum into string sceString
    sce = sceNum; % add one since indices start at 1 (so scenarios will be 1-10 in this case) 

    % Initialize result matrices 
    deaths = zeros(nTimepoints, age+1, 3, nRuns); % time, age (1:17), 3 death data elements, number of parameters, 10 scenarios
    screenTreat = zeros(nTimepointsFut, age+1, 5, nRuns); % time, age (1:16), 5 screen/treat data elements, number of parameters, 10 scenarios
    hpvHealthState = zeros(nTimepoints, age+1, 7, nRuns); % time, age (1:16), 10 HPV/CC health states, number of parameters, 10 scenarios
    ccHealthState = zeros(nTimepoints, age+1, 4, nRuns); 
    hivHealthState = zeros(nTimepoints, age+1, 7, nRuns); % time, age (1:16), 7 HIV health states, number of parameters , 10 scenarios
    totalPerAge = zeros(nTimepoints, age+1, nRuns); % to pull the N per age group at each time point , 10 scenarios
    vax = zeros(nTimepoints, age+1, nRuns); % number of vaccinations is not stratified by age. only time and parameter.    , 10 scenarios  
    nonDisabHealthState = zeros(nTimepoints, age+1, nRuns); 
    scen0CinTx = zeros(nTimepointsFut, endpoints, 2, nRuns); % 2 because only LEEP and cryo (for the 3rd dimension, 1 is LEEP 2 is cryo)
    newCC = zeros(nTimepoints, age+1, nRuns); 
    hivDeath = zeros(nTimepoints, nRuns); 
    womenCount = zeros(nTimepoints, nRuns);
    womenCountAge = zeros(nTimepoints, age, gender, nRuns); 
    hivPrev = zeros(nTimepoints, age, gender, nRuns); 
    hivPrevTotal = zeros(nTimepoints, nRuns); 

    % Feeding in the zeroed result matrix, spitting out the same matrix but with all the counts added in for that scenario
    [hivDeath, womenCount, hivPrev, hivPrevTotal, womenCountAge] = ...
        vaxCEA_multSims_CISNET_UWRuns(1 , sceString , {'0'}, fileInds, hivDeath, womenCount, hivPrev, hivPrevTotal, womenCountAge);  

% turn all the result matrices into 2D 
    for param = 1 : nRuns
        if (param == 1) 
            hivDeathReshape = [transpose(monthlyTimespan), param.*ones(nTimepoints, 1), sce.*ones(nTimepoints,1), hivDeath(:, param)]; 
            womenCountReshape = [transpose(monthlyTimespan), param.*ones(nTimepoints, 1), sce.*ones(nTimepoints,1), womenCount(:, param)];
        else 
            hivDeathReshape = [hivDeathReshape; ...
                        transpose(monthlyTimespan), param.*ones(nTimepoints, 1), sce.*ones(nTimepoints, 1), hivDeath(:, param)];
            womenCountReshape = [womenCountReshape; ...
                        transpose(monthlyTimespan), param.*ones(nTimepoints, 1), sce.*ones(nTimepoints,1), womenCount(:, param)];
        end

        for a = 1 : age
            for g = 1 : gender 
                if (param == 1 && a == 1 && g == 1)
                    hivPrevReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), g.*ones(nTimepoints,1), param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
                        hivPrev(:, a, g, param)]; 
                    womenCountAgeReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), g.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                        sce.*ones(nTimepoints,1), womenCountAge(:, a, g, param)]; 
                else 
                    hivPrevReshape = [hivPrevReshape; ...
                        transpose(monthlyTimespan), a.*ones(nTimepoints,1), g.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                        sce.*ones(nTimepoints,1), hivPrev(:, a, g, param)];
                    womenCountAgeReshape = [womenCountAgeReshape; ...
                        transpose(monthlyTimespan), a.*ones(nTimepoints,1), g.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                        sce.*ones(nTimepoints,1), womenCountAge(:, a, g, param)];
                end 
            end 
        end 

%         for a = 1 : age + 1 
%             % turning deaths into 2D 
%             for var = 1 : 3 
%                 if (param == 1 && a ==1 && var == 1)
%                     deathsReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
%                                         param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), deaths(:, a, var, param)];
%                 else
%                     deathsReshape = [deathsReshape; ...
%                                  transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
%                                  param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), deaths(:, a, var, param)];
%                 end 
%             end 
%             % turning screening/treatment matrix into 2D
%             for var = 1 : 5 
%                 if (param == 1 && a ==1 && var == 1)
%                     screenTreatReshape = [transpose(monthlyTimespanFut), a.*ones(nTimepointsFut,1), var.*ones(nTimepointsFut,1), ...
%                                              param.*ones(nTimepointsFut,1), sce.*ones(nTimepointsFut,1), screenTreat(:, a, var, param)];
%                 else
%                     screenTreatReshape = [screenTreatReshape; ...
%                                             transpose(monthlyTimespanFut), a.*ones(nTimepointsFut,1), var.*ones(nTimepointsFut,1), ...
%                                              param.*ones(nTimepointsFut,1), sce.*ones(nTimepointsFut,1), screenTreat(:, a, var, param)];
%                 end 
%             end 
%             % turning hpv health states matrix into 2D
%             for var = 1 : 7
%                 if (param == 1 && a ==1 && var == 1)
%                     hpvHealthStateReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
%                                              param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), hpvHealthState(:, a, var, param)];
%                 else 
%                     hpvHealthStateReshape = [hpvHealthStateReshape; 
%                                             transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
%                                              param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), hpvHealthState(:, a, var, param)];
%                 end
%             end 
%             % turning cc health states matrix into 2D 
%             for var = 1 : 4
%                 if (param == 1 && a ==1 && var == 1)
%                     ccHealthStateReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
%                                              param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ccHealthState(:, a, var, param)];
%                 else 
%                     ccHealthStateReshape = [ccHealthStateReshape; 
%                                             transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
%                                              param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ccHealthState(:, a, var, param)];
%                 end 
%             end 
%             % turning hiv health states matrix into 2D 
%             for var = 1 : 7
%                 if (param == 1 && a ==1 && var == 1)
%                     hivHealthStateReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
%                                              param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), hivHealthState(:, a, var, param)];
%                 else 
%                     hivHealthStateReshape = [hivHealthStateReshape; 
%                                             transpose(monthlyTimespan), a.*ones(nTimepoints,1), var.*ones(nTimepoints,1), ...
%                                              param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), hivHealthState(:, a, var, param)];
%                 end
%             end 
%             % turning non disability health states matrix into 2D
%             if (param == 1 && a == 1)
%                 nonDisabHealthStateReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
%                     sce.*ones(nTimepoints,1), nonDisabHealthState(:, a, param)]; 
%             else 
%                 nonDisabHealthStateReshape = [nonDisabHealthStateReshape; ...
%                                                     transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
%                                                     sce.*ones(nTimepoints,1), nonDisabHealthState(:, a, param)]; 
%             end 
%             % turning total per age matrix into 2D 
%             if (param == 1 && a == 1)
%                 totalPerAgeReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
%                                         sce.*ones(nTimepoints,1), totalPerAge(:, a, param)]; 
%                 vaxReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
%                                         sce.*ones(nTimepoints,1), vax(:, a, param)];
%                 newCCReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
%                                         newCC(:, a, param)]; 
%             else 
%                 totalPerAgeReshape = [totalPerAgeReshape; 
%                                         transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
%                                         sce.*ones(nTimepoints,1), totalPerAge(:, a, param)]; 
%                 vaxReshape = [vaxReshape; 
%                                         transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
%                                         sce.*ones(nTimepoints,1), vax(:, a, param)]; 
%                 newCCReshape = [newCCReshape; 
%                                         transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
%                                         newCC(:, a, param)]; 
%             end 
%         end
        disp(['Complete Scenario ', num2str(sce), ', Parameter ', num2str(param)])
    end

% turn result matrices for scenario 0 cin analysis into 2D
%     for param = 1 : nRuns
%         for h = 1 : hpvVaxStates
%             % turning the cin treatment scenario 0 analysis into 2D 
%             for var = 1 : 2
%                 if (param == 1 && h ==1 && var == 1)
%                     scen0CinTxReshape = [transpose(monthlyTimespanFut), h.*ones(nTimepointsFut,1), var.*ones(nTimepointsFut,1), ...
%                                          param.*ones(nTimepointsFut,1), sce.*ones(nTimepointsFut,1), scen0CinTx(:, h, var, param)]; 
%                 else 
%                     scen0CinTxReshape = [scen0CinTxReshape; 
%                                             transpose(monthlyTimespanFut), h.*ones(nTimepointsFut,1), var.*ones(nTimepointsFut,1), ...
%                                             param.*ones(nTimepointsFut,1), sce.*ones(nTimepointsFut,1), scen0CinTx(:, h, var, param)]; 
%                 end 
%             end 
%         end
%     end 

% turn into arrays
% deathsReshape1 = array2table(deathsReshape, 'VariableNames', {'year', 'age', 'categ', 'paramNum', ...
%         'sceNum', 'count'}); 
% screenTreatReshape1 = array2table(screenTreatReshape, 'VariableNames', {'year', 'age', 'categ', 'paramNum', ...
%         'sceNum', 'count'}); 
% hpvHealthStateReshape1 = array2table(hpvHealthStateReshape, 'VariableNames', {'year', 'age', 'categ', 'paramNum', ...
%         'sceNum', 'count'}); 
% ccHealthStateReshape1 = array2table(ccHealthStateReshape, 'VariableNames', {'year', 'age', 'categ', 'paramNum', ...
%         'sceNum', 'count'}); 
% hivHealthStateReshape1 = array2table(hivHealthStateReshape, 'VariableNames', {'year', 'age', 'categ', 'paramNum', ...
%         'sceNum', 'count'}); 
% totalPerAgeReshape1 = array2table(totalPerAgeReshape, 'VariableNames', {'year', 'age', 'paramNum', ...
%         'sceNum', 'count'}); 
% vaxReshape1 = array2table(vaxReshape, 'VariableNames', {'year', 'age', 'paramNum', ...
%         'sceNum', 'count'}); 
% nonDisabHealthStateReshape1 = array2table(nonDisabHealthStateReshape, 'VariableNames', {'year', 'age', 'paramNum', ...
%         'sceNum', 'count'}); 
% scen0CinTxReshape1 = array2table(scen0CinTxReshape, 'VariableNames', {'year', 'hpvState', 'categ', 'paramNum', 'sceNum', 'count'}); 
% newCCReshape1 = array2table(newCCReshape, 'VariableNames', {'year', 'age', 'paramNum', 'sceNum', 'count'}); 
hivDeathReshape1 = array2table(hivDeathReshape, 'VariableNames', {'year', 'paramNum', 'sceNum', 'count'});
womenCountReshape1 = array2table(womenCountReshape, 'VariableNames', {'year', 'paramNum', 'sceNum', 'count'}); 
womenCountAgeReshape1 = array2table(womenCountAgeReshape, 'VariableNames', {'year', 'age', 'gender', 'paramNum', 'sceNum', 'count'}); 
hivPrevReshape1 = array2table(hivPrevReshape, 'VariableNames', {'year', 'age', 'gender', 'paramNum', 'sceNum', 'count'}); 


% spit out into CSV 
% writetable(deathsReshape1,[pwd '/SACEA/deaths_S' sceString '.csv']);
% writetable(screenTreatReshape1, [pwd '/SACEA/screenTreat_S' sceString '.csv']); 
% writetable(hpvHealthStateReshape1, [pwd '/SACEA/hpvHealthState_S' sceString '.csv']); 
% writetable(ccHealthStateReshape1, [pwd '/SACEA/ccHealthState_S' sceString '.csv']);
% writetable(hivHealthStateReshape1, [pwd '/SACEA/hivHealthState_S' sceString '.csv']); 
% writetable(totalPerAgeReshape1, [pwd '/SACEA/totalPerAge_S' sceString '.csv']);
% writetable(vaxReshape1, [pwd '/SACEA/vax_S' sceString '.csv']);
% writetable(nonDisabHealthStateReshape1, [pwd '/SACEA/nonDisabHealthState_S' sceString '.csv']);
% writetable(newCCReshape1, [pwd '/SACEA/newCC_S' sceString '.csv']); 
writetable(hivDeathReshape1, [pwd '/CISNET/hivDeath_S' sceString '.csv']); 
writetable(womenCountReshape1, [pwd '/CISNET/womenCount_S' sceString '.csv']); 
writetable(womenCountAgeReshape1, [pwd '/CISNET/womenCountAgeGender_S' sceString '.csv']); 
writetable(hivPrevReshape1, [pwd '/CISNET/hivPrevAgeGender_S' sceString '.csv']); 

% if sceNum == 0 % only write the cin / treatment file if in scenario 0 
%     writetable(scen0CinTxReshape1, [pwd '/SACEA/scen0CinTx_S0.csv']); 
% end 

    end
% end
end