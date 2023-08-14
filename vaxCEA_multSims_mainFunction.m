% DESCRIPTION: Looping processModelResults_IPVC2023 over all the scenario
% numbers, cleaning the output, turning into an array, and exporting to CSV for 
% processing in R. 

function vaxCEA_multSims_mainFunction(username)

clear;

% Initialize variables
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
    c3c2Mults , c2c1Mults , muCC , muCC_ud , muCC_d , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , waning , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear ,circNatStartYear , vaxStartYear , baseline , cisnet , who , whob , ...
    circProtect , condProtect , MTCTRate , hyst , OMEGA , ...
    ccInc2012_dObs , cc_dist_dObs , cin3_dist_dObs , cin1_dist_dObs , ...
    hpv_dist_dObs , cinPos2007_dObs , cin1_2010_dObs , cin2_2010_dObs, ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_all_dObs , hpv_hiv2009_dObs , ...
    hivPrevF_dObs , hivPrevM_dObs , hivPrevAll_dObs, popAgeDist_dObs , totPopSize_dObs , hivCurr , ...
    stageDist_2018_dObs , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ageInd , riskInd , ...
    kSymp , hystMult , ...
    hivNegNonVMMCinds , hivNegVMMCinds , vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3  , dFertMat3, d_partnersMmult, riskAdj, d_riskAdj, ...
    deathMat , deathMat2 , deathMat3 , deathMat4 , deathMat5,...
    dDeathMat , dDeathMat2 , dDeathMat3 , dDeathMat4, dMue , ...
    ccLochpvVaxIndsFrom_treat , ...
    ccReghpvVaxInds_treat , ccDisthpvVaxInds_treat , vaxEff] = loadUp2(1 , 0 , [] , [] , [] , 1);

% Indices of calib runs to plot
% Temporarily commenting out to only run one scenario first to test out
% code
%***SET ME***: this will likely change once we spit out all of the sim results
fileInds = {'1', '2', '3', '4', '5', '6', '7', '8', '9', '10', ...
                '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', ...
                '22', '23', '24', '25'};    % 22Apr20Ph2V11 ***************SET ME****************
% fileInds = {'1'}; % FORTESTING
nRuns = length(fileInds);

lastYear = 2123; % manually set in futureSim
monthlyTimespan = [startYear : timeStep : lastYear]; % list all the timespans in a vector
monthlyTimespan = monthlyTimespan(1 : end-1); % remove the very last date
monthlyTimespanFut = [endYear : timeStep : lastYear]; % screening time span starts at 2021
monthlyTimespanFut = monthlyTimespanFut(1 : end-1); 
nTimepoints = length(monthlyTimespan);
nTimepointsFut = length(monthlyTimespanFut); 
fivYrAgeGrpsOn = 1;

% parallelizing the for loop
loopSegments = {0 , round(2/2) , 5}; % running 10 scenarios ***SET ME***: the number of scenarios will be different
loopSegmentsLength = length(loopSegments);

% scenarios = {'1.1', '1.2', '2.1', '2.2', '3.1'}; ***SET ME***: specify the scenarios to loop through
scenarios = {'0', '1', '2', '3', '4'}; 

for k = 1 : loopSegmentsLength-1 
    parfor j = loopSegments{k}+1 : loopSegments{k+1} % for testing (parfor)
    
% for j = [1] % FORTESTING

    sceNum = j - 1; 
    sceString = scenarios{j}; % turn sceNum into string sceString
    sce = sceNum + 1; % add one since indices start at 1 (so scenarios will be 1-10 in this case) 

    % Initialize result matrices 
    vax = zeros(nTimepoints, age+1, 2, nRuns); % number of vaccinations is not stratified by age. only time and parameter.    , 10 scenarios  
    deaths = zeros(nTimepoints, age+1, 3, nRuns); 
    ccHealthState = zeros(nTimepoints, age+1, endpoints, nRuns); 
    hpvHealthState = zeros(nTimepoints, age+1, 7, nRuns); 
    newCC = zeros(nTimepoints, age+1, nRuns); 
    totalPerAge = zeros(nTimepoints, age+1, nRuns); 
    screenTreat = zeros(nTimepoints, age+1, 3, nRuns); 
    screenSympCCTreat = zeros(nTimepoints, 3, age+1, 3, 2, nRuns); 

    % Feeding in the zeroed result matrix, spitting out the same matrix but with all the counts added in for that scenario
    [vax, deaths, ccHealthState, hpvHealthState, newCC, totalPerAge, screenTreat, screenSympCCTreat] = ...
        vaxCEA_multSims_processResults(1 , sceString , {'0'}, fileInds, vax, deaths, ccHealthState, hpvHealthState, newCC, totalPerAge, screenTreat, screenSympCCTreat);  

% turn all the result matrices into 2D 
    for param = 1 : nRuns
        for a = 1 : (age + 1)
            % turning total per age matrix into 2D 
            if (param == 1 && a == 1)
                newCCReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
                                        newCC(:, a, param)]; 
                totalPerAgeReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
                                        totalPerAge(:, a, param)];
            else 
                newCCReshape = [newCCReshape; 
                                        transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
                                        newCC(:, a, param)]; 
                totalPerAgeReshape = [totalPerAgeReshape; 
                                        transpose(monthlyTimespan), a.*ones(nTimepoints,1), param.*ones(nTimepoints,1), sce.*ones(nTimepoints,1), ...
                                        totalPerAge(:, a, param)]; 
            end 

            for index = 1 : 2
                if (param == 1 && a == 1 && index == 1)
                    vaxReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), index.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                        sce.*ones(nTimepoints,1), vax(:, a, index, param)];
                else 
                    vaxReshape = [vaxReshape; 
                                        transpose(monthlyTimespan), a.*ones(nTimepoints,1), index.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                        sce.*ones(nTimepoints,1), vax(:, a, index, param)]; 
                end 
            end 

            for index = 1 : 3 
                if (param == 1 && a == 1 && index == 1)
                    deathsReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), index.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                        sce.*ones(nTimepoints,1), deaths(:, a, index, param)];
                else 
                    deathsReshape = [deathsReshape; 
                                        transpose(monthlyTimespan), a.*ones(nTimepoints,1), index.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                        sce.*ones(nTimepoints,1), deaths(:, a, index, param)];
                end 
            end 

            for x = 1 : endpoints 
                if (param == 1 && a == 1 && x == 1)
                    ccHealthStateReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                        sce.*ones(nTimepoints,1), ccHealthState(:, a, x, param)];
                else 
                    ccHealthStateReshape = [ccHealthStateReshape; 
                                        transpose(monthlyTimespan), a.*ones(nTimepoints,1), x.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                        sce.*ones(nTimepoints,1), ccHealthState(:, a, x, param)];
                end 
            end 

            for index = 1 : 7 
                if (param == 1 && a == 1 && index == 1)
                    hpvHealthStateReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), index.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                        sce.*ones(nTimepoints,1), hpvHealthState(:, a, index, param)];
                else 
                    hpvHealthStateReshape = [hpvHealthStateReshape; 
                                        transpose(monthlyTimespan), a.*ones(nTimepoints,1), index.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                        sce.*ones(nTimepoints,1), hpvHealthState(:, a, index, param)];
                end 
            end 

            for index = 1 : 3
                if (param == 1 && a == 1 && index == 1)
                    screenTreatReshape = [transpose(monthlyTimespan), a.*ones(nTimepoints,1), index.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                        sce.*ones(nTimepoints,1), screenTreat(:, a, index, param)];
                else 
                    screenTreatReshape = [screenTreatReshape; 
                                        transpose(monthlyTimespan), a.*ones(nTimepoints,1), index.*ones(nTimepoints,1), param.*ones(nTimepoints,1), ...
                                        sce.*ones(nTimepoints,1), screenTreat(:, a, index, param)];
                end 
            end 

            for x = 1 : 3 
                for index = 1 : 2
                    for treat = 1 : 3
                        if (param==1 && a==1 && x==1)
                            screenSympCCTreatReshape = [transpose(monthlyTimespan), x.*ones(nTimepoints,1), a.*ones(nTimepoints,1), treat.*ones(nTimepoints,1), index.*ones(nTimepoints,1) , param.*ones(nTimepoints,1), ...
                                            sce.*ones(nTimepoints,1), screenSympCCTreat(:, x, a, treat, index, param)];

                        else 
                            screenSympCCTreatReshape = [screenSympCCTreatReshape; 
                                            transpose(monthlyTimespan), x.*ones(nTimepoints,1), a.*ones(nTimepoints,1), treat.*ones(nTimepoints,1), index.*ones(nTimepoints,1) , param.*ones(nTimepoints,1), ...
                                            sce.*ones(nTimepoints,1), screenSympCCTreat(:, x, a, treat, index, param)];
                        end 
                    end 
                end 
            end 
        end
        disp(['Complete Scenario ', num2str(sce), ', Parameter ', num2str(param)])
    end

% turn into arrays
vaxReshape1 = array2table(vaxReshape, 'VariableNames', {'year', 'age', 'vaxType', 'paramNum', ...
        'sceNum', 'count'}); 
deathsReshape1 = array2table(deathsReshape, 'VariableNames', {'year', 'age', 'index', 'paramNum', ...
        'sceNum', 'count'}); 
ccHealthStateReshape1 = array2table(ccHealthStateReshape, 'VariableNames', {'year', 'age', 'endpoint', 'paramNum', ...
        'sceNum', 'count'});
hpvHealthStateReshape1 = array2table(hpvHealthStateReshape, 'VariableNames', {'year', 'age', 'healthState', 'paramNum', ...
        'sceNum', 'count'});
newCCReshape1 = array2table(newCCReshape, 'VariableNames', {'year', 'age', 'paramNum', 'sceNum', 'count'}); 
totalPerAgeReshape1 = array2table(totalPerAgeReshape, 'VariableNames', {'year', 'age', 'paramNum', 'sceNum', 'count'}); 
screenTreatReshape1 = array2table(screenTreatReshape, 'VariableNames', {'year', 'age', 'index', 'paramNum', 'sceNum', 'count'}); 
screenSympCCTreatReshape1 = array2table(screenSympCCTreatReshape, 'VariableNames', {'year', 'endpoint', 'age', 'treat', 'index', 'paramNum', 'sceNum', 'count'}); 

% spit out into CSV 
writetable(vaxReshape1, [pwd '/KECEA/vax_S' sceString '.csv']);
writetable(deathsReshape1, [pwd '/KECEA/deaths_S' sceString '.csv']);
writetable(ccHealthStateReshape1, [pwd '/KECEA/ccHealthState_S' sceString '.csv']);
writetable(hpvHealthStateReshape1, [pwd '/KECEA/hpvHealthState_S' sceString '.csv']);
writetable(newCCReshape1, [pwd '/KECEA/newCC_S' sceString '.csv']); 
writetable(totalPerAgeReshape1, [pwd '/KECEA/totalPerAge_S' sceString '.csv']);
writetable(screenTreatReshape1, [pwd '/KECEA/screenTreat_S' sceString '.csv']);
writetable(screenSympCCTreatReshape1, [pwd '/KECEA/screenSympCCTreat_S' sceString '.csv']);

end % closing the parfor loops 
end