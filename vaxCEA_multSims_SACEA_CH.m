function [] = vaxCEA_multSims_SACEA_CH(vaxResultInd , sceNum , fileNameNums)
% sceNum is scenario number
% filNameNums is the parameter numbers? 
% vaxResultInd is the number that often comes after the vaxSimResult matlab
% file. Usually just 1. 
% example: vaxCEA_multSims_SACEA_CH(1 , '0' , {'0'})

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

lastYear = 2122;  % ***SET ME***: last year of simulation (use 2122 for SA screening analysis)

%% TODO: Uncomment all the fileInds
% Indices of calib runs to plot
% Temporarily commenting out to only run one scenario first to test out
% code
% fileInds = {'6_1' , '6_2' , '6_3' , '6_6' , '6_8' , '6_9' , '6_11' , ...
%      '6_12' , '6_13' , '6_15' , '6_20' , '6_21' , '6_22' , '6_26' , ...
%     '6_27' , '6_32' , '6_34' , '6_35' , '6_38' , '6_39' , '6_40' , ...
%     '6_41' , '6_42' , '6_45' , '6_47'};    % 22Apr20Ph2V11 
fileInds = {'6_1', '6_6'};
nRuns = length(fileInds);

% Initialize model output plots
% Timespans
monthlyTimespan = [startYear : timeStep : lastYear]; % list all the timespans in a vector
monthlyTimespan = monthlyTimespan(1 : end-1); % remove the very last date
% Population outputs
popSizeW_hivInds = {1:2, 3, 4, 5, 6, 7, 8};    % {HIV-negative, acute & >500, 500-350, 350-200, <=200, ART}, it combines category 1 and 2 for HIV negative
% popSizeW_hivIndsLength = length(popSizeW_disInds); % what is popSizeW_disInds?
% popSizeW_multSims = zeros(length(monthlyTimespan) , nRuns , popSizeW_disIndsLength , age, 1 ); % ch: note I added 1


resultsDir = [pwd , '\HHCoM_Results\'];
fileKey = {'sim1' , 'sim0'};
fileKeyNums = fileNameNums;
n = vaxResultInd;
baseFileName = ['Vaccine22Apr20Ph2V11_2v57BaseVax_spCytoScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_hivInt2017_SA-S' , sceNum , '_']; % ***SET ME***: name for simulation output file
loopSegments = {0 , round(nRuns/2) , nRuns};
loopSegmentsLength = length(loopSegments);
diseaseVec_vax = {[1:2], 3, 4, 5, 6, 7, 8}; % HIV negative grouped together, and then all the HIV positive states 

%% Start of for loop for running each of the parameters 

% for k = 1 : loopSegmentsLength-1
%     parfor j = loopSegments{k}+1 : loopSegments{k+1}

k=1; % temporarily only running a single parameter to test 
j=1; 

% Load results
pathModifier = [baseFileName , fileInds{j}];
nSims = size(dir([pwd , '\HHCoM_Results\' , pathModifier, '\' , '*.mat']) , 1);
curr = load([pwd , '/HHCoM_Results/toNow_22Apr20Ph2V11_2v57BaseVax_spCytoScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_hivInt2017_' , fileInds{j}]); % ***SET ME***: name for historical run output file

vaxResult = cell(nSims , 1);
resultFileName = [pwd , '\HHCoM_Results\' , pathModifier, '\' , 'vaxSimResult'];

% load results from vaccine run into cell array
vaxResult{n} = load([resultFileName , num2str(n), '.mat']);

% concatenate vectors/matrices of population up to current year to population
% matrices for years past current year
% curr is historical model results 
% vaxResult is future model results
% this section of code combines the historical results with futur
% results
% notice for vaxResult you start at row 2. likely because of
% 2021 being double counted in both (?). 
vaxResult{n}.popVec = [curr.popVec(1 : end  , :); vaxResult{n}.popVec(2 : end , :)]; % consolidating historical population numbers with future
vaxResult{n}.ccDeath = [curr.ccDeath(1 : end , : , : , :) ; vaxResult{n}.ccDeath(2 : end , : , : , :)]; % consolidating historical CC death #s with future... etc.
vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :); vaxResult{n}.newCC(2 : end , : , : , :)];
vaxResult{n}.deaths = [curr.deaths(1 : end, 1); vaxResult{n}.deaths(2 : end, 1)];
vaxResult{n}.newHpvVax = [curr.newHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvVax(2 : end , : , : , : , : , :)]; % infected with vaccine type HPV
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
tVec = vaxResult{n}.tVec; % time vector -- length is # of years * 6 time points per year
tVecYr = tVec(1 : stepsPerYear : end); % calculating the number of years. removing the time points in between the whole number years

% Initialize variables
%         hpvYearVec = hpvYearVec_orig;

               
        %% TOTAL POPULATION SIZE BY AGE, HIV STATE, HPV STATE;
%         for a = 1 : age
%             for dInd = 1 : popSizeW_disIndsLength
%                 d = popSizeW_disInds{dInd};
%                 for % hpv 9v
%                     for % hpv non-9v
%                         for % cc
%                             popSizeWinds = toInd(allcomb(d , 1 : viral , 1 , 1 , ... %ch added 1s
%                                 1 , 1 : intervens , 2 , a , 1 : risk)); %ch added 1
%                             popSizeW_multSims(: , j , dInd , a, 1) = sum(vaxResult{n}.popVec(: , popSizeWinds) , 2); %ch added 1
%                         end
%                     end
%                 end
%             end
%         end
    
        
        
        
        
%% TOTAL NUMBER OF DEATHS??? BY AGE, HIV STATE, HPV STATE;
% ccDeath: same format as ccInc = zeros(year, disease , age , hpvTypeGroups);
% hivDeaths: hivDeaths = zeros(length(s) - 1 , disease , gender , age)
% deaths(: not stratified by age?? Just by year x 1
% vector 
% Only want to filter for females
% The first index of all these death matrices is year as a number
% of 1 to 1182. That is the equivalent of monthlyTimeSpan. 

% 1:1182 time spans (monthlyTimeSpan) 
% Disease states is popSizeW_hivInds = {1:2, 3, 4, 5, 6, 7, 8}; 
% Gender is 1 and 2
% 2 hpvTypeGroups: vaccine-type and non-vaccine-type 

%% Cervical cancer deaths
nTimepoints = length(monthlyTimespan);

% note: ccDeath stratifies by whether the CC was due to vaccine type or
% non vaccine type HPV. we don't care to stratify by this in the shell. but
% i will still pull this for now and can remove when cleaning in R. 
for dInd = 1 : length(diseaseVec_vax)
    d = diseaseVec_vax{dInd}; 
    for a = 1 : age
        for h = 1 : hpvTypeGroups
            % initialize 2D matrix 
            ccDeathReshapeTemp = zeros(nTimepoints, ndims(vaxResult{n}.ccDeath)); 

            % set values of the matrix
            ccDeathReshapeTemp(1 : end, 1 : size(ccDeathReshapeTemp, 2)) = ...
                [transpose(monthlyTimespan), dInd.*ones(nTimepoints,1), ...
                a.*ones(nTimepoints,1), h.*ones(nTimepoints,1), vaxResult{n}.ccDeath(1 : end, d, a, h)]; 

            % append to the end 
            if exist('ccDeathReshape') == 0
                ccDeathReshape = ccDeathReshapeTemp;
            else 
                ccDeathReshape = [ccDeathReshape; ccDeathReshapeTemp];
            end
        end 
    end
end

% Turn into table, add scenario and parameter numbers
ccDeathReshape = array2table(ccDeathReshape, ...
    'VariableNames',{'year','disease','age', 'hpvTypeGroups', 'ccDeaths'});
ccDeathReshape.sceNum(:,1) =  sceNum; 
ccDeathReshape.paramNum(:,1) = fileInds(j); 

% Add it to a running array of all the previous parameter runs.
% If else statement if the starting array exists. If yes, add
% it on. If not, create a new array. 
if exist('ccDeathReshapeAllParam') == 0 
    ccDeathReshapeAllParam = ccDeathReshape; 
else 
    ccDeathReshapeAllParam = [ccDeathReshapeAllParam; ccDeathReshape]; 
end 


%% HIV deaths
for d = 1 : disease
    for a = 1 : age
            % initialize 2D matrix 
            hivDeathReshapeTemp = zeros(nTimepoints, ndims(vaxResult{n}.hivDeaths)); % can take out dimension for gender since we're only looking for gender==2 (female) 

            % set values of the matrix
            hivDeathReshapeTemp(1 : end, 1 : size(hivDeathReshapeTemp, 2)) = ...
                [transpose(monthlyTimespan), d.*ones(nTimepoints,1), ...
                a.*ones(nTimepoints,1), vaxResult{n}.hivDeaths(1 : end, d, 2, a)]; % set 2 for gender female 

            % append to the end 
            if exist('hivDeathReshape') == 0
                hivDeathReshape = hivDeathReshapeTemp;
            else 
                hivDeathReshape = [hivDeathReshape; hivDeathReshapeTemp];
            end
    end
end

% Add scenario and parameter numbers
hivDeathReshape = array2table(hivDeathReshape, ...
'VariableNames',{'year','disease','age', 'hivDeaths'}); 
hivDeathReshape.sceNum(:,1) =  sceNum; 
hivDeathReshape.paramNum(:,1) = fileInds(j); 

% Add it to a running array of all the previous parameter runs.
if exist('hivDeathReshapeAllParam') == 0 
    hivDeathReshapeAllParam = hivDeathReshape; 
else 
    hivDeathReshapeAllParam = [hivDeathReshapeAllParam; hivDeathReshape]; 
end 

%% All cause deaths
% initialize 2D matrix 
        allDeathReshapeTemp = zeros(nTimepoints, ndims(vaxResult{n}.deaths)); % can take out dimension for gender since we're only looking for gender==2 (female) 

        % set values of the matrix
        allDeathReshapeTemp(1 : end, 1 : size(allDeathReshapeTemp, 2)) = ...
            [transpose(monthlyTimespan), vaxResult{n}.deaths(1 : end)]; 

        % append to the end 
        if exist('allDeathReshape') == 0
            allDeathReshape = allDeathReshapeTemp;
        else 
            allDeathReshape = [allDeathReshape; allDeathReshapeTemp];
        end

% Add scenario and parameter numbers
allDeathReshape = array2table(allDeathReshape, ...
'VariableNames',{'year','allDeaths'}); 
allDeathReshape.sceNum(:,1) =  sceNum; 
allDeathReshape.paramNum(:,1) = fileInds(j); 

% Add it to a running array of all the previous parameter runs.
if exist('allDeathReshapeAllParam') == 0 
    allDeathReshapeAllParam = allDeathReshape; 
else 
    allDeathReshapeAllParam = [allDeathReshapeAllParam; allDeathReshape]; 
end 
        
%% Stratify by HPV infection and CC state for population counts 

% Description: CEA doesn't care about vaccine / non vaccine type infections. h and s are compared to see what is the more severe HPV state. 
% If HPV state is cancer, then loop through the CC compartments. 
% The combine the h, s, and x compartments into newHpvCcCateg, one column of the CEA data table shell. 
% Output: ceaPopStates contains count info that is stratified by each
% dimension. 

% Initializing ceaPopStates 
    % dim 1 = time 
    % dim 2 = HIV disease groupings
        % 1 = HIV negative
        % 2 = acute
        % 3 = CD4 > 500
        % 4 = CD4 500-350
        % 5 = CD4 350-200
        % 6 = CD4 <= 200
        % 7 = HIV positive, ART 
    % dim 3 = HPV vax groupings + CC endpoint groupings all consolidated 
        % 11 = immune 
        % 1 = susceptible
        % 2 = infected 
        % 3 = CIN1
        % 4 = CIN2 
        % 5 = CIN3 
        % 6 = CC/hysterectomy 
        % 7 = local CC
        % 8 = regional CC 
        % 9 = distant CC 
        % 10 = hysterectomy 
    % dim 4 = intervens (vaccination/screening state)
    % dim 5 = age 
ceaPopStates = zeros(length(monthlyTimespan), length(diseaseVec_vax), (hpvVaxStates+endpoints), intervens, age); 

% Identify indices for popVec
for dInd = 1 : length(diseaseVec_vax)
    d = diseaseVec_vax{dInd}; 

    for h = 0 : hpvVaxStates % note for loop starting at 0. i did this so that the order of vax states is least severe to most severe.

        for s = 0 : hpvNonVaxStates % note for loop starting at 0

            if h <= s % if hpvNonVaxStates is more severe of a state than hpvVaxStates
                newHpvCcCateg = s; 

                if s == 0
                    s1 = 7; % i arbitrarily made loop start with 0 to make comparison of h and s states easier, but it's not actually an index. s1 and h1 is the actual index. 
                    newHpvCcCateg = 11; % update the index for the newHpvCcCateg to 11
                else 
                    s1 = s; 
                end 

                if h == 0
                    h1 = 7; 
                else 
                    h1 = h; 
                end 

                if newHpvCcCateg == 6 % if 6, then loop through CC compartment 

                    for x = 1 : endpoints

                        for p = 1 : intervens % want to stratify by all p for output

                            for a = 1 : age
                                newHpvCcCateg = x + 5; % create a new category that combines HPV and CC states 
                                vaxInds = toInd(allcomb(d, 1:viral, 1:h1, s1, x, p, 2, a, 1:risk)); % 2 is only for female gender; only stratify by s, not h, since s>h
                                ceaPopStates(1:end, dInd, newHpvCcCateg, p, a) = sum(vaxResult{n}.popVec(:, vaxInds), 2);
                            end
                        end
                    end
                else % if you don't need to loop through the CC compartment 
 
                    for p = 1 : intervens % can skip through looping through x 

                        for a = 1 : age
                            vaxInds = toInd(allcomb(d, 1:viral, 1:h1, s1, 1:endpoints, p, 2, a, 1:risk)); % 2 is only for female gender; only stratify by s, not h, since s>h
                            ceaPopStates(1:end, dInd, newHpvCcCateg, p, a) = sum(vaxResult{n}.popVec(:, vaxInds), 2);
                        end 
                    end 
                end 
            else % if s <= h 
                newHpvCcCateg = h; % use h because h is more severe state than s  

                % repeat everything from above, but this time stratifying by h instead of s

                if s == 0
                    s1 = 7; 
                else 
                    s1 = s; 
                end 

                if h == 0
                    h1 = 7; 
                    newHpvCcCateg = 11; 
                else 
                    h1 = h; 
                end

                if newHpvCcCateg == 6 % if 6, then loop through CC compartment 

                    for x = 1 : endpoints

                        for p = 1 : intervens % want to stratify by all p for output

                            for a = 1 : age
                                newHpvCcCateg = x + 5; % create a new category that combines HPV and CC states 
                                vaxInds = toInd(allcomb(d, 1:viral, h1, 1:s1, x, p, 2, a, 1:risk)); % notice stratify by h and not for s
                                ceaPopStates(1:end, dInd, newHpvCcCateg, p, a) = sum(vaxResult{n}.popVec(:, vaxInds), 2);
                            end
                        end
                    end
                else % if you don't need to loop through the CC compartment 

                    for p = 1 : intervens % can skip through looping through x 

                        for a = 1 : age
                            newHpvCcCateg = h; 
                            vaxInds = toInd(allcomb(d, 1:viral, h1, 1:s1, 1:endpoints, p, 2, a, 1:risk)); 
                            ceaPopStates(1:end, dInd, newHpvCcCateg, p, a) = sum(vaxResult{n}.popVec(:, vaxInds), 2);
                        end
                    end 
                end 
            end 
        end 
    end 
end 

% Read through ceaPopStates into a 2D matrix 
% Columns
    % timepoint
    % HIV state 
    % HPV/CC state 
    % vaccination/screening state 
    % age 
    % count

nTimepoints = length(monthlyTimespan);

for d = 1 : length(diseaseVec_vax)
    for g = 1 : size(ceaPopStates, 3) % g refers to newHpvCcCateg
        for p = 1 : intervens
            for a = 1 : age

                % initialize 2D matrix 
                ceaPopStatesReshapeTemp = zeros(nTimepoints, ndims(ceaPopStates)+1); 

                % set values of the matrix 
                ceaPopStatesReshapeTemp(1 : end, 1 : size(ceaPopStatesReshapeTemp,2)) = ...
                    [transpose(monthlyTimespan), d.*ones(nTimepoints,1), ...
                     g.*ones(nTimepoints,1), p.*ones(nTimepoints,1), ...
                     a.*ones(nTimepoints,1), ceaPopStates(1:end, d, g, p, a)]; 

                % append to the end 
                if exist ('ceaPopStatesReshape') == 0 
                    ceaPopStatesReshape = ceaPopStatesReshapeTemp;
                else 
                    ceaPopStatesReshape = [ceaPopStatesReshape; ceaPopStatesReshapeTemp]; 
                end
            end
        end 
    end
end 

% turn into array, add scenario and parameter numbers, save as csv
ceaPopStatesReshape = array2table(ceaPopStatesReshape, ...
    'VariableNames', {'year', 'hivState', 'hpvCcState', 'vaxScreenState', 'age', 'count'}); 
ceaPopStatesReshape.sceNum(:,1) = sceNum; 
ceaPopStatesReshape.paramNum(:,1) = fileInds(j); 

% Add it to a running array of all the previous parameter runs.
if exist('ceaPopStatesReshapeAllParam') == 0 
    ceaPopStatesReshapeAllParam = ceaPopStatesReshape; 
else 
    ceaPopStatesReshapeAllParam = [ceaPopStatesReshapeAllParam; ceaPopStatesReshape]; 
end 

% end -- end for the for loops dictated by k and j 
% end

%% Exporting to CSV 
% Write all causes of death info into CSVs for cleaning in R...
% this should end up outside the for loop for the parameters. 
% TODO. This should end up outsid ethe for loop. And adjust the
% file name to include the scenario number. 
writetable(ccDeathReshapeAllParam,[pwd '/SACEA/ccDeaths_S' sceNum '.csv']);
writetable(hivDeathReshapeAllParam, [pwd '/SACEA/hivDeaths_S' sceNum '.csv']); 
writetable(allDeathReshapeAllParam, [pwd '/SACEA/allDeaths_S' sceNum '.csv']); 
writetable(ceaPopStatesReshapeAllParam, [pwd '/SACEA/ceaPopStates_S' sceNum '.csv']); 

%% TODO: think about how to incorporate scenario #, and parameter # into these matlab tables as well before export. 
% Need to modify R script accordingly 

% %% Save population size and deaths by gender, age, and HIV/ART status to existing template
% ageLabelVec = {2, 7, 12, 17 , 22 , 27 , 32 , 37 , 42 , 47 , 52 , 57 , 62 , 67 , 72 , 77};    % median age of age group
% hpvStateVec = {1, 2, 3, 4, 5, 6, 7, 8, 9};
% hivStateVec = {1, 3, 4, 5, 6, 2};    % {HIV-negative, CD4>500, CD4350-500, CD4200-350, CD4<200, On ART}
% firstYrInd = ((2021 - startYear)*stepsPerYear +1);
% t2021on = tVec(firstYrInd:end)';    % row vector, times over timespan
% t2021onLen = length(t2021on);    % number of timesteps in timspan
% firstYrAnlInd = ((2021 - startYear) +1);
% outputVec = [];
% 
% % for n = 1 : nRuns
% %     for aInd = 1 : age
% %         for dInd = 1 : popSizeGAD_disIndsLength
% %             for % hpv 9v
% %                 for % hpv non-9v
% %                     for % cc
% %                        outputVec = [outputVec; ...
% %                            [ones(t2021onLen,1).*str2num(sceNum) , ...    % column vector, scenario number
% %                            (ones(t2021onLen,1).*n) , ...    % column vector, parameter set number
% %                            t2021on , ...    % column vector, times
% %                            squeeze(popSizeW_multSims((firstYrInd:end) , n , dInd , aInd , __)) , ...    % column vector, number of women of designated category
% %                            ones(t2021onLen,1).*ageLabelVec{aInd} , ...    % column vector, age group
% %                            ones(t2021onLen,1).*hpvStateVec{dInd} , ...
% %                            ones(t2021onLen,1).*hivStateVec{hInd}]];
% %                     end
% %                 end
% %             end
% %         end
% %     end
% % end
% 
% fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
% 'HIV_HPV_CEAtableshell_S' , num2str(fileInd) , '.xlsx'];
% writematrix(outputVec , fname , 'Range' , 'A2')

