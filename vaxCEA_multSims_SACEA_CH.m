function [] = vaxCEA_multSims_CIs_SACEA(vaxResultInd , sceNum , fileNameNums)
% sceNum is scenario number
% filNameNums is the parameter numbers? 
% example: vaxCEA_multSims_CIs_SACEA(1 , '0' , {'0'})

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
% for k = 1 : loopSegmentsLength-1
%     parfor j = loopSegments{k}+1 : loopSegments{k+1}

k=1;
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
        % ???????? why for vaxResult you start at row 2? Is it because of
        % 2021 being double counted in both? 
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

        % Questions: 
        % - deaths is not stratified by anything. just by year. 
        % - ccDeath and hivDeaths are not stratified by HPV state. Only by
        % vaccine and non vaccine type. Is that okay?

        % Cervical cancer deaths
        nTimepoints = length(monthlyTimespan);

            for d = 1 : disease
                for a = 1 : age
                    for h = 1 : hpvTypeGroups

                        % initialize 2D matrix 
                        ccDeathReshapeTemp = zeros(nTimepoints, ndims(vaxResult{n}.ccDeath)+1); 

                        % set values of the matrix
                        ccDeathReshapeTemp(1 : end, 1 : size(ccDeathReshapeTemp, 2)) = ...
                            [transpose(monthlyTimespan), d.*ones(nTimepoints,1), ...
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


        % HIV deaths
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

        % All cause deaths
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


            % TODO: need to add a section outside the for loop where i
            % append all the reshape tables to append the
            % scenario/parameter data 

        % Write all causes of death info into CSVs for cleaning in R...
        % this should end up outside the for loop
        writetable(ccDeathReshape,[pwd '/SACEA/ccDeaths.csv']);
        writetable(hivDeathReshape, [pwd '/SACEA/hivDeaths.csv']); 
        writetable(allDeathReshape, [pwd '/SACEA/allDeaths.csv']); 

%% TODO: think about how to incorporate scenario #, and parameter # into these matlab tables as well before export. 
% Need to modify R script accordingly 
        
        %% TOTAL NUMBER WHO RECEIVE NONAVALENT HPV VACCINE
        % from vaxCEA_multSims_CIs.m
        % vaxInds = zeros(viral , endpoints , gender , age , risk , disease * 5 * hpvNonVaxStates * intervens); % 5 HPV+ states
        % toInd(allcomb(d,v,h,s,x,p,g,a,r)); note p=vax and screening
        % history
%         vaxTotAge = zeros(nRuns , age , 5 , length(monthlyTimespan));
%         diseaseVec_ccInc = {[1 : disease] , [1 : 2] , [3 : 8] , [3 : 7] , 8}; % creates a cell with these vectors inside 
%         % the above vectors are referring to all the different ways you
%         % could group an HIV disease state (ex: HIV neg, HIV infected no
%         % ART, HIV infected, on ART, etc)
%         diseaseVecLength_ccInc = length(diseaseVec_ccInc); % 5 
% 
%         for dInd = 1 : diseaseVecLength_ccInc
%             d = diseaseVec_ccInc{dInd}; 
%             for a = 1 : age % you are trying to stratify results by d and age -- the two variables you are changing in the for loop
%                 vaxInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : endpoints , [2 , 4] , 2 , a , 1 : risk)); % all the indexes that correspond to vaxInds in popVec
%                 % vaxInds is specifically picking out p=2 and 4, vaccinated
%                 % for HPV
%                 vaxTotAge(j , a , dInd , :) = sum(vaxResult{n}.popVec(: , vaxInds), 2); % popVec vertical dimension is years, horizontal dimension is all the population data that corresponds with vaxInds at those years
%             end
%         end

        %% rewriting the above for this analysis
            % stratify disease in to HIV neg (1,2) and then individual
            % stratifications for 3-8
        diseaseVec_vax = {[1:2], 3, 4, 5, 6, 7, 8}; % 1 and 2 are HIV neg, 3:8 are HIV pos
%         intervensVec_vax = {[3 4], [2 4]}; % 3 is screened, 2 and 4 are vaccinated

        % TODO: i will combine the hpvVax and hpvNonVax states later on

        % initialize
        % time, HIV pos/neg state (2), hpv vax state, hpv non vax state,
        % endpoints (CC/hyst status), age
        % outcome is intervens (screened and non-screened)
        % Will combine HPV vax and non vax state together, and also
        % endpoints if vaxstate
        % hpvVaxState -1 is removing the 6, and then + endpoints is adding
        % CC endpoints 
        vaxTotVaxStates = zeros(length(monthlyTimespan), length(diseaseVec_vax), (hpvVaxStates-1+endpoints), intervens, age); 

        for dInd = 1:length(diseaseVec_vax)
            d = diseaseVec_vax{dInd};

            % no need to stratify by viral load

            % stratify by vaccine type HPV precancer states, but you want
            % to combine the vaccine type and non vaccine type into one
            for h = 1 : hpvVaxStates
            
                for s = 1 : hpvNonVaxStates

                    if h == s & h == 6 % if vax states match up and it's a CC state 

                        for x = 1 : endpoints

                        % no need to stratify by p (intervens: vax and screening
                        % history) because this is the data we are going to pull
        
                        % only want gender = 2 (female) 

                            for pInd = 1:length(intervensVec_vax)
                                p = intervensVec_vax{pInd};
        
                                for a = 1 : age
            
                                % no need to stratify by risk
            
                                    vaxInds = toInd(allcomb(d, 1:viral, h, s, x, ...
                                        p, 2, a, 1 : risk)); %  2 is only for female gender
                
                                    vaxTotVaxStates(:, dInd, h, s, x, pInd, a) = sum(vaxResult{n}.popVec(: , vaxInds), 2);

                            end
                        end
                    end
                end
            end
        end

        % spread vaxTotVaxStates into a 2D matrix
        for d = 1:length(diseaseVec_vax)
            for h = 1:hpvVaxStates
                for s = 1:hpvNonVaxStates
                    for x = 1:endpoints
                        for p = 1:length(intervensVec_vax)
                            for a = 1:age
                                % initalize 2D matrix
                                vaxReshapeTemp = zeros(nTimepoints, ndims(vaxTotVaxStates)); 

                                % set values of the matrix
                                vaxReshapeTemp(:, 1:size(vaxReshapeTemp,2)+1) = ...
                                    [transpose(monthlyTimespan), d.*ones(nTimepoints,1),...
                                    h.*ones(nTimepoints,1), s.*ones(nTimepoints,1), x.*ones(nTimepoints,1), ...
                                    p.*ones(nTimepoints,1), a.*ones(nTimepoints,1), vaxTotVaxStates(:, d, h, s, x, p, a)]; 

                                % append to the end 
                                if exist('vaxReshape') == 0
                                    vaxReshape = vaxReshapeTemp;
                                else 
                                    vaxReshape = [vaxReshape; vaxReshapeTemp];
                                end

                                disp(d)

                            end
                        end
                    end
                end
            end
        end


            % Add scenario and parameter numbers
            vaxReshape = array2table(vaxReshape, ...
            'VariableNames',{'year','disease','vaxTypeHPVState', 'nonvaxTypeHPVState', 'ccState', 'vaxScreenHistory', 'age', 'count'}); 
            vaxReshape.sceNum(:,1) =  sceNum; 
            vaxReshape.paramNum(:,1) = fileInds(j); 
        


        %% TOTAL NUMBER WHO RECEIVE HPV SCREENING
        % newScreen = zeros(length(s) - 1 , disease , hpvVaxStates , hpvNonVaxStates , endpoints , numScreenAge , 2);
        
        
        
%     end
% end

%% Save population size and deaths by gender, age, and HIV/ART status to existing template
ageLabelVec = {2, 7, 12, 17 , 22 , 27 , 32 , 37 , 42 , 47 , 52 , 57 , 62 , 67 , 72 , 77};    % median age of age group
hpvStateVec = {1, 2, 3, 4, 5, 6, 7, 8, 9};
hivStateVec = {1, 3, 4, 5, 6, 2};    % {HIV-negative, CD4>500, CD4350-500, CD4200-350, CD4<200, On ART}
firstYrInd = ((2021 - startYear)*stepsPerYear +1);
t2021on = tVec(firstYrInd:end)';    % row vector, times over timespan
t2021onLen = length(t2021on);    % number of timesteps in timspan
firstYrAnlInd = ((2021 - startYear) +1);
outputVec = [];

% for n = 1 : nRuns
%     for aInd = 1 : age
%         for dInd = 1 : popSizeGAD_disIndsLength
%             for % hpv 9v
%                 for % hpv non-9v
%                     for % cc
%                        outputVec = [outputVec; ...
%                            [ones(t2021onLen,1).*str2num(sceNum) , ...    % column vector, scenario number
%                            (ones(t2021onLen,1).*n) , ...    % column vector, parameter set number
%                            t2021on , ...    % column vector, times
%                            squeeze(popSizeW_multSims((firstYrInd:end) , n , dInd , aInd , __)) , ...    % column vector, number of women of designated category
%                            ones(t2021onLen,1).*ageLabelVec{aInd} , ...    % column vector, age group
%                            ones(t2021onLen,1).*hpvStateVec{dInd} , ...
%                            ones(t2021onLen,1).*hivStateVec{hInd}]];
%                     end
%                 end
%             end
%         end
%     end
% end

fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
'HIV_HPV_CEAtableshell_S' , num2str(fileInd) , '.xlsx'];
writematrix(outputVec , fname , 'Range' , 'A2')

