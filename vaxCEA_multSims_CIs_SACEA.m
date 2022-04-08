function [] = vaxCEA_multSims_CIs_SACEA(vaxResultInd , sceNum , fileNameNums)
% example: vaxCEA_multSims_CIs_SACEA(1 , '0' , {'0' , '__'})

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
fileInds = {'6_1' , '6_2' , '6_3' , '6_6' , '6_8' , '6_9' , '6_11' , ...
     '6_12' , '6_13' , '6_15' , '6_20' , '6_21' , '6_22' , '6_26' , ...
    '6_27' , '6_32' , '6_34' , '6_35' , '6_38' , '6_39' , '6_40' , ...
    '6_41' , '6_42' , '6_45' , '6_47'};    % 22Apr20Ph2V11
nRuns = length(fileInds);

% Initialize model outputs
% Timespans
monthlyTimespan = [startYear : timeStep : lastYear];
monthlyTimespan = monthlyTimespan(1 : end-1);
% Population outputs
popSizeW_hivInds = {1:2, 3:4, 5, 6, 7, 8};    % {HIV-negative, acute & >500, 500-350, 350-200, <=200, ART}
popSizeW_hivIndsLength = length(popSizeW_disInds);
popSizeW_multSims = zeros(length(monthlyTimespan) , nRuns , popSizeW_disIndsLength , age, ___);


% Select screening algorithm based on scenario number

% Treatment retention (proportion who return and comply with treatment)
%     Note: ideally these parameters would be fed into the script via loadUp2.m
cryoRetain = 0.51; % with three-visit algorithm (cytology + colpo + cryotherapy treatment)
leepRetain = 0.80; % LEEP
thrmlRetain = 0.95; % thermal ablation
ccRetain = 0.40; % cancer treatment
eligLeep = [0.0 , 0.1 , 0.3]; % percent referred to/ eligible for LEEP (CIN1 , CIN2 , CIN3)

if ((str2num(sceNum) == 0) || (str2num(sceNum) == 1)) % Scenarios 0 or 1
    % Screening paper cytology algorithm
    screenAlgs = spCyto;
    screenAlgs.genTypBool = 0;
    % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
    screenAlgs.leepRetain = zeros(1,4);
    screenAlgs.cryoRetain = [0.0 , cryoRetain , cryoRetain , ccRetain];     
    screenAlgs.thrmlRetain = zeros(1,4);    
elseif ((str2num(sceNum) == 2) || (str2num(sceNum) == 3)) % Scenarios 2 or 3
    % Screening paper HPV DNA -and-treat algorithm
    screenAlgs = spHpvDna;
    screenAlgs.genTypBool = 0;
    % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
    screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
    screenAlgs.cryoRetain = zeros(1,4);
    screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain];
elseif ((str2num(sceNum) == 4) || (str2num(sceNum) == 5)) % Scenarios 4 or 5
    % Screening paper HPV DNA+genotyping -and-treat algorithm
    screenAlgs = spGentyp;
    screenAlgs.genTypBool = 1;
    % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
    screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
    screenAlgs.cryoRetain = zeros(1,4);
    screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
elseif ((str2num(sceNum) == 6) || (str2num(sceNum) == 7)) % Scenarios 6 or 7
    % Screening paper AVE -and-treat algorithm
    screenAlgs = spAve;
    screenAlgs.genTypBool = 0;
    % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
    screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
    screenAlgs.cryoRetain = zeros(1,4);
    screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
elseif ((str2num(sceNum) == 8) || (str2num(sceNum) == 9)) % Scenarios 8 or 9
    % Screening paper HPV DNA + AVE triage -and-treat algorithm
    screenAlgs = spHpvAve;
    screenAlgs.genTypBool = 0;
    % proportion who return for treatment (susceptible/immune/infected/CIN1 , CIN2 , CIN3 , CC)
    screenAlgs.leepRetain = [[eligLeep.*leepRetain] , ccRetain]; 
    screenAlgs.cryoRetain = zeros(1,4);
    screenAlgs.thrmlRetain = [[(1-eligLeep).*thrmlRetain] , ccRetain]; 
end


resultsDir = [pwd , '\HHCoM_Results\'];
fileKey = {'sim1' , 'sim0'};
fileKeyNums = fileNameNums;
n = vaxResultInd;
baseFileName = ['Vaccine22Apr20Ph2V11_2v57BaseVax_spCytoScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_hivInt2017_SA-S' , sceNum , '_']; % ***SET ME***: name for simulation output file
loopSegments = {0 , round(nRuns/2) , nRuns};
loopSegmentsLength = length(loopSegments);
for k = 1 : loopSegmentsLength-1
    parfor j = loopSegments{k}+1 : loopSegments{k+1}
        % Load results
        pathModifier = [baseFileName , fileInds{j}];
        nSims = size(dir([pwd , '\HHCoM_Results\' , pathModifier, '\' , '*.mat']) , 1);
        curr = load([pwd , '/HHCoM_Results/toNow_22Apr20Ph2V11_2v57BaseVax_spCytoScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_hivInt2017_' , fileInds{j}]); % ***SET ME***: name for historical run output file 

        vaxResult = cell(nSims , 1);
        resultFileName = [pwd , '\HHCoM_Results\' , pathModifier, '\' , 'vaxSimResult'];
        % load results from vaccine run into cell array
        vaxResult{n} = load([resultFileName , num2str(n), '.mat']);
        % concatenate vectors/matrices of population up to current year (prefix curr) to population
        % matrices for years past current year (prefix vaxResult{n})
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
        
        % Initialize variables
        hpvYearVec = hpvYearVec_orig;

               
        %% TOTAL POPULATION SIZE BY AGE, HIV STATE, HPV STATE;
        for a = 1 : age
            for dInd = 1 : popSizeW_disIndsLength
                d = popSizeW_disInds{dInd};
                for % hpv 9v
                    for % hpv non-9v
                        for % cc
                            popSizeWinds = toInd(allcomb(d , 1 : viral , __ , __ , ...
                                __ , 1 : intervens , 2 , a , 1 : risk));
                            popSizeW_multSims(: , j , dInd , a, __) = sum(vaxResult{n}.popVec(: , popSizeWinds) , 2);
                        end
                    end
                end
            end
        end
    
        
        
        
        
        %% TOTAL NUMBER OF DEATHS??? BY AGE, HIV STATE, HPV STATE;
        
        
        
        %% TOTAL NUMBER WHO RECEIVE NONAVALENT HPV VACCINE
        
        
        
        %% TOTAL NUMBER WHO RECEIVE HPV SCREENING/TREATMENT              
        for a = 1 : age % Sum of newScreen(...) should equal 0 for ages not screened in a given HIV disease state
            for dInd = 1 : popSizeW_disIndsLength % HIV disease states in template
                d = popSizeW_disInds{dInd};
                for h = 1 : hpvVaxStates % Vaccine-type HPV precancer state
                    for s = 1 : hpvNonVaxStates % Non-vaccine-type HPV precancer state
                        for x = 1 : endpoints % Cervical cancer or hysterectomy status
                            % Apply selected screening/treatment algorithm
                                % if you're susceptible/immune to both HPV types or have had a hysterectomy
                                if [( ((h==1) || (h==7)) && ((s==1) || (s==7)) && (x==1) )] || (x==4) || [(screenAlgs.genTypBool && ((h==1) || (h==7)) && (((s>=2) && (s<=5)) || ((s==6) && (x<=3))))]
                                    numScreen(: , j , dInd , hInd, a) = sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , h , s , x , aInd , :),2),3),4),5),6),7);
                                    numLEEP(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                    numCryo(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                    numThrml(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                    numHyst(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                % if you're infected with either HPV type
                                elseif [( ((h==2) && ((s<=2) || (s==7))) || (((h<=2) || (h==7)) && (s==2)) ) && (x==1)] && [(~screenAlgs.genTypBool) || (screenAlgs.genTypBool && (h==2))]
                                    numScreen(: , j , dInd , hInd, a) = sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , h , s , x , aInd , :),2),3),4),5),6),7);
                                    numLEEP(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgs.leepRetain(1);
                                    numCryo(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgscryoRetain(1);
                                    numThrml(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgs.thrmlRetain(1);
                                    numHyst(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                % if you have CIN1 of either HPV type
                                elseif [( ((h==3) && ((s<=3) || (s==7))) || (((h<=3) || (h==7)) && (s==3)) ) && (x==1)] && [(~screenAlgs.genTypBool) || (screenAlgs.genTypBool && (h==3))]
                                    numScreen(: , j , dInd , hInd, a) = sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , h , s , x , aInd , :),2),3),4),5),6),7);
                                    numLEEP(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgs.leepRetain(1);
                                    numCryo(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgscryoRetain(1);
                                    numThrml(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,2) * screenAlgs.colpoRetain * screenAlgs.thrmlRetain(1);
                                    numHyst(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                % if you have CIN2 of either HPV type
                                elseif [( ((h==4) && ((s<=4) || (s==7))) || (((h<=4) || (h==7)) && (s==4)) ) && (x==1)] && [(~screenAlgs.genTypBool) || (screenAlgs.genTypBool && (h==4))]
                                    numScreen(: , j , dInd , hInd, a) = sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , h , s , x , aInd , :),2),3),4),5),6),7);
                                    numLEEP(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,3) * screenAlgs.colpoRetain * screenAlgs.leepRetain(2);
                                    numCryo(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,3) * screenAlgs.colpoRetain * screenAlgscryoRetain(2);
                                    numThrml(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,3) * screenAlgs.colpoRetain * screenAlgs.thrmlRetain(2);
                                    numHyst(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                % if you have CIN3 of either HPV type
                                elseif [( ((h==5) && ((s<=5) || (s==7))) || (((h<=5) || (h==7)) && (s==5)) ) && (x==1)] && [(~screenAlgs.genTypBool) || (screenAlgs.genTypBool && (h==5))]
                                    numScreen(: , j , dInd , hInd, a) = sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , h , s , x , aInd , :),2),3),4),5),6),7);
                                    numLEEP(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,4) * screenAlgs.colpoRetain * screenAlgs.leepRetain(3);
                                    numCryo(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,4) * screenAlgs.colpoRetain * screenAlgscryoRetain(3);
                                    numThrml(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,4) * screenAlgs.colpoRetain * screenAlgs.thrmlRetain(3);
                                    numHyst(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                % if you have cervical cancer
                                elseif [( (x==1) && ((h==6) || (s==6)) ) || (x==2) || (x==3)] && [(~screenAlgs.genTypBool) || (screenAlgs.genTypBool && (h==6) && (x<=3))]
                                    numScreen(: , j , dInd , hInd, a) = sum(sum(sum(sum(sum(sum(vaxResult{n}.newScreen(: , d , h , s , x , aInd , :),2),3),4),5),6),7);
                                    numLEEP(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                    numCryo(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                    numThrml(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * 0.0;
                                    numHyst(: , j , dInd , hInd, a) = numScreen(: , j , dInd , hInd, a) * screenAlgs.testSens(d,4) * screenAlgs.colpoRetain * screenAlgs.treatRetain(4);
                                end

                        end
                    end
                end
            end
        end
        
        
        
        
    end
end

%% Save population size and deaths by gender, age, and HIV/ART status to existing template
ageLabelVec = {2, 7, 12, 17 , 22 , 27 , 32 , 37 , 42 , 47 , 52 , 57 , 62 , 67 , 72 , 77};    % median age of age group
hpvStateVec = {1, 2, 3, 4, 5, 6, 7, 8, 9};
hivStateVec = {1, 3, 4, 5, 6, 2};    % {HIV-negative, CD4>500, CD4350-500, CD4200-350, CD4<200, On ART}
firstYrInd = ((2021 - startYear)*stepsPerYear +1);
t2021on = tVec(firstYrInd:end)';    % row vector, times over timespan
t2021onLen = length(t2021on);    % number of timesteps in timspan
firstYrAnlInd = ((2021 - startYear) +1);
outputVec = [];

for n = 1 : nRuns
    for aInd = 1 : age
        for dInd = 1 : popSizeGAD_disIndsLength
            for % hpv 9v
                for % hpv non-9v
                    for % cc
                       outputVec = [outputVec; ...
                           [ones(t2021onLen,1).*str2num(sceNum) , ...    % column vector, scenario number
                           (ones(t2021onLen,1).*n) , ...    % column vector, parameter set number
                           t2021on , ...    % column vector, times
                           squeeze(popSizeW_multSims((firstYrInd:end) , n , dInd , aInd , __)) , ...    % column vector, number of women of designated category
                           ones(t2021onLen,1).*ageLabelVec{aInd} , ...    % column vector, age group
                           ones(t2021onLen,1).*hpvStateVec{dInd} , ...
                           ones(t2021onLen,1).*hivStateVec{hInd}]];
                    end
                end
            end
        end
    end
end

fname = [pwd , '\HHCoM_Results\' , baseFileName , fileInds{1} , '\' , ...
'HIV_HPV_CEAtableshell_S' , num2str(fileInd) , '.xlsx'];
writematrix(outputVec , fname , 'Range' , 'A2')

