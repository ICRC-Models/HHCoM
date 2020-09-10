%function vaxCEA_plotSaveDoArt_071720()

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
    c3c2Mults , c2c1Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , waning , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear , circNatStartYear , vaxStartYear , ...
    baseline , cisnet , who , whob , circProtect , condProtect , MTCTRate , ...
    hyst , OMEGA , ...
    ccInc2012_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , ...
    hivPrevM_dObs , hivPrevF_dObs , popAgeDist_dObs , totPopSize_dObs , ...
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
    infhpvNonVaxInds , ageInd , riskInd , ...
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3 , dFertMat3 , deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2(1 , 0 , [] , [] , []);

%% LOAD SAVED RESULTS

lastYear = 2061;
% Indices of calib runs to plot
fileInds = {'11_1' , '11_2' , '11_3' , '11_4' , '11_5' , '11_6' , '11_7' , '11_8' , ...
    '11_9' , '11_10' , '11_11' , '11_12' , '11_13' , '11_14' , '11_15' , '11_16' , ...
    '11_17' , '11_18' , '11_19' , '11_20' , '11_21' , '11_22' , '11_23' , '11_24' , '11_25'};  % DO ART, 22Apr20Ph2V2, t=11
% fileInds = {'7_1'}; % , '7_2' , '7_3' , '7_4' , '7_5'}; % , '7_6' , '7_7' , '7_8' , ...
%     '7_9' , '7_10' , '7_11' , '7_12' , '7_13' , '7_14' , '7_15' , '7_16' , ...
%     '7_17' , '7_18' , '7_19' , '7_20' , '7_21' , '7_22' , '7_23' , '7_24' , '7_25'};  % DO ART, 22Apr20Ph2V2, t=11
nRuns = length(fileInds);

hivIncF_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
hivIncM_multSims = hivIncF_multSims;
hivIncC_multSims = hivIncF_multSims;
hivPrevF_multSims = zeros(length([startYear : lastYear-1]) , nRuns , 6);
hivPrevM_multSims = hivPrevF_multSims;
hivPrevC_multSims = hivPrevF_multSims;
artPrevF_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
artPrevM_multSims = artPrevF_multSims;
artPrevC_multSims = artPrevF_multSims;
cd4DistF_multSims = zeros(length([startYear : lastYear-1]) , nRuns , 3);
cd4DistM_multSims = cd4DistF_multSims;
cd4DistC_multSims = cd4DistF_multSims;
hivMortF_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
hivMortM_multSims = hivMortF_multSims;
hivMortC_multSims = hivMortF_multSims;
popSizeF_multSims = zeros(length([startYear : lastYear-1]) , nRuns);
popSizeM_multSims = popSizeF_multSims;
popSizeC_multSims = popSizeF_multSims;

resultsDir = [pwd , '\HHCoM_Results\'];

for j = 1 : nRuns
    % Load results
    baseFileNameShort = '22Apr20Ph2V2_baseVax057_baseScreen_baseVMMC_fertDec042-076-052_2020ARTfxd_DoART';
    baseFileName = [baseFileNameShort , '_S1_'];
    pathModifier = [baseFileName , fileInds{j}]; % ***SET ME***: name for simulation output file
    nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
    curr = load([pwd , '/HHCoM_Results/toNow_22Apr20Ph2V2_baseVax057_baseScreen_baseVMMC_fertDec042-076_2020ARTfxd_DoART_S1_' , fileInds{j}]); % ***SET ME***: name for historical run output file 
    
    vaxResult = cell(nSims , 1);
    resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxSimResult'];
    if waning
        resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxWaneSimResult'];
    end
    for n = nSims
        % load results from vaccine run into cell array
        vaxResult{n} = load([resultFileName , num2str(2), '.mat']);
        % concatenate vectors/matrices of population up to current year to population
        % matrices for years past current year
        vaxResult{n}.popVec = [curr.popVec(1 : end  , :); vaxResult{n}.popVec(2 : end , :)];
        vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , : , : , : , :); vaxResult{n}.newHiv(2 : end , : , : , : , : , : , :)];
        vaxResult{n}.hivDeaths = [curr.hivDeaths(1 : end , : , : , :); vaxResult{n}.hivDeaths(2 : end , : , : , :)];
        vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , : , : , : , : , :); vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , :)];
        vaxResult{n}.tVec = [curr.tVec(1 : end), vaxResult{n}.tVec(2 : end)];
    end
    noVaxInd = nSims;
    noV = vaxResult{noVaxInd};
    tVec = noV.tVec;
    tVecYr = tVec(1 : stepsPerYear : end);
    
    %% HIV INCIDENCE
    % Validation data
    hivInc_obs(: , : , 1) = [2005 2.14 1.57 2.93; % AHRI KZN: (Vandormael, 2019)
                         2006 2.24 1.69 2.96;
                         2007 2.30 1.74 3.05;
                         2008 2.35 1.78 3.09;
                         2009 2.45 1.85 3.24;
                         2010 2.45 1.85 3.25;
                         2011 2.30 1.70 3.11;
                         2012 2.49 1.83 3.37;
                         2013 2.22 1.64 3.01;
                         2014 1.83 1.29 2.59;
                         2015 1.39 0.94 2.07;
                         2016 1.24 0.79 1.95;
                         2017 1.01 0.58 1.76];
    hivInc_obs(: , : , 2) = [2005 4.08 3.40 4.90;
                         2006 4.45 3.77 5.27;
                         2007 4.56 3.86 5.39;
                         2008 4.58 3.89 5.40;
                         2009 4.58 3.85 5.44;
                         2010 4.72 3.98 5.61;
                         2011 4.59 3.85 5.47;
                         2012 4.95 4.14 5.92;
                         2013 4.85 4.05 5.81;
                         2014 4.89 4.09 5.84;
                         2015 4.31 3.58 5.20;
                         2016 3.74 3.04 4.61;
                         2017 3.06 2.38 3.94];

    % Calculate female HIV incidence
    hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 2 , 4 : age , 1 : risk));
    hivSus = annlz(sum(noV.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
    hivIncF = annlz(sum(sum(sum(sum(sum(noV.newHiv(: , : , : , : , 2 , 4 : age , ...
        1 : risk), 2), 3), 4), 6), 7)) ./ hivSus * 100;
    hivIncF_multSims(: , j) = hivIncF(1 : end)';
 
    % Plot female HIV incidence  
    if (j == 1)
        fig1 = figure;
%         errorbar(hivInc_obs(: , 1 , 2) , hivInc_obs(: , 2 , 2) , ...
%             hivInc_obs(: , 2 , 2) - hivInc_obs(: , 3 , 2) , hivInc_obs(: , 4 , 2) - hivInc_obs(: , 2 , 2) , ...
%             'rs' , 'LineWidth' , 1.5);
    else
        figure(fig1);
    end
    hold all;
    plot(tVec(1 : stepsPerYear : end)' , hivIncF(1 : end)' , 'b-')
    xlabel('Year'); ylabel('Incidence per 100'); title('Female HIV incidence, ages 15-79');
    xlim([1980 2060]); ylim([0 10]);
    legend('Model: ages 15-79'); %'(Vandormael, 2019) Observed KZN, ages 15-49: 95% CI' , 
    grid on;

    % Calculate male HIV incidence
    hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 1 , 4 : age , 1 : risk));
    hivSus = annlz(sum(noV.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
    hivIncM = annlz(sum(sum(sum(sum(sum(noV.newHiv(: , : , : , : , 1 , 4 : age , ...
        1 : risk), 2), 3), 4), 6), 7)) ./ hivSus * 100;
    hivIncM_multSims(: , j) = hivIncM(1 : end)';

    % Plot male HIV incidence            
    if (j == 1)
        fig2 = figure;
%         errorbar(hivInc_obs(: , 1 , 1) , hivInc_obs(: , 2 , 1) , ...
%             hivInc_obs(: , 2 , 1) - hivInc_obs(: , 3 , 1) , hivInc_obs(: , 4 , 1) - hivInc_obs(: , 2 , 1) , ...
%             'rs' , 'LineWidth' , 1.5);
    else
        figure(fig2);
    end
    hold all;
    plot(tVec(1 : stepsPerYear : end)' , hivIncM(1 : end)' , 'b-')
    xlabel('Year'); ylabel('Incidence per 100'); title('Male HIV incidence, ages 15-79');
    xlim([1980 2060]); ylim([0 10]);
    legend('Model: ages 15-79'); %'(Vandormael, 2019) Observed KZN, ages 15-49: 95% CI' , 
    grid on;
    
    % Calculate combined HIV incidence
    hivSusInds = toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 1 : 2 , 4 : age , 1 : risk));
    hivSus = annlz(sum(noV.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
    hivIncC = annlz(sum(sum(sum(sum(sum(sum(noV.newHiv(: , : , : , : , 1 : 2 , 4 : age , ...
        1 : risk), 2), 3), 4), 5), 6), 7)) ./ hivSus * 100;
    hivIncC_multSims(: , j) = hivIncC(1 : end)';

    % Plot combined HIV incidence            
    if (j == 1)
        fig13 = figure;
%         errorbar(hivInc_obs(: , 1 , 1) , hivInc_obs(: , 2 , 1) , ...
%             hivInc_obs(: , 2 , 1) - hivInc_obs(: , 3 , 1) , hivInc_obs(: , 4 , 1) - hivInc_obs(: , 2 , 1) , ...
%             'rs' , 'LineWidth' , 1.5);
    else
        figure(fig13);
    end
    hold all;
    plot(tVec(1 : stepsPerYear : end)' , hivIncC(1 : end)' , 'b-')
    xlabel('Year'); ylabel('Incidence per 100'); title('Combined HIV incidence, ages 15-79');
    xlim([1980 2060]); ylim([0 10]);
    legend('Model: ages 15-79'); %'(Vandormael, 2019) Observed KZN, ages 15-49: 95% CI' , 
    grid on;
    
    %% HIV PREVALENCE
    cInds = {3 : 8 , 3 : 5 , 6 , 7 , 8};
    cTits = {'all HIV' , 'CD4 350plus noART' , 'CD4 200-350 noART' , 'CD4 below200 noART' , 'onART'};
    
    % Plot female HIV prevalence      
    if (j == 1)
        fig3 = figure;
        %insert observed data
    else
        figure(fig3);
    end
    for c = 1 : length(cInds)
        % Calculate female HIV prevalence
        hivIndsF = toInd(allcomb(cInds{c} , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 4 : age , 1 : risk));
        totIndsF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 2 , 4 : age , 1 : risk));
        hivPopF = sum(noV.popVec(: , hivIndsF) , 2);
        totPopF = sum(noV.popVec(: , totIndsF) , 2);
        hivPrevF = bsxfun(@rdivide , hivPopF , totPopF);
        hivPrevF_multSims(: , j , c) = hivPrevF(1 : stepsPerYear : end);
        
        subplot(2 , 3 , c);
        hold all;
        plot(tVec(1 : stepsPerYear : end)' , hivPrevF(1 : stepsPerYear : end) , 'b-')
        xlabel('Year'); ylabel('Female HIV prevalence, ages 15-79'); title(cTits{c});
        xlim([1980 2060]); ylim([0.0 0.5]);
        grid on;
    end

    % Plot male HIV prevalence      
    if (j == 1)
        fig4 = figure;
        %insert observed data
    else
        figure(fig4);
    end
    for c = 1 : length(cInds)
        % Calculate male HIV prevalence
        hivIndsM = toInd(allcomb(cInds{c} , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 , 4 : age , 1 : risk));
        totIndsM = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 , 4 : age , 1 : risk));
        hivPopM = sum(noV.popVec(: , hivIndsM) , 2);
        totPopM = sum(noV.popVec(: , totIndsM) , 2);
        hivPrevM = bsxfun(@rdivide , hivPopM , totPopM);
        hivPrevM_multSims(: , j , c) = hivPrevM(1 : stepsPerYear : end);
     
        subplot(2 , 3 , c);
        hold all;
        plot(tVec(1 : stepsPerYear : end)' , hivPrevM(1 : stepsPerYear : end) , 'b-')
        xlabel('Year'); ylabel('Male HIV prevalence, ages 15-79'); title(cTits{c});
        xlim([1980 2060]); ylim([0.0 0.5]);
        grid on;
    end
    
    % Plot combined HIV prevalence      
    if (j == 1)
        fig14 = figure;
        %insert observed data
    else
        figure(fig14);
    end
    for c = 1 : length(cInds)
        % Calculate combined HIV prevalence
        hivIndsC = toInd(allcomb(cInds{c} , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : 2 , 4 : age , 1 : risk));
        totIndsC = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
            1 : intervens , 1 : 2 , 4 : age , 1 : risk));
        hivPopC = sum(noV.popVec(: , hivIndsC) , 2);
        totPopC = sum(noV.popVec(: , totIndsC) , 2);
        hivPrevC = bsxfun(@rdivide , hivPopC , totPopC);
        hivPrevC_multSims(: , j , c) = hivPrevC(1 : stepsPerYear : end);
     
        subplot(2 , 3 , c);
        hold all;
        plot(tVec(1 : stepsPerYear : end)' , hivPrevC(1 : stepsPerYear : end) , 'b-')
        xlabel('Year'); ylabel('Combined HIV prevalence, ages 15-79'); title(cTits{c});
        xlim([1980 2060]); ylim([0.0 0.5]);
        grid on;
    end
    
    %% ART PREVALENCE 
    % Calculate female ART prevalence
    artIndsF = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 2 , 4 : age , 1 : risk));
    hivIndsF = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 2 , 4 : age , 1 : risk));
    artPopF = sum(noV.popVec(: , artIndsF) , 2);
    hivPopF = sum(noV.popVec(: , hivIndsF) , 2);
    artPrevF = bsxfun(@rdivide , artPopF , hivPopF);
    artPrevF_multSims(: , j) = artPrevF(1 : stepsPerYear : end);

    % Plot female ART prevalence      
    if (j == 1)
        fig5 = figure;
        %insert observed data
    else
        figure(fig5);
    end
    hold all;
    plot(tVec(1 : stepsPerYear : end)' , artPrevF(1 : stepsPerYear : end) , 'b-')
    xlabel('Year'); ylabel('Proportion WLWHIV on ART, ages 15-79'); title('Proportion WLWHIV on ART, ages 15-79');
    xlim([1980 2060]); ylim([0.0 1.0]);
    legend('Model: ages 15-79'); % Insert observed data description
    grid on;

    % Calculate male ART prevalence
    artIndsM = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 1 , 4 : age , 1 : risk));
    hivIndsM = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 1 , 4 : age , 1 : risk));
    artPopM = sum(noV.popVec(: , artIndsM) , 2);
    hivPopM = sum(noV.popVec(: , hivIndsM) , 2);
    artPrevM = bsxfun(@rdivide , artPopM , hivPopM);
    artPrevM_multSims(: , j) = artPrevM(1 : stepsPerYear : end);

    % Plot male ART prevalence      
    if (j == 1)
        fig6 = figure;
        %insert observed data
    else
        figure(fig6);
    end
    hold all;
    plot(tVec(1 : stepsPerYear : end)' , artPrevM(1 : stepsPerYear : end) , 'b-')
    xlabel('Year'); ylabel('Proportion MLWHIV on ART, ages 15-79'); title('Proportion MLWHIV on ART, ages 15-79');
    xlim([1980 2060]); ylim([0.0 1.0]);
    legend('Model: ages 15-79'); % Insert observed data description
    grid on;
    
    % Calculate combined ART prevalence
    artIndsC = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 1 : 2 , 4 : age , 1 : risk));
    hivIndsC = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 1 : 2 , 4 : age , 1 : risk));
    artPopC = sum(noV.popVec(: , artIndsC) , 2);
    hivPopC = sum(noV.popVec(: , hivIndsC) , 2);
    artPrevC = bsxfun(@rdivide , artPopC , hivPopC);
    artPrevC_multSims(: , j) = artPrevC(1 : stepsPerYear : end);

    % Plot combined ART prevalence      
    if (j == 1)
        fig15 = figure;
        %insert observed data
    else
        figure(fig15);
    end
    hold all;
    plot(tVec(1 : stepsPerYear : end)' , artPrevC(1 : stepsPerYear : end) , 'b-')
    xlabel('Year'); ylabel('Proportion PLWHIV on ART, ages 15-79'); title('Proportion PLWHIV on ART, ages 15-79');
    xlim([1980 2060]); ylim([0.0 1.0]);
    legend('Model: ages 15-79'); % Insert observed data description
    grid on;
    
    %% CD4 DISTRIBUTION AT ART INITIATION
    cInds = {3 : 5 , 6 , 7};
    cTits = {'CD4 350plus' , 'CD4 200-350' , 'CD4 below200'};
    
    % Plot female distribution
    if (j == 1)
        fig7 = figure;
        %insert observed data
    else
        figure(fig7);
    end
    for c = 1 : length(cInds)
        % Calculate female distribution
        init_subF = sum(sum(sum(sum(noV.artTreatTracker(: , cInds{c} , : , 2 , 4 : age , 1 : risk), 2), 3), 5), 6);
        init_allCD4F = sum(sum(sum(sum(noV.artTreatTracker(: , 3 : 7 , : , 2 , 4 : age , 1 : risk), 2), 3), 5), 6);
        cd4DistF = init_subF ./ init_allCD4F;
        cd4DistF_multSims(: , j , c) = cd4DistF(1 : stepsPerYear : end)';
        
        subplot(1 ,3 , c);
        hold all;
        plot(tVec(1 : stepsPerYear : end)' , cd4DistF(1 : stepsPerYear : end)' , 'b-')
        xlabel('Year'); ylabel('Proportion of females initiating ART, ages 15-79'); title(['Proportion CD4: ' , cTits{c}]);
        xlim([1980 2060]); ylim([0.0 1.0]);
        grid on;
    end

    % Plot male distribution
    if (j == 1)
        fig8 = figure;
        %insert observed data
    else
        figure(fig8);
    end
    for c = 1 : length(cInds)
        % Calculate male distribution
        init_subM = sum(sum(sum(sum(noV.artTreatTracker(: , cInds{c} , : , 1 , 4 : age , 1 : risk), 2), 3), 5), 6);
        init_allCD4M = sum(sum(sum(sum(noV.artTreatTracker(: , 3 : 7 , : , 1 , 4 : age , 1 : risk), 2), 3), 5), 6);
        cd4DistM = init_subM ./ init_allCD4M;
        cd4DistM_multSims(: , j , c) = cd4DistM(1 : stepsPerYear : end)';

        subplot(1 , 3 , c);
        hold all;
        plot(tVec(1 : stepsPerYear : end)' , cd4DistM(1 : stepsPerYear : end)' , 'b-')
        xlabel('Year'); ylabel('Proportion of males initiating ART, ages 15-79'); title(['Proportion CD4: ' , cTits{c}]);
        xlim([1980 2060]); ylim([0.0 1.0]);
        grid on;
    end
    
    % Plot combined distribution
    if (j == 1)
        fig16 = figure;
        %insert observed data
    else
        figure(fig16);
    end
    for c = 1 : length(cInds)
        % Calculate combined distribution
        init_subC = sum(sum(sum(sum(sum(noV.artTreatTracker(: , cInds{c} , : , 1 : 2 , 4 : age , 1 : risk), 2), 3), 4), 5), 6);
        init_allCD4C = sum(sum(sum(sum(sum(noV.artTreatTracker(: , 3 : 7 , : , 1 : 2 , 4 : age , 1 : risk), 2), 3), 4), 5), 6);
        cd4DistC = init_subC ./ init_allCD4C;
        cd4DistC_multSims(: , j , c) = cd4DistC(1 : stepsPerYear : end)';

        subplot(1 , 3 , c);
        hold all;
        plot(tVec(1 : stepsPerYear : end)' , cd4DistC(1 : stepsPerYear : end)' , 'b-')
        xlabel('Year'); ylabel('Proportion of persons initiating ART, ages 15-79'); title(['Proportion CD4: ' , cTits{c}]);
        xlim([1980 2060]); ylim([0.0 1.0]);
        grid on;
    end
    
    %% HIV-ASSOCIATED MORTALITY
    
    % Calculate female HIV-associated mortality
    popIndsF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 2 , 4 : age , 1 : risk));
    popTotF = annlz(sum(noV.popVec(: , popIndsF) , 2)) ./ stepsPerYear;
    hivMortF = annlz(sum(sum(noV.hivDeaths(: , : , 2 , 4 : age), 2), 4)) ./ popTotF * 100000;
    hivMortF_multSims(: , j) = hivMortF(1 : end)';

    % Plot female HIV-associated mortality
    if (j == 1)
        fig9 = figure;
        %insert observed data
    else
        figure(fig9);
    end
    hold all;
    plot(tVec(1 : stepsPerYear : end)' , hivMortF(1 : end)' , 'b-')
    xlabel('Year'); ylabel('Mortality per 100K'); title('Female HIV-associated mortality, ages 15-79');
    xlim([1980 2060]); ylim([0 5000]);
    legend('Model: ages 15-79');
    grid on;

    % Calculate male HIV-associated mortality
    popIndsM = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 1 , 4 : age , 1 : risk));
    popTotM = annlz(sum(noV.popVec(: , popIndsM) , 2)) ./ stepsPerYear;
    hivMortM = annlz(sum(sum(noV.hivDeaths(: , : , 1 , 4 : age), 2), 4)) ./ popTotM * 100000;
    hivMortM_multSims(: , j) = hivMortM(1 : end)';

    % Plot male HIV-associated mortality
    if (j == 1)
        fig10 = figure;
        %insert observed data
    else
        figure(fig10);
    end
    hold all;
    plot(tVec(1 : stepsPerYear : end)' , hivMortM(1 : end)' , 'b-')
    xlabel('Year'); ylabel('Mortality per 100K'); title('Male HIV-associated mortality, ages 15-79');
    xlim([1980 2060]); ylim([0 5000]);
    legend('Model: ages 15-79');
    grid on;
    
    % Calculate combined HIV-associated mortality
    popIndsC = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , ...
        1 : intervens , 1 : 2 , 4 : age , 1 : risk));
    popTotC = annlz(sum(noV.popVec(: , popIndsC) , 2)) ./ stepsPerYear;
    hivMortC = annlz(sum(sum(sum(noV.hivDeaths(: , : , 1 : 2 , 4 : age), 2), 3), 4)) ./ popTotC * 100000;
    hivMortC_multSims(: , j) = hivMortC(1 : end)';

    % Plot combined HIV-associated mortality
    if (j == 1)
        fig17 = figure;
        %insert observed data
    else
        figure(fig17);
    end
    hold all;
    plot(tVec(1 : stepsPerYear : end)' , hivMortC(1 : end)' , 'b-')
    xlabel('Year'); ylabel('Mortality per 100K'); title('Combined HIV-associated mortality, ages 15-79');
    xlim([1980 2060]); ylim([0 5000]);
    legend('Model: ages 15-79');
    grid on;
    
    %% TOTAL POPULATION SIZE
    
    % Calculate female population size
    popTotF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 4 : age , 1 : risk));
    popSizeF = sum(noV.popVec(: , popTotF) , 2);
    popSizeF_multSims(: , j) = popSizeF(1 : stepsPerYear : end);

    % Plot female population size
%     if (j == 1)
%         fig11 = figure;
%         %insert observed data
%     else
%         figure(fig11);
%     end
%     hold all;
%     plot(tVec(1 : stepsPerYear : end)' , popSizeF(1 : stepsPerYear : end) , 'b-')
%     xlabel('Year'); ylabel('Individuals'); title('Female population size, ages 15-79');
%     xlim([1980 2060])
%     legend('Model: ages 15-79');
%     grid on;

    % Calculate male population size
    popTotM = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 , 4 : age , 1 : risk));
    popSizeM = sum(noV.popVec(: , popTotM) , 2);
    popSizeM_multSims(: , j) = popSizeM(1 : stepsPerYear : end);

    % Plot female population size
%     if (j == 1)
%         fig12 = figure;
%         %insert observed data
%     else
%         figure(fig12);
%     end
%     hold all;
%     plot(tVec(1 : stepsPerYear : end)' , popSizeM(1 : stepsPerYear : end) , 'b-')
%     xlabel('Year'); ylabel('Individuals'); title('Male population size, ages 15-79');
%     xlim([1980 2060])
%     legend('Model: ages 15-79');
%     grid on;
    
    % Calculate combined population size
    popTotC = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : 2 , 4 : age , 1 : risk));
    popSizeC = sum(noV.popVec(: , popTotC) , 2);
    popSizeC_multSims(: , j) = popSizeC(1 : stepsPerYear : end);

    % Plot combined population size
%     if (j == 1)
%         fig18 = figure;
%         %insert observed data
%     else
%         figure(fig18);
%     end
%     hold all;
%     plot(tVec(1 : stepsPerYear : end)' , popSizeC(1 : stepsPerYear : end) , 'b-')
%     xlabel('Year'); ylabel('Individuals'); title('Combined population size, ages 15-79');
%     xlim([1980 2060])
%     legend('Model: ages 15-79');
%     grid on;

end
    
%% Save HIV incidence
% female
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'HIV_incidence_females_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivIncF_multSims , 2) , ...
    min(hivIncF_multSims , [] , 2) , max(hivIncF_multSims , [] , 2) , ...
    hivIncF_multSims] , fname)
% male
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'HIV_incidence_males_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivIncM_multSims , 2) , ...
    min(hivIncM_multSims , [] , 2) , max(hivIncM_multSims , [] , 2) , ...
    hivIncM_multSims] , fname)
% combined
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'HIV_incidence_combined_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivIncC_multSims , 2) , ...
    min(hivIncC_multSims , [] , 2) , max(hivIncC_multSims , [] , 2) , ...
    hivIncC_multSims] , fname)

%% Save HIV prevalence
cInds = {3 : 8 , 3 : 5 , 6 , 7 , 8};
cTits = {'all_HIV' , 'CD4_350plus_noART' , 'CD4_200-350_noART' , 'CD4_below200_noART' , 'onART'};
for c = 1 : length(cInds)
    % female
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
        'HIV_prevalence_females_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(hivPrevF_multSims(: , : , c) , 2) , ...
        min(hivPrevF_multSims(: , : , c) , [] , 2) , max(hivPrevF_multSims(: , : , c) , [] , 2) , ...
        hivPrevF_multSims(: , : , c)] , fname)
    % male
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
        'HIV_prevalence_males_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(hivPrevM_multSims(: , : , c) , 2) , ...
        min(hivPrevM_multSims(: , : , c) , [] , 2) , max(hivPrevM_multSims(: , : , c) , [] , 2) , ...
        hivPrevM_multSims(: , : , c)] , fname)
    % combined
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
        'HIV_prevalence_combined_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(hivPrevC_multSims(: , : , c) , 2) , ...
        min(hivPrevC_multSims(: , : , c) , [] , 2) , max(hivPrevC_multSims(: , : , c) , [] , 2) , ...
        hivPrevC_multSims(: , : , c)] , fname)
end

%% Save ART prevalence
% female
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'Proportion_WLWHIV_onART_aged15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(artPrevF_multSims , 2) , ...
    min(artPrevF_multSims , [] , 2) , max(artPrevF_multSims , [] , 2) , ...
    artPrevF_multSims] , fname)
% male
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'Proportion_MLWHIV_onART_aged15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(artPrevM_multSims , 2) , ...
    min(artPrevM_multSims , [] , 2) , max(artPrevM_multSims , [] , 2) , ...
    artPrevM_multSims] , fname)
% combined
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'Proportion_PLWHIV_onART_aged15-79.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(artPrevC_multSims , 2) , ...
    min(artPrevC_multSims , [] , 2) , max(artPrevC_multSims , [] , 2) , ...
    artPrevC_multSims] , fname)

%% Save CD4 distribution at ART initiation
cInds = {3 : 5 , 6 , 7};
cTits = {'CD4_350plus' , 'CD4_200-350' , 'CD4_below200'};
for c = 1 : length(cInds)
    % female
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
        'CD4_dist_initART_females_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(cd4DistF_multSims(: , : , c) , 2) , ...
        min(cd4DistF_multSims(: , : , c) , [] , 2) , max(cd4DistF_multSims(: , : , c) , [] , 2) , ...
        cd4DistF_multSims(: , : , c)] , fname)
    % male 
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
        'CD4_dist_initART_males_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(cd4DistM_multSims(: , : , c) , 2) , ...
        min(cd4DistM_multSims(: , : , c) , [] , 2) , max(cd4DistM_multSims(: , : , c) , [] , 2) , ...
        cd4DistM_multSims(: , : , c)] , fname)
    % combined 
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
        'CD4_dist_initART_combined_aged15-79_' , cTits{c} , '.csv'];
    writematrix([tVec(1 : stepsPerYear : end)' , mean(cd4DistC_multSims(: , : , c) , 2) , ...
        min(cd4DistC_multSims(: , : , c) , [] , 2) , max(cd4DistC_multSims(: , : , c) , [] , 2) , ...
        cd4DistC_multSims(: , : , c)] , fname)
end

%% Save HIV mortality
% female
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'HIV_mortality_females_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivMortF_multSims , 2) , ...
    min(hivMortF_multSims , [] , 2) , max(hivMortF_multSims , [] , 2) , ...
    hivMortF_multSims] , fname)
% male
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'HIV_mortality_males_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivMortM_multSims , 2) , ...
    min(hivMortM_multSims , [] , 2) , max(hivMortM_multSims , [] , 2) , ...
    hivMortM_multSims] , fname)
% combined
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'HIV_mortality_combined_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(hivMortC_multSims , 2) , ...
    min(hivMortC_multSims , [] , 2) , max(hivMortC_multSims , [] , 2) , ...
    hivMortC_multSims] , fname)

%% Save population size
% female
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'PopulationSize_females_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(popSizeF_multSims , 2) , ...
    min(popSizeF_multSims , [] , 2) , max(popSizeF_multSims , [] , 2) , ...
    popSizeF_multSims] , fname)
% male 
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'PopulationSize_males_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(popSizeM_multSims , 2) , ...
    min(popSizeM_multSims , [] , 2) , max(popSizeM_multSims , [] , 2) , ...
    popSizeM_multSims] , fname)
% combined
fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileName , '11_1' , '\' , ...
    'PopulationSize_combined_aged15-79' , '.csv'];
writematrix([tVec(1 : stepsPerYear : end)' , mean(popSizeC_multSims , 2) , ...
    min(popSizeC_multSims , [] , 2) , max(popSizeC_multSims , [] , 2) , ...
    popSizeC_multSims] , fname)


%% PLOT SCENARIOS COMPARED %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorList = {'k' , 'b' , 'r'};

%% Plot population size
figure;

% combined
subplot(1 ,3 , 3);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'PopulationSize_combined_aged15-79' , '.csv'];
    popSizeC = xlsread(fname);
    if i == 1
        hold all;
        plot(popSizeC(: , 1) , popSizeC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(popSizeC(: , 1) , popSizeC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(popSizeC(: , 1) , popSizeC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(popSizeC((2019-startYear)+1 : end , 1) , popSizeC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(popSizeC((2019-startYear)+1 : end , 1) , popSizeC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(popSizeC((2019-startYear)+1 : end , 1) , popSizeC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('Individuals'); title('Females + Males, aged 15-79');
xlim([1980 2060]); ylim([0 14*(10^6)]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max');
grid on;
sgtitle('Population size');

%% Plot HIV incidence
figure;

% female
subplot(1 ,3 , 1);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'HIV_incidence_females_aged15-79' , '.csv'];
    hivIncF = xlsread(fname);
    if i == 1
        hold all;
        plot(hivIncF(: , 1) , hivIncF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivIncF(: , 1) , hivIncF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivIncF(: , 1) , hivIncF(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(hivIncF((2019-startYear)+1 : end , 1) , hivIncF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivIncF((2019-startYear)+1 : end , 1) , hivIncF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivIncF((2019-startYear)+1 : end , 1) , hivIncF((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('HIV incidence per 100'); title('Females, aged 15-79');
xlim([1980 2060]); ylim([0 9]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max');
grid on;

% male
subplot(1 ,3 , 2);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'HIV_incidence_males_aged15-79' , '.csv'];
    hivIncM = xlsread(fname);
    if i == 1
        hold all;
        plot(hivIncM(: , 1) , hivIncM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivIncM(: , 1) , hivIncM(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivIncM(: , 1) , hivIncM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(hivIncM((2019-startYear)+1 : end , 1) , hivIncM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivIncM((2019-startYear)+1 : end , 1) , hivIncM((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivIncM((2019-startYear)+1 : end , 1) , hivIncM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('HIV incidence per 100'); title('Males, aged 15-79');
xlim([1980 2060]); ylim([0 9]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max');
grid on;

% combined
subplot(1 ,3 , 3);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'HIV_incidence_combined_aged15-79' , '.csv'];
    hivIncC = xlsread(fname);
    if i == 1
        hold all;
        plot(hivIncC(: , 1) , hivIncC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivIncC(: , 1) , hivIncC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivIncC(: , 1) , hivIncC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(hivIncC((2019-startYear)+1 : end , 1) , hivIncC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivIncC((2019-startYear)+1 : end , 1) , hivIncC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivIncC((2019-startYear)+1 : end , 1) , hivIncC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('HIV incidence per 100'); title('Females + Males, aged 15-79');
xlim([1980 2060]); ylim([0 9]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max');
grid on;
sgtitle('HIV Incidence');

%% Plot HIV prevalence by CD4
cInds = {3 : 8 , 3 : 5 , 6 , 7 , 8};
cTits = {'all_HIV' , 'CD4_350plus_noART' , 'CD4_200-350_noART' , 'CD4_below200_noART' , 'onART'};
cTitsPlot = {'all HIV' , 'CD4 350plus noART' , 'CD4 200-350 noART' , 'CD4 below200 noART' , 'onART'};

% female
figure;
for c = 1 : length(cInds)
    subplot(2 , 3 , c);
    for i = 1 : 3
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
            'HIV_prevalence_females_aged15-79_' , cTits{c} , '.csv'];
        hivPrevF = xlsread(fname);
        if i == 1
            hold all;
            plot(hivPrevF(: , 1) , hivPrevF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(hivPrevF(: , 1) , hivPrevF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(hivPrevF(: , 1) , hivPrevF(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        else
            hold all;
            plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        end
    end
    xlabel('Year'); ylabel('HIV prevalence'); title(cTitsPlot{c}); 
    xlim([1980 2060]); ylim([0.0 0.5]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
        'Scenario 3, mean' , 'min' , 'max');
    grid on;
end
sgtitle('HIV Prevalence: Females, aged 15-79');

% male
figure;
for c = 1 : length(cInds)
    subplot(2 , 3 , c);
    for i = 1 : 3
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
            'HIV_prevalence_males_aged15-79_' , cTits{c} , '.csv'];
        hivPrevM = xlsread(fname);
        if i == 1
            hold all;
            plot(hivPrevM(: , 1) , hivPrevM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(hivPrevM(: , 1) , hivPrevM(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(hivPrevM(: , 1) , hivPrevM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        else
            hold all;
            plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        end
    end
    xlabel('Year'); ylabel('HIV prevalence'); title(cTitsPlot{c}); 
    xlim([1980 2060]); ylim([0.0 0.3]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
        'Scenario 3, mean' , 'min' , 'max');
    grid on;
end   
sgtitle('HIV Prevalence: Males, aged 15-79');
    
% combined
figure;
for c = 1 : length(cInds)
    subplot(2 , 3 , c);
    for i = 1 : 3
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
            'HIV_prevalence_combined_aged15-79_' , cTits{c} , '.csv'];
        hivPrevC = xlsread(fname);
        if i == 1
            hold all;
            plot(hivPrevC(: , 1) , hivPrevC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(hivPrevC(: , 1) , hivPrevC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(hivPrevC(: , 1) , hivPrevC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        else
            hold all;
            plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        end
    end
    xlabel('Year'); ylabel('HIV prevalence'); title(cTitsPlot{c}); 
    xlim([1980 2060]); ylim([0.0 0.4]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
        'Scenario 3, mean' , 'min' , 'max');
    grid on;
end   
sgtitle('HIV Prevalence: Females + Males, aged 15-79');

%% Plot HIV prevalence overall
figure;

% female
subplot(1 , 3 , 1);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'HIV_prevalence_females_aged15-79_all_HIV' , '.csv'];
    hivPrevF = xlsread(fname);
    if i == 1
        hold all;
        plot(hivPrevF(: , 1) , hivPrevF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivPrevF(: , 1) , hivPrevF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivPrevF(: , 1) , hivPrevF(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivPrevF((2019-startYear)+1 : end , 1) , hivPrevF((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('HIV prevalence'); title('Females, aged 15-79'); 
xlim([1980 2060]); ylim([0.0 0.5]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max');
grid on;

% male
subplot(1 , 3 , 2);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'HIV_prevalence_males_aged15-79_all_HIV' , '.csv'];
    hivPrevM = xlsread(fname);
    if i == 1
        hold all;
        plot(hivPrevM(: , 1) , hivPrevM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivPrevM(: , 1) , hivPrevM(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivPrevM(: , 1) , hivPrevM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivPrevM((2019-startYear)+1 : end , 1) , hivPrevM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('HIV prevalence'); title('Males, aged 15-79'); 
xlim([1980 2060]); ylim([0.0 0.5]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max');
grid on;

% combined
subplot(1 , 3 , 3);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'HIV_prevalence_combined_aged15-79_all_HIV' , '.csv'];
    hivPrevC = xlsread(fname);
    if i == 1
        hold all;
        plot(hivPrevC(: , 1) , hivPrevC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivPrevC(: , 1) , hivPrevC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivPrevC(: , 1) , hivPrevC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivPrevC((2019-startYear)+1 : end , 1) , hivPrevC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('HIV prevalence'); title('Females + Males, aged 15-79'); 
xlim([1980 2060]); ylim([0.0 0.5]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max');
grid on;
sgtitle('HIV Prevalence');

%% Plot ART prevalence
figure;

% female
subplot(1 , 3 , 1);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'Proportion_WLWHIV_onART_aged15-79.csv'];
    artPrevF = xlsread(fname);
    if i == 1
        hold all;
        plot(artPrevF(: , 1) , artPrevF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(artPrevF(: , 1) , artPrevF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(artPrevF(: , 1) , artPrevF(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(artPrevF((2019-startYear)+1 : end , 1) , artPrevF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(artPrevF((2019-startYear)+1 : end , 1) , artPrevF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(artPrevF((2019-startYear)+1 : end , 1) , artPrevF((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('Proportion WLWHIV VS on ART'); title('WLWHIV on ART, ages 15-79');
xlim([2000 2060]); ylim([0.0 1.0]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max' , 'Location' , 'SouthEast');
grid on;

% male
subplot(1 , 3 , 2);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'Proportion_MLWHIV_onART_aged15-79.csv'];
    artPrevM = xlsread(fname);
    if i == 1
        hold all;
        plot(artPrevM(: , 1) , artPrevM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(artPrevM(: , 1) , artPrevM(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(artPrevM(: , 1) , artPrevM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(artPrevM((2019-startYear)+1 : end , 1) , artPrevM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(artPrevM((2019-startYear)+1 : end , 1) , artPrevM((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(artPrevM((2019-startYear)+1 : end , 1) , artPrevM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('Proportion MLWHIV VS on ART'); title('MLWHIV on ART, ages 15-79');
xlim([2000 2060]); ylim([0.0 1.0]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max' , 'Location' , 'SouthEast');
grid on;

% combined
subplot(1 , 3 , 3);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'Proportion_PLWHIV_onART_aged15-79.csv'];
    artPrevC = xlsread(fname);
    if i == 1
        hold all;
        plot(artPrevC(: , 1) , artPrevC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(artPrevC(: , 1) , artPrevC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(artPrevC(: , 1) , artPrevC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(artPrevC((2019-startYear)+1 : end , 1) , artPrevC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(artPrevC((2019-startYear)+1 : end , 1) , artPrevC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(artPrevC((2019-startYear)+1 : end , 1) , artPrevC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('Proportion PLWHIV VS on ART'); title('PLWHIV on ART, ages 15-79');
xlim([2000 2060]); ylim([0.0 1.0]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max' , 'Location' , 'SouthEast');
grid on;
sgtitle('Coverage of PLWHIV on ART and virally suppressed');

%% Plot CD4 distribution at ART initiation
cInds = {3 : 5 , 6 , 7};
cTits = {'CD4_350plus' , 'CD4_200-350' , 'CD4_below200'};
cTitsPlot = {'CD4 350plus' , 'CD4 200-350' , 'CD4 below200'};

% female
figure;
for c = 1 : length(cInds)
    subplot(1 , 3 , c);
    for i = 1 : 3
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
            'CD4_dist_initART_females_aged15-79_' , cTits{c} , '.csv'];
        cd4DistF = xlsread(fname);
        if i == 1
            hold all;
            plot(cd4DistF(: , 1) , cd4DistF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(cd4DistF(: , 1) , cd4DistF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i})'
            plot(cd4DistF(: , 1) , cd4DistF(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        else
            hold all;
            plot(cd4DistF((2019-startYear)+1 : end , 1) , cd4DistF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(cd4DistF((2019-startYear)+1 : end , 1) , cd4DistF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(cd4DistF((2019-startYear)+1 : end , 1) , cd4DistF((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        end
    end
    xlabel('Year'); ylabel('Proportion of females initiating ART'); title(['Proportion CD4: ' , cTitsPlot{c}]);
    xlim([1980 2060]); ylim([0.0 1.0]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
        'Scenario 3, mean' , 'min' , 'max');
    grid on;
end
sgtitle('Proportion of persons initiating ART: Females, aged 15-79');
   
% male 
figure;
for c = 1 : length(cInds)
    subplot(1 , 3 , c);
    for i = 1 : 3
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
            'CD4_dist_initART_males_aged15-79_' , cTits{c} , '.csv'];
        cd4DistM = xlsread(fname);
        if i == 1 
            hold all;
            plot(cd4DistM(: , 1) , cd4DistM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(cd4DistM(: , 1) , cd4DistM(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(cd4DistM(: , 1) , cd4DistM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        else
            hold all;
            plot(cd4DistM((2019-startYear)+1 : end , 1) , cd4DistM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(cd4DistM((2019-startYear)+1 : end , 1) , cd4DistM((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(cd4DistM((2019-startYear)+1 : end , 1) , cd4DistM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        end
    end
    xlabel('Year'); ylabel('Proportion of males initiating ART'); title(['Proportion CD4: ' , cTitsPlot{c}]);
    xlim([1980 2060]); ylim([0.0 1.0]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
        'Scenario 3, mean' , 'min' , 'max');
    grid on;
end
sgtitle('Proportion of persons initiating ART: Males, aged 15-79');
   
% combined 
figure;
for c = 1 : length(cInds)
    subplot(1 , 3 , c);
    for i = 1 : 3
        fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
            'CD4_dist_initART_combined_aged15-79_' , cTits{c} , '.csv'];
        cd4DistC = xlsread(fname);
        if i == 1
            hold all;
            plot(cd4DistC(: , 1) , cd4DistC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(cd4DistC(: , 1) , cd4DistC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(cd4DistC(: , 1) , cd4DistC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        else
            hold all;
            plot(cd4DistC((2019-startYear)+1 : end , 1) , cd4DistC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
            plot(cd4DistC((2019-startYear)+1 : end , 1) , cd4DistC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
            plot(cd4DistC((2019-startYear)+1 : end , 1) , cd4DistC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
        end
    end
    xlabel('Year'); ylabel('Proportion of persons initiating ART'); title(['Proportion CD4: ' , cTitsPlot{c}]);
    xlim([1980 2060]); ylim([0.0 1.0]);
    legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
        'Scenario 3, mean' , 'min' , 'max');
    grid on;
end
sgtitle('Proportion of persons initiating ART: Females + Males, aged 15-79');

%% Plot HIV mortality
figure;

% female
subplot(1 , 3 , 1);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'HIV_mortality_females_aged15-79' , '.csv'];
    hivMortF = xlsread(fname);
    if i == 1
        hold all;
        plot(hivMortF(: , 1) , hivMortF(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivMortF(: , 1) , hivMortF(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivMortF(: , 1) , hivMortF(: , 4) , 'LineStyle' ,  '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(hivMortF((2019-startYear)+1 : end , 1) , hivMortF((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivMortF((2019-startYear)+1 : end , 1) , hivMortF((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivMortF((2019-startYear)+1 : end , 1) , hivMortF((2019-startYear)+1 : end , 4) , 'LineStyle' ,  '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('Mortality per 100K'); title('Females, aged 15-79');
xlim([1980 2060]); ylim([0 4000]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max');
grid on;

% male
subplot(1 , 3 , 2);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'HIV_mortality_males_aged15-79' , '.csv'];
    hivMortM = xlsread(fname);
    if i == 1
        hold all;
        plot(hivMortM(: , 1) , hivMortM(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivMortM(: , 1) , hivMortM(: , 3) , 'LineStyle' ,  '--' , 'Color' , colorList{i});
        plot(hivMortM(: , 1) , hivMortM(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(hivMortM((2019-startYear)+1 : end , 1) , hivMortM((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivMortM((2019-startYear)+1 : end , 1) , hivMortM((2019-startYear)+1 : end , 3) , 'LineStyle' ,  '--' , 'Color' , colorList{i});
        plot(hivMortM((2019-startYear)+1 : end , 1) , hivMortM((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('Mortality per 100K'); title('Males, aged 15-79');
xlim([1980 2060]); ylim([0 4000]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max');
grid on;

% combined
subplot(1 , 3 , 3);
for i = 1 : 3
    fname = [pwd , '\HHCoM_Results\Vaccine' , baseFileNameShort , '_S' , num2str(i) , '_11_1' , '\' , ...
        'HIV_mortality_combined_aged15-79' , '.csv'];
    hivMortC = xlsread(fname);
    if i == 1
        hold all;
        plot(hivMortC(: , 1) , hivMortC(: , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivMortC(: , 1) , hivMortC(: , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivMortC(: , 1) , hivMortC(: , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    else
        hold all;
        plot(hivMortC((2019-startYear)+1 : end , 1) , hivMortC((2019-startYear)+1 : end , 2) , 'LineStyle' , '-' , 'Color' , colorList{i} , 'LineWidth' , 2);
        plot(hivMortC((2019-startYear)+1 : end , 1) , hivMortC((2019-startYear)+1 : end , 3) , 'LineStyle' , '--' , 'Color' , colorList{i});
        plot(hivMortC((2019-startYear)+1 : end , 1) , hivMortC((2019-startYear)+1 : end , 4) , 'LineStyle' , '--' , 'Color' , colorList{i});
    end
end
xlabel('Year'); ylabel('Mortality per 100K'); title('Females + Males, aged 15-79');
xlim([1980 2060]); ylim([0 4000]);
legend('Scenario 1, mean' , 'min' , 'max' , 'Scenario 2, mean' , 'min' , 'max' , ...
    'Scenario 3, mean' , 'min' , 'max');
grid on;
sgtitle('HIV-Associated Mortality');
