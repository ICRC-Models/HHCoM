function[] = showResults_22Apr20calib_bw()

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

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1)

% Indices of calib runs to plot
fileInds = {'12_3346' , '12_2618' , '11_932' , '16_3038' , '8_597' , '12_2550' , ...  % 22Apr20Ph2
    '12_3895' , '22_487' , '8_2705' , '22_3250' , '15_2550' , '14_563' , ...
    '4_1887' , '10_688' , '18_3391' , '14_2659' , '19_2814' , '18_903' , ...
    '22_2697' , '4_1676' , '4_2471' , '15_2517' , '16_1709' , '12_2481' , '16_3992'};
% {'2_846' , '14_947' , '16_3127' , '12_689' , ...  % 22Apr20Ph1, top 50 sets
%     '10_2727' , '17_3986' , '16_2194' , '15_3850' , '9_334' , '0_6657' , ...
%     '16_2364' , '4_711' , '19_1017' , '4_2361' , '15_2155' , '17_594' , ...
%     '19_1779' , '11_1541' , '12_3055' , '6_746' , '20_944' , '13_3012' , ...
%     '18_387' , '17_1649' , '17_3242' , '14_2649' , '16_2701' , '20_1864' , ...
%     '19_3788' , '10_425' , '11_3176' , '18_45' , '15_532' , '20_3201' , ...
%     '5_1822' , '7_2249' , '15_2965' , '19_2304' , '10_2162' , '11_301' , ...
%     '14_1533' , '18_1043' , '0_8616' , '5_4018' , '11_3613' , '18_2578' , ...
%     '7_1869' , '15_1004' , '12_2230' , '8_3465'};
nRuns = length(fileInds);

% Initialize model output plots
% Total population size
popYearVec = unique(totPopSize_dObs(: ,1));
popSize = zeros(nRuns , length(popYearVec));
% Population age distribution
popYearVec = unique(popAgeDist_dObs(: ,1));
popProp = zeros(nRuns , length(popYearVec) , age);
% HIV prevalence
hivYearVec = unique(hivPrevM_dObs(: ,1));
hivAgeM = zeros(nRuns , 7 , length(hivYearVec));
hivAgeF = hivAgeM;
% Female HPV prevalence
hpv_hiv = zeros(nRuns , 9);
hpv_hivNeg = hpv_hiv;
% Male HPV prevalence
hpv_hivM2008 = zeros(nRuns , 4);
hpv_hivMNeg2008 = hpv_hivM2008;
% CIN2/3 prevalence
cinPos2002 = zeros(nRuns , 10);
cinNeg2002 = cinPos2002;
% CC incidence
ccInc2012 = zeros(nRuns , age);
% HPV/CIN/CC type distribution
typeDistYearVec = [2012];
cc_vax = zeros(nRuns , length(typeDistYearVec));
cc_nonVax = cc_vax;
cin3_vax = cc_vax;
cin3_nonVax = cc_vax;
cin1_vax = cc_vax;
cin1_nonVax = cc_vax;
hpv_vax = cc_vax;
hpv_nonVax = cc_vax;


resultsDir = [pwd , '\HHCoM_Results\'];
for j = 1 : nRuns
    % Load results
    pathModifier = ['toNow_22Apr20Ph2_noBaseVax_baseScreen_hpvHIVcalib_' , fileInds{j}];
    load([resultsDir , pathModifier])
   

    %% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

    %% Population size over time vs. UN/SSA data
    popYearVec = unique(totPopSize_dObs(: ,1));
    for t = 1 : length(popYearVec)
        popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
        popSize(j , t) = sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTot),2);
    end

    %% Population size by 5-year age groups over time vs. SSA data
    popYearVec = unique(popAgeDist_dObs(: ,1));
    for t = 1 : length(popYearVec)
        for a = 1 : age
            popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
            popProp(j , t , a) = sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popAge),2) ./ ...
                sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTot),2);
        end
    end

    %% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

    %% HIV prevalence by age over time vs. Africa Centre data
    hivYearVec = unique(hivPrevM_dObs(: ,1));
    
    for t = 1 : length(hivYearVec)
        for a = 4 : 10
            hivMInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            artMInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            totMInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            hivFInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            artFInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            totFInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            hivAgeM(j , a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivMInds)) ...
                + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artMInds))) ...
                / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totMInds));
            hivAgeF(j , a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivFInds)) ...
                + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artFInds))) ...
                / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totFInds));
        end
    end

    %% ********************************** HPV FIGURES **********************************************************************************************

    %% HPV Prevalence by age in 2002 vs. McDonald 2014 data
    yr = 2002;
    for a = 4 : 12 % 15-19 -> 55-65
        hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpv_hiv(j , a - 3) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
            ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
        
        hpvInds_hivNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds_hivNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpv_hivNeg(j , a - 3) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds_hivNeg))...
            ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds_hivNeg));
    end

    %% HPV prevalence by age and HIV status in 2008 vs. Mbulawa data
    ageVec = {[4:5],[6:7],[8:9],[10:13]};
    for aV = 1 : length(ageVec)
        a = ageVec{aV};
        hpvInds_hivM = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 , [1 : 2 , 7] , ...
            1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
        ageInds_hivM = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        hpv_hivM2008(j , aV) = sum(popVec((2008 - startYear) * stepsPerYear +1 , hpvInds_hivM))...
            ./ sum(popVec((2008 - startYear) * stepsPerYear+1 , ageInds_hivM));
        
        hpvInds_hivMNeg = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 , [1 : 2 , 7] , ...
            1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
        ageInds_hivMNeg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        hpv_hivMNeg2008(j , aV) = sum(popVec((2008 - startYear) * stepsPerYear +1 , hpvInds_hivMNeg))...
        ./ sum(popVec((2008 - startYear) * stepsPerYear +1 , ageInds_hivMNeg));
    end

    %% ********************************** CIN FIGURES *********************************************************************************************

    %% CIN2/3 prevalence for All HR HPV types combined by HIV status and age in 2002 vs. McDonald 2014 data
    for a = 4 : 13 % 15-19 -> 60-64
        % HIV-positive (on and not on ART)
        cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        cinPos2002(j , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cinInds)))...
            ./ sum(popVec((2002 - startYear) * stepsPerYear +1 , ageInds));
        % HIV-negative
        cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        cinNeg2002(j , a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cinNegInds)))...
            ./ (sum(popVec((2002 - startYear) * stepsPerYear +1 , ageNegInds)));
    end
    
    %% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

    %% Cervical cancer incidence in 2011 by age vs. Globocan 2012 data and other sources
    incTimeSpan = [((2012 - startYear) * stepsPerYear +1) : ((2012 - startYear) * stepsPerYear +6)];
    fac = 10 ^ 5;

    for a = 1 : age
        allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012(j , a) = ...
            (annlz(sum(sum(sum(newCC(incTimeSpan , : , a , :),2),3),4)) ./ ...
            (annlz(sum(popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);
    end
   
    %% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************
    
    %% HPV type distribution by state over time (coinfections grouped as 9v-type HPV)
    typeDistYearVec = [2012];
    ccInds_vax = toInd(allcomb(1 : disease , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
        1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
    ccInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , 6 , ...
        1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
    ccInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
            1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
        
    cin3Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin3Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 4 , 7] , 5 , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin3Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
        
    cin1Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin1Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 2 , 7] , 3 , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin1Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
   
    hpvInds_vax = toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    hpvInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 2 , 7] , 3 , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    hpvInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

    for i = 1 : length(typeDistYearVec)
        cc_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);
        cc_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);

        cin3_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);
        cin3_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);

        cin1_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_tot) , 2);
        cin1_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_tot) , 2);

        hpv_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_tot) , 2);
        hpv_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_tot) , 2);
    end
    
end

%% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

%% Population size over time vs. UN/SSA data
popYearVec = unique(totPopSize_dObs(: ,1));
% Load calibration data from Excel (years, values)
file = [pwd , '/Config/Population_validation_targets.xlsx'];
historicalPop0_69 = zeros(5,2);
futurePop0_69 = zeros(16,2);
historicalPop0_69(:,1) = xlsread(file , 'Demographics' , 'B91:F91'); % years
historicalPop0_69(:,2) = xlsread(file , 'Demographics' , 'B130:F130') .* 1000; % estimates
futurePop0_69(:,1) = xlsread(file , 'Demographics' , 'C144:R144'); % years
futurePop0_69(:,2) = xlsread(file , 'Demographics' , 'C146:R146') .* 1000; % projections

% Calibration error bars
mean = totPopSize_dObs(: , 2);
sdev = (totPopSize_dObs(: , 3).^(1/2));

figure;
errorbar(totPopSize_dObs(: , 1) , mean , sdev , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot(popSize , 'Positions' , popYearVec , 'Labels' , popYearVec , 'Color' , 'k' , 'Whisker' , 5);
title('KZN Population Size Ages 0-79')
xlabel('Year'); ylabel('Individuals')
xlim([2000 2020]); ylim([0 (14*10^6)]);
legend('Observed KZN ages 0-79: mean, SD');
grid on;

%% Population size by 5-year age groups over time vs. SSA data
% % Load calibration data from Excel
% file = [pwd , '/Config/Population_validation_targets.xlsx'];
% years = xlsread(file , 'Demographics' , 'B91:F91');    % years
% kzn_popByage_yrs(: , :) = xlsread(file , 'Demographics' , 'M92:Q107').*1000;    % males and females by age in 1996-2019
% 
% % Calibration error bars
% mean = [popAgeDist_dObs(1:16 , 2) , popAgeDist_dObs(17:32 , 2) , popAgeDist_dObs(33:48 , 2)]';
% sdev = ([popAgeDist_dObs(1:16 , 3) , popAgeDist_dObs(17:32 , 3) , popAgeDist_dObs(33:48 , 3)]'.^(1/2));
% 
% figure;
% subplot(1,3,1);
% set(gca,'ColorOrderIndex',1)
% plot(years , popPropYrs_obs(: , 1:7) , 'o');
% hold on;
% calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , ...
%     unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) ,];
% errorbar(calibYrs , mean(: , 1:7) , sdev(: , 1:7) , ...
%     's' , 'Color' , [0.9290, 0.6940, 0.1250])
% hold on;
% % insert boxplot
% ylim([0.05 0.15]);
% ylabel('Population proportion by age'); xlabel('Year');
% legend('0-4, Observed' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' , '30-34' , ...
%     'Calibration SD' , 'Location' , 'EastOutside');
% 
% subplot(1,3,2);
% set(gca,'ColorOrderIndex',1)
% hold on;
% plot(years , popPropYrs_obs(: , 8:14) , 'o');
% hold on;
% calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , ...
%     unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) ,];
% errorbar(calibYrs , mean(: , 8:14) , sdev(: , 8:14) , ...
%     's' , 'Color' , [0.9290, 0.6940, 0.1250])
% hold on;
% % insert boxplot
% ylim([0.0 0.1]);
% ylabel('Population proportion by age'); xlabel('Year');
% legend('35-39, Observed' , '40-44' , '45-49' , '50-54' , '55-59' , '60-64' , '65-69' , ...
%     'Calibration SD' , 'Location' , 'EastOutside');
% 
% subplot(1,3,3);
% set(gca,'ColorOrderIndex',1)
% plot(years , popPropYrs_obs(: , 15:16) , 'o');
% hold on;
% calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1))];
% errorbar(calibYrs , mean(: , 15:16) , sdev(: , 15:16) , ...
%     's' , 'Color' , [0.9290, 0.6940, 0.1250])
% hold on;
% % insert boxplot
% ylim([0.0 0.02]);
% ylabel('Population proportion by age'); xlabel('Year'); %title('KZN age distribution in 5-year groups');
% legend('70-74, Observed' , '75-79' , 'Calibration SD' , 'Location' , 'EastOutside');

%% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

%% HIV prevalence by age over time vs. Africa Centre data
ageGroup = {'15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74'};

% 2010-2016 AC data
hivPrevF_val = [9.29	9.02	10.45	9.33	10.37	11.00	9.35
31.41	31.68	30.64	33.95	34.56	34.12	33.42
53.27	51.72	50.80	51.33	51.94	53.98	52.41
59.18	61.35	58.66	64.90	62.57	64.71	63.09
53.97	54.08	58.77	65.12	65.28	64.66	66.95
42.69	43.27	45.29	49.16	54.25	56.37	61.28
32.34	34.30	39.18	41.47	48.21	49.57	50.23];

hivPrevM_val = [1.60	1.85	2.75	3.46	2.87	3.95	4.50
9.56	8.02	9.87	9.65	11.86	7.19	8.02
28.99	21.92	24.88	29.84	35.40	27.65	27.31
46.47	44.51	39.49	47.22	46.35	41.64	42.08
52.03	44.30	49.61	63.33	51.41	52.05	51.35
41.73	41.53	51.55	51.64	59.40	52.69	51.18
36.64	37.12	33.01	40.00	40.54	44.52	52.17];

% Calibration error bars
hivM(: , 1) = hivPrevM_dObs(: , 2) .* 100; % mean
hivM(: , 2) = (hivPrevM_dObs(: , 3).^(1/2)) .* 100; % calibration SD
hivF(: , 1) = hivPrevF_dObs(: , 2) .* 100; % mean
hivF(: , 2) = (hivPrevF_dObs(: , 3).^(1/2)) .* 100; % calibration SD

prevYears = unique(hivPrevF_dObs(: , 1));
prevYears2 = [2010 : 2016];

gen = {'Male' , 'Female'};
for g = 1 : gender
    hivPrevs = hivM;
    hivPrevs2 = hivPrevM_val;
    hivModel = hivAgeM;
    if g == 2
        hivPrevs = hivF;
        hivPrevs2 = hivPrevF_val;
        hivModel = hivAgeF;
    end

    figure; 
    for a = 4 : 10
        subplot(3 , 3 , a-3)
        hold all;
        if a <= 11            
            errorbar(prevYears , hivPrevs(((a-3) : 7 : end) , 1) , hivPrevs(((a-3) : 7 : end) , 2) , ...
                'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
            hold on;
            boxplot((squeeze(hivModel(: , a-3 , :)) .* 100) , ...
                'Positions' , prevYears , 'Labels' , prevYears , 'Color' , 'k' , 'Whisker' , 5);
        else
            boxplot((squeeze(hivModel(: , a-3 , :)) .* 100) , ...
                'Positions' , prevYears , 'Labels' , prevYears , 'Color' , 'k' , 'Whisker' , 5);
        end
        xlabel('Year'); ylabel('Prevalence (%)'); title([gen{g} , 's ages ' , ageGroup{a-3}]) % , ' HIV Prevalence'])
        xlim([2000 2010]); ylim([0 60]);
        grid on;
    end
    legend('Observed KZN: mean, SD')
end

%% ********************************** HPV FIGURES **********************************************************************************************

%% HPV Prevalence by age in 2002 vs. McDonald 2014 data
ageGroup = {'17 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};

% Calibration error bars
mean = hpv_hiv_dObs(: , 2) .* 100;
sdev = (hpv_hiv_dObs(: , 3).^(1/2)) .* 100;
meanNeg = hpv_hivNeg_dObs(: , 2) .* 100;
sdevNeg = (hpv_hivNeg_dObs(: , 3).^(1/2)) .* 100;

figure;
subplot(2,1,1);
errorbar(1 : length(mean) , mean , sdev , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot((hpv_hiv .* 100) , 'Color' , 'k' , 'Whisker' , 5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('Observed HIV-Positive Females: mean, SD');
title('HPV Prevalence in 2002 - Females, HIV+');

subplot(2,1,2);
errorbar(1 : length(mean) , meanNeg , sdevNeg , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot((hpv_hivNeg .* 100) , 'Color' , 'k' , 'Whisker' , 5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('Observed HIV-Negative Females: mean, SD');
title('HPV Prevalence in 2002 - Females, HIV-');
grid on;

%% HPV prevalence by age and HIV status in 2008 vs. Mbulawa data
ageGroup = {'15-24' , '25-34' , '35-44' , '45-64'};

% Calibration error bars
mean = hpv_hivM2008_dObs(: , 2) .* 100;
sdev = (hpv_hivM2008_dObs(: , 3).^(1/2)) .* 100;
meanNeg = hpv_hivMNeg2008_dObs(: , 2) .* 100;
sdevNeg = (hpv_hivMNeg2008_dObs(: , 3).^(1/2)) .* 100;

figure;
subplot(2,1,1)
errorbar(1 : length(mean) , mean , sdev , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot((hpv_hivM2008 .* 100) , 'Color' , 'k' , 'Whisker' , 5)
set(gca , 'xtick' , [1 : length(ageGroup)] , 'xtickLabel' , ageGroup);
legend('Observed HIV-Positive Males: mean, SD');
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)'); ylim([0 100]);
title('HPV Prevalence in 2008 - Males, HIV+');
grid on;

subplot(2,1,2)
errorbar(1 : length(mean) , meanNeg , sdevNeg , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot((hpv_hivMNeg2008 .* 100) , 'Color' , 'k' , 'Whisker' , 5)
set(gca , 'xtick' , [1 : length(ageGroup)] , 'xtickLabel' , ageGroup);
legend('Observed HIV-Negative Males: mean, SD');
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)'); ylim([0 100]);
title('HPV Prevalence in 2008 - Males, HIV-');
grid on;

%% ********************************** CIN FIGURES *********************************************************************************************

%% CIN2/3 prevalence for All HR HPV types combined by HIV status and age in 2002 vs. McDonald 2014 data
ageGroup = {'17-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

% Calibration error bars
cinPos_mean = cinPos2002_dObs(: , 2) .* 100;
cinPos_sdev = (cinPos2002_dObs(: , 3).^(1/2)) .* 100;
cinNeg_mean = cinNeg2002_dObs(: , 2) .* 100;
cinNeg_sdev = (cinNeg2002_dObs(: , 3).^(1/2)) .* 100;

figure;
subplot(2 , 1 , 1);
errorbar(1 : length(cinPos_mean) , cinPos_mean , cinPos_sdev , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot((cinPos2002 .* 100) , 'Color' , 'k' , 'Whisker' , 5)
legend('Observed HIV-positive Females: mean, SD');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('CIN 2/3 Prevalence in 2002 - HIV+')
ylim([0 25])
grid on;

subplot(2 , 1 , 2)
errorbar(1 : length(cinNeg_mean) , cinNeg_mean , cinNeg_sdev , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot((cinNeg2002 .* 100) , 'Color' , 'k' , 'Whisker' , 5);
legend('Observed HIV-negative Females: mean, SD');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('CIN 2/3 Prevalence in 2002 - HIV-')
ylim([0 10])
grid on;

%% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

%% Cervical cancer incidence in 2011 by age vs. Globocan 2012 data and other sources
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

globocan = [0.00
2.646467154
8.848389036
45.1937379
53.40682334
63.4
68.3
70.7
73
77.4
82.7
88.6
95.2];

combined_ub = [0.00
2.65
8.85
45.19
53.41
67.05
80.83
78.97
128.87
105.27
118.70
111.81
95.20];

combined_lb = [0.00
0.00
0.41
9.97
12.61
25.00
45.69
36.31
50.55
57.08
62.69
42.43
52.01];

% Calibration error bars
mean = ccInc2012_dObs(: , 2);
sdev = (ccInc2012_dObs(: , 3).^(1/2));

figure;    
% Plot observed data
plot(4 : age , combined_ub , 'r-' , 4 : age , combined_lb , 'r-'); % , ...
    %'Color' , [0.9290, 0.6940, 0.1250]);
hold on; 
errorbar(4 : age , mean , sdev , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot(ccInc2012 , 'Color' , 'k' , 'Whisker' , 5);
xlabel('Age Group'); ylabel('Incidence per 100,000');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
ylim([0 160]);
title(['Cervical Cancer Incidence in 2012']);
legend('Combined SA: upper bound' , 'Combined SA: lower bound' , ...
    'Observed Females: mean, SD');
grid on;

%% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************

%% HPV type distribution by state over time (coinfections grouped as 9v-type HPV)

% HPV infected
% Calibration error bars
mean = hpv_dist_dObs(: , 2) .* 100;
sdev = (hpv_dist_dObs(: , 3).^(1/2)) .* 100;

figure;
subplot(2,3,1)
plot(2012 , 46.82 , 'k*');
hold on;
plot(2012 , 53.18 , 'b*');
hold on;
errorbar([2012, 2012] , mean , sdev , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot(hpv_vax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
    'Colors' , 'k' , 'Whisker' , 5);
hold on;
boxplot(hpv_nonVax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
    'Colors' , 'b' , 'Whisker' , 5);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('HPV');
ylim([0 100]);
xlim([2010 2015]);
grid on;

% CIN1
% Calibration error bars
mean = cin1_dist_dObs(: , 2) .* 100;
sdev = (cin1_dist_dObs(: , 3).^(1/2)) .* 100;

subplot(2,3,2)
plot(2012 , 51.92 , 'k*');
hold on;
plot(2012 , 48.08 , 'b*');
hold on;
errorbar([2012, 2012] , mean , sdev , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot(cin1_vax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
    'Colors' , 'k' , 'Whisker' , 5);
hold on;
boxplot(cin1_nonVax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
    'Colors' , 'b' , 'Whisker' , 5);
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('CIN1');
grid on;

% CIN2
subplot(2,3,3);
plot(2012 , 62.81 , 'k*');
hold on;
plot(2012 , 37.19 , 'b*');
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('CIN2');
grid on;

% CIN3
% Calibration error bars
mean = cin3_dist_dObs(: , 2) .* 100;
sdev = (cin3_dist_dObs(: , 3).^(1/2)) .* 100;

subplot(2,3,4)
plot(2012 , 73.71 , 'k*');
hold on;
plot(2012 , 26.29 , 'b*');
hold on;
errorbar([2012, 2012] , mean , sdev , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot(cin3_vax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
    'Colors' , 'k' , 'Whisker' , 5);
hold on;
boxplot(cin3_nonVax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
    'Colors' , 'b' , 'Whisker' , 5);
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('CIN3');
grid on;

% CC
% Calibration error bars
mean = cc_dist_dObs(: , 2) .* 100;
sdev = (cc_dist_dObs(: , 3).^(1/2)) .* 100;

subplot(2,3,5)
plot(2012 , 85.78 , 'k*');
hold on;
plot(2012 , 14.22 , 'b*');
hold on;
errorbar([2012, 2012] , mean , sdev , ...
    'rs'); % , 'Color' , [0.9290, 0.6940, 0.1250])
hold on;
boxplot(cc_vax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
    'Colors' , 'k' , 'Whisker' , 5);
hold on;
boxplot(cc_nonVax .* 100 , 'Positions' , 2012 , 'Labels' , num2str(2012) , ...
    'Colors' , 'b' , 'Whisker' , 5);
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('Cervical Cancer')
legend('9v-type HPV' , 'Non-9v-type HPV' , 'Observed 2011: 9v' , 'Observed 2011: non-9v');
grid on;
