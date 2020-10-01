% function[] = showResults_multSims_CIs()

%% Load parameters and results
paramDir = [pwd , '\Params\'];

[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd, annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , partnersMmult, maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4,...
    mue , mue2 , mue3 , mue4 , epsA_vec , epsR_vec , yr , ...
    hivOn , betaHIV_mod ,  muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , muCC , kRL , kDR , artHpvMult , ...
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
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ageInd , riskInd , ...
    hivNegNonVMMCinds , hivNegVMMCinds , vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3  , dFertMat3, d_partnersMmult, riskAdj, d_riskAdj, ...
    deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2(1 , 0 , [] , [] , []);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)

% Indices of calib runs to plot
fileInds = {''};    % 22Apr20Ph2V11
nRuns = length(fileInds);

% Initialize model output plots
% Total population size
popYearVec = unique(totPopSize_dObs(: ,1));
popSize = zeros(nRuns , length(popYearVec));
% Population age distribution
popYearVecAge = unique(popAgeDist_dObs(: ,1));
popProp = zeros(nRuns , length(popYearVecAge) , age);
% HIV prevalence
hivYearVec = unique(hivPrevM_dObs(: ,1));
hivAgeM = zeros(nRuns , 8 , length(hivYearVec));
hivAgeF = zeros(nRuns , 7 , length(hivYearVec));
hivAgeAll = zeros(nRuns , 10 , 1);
% HPV Prevalence in high risk women
hpv_hiv = zeros(nRuns , 4);
hpv_hivNeg = zeros(nRuns , 4);
% HPV prevalence in all women (not including CIN2+)
hpv_all = zeros(nRuns , 6);
% HPV prevalence in HIV+ women (including CIN)
hpv_hiv2009 = zeros(nRuns , 5);
% CIN1 prevalence among women aged 20-50 in 2010 by HIV status 
cin1_2010 = zeros(nRuns , 2);
% CIN2/3 prevalence among women aged 20-50 in 2010 by HIV status 
cin2_2010 = zeros(nRuns , 2);
% CIN2/3 prevalence among HIV+ women aged 20-50 in 2007
cinPos2007 = zeros(nRuns , 3);
% CC incidence
ccInc2012 = zeros(nRuns , age);
ccInc2012neg = ccInc2012;
ccInc2012pos = ccInc2012;
ccInc2012art = ccInc2012;
% HPV/CIN/CC type distribution
typeDistYearVec = [2010 , 2011 , 2012 , 2013 , 2014 , 2015];
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
    pathModifier = ['toNow_24Aug20_K_HPVprogAge' , fileInds{j}];
    load([resultsDir , pathModifier])
   
    %% ***************************** DEMOGRAPHY FIGURES **********************************************************************************************

    %% Population size over time vs. Kenya Population Census data (calibration)
    for t = 1 : length(popYearVec)
        popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
        popSize(j , t) = sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTot),2);
    end

    %% Population size by 5-year age groups over time vs. Kenya Population Census data (calibration)
    for t = 1 : length(popYearVecAge)
        for a = 1 : age
            popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
            popProp(j , t , a) = sum(popVec(((popYearVecAge(t) - startYear) * stepsPerYear +1) , popAge),2) ./ ...
                sum(popVec(((popYearVecAge(t) - startYear) * stepsPerYear +1) , popTot),2);
        end
    end

    %% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

    %% HIV prevalence by age and gender over time vs. KAIS and DHS (calibration)
    for t = 1 : length(hivYearVec)
        for a = 4 : 10 
            hivFInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            artFInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            totFInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
            hivAgeF(j , a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivFInds)) ...
                + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artFInds))) ...
                / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totFInds));
        end
    end
    
    for t = 1 : length(hivYearVec)
        for a = 4 : 11
            hivMInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            artMInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            totMInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
            hivAgeM(j , a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivMInds)) ...
                + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artMInds))) ...
                / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totMInds));
        end
    end
    
    for t = 1 : length(hivYearVec)
        for a = 4 : 13
            hivAllInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            artAllInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            totAllInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
            hivAgeAll(j , a - 3 , 1) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivAllInds)) ...
                + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artAllInds))) ...
                / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totAllInds));
        end
    end
    
    %% ********************************** HPV FIGURES **********************************************************************************************
    %% HPV Prevalence in high risk women
    yr = 2006;
    ageVec = {[4:5],[6],[7:8],[9:10]};
    for aV = 1 : length(ageVec)
        a = ageVec{aV};
        %HIV-positive
        hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 1 , a , 3)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 1 , a , 3))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 3));
        hpv_hiv(j , aV) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
            ./ sum(popVec((yr - startYear) * stepsPerYear+1 , ageInds));
        %HIV-negative
        hpvInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 1 , a , 3)); toInd(allcomb(1 : 2 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5, 1 , 1 : intervens , 1 , a , 3))]);
        ageInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 3));
        hpv_hivNeg(j , aV) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
            ./ sum(popVec((yr - startYear) * stepsPerYear+1 , ageInds));
    end

    %% HPV prevalence in all women (not including CIN2+) by age
    yr = 2000;
    for a = 6 : 11
        hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 3 , 7] , 2 : 3 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(1 : disease, 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpv_all(j , a - 5) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
            ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    end

    %% HPV prevalence in HIV+ women (including CIN)
    yr = 2009;
    ageVec = {[4:6], 7, 8, 9, [10:11]};
    for aV = 1 : length(ageVec)
        a = ageVec{aV}; % 15-19 ->  55-65
        hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hpv_hiv2009(j , aV) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
            ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    end

    %% ********************************** CIN FIGURES *********************************************************************************************
    %% CIN2/3 prevalence among HIV+ women aged 20-50 in 2007
    yr = 2007;
    ageVec = {[5:6],[7:8],[9:10]};
    for aV = 1 : length(ageVec)
        a = ageVec{aV};
        % HIV-positive (on and not on ART)
        cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
        ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        cinPos2007(j , aV) = (sum(popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
            ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    end

    %% CIN1 prevalence among women aged 20-50 in 2010 by HIV status 
    yr = 2010;
    dVec = {[1:2],[3:8]};
    for dV = 1 : length(dVec)
        d = dVec{dV};
        cinInds = unique([toInd(allcomb(d , 1 : viral , 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , 5 : 10 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
            [1 : 3 , 7] , 3 , 1 , 1 : intervens , 2 , 5 : 10 , 1 : risk))]);
        ageInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 5 : 10 , 1 : risk));
        cin1_2010(j , dV) = (sum(popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
            ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    end

    %% CIN2/3 prevalence among women aged 20-50 in 2010 by HIV status 
    yr = 2010;
    dVec = {[1:2],[3:8]};
    for dV = 1 : length(dVec)
        d = dVec{dV};
        cinInds = unique([toInd(allcomb(d , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , 5 : 10 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
            [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , 5 : 10 , 1 : risk))]);
        ageInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , 5 : 10 , 1 : risk));
        cin2_2010(j , dV) = (sum(popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
            ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
    end
    

    %% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

    %% Cervical cancer incidence in 2012 by age vs. Globocan 2012 data and other sources (calibration)
    incTimeSpan = [((2012 - startYear) * stepsPerYear +1) : ((2012 - startYear) * stepsPerYear +6)];
    fac = 10 ^ 5;

    for a = 1 : age
        % General population
        allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012(j , a) = ...
            (annlz(sum(sum(sum(newCC(incTimeSpan , : , a , :),2),3),4)) ./ ...
            (annlz(sum(popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);
        
        % HIV-negative
        allFneg = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012neg(j , a) = ...
            (annlz(sum(sum(sum(newCC(incTimeSpan , 1 : 2 , a , :),2),3),4)) ./ ...
            (annlz(sum(popVec(incTimeSpan , allFneg) , 2) ./ stepsPerYear)) * fac);
        
        % HIV-positive untreated
        allFpos = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012pos(j , a) = ...
            (annlz(sum(sum(sum(newCC(incTimeSpan , 3 : 7 , a , :),2),3),4)) ./ ...
            (annlz(sum(popVec(incTimeSpan , allFpos) , 2) ./ stepsPerYear)) * fac);
        
        % HIV-positive on ART
        allFart = toInd(allcomb(8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        % Calculate incidence
        ccInc2012art(j , a) = ...
            (annlz(sum(sum(sum(newCC(incTimeSpan , 8 , a , :),2),3),4)) ./ ...
            (annlz(sum(popVec(incTimeSpan , allFart) , 2) ./ stepsPerYear)) * fac);
    end

   
    %% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************
    
    %% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)
    typeDistYearVec = [2010 , 2011 , 2012 , 2013 , 2014 , 2015];
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
        
    cin2Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 4 , [1 : 4 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin2Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 3 , 7] , 4 , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin2Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 4 , [1 : 4 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 3 , 7] , 4 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
        
    cin1Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin1Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 2 , 7] , 3 , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    cin1Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
   
    hpvInds_vax = toInd(allcomb(1 : disease , 1 : viral , 2 , [1 : 2 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    hpvInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 , 7] , 2 , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk));
    hpvInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 2 , [1 : 2 , 7] , ...
            1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
            [1 , 7] , 2 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

    for i = 1 : length(typeDistYearVec)
        cc_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);
        cc_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);

        cin3_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);
        cin3_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);
        
        cin2_vax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin2Inds_vax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin2Inds_tot) , 2);
        cin2_nonVax(j , i) = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin2Inds_nonVax) , 2)...
            ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin2Inds_tot) , 2);

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

%% Population size over time vs. Kenya Population Census data (calibration)
% Calibration error bars
meanObs = totPopSize_dObs(: , 2);
sdevObs = (totPopSize_dObs(: , 3).^(1/2)).*2;

figure;
errorbar(totPopSize_dObs(: , 1) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold all;
plot(totPopSize_dObs(: , 1) , mean(popSize,1)' , 'k-' , ...
    totPopSize_dObs(: , 1) , min(popSize,[],1)' , 'k--' , ...
    totPopSize_dObs(: , 1) , max(popSize,[],1)' , 'k--' , 'LineWidth' , 1.5);
title('Kenya Population Size Ages 0-79')
xlabel('Year'); ylabel('Individuals')
xlim([1970 2020]); ylim([0 (4*10^7)]);
legend('(Kenya Population Census) Observed Kenya, ages 0-79: mean, 2SD' , ...
    ['Model, ages 0-79: ' , num2str(nRuns) , '-sets mean'] , ['Model, ages 0-79: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model, ages 0-79: ' , num2str(nRuns) , '-sets maximum'] , 'Location' , 'NorthWest');
grid on;

%% Population size by 5-year age groups over time vs. Kenya Population Census data (calibration)
% Calibration error bars
meanObs = [popAgeDist_dObs(1:16 , 2) , popAgeDist_dObs(17:32 , 2) , popAgeDist_dObs(33:48 , 2) , popAgeDist_dObs(49:64 , 2)]';
sdevObs = ([popAgeDist_dObs(1:16 , 3) , popAgeDist_dObs(17:32 , 3) , popAgeDist_dObs(33:48 , 3) , popAgeDist_dObs(49:64 , 3)]'.^(1/2)).*2;

figure;
subplot(1,3,1);
set(gca,'ColorOrderIndex',1)
calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , ...
    unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) ,];
errorbar(calibYrs , meanObs(: , 1:7) , sdevObs(: , 1:7) , ...
    's' , 'LineWidth' , 1.5);
hold on;
for a = 1 : 7
    set(gca,'ColorOrderIndex',a)
    plot(popYearVecAge , mean(popProp(: , : , a),1)' , '-' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a)
    plot(popYearVecAge , min(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a)
    plot(popYearVecAge , max(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
end
ylim([0.05 0.25]);
ylabel('Population proportion by age'); xlabel('Year');
legend('(Kenya Population Census) Observed Kenya, ages 0-4: mean, 2SD' , ...
    'ages 5-9: mean, 2SD' , 'ages 10-14: mean, 2SD' , 'ages 15-19: mean, 2SD' , ...
    'ages 20-24: mean, 2SD' , 'ages 25-29: mean, 2SD' , 'ages 30-34: mean, 2SD' , ...
    ['Model, ages 0-4: ' , num2str(nRuns) , '-sets mean'] , ...
    ['ages 0-4: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['ages 0-4: ' , num2str(nRuns) , '-sets maximum'] , ...
    '...' , 'Location' , 'north');

subplot(1,3,2);
set(gca,'ColorOrderIndex',1)
calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , ...
    unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1)) ,];
errorbar(calibYrs , meanObs(: , 8:14) , sdevObs(: , 8:14) , ...
    's' , 'LineWidth' , 1.5);
hold on;
for a = 8 : 14
    set(gca,'ColorOrderIndex',a-7)
    plot(popYearVecAge , mean(popProp(: , : , a),1)' , '-' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-7)
    plot(popYearVecAge , min(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-7)
    plot(popYearVecAge , max(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
end
ylim([0.0 0.12]);
ylabel('Population proportion by age'); xlabel('Year');
legend('(Kenya Population Census) Observed Kenya, ages 35-39: mean, 2SD' , ...
    'ages 40-44: mean, 2SD' , 'ages 45-49: mean, 2SD' , 'ages 50-54: mean, 2SD' , ...
    'ages 55-59: mean, 2SD' , 'ages 60-64: mean, 2SD' , 'ages 65-69: mean, 2SD' , ...
    ['Model, ages 35-39: ' , num2str(nRuns) , '-sets mean'] , ...
    ['ages 35-39: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['ages 35-39: ' , num2str(nRuns) , '-sets maximum'] , '...' , 'Location' , 'north');


subplot(1,3,3);
set(gca,'ColorOrderIndex',1)
calibYrs = [unique(popAgeDist_dObs(: , 1)) , unique(popAgeDist_dObs(: , 1))];
errorbar(calibYrs , meanObs(: , 15:16) , sdevObs(: , 15:16) , ...
    's' , 'LineWidth' , 1.5);
hold on;
for a = 15 : 16
    set(gca,'ColorOrderIndex',a-14)
    plot(popYearVecAge , mean(popProp(: , : , a),1)' , '-' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-14)
    plot(popYearVecAge , min(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
    set(gca,'ColorOrderIndex',a-14)
    plot(popYearVecAge , max(popProp(: , : , a),[],1)' , '--' , 'LineWidth' , 1.5);
    hold all;
end
ylim([0.0 0.04]);
ylabel('Population proportion by age'); xlabel('Year'); %title('KZN age distribution in 5-year groups');
legend('(Kenya Population Census) Observed Kenya, ages 70-74: mean, 2SD' , ...
    'ages 75-79: mean, 2SD' , ['Model, ages 70-74: ' , num2str(nRuns) , '-sets mean'] , ...
    ['ages 70-74: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['ages 70-74: ' , num2str(nRuns) , '-sets maximum'] , '...' , 'Location' , 'north');

%% ***************************** HIV AND HIV TREATMENT FIGURES ******************************************************************************

%% HIV prevalence by age and gender over time vs. DHS and KAIS (calibration)
% Calibration error bars
hivM(: , 1) = hivPrevM_dObs(: , 2) .* 100; % mean
hivM(: , 2) = (hivPrevM_dObs(: , 3).^(1/2)).*2 .* 100; % calibration SD
hivF(: , 1) = hivPrevF_dObs(: , 2) .* 100; % mean
hivF(: , 2) = (hivPrevF_dObs(: , 3).^(1/2)).*2 .* 100; % calibration SD
hivAll(: , 1) = hivPrevAll_dObs(: , 2) .* 100; % mean
hivAll(: , 2) = (hivPrevAll_dObs(: , 3).^(1/2)).*2 .* 100; % calibration SD

ageGroup = {'15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74'};
ageEndVec = {11-3 , 10-3 , 13-3};
genVecLabels = {'Male' , 'Female' , 'Male+Female Combined'};
genVec = {1 , 2 , 1:2};

figure;
for gInd = 1 : length(genVec)
    g = genVec{gInd};
    if gInd == 1
        hivPrevs = hivM;
        hivModel = hivAgeM;
    elseif gInd == 2
        hivPrevs = hivF;
        hivModel = hivAgeF;
    elseif gInd == 3
        hivPrevs = hivAll;
        hivModel = hivAgeAll;
    end

    subplot(3 , 1 , gInd);
    for i = 1 : size(hivModel , 3)
        if i == 1
            hold all;       
            errorbar(1 : ageEndVec{gInd} , hivPrevs(((i-1)*ageEndVec{gInd})+1 : ageEndVec{gInd}*i , 1) , ...
                hivPrevs(((i-1)*ageEndVec{gInd})+1 : ageEndVec{gInd}*i , 2) , ...
                'rs' , 'LineWidth' , 1.5);
        end
        hold all;
        set(gca,'ColorOrderIndex',i)
        plot(1 : ageEndVec{gInd} , mean((squeeze(hivModel(: , : , i)) .* 100),1)' , '-' , 'LineWidth' , 1.5);
        set(gca,'ColorOrderIndex',i)
        plot(1 : ageEndVec{gInd} , min((squeeze(hivModel(: , : , i)) .* 100),[],1)' , '--' , 'LineWidth' , 1.5);
        set(gca,'ColorOrderIndex',i)
        plot(1 : ageEndVec{gInd} , max((squeeze(hivModel(: , : , i)) .* 100),[],1)' , '--' , 'LineWidth' , 1.5);
        set(gca , 'xtickLabel' , ageGroup); set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
    end
    xlabel('Age Group'); ylabel('HIV Prevalence (%)'); title(genVecLabels{gInd});
    ylim([0 100]);
    grid on;
    if (gInd == 1) || (gInd == 2)
        legend('(DHS/KAIS) Observed Kenya: mean, 2SD' , ...
            ['Model, 2003: ' , num2str(nRuns) , '-sets mean'] , ...
            ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
            ['Model: ' , num2str(nRuns) , '-sets maximum'] , ...
            ['Model, 2007: ' , num2str(nRuns) , '-sets mean'] , ...
            ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
            ['Model: ' , num2str(nRuns) , '-sets maximum'] , ...
            ['Model, 2009: ' , num2str(nRuns) , '-sets mean'] , ...
            ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
            ['Model: ' , num2str(nRuns) , '-sets maximum']);
    elseif gInd == 3
        legend('(DHS/KAIS) Observed Kenya: mean, 2SD' , ...
            ['Model, 2012: ' , num2str(nRuns) , '-sets mean'] , ...
            ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
            ['Model: ' , num2str(nRuns) , '-sets maximum']);
    end
end


%% ********************************** HPV FIGURES **********************************************************************************************

%% HPV Prevalence in high risk women
ageGroup = {'15 -24' , '25 - 29' , '30 - 39' , '40 - 49'};

% Calibration error bars
meanObs = hpv_hiv_dObs(: , 2) .* 100;
sdevObs = (hpv_hiv_dObs(: , 3).^(1/2)).*2 .* 100;
meanNeg = hpv_hivNeg_dObs(: , 2) .* 100;
sdevNeg = (hpv_hivNeg_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
subplot(2,1,1);
errorbar(1 : length(meanObs) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(meanObs) , mean((hpv_hiv .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_hiv .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_hiv .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('(Luchter, 2010) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
title('hrHPV Prevalence in 2006 - High-risk females, HIV+');
grid on;

subplot(2,1,2);
errorbar(1 : length(meanObs) , meanNeg , sdevNeg , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(meanObs) , mean((hpv_hivNeg .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_hivNeg .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_hivNeg .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('(Luchter, 2010) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
title('hrHPV Prevalence in 2006 - High-risk females, HIV-');
grid on;

%% HPV prevalence in all women (not including CIN2+) by age
ageGroup = {'25 - 29' , '30 - 34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54'};

% Calibration error bars
meanObs = hpv_all_dObs(: , 2) .* 100;
sdevObs = (hpv_all_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
errorbar(1 : length(meanObs) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(meanObs) , mean((hpv_all .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_all .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_all .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('(DeVuyst, 2003) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
title('hrHPV Prevalence in 2000 - All females, not including CIN2+');
grid on;

%% HPV prevalence in HIV+ women (including CIN)
ageGroup = {'15 - 29' , '30 - 34' , '35 - 39' , '40 - 44' , '45 - 54'};

% Calibration error bars
meanObs = hpv_hiv2009_dObs(: , 2) .* 100;
sdevObs = (hpv_hiv2009_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
errorbar(1 : length(meanObs) , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(meanObs) , mean((hpv_hiv2009 .* 100),1)' , 'k-' , ...
    1 : length(meanObs) , min((hpv_hiv2009 .* 100),[],1)' , 'k--' , ...
    1 : length(meanObs) , max((hpv_hiv2009 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
set(gca , 'xtickLabel' , ageGroup);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('hrHPV Prevalence (%)');
ylim([0 100]);
legend('(DeVuyst, 2012; BJC) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
title('hrHPV Prevalence in 2009 - Females, HIV+');
grid on;


%% ********************************** CIN FIGURES *********************************************************************************************

%% CIN2/3 prevalence among HIV+ women aged 20-50 in 2007
ageGroup = {'20 - 29' , '30 - 39' , '40 - 49'};

% Calibration error bars
cinPos_mean = cinPos2007_dObs(: , 2) .* 100;
cinPos_sdev = (cinPos2007_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
errorbar(1 : length(cinPos_mean) , cinPos_mean , cinPos_sdev , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(cinPos_mean) , mean((cinPos2007 .* 100),1)' , 'k-' , ...
    1 : length(cinPos_mean) , min((cinPos2007 .* 100),[],1)' , 'k--' , ...
    1 : length(cinPos_mean) , max((cinPos2007 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
legend('(Hutchko) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('CIN 2/3 Prevalence (%)')
title('CIN 2/3 Prevalence in 2007 - Females, HIV+')
ylim([0 25])
grid on;

%% CIN1 prevalence among women aged 20-50 in 2010 by HIV status 
diseaseGroup = {'HIV-negative' , 'HIV-positive'};

% Calibration error bars
cin1_mean = cin1_2010_dObs(: , 2) .* 100;
cin1_sdev = (cin1_2010_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
errorbar(1 : length(cin1_mean) , cin1_mean , cin1_sdev , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(cin1_mean) , mean((cin1_2010 .* 100),1)' , 'k-' , ...
    1 : length(cin1_mean) , min((cin1_2010 .* 100),[],1)' , 'k--' , ...
    1 : length(cin1_mean) , max((cin1_2010 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
legend('(Sweet, 2020) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
set(gca , 'xtick' , 1 : length(diseaseGroup) , 'xtickLabel' , diseaseGroup);
xlabel('HIV Status'); ylabel('CIN1 Prevalence (%)')
title('CIN1 Prevalence in 2010 - Females aged 20-50')
ylim([0 25])
grid on;

%% CIN2/3 prevalence among women aged 20-50 in 2010 by HIV status 
diseaseGroup = {'HIV-negative' , 'HIV-positive'};

% Calibration error bars
cin2_mean = cin2_2010_dObs(: , 2) .* 100;
cin2_sdev = (cin2_2010_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
errorbar(1 : length(cin2_mean) , cin2_mean , cin2_sdev , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(1 : length(cin2_mean) , mean((cin2_2010 .* 100),1)' , 'k-' , ...
    1 : length(cin2_mean) , min((cin2_2010 .* 100),[],1)' , 'k--' , ...
    1 : length(cin2_mean) , max((cin2_2010 .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
legend('(Sweet, 2020) Observed Kenya: mean, 2SD' , ...
    ['Model: ' , num2str(nRuns) , '-sets mean'] , ...
    ['Model: ' , num2str(nRuns) , '-sets minimum'] , ...
    ['Model: ' , num2str(nRuns) , '-sets maximum']);
set(gca , 'xtick' , 1 : length(diseaseGroup) , 'xtickLabel' , diseaseGroup);
xlabel('HIV Status'); ylabel('CIN2/3 Prevalence (%)')
title('CIN2/3 Prevalence in 2010 - Females aged 20-50')
ylim([0 25])
grid on;


%% ****************************** CERVICAL CANCER FIGURES ****************************************************************************************

%% Cervical cancer incidence in 2012 by age vs. Globocan 2012 data and other sources (calibration)
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79'};

% Calibration error bars
meanObs = ccInc2012_dObs(: , 2);
sdevObs = (ccInc2012_dObs(: , 3).^(1/2)).*2;

figure;    
% Plot observed data
errorbar(4 : age , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
% General
plot(1 : age , mean(ccInc2012,1)' , 'k-' , ...
    1 : age , min(ccInc2012,[],1)' , 'k--' , ...
    1 : age , max(ccInc2012,[],1)' , 'k--' , 'LineWidth' , 1.5);
hold all;
% HIV-negative
plot(1 : age , mean(ccInc2012neg,1)' , 'g-' , ...
    1 : age , min(ccInc2012neg,[],1)' , 'g--' , ...
    1 : age , max(ccInc2012neg,[],1)' , 'g--' , 'LineWidth' , 1.5);
hold all;
% HIV-positive untreated
plot(1 : age , mean(ccInc2012pos,1)' , 'm-' , ...
    1 : age , min(ccInc2012pos,[],1)' , 'm--' , ...
    1 : age , max(ccInc2012pos,[],1)' , 'm--' , 'LineWidth' , 1.5);
hold all;
% HIV-positive on ART
plot(1 : age , mean(ccInc2012art,1)' , 'c-' , ...
    1 : age , min(ccInc2012art,[],1)' , 'c--' , ...
    1 : age , max(ccInc2012art,[],1)' , 'c--' , 'LineWidth' , 1.5);
hold all;
xlabel('Age Group'); ylabel('Cervical cancer incidence per 100K');
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
ylim([0 450]);
title(['Cervical Cancer Incidence in 2012']);
legend('Combined SA: upper bound' , 'Combined SA: lower bound' , ...
    '(Globocan, 2012) Observed SA: mean, 2SD' , 'Model, general: ' , num2str(nRuns) , '-sets mean' , 'Model: ' , num2str(nRuns) , '-sets minimum' , 'Model: ' , num2str(nRuns) , '-sets maximum' , ...
    'Model, HIV-negative: ' , num2str(nRuns) , '-sets mean' , 'Model: ' , num2str(nRuns) , '-sets minimum' , 'Model: ' , num2str(nRuns) , '-sets maximum' , ...
    'Model, WLWHIV untreated: ' , num2str(nRuns) , '-sets mean' , 'Model: ' , num2str(nRuns) , '-sets minimum' , 'Model: ' , num2str(nRuns) , '-sets maximum' , ...
    'Model, WLWHIV on ART: ' , num2str(nRuns) , '-sets mean' , 'Model: ' , num2str(nRuns) , '-sets minimum' , 'Model: ' , num2str(nRuns) , '-sets maximum' , ...
    'Location' , 'Northwest');
grid on;


%% ************************** HPV/CIN/CC TYPE DISTRIBUTION FIGURES *******************************************************************************

%% HPV type distribution by state over time (coinfections grouped as 9v-type HPV) (calibration)

% HPV infected
% Calibration error bars
meanObs = hpv_dist_dObs(: , 2) .* 100;
sdevObs = (hpv_dist_dObs(: , 3).^(1/2)).*2 .* 100;

figure;
subplot(2,3,1)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 46.82 , 'k*');
hold on;
plot(2012 , 53.18 , 'b*');
hold on;
plot(typeDistYearVec , mean((hpv_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((hpv_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((hpv_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(typeDistYearVec , mean((hpv_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((hpv_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((hpv_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('HPV');
ylim([0 100]);
xlim([2010 2015]);
grid on;

% CIN1
% Calibration error bars
meanObs = cin1_dist_dObs(: , 2) .* 100;
sdevObs = (cin1_dist_dObs(: , 3).^(1/2)).*2 .* 100;

subplot(2,3,2)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 51.92 , 'k*');
hold on;
plot(2012 , 48.08 , 'b*');
hold on;
plot(typeDistYearVec , mean((cin1_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cin1_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cin1_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(typeDistYearVec , mean((cin1_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((cin1_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((cin1_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
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
hold on;
plot(typeDistYearVec , mean((cin2_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cin2_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cin2_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(typeDistYearVec , mean((cin2_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((cin2_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((cin2_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('CIN2');
grid on;

% CIN3
% Calibration error bars
meanObs = cin3_dist_dObs(: , 2) .* 100;
sdevObs = (cin3_dist_dObs(: , 3).^(1/2)).*2 .* 100;

subplot(2,3,4)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 73.71 , 'k*');
hold on;
plot(2012 , 26.29 , 'b*');
hold on;
plot(typeDistYearVec , mean((cin3_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cin3_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cin3_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(typeDistYearVec , mean((cin3_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((cin3_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((cin3_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)');
title('CIN3');
grid on;

% CC
% Calibration error bars
meanObs = cc_dist_dObs(: , 2) .* 100;
sdevObs = (cc_dist_dObs(: , 3).^(1/2)).*2 .* 100;

subplot(2,3,5)
errorbar([2012, 2012] , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 1.5);
hold on;
plot(2012 , 85.78 , 'k*');
hold on;
plot(2012 , 14.22 , 'b*');
hold on;
plot(typeDistYearVec , mean((cc_vax .* 100),1)' , 'k-' , ...
    typeDistYearVec , min((cc_vax .* 100),[],1)' , 'k--' , ...
    typeDistYearVec , max((cc_vax .* 100),[],1)' , 'k--' , 'LineWidth' , 1.5);
hold on;
plot(typeDistYearVec , mean((cc_nonVax .* 100),1)' , 'b-' , ...
    typeDistYearVec , min((cc_nonVax .* 100),[],1)' , 'b--' , ...
    typeDistYearVec , max((cc_nonVax .* 100),[],1)' , 'b--' , 'LineWidth' , 1.5);
ylim([0 100]);
xlim([2010 2015]);
xlabel('Year'); ylabel('Prevalence Proportion by Type (%)')
title('Cervical Cancer')
legend('Observed 2012: mean, 2SD' , 'Observed 2012- 9v' , 'Observed 2012- non-9v' , ...
    'Model- 9v: ' , num2str(nRuns) , '-sets mean' , 'Model- 9v: ' , num2str(nRuns) , '-sets minimum' , 'Model- 9v: ' , num2str(nRuns) , '-sets maximum' , ...
    'Model- non-9v: ' , num2str(nRuns) , '-sets mean' , 'Model- non-9v: ' , num2str(nRuns) , '-sets minimum' , 'Model- non-9v: ' , num2str(nRuns) , '-sets maximum');
grid on;
