%function vaxCEA_multSims_IPVCsces()

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

%% LOAD SAVED RESULTS
% ***SET ME***: save names of potential scenarios to analyze as variables
baseDirName = 'Vaccine22Apr20Ph2V11_baseVax057_baseScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_8ts-2021_';
dirName_reductBaseline = [baseDirName , 'WHO-SCES0b_6_1'];
dirName_SCE2_057 = [baseDirName , 'WHO-SCES2_057_6_1'];
dirName_SCE2_057gbV = [baseDirName , 'WHO-SCES2_057gbV_6_1'];
dirName_SCE7_057 = [baseDirName , 'WHO-SCES7_057_6_1']; % Note that these really should be named WHO-SCES8...

% ***SET ME***: choose which scenarios you want to save data in Excel for
simVec = {dirName_reductBaseline , dirName_SCE2_057 , ...
    dirName_SCE2_057gbV , dirName_SCE7_057};
fileVec = {'sim0' , 'sim2' , 'sim2' , 'sim2'};
fileVec2 = {'0b' , '2_057' , '2_057gbV' ,'7_057'};
% ***SET ME***: make sure the names here correspond to scenarios in simVec above
fileTits = {'Baseline with 57% girls vax' , ...
    '57% girls multi-cohort vax' , ...
    '57% girls + boys multi-cohort vax' , ...
    '57% girls multi-cohort + catch-up vax'};
fileTits2 = {'Girls-only multi-cohort, Median' , '  Uncertainty range' , ...
    'Gender-neutral multi-cohort, Median' , '  Uncertainty range' , ...
    'Girls-only multi-cohort + catch-up, Median' , '  Uncertainty range'};

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2.0)

nResults = length(simVec);

for j = 1 : nResults 
    %% CC INCIDENCE - age standardized
    % Note: the age-standardization process shifts the incidence rate of the
    % last modelled age group to the next age group in the following year.
    % However, CC incidence is NaN prior to HIV introduction in the
    % HIV-positive no ART group, and NaN prior to ART introduction in the
    % HIV-positive ART group. Since we have four age groups past the 16 we
    % model, a NaN value is present for four years past the introduction of
    % HIV/ART, leading to a NaN value for summed incidence during these 
    % years. We therefore lack data in this four-year interval in the
    % saved/plotted results.
    
    fac = 10 ^ 5;
    worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
        247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
        23477 9261 2155];
    
    diseaseLabels = {'Pop(All) ICC' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'};    
    diseaseTits = {'General' , 'HIV-negative' , 'All WLHIV' , 'WLHIV-untreated' , 'WLHIV-on ART'};    
    
    for dInd = 1 : 1 %length(diseaseLabels)
        % Load results
        fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
            'UofW_Impact_CC_IncidenceRates-standardised-(2020-2120)_S' , fileVec2{j} , '.xlsx'];
        ccIncHivAgeTime = readmatrix(fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , ['C8:CZ21']);
        ccIncHivAgeTime = [zeros(2 , size(ccIncHivAgeTime , 2)) ; ccIncHivAgeTime];

        ccIncRefTot = zeros(1 , size(ccIncHivAgeTime,2));       
        for aInd = 1:age+4
            a = aInd;
            if aInd >= age
                a = age;
            end
            if aInd <= age    
                ccIncRef = ccIncHivAgeTime(a , :) .* worldStandard_WP2015(aInd);
                if (dInd == 5) && (a < 3)
                    ccIncRef = zeros(1 , size(ccIncHivAgeTime,2));
                end
            elseif aInd > age
                ccIncRef = ccIncHivAgeTime(a , :);
                ccIncRef = [(ones(1,aInd-a).*ccIncRef(1,1)) , ccIncRef(1,1:end-(aInd-a))];
                ccIncRef = ccIncRef .* worldStandard_WP2015(aInd);
            end
            ccIncRefTot = ccIncRefTot + ccIncRef;
        end
        ccInc = ccIncRefTot ./ (sum(worldStandard_WP2015(1:age+4)));
             
        % Plot baseline incidence
        if (j == 1) && (dInd == 1)
            fig = figure;
            set(fig,'DefaultAxesFontSize' , 18);
        end
        %subplot(1,3,dInd);
        hold all;
        firstYrInd = (currYear-1-startYear)+1;
        firstYrInd2 = (currYear-1-startYear)+1;
        p = plot([2020:2121] , ccInc(1 , :) , '-');
        axis([2020 2121 0 80])
        grid on;
        xlabel('Year'); ylabel('AS ICC per 100K');
        title(diseaseTits{dInd});
        
        if j == nResults
            hold all;
            plot ([2020:2121] , ones(1,size(ccIncHivAgeTime,2)).*4.0 , 'k:')
            hold all;
            plot ([2020:2121] , ones(1,size(ccIncHivAgeTime,2)).*10.0 , 'k--')
        end  
    end
    legend(fileTits{:} , 'Elimination: <4/100K' , 'Benchmark: <10/100K');
end

for j = 1 : nResults 
    %% Age-standardized HPV prevalence over time
    % Note: the age-standardization process shifts the prevalence of the
    % last modelled age group to the next age group in the following year.
    % However, prevalence is NaN prior to HIV introduction in the
    % HIV-positive no ART group, and NaN prior to ART introduction in the
    % HIV-positive ART group. Since we have four age groups past the 16 we
    % model, a NaN value is present for four years past the introduction of
    % HIV/ART, leading to a NaN value for summed HPV infected during these 
    % years. We therefore lack data in this four-year interval in the
    % saved/plotted results.
    worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
        247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
        23477 9261 2155];
    
    %diseaseLabels = {'Pop(All) (hrHPV)' , 'HIV-  (hrHPV)' , 'HIV+   (hrHPV)' , 'HIV+ no ART  (hrHPV)' , 'HIV+ ART  (hrHPV)'};
    diseaseLabels = {'Pop(All) (vtHPV)' , 'HIV-  (vtHPV)' , 'HIV+   (vtHPV)' , 'HIV+ no ART  (vtHPV)' , 'HIV+ ART  (vtHPV)'};

    for dInd = 1 : 1 %length(diseaseLabels)
        % Load results
        fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
            'UofW_Impact_VThrHPV-MWcombined_Prevalence-standardised-(2020-2120)_S' , fileVec2{j} , '.xlsx'];
        hpv_hivAgeW_dis = readmatrix(fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , ['B4:CY4'])' ./ sum(worldStandard_WP2015);
        
        % Plot hrHPV prevalence
        if (j == 1) && (dInd == 1)
            fig = figure;
            set(fig,'DefaultAxesFontSize' , 18);
        end
        %subplot(1,3,dInd);
        hold all;
        p = plot([2020 : 2121] , hpv_hivAgeW_dis , '-');
        axis([2020 2121 0 0.07])
        grid on;
        xlabel('Year'); ylabel('9v-type hrHPV population prevalence');
        title(diseaseTits{dInd});
    end
    legend(fileTits{:});
end

%% Calculate number vaccine doses per cervical cancer case averted
fname = [pwd , '\HHCoM_Results\' , simVec{1} , '\' , ...
    'CumulativeImpact_CC-standardised_wUncert-(2020-2120)_S' , fileVec2{1} , '.xlsx'];
ccCumCasesTot_base = readmatrix(fname , 'Sheet' , 'Pop(All) CCC' , 'Range' , ['CW5:CW29']);

for j = 1 : nResults 
    fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
        'CumulativeImpact_CC-standardised_wUncert-(2020-2120)_S' , fileVec2{j} , '.xlsx'];
    ccCumCasesTot = readmatrix(fname , 'Sheet' , 'Pop(All) CCC' , 'Range' , ['CW5:CW29']);
    ccCumCasesAverted = ccCumCasesTot_base - ccCumCasesTot;

    fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
        'RoutineSchoolVax_wUncert-(2021-2121)_S' , fileVec2{j} , '.xlsx'];
    routineVaxTot = readmatrix(fname , 'Sheet' , 'Vaccinated (HIV+ on ART) (N)' , 'Range' , ['A4:A28']); % Note: sheet name doesn't make sense; accidently autopopulated variable from previous section
    fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
        'CUVax_wUncert-(2021-2121)_S' , fileVec2{j} , '.xlsx'];
    cuVaxTot = readmatrix(fname , 'Sheet' , 'Vaccinated (HIV+ on ART) (N)' , 'Range' , ['A4:A28']); % Note: sheet name doesn't make sense; accidently autopopulated variable from previous section
    vaxTot = routineVaxTot + cuVaxTot;
    personsVaxTot = (routineVaxTot./2) + (cuVaxTot./3);
    
    vaxPerCase = vaxTot ./ ccCumCasesAverted;
    vaxPersonsPerCase = personsVaxTot ./ ccCumCasesAverted;
    
    disp(['Scenario vaccine efficiency ' , num2str(j) , ': '])
    disp(['Median: ' , num2str(median(vaxPerCase , 1))])
    disp([' Minimum: ' , num2str(min(vaxPerCase , [] , 1))])
    disp([' Maximum: ' , num2str(max(vaxPerCase , [] , 1))])
    
    disp(['Scenario vaccine efficiency per person ' , num2str(j) , ': '])
    disp(['Median: ' , num2str(median(vaxPersonsPerCase , 1))])
    disp([' Minimum: ' , num2str(min(vaxPersonsPerCase , [] , 1))])
    disp([' Maximum: ' , num2str(max(vaxPersonsPerCase , [] , 1))])
    
    disp(['Scenario vaccine doses ' , num2str(j) , ': '])
    disp(['Median: ' , num2str(median(vaxTot , 1))])
    disp([' Minimum: ' , num2str(min(vaxTot , [] , 1))])
    disp([' Maximum: ' , num2str(max(vaxTot , [] , 1))])
    
    disp(['Scenario # vaccinated ' , num2str(j) , ': '])
    disp(['Median: ' , num2str(median(personsVaxTot , 1))])
    disp([' Minimum: ' , num2str(min(personsVaxTot , [] , 1))])
    disp([' Maximum: ' , num2str(max(personsVaxTot , [] , 1))])
    
    disp(['Scenario CC averted ' , num2str(j) , ': '])
    disp(['Median: ' , num2str(median(ccCumCasesAverted , 1))])
    disp([' Minimum: ' , num2str(min(ccCumCasesAverted , [] , 1))])
    disp([' Maximum: ' , num2str(max(ccCumCasesAverted , [] , 1))])
end

%% Vaccine efficiency over time
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
colorVec = {[0.5 , 0.5, 0.5] , 'b' , 'r'};

fname = [pwd , '\HHCoM_Results\' , simVec{1} , '\' , ...
    'CumulativeImpact_CC-standardised_wUncert-(2020-2120)_S' , fileVec2{1} , '.xlsx'];
ccCumCasesTot_base = readmatrix(fname , 'Sheet' , 'Pop(All) CCC' , 'Range' , ['A5:CW29']);

for j = 2 : nResults 
    fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
        'CumulativeImpact_CC-standardised_wUncert-(2020-2120)_S' , fileVec2{j} , '.xlsx'];
    ccCumCasesTot = readmatrix(fname , 'Sheet' , 'Pop(All) CCC' , 'Range' , ['A5:CW29']);
    ccCumCasesAverted = ccCumCasesTot_base - ccCumCasesTot;

    fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
        'RoutineSchoolVax_wUncert-(byYear)_S' , fileVec2{j} , '.xlsx'];
    routineVaxTot = readmatrix(fname , 'Sheet' , 'HIV+ on ART (CCC)' , 'Range' , ['A4:CW28']); % Note: sheet name doesn't make sense; accidently autopopulated variable from previous section
    fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
        'CUVax_wUncert-(byYear)_S' , fileVec2{j} , '.xlsx'];
    cuVaxTot = readmatrix(fname , 'Sheet' , 'HIV+ on ART (CCC)' , 'Range' , ['A4:CW28']); % Note: sheet name doesn't make sense; accidently autopopulated variable from previous section
    personsVaxTot = (routineVaxTot./2) + (cuVaxTot./3);
    
    vaxPersonsPerCase = cumsum(personsVaxTot,2) ./ ccCumCasesAverted;
    
    hold all;
    p = plot([2021:2121] , median(vaxPersonsPerCase , 1) , 'Color' , colorVec{j-1})
    hold all;
    x2 = [[2021:2121] , fliplr([2021:2121])];
    inBetween = [max(vaxPersonsPerCase , [] , 1) , fliplr(min(vaxPersonsPerCase , [] , 1))];
    colorP = get(p,'Color');
    h = fill(x2 , inBetween , colorP);
    h.FaceAlpha = 0.2;
    h.LineStyle = '--';
    set(h,'EdgeColor', colorP);
    set(gca, 'YScale', 'log'); ax = gca; ax.XGrid = 'on';
    xlim([2020 2121]); ylim([0 10^12]);
    ylabel('Vaccine efficiency (log scale)'); xlabel('Year');
end
legend(fileTits2{1:end})
