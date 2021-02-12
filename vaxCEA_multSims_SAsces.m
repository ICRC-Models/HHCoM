%function vaxCEA_multSims_SAsces()

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
    baseline , who , spCyto , spHpvDna , spGentyp , spAve , ...
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
baseDirName = 'Vaccine22Apr20Ph2V11_2v57BaseVax_spCytoScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_SA-S';
dirName_S0 = [baseDirName , '0_6_1'];
dirName_S1 = [baseDirName , '1_6_1'];
dirName_S2 = [baseDirName , '2_6_1'];
dirName_S3 = [baseDirName , '3_6_1'];
dirName_S4 = [baseDirName , '4_6_1'];
dirName_S5 = [baseDirName , '5_6_1'];
dirName_S6 = [baseDirName , '6_6_1'];
dirName_S7 = [baseDirName , '7_6_1'];
dirName_S8 = [baseDirName , '8_6_1'];

% ***SET ME***: choose which scenarios you want to save data in Excel for
simVec = {dirName_S0 , dirName_S1 , dirName_S2 , dirName_S3 , dirName_S4 , ...
    dirName_S5 , dirName_S6 , dirName_S7 , dirName_S8};
fileVec = {'0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8'}; 
colorVec = {[0.5 , 0.5 , 0.5] , [0.5 , 0.5 , 0.5] , [0.5 , 0.5 , 0.5] , 'b' , 'b' , [0.0 , 0.5 , 0.0] , [0.0 , 0.5 , 0.0] , 'r' , 'r'}; 
styleVec = {'-.' , ':' , '-' , '-' , '--' , '-' , '--' , '-' , '--'}; 
% ***SET ME***: make sure the names here correspond to scenarios in simVec above
% fileTits = {'S0 (Baseline vaccination): median', 'range' , ...
%     'S1 (Vaccine scale-up): median' , 'range' , ...
%     'S2 (Repeat screening by HIV status): median' , 'range' , ...
%     'S3 (HPV-and-treat at baseline levels): median' , 'range' , ... 
%     'S4 (HPV-and-treat with scale-up): median' , 'range' , ...
%     'S5 (HPV genotyping at baseline levels): median' , 'range' , ...
%     'S6 (HPV genotyping at baseline levels): median' , 'range' , ...
%     'S7 (AVE-and-treat at current levels): median' , 'range' , ...
%     'S8 (AVE-and-treat with scale-up): median' , 'range'};
fileTits = {'S0 (Baseline vaccination): median', ...
    'S1 (Vaccine scale-up): median' , ...
    'S2 (Repeat screening by HIV status): median' , ...
    'S3 (HPV-and-treat at baseline levels): median' , ... 
    'S4 (HPV-and-treat with scale-up): median' , ...
    'S5 (HPV genotyping at baseline levels): median' , ...
    'S6 (HPV genotyping at baseline levels): median' , ...
    'S7 (AVE-and-treat at current levels): median' , ...
    'S8 (AVE-and-treat with scale-up): median'};
fileTits2 = {'Cytology' , 'HPV DNA test + treat' , 'HPV DNA test w/ genotyping + treat' , ...
    'AVE + treat' , '57% 9v vaccine; 48% screening at age 35' , ...
    '90% 9v vaccine; 48% screening at age 35' , ...
    '90% 9v vaccine; 48% repeat screening by HIV status' , ...
    '90% 9v vaccine; scaled-up repeat screening by HIV status'};
nResults = length(simVec);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)

%% Cumulative cervical cancer cases over time
diseaseLabels = {'Tot (CCC)' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)'}; 
figure;

for j = 1 : nResults  
    for dInd = 1 : 1 %length(diseaseLabels);
        % Load results
        fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
            'SA_screening_S' , fileVec{j} , '.xlsx'];
        ccCumHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'A3:D103');
        hold all;
        p = plot(ccCumHivTime(: , 1) , ccCumHivTime(: , 2) , 'Color' , colorVec{j} , 'LineStyle' , styleVec{j});
%         hold all;
%         x2 = [ccCumHivTime(: , 1)' , fliplr(ccCumHivTime(: , 1)')];
%         inBetween = [ccCumHivTime(: , 4)' , fliplr(ccCumHivTime(: , 3)')];
%         colorP = get(p,'Color');
%         h = fill(x2 , inBetween , colorP);
%         h.FaceAlpha = 0.3;
%         h.LineStyle = '--';
    end
end
f(1) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{1}); 
f(2) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{4});
f(3) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{6});
f(4) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{8});
f(5) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{1});
f(6) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{2});
f(7) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{3}); 
f(8) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{5});

legend(f , fileTits2{:} , 'Location' , 'northwest');
xlim([2020 2120]);
grid on;
xlabel('Year'); ylabel('Cumulative cervical cancer cases'); 
title('Cumulative cervical cancer cases');

%% Cervical cancer cases averted over time
diseaseLabels = {'Tot (CCC)' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)'};   
figure;
% Load baseline results
fname = [pwd , '\HHCoM_Results\' , simVec{1} , '\' , ...
    'SA_screening_S0.xlsx'];
ccCumHivTime_baseline = readmatrix(fname , 'Sheet' , diseaseLabels{1} , 'Range' , 'E3:AC103');
for j = 1 : nResults 
        % Load results
        fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
            'SA_screening_S' , fileVec{j} , '.xlsx'];
        ccCumHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{1} , 'Range' , 'E3:AC103');

        reductCases = ccCumHivTime_baseline - ccCumHivTime;
        medReduct = median(reductCases , 2);
        minReduct = min(reductCases , [] , 2);
        maxReduct = max(reductCases , [] , 2);
        
        hold all;
        p = plot([2020:2120]' , medReduct , 'Color' , colorVec{j} , 'LineStyle' , styleVec{j});
%         hold all;
%         x2 = [[2020:2120] , fliplr([2020:2120])];
%         inBetween = [maxReduct' , fliplr(minReduct')];
%         colorP = get(p,'Color');
%         h = fill(x2 , inBetween , colorP);
%         h.FaceAlpha = 0.3;
%         h.LineStyle = '--';
end
f(1) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{1}); 
f(2) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{4});
f(3) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{6});
f(4) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{8});
f(5) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{1});
f(6) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{2});
f(7) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{3}); 
f(8) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{5});

legend(f , fileTits2{:} , 'Location' , 'northwest');
xlim([2020 2120])
grid on;
xlabel('Year'); ylabel('Cervical cancer cases averted'); 
title('Cervical cancer cases averted');

%% Age-standardized cervical cancer incidence over time
diseaseLabels = {'Tot (ICC)' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'}; 
figure;
for j = 1 : nResults  
    for dInd = 1 : 1 %length(diseaseLabels);
        % Load results
        fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
            'SA_screening_S' , fileVec{j} , '.xlsx'];
        ccIncHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'A3:D103');
        hold all;
        p = plot(ccIncHivTime(: , 1) , ccIncHivTime(: , 2) , 'Color' , colorVec{j} , 'LineStyle' , styleVec{j});
%         hold all;
%         x2 = [ccIncHivTime(: , 1)' , fliplr(ccIncHivTime(: , 1)')];
%         inBetween = [ccIncHivTime(: , 4)' , fliplr(ccIncHivTime(: , 3)')];
%         colorP = get(p,'Color');
%         h = fill(x2 , inBetween , colorP);
%         h.FaceAlpha = 0.3;
%         h.LineStyle = '--';
    end
end
f(1) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{1}); 
f(2) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{4});
f(3) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{6});
f(4) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{8});
f(5) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{1});
f(6) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{2});
f(7) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{3}); 
f(8) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{5});

hold all;
plot (ccIncHivTime(: , 1) , ones(size(ccIncHivTime,1),1).*4.0 , 'k')
hold all;
plot (ccIncHivTime(: , 1) , ones(size(ccIncHivTime,1),1).*10.0 , 'k')
legend(f , fileTits2{:} , 'Elimination: <4/100K' , 'Benchmark: <10/100K' , 'Location' , 'northeast');
xlim([2020 2120])
grid on;
xlabel('Year'); ylabel('Cumulative cervical cancer cases'); 
title('Age-standardized cervical cancer incidence');

%% Percent reduction in age-standardized cervical cancer incidence over time
diseaseLabels = {'Tot (ICC)' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'};
figure;
% Load baseline results
fname = [pwd , '\HHCoM_Results\' , simVec{1} , '\' , ...
    'SA_screening_S0.xlsx'];
ccIncHivTime_baseline = readmatrix(fname , 'Sheet' , diseaseLabels{1} , 'Range' , 'E3:AC103');
for j = 1 : nResults 
    % Load results
    fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
        'SA_screening_S' , fileVec{j} , '.xlsx'];
    ccIncHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{1} , 'Range' , 'E3:AC103');

    perReductInc = ((ccIncHivTime_baseline - ccIncHivTime) ./ ccIncHivTime_baseline) .* 100;
    medReduct = median(perReductInc , 2);
    minReduct = min(perReductInc , [] , 2);
    maxReduct = max(perReductInc , [] , 2);

    hold all;
    p = plot([2020:2120]' , medReduct , 'Color' , colorVec{j} , 'LineStyle' , styleVec{j});
%     hold all;
%     x2 = [[2020:2120] , fliplr([2020:2120])];
%     inBetween = [maxReduct' , fliplr(minReduct')];
%     colorP = get(p,'Color');
%     h = fill(x2 , inBetween , colorP);
%     h.FaceAlpha = 0.3;
%     h.LineStyle = '--';
end
f(1) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{1}); 
f(2) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{4});
f(3) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{6});
f(4) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{8});
f(5) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{1});
f(6) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{2});
f(7) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{3}); 
f(8) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{5});

hold all;
plot ([2020:2120]' , ones(size(ccIncHivTime,1),1).*85.0 , 'k')
legend(f , fileTits2{:} , 'Elimination: >85% reduction' , 'Location' , 'southeast');
xlim([2020 2120]); ylim([0 100]);
grid on;
xlabel('Year'); ylabel('Percent reduction in AS-ICC'); 
title('Percent reduction in age-standardized cervical cancer incidence');
