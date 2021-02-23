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
dirName_S9 = [baseDirName , '9_6_1']; 
dirName_S10 = [baseDirName , '10_6_1'];

% ***SET ME***: choose which scenarios you want to save data in Excel for
simVec = {dirName_S0 , dirName_S1 , dirName_S2 , dirName_S3 , dirName_S4 , ...
    dirName_S5 , dirName_S6 , dirName_S7 , dirName_S8 , dirName_S9 , dirName_S10};
fileVec = {'0' , '1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10'}; 
colorVec = {[0.5 , 0.5 , 0.5] , [0.5 , 0.5 , 0.5] , [0.5 , 0.5 , 0.5] , [0 0.4470 0.7410] , ...
    [0 0.4470 0.7410] , [0.4660 , 0.6740 , 0.1880] , [0.4660 , 0.6740 , 0.1880] , ...
    [0.8500 0.3250 0.0980] , [0.8500 0.3250 0.0980] , [0.4940 , 0.1840 , 0.5560] , [0.4940 , 0.1840 , 0.5560]}; 
styleVec = {'-.' , ':' , '-' , '-' , '--' , '-' , '--' , '-' , '--' , '-' , '--'}; 
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
fileTits = {'Cytology', ...
    'HPV DNA' , ... 
    'HPV DNA + genotyping' , ...
    'AVE' , 'HPV DNA + AVE triage'};
fileTits2 = {'Three-visit cytology + colposcopy + treat' , 'One-visit HPV DNA test + treat' , 'One-visit HPV DNA test w/ genotyping + treat' , ...
    'One-visit AVE + treat' , 'One-visit HPV DNA test w/ AVE triage + treat' , '57% 9v vaccine; 48% screening at age 35' , ...
    '90% 9v vaccine; 48% screening at age 35' , ...
    '90% 9v vaccine; 48% repeat screening by HIV status' , ...
    '90% 9v vaccine; scaled-up repeat screening by HIV status'};
fileTits3 = {'CIN1' , 'HPV+ (no CIN)' , 'HPV-susceptible' , 'HIV-negative' , 'HIV-positive, untreated' , 'HIV-positive on ART'};
nResults = length(simVec);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)

%% Cumulative cervical cancer cases over time
diseaseLabels = {'Tot (CCC)' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)'}; 
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);

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
f(5) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{10});
f(6) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{1});
f(7) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{2});
f(8) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{3}); 
f(9) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{5});

legend(f , fileTits2{:} , 'Location' , 'northwest');
xlim([2020 2120]);
grid on;
xlabel('Year'); ylabel('Cumulative cervical cancer cases'); 
title('Cumulative cervical cancer cases');

%% Cervical cancer cases averted over time
diseaseLabels = {'Tot (CCC)' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)'};   
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
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
f(5) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{10});
f(6) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{1});
f(7) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{2});
f(8) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{3}); 
f(9) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{5});

legend(f , fileTits2{:} , 'Location' , 'northwest');
xlim([2020 2120])
grid on;
xlabel('Year'); ylabel('Cervical cancer cases averted'); 
title('Cervical cancer cases averted');

%% Percent cervical cancer cases averted over time
diseaseLabels = {'Tot (CCC)' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)'};   
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
% Load baseline results
fname = [pwd , '\HHCoM_Results\' , simVec{1} , '\' , ...
    'SA_screening_S0.xlsx'];
ccCumHivTime_baseline = readmatrix(fname , 'Sheet' , diseaseLabels{1} , 'Range' , 'E3:AC103');
for j = 1 : nResults 
        % Load results
        fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
            'SA_screening_S' , fileVec{j} , '.xlsx'];
        ccCumHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{1} , 'Range' , 'E3:AC103');

        reductCases = (ccCumHivTime_baseline - ccCumHivTime) ./ ccCumHivTime_baseline .* 100;
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
f(5) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{10});
f(6) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{1});
f(7) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{2});
f(8) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{3}); 
f(9) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{5});

legend(f , fileTits2{:} , 'Location' , 'northwest');
xlim([2020 2120]); ylim([0 100]);
grid on;
xlabel('Year'); ylabel('Percent reduction'); 
title('Percent reduction in cervical cancer cases');

%% Percent cervical cancer cases averted over time - cytology; 90% 9v vaccine; 48% repeat screening by HIV status baseline
% diseaseLabels = {'Tot (CCC)' , 'HIV- (CCC)' , 'HIV+ (CCC)' , 'HIV+ no ART (CCC)' , 'HIV+ ART (CCC)'};   
% fig = figure;
% set(fig,'DefaultAxesFontSize' , 18);
% % Load baseline results
% fname = [pwd , '\HHCoM_Results\' , simVec{3} , '\' , ...
%     'SA_screening_S2.xlsx'];
% ccCumHivTime_baseline = readmatrix(fname , 'Sheet' , diseaseLabels{1} , 'Range' , 'E3:AC103');
% for j = 3 : nResults 
%         % Load results
%         fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
%             'SA_screening_S' , fileVec{j} , '.xlsx'];
%         ccCumHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{1} , 'Range' , 'E3:AC103');
% 
%         reductCases = (ccCumHivTime_baseline - ccCumHivTime) ./ ccCumHivTime_baseline .* 100;
%         medReduct = median(reductCases , 2);
%         minReduct = min(reductCases , [] , 2);
%         maxReduct = max(reductCases , [] , 2);
%         
%         hold all;
%         p = plot([2020:2120]' , medReduct , 'Color' , colorVec{j} , 'LineStyle' , styleVec{j});
% %         hold all;
% %         x2 = [[2020:2120] , fliplr([2020:2120])];
% %         inBetween = [maxReduct' , fliplr(minReduct')];
% %         colorP = get(p,'Color');
% %         h = fill(x2 , inBetween , colorP);
% %         h.FaceAlpha = 0.3;
% %         h.LineStyle = '--';
% end
% f(1) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{4});
% f(2) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{6});
% f(3) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{8});
% f(4) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{4}); 
% f(5) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{5});
% 
% legend(f , fileTits2{[2,3,4,7,8]} , 'Location' , 'northwest');
% xlim([2020 2120]); ylim([0 100]);
% grid on;
% xlabel('Year'); ylabel('Percent reduction'); 
% title('Percent reduction in cervical cancer cases');

%% Age-standardized cervical cancer incidence over time
diseaseLabels = {'Tot (ICC)' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'}; 
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
for j = 1 : nResults  
    for dInd = 1 : 1 %length(diseaseLabels)
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
f(5) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{10});
f(6) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{1});
f(7) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{2});
f(8) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{3}); 
f(9) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{5});
f(10) = plot(NaN,NaN,'k--');
f(11) = plot(NaN,NaN,'k:');

hold all;
plot (ccIncHivTime(: , 1) , ones(size(ccIncHivTime,1),1).*4.0 , 'k--')
hold all;
plot (ccIncHivTime(: , 1) , ones(size(ccIncHivTime,1),1).*10.0 , 'k:')
legend(f , fileTits2{:} , 'Elimination: <4/100K' , 'Benchmark: <10/100K' , 'Location' , 'northeast');
xlim([2020 2120])
grid on;
xlabel('Year'); ylabel('AS-ICC (per 100K)'); 
title('Age-standardized cervical cancer incidence');

%% Percent reduction in age-standardized cervical cancer incidence over time
diseaseLabels = {'Tot (ICC)' , 'HIV- (ICC)' , 'HIV+ (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'}; %{'HIV- (ICC)' , 'HIV+ no ART (ICC)' , 'HIV+ ART (ICC)'};
widthVec = {2}; %{7 3 1};
incRedSceVec = {}; %{4 , 6 , 8};
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
for jInd = 1 : nResults  %length(incRedSceVec) 
    j = jInd; %incRedSceVec{jInd};
    for dInd =  1 : 1 %length(diseaseLabels)
        % Load baseline results
        fname = [pwd , '\HHCoM_Results\' , simVec{1} , '\' , ...
            'SA_screening_S0.xlsx'];
        ccIncHivTime_baseline = readmatrix(fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'E3:AC103');
        % Load results
        fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
            'SA_screening_S' , fileVec{j} , '.xlsx'];
        ccIncHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{dInd} , 'Range' , 'E3:AC103');

        perReductInc = ((ccIncHivTime_baseline - ccIncHivTime) ./ ccIncHivTime_baseline) .* 100;
        medReduct = median(perReductInc , 2);
        minReduct = min(perReductInc , [] , 2);
        maxReduct = max(perReductInc , [] , 2);

        hold all;
        p = plot([2020:2120]' , medReduct , 'Color' , colorVec{j} , 'LineStyle' , styleVec{j} , 'LineWidth' , widthVec{dInd});
    %     hold all;
    %     x2 = [[2020:2120] , fliplr([2020:2120])];
    %     inBetween = [maxReduct' , fliplr(minReduct')];
    %     colorP = get(p,'Color');
    %     h = fill(x2 , inBetween , colorP);
    %     h.FaceAlpha = 0.3;
    %     h.LineStyle = '--';
    end
end
f(1) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{1}); 
f(2) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{4});
f(3) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{6});
f(4) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{8});
f(5) = plot(NaN,NaN,'o','MarkerEdgeColor' , colorVec{10});
f(6) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{1});
f(7) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{2});
f(8) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{3}); 
f(9) = plot(NaN,NaN,'Color' , colorVec{1},'LineStyle' , styleVec{5});
f(10) = plot(NaN,NaN,'k--');

hold all;
plot ([2020:2120]' , ones(size(ccIncHivTime,1),1).*85.0 , 'k--')
legend(f , fileTits2{:} , 'Elimination: >85% reduction' , 'Location' , 'southeast');
xlim([2020 2120]); ylim([0 100]);
grid on;
xlabel('Year'); ylabel('Percent reduction in AS-ICC'); 
title('Percent reduction in age-standardized cervical cancer incidence');

%% Number of women (over)screened in states <CIN2 over time
diseaseLabels = {'Tot (OS)' , 'HIV- (OS)' , 'HIV+ (OS)' , 'HIV+ no ART (OS)' , 'HIV+ ART (OS)'}; 
diseaseSheetInds = {2 , 4 , 5};
sceInds = [3 , 4 , 6 , 8 , 10];
fig = figure;
set(fig,'DefaultAxesFontSize' , 18)

ovrScrnVec = zeros(length(sceInds) , 3 , length(diseaseSheetInds));
for jInd = 1 : length(sceInds)
    j = sceInds(jInd);
    for dInd = 1 : length(diseaseSheetInds)
        d = diseaseSheetInds{dInd};
        % Load results
        fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
            'SA_screening_S' , fileVec{j} , '.xlsx'];
        scrnSusHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C3:C103');
        scrnHpvHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C104:C204') + readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C205:C305');
        scrnCin1HivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C306:C406') + readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C407:C507');
        ovrScrnVec(jInd , 3 , dInd) = scrnSusHivTime(end,1);
        ovrScrnVec(jInd , 2 , dInd) = scrnHpvHivTime(end,1);
        ovrScrnVec(jInd , 1 , dInd) = scrnCin1HivTime(end,1);
    end
end
f(1) = bar(NaN,NaN,'FaceColor' , 'r'); hold all;
f(2) = bar(NaN,NaN,'FaceColor' , 'b'); hold all;
f(3) = bar(NaN,NaN,'FaceColor' , [0 0.5 0]); hold all;
f(4) = bar(NaN,NaN,'FaceColor' , 'k' , 'FaceAlpha' , 1.0); hold all; 
f(5) = bar(NaN,NaN,'FaceColor' , 'k' , 'FaceAlpha' , 0.6); hold all;
f(6) = bar(NaN,NaN,'FaceColor' , 'k' , 'FaceAlpha' , 0.2);

b = bar([1:5],reshape(permute(ovrScrnVec,[1,3,2]),[5,9]) , 'stacked');
set(b(1:3),'FaceColor','r');
set(b(4:6),'FaceColor','b');
set(b(7:9),'FaceColor',[0 0.5 0]);
set(b([1,4,7]),'FaceAlpha',1.0);
set(b([2,5,8]),'FaceAlpha',0.6);
set(b([3,6,9]),'FaceAlpha',0.2);
grid on;
xlabel('Scenario'); ylabel('Number of women'); 
set(gca,'xticklabel',fileTits);
set(gca,'xtick', 1 : 1 : 5); xlim([0,6]);
title('Number of women screened in states <CIN2');
legend(f , fileTits3{:} , 'Location' , 'northwest');


%% Number of women false positive screen results
diseaseLabels = {'Tot (OS)' , 'HIV- (OS)' , 'HIV+ (OS)' , 'HIV+ no ART (OS)' , 'HIV+ ART (OS)'};
diseaseSheetInds = {2 , 4 , 5};
diseaseScrnInds = {1 , 3 , 8};
sceScrnPos = {spCyto.testSens , spCyto.testSens , spCyto.testSens , ...
    spHpvDna.testSens , spHpvDna.testSens , ...
    spGentyp.testSens , spGentyp.testSens , ...
    spAve.testSens , spAve.testSens , ...
    spHpvAve.testSens , spHpvAve.testSens};
sceInds = [3 , 4 , 6 , 8 , 10];
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);

falsePosVec = zeros(length(sceInds) , 3 , length(diseaseSheetInds));
for jInd = 1 : length(sceInds)
    j = sceInds(jInd);
    if (j == 6 || j == 7)
        for dInd = 1 : length(diseaseSheetInds)
            d = diseaseSheetInds{dInd};
            % Load results
            fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
                'SA_screening_S' , fileVec{j} , '.xlsx'];
            scrnSusHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C3:C103') .* sceScrnPos{j}(diseaseScrnInds{dInd},1);
            scrnHpvHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C104:C204') .* sceScrnPos{j}(diseaseScrnInds{dInd},2);
            scrnCin1HivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C306:C406') .* sceScrnPos{j}(diseaseScrnInds{dInd},2);
            falsePosVec(jInd , 3 , dInd) = scrnSusHivTime(end,1);
            falsePosVec(jInd , 2 , dInd) = scrnHpvHivTime(end,1);
            falsePosVec(jInd , 1 , dInd) = scrnCin1HivTime(end,1);
        end
    else
        for dInd = 1 : length(diseaseSheetInds)
            d = diseaseSheetInds{dInd};
            % Load results
            fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
                'SA_screening_S' , fileVec{j} , '.xlsx'];
            scrnSusHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C3:C103') .* sceScrnPos{j}(diseaseScrnInds{dInd},1);
            scrnHpvHivTime = (readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C104:C204') + readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C205:C305')) .* sceScrnPos{j}(diseaseScrnInds{dInd},2);
            scrnCin1HivTime = (readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C306:C406') + readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C407:C507')) .* sceScrnPos{j}(diseaseScrnInds{dInd},2);
            falsePosVec(jInd , 3 , dInd) = scrnSusHivTime(end,1);
            falsePosVec(jInd , 2 , dInd) = scrnHpvHivTime(end,1);
            falsePosVec(jInd , 1 , dInd) = scrnCin1HivTime(end,1);
        end
    end
end
f(1) = bar(NaN,NaN,'FaceColor' , 'r'); hold all;
f(2) = bar(NaN,NaN,'FaceColor' , 'b'); hold all;
f(3) = bar(NaN,NaN,'FaceColor' , [0 0.5 0]); hold all;
f(4) = bar(NaN,NaN,'FaceColor' , 'k' , 'FaceAlpha' , 1.0); hold all; 
f(5) = bar(NaN,NaN,'FaceColor' , 'k' , 'FaceAlpha' , 0.6); hold all;
f(6) = bar(NaN,NaN,'FaceColor' , 'k' , 'FaceAlpha' , 0.2);

b = bar([1:5],reshape(permute(falsePosVec,[1,3,2]),[5,9]) , 'stacked');
set(b(1:3),'FaceColor','r');
set(b(4:6),'FaceColor','b');
set(b(7:9),'FaceColor',[0 0.5 0]);
set(b([1,4,7]),'FaceAlpha',1.0);
set(b([2,5,8]),'FaceAlpha',0.6);
set(b([3,6,9]),'FaceAlpha',0.2);
grid on;
xlabel('Scenario'); ylabel('Number of women'); 
set(gca,'xticklabel',fileTits);
set(gca,'xtick', 1 : 1 : 5); xlim([0,6]);
title('Number of women false-positive screen');
legend(f , fileTits3{:} , 'Location' , 'northwest');

%% Number of women overtreated
diseaseLabels = {'Tot (OS)' , 'HIV- (OS)' , 'HIV+ (OS)' , 'HIV+ no ART (OS)' , 'HIV+ ART (OS)'};
diseaseSheetInds = {2 , 4 , 5};
diseaseScrnInds = {1 , 3 , 8};
sceScrnPos = {spCyto.testSens , spCyto.testSens , spCyto.testSens , ...
    spHpvDna.testSens , spHpvDna.testSens , ...
    spGentyp.testSens , spGentyp.testSens , ...
    spAve.testSens , spAve.testSens , ...
    spHpvAve.testSens , spHpvAve.testSens};
sceTreat = {0 , 0 , 0 , 0.95 , 0.95 , 0.95 , 0.95 , 0.95 , 0.95 , 0.95 , 0.95};
sceInds = [3 , 4 , 6 , 8 , 10];
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);

overTreatVec = zeros(length(sceInds) , length(diseaseSheetInds));
for jInd = 1 : length(sceInds)
    j = sceInds(jInd);
    if (j == 6 || j == 7)
        for dInd = 1 : length(diseaseSheetInds)
            d = diseaseSheetInds{dInd};
            % Load results
            fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
                'SA_screening_S' , fileVec{j} , '.xlsx'];
            scrnSusHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C3:C103') .* sceScrnPos{j}(diseaseScrnInds{dInd},1) .* sceTreat{j};
            scrnHpvHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C104:C204') .* sceScrnPos{j}(diseaseScrnInds{dInd},2) .* sceTreat{j};
            scrnCin1HivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C306:C406') .* sceScrnPos{j}(diseaseScrnInds{dInd},2) .* sceTreat{j};
            overTreatVec(jInd , 3 , dInd) = scrnSusHivTime(end,1);
            overTreatVec(jInd , 2 , dInd) = scrnHpvHivTime(end,1);
            overTreatVec(jInd , 1 , dInd) = scrnCin1HivTime(end,1);
        end
    else
        for dInd = 1 : length(diseaseSheetInds)
            d = diseaseSheetInds{dInd};
            % Load results
            fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
                'SA_screening_S' , fileVec{j} , '.xlsx'];
            scrnSusHivTime = readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C3:C103') .* sceScrnPos{j}(diseaseScrnInds{dInd},1) .* sceTreat{j};
            scrnHpvHivTime = (readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C104:C204') + ...
                readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C205:C305')) .* sceScrnPos{j}(diseaseScrnInds{dInd},2) .* sceTreat{j};
            scrnCin1HivTime = (readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C306:C406') + ...
                readmatrix(fname , 'Sheet' , diseaseLabels{d} , 'Range' , 'C407:C507')) .* sceScrnPos{j}(diseaseScrnInds{dInd},2) .* sceTreat{j};
            overTreatVec(jInd , 3 , dInd) = scrnSusHivTime(end,1);
            overTreatVec(jInd , 2 , dInd) = scrnHpvHivTime(end,1);
            overTreatVec(jInd , 1 , dInd) = scrnCin1HivTime(end,1);
        end
    end
end
f(1) = bar(NaN,NaN,'FaceColor' , 'r'); hold all;
f(2) = bar(NaN,NaN,'FaceColor' , 'b'); hold all;
f(3) = bar(NaN,NaN,'FaceColor' , [0 0.5 0]); hold all;
f(4) = bar(NaN,NaN,'FaceColor' , 'k' , 'FaceAlpha' , 1.0); hold all; 
f(5) = bar(NaN,NaN,'FaceColor' , 'k' , 'FaceAlpha' , 0.6); hold all;
f(6) = bar(NaN,NaN,'FaceColor' , 'k' , 'FaceAlpha' , 0.2);

b = bar([1:5],reshape(permute(overTreatVec,[1,3,2]),[5,9]) , 'stacked');
set(b(1:3),'FaceColor','r');
set(b(4:6),'FaceColor','b');
set(b(7:9),'FaceColor',[0 0.5 0]);
set(b([1,4,7]),'FaceAlpha',1.0);
set(b([2,5,8]),'FaceAlpha',0.6);
set(b([3,6,9]),'FaceAlpha',0.2);
grid on;
xlabel('Scenario'); ylabel('Number of women'); 
set(gca,'xticklabel',fileTits);
set(gca,'xtick', 1 : 1 : 5); xlim([0,6]);
title('Number of women overtreated');
legend(f , fileTits3{:} , 'Location' , 'northwest');
