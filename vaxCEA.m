function vaxCEA(pathModifier)

%% Load parameters and results
paramDir = [pwd , '\Params\'];

[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd , ...
    sumall , annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , mue , epsA_vec , epsR_vec , ...
    yr , ...
    hivOn , betaHIVM2F , betaHIVF2M , muHIV , kVl , kCD4 , ...
    hpvOn , perPartnerHpv_vax , perPartnerHpv_nonV , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , ...
    circ , condUse , screenYrs , hpvScreenStartYear , waning , ...
    maxRateM1 , maxRateF1 , maxRateM2 , maxRateF2 , hivStartYear , ...
    circStartYear , vaxStartYear , baseline , cisnet , who , whob , ...
    circProtect , condProtect , MTCTRate , hyst , ...
    OMEGA , ...
    ccInc2011_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , ...
    hivPrevM_dObs , hivPrevF_dObs , ...
    mCurr , fCurr , mCurrArt , fCurrArt , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , ...
    hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ...
    ageInd , riskInd , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , deathMat , circMat , circMat2] = loadUp2(1 , 0 , [] , [] , []);

% Helper functions
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

% Load results
nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
curr = load([pwd , '\HHCoM_Results\toNow_8Nov19_sameAssum_adjNonVtrans050_125vClr075nvClr_10-20-20-30_40-30-30_evenHIVInit_0_0_fixAcutART-incInitHIV_savIncHivArtChar_fixHIVMult']); % Population up to current year

vaxResult = cell(nSims , 1);
resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxSimResult'];
if waning
    resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxWaneSimResult'];
end
parfor n = 1 : nSims
    % load results from vaccine run into cell array
    vaxResult{n} = load([resultFileName , num2str(n), '.mat']);
    % concatenate vectors/matrices of population up to current year to population
    % matrices for years past current year
    vaxResult{n}.popVec = [curr.popVec(1 : end  , :); vaxResult{n}.popVec(2 : end , :)];
    vaxResult{n}.newHpvVax = [curr.newHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvVax(2 : end , : , : , : , : , :)];
    vaxResult{n}.newImmHpvVax = [curr.newImmHpvVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvVax(2 : end , : , : , : , : , :)];
    vaxResult{n}.newHpvNonVax = [curr.newHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newHpvNonVax(2 : end , : , : , : , : , :)];
    vaxResult{n}.newImmHpvNonVax = [curr.newImmHpvNonVax(1 : end , : , : , : , : , :); vaxResult{n}.newImmHpvNonVax(2 : end , : , : , : , : , :)];
    vaxResult{n}.ccDeath = [curr.ccDeath(1 : end , : , : , :) ; vaxResult{n}.ccDeath(2 : end , : , : , :)];
    vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :); vaxResult{n}.newCC(2 : end , : , : ,:)];
    vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , :); vaxResult{n}.newHiv(2 : end , : , : ,:)];
    vaxResult{n}.hivDeaths = [curr.hivDeaths(1 : end , : , :); vaxResult{n}.hivDeaths(2 : end , : , :)];
    vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , : , : , : , :); vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , : , : , : , :)];
    vaxResult{n}.tVec = [curr.tVec(1 : end), vaxResult{n}.tVec(2 : end)];
    %vaxResult{n}.ccTreated = [curr.ccTreated(1 : end) , vaxResult{n}.ccTreated(2 : end)];
end
noVaxInd = nSims;
noV = vaxResult{noVaxInd};
tVec = noV.tVec;
tVecYr = tVec(1 : stepsPerYear : end);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)

%% CC incidence and incidence reduction (age-standardized)
inds = {':' , 1 : 2 , 3 : 7 , 8 , 3 : 8};
fac = 10 ^ 5;
files = {'CC_General_Hpv_VaxCover' , ...
     'CC_HivNeg_Hpv_VaxCover' , 'CC_HivNoART_Hpv_VaxCover' ,...
     'CC_ART_HPV_VaxCover' , 'CC_HivAll_HPV_VaxCover'};
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
plotTits2 = {'80% coverage: Total female population' , '90% coverage'  , ...
    '80% coverage: HIV-Negative' , '90% coverage' , ...
    '80% coverage: HIV-Positive no ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive all' , '90% coverage'};
linStyle = {'-' , '--'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g'};
% worldStandard_WP2015 = [325428 (311262/5.0) 295693 287187 291738 299655 272348 ...
%         247167 240167 226750 201603 171975 150562 113118 82266 64484];
worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
        247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
        23477 9261 2155];

figure();

for i = 1 : length(inds)
    % % plotTits = {plotTits2{(i*2-1):(i*2)}};
    
    for n = 1 : nSims
        vaxResult{n}.ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
    end
     
    % % % General, all ages
    % % allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
    % %     1 : intervens , 2 , 3 : age , 1 : risk)); ...
    % %     toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
    % %     1 : intervens , 2 , 3 : age , 1 : risk))];
    % % allhivNegFAge = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : intervens , ...
    % %     2 , 3 : age , 1 : risk)); ...
    % %     toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : intervens , ...
    % %     2 , 3 : age , 1 : risk))];
    
    for aInd = 1 : age + 4
        if aInd >= age
            a = age;
        else
            a = aInd;
        end
        
        % General
        allF = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk));
        % All HIV-negative women
        hivNeg = toInd(allcomb(1 : 2 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk));
        % HIV-positive women not on ART
        hivNoARTF = toInd(allcomb(3 : 7 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk));
        % Women on ART
        artF = toInd(allcomb(8 , 6 , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk));
        % All HIV-positive women
        hivAllF = toInd(allcomb(3 : 8 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
            1 , 1 : intervens , 2 , a , 1 : risk));
        genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};
        
        % Calculate incidence
        for n = 1 : nSims
            if aInd <= age
                ccIncRef = ...
                    (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , a , :),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac) ...
                    .* (worldStandard_WP2015(a));
                    %.* (annlz(sum(noV.popVec(: , genArray{3}) , 2) ./ stepsPerYear));
                if (i == 4) && (a < 3) && (max(annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , a , :),2),3),4))) == 0.0)
                    ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
                end
            elseif aInd > age
                ccIncRef = ...
                    (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , a , :),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
                ccIncRef = [(ones(1,aInd-a).*ccIncRef(1,1)) , ccIncRef(1,1:end-(aInd-a))];
                ccIncRef = ccIncRef .* worldStandard_WP2015(aInd);
            end
            vaxResult{n}.ccIncRef = vaxResult{n}.ccIncRef + ccIncRef;
            % % vaxResult{n}.ccIncRef = ccIncRef;
        end
    end
    for n = 1 : nSims
        vaxResult{n}.ccInc = vaxResult{n}.ccIncRef ./ (sum(worldStandard_WP2015(1:age+4)));
    end
        
%     plot(tVec(1 : stepsPerYear : end) , vaxResult{noVaxInd}.ccInc ,'DisplayName' , ...
%          [plotTits1{i} , ': Efficacy: ' , num2str(round(noV.vaxEff * 100)) '% ,', ...
%          'Coverage: ' , num2str(round(noV.vaxRate * 100)) , '%']);
%     legend('-DynamicLegend');
%     grid on;
%     xlim([1980 2120]);
%     ylim([0 150]);
%     xticks([1980 : 20 : 2120]);
%     hold all;

    for n = 2 : length(vaxResult) - 1
%         plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccInc , 'DisplayName' , ...
%             [plotTits1{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
        
        % Reduction
        vaxResult{n}.ccRed = (vaxResult{n}.ccInc - vaxResult{noVaxInd}.ccInc) ./ vaxResult{noVaxInd}.ccInc * 100;
        plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccRed , ...
            'LineStyle' , linStyle{1} , 'Color' , linColor{i} , 'DisplayName' , plotTits1{i})
%             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,' , ...
%             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
        grid on
        legend('-DynamicLegend')
        xlim([2020 2120]);
        %ylim([-100 0]);
        xticks([2020 : 20 : 2120]);
        hold all
%               
%         % Save reduction results
% %         fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
% %             'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
% %         sname = [plotTits1{i} , '_IncRed'];
% %         if exist(fname , 'file') == 2
% %             M = xlsread(fname);
% %             M = catpad(2 , [tVec(1 : stepsPerYear : end)' , noV.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , M);
% %             xlswrite(fname , M , sname)
% %         else
% %             xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , sname)
% %         end
%         
     end
% %     title(' Cervical Cancer Incidence')
% %     xlabel('Year'); ylabel('Incidence per 100,000')
% %     hold all;
% 
%     title('Figure 1: Percent reduction in cervical cancer cases, by HIV status')
%     xlabel('Year'); ylabel('Percent change')
end        

%% CC incidence and incidence reduction (not age-standardized)
inds = {':' , 1 : 2 , 3 : 7 , 8 , 3 : 8};
fac = 10 ^ 5;
files = {'CC_General_Hpv_VaxCover' , ...
     'CC_HivNeg_Hpv_VaxCover' , 'CC_HivNoART_Hpv_VaxCover' ,...
     'CC_ART_HPV_VaxCover' , 'CC_HivAll_HPV_VaxCover'};
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
plotTits2 = {'80% coverage: Total female population' , '90% coverage'  , ...
    '80% coverage: HIV-Negative' , '90% coverage' , ...
    '80% coverage: HIV-Positive no ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive all' , '90% coverage'};
linStyle = {'-' , '--'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g'};
% worldStandard_WP2015 = [325428 (311262/5.0) 295693 287187 291738 299655 272348 ...
%         247167 240167 226750 201603 171975 150562 113118 82266 64484];
worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
        247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
        23477 9261 2155];

figure();

for i = 1 : length(inds)
    % % plotTits = {plotTits2{(i*2-1):(i*2)}};
    
    for n = 1 : nSims
        vaxResult{n}.ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
    end
        
    % General
    allF = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 3 : age , 1 : risk));
    % All HIV-negative women
    hivNeg = toInd(allcomb(1 : 2 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 3 : age , 1 : risk));
    % HIV-positive women not on ART
    hivNoARTF = toInd(allcomb(3 : 7 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 3 : age , 1 : risk));
    % Women on ART
    artF = toInd(allcomb(8 , 6 , [1 : 5 , 7] , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 3 : age , 1 : risk));
    % All HIV-positive women
    hivAllF = toInd(allcomb(3 : 8 , 1 : viral , [1 : 5 , 7] , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 3 : age , 1 : risk));
    genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};

    % Calculate incidence
    for n = 1 : nSims
        ccIncRef = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , 3 : age , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac);
        vaxResult{n}.ccIncRef = ccIncRef;
    end

    for n = 1 : nSims
        vaxResult{n}.ccInc = vaxResult{n}.ccIncRef;
    end
        
    plot(tVec(1 : stepsPerYear : end) , vaxResult{noVaxInd}.ccInc ,'DisplayName' , ...
         [plotTits1{i} , ': Efficacy: ' , num2str(round(noV.vaxEff * 100)) '% ,', ...
         'Coverage: ' , num2str(round(noV.vaxRate * 100)) , '%']);
    legend('-DynamicLegend');
    grid on;
    xlim([1980 2120]);
    ylim([0 250]);
    xticks([1980 : 20 : 2120]);
    hold all;

%     for n = 1 : length(vaxResult)-2
% %         plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccInc , 'DisplayName' , ...
% %             [plotTits1{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
% %             'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
%         
%         % Reduction
%         vaxResult{n}.ccRed = (vaxResult{n}.ccInc - noV.ccInc) ./ noV.ccInc * 100;
%         plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccRed , 'LineStyle' , linStyle{n} , 'Color' , linColor{i} , 'DisplayName' , plotTits1{i})
% %             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
% %             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
%         grid on
%         legend('-DynamicLegend')
%         xlim([2020 2120]);
%         %ylim([-100 0]);
%         xticks([2020 : 20 : 2120]);
%         hold all
%               
%         % Save reduction results
% %         fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
% %             'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
% %         sname = [plotTits1{i} , '_IncRed'];
% %         if exist(fname , 'file') == 2
% %             M = xlsread(fname);
% %             M = catpad(2 , [tVec(1 : stepsPerYear : end)' , noV.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , M);
% %             xlswrite(fname , M , sname)
% %         else
% %             xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , sname)
% %         end
%         
%     end
% %     title(' Cervical Cancer Incidence')
% %     xlabel('Year'); ylabel('Incidence per 100,000')
% %     hold all;
% 
%     title('Figure 1: Percent reduction in cervical cancer cases, by HIV status')
%     xlabel('Year'); ylabel('Percent change')
end        


%% CC incidence by risk
% inds = {':' , [2 : 6] , 1 , 10};
% files = {'CC_General_Hpv_VaxCover' , ...
%      'CC_HivNoART_Hpv_VaxCover' , 'CC_HivNeg_Hpv_VaxCover' ,...
%      'CC_ART_HPV_VaxCover'};
% plotTits = {'General' , 'HIV-Positive (No ART)' , ....
%      'HIV-Negative' , 'HIV-Positive on ART'};
% fac = 10 ^ 5;
% 
% figure();
% for i = 2 : length(inds)
% %     figure();
%     % General, all risk groups
%     allFRisk = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
%         1 : intervens , 2 , 8 : age , 1 : risk)); ...
%         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
%         1 : intervens , 2 , 8 : age , 1 : risk))];
%     % HIV-positive women not on ART
%     hivNoARTFRisk = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
%         1 : intervens , 2 , 8 : age , 1 : risk)); ...
%         toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
%         1 : intervens , 2 , 8 : age , 1 : risk))];
%     % All HIV-negative women
%     hivNegRisk = [toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : intervens , ...
%         2 , 8 : age , 1 : risk)); ...
%         toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : intervens , ...
%         2 , 8 : age , 1 : risk))];
%     % Women on ART
%     artFRisk = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : 4 , ...
%         1 : intervens , 2 , 8 : age , 1 : risk)); ...
%         toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 : 10 , ...
%         1 : intervens , 2 , 8 : age , 1 : risk))];
%     genArrayRisk = {allFRisk , hivNoARTFRisk , hivNegRisk , artFRisk};
%     
%     for r = 1 : risk
%         % General
%         allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
%             1 : intervens , 2 , 8 : age , r)); ...
%             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , 8 : age , r))];
%         % HIV-positive women not on ART
%         hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
%             1 : intervens , 2 , 8 : age , r)); ...
%             toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , 8 : age , r))];
%         % All HIV-negative women
%         hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : intervens , ...
%             2 , 8 : age , r)); ...
%             toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : intervens , ...
%             2 , 8 : age , r))];
%         % Women on ART
%         artF = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : 4 , ...
%             1 : intervens , 2 , 8 : age , r)); ...
%             toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , 8 : age , r))];
%         genArray = {allF , hivNoARTF , hivNeg , artF};
% 
%         prop = ((annlz(sum(noV.popVec(: , genArray{i}) , 2))) ./ (annlz(sum(noV.popVec(: , genArrayRisk{i}) , 2)))) .* 100;
%             
%         
%         plot(tVec(1 : stepsPerYear : end) , prop ,'DisplayName' , ...
%              ['Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%              'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
%         ylim([0 100]);
%         hold all;
%         legend('Low risk' , 'Medium risk' , 'High risk');
%         
%     end
%     title([plotTits{i} , ' Proportion by Risk']);
%     xlabel('Year'); ylabel('Percent');
% %     hold off;
% end        

%% CC mortality and mortality reduction
% inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]};
% files = {'CC_General_Hpv_VaxCover' , ...
%      'CC_HivNeg_Hpv_VaxCover' , 'CC_HivNoART_Hpv_VaxCover' ,...
%      'CC_ART_HPV_VaxCover' , 'CC_HivAll_HPV_VaxCover'};
% plotTits1 = {'General' , 'HIV-' , ...
%     'HIV+ no ART' , 'HIV+ ART' , 'HIV all'};
% plotTits2 = {'80% coverage: Total female population' , '90% coverage'  , ...
%     '80% coverage: HIV-Negative' , '90% coverage' , ...
%     '80% coverage: HIV-Positive no ART' , '90% coverage' , ...
%     '80% coverage: HIV-Positive ART' , '90% coverage' , ...
%     '80% coverage: HIV-Positive all' , '90% coverage'};
% fac = 10 ^ 5;
% linStyle = {'-' , '--'};
% linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g'};
% 
% figure();
% 
% for i = 1 : length(inds)-1
%     plotTits = {plotTits2{(i*2-1):(i*2)}};
% % %     figure();
% % %     noV.ccMortRef = zeros(length(tVec(length(curr.tVec) + 1 : stepsPerYear : end)),1)';
% % %     for n = 1 : length(vaxResult)-1
% % %         vaxResult{n}.ccMortRef = zeros(length(tVec(length(curr.tVec) + 1 : stepsPerYear : end)),1)';
% % %     end
%     
%     % General, all ages
% % %     allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
% % %         1 : intervens , 2 , 3 : age , 1 : risk)); ...
% % %         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
% % %         1 : intervens , 2 , 3 : age , 1 : risk))];
% % %     allhivNegFAge = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : intervens , ...
% % %             2 , 3 : age , 1 : risk)); ...
% % %             toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : intervens , ...
% % %             2 , 3 : age , 1 : risk))];
%       
% % %     for a = 3 : age
%         % General
%         allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
%             1 : intervens , 2 , 3 : age , 1 : risk)); ...
%             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , 3 : age , 1 : risk))];
%         % All HIV-negative women
%         hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : intervens , ...
%             2 , 3 : age , 1 : risk)); ...
%             toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : intervens , ...
%             2 , 3 : age , 1 : risk))];
%         % HIV-positive women not on ART
%         hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
%             1 : intervens , 2 , 3 : age , 1 : risk)); ...
%             toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , 3 : age , 1 : risk))];
%         % Women on ART
%         artF = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : 4 , ...
%             1 : intervens , 2 , 3 : age , 1 : risk)); ...
%             toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , 3 : age , 1 : risk))];
%         % All HIV-positive women
%         hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvVaxStates , 1 : 4 , ...
%             1 : intervens , 2 , 3 : age , 1 : risk)); ...
%             toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , 3 : age , 1 : risk))];
%         genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};
% 
%         ccMortRef = ...
%             (annlz(sum(sum(sum(noV.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
%             (annlz(sum(noV.popVec(length(curr.tVec) : end , genArray{i}) , 2) ./ stepsPerYear))* fac);
% % %             .* (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , genArray{1}) , 2) ./ stepsPerYear));
% % %         noV.ccMortRef = noV.ccMortRef + ccMortRef;
%         noV.ccMortRef = ccMortRef;
%                 
%         for n = 1 : length(vaxResult)-1
%             ccMortRef = ...
%                 (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
%                 (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) : end  , genArray{i}) , 2) ./ stepsPerYear)) * fac);
% % %                 .* (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end  , genArray{1}) , 2) ./ stepsPerYear));
% % %             vaxResult{n}.ccMortRef = vaxResult{n}.ccMortRef + ccMortRef;
%             vaxResult{n}.ccMortRef = ccMortRef;
%         end
%         
% % %     end
% % %     noV.ccMort = noV.ccMortRef ./ (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , allFAge) , 2) ./ stepsPerYear));
%       noV.ccMort = noV.ccMortRef;
% %     plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , noV.ccMort ,'DisplayName' , ...
% %          [plotTits1{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
% %          'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
% %     legend('-DynamicLegend');
% %     hold all;
%     for n = 1 : length(vaxResult)-1
% % %         vaxResult{n}.ccMort = vaxResult{n}.ccMortRef ./ (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end , allFAge) , 2) ./ stepsPerYear));
%         vaxResult{n}.ccMort = vaxResult{n}.ccMortRef;
% %         plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , vaxResult{n}.ccMort , 'DisplayName' , ...
% %             [plotTits1{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
% %             'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
%         
%         % Reduction
%         vaxResult{n}.ccRed = (vaxResult{n}.ccMort - noV.ccMort) ./ noV.ccMort * 100;
%         plot(tVec(length(curr.tVec) : stepsPerYear : end) , vaxResult{n}.ccRed , 'LineStyle' , linStyle{n} , 'Color' , linColor{i} , 'DisplayName' , plotTits{n}) 
% %             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
% %             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
%         grid on        
%         legend('-DynamicLegend')
%         xlim([2019 2099]);
%         xticks([2019 : 10 : 2099]);
%         hold all
%         
%         % Save reduction results
% %         fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
% %             'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_Mort' , '.xlsx'];
% %         sname = [plotTits1{i} , '_MortRed'];
% %         if exist(fname , 'file') == 2
% %             M = xlsread(fname);
% %             M = catpad(2 , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , noV.ccMort' , vaxResult{n}.ccMort' , vaxResult{n}.ccRed'] , M);
% %             xlswrite(fname , M , sname)
% %         else
% %             xlswrite(fname , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , noV.ccMort' , vaxResult{n}.ccMort' , vaxResult{n}.ccRed'] , sname)
% %         end
%         
%     end
% %     title('Cervical Cancer Mortality')
% %     xlabel('Year'); ylabel('Incidence per 100,000')
% %     hold all;
% 
%     title('Figure 2: Percent reduction in cervical cancer mortality, by HIV status')
%     xlabel('Year'); ylabel('Percent change')
% end

%% CC Mortality by age
% inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
% aVec = {11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80}; %{10:15,16:25,26:35,36:50,51:75};
% ageGroup = {'10-14' , '15-19' , '20-24' , '25-29' ,...
%      '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
%      '60-64' , '65-69' , '70-74' , '75-79'};
% aMatrix = zeros(length(inds) , length(aVec));
% for i = 1 : length(inds)
%     for aInd = 1 : length(aVec)
%         a = aVec{aInd};
% 
%         % General
%         allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : 7 , ...
%             1 : intervens , 2 , a , 1 : risk)); ...
%             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , a , 1 : risk))];
%         % All HIV-negative women
%         hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : 7 , 1 : intervens , ...
%             2 , a , 1 : risk)); ...
%             toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : intervens , ...
%             2 , a , 1 : risk))];
%         % HIV-positive women not on ART
%         hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : 7 , ...
%             1 : intervens , 2 , a , 1 : risk)); ...
%             toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , a , 1 : risk))];
%         % Women on ART
%         artF = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : 7 , ...
%             1 : intervens , 2 , a , 1 : risk)); ...
%             toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , a , 1 : risk))];
%         % All HIV-positive women
%         hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvVaxStates , 1 : 7 , ...
%             1 : intervens , 2 , a , 1 : risk)); ...
%             toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvVaxStates , 9 : 10 , ...
%             1 : intervens , 2 , a , 1 : risk))];
%         genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};
% 
%         % Calculate mortality
%         for n = nSims %1 : nSims
%             ccMortRef = ...
%                 (annlz(sum(sum(sum(vaxResult{n}.ccDeath(end-5:end , inds{i} , : , a),2),3),4)) ./ ...
%                 (annlz(sum(vaxResult{n}.popVec(end-5:end , genArray{i}) , 2) ./ stepsPerYear)) * fac);
%             
%             aMatrix(i,aInd) = ccMortRef;
%         end
%     end
%             
%     plot([1:length(aVec)] , aMatrix(i,:))
%     hold all;
%     ylabel('Mortality rates (per 100K)'); 
%     legend('General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all') % , 'General-ARToff' , 'HIV-negative-ARToff' , 'HIV-positive no ART-ARToff' , 'HIV-positive ART-ARToff' , 'HIV all-ARToff');
%     set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
% 
% end

%% Population Size
% % HIV-positive women not on ART
% hivNoART = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%     1 : intervens , 1 : gender , 10 : 75 , 1 : risk));
% % All HIV-negative women
% hivNeg = toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%     1 : gender , 10 : 75 , 1 : risk));
% % Women on ART
% art = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%     1 : intervens , 1 : gender , 10 : 75 , 1 : risk));
% genArray = {hivNeg , hivNoART , art};
% 
% %figure()
% subplot(2,1,1);
% plot(tVec , sum(noV.popVec(:,genArray{1}),2) + sum(noV.popVec(:,genArray{2}),2) + sum(noV.popVec(:,genArray{3}),2))
% % hold all
% % for i = 1 : length(genArray)
% %     plot(tVec , sum(noV.popVec(: , genArray{i}) , 2))
% %     hold all;
% % end
% title('Population Size')
% xlabel('Year'); ylabel('Individuals')
% %xlim([1950 2120]);
% legend('Total Population' , 'HIV-Negative' , 'HIV-Positive no ART' , 'HIV-Positive ART');
% hold off
% 
% %figure;
% subplot(2,1,2);
% aVec = {10:15,16:25,26:35,36:50,51:75};
% popProp = zeros(length(tVec),length(aVec));
% for aInd = 1:length(aVec)
%     a = aVec{aInd};
%     popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , 1 : gender , a , 1 : risk));
%     popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , 1 : gender , 10 : 75 , 1 : risk));
%     popProp(:,aInd) = sum(noV.popVec(:,popAge),2) ./ sum(noV.popVec(:,popTot),2);
% end
% plot(tVec , popProp);
% ylabel('Proportion'); xlabel('Year');
% %xlim([1950 2120]);
% legend('9-14' , '15-24' , '25-34' , '35-49' , '50-74');

%% Proportion HIV population on ART
% gender = {'Females (ages 15-74)','Males (ages 15-74)'};
% gList = [2,1];
% 
% for gInt = 1 : 2
%     g = gList(gInt);
%     artInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , g , 16 : 75 , 1 : risk));
%     artPop = sum(noV.popVec(: , artInds) , 2);
%     hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%         1 : intervens , g , 16 : 75 , 1 : risk));
%     hivPop = sum(noV.popVec(: , hivInds) , 2);
%     subplot(2,1,gInt);
%     plot(tVec , 100 * artPop ./ (hivPop + artPop))
%     hold on
%     xlabel('Year')
%     %ylabel('Proportion of HIV Population')
%     ylabel('Coverage (%)');
%     title(gender{gInt});
%     %legend('Model (Male)' , 'Model (Female)')
%     axis([2000 2050 0 100])
% 
% end
% % xlabel('Year')
% % ylabel('Proportion of HIV Population')
% % title('Proportion on ART')
% % legend('Model (Male)' , 'Model (Female)')

%% On ART by age
% ageGroup = {'0-4' , '5-9' ,'10-14' , '15-19' , '20-24' , '25-29' ,...
%      '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
%      '60-64' , '65-69' , '70-74' , '75-79'};
% aMatrix = zeros(1 , age);
% for a = 1 : age
%     artInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , 2 , a , 1 : risk));
%     artPop = sum(noV.popVec(end , artInds) , 2); %end-605
%     hivInds = toInd(allcomb([2 : 6,10] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%         1 : intervens , 2 , a , 1 : risk));
%     hivPop = sum(noV.popVec(end , hivInds) , 2);
%     hiv_art = [100 * artPop ./ hivPop];
% 
%     aMatrix(a) = hiv_art;
% end
% 
% figure;
% hold all;
% plot([1:age] , aMatrix(1,:) , ':')
% hold all;
% ylabel('Percent females on ART in 2120 by age');
% ylim([0 110])
% set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
% legend('Without ART dropout' , 'With ART dropout: 6.19%' , 'With ART dropout: 11.8%' , 'With ART dropout: 11.8%, HIV mort on ART');

%% HIV prevalance
% figure()
% for g = 1 : 2
%     artInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , g , 1 : age , 1 : risk));
%     artPop = sum(noV.popVec(: , artInds) , 2);
%     hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%         1 : intervens , g , 1 : age , 1 : risk));
%     allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%         1 : intervens , g , 1 : age , 1 : risk)); 
%     hivPop = sum(noV.popVec(: , hivInds) , 2);
%     allPop = sum(noV.popVec(: , allInds) , 2);
%     plot(tVec , 100 * (hivPop + artPop) ./ allPop)
%     hold on
% end
% xlabel('Year')
% ylabel('Prevalence')
% title('HIV Prevalence')
% 
% figure()
% artInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%     1 : intervens , 1 : gender , 1 : age , 1 : risk));
% artPop = sum(noV.popVec(: , artInds) , 2);
% hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%     1 : intervens , 1 : gender , 1 : age , 1 : risk));
% allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%     1 : intervens , 1 : gender , 1 : age , 1 : risk)); 
% hivPop = sum(noV.popVec(: , hivInds) , 2);
% allPop = sum(noV.popVec(: , allInds) , 2);
% plot(tVec , 100 * (hivPop + artPop) ./ allPop)
% 
% xlabel('Year')
% ylabel('Prevalence')
% title('HIV Prevalence All')
% 
% % Total HIV positive
% hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%     1 : intervens , 1 : 2 , 16 : 50 , 1 : risk));
% hivPop = sum(noV.popVec(: , hivInds) , 2);
% artInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%     1 : intervens , 1 : 2 , 16 : 50 , 1 : risk));
% art = sum(noV.popVec(: , artInds) , 2);
% % Compared to ANC data
% overallHivPrev_KZN_AC(1 , :) = 1990 : 2009;
% overallHivPrev_KZN_AC(2 , :) = [0.464072571
%     0.985438052
%     1.506803533
%     5.576509907
%     8.126044402
%     13.04177608
%     12.54905705
%     16.61876343
%     19.50632609
%     19.52064932
%     22.2391979
%     20.22439723
%     22.09787539
%     22.78825495
%     25.16877536
%     26.19622822
%     25.36548102
%     27.2380043
%     27.42134161
%     28.44974934];
% 
% % Compared to Africa Center (validation years)
% prevValYrs = 2010 : 2016;
% prevVal = [0.290
% 0.290
% 0.293
% 0.312
% 0.338
% 0.338
% 0.344
% ] .* 100;
% 
% upper_prevVal = [0.30
% 0.30
% 0.31
% 0.32
% 0.35
% 0.35
% 0.36
% ] .* 100;
% 
% lower_prevVal = [0.27
% 0.27
% 0.27
% 0.29
% 0.32
% 0.32
% 0.33] .* 100;
% 
% %figure()
% popTot = noV.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%     1 : 2 , 16 : 50, 1 : risk)));
% plot(tVec , (hivPop + art) ./ sum(popTot , 2) * 100 , overallHivPrev_KZN_AC(1 , :) , overallHivPrev_KZN_AC(2 , :) , '*')
% hold on 
% yPosError = abs(upper_prevVal - prevVal);
% yNegError = abs(lower_prevVal - prevVal);
% errorbar(prevValYrs , prevVal , yNegError , yPosError , 'ms')
% xlabel('Year'); ylabel('Proportion of Population (%)'); title('HIV Prevalence (Ages 15-49)')
% legend('Model' , 'National ANC Data' , 'Validation set: KZN Actual (Africa Center Data)' , 'Model Western Kenya' , 'Model Western Kenya' , 'Model Western Kenya')

%% HIV prevalence by gender vs. AC data
% prevYears = unique(hivPrevF_obs(: , 1));
% hivRaw(:,:,1) = hivPrevM_obs(: , 2:3);
% hivRaw(:,:,2) = hivPrevF_obs(: , 2:3);
% 
% hivData(: , : , 1) = zeros(length(prevYears) , 1);
% hivData(: , : , 2) = zeros(length(prevYears) , 1);
% 
% for i = 1 : length(prevYears)
%     for g = 1 : gender
%         hivData(i,1,g) = sumall(hivRaw(((i-1)*7+1):(i*7) , 1 , g)) ./ sumall(hivRaw(((i-1)*7+1):(i*7) , 2 , g)) .* 100;
%     end
% end
% 
% gen = {'Male' , 'Female'};
% for g = 1 : gender
%     hivInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         g , 16 : 50 , 1 : risk)); toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         g , 16 : 50 , 1 : risk))];
%     totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         g , 16 : 50 , 1 : risk));
%     hivPop = sum(noV.popVec(: , hivInds) , 2);
%     hivPopPrev = bsxfun(@rdivide , hivPop , sum(noV.popVec(: , totInds) , 2)) * 100;
%     subplot(1,2,g)
%     plot(tVec' , hivPopPrev);
%     hold all;
%     plot(prevYears , hivData(:,:,g) , 'ro');
%     xlabel('Year'); ylabel('Prevalence (%)'); title(gen{g});
%     %xlim([1980 2120])
%     legend('Model' , 'Africa Center Data (Calibration)')
% end

%% HIV prevalence by age on x-axis
% genderVec = {'Males (on and off ART)' , 'Females (on and off ART)'};
% hivAge = zeros(16 , 2);
% ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
%     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
%     '60-64' , '65-69' , '70-74' , '75-79'};
% 
% %figure;
% for g = 1 : gender
%     for a = 1 : age
%         ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : intervens , g , a , 1 : risk));
%         hivInds = toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : intervens , g , a , 1 : risk));
%         hivAge(a , g) = (sum(noV.popVec(end,hivInds),2)/sum(noV.popVec(end,ageInds),2))*100; %end-605
%     end
%     hold all;
%     subplot(2,2,g+2);
%     plot(1 : size(hivAge , 1) , hivAge(: , g) , ':');
%     hold all;
%     xlabel('Age Group'); ylabel('HIV Prevalence')
%     set(gca , 'xtick' , 1 : length(hivAge) , 'xtickLabel' , ageGroup);
%     title(genderVec{g})
%     legend('With ART dropout: 6.19%' , 'Without ART dropout' , 'With ART dropout: 11.8%' , 'With ART dropout: 11.8%, HIV mort on ART' , 'With ART dropout: 11.8%, 90% VS');
% end

%% CD4 count and VL proportions over time
% %figure;
% hold on;
% subplot(1,2,1);
% hold on;
% for d = 2 : 6
%     cd4 = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : intervens , 2 , 1 : age , 1 : risk));
%     hivAll = toInd(allcomb(2:6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : intervens , 2 , 1 : age , 1 : risk));
%         
%     hold all;
%     prop = (sum(noV.popVec(:,cd4),2)./sum(noV.popVec(:,hivAll),2))*100;
%     plot(tVec(:) , prop)
%     legend('Acute- w/DO' , 'CD4>500- w/DO' , 'CD4 500-350- w/DO' , 'CD4 350-200- w/DO' , 'CD4 <=200- w/DO' , ...
%         'Acute- noDO' , 'CD4>500- noDO' , 'CD4 500-350- noDO' , 'CD4 350-200- noDO' , 'CD4 <=200- noDO');
%     ylabel('Percentage (%)');
%     xlabel('Year');
%     grid on;
% end
% hold on;
% subplot(1,2,2);
% hold on;
% for v = 1 : 5
%     vL = toInd(allcomb(2 : 6 , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : intervens , 2 , 1 : age , 1 : risk));
%     hivAll = toInd(allcomb(2:6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : intervens , 2 , 1 : age , 1 : risk));
%         
%     hold all;
%     prop = (sum(noV.popVec(:,vL),2)./sum(noV.popVec(:,hivAll),2))*100;
%     plot(tVec(:) , prop)
%     legend('Acute- w/DO' , 'VL<1000- w/DO' , 'VL 1000-10000- w/DO' , 'VL 10000-50000- w/DO' , 'VL>50000- w/DO' , ...
%         'Acute- noDO' , 'VL<1000- noDO' , 'VL 1000-10000- noDO' , 'VL 10000-50000- noDO' , 'VL>50000- noDO');
%     ylabel('Percentage (%)');
%     xlabel('Year');
%     grid on;
% end

%% VL proportion by age
% genderVec = {'Males (on and off ART)' , 'Females (on and off ART)'};
% hivAge = zeros(16 , 2 , 5);
% ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
%     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
%     '60-64' , '65-69' , '70-74' , '75-79'};
% 
% %figure;
% for g = 1 : gender
%     for a = 1 : age
%         for v = 1 : 5
%             vL = toInd(allcomb(2 : 6 , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : intervens , 2 , a , 1 : risk));
%             hivAll = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%                     1 : intervens , 2 , a , 1 : risk));
% 
%             prop = (sum(noV.popVec(end-713,vL),2)./sum(noV.popVec(end-713,hivAll),2))*100;
%             hivAge(a , g , v) = prop;
%         end
%     end
%     hold all;
%     subplot(2,2,g);
%     for v = 1 : 5
%         plot(1 : size(hivAge , 1) , hivAge(: , g , v) , ':');
%         hold all;
%     end
%     hold all;
%     xlabel('Age Group'); ylabel('VL')
%     set(gca , 'xtick' , 1 : length(hivAge) , 'xtickLabel' , ageGroup);
%     title(genderVec{g})
%     legend('With 11.8% DO: Acute' , 'VL<1000' , 'VL 1000-10000' , 'VL 10000-50000' , 'VL>50000');
%     % 'WithOUT DO: Acute' , 'VL<1000' , 'VL 1000-10000' , 'VL 10000-50000' , 'VL>50000' , ...
% end

%% HIV Mortality by age
% aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80}; %{10:15,16:25,26:35,36:50,51:75};
% ageGroup = {'0-4','5-9','10-14' , '15-19' , '20-24' , '25-29' ,...
%      '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
%      '60-64' , '65-69' , '70-74' , '75-79'};
% aMatrix = zeros(1 , length(aVec));
% fac = 10 ^ 5;
% 
%     for aInd = 1 : length(aVec)
%         a = aVec{aInd};
% 
%         % General
%         allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : intervens , 2 , a , 1 : risk));
%         
%         % Calculate mortality
%         for n = nSims %1 : nSims
%             hivMortRef = ...
%                 (annlz(sum(vaxResult{n}.hivDeaths(end-5:end , 2 , a),3)) ./ ... % end-5:end
%                 (annlz(sum(vaxResult{n}.popVec(end-5:end , allF) , 2) ./ stepsPerYear)) * fac); %553:553+5  , 661:661+5
%             
%             aMatrix(1,aInd) = hivMortRef;
%         end
%     end
%     hold all;    
%     plot([1:length(aVec)] , aMatrix(1,:))
%     hold all;
%     ylabel('HIV Mortality rates in women (per 100K)'); 
%     set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
%     legend('Year: 2002 (before ART)' , 'Year: 2020 (after ART, current year)' , 'Year: 2120 (end year)');

%% On ART by age
% aVec = {1:5,6:1011:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80}; %{10:15,16:25,26:35,36:50,51:75};
% ageGroup = {'0-4','5-9' ,'10-14' , '15-19' , '20-24' , '25-29' ,...
%      '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
%      '60-64' , '65-69' , '70-74' , '75-79'};
% aMatrix = zeros(1 , length(aVec));
% for aInd = 1 : length(aVec)
%     a = aVec{aInd};
% 
%     artInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , 2 , a , 1 : risk));
%     artPop = sum(vaxResult{noVaxInd}.popVec(end , artInds) , 2);
%     hivInds = toInd(allcomb([2 : 6,10] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%         1 : intervens , 2 , a , 1 : risk));
%     hivPop = sum(vaxResult{noVaxInd}.popVec(end , hivInds) , 2);
%     hiv_art = [100 * artPop ./ hivPop];
% 
%     aMatrix(aInd) = hiv_art;
% end
% 
% figure;
% plot([1:length(aVec)] , aMatrix(1,:))
% hold all;
% ylabel('Percent females on ART');
% ylim([0 110])
% set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);

%% HIV incidence
% fac = 10 ^ 5;
% 
% figure();
% 
% for a = 1 : age
%     % All HIV-negative women
%     hivNeg = toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         2 , a , 1 : risk));
% 
%     noV.hivInc = ...
%         (annlz(sum(noV.newHiv(: , 2 , a , 1 : risk),4)) ./ (annlz(sum(noV.popVec(: , hivNeg) , 2) ./ stepsPerYear))* fac);
% 
%     subplot(4,4,a)
%     plot(tVec(1 : stepsPerYear : end) , noV.hivInc);
%     title(' HIV Incidence')
%     xlabel('Year'); ylabel('Incidence per 100,000')
%     %xlim([2018 2050])
%     hold all;
% end
% 
% % legend('lr','mr','hr')

%% HIV Mortality by age
% aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80}; %{10:15,16:25,26:35,36:50,51:75};
% ageGroup = {'0-4','5-9','10-14' , '15-19' , '20-24' , '25-29' ,...
%      '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
%      '60-64' , '65-69' , '70-74' , '75-79'};
% aMatrix = zeros(1 , length(aVec));
% fac = 10 ^ 5;
% 
%     for a = 1 : age
%         % General
%         allF = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : intervens , 2 , a , 1 : risk));
%         
%         % Calculate mortality
%         hivMortRef = ...
%             (annlz(sum(noV.hivDeaths(end-5:end , 2 , a),3)) ./ ... % end-5:end
%             (annlz(sum(noV.popVec(end-5:end , allF) , 2) ./ stepsPerYear)) * fac); %553:553+5  , 661:661+5
% 
%         aMatrix(1,a) = hivMortRef;
%     end
%     hold all;    
%     plot([1:length(aVec)] , aMatrix(1,:))
%     hold all;
%     ylabel('HIV Mortality rates in women (per 100K)'); 
%     set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
%     %legend('Year: 2002 (before ART)' , 'Year: 2020 (after ART, current year)' , 'Year: 2120 (end year)');
%     legend('Year: 2020 (end year)');
    
%% ART 'incidence'
% fac = 10 ^ 5;
% ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
%     '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
%     '60-64' , '65-69' , '70-74' , '75-79'};
% 
% %figure();
% hold all;
% 
% for a = 1 : age
%     % All HIV-positive women
%     hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         2 , a , 1 : risk));
% 
%     noV.artInc = ...
%         annlz(sum(sum(sum(noV.artTreatTracker(: , 2 : 6 , 1 : viral , 2 , a , 1 : risk),2),3),6));% ./ (annlz(sum(noV.popVec(: , hivPos) , 2) ./ stepsPerYear))* fac);
% 
%     subplot(4,4,a)
%     plot(tVec(1 : stepsPerYear : end) , noV.artInc , '-');
%     title(['Ages: ' , ageGroup{a}])
%     xlabel('Year'); ylabel('New females on ART')
%     %xlim([2000 2120])
%     hold all;
% end

%% ART coverage
% figure()
% popTot = noV.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%     1 : gender , 3 : age , 1 : risk)));
% artInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%     1 : intervens , 1 : gender , 3 : age , 1 : risk));
% artPop = sum(noV.popVec(: , artInds) , 2);
% hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%     1 : intervens , 1 : gender , 3 : age , 1 : risk));
% hivPop = sum(noV.popVec(: , hivInds) , 2);
% hiv_art = [100 * hivPop ./ sum(popTot , 2), 100 * artPop ./ sum(popTot , 2)];
% % area(tVec , hiv_art); %art ./ sum(popVec , 2) , tVec , hiv ./ sum(popVec , 2))
% % xlabel('Year')
% % ylabel('Proportion of Population (%)')
% % title('Relative HIV Prevalence')
% % legend('Untreated', 'On ART' , 'Location' , 'NorthWest')
% 
% %figure()
% hold all;
% artActual = [0	0	0	0	1	2	3	6, ...
%     9	14	19	27	34	40	45	48];
% yrsArtActual = [2000	2001	2002	2003	2004	2005 ...
%     2006	2007	2008	2009	2010	2011	2012	2013 ...
%     2014	2015];
% plot(tVec , 100 * artPop ./ (hivPop + artPop) , '--'); % , yrsArtActual , artActual , '*')
% xlabel('Year')
% ylabel('Proportion of HIV Population')
% title('Proportion on ART')
% legend('Model' , 'Observed')
% xlim([2000 2120])

%% HIV by age group
% hivAge = zeros(length(tVec) , 12);
% for a = 1 : age
%     hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , 1 : gender , a , 1 : risk));
%     hivArt = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , 1 : gender , a , 1 : risk));
%     hivNeg = toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , 1 : gender , a , 1 : risk));
%     hivAge(: , a) = sum(noV.popVec(: , hivPos) , 2) + sum(noV.popVec(: , hivArt) , 2);
% end
% hivPosAllInd = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%     1 : gender , 1 : age , 1 : risk));
% hivArtAllInd = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%     1 : gender , 1 : age , 1 : risk));
% hivNegAllInd = toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , 1 : gender , 1 : age , 1 : risk));
% hivPosAll = sum(noV.popVec(: , hivPosAllInd) , 2 ) + sum(noV.popVec(: , hivArtAllInd),2);
% 
% figure()
% subplot(1 , 2 , 1)
% area(tVec , bsxfun(@rdivide , hivAge , hivPosAll));
% title('HIV Status by Age Group')
% xlabel('Year')
% ylabel('Proportion of HIV Positive')
% legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
%     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
%     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')
% 
% % HIV by risk group
% hivRisk = zeros(length(tVec) , risk);
% for r = 1 : risk
%     hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , 1 : gender , 1 : age , r));
%     hivArt = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , 1 : gender , 1 : age , r));
%     hivRisk(: , r) = sum(noV.popVec(: , hivPos) , 2) + sum(noV.popVec(: , hivArt) , 2);
% 
% end
% subplot(1 , 2 , 2)
% area(tVec , bsxfun(@rdivide , hivRisk , hivPosAll));
% title('HIV Status by Risk Group')
% xlabel('Year')
% ylabel('Proportion of HIV Positive')
% legend('Low' , 'Medium' , 'High' , 'Location' , 'NorthEastOutside')

%% ART treatment tracker- CD4/risk/age
% cd4ARTFrac = zeros(length(tVec) , risk);
% for i = 1 : length(tVec)
%     %currTot = sumall(noV.artTreatTracker(i , 2 : 6 , 1 : 5 , 2 , 3 , 1 : risk));
%     for c = 2 : 6
%         curr = sumall(noV.artTreatTracker(i , c , 1 : viral , 2 , 3 : age , 1 : risk));
%         %hivArt = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         %    1 : intervens , 2 , 3 , r));
%         %curr = sum(noV.popVec(i , hivArt));
% %         hivPos = toInd(allcomb(c , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %             1 : intervens , 2 , 4 , 1 : risk));
%         hivPos = sumall(noV.artTreatTracker(i , 2 : 6 , 1 : viral , 2 , 3 : age , 1 : risk));
%         %cd4ARTFrac(i , r) = (curr / (sum(noV.popVec(i , hivPos))+curr)) .* 100;
% %         cd4ARTFrac(i , c-1) = (curr / (sum(noV.popVec(i , hivPos)))).* 100;
%         cd4ARTFrac(i , c-1) = (curr / hivPos).* 100;
%     end
% end
% 
% figure()
% plot(tVec , cd4ARTFrac)
% legend('Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
%    'CD4 <= 200 cells/uL' , 'Location' , 'NorthEastOutside')
% %legend('1' , '2' , '3' , '4' , '5' , '6' , '7' , '8' , '9' , '10' , 'Location' , 'NorthEastOutside')
% xlabel('Year')
% ylabel('Initiated ART')

%% HPV prevalence by HIV group and HPV type
% Vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , 2 , [1 : 1 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(vaxResult{2}.popVec(: , ccHivInds) , 2);
popHivTot = vaxResult{2}.popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , 2 , [1 : 1 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(vaxResult{2}.popVec(: , ccArtInds) , 2);
popArtTot = vaxResult{2}.popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , 2 , [1 : 1 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(vaxResult{2}.popVec(: , ccHivNegInds) , 2);
popHivNegTot = vaxResult{2}.popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

figure();
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'-')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'-')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'-')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' HPV Prevalence')
%legend('HIV-' , 'HIV+ noART' , 'ART')
hold all;

% Non-vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , ...
     [1 : 1 , 7] , 2 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(vaxResult{2}.popVec(: , ccHivInds) , 2);
popHivTot = vaxResult{2}.popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     [1 : 1 , 7] , 2 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(vaxResult{2}.popVec(: , ccArtInds) , 2);
popArtTot = vaxResult{2}.popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     [1 : 1 , 7] , 2 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(vaxResult{2}.popVec(: , ccHivNegInds) , 2);
popHivNegTot = vaxResult{2}.popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

hold all;
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'--')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'--')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'--')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' HPV Prevalence')
legend('HIV-' , 'HIV+ noART' , 'ART' , 'HIV-, nonVax' , 'HIV+ noART, nonVax' , 'ART, nonVax')
hold all;

%% CIN1 prevalence by HIV group and HPV type
% Vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , 3 , [1 : 3 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(vaxResult{2}.popVec(: , ccHivInds) , 2);
popHivTot = vaxResult{2}.popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , 3 , [1 : 3 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(vaxResult{2}.popVec(: , ccArtInds) , 2);
popArtTot = vaxResult{2}.popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , 3 , [1 : 3 , 7] , ...
     1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(vaxResult{2}.popVec(: , ccHivNegInds) , 2);
popHivNegTot = vaxResult{2}.popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

figure();
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'-')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'-')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'-')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' CIN1 Prevalence')
%legend('HIV-' , 'HIV+ noART' , 'ART')
hold all;

% Non-vaccine-type HPV
% HIV+
ccHivInds = toInd(allcomb(3 : 7 , 1 : viral , ...
     [1 : 3 , 7] , 3 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivPop = sum(vaxResult{2}.popVec(: , ccHivInds) , 2);
popHivTot = vaxResult{2}.popVec(: , toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%ART
ccArtInds = toInd(allcomb(8 , 6 , ...
     [1 : 3 , 7] , 3 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccArtPop = sum(vaxResult{2}.popVec(: , ccArtInds) , 2);
popArtTot = vaxResult{2}.popVec(: , toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));
%HIV-
ccHivNegInds = toInd(allcomb(1 : 2 , 1 , ...
     [1 : 3 , 7] , 3 , 1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk));
ccHivNegPop = sum(vaxResult{2}.popVec(: , ccHivNegInds) , 2);
popHivNegTot = vaxResult{2}.popVec(: , toInd(allcomb(1 : 2 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 3 : age , 1 : risk)));

hold all;
% plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'--')
hold all
plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'--')
hold all
plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'--')
%axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Prevalence (%)'); title(' CIN1 Prevalence')
legend('HIV-' , 'HIV+ noART' , 'ART' , 'HIV-, nonVax' , 'HIV+ noART, nonVax' , 'ART, nonVax')
hold all;

%% HPV prevalence by HIV group
% %figure()
% linStyle = {'--' , '-' , ':'};
% %for a = 1 : age
% %    for r = 1 : risk
%         % HIV+
%         hpvHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 2:4 , 4 , ...
%              1 : intervens , 2 : gender , 10:75 , 1:risk));
%         hpvHivPop = sum(vaxResult{3}.popVec(: , hpvHivInds) , 2);
%         popHivTot = vaxResult{3}.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates ,  1 : intervens , ...
%             2 : gender , 10:75 , 1:risk)));
%         %ART
%         hpvArtInds = toInd(allcomb(10 , 6 , 2:4 , 4 , ...
%              1 : intervens , 2 : gender , 10:75 , 1:risk));
%         hpvArtPop = sum(vaxResult{3}.popVec(: , hpvArtInds) , 2);
%         popArtTot = vaxResult{3}.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%             2 : gender , 10:75 , 1:risk)));
%         %HIV ALL
%         hpvAllInds = toInd(allcomb([2:6,10] , 1 : viral , 2:4 , 4 , ...
%              1 : intervens , 2 : gender , 10:75 , 1:risk));
%         hpvAllPop = sum(vaxResult{3}.popVec(: , hpvAllInds) , 2);
%         popAllTot = vaxResult{3}.popVec(: , toInd(allcomb([2:6,10] , 1:viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%             2 : gender , 10:75 , 1:risk)));
%         %HIV-
%         hpvHivNegInds = toInd(allcomb([1,7:9] , 1 , 2:4 , 4 , ...
%              1 : intervens , 2 : gender , 10:75 , 1:risk));
%         hpvHivNegPop = sum(vaxResult{3}.popVec(: , hpvHivNegInds) , 2);
%         popHivNegTot = vaxResult{3}.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates ,  1 : intervens , ...
%             2 : gender , 10:75 , 1:risk)));
%         
%         %General
%         genInds = toInd(allcomb(1 : disease , 1 : viral , 2:4 , 4 , ...
%              1 : intervens , 2 : gender , 10:75 , 1:risk));
%         genPop = sum(vaxResult{3}.popVec(: , genInds) , 2);
%         genTot = vaxResult{3}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates ,  1 : intervens , ...
%             2 : gender , 10:75 , 1:risk)));
% 
%         %subplot(4,4,a)
%         hold all;
%         plot(tVec , 100 * genPop ./ sum(genTot , 2) , '-')
%         hold all;
%         plot(tVec , 100 * hpvHivNegPop ./ sum(popHivNegTot , 2) , '-')
%         hold all;
%         plot(tVec , 100 * hpvHivPop ./ sum(popHivTot , 2) , '-')
%         hold all;
%         plot(tVec , 100 * hpvArtPop ./ sum(popArtTot , 2) , '-')
%         hold all;
%         plot(tVec , 100 * hpvAllPop ./ sum(popAllTot , 2) , '-')
%         %axis([tVec(1) tVec(end) 0 100])
%         xlim([1950 2120])
%         ylim([0 20])
%         xlabel('Year'); ylabel('Prevalence (%)'); title(' CIN3 Prevalence, ages 9-74')
%         legend('General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all');
% %    end
% %end
% %legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')

%% HPV/CIN/CC prevalence by HIV group
% figure()
% linStyle = {'--' , '-' , ':'};
% for a = 1 : age
%     for r = 1 : risk
%         % HIV+
%         hpvHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 2 : 4 , 1:7, ...
%              1 : intervens , 2 , a , r));
%         hpvHivPop = sum(noV.popVec(: , hpvHivInds) , 2);
%         popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates ,  1 : intervens , ...
%             2 , a , r)));
%         %ART
%         hpvArtInds = toInd(allcomb(10 , 6 , 2 : 4 , 1:7, ...
%              1 : intervens , 2 , a , r));
%         hpvArtPop = sum(noV.popVec(: , hpvArtInds) , 2);
%         popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%             2 , a , r)));
%         %HIV-
%         hpvHivNegInds = toInd(allcomb([1,7:9] , 1 , 2 : 4 , 1:7 , ...
%              1 : intervens , 2 , a , r));
%         hpvHivNegPop = sum(noV.popVec(: , hpvHivNegInds) , 2);
%         popHivNegTot = noV.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates ,  1 : intervens , ...
%             2 , a , r)));
% 
%         subplot(4,4,a)
%         plot(tVec , 100 * hpvHivNegPop ./ sum(popHivNegTot , 2),linStyle{r})
%         hold all
%         plot(tVec , 100 * hpvHivPop ./ sum(popHivTot , 2),linStyle{r})
%         hold all
%         plot(tVec , 100 * hpvArtPop ./ sum(popArtTot , 2),linStyle{r})
%         %axis([tVec(1) tVec(end) 0 100])
%         xlim([1980 2100])
%         xlabel('Year'); ylabel('Prevalence (%)'); title(' CC Prevalence all ages')
%     end
% end
% legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')

%% CC prevalence by HIV group
% % HIV+
% ccHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 2 : 4 , 5 : 7, ...
%      [1,6] , 2 , 3 : age , 1 : risk));
% ccHivPop = sum(noV.popVec(: , ccHivInds) , 2);
% popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , [1:8,10] ,  [1,6] , ...
%     2 , 3 : age , 1 : risk)));
% %ART
% ccArtInds = toInd(allcomb(10, 6 , 2 : 4 , 5 : 7, ...
%      [1,6] , 2 , 3 : age , 1 : risk));
% ccArtPop = sum(noV.popVec(: , ccArtInds) , 2);
% popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , [1:8,10] ,  [1,6] , ...
%     2 , 3 : age , 1 : risk)));
% %HIV-
% ccHivNegInds = toInd(allcomb(1 , 1 , 2 : 4 , 5 : 7, ...
%      [1,6] , 2 , 3 : age , 1 : risk));
% ccHivNegPop = sum(noV.popVec(: , ccHivNegInds) , 2);
% popHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvVaxStates , [1:8,10] ,  [1,6] , ...
%     2 , 3 : age , 1 : risk)));
% 
% figure();
% plot(tVec , 100 * ccHivNegPop ./ sum(popHivNegTot , 2),'o')
% hold all
% plot(tVec , 100 * ccHivPop ./ sum(popHivTot , 2),'o')
% hold all
% plot(tVec , 100 * ccArtPop ./ sum(popArtTot , 2),'o')
% %axis([tVec(1) tVec(end) 0 100])
% xlabel('Year'); ylabel('Prevalence (%)'); title(' CC Prevalence')
% legend('HIV-' , 'HIV+ noART' , 'ART')

%% Vaccinated proportion

% By scenario
figure;
vaxInds = [toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , [1,6] , 2 , 10 : 75 , 1 : risk)); ...
    toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , 10 : 75 , 1 : risk))];
vaxPop = sum(vaxResult{1}.popVec(: , vaxInds) , 2);
popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
    2 , 10 : 75 , 1 : risk)));

hold all;
plot(tVec , 100 * vaxPop ./ sum(popTot , 2))

axis([2020 2120 0 100])
xlabel('Year'); ylabel('Coverage (%)');
legend('Phase I, Scenario 1' , 'Phase 2, Scenario 5');

% By age
figure;
aVec = {10:15,16:25,26:35,36:50,51:75};
saveResults = [];
for aInd = 1:length(aVec)
    a = aVec{aInd};
    vaxInds = [toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , [1,6] , 2 , a , 1 : risk)); ...
    toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , a , 1 : risk))];
    vaxPop = sum(vaxResult{1}.popVec(: , vaxInds) , 2);
    popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , 1 : risk)));
    
    covAge = 100 * vaxPop ./ sum(popTot , 2);
    saveResults = [saveResults , covAge(1 : stepsPerYear : end)];
    hold all;
    plot(tVec , 100 * vaxPop ./ sum(popTot , 2))
    axis([2010 2120 0 100])
    xlabel('Year'); ylabel('Coverage (%)');
    %legend('Phase I, Scenario 1' , 'Phase 2, Scenario 5');
    legend('9-14' , '15-24' , '25-34' , '35-49' , '50-74');
end
saveResults = [tVec(1 : stepsPerYear : end)' , saveResults];
fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\' , 'vaxCoverage' , '.xlsx'];
xlswrite(fname , saveResults)

% Age-standardized
figure;
inds = {[1 : disease] , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
worldStandard_WP2015 = [325428 (311262/5.0) 295693 287187 291738 299655 272348 ...
        247167 240167 226750 201603 171975 150562 113118 82266 64484];
aVec = {1:5,10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};

for i = 1 : length(inds)    
    vaxResult{1}.propVaxRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
    
    for aInd = 2:(age/5)-1
        a = aVec{aInd};
        
        % General
        allFVax = [toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , [1,6] , 2 , a , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , a , 1 : risk))];
        % All HIV-negative women
        hivNegVax = [toInd(allcomb([1,7:9] , 1 : viral , 1 , 9 , [1,6] , 2 , a , 1 : risk)); ...
            toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , a , 1 : risk))];
        % HIV-positive women not on ART
        hivNoARTFVax = [toInd(allcomb([2 : 6] , 1 : viral , 1 , 9 , [1,6] , 2 , a , 1 : risk)); ...
            toInd(allcomb([2 : 6] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , a , 1 : risk))];
        % Women on ART
        artFVax = [toInd(allcomb(10 , 1 : viral , 1 , 9 , [1,6] , 2 , a , 1 : risk)); ...
            toInd(allcomb(10 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , a , 1 : risk))];
        % All HIV-positive women
        hivAllFVax = [toInd(allcomb([2:6,10] , 1 : viral , 1 , 9 , [1,6] , 2 , a , 1 : risk)); ...
            toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , a , 1 : risk))];
        genArrayVax = {allFVax , hivNegVax , hivNoARTFVax , artFVax , hivAllFVax};
        
        vaxPop = sum(vaxResult{1}.popVec(: , genArrayVax{i}) , 2);
        popTot = sum(vaxResult{1}.popVec(: , toInd(allcomb(inds{i} , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
            2 , a , 1 : risk))) , 2);
        
        vaxResult{1}.propVaxRef = vaxResult{1}.propVaxRef + (vaxPop ./ popTot) .* worldStandard_WP2015(aInd);
    end
    
    vaxResult{1}.propVax = (vaxResult{1}.propVaxRef ./ (sum(worldStandard_WP2015(2:(age/5)-1)))) .* 100;
    
    hold all;
    plot(tVec , vaxResult{1}.propVax)
    axis([2010 2120 0 90])
    xlabel('Year'); ylabel('Coverage (%)'); 
    legend('General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all');
end
    
%% Vaccinated proportion by HIV group
figure();
linStyle = {'--' , '-' , ':'};
for a = 1 : age
    for r = 1 : risk
    % HIV+
    vaxHivInds = [toInd(allcomb(2 : 6 , 1 : 5 , 1 , 9 , [1,6] , 2 , a , r)); ...
        toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , a , r))];
    vaxHivPop = sum(vaxResult{2}.popVec(: , vaxHivInds) , 2);
    popHivTot = vaxResult{2}.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , r)));
    %ART
    vaxArtInds = [toInd(allcomb(10 , 6 , 1 , 9 , [1,6] , 2 , a , r)); ...
        toInd(allcomb(10, 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , a , r))];
    vaxArtPop = sum(vaxResult{2}.popVec(: , vaxArtInds) , 2);
    popArtTot = vaxResult{2}.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , r)));
    %HIV-
    vaxHivNegInds = [toInd(allcomb([1,7:9] , 1 , 1 , 9 , [1,6] , 2 , a , r)); ...
    toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , a , r))];
    vaxHivNegPop = sum(vaxResult{2}.popVec(: , vaxHivNegInds) , 2);
    popHivNegTot = vaxResult{2}.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , r)));

    subplot(4,4,a)
    plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2) , linStyle{r})
    hold all
    plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2) , linStyle{r})
    hold all
    plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2) , linStyle{r})
    %axis([tVec(1) tVec(end) 0 100])
    xlabel('Year'); ylabel('Proportion (%)'); title('Vaccinated Proportion')
    xlim([2018 2100])
    
    hold all;
    end
end
legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')

%% Vaccinated proportion by CD4 count
figure();
linStyle = {'--' , '-' , ':'};
for c = 2 : 6
% HIV+
vaxHivInds = [toInd(allcomb(c , 1 : 5 , 1 , 9 , [1,6] , 2 , 7 , 1 : risk)); ...
    toInd(allcomb(c , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , 7 , 1 : risk))];
vaxHivPop = sum(vaxResult{2}.popVec(: , vaxHivInds) , 2);
popHivTot = vaxResult{2}.popVec(: , toInd(allcomb(c , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
    2 , 7 , 1 : risk)));
%ART
vaxArtInds = [toInd(allcomb(10 , 6 , 1 , 9 , [1,6] , 2 , 7 , 1 : risk)); ...
   toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , 7 , 1 : risk))];
vaxArtPop = sum(vaxResult{2}.popVec(: , vaxArtInds) , 2);
popArtTot = vaxResult{2}.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
    2 , 7 , 1 : risk)));
%HIV-
vaxHivNegInds = [toInd(allcomb([1,7:9] , 1 , 1 , 9 , [1,6] , 2 , 7 , 1 : risk)); ...
    toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , 7 , 1 : risk))];
vaxHivNegPop = sum(vaxResult{2}.popVec(: , vaxHivNegInds) , 2);
popHivNegTot = vaxResult{2}.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
    2 , 7 , 1 : risk)));

% hold all
% plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2),'k--')
% hold all
plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2),'m:')
hold all
% plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2),':')
% hold all
axis([tVec(1) tVec(end) 0 100])
xlabel('Year'); ylabel('Proportion (%)'); title('Vaccinated Proportion- age group 7')
xlim([2018 2100])
end

% legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')
legend('HIV-','ART','Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
   'CD4 <= 200 cells/uL')

%% Vaccinated proportion by risk
figure()
for r = 1 : risk
    %HIV-
    vaxHivNegInds = [toInd(allcomb([1,7:9] , 1 , 1 , 9 , [1,6] , 2 , 3 , r)); ...
    toInd(allcomb([1,7:9] , 1 , 2 : 4 , 1 : hpvNonVaxStates , [2,4] , 2 , 3 , r))];
    vaxHivNegPop = sum(vaxResult{1}.popVec(: , vaxHivNegInds) , 2);
    popHivNegTot = vaxResult{1}.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , 3 , r)));
    % HIV+
    vaxHivInds = [toInd(allcomb(2 : 6 , 1 : 5 , 1 , 9 , [1,6] , 2 , 3 , r)); ...
        toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , 2 , 3 , r))];
    vaxHivPop = sum(vaxResult{1}.popVec(: , vaxHivInds) , 2);
    popHivTot = vaxResult{1}.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , 3 , r)));
    %ART
    vaxArtInds = [toInd(allcomb(10, 6 , 1 , 9 , [1,6] , 2 , 3 , r)); ...
        toInd(allcomb(10, 6 , 2 : 4 , 1 : hpvNonVaxStates , [2,4] , 2 , 3 , r))];
    vaxArtPop = sum(vaxResult{1}.popVec(: , vaxArtInds) , 2);
    popArtTot = vaxResult{1}.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , 3 , r)));

    plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2),'o')
    hold all
    plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2),'o')
    hold all
    plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2),'o')
    %axis([tVec(1) tVec(end) 0 100])
    xlabel('Year'); ylabel('Proportion (%)'); title('Vaccinated Proportion')
    legend('lr HIV-' , 'lr HIV+' , 'lr ART' , 'mr HIV-' , 'mr HIV+' , 'mr ART' , 'hr HIV-' , 'hr HIV+' , 'hr ART')
end

%% Vaccinateable prevalence by HIV group
% figure();
% linStyle = {'--' , '-' , ':'};
% for a = 1 : age
%     for r = 1 : risk
%         % HIV+
%         vaxHivInds = [toInd(allcomb(2 : 6 , 1 : 5 , 1 , 1 , [1,6] , 2 , a , r)); ...
%             toInd(allcomb(2 : 6 , 1 : 5 , 2:4 , 10 , [1,6] , 2 , a , r))];
%         vaxHivPop = sum(noV.popVec(: , vaxHivInds) , 2);
%         popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%             2 , a , r)));
%         %ART
%         vaxArtInds = [toInd(allcomb(10, 6 , 1 , 1 , [1,6] , 2 , a , r)); ...
%             toInd(allcomb(10, 6 , 2:4 , 10 , [1,6] , 2 , a , r))];
%         vaxArtPop = sum(noV.popVec(: , vaxArtInds) , 2);
%         popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%             2 , a , r)));
%         %HIV-
%         vaxHivNegInds = [toInd(allcomb([1,7:9] , 1 , 1 , 1 , [1,6] , 2 , a , r)); ...
%             toInd(allcomb([1,7:9] , 1 , 2:4 , 10 , [1,6] , 2 , a , r))];
%         vaxHivNegPop = sum(noV.popVec(: , vaxHivNegInds) , 2);
%         popHivNegTot = noV.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%             2 , a , r)));
% 
%         subplot(4,4,a)
%         plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2),linStyle{r})
%         hold all
%         plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2),linStyle{r})
%         hold all
%         plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2),linStyle{r})
%         %axis([tVec(1) tVec(end) 0 100])
%         xlabel('Year'); ylabel('Prevalence (%)'); title(' Vaccinateable Prevalence')
%         ylim([0 100])
%     end
% end
% legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')

%% Population by risk and HIV group
figure();
for a = 1 : age
for r = 1 : risk
    % HIV+
    rpopHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , r)));
    popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , 1 : risk)));
    %ART
    rpopArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , r)));
    popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , 1 : risk)));
    %HIV-
    rpopHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , r)));
    popHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , 1 : risk)));

    subplot(4,4,a)
    plot(tVec , 100 * sum(rpopHivNegTot , 2) ./ sum(popHivNegTot , 2))
    hold all
    plot(tVec , 100 * sum(rpopHivTot , 2) ./ sum(popHivTot , 2))
    hold all
    plot(tVec , 100 * sum(rpopArtTot , 2) ./ sum(popArtTot , 2))
    xlabel('Year'); ylabel('Prevalence (%)'); title(' Risk Proportion')
    hold all;
end
end
legend('HIV- : lr' , 'HIV+ noART : lr' , 'ART : lr' , 'HIV- : mr' , 'HIV+ noART : mr' , 'ART : mr' , 'HIV- : hr' , 'HIV+ noART : hr'  , 'ART : hr')

%% Population by age and HIV group
% for a = 3 : 6
%     % HIV+
%     rpopHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         2 , a , 1 : risk)));
%     popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         2 , 3 : age , 1 : risk)));
%     %ART
%     rpopArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         2 , a , 1 : risk)));
%     popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         2 , 3 : age , 1 : risk)));
%     %HIV-
%     rpopHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         2 , a , 1 : risk)));
%     popHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
%         2 , 3 : age , 1 : risk)));
% 
%     figure();
%     % plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
%     plot(tVec , 100 * sum(rpopHivNegTot , 2) ./ sum(popHivNegTot , 2))
%     hold all
%     plot(tVec , 100 * sum(rpopHivTot , 2) ./ sum(popHivTot , 2))
%     hold all
%     plot(tVec , 100 * sum(rpopArtTot , 2) ./ sum(popArtTot , 2))
%     xlabel('Year'); ylabel('Prevalence (%)'); title([num2str(a) , ': Age Distribution'])
%     legend('HIV-' , 'HIV+ noART' , 'ART' )
%     hold all;
% end
    
%% Population by "p"
figure();
subplot(2,2,1);
for p = 1 : intervens
    % General
    inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
         p , 1 : gender , 8 , 1 : risk));
    pop = sum(vaxResult{1}.popVec(: , inds) , 2);
    popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
         1 : intervens , 1 : gender , 8 , 1 : risk)));
    plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
    xlabel('Year'); ylabel('Proportion (%)'); title(' p Proportion')
    legend('1' , '2' , '3' , '4' , '5' ,'6')
    hold all;
end

subplot(2,2,3);
inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
     [1,6] , 1 : gender , 8 , 1 : risk));
pop = sum(vaxResult{1}.popVec(: , inds) , 2);
popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
     1 : intervens , 1 : gender , 8 , 1 : risk)));
plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
xlabel('Year'); ylabel('Proportion (%)'); title(' Proportion unvaccinated or vaccinated and not reinfected ')

subplot(2,2,4);
inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
     [2,4] , 1 : gender , 8 , 1 : risk));
pop = sum(vaxResult{1}.popVec(: , inds) , 2);
popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
     1 : intervens , 1 : gender , 8 , 1 : risk)));
plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
xlabel('Year'); ylabel('Proportion (%)'); title(' Proportion vaccinated and reinfected ')

subplot(2,2,2);
inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
     [4,6] , 1 : gender , 8 , 1 : risk));
pop = sum(vaxResult{1}.popVec(: , inds) , 2);
popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
     1 : intervens , 1 : gender , 8 , 1 : risk)));
plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
xlabel('Year'); ylabel('Proportion (%)'); title(' Proportion screened')

%% Screened proportion by HIV group
figure();
linStyle = {'--' , '-' , ':'};
for a = 1 : age
    for r = 1 : risk
    % HIV+
    vaxHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 4 : 6 , 2 , a , r));
    vaxHivPop = sum(noV.popVec(: , vaxHivInds) , 2);
    popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , r)));
    %ART
    vaxArtInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 4 : 6 , 2 , a , r));
    vaxArtPop = sum(noV.popVec(: , vaxArtInds) , 2);
    popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , r)));
    %HIV-
    vaxHivNegInds = toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 4 : 6 , 2 , a , r));
    vaxHivNegPop = sum(noV.popVec(: , vaxHivNegInds) , 2);
    popHivNegTot = noV.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , r)));

    subplot(4,4,a)
    plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2) , linStyle{r})
    hold all
    plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2) , linStyle{r})
    hold all
    plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2) , linStyle{r})
    %axis([tVec(1) tVec(end) 0 100])
    xlabel('Year'); ylabel('Proportion (%)'); title('Screened Proportion')
%     xlim([2018 2100])
    
    hold all;
    end
end
legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')

%%
% hold on
% for g = 1 : 2
%     artInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , g , 4 : 10 , 1 : risk));
%     artPop = sum(c90_2vFull.popVec(: , artInds) , 2);
%     hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%         1 : intervens , g , 4 : 10 , 1 : risk));
%     allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates, ...
%         1 : intervens , g , 4 : 10 , 1 : risk)); 
%     hivPop = sum(c90_2vFull.popVec(: , hivInds) , 2);
%     allPop = sum(c90_2vFull.popVec(: , allInds) , 2);
%     plot(tVec , 100 * (hivPop + artPop) ./ allPop)
%     hold on
% end
% legend('Male' , 'Female' , 'Male Vax' , 'Female Vax')

%%
% figure()    
% for g = 1 : 2
%     hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , g , 4 : 10 , 1 : risk)); ...
%         toInd(allcomb(7 : 9 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , g , 4 : 10 , 1 : risk))];
%     hivSus = annlz(sum(c90_2vFull.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;    
%     plot(tVec(1 : stepsPerYear : end) , ...
%         annlz(sum(sum(c90_2vFull.newHiv(: , g , 4 : 10 , :) ...
%         , 3) , 4)) ./ hivSus * 100)
%     hold on
% end
% 
% xlabel('Year'); ylabel('Rate Per 100'); title('HIV Incidence')
% hold on
% 
% for g = 1 : 2
%     hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , g , 4 : 10 , 1 : risk)); ...
%         toInd(allcomb(7 : 9 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         1 : intervens , g , 4 : 10 , 1 : risk))];
%     hivSusNo = annlz(sum(noV.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
%     plot(tVec(1 : stepsPerYear : end) , ...
%         annlz(sum(sum(noV.newHiv(: , g , 4 : 10 , :) ...
%         , 3) , 4)) ./ hivSusNo * 100 )
% end
% legend('Male' , 'Female' , 'Male No Vax' , 'Female No vax')

%% HPV incidence
inds = {':' , [2 : 6] , [1,7:9] , 10};
files = {'CC_General_Hpv_VaxCover' , ...
     'CC_HivNoART_Hpv_VaxCover' , 'CC_HivNeg_Hpv_VaxCover' ,...
     'CC_ART_HPV_VaxCover'};
plotTits = {'General' , 'HIV-Positive (No ART)' , ....
     'HIV-Negative' , 'HIV-Positive on ART'};
fac = 10 ^ 5;

figure();

for i = 2 : length(inds)
%     figure();
    noV.hpvIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
    for n = 1 : length(vaxResult)-1
        vaxResult{n}.hpvIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
    end
    
    % General, all ages
    allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , ...
        1 : intervens , 2 , 3 : age , 1 : risk)); ...
        toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 10 , ...
        1 : intervens , 2 , 3 : age , 1 : risk))];
    allhivNegFAge = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 1 : 4 , 1 : intervens , ...
            2 , 3 : age , 1 : risk)); ...
            toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 9 : 10 , 1 : intervens , ...
            2 , 3 : age , 1 : risk))];
    
    for a = 3 : age
        % General
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , ...
            1 : intervens , 2 , a , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 10 , ...
            1 : intervens , 2 , a , 1 : risk))];
        % HIV-positive women not on ART
        hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 , ...
            1 : intervens , 2 , a , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 10 , ...
            1 : intervens , 2 , a , 1 : risk))];
        % All HIV-negative women
        hivNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 , 1 , 1 : intervens , ...
            2 , a, 1 : risk)); ...
            toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvVaxStates , 10 , 1 : intervens , ...
            2 , a , 1 : risk))];
        % Women on ART
        artF = [toInd(allcomb(10 , 6 , 1 , 1 , ...
            1 : intervens , 2 , a , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 10 , ...
            1 : intervens , 2 , a , 1 : risk))];
        genArray = {allF , hivNoARTF , hivNeg , artF};

        hpvIncRef = ...
            ((annlz(sum(sum(noV.newHpv(: , 2 , inds{i} , a , :),3),5))+annlz(sum(sum(noV.newImmHpv(: , 2 , inds{i} , a , :),3),5))+annlz(sum(sum(noV.newVaxHpv(: , 2 , inds{i} , a , :),3),5))) ./ ...
            (annlz(sum(noV.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac) ...
            .* (annlz(sum(noV.popVec(: , genArray{3}) , 2) ./ stepsPerYear));
        noV.hpvIncRef = noV.hpvIncRef + hpvIncRef;
                
%         for n = 1 : length(vaxResult)-1
%             hpvIncRef = ...
%                 ((annlz(sum(sum(vaxResult{n}.newHpv(: , 2 , inds{i} , a , :),3),5))+annlz(sum(sum(vaxResult{n}.newImmHpv(: , 2 , inds{i} , a , :),3),5))+annlz(sum(sum(vaxResult{n}.newVaxHpv(: , 2 , inds{i} , a , :),3),5))) ./ ...
%                 (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac) ...
%                 .* (annlz(sum(vaxResult{n}.popVec(: , genArray{3}) , 2) ./ stepsPerYear));
%             vaxResult{n}.hpvIncRef = vaxResult{n}.hpvIncRef + hpvIncRef;
%         end
        
    end
    noV.hpvInc = noV.hpvIncRef ./ (annlz(sum(noV.popVec(: , allhivNegFAge) , 2) ./ stepsPerYear));
    plot(tVec(1 : stepsPerYear : end) , noV.hpvInc ,'DisplayName' , ...
         [plotTits{i} , ': Efficacy: ' , num2str(round(noV.vaxEff * 100)) '% ,', ...
         'Coverage: ' , num2str(round(noV.vaxRate * 100)) , '%']);
    legend('-DynamicLegend');
    hold all;
%     for n = 1 : length(vaxResult)-1
%         vaxResult{n}.hpvInc = vaxResult{n}.hpvIncRef ./ (annlz(sum(vaxResult{n}.popVec(: , allhivNegFAge) , 2) ./ stepsPerYear));
%         plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.hpvInc , 'DisplayName' , ...
%             [plotTits{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
%     end
    title(' HPV Incidence')
    xlabel('Year'); ylabel('Incidence per 100,000')
    hold all;
end       

%% HPV incidence by age, HIV+ on ART
inds = {':' , [2 : 6] , 1 , 10};
fac = 10 ^ 5;

figure();
    
for a = 3 : 10
    for r = 1 : risk
    % Women on ART
    artF = [toInd(allcomb(10 , 6 , 1 , 1 , ...
        1 : intervens , 2 , a , r)); ...
        toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 10 , ...
        1 : intervens , 2 , a , r))];

    noV.hpvInc = ...
        ((annlz(sum(noV.newHpv(: , 2 , inds{4} , a , r),3))+annlz(sum(noV.newImmHpv(: , 2 , inds{4} , a , r),3))+annlz(sum(noV.newVaxHpv(: , 2 , inds{4} , a , r),3))) ./ ...
        (annlz(sum(noV.popVec(: , artF) , 2) ./ stepsPerYear))* fac);
    subplot(1,3,r)
    plot(tVec(1 : stepsPerYear : end) , noV.hpvInc ,'DisplayName' , ...
     num2str(a));
    ylim([0 12*10^4])
    title(' HPV Incidence: HIV+ on ART')
    xlabel('Year'); ylabel('Incidence per 100,000')    
    hold all;
    end
end
legend('-DynamicLegend');

%% HPV incidence by age, HIV+
inds = {':' , [2 : 6] , 1 , 10};
fac = 10 ^ 5;

figure();
    
for a = 3 : 10
    for r = 1 : risk
    % HIV-positive women not on ART
    hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 , ...
        1 : intervens , 2 , a , r)); ...
        toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 10 , ...
        1 : intervens , 2 , a , r))];

    noV.hpvInc = ...
        ((annlz(sum(noV.newHpv(: , 2 , inds{2} , a , r),3))+annlz(sum(noV.newImmHpv(: , 2 , inds{2} , a , r),3))+annlz(sum(noV.newVaxHpv(: , 2 , inds{2} , a , r),3)))./ ...
        (annlz(sum(noV.popVec(: , hivNoARTF) , 2) ./ stepsPerYear))* fac);

    subplot(1,3,r)
    plot(tVec(1 : stepsPerYear : end) , noV.hpvInc ,'DisplayName' , ...
     num2str(a));
    ylim([0 12*10^4])
    title(' HPV Incidence: HIV+')
    xlabel('Year'); ylabel('Incidence per 100,000')
    hold all;
    end
end   
legend('-DynamicLegend');

%% HPV incidence by age, HIV-
inds = {':' , [2 : 6] , 1 , 10};
fac = 10 ^ 5;

figure();
    
for a = 3 : 10
    for r = 1 : risk
    % HIV-negative women
    hivNeg = [toInd(allcomb(1 , 1 : viral , 1 , 1 , 1 : intervens , ...
            2 , a, r)); ...
            toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 10 , 1 : intervens , ...
            2 , a , r))];

    noV.hpvInc = ...
        (annlz(sum(noV.newHpv(: , 2 , inds{3} , a , r),3)) ./ ...
        (annlz(sum(noV.popVec(: , hivNeg) , 2) ./ stepsPerYear))* fac);

    subplot(1,3,r)
    plot(tVec(1 : stepsPerYear : end) , noV.hpvInc ,'DisplayName' , ...
     num2str(a));
    title(' HPV Incidence: HIV-')
    xlabel('Year'); ylabel('Incidence per 100,000')
    ylim([0 5*10^4])
    hold all;
    end
end
legend('-DynamicLegend');

%% Population Size by age: vax vs. non-vax
figure()
for a = 1: age  
    subplot(4,4,a)
    % HIV-positive women not on ART
    hivNoARTnoV = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , [1:7,10] , ...
        [1,6] , 2 , a , 1 : risk));
    hivNoARTV = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        [2,4] , 2 , a , 1 : risk)); ...
        toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 , ...
        [1,6] , 2 , a , 1 : risk))];
    % All HIV-negative women
    hivNegnoV = toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , [1:7,10] , [1,6] , ...
        2 , a , 1 : risk));
    hivNegV = [toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , [2,4] , ...
        2 , a , 1 : risk)); ...
        toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 9 , [1,6] , ...
        2 , a , 1 : risk))];
    % Women on ART
    artnoV = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , [1:7,10] , ...
        [1,6] , 2 , a , 1 : risk));
    artV = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        [2,4] , 2 , a , 1 : risk)); ...
        toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 , ...
        [1,6] , 2 , a , 1 : risk))];
    genArraynoV = {hivNoARTnoV , hivNegnoV , artnoV};
    genArrayV = {hivNoARTV , hivNegV , artV};

    for i = 1 : length(genArraynoV)
        plot(tVec , sum(vaxResult{2}.popVec(: , genArraynoV{i}) , 2),'-')
        hold all;
    end
    set(gca,'ColorOrderIndex',1)
    for i = 1 : length(genArrayV)
        plot(tVec , sum(vaxResult{2}.popVec(: , genArrayV{i}) , 2),'--')
        hold all;
    end
    title('Population Size')
    xlabel('Year'); ylabel('Individuals')
    xlim([1910 2200]);
end
legend('HIV+ , no ART: no vax' , 'HIV-: no vax' , 'HIV+ , ART: no vax');

%% Population Size by age
figure()
for a = 1: age  
    subplot(4,4,a)
    % HIV-positive women not on ART
    hivNoART = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : intervens , 2 , a , 1 : risk));
    % All HIV-negative women
    hivNeg = toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
        2 , a , 1 : risk));
    % Women on ART
    art = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : intervens , 2 , a , 1 : risk));
    genArray = {hivNoART , hivNeg , art};

    for i = 1 : length(genArray)
        plot(tVec , sum(vaxResult{2}.popVec(: , genArray{i}) , 2),'-')
        hold all;
    end
    title('Population Size')
    xlabel('Year'); ylabel('Individuals')
    xlim([1910 2200]);
end
legend('HIV+ , no ART' , 'HIV-' , 'HIV+ , ART');

%% Population age distribution by vaccination status
% figure()
% linStyle = {'o' , '*'};
% linColor = {'[0, 0.4470, 0.7410]' , '[0.8500, 0.3250, 0.0980]' , '[0.9290, 0.6940, 0.1250]'};
% 
% for a = 3: age  
%     % HIV-positive women not on ART
%     hivNoARTnoV = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , [1:7,10] , ...
%         1 , 2 , a , 1 : risk));
%     hivNoARTV = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         2 , 2 , a , 1 : risk)); ...
%         toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 , ...
%         1 , 2 , a , 1 : risk))];
%     hivNoARTnoVTot = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , [1:7,10] , ...
%         1 , 2 , 3 : age , 1 : risk));
%     hivNoARTVTot = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         2 , 2 , 3 : age , 1 : risk)); ...
%         toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvVaxStates , 9 , ...
%         1 , 2 , 3 : age , 1 : risk))];
%     % All HIV-negative women
%     hivNegnoV = toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , [1:7,10] , 1 , ...
%         2 , a , 1 : risk));
%     hivNegV = [toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 2 , ...
%         2 , a , 1 : risk)); ...
%         toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 9 , 1 , ...
%         2 , a , 1 : risk))];
%     hivNegnoVTot = toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , [1:7,10] , 1 , ...
%         2 , 3 : age , 1 : risk));
%     hivNegVTot = [toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 2 , ...
%         2 , 3 : age , 1 : risk)); ...
%         toInd(allcomb(1 , 1 : viral , 1 : hpvVaxStates , 9 , 1 , ...
%         2 , 3 : age , 1 : risk))];
%     % Women on ART
%     artnoV = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , [1:7,10] , ...
%         1 , 2 , a , 1 : risk));
%     artV = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         2 , 2 , a , 1 : risk)); ...
%         toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 , ...
%         1 , 2 , a , 1 : risk))];
%     artnoVTot = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , [1:7,10] , ...
%         1 , 2 , 3 : age , 1 : risk));
%     artVTot = [toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%         2 , 2 , 3 : age , 1 : risk)); ...
%         toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 9 , ...
%         1 , 2 , 3 : age , 1 : risk))];
%     genArraynoV = {hivNoARTnoV , hivNegnoV , artnoV};
%     genArrayV = {hivNoARTV , hivNegV , artV};
%     genArraynoVTot = {hivNoARTnoVTot , hivNegnoVTot , artnoVTot};
%     genArrayVTot = {hivNoARTVTot , hivNegVTot , artVTot};
% 
%     for i = 1 : length(genArraynoV)
%         scatter(a , sum(vaxResult{2}.popVec(end-40*stepsPerYear , genArraynoV{i}) , 2)/sum(vaxResult{2}.popVec(end-70*stepsPerYear , genArraynoVTot{i}) , 2),'MarkerEdgeColor',linColor{i},'MarkerFaceColor','k')
%         hold all;
%     end
%     for i = 1 : length(genArrayV)
%         scatter(a , sum(vaxResult{2}.popVec(end-40*stepsPerYear , genArrayV{i}) , 2)/sum(vaxResult{2}.popVec(end-70*stepsPerYear , genArrayVTot{i}) , 2),'MarkerEdgeColor',linColor{i})
%         hold all;
%     end
% 
%     title('Population Size')
%     xlabel('Age'); ylabel('Age Proportion')
% end
% legend('HIV+ , no ART' , 'HIV-' , 'HIV+ , ART');

%% Screened by risk and HIV group
% figure();
% % for r = 1 : risk
%     for v = 1 : 2
%         for i = 1 : (length(vaxResult{1}.tVec) - length(curr.tVec))
%         %HIV-
%         rpopHivNegTot(i,1) = sumall(vaxResult{1}.newScreen(i , [1,7:9] , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : risk , v ));
% 
%         %HIV+
%         rpopHivTot(i,1) = sumall(vaxResult{1}.newScreen(i , 2 : 6 , 1 : 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : risk , v));
% 
%         %ART
%         rpopArtTot(i,1) = sumall(vaxResult{1}.newScreen(i , 10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : risk , v));
%         end
% 
%         subplot(1,3,1)
%         plot([1:(length(vaxResult{1}.tVec)-length(curr.tVec))]' , rpopHivNegTot)
%         hold all
%         subplot(1,3,2)
%         plot([1:(length(vaxResult{1}.tVec)-length(curr.tVec))]' , rpopHivTot)
%         hold all
%         subplot(1,3,3)
%         plot([1:(length(vaxResult{1}.tVec)-length(curr.tVec))]' , rpopArtTot)
%         xlabel('Year'); ylabel('Count'); title('Number')
%         hold all;
%     end
% % end
% legend('1','2','3','4','5','6','7','8','9','10')

%% Calculate life years saved
% yrIntStart = 2018;
% for n = 1 : length(vaxResult)
%     vaxResult{n}.ly = zeros((length(tVec) - length(curr.tVec)) , 1);
%     vaxResult{n}.daly = zeros((length(tVec) - length(curr.tVec)) , 1);
% end
% noV.ly = zeros((length(tVec) - length(curr.tVec)) , 1);
% noV.daly = zeros((length(tVec) - length(curr.tVec)) , 1);

%% CC Costs
% ccCost = [2617 , 8533 , 8570]; % local, regional, distant
% ccDalyWeight = 1 - [0.288 , 0.288 , 0.288]; % corresponds to local, regional, distant CC
% 
% for i = 1 : (length(tVec) - length(curr.tVec))
%     % If y = current year, count benefits and CC treatment costs for women aged
%         % >= y - B, where B = last year eligible for inclusion
%         % Since 5 year age groups are being used, at each year y, count benefits
%         % for women in age group (round((y-B)/5)) and above.
%         a = min(max(round((tVec(i + length(curr.tVec) - 1) - yrIntStart) / 5) , 1) , age);
%         ageCounted = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
%             1 : intervens , 2 , a : age , 1 : risk));
%     for n = 1 : length(vaxResult)    
%         % CC Indices
%         localInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvVaxStates , 5 , 1 : intervens , ...
%             1 : gender , a : age , 1 : risk));
%         regionalInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvVaxStates , 6 , 1 : intervens , ...
%             1 : gender , a : age , 1 : risk));
%         distantInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvVaxStates , 7 , 1 : intervens , ...
%             1 : gender , a : age , 1 : risk));
%         
%         % Count life years
%         vaxResult{n}.ly(i) = sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1  , ageCounted) , 2);
%         
%         % Count DALYs
%         % Adjust life years for CC by region according to disability
%         % weights
%         % calculate CC DALYs for each time step
%         cc_dalys = ...
%             sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , localInds) , 2) .* ccDalyWeight(1) ...
%             + sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , regionalInds) , 2) .* ccDalyWeight(2)...
%             + sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , distantInds) , 2) .* ccDalyWeight(3); 
%         
%         % DALYs obtained by subtracting full life years corresponding to CC
%         % then adding DALYs corresponding to CC (done since LY already
%         % calculated above)
%         vaxResult{n}.daly(i) = vaxResult{n}.ly(i) ...
%             - (sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , localInds) , 2) ... % subtract full LY corresponding to CC 
%             + sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , regionalInds) , 2) ...
%             + sum(vaxResult{n}.popVec(i + length(curr.tVec) - 1 , distantInds) , 2)) ...
%             + cc_dalys; % Add CC DALYs
% 
%         % Cervical cancer costs
%         vaxResult{n}.ccCosts(i) = ...
%             sum(sum(sum(vaxResult{n}.ccTreated(i , : , : , a : age , 1) , 2) , 3) , 4) .* ccCost(1) + ...
%             sum(sum(sum(vaxResult{n}.ccTreated(i , : , : , a : age , 2) , 2) , 3) , 4) .* ccCost(2) + ...
%             sum(sum(sum(vaxResult{n}.ccTreated(i , : , : , a : age , 3) , 2) , 3) , 4) .* ccCost(3);
%     end
%     
%     % no vaccine scenario
%     % Count life years
%     noV.ly(i) = sum(noV.popVec(i + length(curr.tVec) - 1 , ageCounted) , 2);
%     
%     cc_dalys = ...
%         sum(noV.popVec(i + length(curr.tVec) - 1 , localInds) , 2) .* ccDalyWeight(1) ...
%         + sum(noV.popVec(i + length(curr.tVec) - 1 , regionalInds) , 2) .* ccDalyWeight(2)...
%         + sum(noV.popVec(i + length(curr.tVec) - 1 , distantInds) , 2) .* ccDalyWeight(3);
%     
%     % DALYs obtained by subtracting full life years corresponding to CC
%     % then adding DALYs corresponding to CC (done since LY already
%     % calculated above)
%     noV.daly(i) = noV.ly(i) ...
%         - (sum(noV.popVec(i + length(curr.tVec) - 1 , localInds) , 2) ... % subtract full LY corresponding to CC
%         + sum(noV.popVec(i + length(curr.tVec) - 1 , regionalInds) , 2) ...
%         + sum(noV.popVec(i + length(curr.tVec) - 1 , distantInds) , 2)) ...
%         + cc_dalys; % Add CC DALYs
%     
%     % Cervical cancer costs
%     noV.ccCosts(i) = sum(sum(sum(noV.ccTreated(i, : , : , a : age , 1) , 2) , 3) , 4) .* ccCost(1) + ...
%         sum(sum(sum(noV.ccTreated(i , : , : , a : age , 2) , 2) , 3) , 4) .* ccCost(2) + ...
%         sum(sum(sum(noV.ccTreated(i , : , : , a : age , 3) , 2) , 3) , 4) .* ccCost(3);
% end

%%
% for n = 1 : length(vaxResult)
%     vaxResult{n}.lys = vaxResult{n}.ly - noV.ly;
% end
% 
% figure()
% for n = 1 : length(vaxResult)
%     plot(tVec(length(curr.tVec) + 1 : end) , vaxResult{n}.lys , ...
%         'DisplayName' , ['Vaccine Efficacy: ' , ...
%         num2str(round(vaxResult{n}.vaxEff * 100)) ,'%, ' , ...
%         'Vaccine Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) ,'%'])
%     hold on
%     legend('-DynamicLegend' , 'Location' , 'NorthwestOutside')
% end
% 
% figure()
% for n = 1 : length(vaxResult)
%     plot(tVec(length(curr.tVec) + 1 : end) , sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end , :) , 2)...
%         -sum(noV.popVec(length(curr.tVec) + 1 : end , :),2), ...
%         'DisplayName' , ['Vaccine Efficacy: ' , ...
%         num2str(round(vaxResult{n}.vaxEff * 100)) ,'%, ' , ...
%         'Vaccine Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) ,'%'])
%     hold on
%     legend('-DynamicLegend' , 'Location' , 'NorthwestOutside')
% end
% 
% for n = 1 : length(vaxResult)
%     vaxResult{n}.dalyPlus = vaxResult{n}.daly - noV.daly;
% end

%% Calculate annual costs
% % HIV Costs (Not used in ICER calculation. Included for completeness)
% hospCost = [117 , 56 , 38 , 38]; % <200 , 200-350, >350 , on ART
% artCost = 260; 
% 
% % CC Costs (Incurred once per person at time of cervical cancer detection)
% ccCost = [2617 , 8533 ,8570]; % local, regional, distant
% 
% % HIV Indices (Not used)
% % above350Inds = toInd(allcomb(4 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
% %     1 : gender , 1 : age , 1 : risk));
% % cd200_350Inds = toInd(allcomb(5 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
% %     1 : gender , 1 : age , 1 : risk));
% % under200Inds = toInd(allcomb(6 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
% %     1 : gender , 1 : age , 1 : risk));
% % artInds = toInd(allcomb(10 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : intervens , ...
% %     1 : gender , 1 : age , 1 : risk));

%% Vaccination 
% discountRate = 0.03; % discount rate of 3% per annum
% cost2v = 27; % cost of 2 doses of bivalent vaccine
% import java.util.LinkedList
% sims2v = LinkedList();
% full9v = LinkedList();
% full2v = LinkedList();
% maxRate9vSim = -1;
% maxRate2vSim = -1;
% maxRate9v = -1;
% maxRate2v = -1;
% for n = 1 : length(vaxResult)
%     if vaxResult{n}.vaxEff <= 0.7 && vaxResult{n}.vaxEff > 0 % Vaccine efficacy <= 70% corresponds to 2v
%         sims2v.add(n); % save index of 2v scenarios
%         vaxResult{n}.vax2vCost = annlz(vaxResult{n}.vaxd) * cost2v;
%         vaxResult{n}.vax2vCostNPV = pvvar(vaxResult{n}.vax2vCost , discountRate); % NPV of vaccination cost
%     end
%     if vaxResult{n}.vaxEff == 0.85
%         full9v.add(n);
%         if vaxResult{n}.vaxRate > maxRate9v
%             maxRate9vSim = n;
%             maxRate9v = vaxResult{n}.vaxRate;
%         end
%     elseif vaxResult{n}.vaxEff == 0.65
%         full2v.add(n);
%         if vaxResult{n}.vaxRate > maxRate2v
%             maxRate2vSim = n;
%             maxRate2v = vaxResult{n}.vaxRate; 
%         end
%     end
% end

%% Find price threshold for 9v (CC costs only)
% % 3 thresholds: 0.5x GDP , 1x GDP , 500 USD per LYS (BMGF)
% ceThreshold = 1540; % USD per LYS
% ceThresholds = [0.5 * ceThreshold , ceThreshold , 500];
% 
% % Using Life years
% % High coverage scenario (9v vs 2v)
% for i = 1 : length(ceThresholds)
%     priceGuess = 100; % Enter a price guess for 9v to seed the search process
%     % ce9v is an anonymous function that finds the vaccine price that
%     % places the 9v vaccine right at the cost-effectiveness threshold
%     % specified by ceThresholds(i)
%     ce9v = @(x) abs(pvvar(annlz(vaxResult{maxRate9vSim}.vaxd) * x - annlz(vaxResult{maxRate2vSim}.vaxd) .* cost2v ... % difference in vaccine cost for 9v vs 2v 
%         + annlz(vaxResult{maxRate9vSim}.ccCosts') - annlz(vaxResult{maxRate2vSim}.ccCosts') , discountRate) ... % difference in CC cost for 9v vs 2v scenario
%         / pvvar(annAvg(vaxResult{maxRate9vSim}.lys) - annAvg(vaxResult{maxRate2vSim}.lys) , discountRate) - ceThresholds(i)); % difference in LYS for 9v vs 2v scenario
%     priceThreshold_9v = fminsearch(ce9v , priceGuess);
%     fprintf(['\n 9v vs 2v: Considering only CC costs, with a cost-effectiveness \n' , ...
%         ' threshold of ' , num2str(ceThresholds(i)) , ' USD per LYS, ' ,...
%         ' the unit cost of 9v vaccine must be less than or equal to \n ' , ...
%         num2str(round(priceThreshold_9v , 2)),' USD. \n']) 
% end
% 
% disp(' ')
% % Using DALYs
% % High coverage scenario (9v vs 2v)
% for i = 1 : length(ceThresholds)
%     priceGuess = 100; % Enter a price guess for 9v to seed the search process
%     % ce9v is an anonymous function that finds the vaccine price that
%     % places the 9v vaccine right at the cost-effectiveness threshold
%     % specified by ceThresholds(i)
%     ce9v = @(x) abs(pvvar(annlz(vaxResult{maxRate9vSim}.vaxd) * x - annlz(vaxResult{maxRate2vSim}.vaxd) .* cost2v ... % difference in vaccine cost for 9v vs 2v 
%         + annlz(vaxResult{maxRate9vSim}.ccCosts') - annlz(vaxResult{maxRate2vSim}.ccCosts') , discountRate) ... % difference in CC cost for 9v vs 2v scenario
%         / pvvar(annAvg(vaxResult{maxRate9vSim}.daly) - annAvg(vaxResult{maxRate2vSim}.daly) , discountRate) - ceThresholds(i)); % difference in DALYs for 9v vs 2v scenario
%     priceThreshold_9v = fminsearch(ce9v , priceGuess);
%     fprintf(['\n 9v vs 2v: Considering only CC costs, with a cost-effectiveness \n' , ...
%         ' threshold of ' , num2str(ceThresholds(i)) , ' USD per DALY, ' ,...
%         ' the unit cost of 9v vaccine must be less than or equal to \n ' , ...
%         num2str(round(priceThreshold_9v , 2)),' USD.\n']) 
%end
% end

%% YLS
% figure()
% plot(tVec(1 : stepsPerYear : end) , annlz(c90_9vFull.vaxd))
% title('Vaccinated with 9v'); xlabel('Year'); ylabel('Number vaccinated')

