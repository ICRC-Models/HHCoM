function vaxCEA_screen041819(pathModifier)

waning = 0;    % turn waning on or off

%% Load parameters
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'],'stepsPerYear','circ','condUse','disease','viral',...
    'hpvTypes','hpvStates','periods','gender','age','risk','dim','k','toInd','sumall','modelYr1')

sumall = @(x) sum(x(:));

% Load results
nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
curr = load([pwd , '\HHCoM_Results\toNow_052219']); % Population up to current year

% Helper functions
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

% Time
c = fix(clock); % get time
currYear = c(1); % get the current year from time

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
    vaxResult{n}.popVec = [curr.popVec(1 : end  , :) ; vaxResult{n}.popVec];
    vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :) ; vaxResult{n}.newCC];
    vaxResult{n}.tVec = [curr.tVec(1 : end) , vaxResult{n}.tVec];
end

% Find no vaccine scenario
% noVaxInd = -1;
% for n = 1 : nSims
%     if vaxResult{n}.vaxEff == 0
%         noVaxInd = n;
%     end
% end
noVaxInd = nSims;

noV = vaxResult{noVaxInd};
tVec = noV.tVec;
tVecYr = tVec(1 : stepsPerYear : end);

%% Run this section to plot percent reduction comparison to standard baseline
resultFileNameBaseline = [pwd , '\HHCoM_Results\Vaccine' , '051519_CISNET_baseline_9yoVax', '\' , 'vaxSimResult'];
for n = nSims
    % load results from vaccine run into cell array
    noV = load([resultFileNameBaseline , num2str(n), '.mat']);
    % concatenate vectors/matrices of population up to current year to population
    % matrices for years past current year
    noV.popVec = [curr.popVec(1 : end  , :) ; noV.popVec];
    noV.newCC = [curr.newCC(1 : end , : , : , :) ; noV.newCC];
    noV.tVec = [curr.tVec(1 : end) , noV.tVec];
end
tVecYr = tVec(1 : stepsPerYear : end);

%% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)

%% CC incidence and incidence reduction
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]};
files = {'CC_General_Hpv_VaxCover' , ...
     'CC_HivNeg_Hpv_VaxCover' , 'CC_HivNoART_Hpv_VaxCover' ,...
     'CC_ART_HPV_VaxCover' , 'CC_HivAll_HPV_VaxCover'};
plotTits1 = {'General' , 'HIV-' , ...
    'HIV+ no ART' , 'HIV+ ART' , 'HIV all'};
plotTits2 = {'80% coverage: Total female population' , '90% coverage'  , ...
    '80% coverage: HIV-Negative' , '90% coverage' , ...
    '80% coverage: HIV-Positive no ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive all' , '90% coverage'};
fac = 10 ^ 5;
%plotTitsS = {'No Screening' , 'Screen- NH baseline' , 'Screen- WHO35' , 'Screen- WHO3545' , 'Screen- WHO45'};
plotTitsS = {'Baseline Screening' , 'Increased Screening- Age 35' , 'Increased Screening- Ages 35,45'};
%linStyle = {'-' , '--' , '-.' , ':' , '--'};
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};

%figure();
subplot(1,2,1);
set(gca,'ColorOrderIndex',1)

for i = 1 %2 : length(inds)-1
     plotTits = {plotTits2{(i*2-1):(i*2)}}; % Turn on for percent reduction plots
% %     figure();
% %     noV.ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
% %     for n = 1 : length(vaxResult)-1
% %         vaxResult{n}.ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
% %     end
%     
% %     % General, all ages
% %     allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
% %         1 : periods , 2 , 3 : age , 1 : risk)); ...
% %         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
% %         1 : periods , 2 , 3 : age , 1 : risk))];
% %     allhivNegFAge = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
% %             2 , 3 : age , 1 : risk)); ...
% %             toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
% %             2 , 3 : age , 1 : risk))];
    
% %     for a = 3 : age
        % General
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % All HIV-negative women
        hivNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , 3 : age , 1 : risk)); ...
            toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
            2 , 3 : age , 1 : risk))];
        % HIV-positive women not on ART
        hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % Women on ART
        artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % All HIV-positive women
        hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};

        ccIncRef = ...
            (annlz(sum(sum(sum(noV.newCC(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(noV.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac);
% %             .* (annlz(sum(noV.popVec(: , genArray{3}) , 2) ./ stepsPerYear));
% %         noV.ccIncRef = noV.ccIncRef + ccIncRef; 
        noV.ccIncRef = ccIncRef; 
                
        for n = 1 : length(vaxResult)-1
            ccIncRef = ...
                (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , 3 : age),2),3),4)) ./ ...
                (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
% %                 .* (annlz(sum(vaxResult{n}.popVec(: , genArray{3}) , 2) ./ stepsPerYear));
% %             vaxResult{n}.ccIncRef = vaxResult{n}.ccIncRef + ccIncRef;
            vaxResult{n}.ccIncRef = ccIncRef;
        end
        
% %     end
    
% %     noV.ccInc = noV.ccIncRef ./ (annlz(sum(noV.popVec(: , allhivNegFAge) , 2) ./ stepsPerYear));
    noV.ccInc = noV.ccIncRef;
%     plot(tVec(1 : stepsPerYear : end) , noV.ccInc , 'LineStyle' , linStyle{2} , 'DisplayName' , ...
%          [plotTits1{i} , ': Coverage: ' , num2str(round(noV.vaxRate * 100)) , '%']);    %'LineWidth' ,0.5,
%     %legend('-DynamicLegend');
%     hold all;
    for n = 1 : 1 %length(vaxResult)-1
%         vaxResult{n}.ccInc = vaxResult{n}.ccIncRef ./ (annlz(sum(vaxResult{n}.popVec(: , allhivNegFAge) , 2) ./ stepsPerYear));
        vaxResult{n}.ccInc = vaxResult{n}.ccIncRef;
%         plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccInc , 'LineStyle' , linStyle{2} , 'DisplayName' , ...
%             [plotTits1{i} , ': Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);    %'LineWidth' ,0.5,
        
        % Reduction
        vaxResult{n}.ccRed = (vaxResult{n}.ccInc - noV.ccInc) ./ noV.ccInc * 100;
        plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccRed , 'LineStyle' , linStyle{n} , 'Color' , linColor{3} , 'DisplayName' , plotTits{n})
%             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
        grid on
        legend('-DynamicLegend')
        xlim([2019 2099]);
        ylim([-100 0]);
        xticks([2019 : 10 : 2099]);
        hold all
              
        % Save reduction results
%         fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
%             'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
%         sname = [plotTits1{i} , '_IncRed'];
%         if exist(fname , 'file') == 2
%             M = xlsread(fname);
%             M = catpad(2 , [tVec(1 : stepsPerYear : end)' , noV.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , M);
%             xlswrite(fname , M , sname)
%         else
%             xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , sname)
%         end
        
    end
%     title(' Cervical Cancer Incidence')
%     xlabel('Year'); ylabel('Incidence per 100,000')
%     hold all;

    title('Figure 1: Percent reduction in cervical cancer cases, by HIV status')
    xlabel('Year'); ylabel('Percent change')
end        

%% CC mortality and mortality reduction
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]};
files = {'CC_General_Hpv_VaxCover' , ...
     'CC_HivNeg_Hpv_VaxCover' , 'CC_HivNoART_Hpv_VaxCover' ,...
     'CC_ART_HPV_VaxCover' , 'CC_HivAll_HPV_VaxCover'};
plotTits1 = {'General' , 'HIV-' , ...
    'HIV+ no ART' , 'HIV+ ART' , 'HIV all'};
plotTits2 = {'80% coverage: Total female population' , '90% coverage'  , ...
    '80% coverage: HIV-Negative' , '90% coverage' , ...
    '80% coverage: HIV-Positive no ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive ART' , '90% coverage' , ...
    '80% coverage: HIV-Positive all' , '90% coverage'};
fac = 10 ^ 5;
linStyle = {'-' , '--'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};

%figure();
subplot(1,2,2);

for i = 1 %1 : length(inds)-1
    plotTits = {plotTits2{(i*2-1):(i*2)}}; % Turn on for percent reduction plots
% %     figure();
% %     noV.ccMortRef = zeros(length(tVec(length(curr.tVec) + 1 : stepsPerYear : end)),1)';
% %     for n = 1 : length(vaxResult)-1
% %         vaxResult{n}.ccMortRef = zeros(length(tVec(length(curr.tVec) + 1 : stepsPerYear : end)),1)';
% %     end
    
    % General, all ages
% %     allFAge = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
% %         1 : periods , 2 , 3 : age , 1 : risk)); ...
% %         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
% %         1 : periods , 2 , 3 : age , 1 : risk))];
% %     allhivNegFAge = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
% %             2 , 3 : age , 1 : risk)); ...
% %             toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
% %             2 , 3 : age , 1 : risk))];
      
% %     for a = 3 : age
        % General
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % All HIV-negative women
        hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
            2 , 3 : age , 1 : risk))];
        % HIV-positive women not on ART
        hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % Women on ART
        artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        % All HIV-positive women
        hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 3 : age , 1 : risk)); ...
            toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 3 : age , 1 : risk))];
        genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};

        ccMortRef = ...
            (annlz(sum(sum(sum(noV.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , genArray{i}) , 2) ./ stepsPerYear))* fac);
% %             .* (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , genArray{1}) , 2) ./ stepsPerYear));
% %         noV.ccMortRef = noV.ccMortRef + ccMortRef;
        noV.ccMortRef = ccMortRef;
                
        for n = 1 : length(vaxResult)-1
            ccMortRef = ...
                (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
                (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end  , genArray{i}) , 2) ./ stepsPerYear)) * fac);
% %                 .* (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end  , genArray{1}) , 2) ./ stepsPerYear));
% %             vaxResult{n}.ccMortRef = vaxResult{n}.ccMortRef + ccMortRef;
            vaxResult{n}.ccMortRef = ccMortRef;
        end
        
% %     end
% %     noV.ccMort = noV.ccMortRef ./ (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , allFAge) , 2) ./ stepsPerYear));
      noV.ccMort = noV.ccMortRef;
%     plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , noV.ccMort ,'DisplayName' , ...
%          [plotTits1{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%          'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
%     legend('-DynamicLegend');
%     hold all;
    for n = 1 : 1 %length(vaxResult)-1
% %         vaxResult{n}.ccMort = vaxResult{n}.ccMortRef ./ (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end , allFAge) , 2) ./ stepsPerYear));
        vaxResult{n}.ccMort = vaxResult{n}.ccMortRef;
%         plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , vaxResult{n}.ccMort , 'DisplayName' , ...
%             [plotTits1{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
        
        % Reduction
        vaxResult{n}.ccRed = (vaxResult{n}.ccMort - noV.ccMort) ./ noV.ccMort * 100;
        plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , vaxResult{n}.ccRed , ...
            'LineStyle' , linStyle{n} , 'Color' , linColor{3} , 'DisplayName' , plotTits{n}) 
%             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
        grid on        
        legend('-DynamicLegend')
        xlim([2019 2099]);
        ylim([-100 0]);
        xticks([2019 : 10 : 2099]);
        hold all       
        
        % Save reduction results
%         fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
%             'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_Mort' , '.xlsx'];
%         sname = [plotTits1{i} , '_MortRed'];
%         if exist(fname , 'file') == 2
%             M = xlsread(fname);
%             M = catpad(2 , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , noV.ccMort' , vaxResult{n}.ccMort' , vaxResult{n}.ccRed'] , M);
%             xlswrite(fname , M , sname)
%         else
%             xlswrite(fname , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , noV.ccMort' , vaxResult{n}.ccMort' , vaxResult{n}.ccRed'] , sname)
%         end
        
    end
%     title('Cervical Cancer Mortality')
%     xlabel('Year'); ylabel('Incidence per 100,000')
%     hold all;

    title('Figure 2: Percent reduction in cervical cancer mortality, by HIV status')
    xlabel('Year'); ylabel('Percent change')
end

%% Population Size
% % HIV-positive women not on ART
% hivNoART = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 7 , ...
%     1 : periods , 1 : gender , 1 : age , 1 : risk)); ...
%     toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%     1 : periods , 1 : gender , 1 : age , 1 : risk))];
% % All HIV-negative women
% hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 7 , 1 : periods , ...
%     1 : gender , 1 : age , 1 : risk)); ...
%     toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
%     1 : gender , 1 : age , 1 : risk))];
% % Women on ART
% art = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 7 , ...
%     1 : periods , 1 : gender , 1 : age , 1 : risk)); ...
%     toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
%     1 : periods , 1 : gender , 1 : age , 1 : risk))];
% genArray = {hivNoART , hivNeg , art};
% 
% figure()
% subplot(1,2,1)
% for i = 1 : length(genArray)
%     plot(tVec , sum(noV.popVec(: , genArray{i}) , 2))
%     hold all;
% end
% subplot(1,2,2)
% plot(tVec , sum(noV.popVec(: , genArray{1}) , 2)+sum(noV.popVec(: , genArray{2}) , 2)+sum(noV.popVec(: , genArray{3}) , 2))
% title('Population Size')
% xlabel('Year'); ylabel('Individuals')
% xlim([1910 2099]);
% legend('HIV+ , no ART' , 'HIV-' , 'HIV+ , ART');
% hold off

%% CIN2+ prevalence by HIV group
%figure()
%subplot(1,2,2);

%plotTitsS = {'No Screening' , 'Screen- NH baseline' , 'Screen- WHO35' , 'Screen- WHO3545'};
plotTitsS = {'Baseline Screening' , 'Increased Screening- Age 35' , 'Increased Screening- Ages 35,45'};
%linStyle = {'-' , '--' , '-.' , ':' , '--'};
linStyle = {'-' , '--' , ':'};

% HIV+
hpvHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 2 , 3 : 4 , ...
     1 : periods , 2 , 3 : age , 1 : risk));
hpvHivPop = sum(noV.popVec(: , hpvHivInds) , 2);
popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates ,  1 : periods , ...
    2 , 3 : age , 1 : risk)));
%ART
hpvArtInds = toInd(allcomb(10 , 6 , 2  , 3 : 4 , ...
     1 : periods , 2 , 3 : age , 1 : risk));
hpvArtPop = sum(noV.popVec(: , hpvArtInds) , 2);
popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    2 , 3 : age , 1 : risk)));
%HIV-
hpvHivNegInds = toInd(allcomb([1,7:9] , 1 , 2 , 3 : 4 , ...
     1 : periods , 2 , 3 : age , 1 : risk));
hpvHivNegPop = sum(noV.popVec(: , hpvHivNegInds) , 2);
popHivNegTot = noV.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates ,  1 : periods , ...
    2 , 3 : age , 1 : risk)));

set(gca,'ColorOrderIndex',1)
plot(tVec , 100 * hpvHivNegPop ./ sum(popHivNegTot , 2),linStyle{3})    %'LineWidth' ,0.5
hold all
plot(tVec , 100 * hpvHivPop ./ sum(popHivTot , 2),linStyle{3})    %'LineWidth' ,0.5
hold all
plot(tVec , 100 * hpvArtPop ./ sum(popArtTot , 2),linStyle{3})    %'LineWidth' ,0.5
xlabel('Year'); ylabel('Prevalence (%)'); title(' CIN2+ Prevalence')
%legend('HIV-' , 'HIV+ noART' , 'HIV+ ART')

%% CC prevalence by HIV group
% % HIV+
% ccHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 2 : 4 , 5 : 7, ...
%      [1,6] , 2 , 3 : age , 1 : risk));
% ccHivPop = sum(noV.popVec(: , ccHivInds) , 2);
% popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , [1:8,10] ,  [1,6] , ...
%     2 , 3 : age , 1 : risk)));
% %ART
% ccArtInds = toInd(allcomb(10, 6 , 2 : 4 , 5 : 7, ...
%      [1,6] , 2 , 3 : age , 1 : risk));
% ccArtPop = sum(noV.popVec(: , ccArtInds) , 2);
% popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , [1:8,10] ,  [1,6] , ...
%     2 , 3 : age , 1 : risk)));
% %HIV-
% ccHivNegInds = toInd(allcomb(1 , 1 , 2 : 4 , 5 : 7, ...
%      [1,6] , 2 , 3 : age , 1 : risk));
% ccHivNegPop = sum(noV.popVec(: , ccHivNegInds) , 2);
% popHivNegTot = noV.popVec(: , toInd(allcomb(1 , 1 , 1 : hpvTypes , [1:8,10] ,  [1,6] , ...
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

%% Population by "p"
% figure();
% subplot(2,2,1);
% for p = 1 : periods
%     % General
%     inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%          p , 1 : gender , 8 , 1 : risk));
%     pop = sum(vaxResult{1}.popVec(: , inds) , 2);
%     popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%          1 : periods , 1 : gender , 8 , 1 : risk)));
%     plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
%     xlabel('Year'); ylabel('Proportion (%)'); title(' p Proportion')
%     legend('1' , '2' , '3' , '4' , '5' ,'6')
%     hold all;
% end
% 
% subplot(2,2,3);
% inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%      [1,6] , 1 : gender , 8 , 1 : risk));
% pop = sum(vaxResult{1}.popVec(: , inds) , 2);
% popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%      1 : periods , 1 : gender , 8 , 1 : risk)));
% plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
% xlabel('Year'); ylabel('Proportion (%)'); title(' Proportion unvaccinated or vaccinated and not reinfected ')
% 
% subplot(2,2,4);
% inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%      [2,4] , 1 : gender , 8 , 1 : risk));
% pop = sum(vaxResult{1}.popVec(: , inds) , 2);
% popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%      1 : periods , 1 : gender , 8 , 1 : risk)));
% plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
% xlabel('Year'); ylabel('Proportion (%)'); title(' Proportion vaccinated and reinfected ')
% 
% subplot(2,2,2);
% inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%      [4,6] , 1 : gender , 8 , 1 : risk));
% pop = sum(vaxResult{1}.popVec(: , inds) , 2);
% popTot = vaxResult{1}.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%      1 : periods , 1 : gender , 8 , 1 : risk)));
% plot(tVec , 100 * pop ./ sum(popTot , 2),'o')
% xlabel('Year'); ylabel('Proportion (%)'); title(' Proportion screened')

%% Screened proportion by HIV group
% figure();
% linStyle = {'--' , '-' , ':'};
% for a = 1 : age
%     for r = 1 : risk
%     % HIV+
%     vaxHivInds = toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 4 : 6 , 2 , a , r));
%     vaxHivPop = sum(noV.popVec(: , vaxHivInds) , 2);
%     popHivTot = noV.popVec(: , toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , a , r)));
%     %ART
%     vaxArtInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 4 : 6 , 2 , a , r));
%     vaxArtPop = sum(noV.popVec(: , vaxArtInds) , 2);
%     popArtTot = noV.popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , a , r)));
%     %HIV-
%     vaxHivNegInds = toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , 4 : 6 , 2 , a , r));
%     vaxHivNegPop = sum(noV.popVec(: , vaxHivNegInds) , 2);
%     popHivNegTot = noV.popVec(: , toInd(allcomb([1,7:9] , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , a , r)));
% 
%     subplot(4,4,a)
%     plot(tVec , 100 * vaxHivNegPop ./ sum(popHivNegTot , 2) , linStyle{r})
%     hold all
%     plot(tVec , 100 * vaxHivPop ./ sum(popHivTot , 2) , linStyle{r})
%     hold all
%     plot(tVec , 100 * vaxArtPop ./ sum(popArtTot , 2) , linStyle{r})
%     %axis([tVec(1) tVec(end) 0 100])
%     xlabel('Year'); ylabel('Proportion (%)'); title('Screened Proportion')
% %     xlim([2018 2100])
%     
%     hold all;
%     end
% end
% legend('HIV- lr' , 'HIV+ noART lr' , 'ART lr' , 'HIV- mr' , 'HIV+ noART mr' , 'ART mr' , 'HIV- hr' , 'HIV+ noART hr' , 'ART hr')
