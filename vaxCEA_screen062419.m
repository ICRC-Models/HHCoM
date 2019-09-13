function vaxCEA_screen062419()

waning = 0;    % turn waning on or off

%% LOAD PARAMETERS
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'],'stepsPerYear','circ','condUse','disease','viral',...
    'hpvTypes','hpvStates','periods','gender','age','risk','dim','k','toInd','sumall')

% Helper functions
sumall = @(x) sum(x(:));
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

% Time
c = fix(clock); % get time
currYear = c(1); % get the current year from time

%% LOAD SAVED RESULTS
curr = load([pwd , '\HHCoM_Results\toNow_061319_WHOp2']); % ***SET ME***: name for historical run file
dirName_baselineScreen = '062519_WHOp1_noBaseVax_baseScreen'; % ***SET ME***: name for baseline screening file
dirName_35Screen = '062519_WHOp1_noBaseVax_35Screen'; % ***SET ME***: name for screening @ 35 file
dirName_3545Screen = '062519_WHOp1_noBaseVax_3545Screen'; % ***SET ME***: name for screening @ 35,45 file
dirName_reductBaseline = '062519_WHOp1_noBaseVax_baseScreen'; % ***SET ME***: name for file used as baseline in percent reduction figs
simVec = {dirName_baselineScreen , dirName_35Screen , dirName_3545Screen , dirName_reductBaseline};
% % dirName_baselineScreen = '061319_WHOp2_baselineScreening'; % For Phase 2
% % dirName_SCE1 = '061319_WHOp2_SCE1';
% % dirName_SCE2 = '061319_WHOp2_SCE2';
% % dirName_SCE3 = '061319_WHOp2_SCE3';
% % dirName_SCE4 = '061319_WHOp2_SCE4';
% % dirName_SCE5 = '061319_WHOp2_SCE5';
% % dirName_reductBaseline = '061319_WHOp2_baselineScreening';
% % simVec = {dirName_baselineScreen , dirName_SCE1 , dirName_SCE2 , dirName_SCE3 , dirName_SCE4 , dirName_SCE5 , dirName_reductBaseline};

nResults = length(simVec);

allResult = cell(nResults , 1);
for i = 1 : nResults
    % Load results
    pathModifier = simVec{i};
    nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
    
    % ID correct file naming scheme for waning or no waning
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
        vaxResult{n}.popVec = [curr.popVec(1 : end  , :) ; vaxResult{n}.popVec(2 : end , :)];
        vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :) ; vaxResult{n}.newCC(2 : end , : , : , :)];
        vaxResult{n}.tVec = [curr.tVec(1 : end) , vaxResult{n}.tVec(2 : end)];
    end
    allResult{i} = vaxResult;
end
nSims = length(allResult{1,1});
noVaxInd = nSims;
tVec = allResult{nResults,1}{noVaxInd}.tVec;
tVecYr = tVec(1 : stepsPerYear : end);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)

%% CC INCIDENCE BY HIV STATUS
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
fac = 10 ^ 5;
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};

set(gca,'ColorOrderIndex',1)

for i = 1 : length(inds)-1
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

    % Calculate incidence
    % Baseline
    ccIncRef = ...
        (annlz(sum(sum(sum(allResult{1,1}{noVaxInd}.newCC(: , inds{i} , : , :),2),3),4)) ./ ...
        (annlz(sum(allResult{1,1}{noVaxInd}.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac);
    allResult{1,1}{noVaxInd}.ccInc = ccIncRef; 
    % Increased vaccination scenarios
    for n = 1 : nSims-1
        ccIncRef = ...
            (annlz(sum(sum(sum(allResult{1,1}{n}.newCC(: , inds{i} , : , 3 : age),2),3),4)) ./ ...
            (annlz(sum(allResult{1,1}{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
        allResult{1,1}{n}.ccInc = ccIncRef;
    end
    
    % Plot incidence
    % Baseline
    plot(tVec(1 : stepsPerYear : end) , allResult{1,1}{noVaxInd}.ccInc , 'LineStyle' , linStyle{1} , 'DisplayName' , ...
         [plotTits1{i} ]); %, ': Coverage: ' , num2str(round(allResult{1,1}{n}.vaxRate * 100)) , '%']);    %'LineWidth' ,0.5,
    legend('-DynamicLegend');
    hold all;
    % Increased vavccination scenarios
%     for n = 1 : nSims-1
%         plot(tVec(1 : stepsPerYear : end) , allResult{1,1}{n}.ccInc , 'LineStyle' , linStyle{2} , 'DisplayName' , ...
%             [plotTits1{i} , ': Coverage: ' , num2str(round(allResult{1,1}{n}.vaxRate * 100)) , '%']);    %'LineWidth' ,0.5,      
%     end
    
    %title(' Cervical Cancer Incidence')
    xlabel('Year'); ylabel('Incidence per 100,000')
    xlim([2019 2099]);
    xticks([2019 : 10 : 2099]);
    set(gca,'FontSize',18)
    hold all;

end        

%% CC MORTALITY BY HIV STATUS
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
fac = 10 ^ 5;
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};

set(gca,'ColorOrderIndex',1)

for i = 1 : length(inds)-1
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

    % Calculate mortality
    % Baseline
    ccMortRef = ...
        (annlz(sum(sum(sum(allResult{1,1}{noVaxInd}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
        (annlz(sum(allResult{1,1}{noVaxInd}.popVec(length(curr.tVec) : end , genArray{i}) , 2) ./ stepsPerYear))* fac);
    allResult{1,1}{noVaxInd}.ccMort = ccMortRef;
    % Increased vaccination scenarios
    for n = 1 : nSims-1
        ccMortRef = ...
            (annlz(sum(sum(sum(allResult{1,1}{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(allResult{1,1}{n}.popVec(length(curr.tVec) : end  , genArray{i}) , 2) ./ stepsPerYear)) * fac);
        allResult{1,1}{n}.ccMort = ccMortRef;
    end
        
    % Plot mortality
    % Baseline
    plot(tVec(length(curr.tVec) : stepsPerYear : end) , allResult{1,1}{noVaxInd}.ccMort ,'DisplayName' , ...
         [plotTits1{i}]); % , ': Efficacy: ' , num2str(round(allResult{1,1}{n}.vaxEff * 100)) '% ,', ...
         %'Coverage: ' , num2str(round(allResult{1,1}{n}.vaxRate * 100)) , '%']);
    legend('-DynamicLegend');
    hold all;
    
    for n = 1 : nSims-1
%         plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , allResult{1,1}{n}.ccMort , 'DisplayName' , ...
%             [plotTits1{i} , ': Efficacy: ' , num2str(round(allResult{1,1}{n}.vaxEff * 100)) '% ,', ...
%             'Coverage: ' , num2str(round(allResult{1,1}{n}.vaxRate * 100)) , '%']);
    end
    
    %title('Cervical Cancer Mortality')
    xlabel('Year'); ylabel('Mortality per 100,000')
    xlim([2019 2099]);
    ylim([0 220])
    xticks([2019 : 10 : 2099]);
    yticks([0:20:220])
    set(gca,'FontSize',18)
    hold all;
end

%% CC INCIDENCE REDUCTION- WITH VACCINATION
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
plotTits2 = {'80% coverage: Total female population' , '90% coverage'  , ...
    '80% coverage: HIV-negative' , '90% coverage' , ...
    '80% coverage: HIV-positive no ART' , '90% coverage' , ...
    '80% coverage: HIV-positive ART' , '90% coverage' , ...
    '80% coverage: HIV-positive all' , '90% coverage'};
fileTits = {'baselineScreen' , '35Screen' , '3545Screen'};
fac = 10 ^ 5;
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};
set(gca,'ColorOrderIndex',1)
    
for j = 1 : 1
    for i = 1 : length(inds)-1
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

        % Calculate incidence
        % Baseline
        ccIncRef = ...
            (annlz(sum(sum(sum(allResult{nResults,1}{noVaxInd}.newCC(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(allResult{nResults,1}{noVaxInd}.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac);
        allResult{nResults,1}{noVaxInd}.ccInc = ccIncRef; 
        %Increased vaccination scenarios
        for n = 1 : nSims-1
            ccIncRef = ...
                (annlz(sum(sum(sum(allResult{j,1}{n}.newCC(: , inds{i} , : , 3 : age),2),3),4)) ./ ...
                (annlz(sum(allResult{j,1}{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
            allResult{j,1}{n}.ccInc = ccIncRef;
        end

        % Plot reduction in incidence
        for n = 1 : nSims-1
            % Calculate reduction
            allResult{j,1}{n}.ccRed = (allResult{j,1}{n}.ccInc - allResult{nResults,1}{noVaxInd}.ccInc) ./ allResult{nResults,1}{noVaxInd}.ccInc * 100;
            plot(tVec(1 : stepsPerYear : end) , allResult{j,1}{n}.ccRed , 'LineStyle' , linStyle{n} , 'DisplayName' , plotTits2{(i*2-1)+(n-1)} , 'Color' , linColor{i})
    %             ': Efficacy ' , num2str(round(allResult{j,1}{n}.vaxEff * 100)) '% ,', ...
    %             'Coverage ' , num2str(round(allResult{j,1}{n}.vaxRate * 100)) , '%'])
            grid on
            legend('-DynamicLegend')
            xlim([2019 2099]);
            ylim([-100 0]);
            xticks([2019 : 10 : 2099]);
            hold all

            % Save reduction results
%             fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_baselineScreen, '\' , fileTits{j} , 'Efficacy' , num2str(round(allResult{j,1}{n}.vaxEff * 100)) , ...
%                 'Coverage' , num2str(round(allResult{j,1}{n}.vaxRate * 100)) , '.xlsx'];
%             sname = [plotTits1{i} , '_IncRed']
%             if exist(fname , 'file') == 2
%                 M = xlsread(fname);
%                 M = catpad(2 , [tVec(1 : stepsPerYear : end)' , allResult{nResults,1}{noVaxInd}.ccInc' , allResult{j,1}{n}.ccInc' , allResult{j,1}{n}.ccRed'] , M);
%                 xlswrite(fname , M , sname)
%             else
%                 xlswrite(fname , [tVec(1 : stepsPerYear : end)' , allResult{nResults,1}{noVaxInd}.ccInc' , allResult{j,1}{n}.ccInc' , allResult{j,1}{n}.ccRed'] , sname)
%             end

        end
    end     
    %title('Percent Reduction in Incidence')
    xlabel('Year'); ylabel('Percent change')
    set(gca,'FontSize',18)
end

%% CC INCIDENCE REDUCTION- WITH VACCINATION & SCREENING
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
plotTits2 = {'Total female population: 80% coverage, Baseline screening' , ...
    'HIV-negative: 80% coverage, Baseline screening' , ...
    'HIV-positive no ART' , ...
    'HIV-positive ART' , ...
    'HIV-positive all' , ...
    'Total female population: 80% coverage, Increased screening: Age 35'  , ...
    'HIV-negative: 80% coverage, Increased screening: Age 35'  , ...
    'HIV-positive no ART'  , ...
    'HIV-positive ART'  , ...
    'HIV-positive all'  , ...
    'Total female population: 80% coverage, Increased screening: Ages 35,45' , ...
    'HIV-negative: 80% coverage, Increased screening: Ages 35,45' , ...
    'HIV-positive no ART' , ...
    'HIV-positive ART' , ...
    'HIV-positive all'};
plotTits3 = {'3x screening for HIV-positive' , '3x screening for HIV-positive, 50% catch-up vax ages 15-26' , ...
    '3x screening for HIV-positive, 80% catch-up vax ages 15-26' , '3x screening for HIV-positive, 50% catch-up vax all ages' , ...
    '3x screening for HIV-positive, 80% catch-up vax all ages'};
fileTits = {'baselineScreen' , '35Screen' , '3545Screen'};
fac = 10 ^ 5;
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};
set(gca,'ColorOrderIndex',1)
    
for j = 1 : nResults-1
    for i = 2 : length(inds)-1
        %subplot(3,2,i); % Turn on for Phase 2
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

        % Calculate incidence
        % Baseline
        ccIncRef = ...
            (annlz(sum(sum(sum(allResult{nResults,1}{noVaxInd}.newCC(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(allResult{nResults,1}{noVaxInd}.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac);
        allResult{nResults,1}{noVaxInd}.ccInc = ccIncRef; 
        %Increased vaccination scenarios
        for n = 1 : nSims-1
            ccIncRef = ...
                (annlz(sum(sum(sum(allResult{j,1}{n}.newCC(: , inds{i} , : , 3 : age),2),3),4)) ./ ...
                (annlz(sum(allResult{j,1}{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
            allResult{j,1}{n}.ccInc = ccIncRef;
        end

        % Plot reduction in incidence
        for n = 1 : nSims-2 % For 80% vaccination
        %for n = 2 : nSims-1 % For 90% vaccination
            % Calculate reduction
            allResult{j,1}{n}.ccRed = (allResult{j,1}{n}.ccInc - allResult{nResults,1}{noVaxInd}.ccInc) ./ allResult{nResults,1}{noVaxInd}.ccInc * 100;
            plot(tVec(1 : stepsPerYear : end) , allResult{j,1}{n}.ccRed , 'LineStyle' , linStyle{j} , 'DisplayName' , plotTits2{(i-1)+(j*5-4)} , 'Color' , linColor{i})
    %             ': Efficacy ' , num2str(round(allResult{j,1}{n}.vaxEff * 100)) '% ,', ...
    %             'Coverage ' , num2str(round(allResult{j,1}{n}.vaxRate * 100)) , '%'])
            grid on
            legend('-DynamicLegend')
            xlim([2019 2099]);
            ylim([-100 0]);
            xticks([2019 : 10 : 2099]);
            hold all

            % Save reduction results
%             fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_baselineScreen, '\' , fileTits{j} , '_Screen_' , 'Efficacy' , num2str(round(allResult{j,1}{n}.vaxEff * 100)) , ...
%                 'Coverage' , num2str(round(allResult{j,1}{n}.vaxRate * 100)) , '.xlsx'];
%             sname = [plotTits1{i} , '_IncRed']
%             if exist(fname , 'file') == 2
%                 M = xlsread(fname);
%                 M = catpad(2 , [tVec(1 : stepsPerYear : end)' , allResult{nResults,1}{noVaxInd}.ccInc' , allResult{j,1}{n}.ccInc' , allResult{j,1}{n}.ccRed'] , M);
%                 xlswrite(fname , M , sname)
%             else
%                 xlswrite(fname , [tVec(1 : stepsPerYear : end)' , allResult{nResults,1}{noVaxInd}.ccInc' , allResult{j,1}{n}.ccInc' , allResult{j,1}{n}.ccRed'] , sname)
%             end

        end
    end     
    %title(plotTits1{i}) % Turn on for Phase 2
    %title('Percent Reduction in Incidence')
    xlabel('Year'); ylabel('Percent change')
    set(gca,'FontSize',18)
end

%% CC MORTALITY REDUCTION- WITH VACCINATION
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
plotTits2 = {'80% coverage: Total female population' , '90% coverage'  , ...
    '80% coverage: HIV-negative' , '90% coverage' , ...
    '80% coverage: HIV-positive no ART' , '90% coverage' , ...
    '80% coverage: HIV-positive ART' , '90% coverage' , ...
    '80% coverage: HIV-positive all' , '90% coverage'};
fileTits = {'baselineScreen' , '35Screen' , '3545Screen'};
fac = 10 ^ 5;
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};
set(gca,'ColorOrderIndex',1)

for j = 1 : 1
    for i = 1 : length(inds)-1
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

        % Calculate mortality
        % Baseline
        ccMortRef = ...
            (annlz(sum(sum(sum(allResult{nResults,1}{noVaxInd}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(allResult{nResults,1}{noVaxInd}.popVec(length(curr.tVec) : end , genArray{i}) , 2) ./ stepsPerYear))* fac);
        allResult{nResults,1}{noVaxInd}.ccMort = ccMortRef;
        % Increased vaccination scenarios
        for n = 1 : length(vaxResult)-1
            ccMortRef = ...
                (annlz(sum(sum(sum(allResult{j,1}{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
                (annlz(sum(allResult{j,1}{n}.popVec(length(curr.tVec) : end  , genArray{i}) , 2) ./ stepsPerYear)) * fac);
            allResult{j,1}{n}.ccMort = ccMortRef;
        end

        % Plot reduction in mortality
        for n = 1 : nSims-1       
            % Calculate reduction
            allResult{j,1}{n}.ccRed = (allResult{j,1}{n}.ccMort - allResult{nResults,1}{noVaxInd}.ccMort) ./ allResult{nResults,1}{noVaxInd}.ccMort * 100;
            plot(tVec(length(curr.tVec) : stepsPerYear : end) , allResult{j,1}{n}.ccRed , ...
                'LineStyle' , linStyle{n} , 'DisplayName' , plotTits2{(i*2-1)+(n-1)} , 'Color' , linColor{i})
    %             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
    %             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
            grid on        
            legend('-DynamicLegend')
            xlim([2019 2099]);
            ylim([-100 0]);
            xticks([2019 : 10 : 2099]);
            hold all       

            % Save reduction results
%             fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_baselineScreen, '\' , fileTits{j} , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
%                 'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_Mort' , '.xlsx'];
%             sname = [plotTits1{i} , '_MortRed'];
%             if exist(fname , 'file') == 2
%                 M = xlsread(fname);
%                 M = catpad(2 , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , allResult{nResults,1}{noVaxInd}.ccMort' , allResult{j,1}{n}.ccMort' , allResult{j,1}{n}.ccRed'] , M);
%                 xlswrite(fname , M , sname)
%             else
%                 xlswrite(fname , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , allResult{nResults,1}{noVaxInd}.ccMort' , allResult{j,1}{n}.ccMort' , allResult{j,1}{n}.ccRed'] , sname)
%             end

        end
    end
    %title('Percent Reduction in Mortality')
    xlabel('Year'); ylabel('Percent change')
    set(gca,'FontSize',18)
end

%% CC MORTALITY REDUCTION- WITH VACCINATION & SCREENING
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
plotTits2 = {'Total female population: 80% coverage, Baseline screening' , ...
    'HIV-negative: 80% coverage, Baseline screening' , ...
    'HIV-positive no ART' , ...
    'HIV-positive ART' , ...
    'HIV-positive all' , ...
    'Total female population: 80% coverage, Increased screening: Age 35'  , ...
    'HIV-negative: 80% coverage, Increased screening: Age 35'  , ...
    'HIV-positive no ART'  , ...
    'HIV-positive ART'  , ...
    'HIV-positive all'  , ...
    'Total female population: 80% coverage, Increased screening: Ages 35,45' , ...
    'HIV-negative: 80% coverage, Increased screening: Ages 35,45' , ...
    'HIV-positive no ART' , ...
    'HIV-positive ART' , ...
    'HIV-positive all'};
% plotTits3 = {'2x screen' , '3x screen HIV-positive' , '3x screen HIV-positive, 50% catch-up vax 15-26' , ...
%     '3x screen HIV-positive, 80% catch-up vax 15-26' , '3x screen HIV-positive, 50% catch-up vax all' , ...
%     '3x screen HIV-positive, 80% catch-up vax all'};
fileTits = {'baselineScreen' , '35Screen' , '3545Screen'};
fac = 10 ^ 5;
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};
set(gca,'ColorOrderIndex',1)

for j = 1 : nResults-1
    for i = 2 : length(inds)-1
        %subplot(3,2,i); % Turn on for Phase 2
        % plotTits = {plotTits3{1}}; % Turn on for Phase 2

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

        % Calculate mortality
        % Baseline
        ccMortRef = ...
            (annlz(sum(sum(sum(allResult{nResults,1}{noVaxInd}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(allResult{nResults,1}{noVaxInd}.popVec(length(curr.tVec) : end , genArray{i}) , 2) ./ stepsPerYear))* fac);
        allResult{nResults,1}{noVaxInd}.ccMort = ccMortRef;
        % Increased vaccination scenarios
        for n = 1 : nSims-1
            ccMortRef = ...
                (annlz(sum(sum(sum(allResult{j,1}{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
                (annlz(sum(allResult{j,1}{n}.popVec(length(curr.tVec) : end  , genArray{i}) , 2) ./ stepsPerYear)) * fac);
            allResult{j,1}{n}.ccMort = ccMortRef;
        end

        % Plot reduction in mortality
        for n = 1 : nSims-2 % For 80% vaccination
        %for n = 2 : nSims-1 % For 90% vaccination 
            % Calculate reduction
            allResult{j,1}{n}.ccRed = (allResult{j,1}{n}.ccMort - allResult{nResults,1}{noVaxInd}.ccMort) ./ allResult{nResults,1}{noVaxInd}.ccMort * 100;
            plot(tVec(length(curr.tVec) : stepsPerYear : end) , allResult{j,1}{n}.ccRed , ...
                'LineStyle' , linStyle{j} , 'DisplayName' , plotTits2{(i-1)+(j*5-4)} , 'Color' , linColor{i})
    %             ': Efficacy ' , num2str(round(allResult{j,1}{n}.vaxEff * 100)) '% ,', ...
    %             'Coverage ' , num2str(round(allResult{j,1}{n}.vaxRate * 100)) , '%'])
            grid on        
            legend('-DynamicLegend')
            xlim([2019 2099]);
            ylim([-100 0]);
            xticks([2019 : 10 : 2099]);
            hold all       

            % Save reduction results
%             fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_baselineScreen, '\' , fileTits{j} , '_Screen_' , 'Efficacy' , num2str(round(allResult{j,1}{n}.vaxEff * 100)) , ...
%                 'Coverage' , num2str(round(allResult{j,1}{n}.vaxRate * 100)) , '_Mort' , '.xlsx'];
%             sname = [plotTits1{i} , '_MortRed'];
%             if exist(fname , 'file') == 2
%                 M = xlsread(fname);
%                 M = catpad(2 , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , allResult{nResults,1}{noVaxInd}.ccMort' , allResult{j,1}{n}.ccMort' , allResult{j,1}{n}.ccRed'] , M);
%                 xlswrite(fname , M , sname)
%             else
%                 xlswrite(fname , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , allResult{nResults,1}{noVaxInd}.ccMort' , allResult{j,1}{n}.ccMort' , allResult{j,1}{n}.ccRed'] , sname)
%             end

        end
    end
    %title(plotTits1{i}) % Turn on for Phase 2
    %title('Percent Reduction in Mortality')
    xlabel('Year'); ylabel('Percent change')
    set(gca,'FontSize',18)
end
