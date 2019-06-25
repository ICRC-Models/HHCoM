function vaxCEA_screen062419(pathModifier)

waning = 0;    % turn waning on or off

%% Load parameters
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'],'stepsPerYear','circ','condUse','disease','viral',...
    'hpvTypes','hpvStates','periods','gender','age','risk','dim','k','toInd','sumall','modelYr1')

% Load results
nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
curr = load([pwd , '\HHCoM_Results\toNow_052919_WHO_noVax']); % Population up to current year

% Helper functions
sumall = @(x) sum(x(:));
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

% Time
c = fix(clock); % get time
currYear = c(1); % get the current year from time

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
    vaxResult{n}.popVec = [curr.popVec(1 : end  , :) ; vaxResult{n}.popVec];
    vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :) ; vaxResult{n}.newCC];
    vaxResult{n}.tVec = [curr.tVec(1 : end) , vaxResult{n}.tVec];
end

noVaxInd = nSims;
noV = vaxResult{noVaxInd};
tVec = noV.tVec;
tVecYr = tVec(1 : stepsPerYear : end);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)

%% CC incidence
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
fac = 10 ^ 5;
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};

subplot(1,2,1)
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
        (annlz(sum(sum(sum(noV.newCC(: , inds{i} , : , :),2),3),4)) ./ ...
        (annlz(sum(noV.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac);
    noV.ccInc = ccIncRef; 
    % Increased vaccination scenarios
    for n = 1 : length(vaxResult)-1
        ccIncRef = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , 3 : age),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
        vaxResult{n}.ccInc = ccIncRef;
    end
    
    % Plot incidence
    % Baseline
    plot(tVec(1 : stepsPerYear : end) , noV.ccInc , 'LineStyle' , linStyle{1} , 'DisplayName' , ...
         [plotTits1{i} ]); %, ': Coverage: ' , num2str(round(noV.vaxRate * 100)) , '%']);    %'LineWidth' ,0.5,
    legend('-DynamicLegend');
    hold all;
    % Increased vavccination scenarios
%     for n = 1 : length(vaxResult)-1
%         plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccInc , 'LineStyle' , linStyle{2} , 'DisplayName' , ...
%             [plotTits1{i} , ': Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);    %'LineWidth' ,0.5,      
%     end
    
    title(' Cervical Cancer Incidence')
    xlabel('Year'); ylabel('Incidence per 100,000')
    xlim([2019 2099]);
    xticks([2019 : 10 : 2099]);
    hold all;

end        

%% CC mortality
inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
fac = 10 ^ 5;
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};

subplot(1,2,2)
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
        (annlz(sum(sum(sum(noV.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
        (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , genArray{i}) , 2) ./ stepsPerYear))* fac);
    noV.ccMort = ccMortRef;
    % Increased vaccination scenarios
    for n = 1 : length(vaxResult)-1
        ccMortRef = ...
            (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end  , genArray{i}) , 2) ./ stepsPerYear)) * fac);
        vaxResult{n}.ccMort = ccMortRef;
    end
        
    % Plot mortality
    % Baseline
    plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , noV.ccMort ,'DisplayName' , ...
         [plotTits1{i}]); % , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
         %'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
    legend('-DynamicLegend');
    hold all;
    
    for n = 1 : length(vaxResult)-1
%         plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , vaxResult{n}.ccMort , 'DisplayName' , ...
%             [plotTits1{i} , ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%']);
    end
    
    title('Cervical Cancer Mortality')
    xlabel('Year'); ylabel('Mortality per 100,000')
    xlim([2019 2099]);
    xticks([2019 : 10 : 2099]);
    hold all;
end

%% CC incidence reduction

% Run this section to plot percent reduction comparison to standard baseline
resultFileNameBaseline = [pwd , '\HHCoM_Results\Vaccine' , '052919_WHO_baselineScreen', '\' , 'vaxSimResult'];
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

inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
plotTits2 = {'General:    80% Vax Coverage' , '90% coverage'  , ...
    'HIV-negative:    80% Vax Coverage' , '90% coverage' , ...
    'HIV-positive no ART:    80% Vax Coverage' , '90% coverage' , ...
    'HIV-positive ART:    80% Vax Coverage' , '90% coverage' , ...
    'HIV-positive all:    80% Vax Coverage' , '90% coverage'};
% plotTits3 = {'2x screen' , '3x screen HIV-positive' , '3x screen HIV-positive, 50% catch-up vax 15-26' , ...
%     '3x screen HIV-positive, 80% catch-up vax 15-26' , '3x screen HIV-positive, 50% catch-up vax all' , ...
%     '3x screen HIV-positive, 80% catch-up vax all'};
fac = 10 ^ 5;
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};

subplot(2,1,2) % Turn on for Phase 1
set(gca,'ColorOrderIndex',1)

for i = 1 : length(inds) %-1
    %subplot(3,2,i); % Turn on for Phase 2
    plotTits = {plotTits2{(i*2-1):(i*2)}}; % Turn on for Phase 1
    % plotTits = {plotTits3{1}}; % Turn on for Phase 2

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
        (annlz(sum(sum(sum(noV.newCC(: , inds{i} , : , :),2),3),4)) ./ ...
        (annlz(sum(noV.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac);
    noV.ccInc = ccIncRef; 
    %Increased vaccination scenarios
    for n = 1 : length(vaxResult)-1
        ccIncRef = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , 3 : age),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
        vaxResult{n}.ccInc = ccIncRef;
    end
        
    % Plot reduction in incidence
    for n = 1 : length(vaxResult)-2
        % Calculate reduction
        vaxResult{n}.ccRed = (vaxResult{n}.ccInc - noV.ccInc) ./ noV.ccInc * 100;
        plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccRed , 'LineStyle' , linStyle{1} , 'DisplayName' , plotTits{1})  %'Color' , linColor{1} ,
%             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
        grid on
        legend('-DynamicLegend')
        xlim([2019 2099]);
        ylim([-100 0]);
        xticks([2019 : 10 : 2099]);
        hold all
              
        % Save reduction results
%         fname = [pwd , '\HHCoM_Results\Vaccine' , '052919_WHO_baselineScreen', '\' , 'baseline' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
%             'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
%         sname = [plotTits1{i} , '_IncRed']
%         if exist(fname , 'file') == 2
%             M = xlsread(fname);
%             M = catpad(2 , [tVec(1 : stepsPerYear : end)' , noV.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , M);
%             xlswrite(fname , M , sname)
%         else
%             xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , sname)
%         end
    end

    %title(plotTits1{i}) % Turn on for Phase 2
    title('Percent Reduction in Incidence')
    xlabel('Year'); ylabel('Percent change')
end        

%% CC mortality reduction

% Run this section to plot percent reduction comparison to standard baseline
resultFileNameBaseline = [pwd , '\HHCoM_Results\Vaccine' , '052919_WHO_baselineScreen', '\' , 'vaxSimResult'];
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

inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
plotTits2 = {'General:    80% Vax Coverage' , '90% coverage'  , ...
    'HIV-negative:    80% Vax Coverage' , '90% coverage' , ...
    'HIV-positive no ART:    80% Vax Coverage' , '90% coverage' , ...
    'HIV-positive ART:    80% Vax Coverage' , '90% coverage' , ...
    'HIV-positive all:    80% Vax Coverage' , '90% coverage'};
% plotTits3 = {'2x screen' , '3x screen HIV-positive' , '3x screen HIV-positive, 50% catch-up vax 15-26' , ...
%     '3x screen HIV-positive, 80% catch-up vax 15-26' , '3x screen HIV-positive, 50% catch-up vax all' , ...
%     '3x screen HIV-positive, 80% catch-up vax all'};
fac = 10 ^ 5;
linStyle = {'-' , '--' , ':'};
linColor = {'k' , '[0.8500, 0.3250, 0.0980]' , '[0, 0.4470, 0.7410]' , '[0.9290, 0.6940, 0.1250]' , 'g' , 'c'};

subplot(1,2,2);
set(gca,'ColorOrderIndex',1)

for i = 1 : length(inds) %-1
    %subplot(3,2,i); % Turn on for Phase 2
    plotTits = {plotTits2{(i*2-1):(i*2)}}; % Turn on for Phase 1
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
        (annlz(sum(sum(sum(noV.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
        (annlz(sum(noV.popVec(length(curr.tVec) + 1 : end , genArray{i}) , 2) ./ stepsPerYear))* fac);
    noV.ccMort = ccMortRef;
    % Increased vaccination scenarios
    for n = 1 : length(vaxResult)-1
        ccMortRef = ...
            (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) + 1 : end  , genArray{i}) , 2) ./ stepsPerYear)) * fac);
        vaxResult{n}.ccMort = ccMortRef;
    end
        
    % Plot reduction in mortality
    for n = 1 : length(vaxResult)-2       
        % Calculate reduction
        vaxResult{n}.ccRed = (vaxResult{n}.ccMort - noV.ccMort) ./ noV.ccMort * 100;
        plot(tVec(length(curr.tVec) + 1 : stepsPerYear : end) , vaxResult{n}.ccRed , ...
            'LineStyle' , linStyle{1} , 'DisplayName' , plotTits{n})  %'Color' , linColor{3} ,
%             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
%             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
        grid on        
        legend('-DynamicLegend')
        xlim([2019 2099]);
        ylim([-100 0]);
        xticks([2019 : 10 : 2099]);
        hold all       
        
        % Save reduction results
%         fname = [pwd , '\HHCoM_Results\Vaccine' , '052919_WHO_baselineScreen', '\' , 'baseline' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
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
    %title(plotTits1{i}) % Turn on for Phase 2
    title('Percent Reduction in Mortality')
    xlabel('Year'); ylabel('Percent change')
end
