function vaxCEA_calibReduct090319()

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
curr = load([pwd , '\HHCoM_Results\toNow_090319calib_22Aug19_baseline']); % ***SET ME***: name for historical run file
dirName_calibSim = '090319calib_22Aug19_baseline'; % ***SET ME***: name for baseline screening file
simVec = {dirName_calibSim};


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
    
noVaxInd = nSims;
tVec = vaxResult{noVaxInd}.tVec;
tVecYr = tVec(1 : stepsPerYear : end);

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)


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
    
for i = 1 : length(inds)
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
        (annlz(sum(sum(sum(vaxResult{noVaxInd}.newCC(: , inds{i} , : , :),2),3),4)) ./ ...
        (annlz(sum(vaxResult{noVaxInd}.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac);
    vaxResult{noVaxInd}.ccInc = ccIncRef; 
    %Increased vaccination scenarios
    for n = 1 : nSims-1
        ccIncRef = ...
            (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , 3 : age),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
        vaxResult{n}.ccInc = ccIncRef;
    end

    % Plot reduction in incidence
    for n = 1 : nSims-1
        % Calculate reduction
        allResult{j,1}{n}.ccRed = (vaxResult{n}.ccInc - vaxResult{noVaxInd}.ccInc) ./ vaxResult{noVaxInd}.ccInc * 100;
%             plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccRed , 'LineStyle' , linStyle{n} , 'DisplayName' , plotTits2{(i*2-1)+(n-1)} , 'Color' , linColor{i})
%     %             ': Efficacy ' , num2str(round(allResult{j,1}{n}.vaxEff * 100)) '% ,', ...
%     %             'Coverage ' , num2str(round(allResult{j,1}{n}.vaxRate * 100)) , '%'])
%             grid on
%             legend('-DynamicLegend')
%             xlim([2019 2099]);
%             ylim([-100 0]);
%             xticks([2019 : 10 : 2099]);
%             hold all

        % Save reduction results
        fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_calibSim, '\' , fileTits{j} , 'Efficacy' , num2str(round(allResult{j,1}{n}.vaxEff * 100)) , ...
            'Coverage' , num2str(round(allResult{j,1}{n}.vaxRate * 100)) , '.xlsx'];
        sname = [plotTits1{i} , '_IncRed']
        if exist(fname , 'file') == 2
            M = xlsread(fname);
            M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{noVaxInd}.ccInc' , allResult{j,1}{n}.ccInc' , vaxResult{n}.ccRed'] , M);
            xlswrite(fname , M , sname)
        else
            xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{noVaxInd}.ccInc' , allResult{j,1}{n}.ccInc' , vaxResult{n}.ccRed'] , sname)
        end

    end
end     
%     %title('Percent Reduction in Incidence')
%     xlabel('Year'); ylabel('Percent change')
%     set(gca,'FontSize',18)


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

for i = 1 : length(inds)
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
        (annlz(sum(sum(sum(vaxResult{noVaxInd}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
        (annlz(sum(vaxResult{noVaxInd}.popVec(length(curr.tVec) : end , genArray{i}) , 2) ./ stepsPerYear))* fac);
    vaxResult{noVaxInd}.ccMort = ccMortRef;
    % Increased vaccination scenarios
    for n = 1 : length(vaxResult)-1
        ccMortRef = ...
            (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
            (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) : end  , genArray{i}) , 2) ./ stepsPerYear)) * fac);
        vaxResult{n}.ccMort = ccMortRef;
    end

    % Plot reduction in mortality
    for n = 1 : nSims-1       
        % Calculate reduction
        vaxResult{n}.ccRed = (vaxResult{n}.ccMort - vaxResult{noVaxInd}.ccMort) ./ vaxResult{noVaxInd}.ccMort * 100;
%         plot(tVec(length(curr.tVec) : stepsPerYear : end) , allResult{j,1}{n}.ccRed , ...
%             'LineStyle' , linStyle{n} , 'DisplayName' , plotTits2{(i*2-1)+(n-1)} , 'Color' , linColor{i})
% %             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
% %             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
%         grid on        
%         legend('-DynamicLegend')
%         xlim([2019 2099]);
%         ylim([-100 0]);
%         xticks([2019 : 10 : 2099]);
%         hold all       

        % Save reduction results
            fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_calibSim, '\' , fileTits{j} , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_Mort' , '.xlsx'];
            sname = [plotTits1{i} , '_MortRed'];
            if exist(fname , 'file') == 2
                M = xlsread(fname);
                M = catpad(2 , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , vaxResult{noVaxInd}.ccMort' , vaxResult{n}.ccMort' , vaxResult{n}.ccRed'] , M);
                xlswrite(fname , M , sname)
            else
                xlswrite(fname , [tVec(length(curr.tVec) + 1 : stepsPerYear : end)' , vaxResult{noVaxInd}.ccMort' , vaxResult{n}.ccMort' , vaxResult{n}.ccRed'] , sname)
            end

    end
end
% %title('Percent Reduction in Mortality')
% xlabel('Year'); ylabel('Percent change')
% set(gca,'FontSize',18)
