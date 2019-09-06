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
dirName_calibCurr1 = 'toNow_090319calib_22Aug19_6_6';
dirName_calibCurr2 = 'toNow_090319calib_22Aug19_6_7';
dirName_calibCurr3 = 'toNow_090319calib_22Aug19_6_8';
dirName_calibCurr4 = 'toNow_090319calib_22Aug19_6_9';
dirName_calibCurr5 = 'toNow_090319calib_22Aug19_6_10';
dirName_calibCurr6 = 'toNow_090319calib_22Aug19_6_1';
dirName_calibCurr7 = 'toNow_090319calib_22Aug19_6_2';
dirName_calibCurr8 = 'toNow_090319calib_22Aug19_6_3';
dirName_calibCurr9 = 'toNow_090319calib_22Aug19_6_4';
dirName_calibCurr10 = 'toNow_090319calib_22Aug19_6_5';
dirName_calibCurr11 = 'toNow_090319calib_22Aug19_6_6';
dirName_calibCurr12 = 'toNow_090319calib_22Aug19_6_7';
dirName_calibCurr13 = 'toNow_090319calib_22Aug19_6_8';
dirName_calibCurr14 = 'toNow_090319calib_22Aug19_6_9';
dirName_calibCurr15 = 'toNow_090319calib_22Aug19_6_10';
dirName_calibCurr16 = 'toNow_090319calib_22Aug19_6_1';
dirName_calibCurr17 = 'toNow_090319calib_22Aug19_6_2';
dirName_calibCurr18 = 'toNow_090319calib_22Aug19_6_3';
dirName_calibCurr19 = 'toNow_090319calib_22Aug19_6_4';
dirName_calibCurr20 = 'toNow_090319calib_22Aug19_6_5';
dirName_calibCurr21 = 'toNow_090319calib_22Aug19_6_6';
dirName_calibCurr22 = 'toNow_090319calib_22Aug19_6_7';
dirName_calibCurr23 = 'toNow_090319calib_22Aug19_6_8';
dirName_calibCurr24 = 'toNow_090319calib_22Aug19_6_9';
dirName_calibCurr25 = 'toNow_090319calib_22Aug19_6_10';
dirName_calibCurr26 = 'toNow_090319calib_22Aug19_baseline';
dirName_calibCurr27 = 'toNow_090319calib_22Aug19_baseline';
dirName_calibCurr28 = 'toNow_090319calib_22Aug19_baseline';
dirName_calibCurr29 = 'toNow_090319calib_22Aug19_6_1';
dirName_calibCurr30 = 'toNow_090319calib_22Aug19_6_2';
dirName_calibCurr31 = 'toNow_090319calib_22Aug19_6_3';
dirName_calibCurr32 = 'toNow_090319calib_22Aug19_6_4';
dirName_calibCurr33 = 'toNow_090319calib_22Aug19_6_5';
currVec = {dirName_calibCurr26 , dirName_calibCurr27 , dirName_calibCurr28 , ...
    dirName_calibCurr29 , dirName_calibCurr30 , dirName_calibCurr31 , ...
    dirName_calibCurr32 , dirName_calibCurr33}; 
% dirName_calibCurr2 , ... % dirName_calibCurr1 , dirName_calibCurr3
%     dirName_calibCurr6 , dirName_calibCurr7 , ... %dirName_calibCurr4 , dirName_calibCurr5 , 
%     dirName_calibCurr8 , dirName_calibCurr9 , dirName_calibCurr10 , ...
%     dirName_calibCurr11 , dirName_calibCurr12 , dirName_calibCurr13 ,  ... %dirName_calibCurr14 ,
%     dirName_calibCurr15 , 
% dirName_calibCurr16 , dirName_calibCurr17 , dirName_calibCurr18 , ...
%     dirName_calibCurr19 , dirName_calibCurr20 , dirName_calibCurr21 ,  ... %dirName_calibCurr22 ,
%     dirName_calibCurr23
%dirName_calibCurr14 , dirName_calibCurr15 , dirName_calibCurr23 , dirName_calibCurr24 ,  dirName_calibCurr25
%dirName_calibCurr24 ,  dirName_calibCurr25

dirName_calibSim1 = '090319calib_22Aug19_6_6';
dirName_calibSim2 = '090319calib_22Aug19_6_7';
dirName_calibSim3 = '090319calib_22Aug19_6_8';
dirName_calibSim4 = '090319calib_22Aug19_6_9';
dirName_calibSim5 = '090319calib_22Aug19_6_10';
dirName_calibSim6 = '090319calib_22Aug19_6_1_CU50';
dirName_calibSim7 = '090319calib_22Aug19_6_2_CU50';
dirName_calibSim8 = '090319calib_22Aug19_6_3_CU50';
dirName_calibSim9 = '090319calib_22Aug19_6_4_CU50';
dirName_calibSim10 = '090319calib_22Aug19_6_5_CU50';
dirName_calibSim11 = '090319calib_22Aug19_6_6_CU50';
dirName_calibSim12 = '090319calib_22Aug19_6_7_CU50';
dirName_calibSim13 = '090319calib_22Aug19_6_8_CU50';
dirName_calibSim14 = '090319calib_22Aug19_6_9_CU50';
dirName_calibSim15 = '090319calib_22Aug19_6_10_CU50';
dirName_calibSim16 = '090319calib_22Aug19_6_1_CU80';
dirName_calibSim17 = '090319calib_22Aug19_6_2_CU80';
dirName_calibSim18 = '090319calib_22Aug19_6_3_CU80';
dirName_calibSim19 = '090319calib_22Aug19_6_4_CU80';
dirName_calibSim20 = '090319calib_22Aug19_6_5_CU80';
dirName_calibSim21 = '090319calib_22Aug19_6_6_CU80';
dirName_calibSim22 = '090319calib_22Aug19_6_7_CU80';
dirName_calibSim23 = '090319calib_22Aug19_6_8_CU80';
dirName_calibSim24 = '090319calib_22Aug19_6_9_CU80';
dirName_calibSim25 = '090319calib_22Aug19_6_10_CU80';
dirName_calibSim26 = '090319calib_22Aug19_baseline';
dirName_calibSim27 = '090319calib_22Aug19_baseline_CU50';
dirName_calibSim28 = '090319calib_22Aug19_baseline_CU80';
dirName_calibSim29 = '090319calib_22Aug19_6_1';
dirName_calibSim30 = '090319calib_22Aug19_6_2';
dirName_calibSim31 = '090319calib_22Aug19_6_3';
dirName_calibSim32 = '090319calib_22Aug19_6_4';
dirName_calibSim33 = '090319calib_22Aug19_6_5';
simVec = {dirName_calibSim26 , dirName_calibSim27 , dirName_calibSim28 , ...
    dirName_calibSim29 , dirName_calibSim30 , dirName_calibSim31 , dirName_calibSim32 , ...
    dirName_calibSim33}; 
% dirName_calibSim2 , ... %dirName_calibSim1 ,  dirName_calibSim3 ,
%     dirName_calibSim6 , dirName_calibSim7 , ... %dirName_calibSim4 , dirName_calibSim5 , 
%     dirName_calibSim8 , dirName_calibSim9 , dirName_calibSim10 , ...
%     dirName_calibSim11 , dirName_calibSim12 , dirName_calibSim13 ,  ... %dirName_calibSim14 ,
%     dirName_calibSim15 , 
% dirName_calibSim16 , dirName_calibSim17 , dirName_calibSim18 , ...
%     dirName_calibSim19 , dirName_calibSim20 , dirName_calibSim21 , ... %dirName_calibSim22 , 
%     dirName_calibSim23 
%dirName_calibSim14 , dirName_calibSim15 , dirName_calibSim23 ,
%dirName_calibSim24 , dirName_calibSim25 %, dirName_calibSim24 , dirName_calibSim25

r80inc = [];
r90inc = [];
r80mort = [];
r90mort = [];

for j = 1 : length(simVec)
    % Load results
    currModifier = currVec{j};
    curr = load([pwd , '\HHCoM_Results\' , currModifier]); % ***SET ME***: name for historical run file
    pathModifier = simVec{j};
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
        vaxResult{n}.popVec = [curr.popVec(1 : end  , :) ; vaxResult{n}.popVec(2 : end-5 , :)];
        vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :) ; vaxResult{n}.newCC(2 : end-5 , : , : , :)];
        vaxResult{n}.ccDeath = [vaxResult{n}.ccDeath(2 : end-5 , : , : , :)];
        vaxResult{n}.tVec = [curr.tVec(1 : end) , vaxResult{n}.tVec(2 : end-5)];
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

    for i = 1 : 1 %length(inds)
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
            vaxResult{n}.ccRed = (vaxResult{n}.ccInc - vaxResult{noVaxInd}.ccInc) ./ vaxResult{noVaxInd}.ccInc * 100;
    %             plot(tVec(1 : stepsPerYear : end) , vaxResult{n}.ccRed , 'LineStyle' , linStyle{n} , 'DisplayName' , plotTits2{(i*2-1)+(n-1)} , 'Color' , linColor{i})
    %     %             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
    %     %             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
    %             grid on
    %             legend('-DynamicLegend')
    %             xlim([2019 2099]);
    %             ylim([-100 0]);
    %             xticks([2019 : 10 : 2099]);
    %             hold all

            if n == 1
                r80inc = [r80inc , vaxResult{n}.ccRed(end)];
            elseif n == 2
                r90inc = [r90inc , vaxResult{n}.ccRed(end)];
            end
    
            % Save reduction results
            fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
            sname = [plotTits1{i} , '_IncRed'];
            if exist(fname , 'file') == 2
                M = xlsread(fname);
                M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{noVaxInd}.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , M);
                xlswrite(fname , M , sname)
            else
                xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{noVaxInd}.ccInc' , vaxResult{n}.ccInc' , vaxResult{n}.ccRed'] , sname)
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

    for i = 1 : 1 %length(inds)
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
            (annlz(sum(vaxResult{noVaxInd}.popVec(length(curr.tVec) +1 : end , genArray{i}) , 2) ./ stepsPerYear))* fac);
        vaxResult{noVaxInd}.ccMort = ccMortRef;
        % Increased vaccination scenarios
        for n = 1 : length(vaxResult)-1
            ccMortRef = ...
                (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
                (annlz(sum(vaxResult{n}.popVec(length(curr.tVec) +1 : end  , genArray{i}) , 2) ./ stepsPerYear)) * fac);
            vaxResult{n}.ccMort = ccMortRef;
        end

        % Plot reduction in mortality
        for n = 1 : nSims-1       
            % Calculate reduction
            vaxResult{n}.ccRed = (vaxResult{n}.ccMort - vaxResult{noVaxInd}.ccMort) ./ vaxResult{noVaxInd}.ccMort * 100;
    %         plot(tVec(length(curr.tVec) : stepsPerYear : end) , vaxResult{n}.ccRed , ...
    %             'LineStyle' , linStyle{n} , 'DisplayName' , plotTits2{(i*2-1)+(n-1)} , 'Color' , linColor{i})
    % %             ': Efficacy ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
    % %             'Coverage ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
    %         grid on        
    %         legend('-DynamicLegend')
    %         xlim([2019 2099]);
    %         ylim([-100 0]);
    %         xticks([2019 : 10 : 2099]);
    %         hold all       

            if n == 1
                r80mort = [r80mort , vaxResult{n}.ccRed(end)];
            elseif n == 2
                r90mort = [r90mort , vaxResult{n}.ccRed(end)];
            end
    
            % Save reduction results
                fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\' , 'Efficacy' , num2str(round(vaxResult{n}.vaxEff * 100)) , ...
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
    
    r80inc
    r90inc
    r80mort
    r90mort
end

r80inc
r90inc
r80mort
r90mort
