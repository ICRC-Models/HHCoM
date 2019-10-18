function vaxCEA_whoReport_101819()

waning = 0;    % turn waning on or off

%% LOAD PARAMETERS
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'],'stepsPerYear','disease','viral',...
    'hpvTypes','hpvStates','periods','gender','age','risk','dim','k','toInd','sumall')

% Helper functions
sumall = @(x) sum(x(:));
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

% Time
c = fix(clock); % get time
currYear = 2020; %c(1); % get the current year from time

%% LOAD SAVED RESULTS
curr = load([pwd , '\HHCoM_Results\toNow_090519_singleAge_baseScreen_noBaseVax_2020']); % ***SET ME***: name for historical run file
% ***SET ME***: save names of potential scenarios to analyze as variables
dirName_reductBaseline = '090519_WHOP1_SCES12';
dirName_P1_SCE3 = '090519_WHOP1_SCES34';
dirName_P1_SCE5 = '090519_WHOP1_SCES56';
dirName_P2_SCE1 = '090519_WHOP2_SCE1';
dirName_P2_SCE2 = '090519_WHOP2_SCE2';
dirName_P2_SCE3 = '090519_WHOP2_SCE3';
dirName_P2_SCE4 = '090519_WHOP2_SCE4';
dirName_P2_SCE5 = '090519_WHOP2_SCE5';

% ***SET ME***: choose which scenarios you want to save data in Excel for
simVec = {dirName_reductBaseline , dirName_P1_SCE3 , dirName_P1_SCE5 , dirName_P2_SCE1 , dirName_P2_SCE2 , dirName_P2_SCE3 , dirName_P2_SCE4 , dirName_P2_SCE5};
% ***SET ME***: make sure the names here correspond to scenarios in simVec above
fileTits = {'P1_SCES12' , 'P1_SCES34' , 'P1_SCES56' , 'P2_SCE1' , 'P2_SCE2' , 'P2_SCE3' , 'P2_SCE4' , 'P2_SCE5'};

nResults = length(simVec);

%% SAVE INCIDENCE AND MORTALITY
for j = 1 : nResults 
    % Load results
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
        vaxResult{n}.popVec = [curr.popVec(end , :) ; vaxResult{n}.popVec(2 : end , :)];
        vaxResult{n}.newCC = [curr.newCC(end , : , : , :) ; vaxResult{n}.newCC(2 : end , : , : , :)];
        vaxResult{n}.ccDeath = [curr.ccDeath(end , : , : ,:) ; vaxResult{n}.ccDeath(2 : end , : , : ,:)];
        vaxResult{n}.tVec = [curr.tVec(end) , vaxResult{n}.tVec(2 : end)];
    end
    
    noVaxInd = nSims;
    tVec = vaxResult{noVaxInd}.tVec;
    tVecYr = tVec(1 : stepsPerYear : end);
    
    %% CC INCIDENCE - not age standardized
%     inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
%     plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
%     fac = 10 ^ 5;
%     for i = 1 : length(inds)
%         % General
%         allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk)); ...
%             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk))];
%         % All HIV-negative women
%         hivNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
%             2 , 10 : 75 , 1 : risk)); ...
%             toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
%             2 , 10 : 75 , 1 : risk))];
%         % HIV-positive women not on ART
%         hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk)); ...
%             toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk))];
%         % Women on ART
%         artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk)); ...
%             toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk))];
%         % All HIV-positive women
%         hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk)); ...
%             toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk))];
%         genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};
% 
%         % Calculate incidence
%         %Increased vaccination scenarios
%         for n = 1 : nSims
%             ccIncRef = ...
%                 (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , 10 : 75),2),3),4)) ./ ...
%                 (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
%             vaxResult{n}.ccInc = ccIncRef;
%         end
% 
%         % Save incidence results
%         for n = 1 : nSims
%             fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_reductBaseline , '\' , fileTits{j} , ...
%                 'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_RawInc_ages9-74' , '.xlsx'];
%             sname = plotTits1{i};
%             if exist(fname , 'file') == 2
%                 M = xlsread(fname);
%                 M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , M);
%                 xlswrite(fname , M , sname)
%             else
%                 xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , sname)
%             end
%         end
% %         hold all;
% %         plot(tVec(1 : stepsPerYear : end)' , vaxResult{noVaxInd}.ccInc')
% %         hold all;
% %         axis([1980 2120 0 275])
% %         xlabel('Year'); ylabel('Incidence rates (per 100K)'); 
% %         legend('General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all');
% 
% %         for n = 1 : 2
% %             if (n == 1)
% %                 s = 2;
% %             end
% %             if (n == 2)
% %                 s = 1;
% %             end
% %             if ((j == 1) && (i == 1)) || ((i == 1) && (s == 1))
% %                 hold all;
% %                 plot(tVec(1 : stepsPerYear : end)' , vaxResult{s}.ccInc')
% %                
% %                 hold all;
% %                 axis([2020 2120 0 70])
% %                 xlabel('Year'); ylabel('Cervical cancer incidence rate (per 100K)'); 
% %                 legend('Scenario 0' , 'Scenario 1' , 'Scenario 2' , 'Scenario 3' , 'Scenario 4' , 'Scenario 5');
% %                 %legend('General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all');
% %             end
% %         end
%         
%     end     
    
    %% CC INCIDENCE - age standardized
    inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
    plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
    fac = 10 ^ 5;
    
    %worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
    %    247167 240167 226750 201603 171975 150562 113118 82266 64484];
    worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
        247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
        23477 9261 2155];
    
    for i = 1 : length(inds)
        for n = 1 : nSims
            vaxResult{n}.ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
        end
        
        for aInd = 1:age+4
            if aInd >= age
                a = age;
            end
            % General
            allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            % All HIV-negative women
            hivNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
                2 , a , 1 : risk)); ...
                toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
                2 , a , 1 : risk))];
            % HIV-positive women not on ART
            hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            % Women on ART
            artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            % All HIV-positive women
            hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};

            % Calculate incidence
            for n = 1 : nSims
                if aInd <= age    
                    ccIncRef = ...
                        (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , a),2),3),4)) ./ ...
                        (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac) ...
                        .* (worldStandard_WP2015(aInd));
                    if (i == 4) && (a < 3) && (max(annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , a),2),3),4))) == 0.0)
                        ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
                    end
                elseif aInd > age
                    ccIncRef = ...
                        (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , a),2),3),4)) ./ ...
                        (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
                    ccIncRef = [(ones(1,aInd-a).*ccIncRef(1,1)) , ccIncRef(1,1:end-(aInd-a))];
                    ccIncRef = ccIncRef .* worldStandard_WP2015(aInd);
                end
                vaxResult{n}.ccIncRef = vaxResult{n}.ccIncRef + ccIncRef;
            end
        end
        for n = 1 : nSims
            vaxResult{n}.ccInc = vaxResult{n}.ccIncRef ./ (sum(worldStandard_WP2015(1:age+4)));
        end
        
        % Save incidence results
        for n = 1 : nSims
            fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_reductBaseline , '\' , fileTits{j} , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_AgeStandInc_ages0-99' , '.xlsx'];
            sname = plotTits1{i};
            if exist(fname , 'file') == 2
                M = xlsread(fname);
                M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , M);
                xlswrite(fname , M , sname)
            else
                xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , sname)
            end
        end
%         hold all;
%         plot(tVec(1 : stepsPerYear : end)' , vaxResult{noVaxInd}.ccInc' , '--')
%         hold all;
%         axis([1980 2120 0 150])
%         xlabel('Year'); ylabel('Incidence rates (per 100K)'); 
%         legend('General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all' , 'General-ARToff' , 'HIV-negative-ARToff' , 'HIV-positive no ART-ARToff' , 'HIV-positive ART-ARToff' , 'HIV all-ARToff');
        
    end  

    %% CC MORTALITY - not age standardized
%     inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
%     plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
%     fac = 10 ^ 5;
%     for i = 1 : length(inds)
%         % General
%         allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 7 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk)); ...
%             toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk))];
%         % All HIV-negative women
%         hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 7 , 1 : periods , ...
%             2 , 10 : 75 , 1 : risk)); ...
%             toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
%             2 , 10 : 75 , 1 : risk))];
%         % HIV-positive women not on ART
%         hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 7 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk)); ...
%             toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk))];
%         % Women on ART
%         artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk)); ...
%             toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk))];
%         % All HIV-positive women
%         hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 7 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk)); ...
%             toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%             1 : periods , 2 , 10 : 75 , 1 : risk))];
%         genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};
% 
%         % Calculate mortality
%         % Increased vaccination scenarios
%         for n = 1 : nSims
%             ccMortRef = ...
%                 (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , 10 : 75),2),3),4)) ./ ...
%                 (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
%             vaxResult{n}.ccMort = ccMortRef;
%         end
% 
%         % Save mortality results
%         for n = 1 : nSims
%             fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_reductBaseline , '\' , fileTits{j} , ...
%                 'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_RawMort_ages9-74' , '.xlsx'];
%             sname = plotTits1{i};
%             if exist(fname , 'file') == 2
%                 M = xlsread(fname);
%                 M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccMort'] , M);
%                 xlswrite(fname , M , sname)
%             else
%                 xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccMort'] , sname)
%             end
%         end
%         
%     end

%% CC MORTALITY - age standardized
    inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
    plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
    fac = 10 ^ 5;
    
    %worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
    %    247167 240167 226750 201603 171975 150562 113118 82266 64484];
    worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
        247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
        23477 9261 2155];
    
    for i = 1 : length(inds)
        for n = 1 : nSims
            vaxResult{n}.ccMortRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
        end
        
        for aInd = 1:age+4
            if aInd >= age
                a = age;
            end
            % General
            allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 7 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            % All HIV-negative women
            hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 7 , 1 : periods , ...
                2 , a , 1 : risk)); ...
                toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
                2 , a , 1 : risk))];
            % HIV-positive women not on ART
            hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 7 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            % Women on ART
            artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 7 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            % All HIV-positive women
            hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 7 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};
            
            % Calculate mortality
            for n = 1 : nSims
                if aInd <= age    
                    ccMortRef = ...
                        (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , a),2),3),4)) ./ ...
                        (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac) ...
                        .* (worldStandard_WP2015(aInd));
                    if (i == 4) && (a < 3) && (max(annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , a),2),3),4))) == 0.0)
                        ccMortRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
                    end
                elseif aInd > age
                    ccMortRef = ...
                        (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , a),2),3),4)) ./ ...
                        (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
                    ccMortRef = [(ones(1,aInd-a).*ccMortRef(1,1)) , ccMortRef(1,1:end-(aInd-a))];
                    ccMortRef = ccMortRef .* worldStandard_WP2015(aInd);
                end
                vaxResult{n}.ccMortRef = vaxResult{n}.ccMortRef + ccMortRef;
            end
        end
        
        for n = 1 : nSims
            vaxResult{n}.ccMort = vaxResult{n}.ccMortRef ./ (sum(worldStandard_WP2015(1:age+4)));
        end

        % Save mortality results
        for n = 1 : nSims
            fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_reductBaseline , '\' , fileTits{j} , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_AgeStandMort_ages0-99' , '.xlsx'];
            sname = plotTits1{i};
            if exist(fname , 'file') == 2
                M = xlsread(fname);
                M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccMort'] , M);
                xlswrite(fname , M , sname)
            else
                xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccMort'] , sname)
            end
        end
            
    end

    clear vaxResult; 
end
