function vaxCEA_whoCompare_090619()

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
dirName_P1_SCE1 = '090519_WHOP1_SCES12';
dirName_P1_SCE3 = '090519_WHOP1_SCES34';
dirName_P1_SCE5 = '090519_WHOP1_SCES56';

% ***SET ME***: choose which scenarios you want to save data in Excel for
simVec = {dirName_P1_SCE1 , dirName_P1_SCE3 , dirName_P1_SCE5}; %{dirName_P2_SCE1_2xHIVneg , dirName_P2_SCE2_2xHIVneg};
% ***SET ME***: make sure the names here correspond to scenarios in simVec above
fileTits = {'P1_SCE1' , 'P1_SCE3' , 'P1_SCE5'};

nResults = length(simVec);

%% SAVE INCIDENCE
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
        vaxResult{n}.tVec = [curr.tVec(end) , vaxResult{n}.tVec(2 : end)];
    end
    
    noVaxInd = nSims;
    tVec = vaxResult{noVaxInd}.tVec;
    tVecYr = tVec(1 : stepsPerYear : end);
    
    %% CC INCIDENCE - not age standardized
    inds = {':' , [1,7:9] , [2:6,10] , 10 , [2 : 6]}; % HIV state inds
    plotTits1 = {'Pop' , 'HIV-' , 'HIV+' , 'HIV+ ART' , 'HIV+ no ART'};
    fac = 10 ^ 5;
    for i = 1 : length(inds)
        
        for n = 1 : nSims
            vaxResult{n}.incAgeVec = zeros(age , length(tVec(1 : stepsPerYear : end)));
            vaxResult{n}.denomAgeVec = zeros(age , length(tVec(1 : stepsPerYear : end)));
        end
        
        for a = 1 : age
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
            % All HIV-positive women
            hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            % Women on ART
            artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            % HIV-positive women not on ART
            hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];

            genArray = {allF , hivNeg , hivAllF , artF , hivNoARTF};

            % Calculate incidence
            %Increased vaccination scenarios
            for n = 1 : nSims
                ccIncRef = ...
                    (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , a),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
                vaxResult{n}.ccInc = ccIncRef;
            
                vaxResult{n}.incAgeVec(a , :) = vaxResult{n}.ccInc;
                vaxResult{n}.denomAgeVec(a , :) = (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear));

            end  
        end

        % Save incidence results
        for n = 1 : nSims
            fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_P1_SCE1 , '\' , fileTits{j} , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_RawInc_byAge' , '.xlsx'];
            sname = plotTits1{i};
            if exist(fname , 'file') == 2
                M = xlsread(fname);
                M = catpad(2 , [tVec(1 : stepsPerYear : end) ; vaxResult{n}.incAgeVec] , M);
                xlswrite(fname , M , sname)
            else
                xlswrite(fname , [tVec(1 : stepsPerYear : end) ; vaxResult{n}.incAgeVec] , sname)
            end

        end       
        
        % Save denominator results
        for n = 1 : nSims
            fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_P1_SCE1 , '\' , fileTits{j} , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_Denom_byAge' , '.xlsx'];
            sname = plotTits1{i};
            if exist(fname , 'file') == 2
                M = xlsread(fname);
                M = catpad(2 , [tVec(1 : stepsPerYear : end) ; vaxResult{n}.denomAgeVec] , M);
                xlswrite(fname , M , sname)
            else
                xlswrite(fname , [tVec(1 : stepsPerYear : end) ; vaxResult{n}.denomAgeVec] , sname)
            end

        end       
    end     
    
    clear vaxResult; 
end

