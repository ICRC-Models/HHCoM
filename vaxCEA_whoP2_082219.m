function vaxCEA_whoP2_082219()

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
currYear = 2019; %c(1); % get the current year from time

%% LOAD SAVED RESULTS
curr = load([pwd , '\HHCoM_Results\toNow_082119_singleAge_baseScreen_noBaseVax_2020']); % ***SET ME***: name for historical run file
dirName_reductBaseline = '082119_WHOP1_SCE1';
dirName_P1_SCE3 = '082119_WHOP1_SCE3';
dirName_P1_SCE5 = '082119_WHOP1_SCE5';
dirName_P2_SCE1 = '082119_WHOP2_SCE1';
dirName_P2_SCE2 = '082119_WHOP2_SCE2';
dirName_P2_SCE3 = '082119_WHOP2_SCE3';
dirName_P2_SCE4 = '082119_WHOP2_SCE4';
dirName_P2_SCE5 = '082119_WHOP2_SCE5';
simVec = {dirName_reductBaseline , dirName_P1_SCE3 , dirName_P1_SCE5 , dirName_P2_SCE1 , dirName_P2_SCE2 , dirName_P2_SCE3 , dirName_P2_SCE4 , dirName_P2_SCE5};

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
    inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
    plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
    fileTits = {'P1_SCE1' , 'P1_SCE3' , 'P1_SCE5' , 'P2_SCE1' , 'P2_SCE2' , 'P2_SCE3' , 'P2_SCE4' , 'P2_SCE5'};
    fac = 10 ^ 5;
    for i = 1 : length(inds)
        % General
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 1 : age , 1 : risk))];
        % All HIV-negative women
        hivNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , 1 : age , 1 : risk)); ...
            toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
            2 , 1 : age , 1 : risk))];
        % HIV-positive women not on ART
        hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 1 : age , 1 : risk))];
        % Women on ART
        artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 1 : age , 1 : risk))];
        % All HIV-positive women
        hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 1 : age , 1 : risk))];
        genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};

        % Calculate incidence
        %Increased vaccination scenarios
        for n = 1 : nSims
            ccIncRef = ...
                (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , 1 : age),2),3),4)) ./ ...
                (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
            vaxResult{n}.ccInc = ccIncRef;
        end

        % Save incidence results
        for n = 1 : nSims
            fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_reductBaseline , '\' , fileTits{j} , '_Screen_' , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_RawInc' , '.xlsx'];
            sname = plotTits1{i};
            if exist(fname , 'file') == 2
                M = xlsread(fname);
                M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , M);
                xlswrite(fname , M , sname)
            else
                xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , sname)
            end

        end
    end     
    
    
    %% CC INCIDENCE - age standardized
    inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
    plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
    fileTits = {'P1_SCE1' , 'P1_SCE3' , 'P1_SCE5' , 'P2_SCE1' , 'P2_SCE2' , 'P2_SCE3' , 'P2_SCE4' , 'P2_SCE5'};
    fac = 10 ^ 5;
    
    worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
        247167 240167 226750 201603 171975 150562 113118 82266 64484];
    
    aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
    
    for i = 1 : length(inds)
        for n = 1 : nSims
            vaxResult{n}.ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
        end
        
        for aInd = 1 : (age/5)
            a = aVec{aInd};
            
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
                ccIncRef = ...
                    (annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , a),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac) ...
                    .* (worldStandard_WP2015(aInd));
                if (i == 4) && (aInd < 3) && (max(annlz(sum(sum(sum(vaxResult{n}.newCC(: , inds{i} , : , a),2),3),4))) == 0.0)
                    ccIncRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
                end
                vaxResult{n}.ccIncRef = vaxResult{n}.ccIncRef + ccIncRef;
            end
        end
        for n = 1 : nSims
            vaxResult{n}.ccInc = vaxResult{n}.ccIncRef ./ (sum(worldStandard_WP2015(1:(age/5))));
        end
        
        % Save incidence results
        for n = 1 : nSims
            fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_reductBaseline , '\' , fileTits{j} , '_Screen_' , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_AgeStandInc' , '.xlsx'];
            sname = plotTits1{i};
            if exist(fname , 'file') == 2
                M = xlsread(fname);
                M = catpad(2 , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , M);
                xlswrite(fname , M , sname)
            else
                xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.ccInc'] , sname)
            end
        end
        
    end  

    %% CC MORTALITY - not age standardized
    inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
    plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
    fileTits = {'P1_SCE1' , 'P1_SCE3' , 'P1_SCE5' , 'P2_SCE1' , 'P2_SCE2' , 'P2_SCE3' , 'P2_SCE4' , 'P2_SCE5'};
    fac = 10 ^ 5;
    for i = 1 : length(inds)
        % General
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 1 : age , 1 : risk))];
        % All HIV-negative women
        hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
            2 , 1 : age , 1 : risk))];
        % HIV-positive women not on ART
        hivNoARTF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 1 : age , 1 : risk))];
        % Women on ART
        artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 1 : age , 1 : risk))];
        % All HIV-positive women
        hivAllF = [toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 1 : age , 1 : risk)); ...
            toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 1 : age , 1 : risk))];
        genArray = {allF , hivNeg , hivNoARTF , artF , hivAllF};

        % Calculate mortality
        % Increased vaccination scenarios
        for n = 1 : nSims
            ccMortRef = ...
                (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , :),2),3),4)) ./ ...
                (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
            vaxResult{n}.ccMort = ccMortRef;
        end

        % Save mortality results
        for n = 1 : nSims
            fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_reductBaseline , '\' , fileTits{j} , '_Screen_' , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_RawMort' , '.xlsx'];
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

%% CC MORTALITY - age standardized
    inds = {':' , [1,7:9] , [2 : 6] , 10 , [2:6,10]}; % HIV state inds
    plotTits1 = {'General' , 'HIV-negative' , 'HIV-positive no ART' , 'HIV-positive ART' , 'HIV all'};
    fileTits = {'P1_SCE1' , 'P1_SCE3' , 'P1_SCE5' , 'P2_SCE1' , 'P2_SCE2' , 'P2_SCE3' , 'P2_SCE4' , 'P2_SCE5'};
    fac = 10 ^ 5;
    
    worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
        247167 240167 226750 201603 171975 150562 113118 82266 64484];
    
    aVec = {1:5,6:10,11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,51:55,56:60,61:65,66:70,71:75,76:80};
    
    for i = 1 : length(inds)
        for n = 1 : nSims
            vaxResult{n}.ccMortRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
        end
        
        for aInd = 1 : (age/5)
            a = aVec{aInd};
    
            % General
            allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
                1 : periods , 2 , a , 1 : risk)); ...
                toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
                1 : periods , 2 , a , 1 : risk))];
            % All HIV-negative women
            hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
                2 , a , 1 : risk)); ...
                toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
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
            
            % Calculate mortality
            for n = 1 : nSims
                ccMortRef = ...
                    (annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , a),2),3),4)) ./ ...
                    (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac) ...
                    .* (worldStandard_WP2015(aInd));
                if (i == 4) && (aInd < 3) && (max(annlz(sum(sum(sum(vaxResult{n}.ccDeath(: , inds{i} , : , a),2),3),4))) == 0.0)
                    ccMortRef = zeros(length(tVec(1 : stepsPerYear : end)),1)';
                end
                vaxResult{n}.ccMortRef = vaxResult{n}.ccMortRef + ccMortRef;
            end
        end
        
        for n = 1 : nSims
            vaxResult{n}.ccMort = vaxResult{n}.ccMortRef ./ (sum(worldStandard_WP2015(1:(age/5))));
        end

        % Save mortality results
        for n = 1 : nSims
            fname = [pwd , '\HHCoM_Results\Vaccine' , dirName_reductBaseline , '\' , fileTits{j} , '_Screen_' , ...
                'Coverage' , num2str(round(vaxResult{n}.vaxRate * 100)) , '_AgeStandMort' , '.xlsx'];
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
