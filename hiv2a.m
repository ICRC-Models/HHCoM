% HIV Progression
% Simulates progression through HIV states and takes treatments such as
% ART, PrEP, etc. into account (currently no PrEP)
% Accepts:
% 1) Population matrix (pop)
% 2) Matrix describing the distribution of individuals (by disease status,
% viral load) that went on ART over the last stepsPerYear*2 time steps as input (artDist).
% Returns:
% 1) dPop, a matrix of derivatives that describes the change in the
% population's subgroups due to HIV progression.
% 2) artTreat, a matrix describing the distribution of individuals who went
% on ART according to their disease and viral load status at the time they
% went on treatment.
function[dPop , extraOuts] = hiv2a(t , pop , vlAdvancer , artDist , muHIV , ...
    kCD4 ,  maxRateM , maxRateF , minLim , maxLim , disease , ...
    viral , gender , age , risk , k , hivInds , ...
    stepsPerYear , year)

%% Load constants and parameters
toInd = @(x) (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + (x(: , 6) - 1) * k(5) ...
    + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + (x(: , 3) - 1) * k(2) ...
    + (x(: , 2) - 1) * k(1) + x(: , 1);
sumall = @(x) sum(x(:));

%% Calculate ART treatment coverage
artOut = zeros(gender , age ,risk);
artDist = reshape(artDist, [disease , viral , gender , age , risk]);
treat = zeros(disease , viral , gender , age ,risk);

% CD4 <= 200, from 2004 to 2010
if year >= 2004 && year < 2010
    if year >= 2004 && year < 2007
        yrs = 2004 : 1 / stepsPerYear : 2007;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.0 , maxRateM(1) , length(yrs)) ,...
            linspace(0.0 , maxRateF(1) , length(yrs))};
    elseif year >= 2007 && year < 2010
        yrs = 2007 : 1 / stepsPerYear : 2011;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(maxRateM(1) , maxRateM(2) , length(yrs)) ,...
            linspace(maxRateF(1) , maxRateF(2) , length(yrs))};
    end
    ageVec = [1 : age];
    dRange = [6];
    for g = 1 : gender
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVSubTots = zeros(1 , age); % number HIV-positives by age
        for a = 11 : age
            onArtAge = sumall(pop(hivInds(10 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            for v = 1 : 5
                totHivPosAge = totHivPosAge + sumall(pop(hivInds(6 , v , g , a , : , :)));
            end
            ageARTSubTots(1 , a) = onArtAge;
            ageHIVSubTots(1 , a) = totHivPosAge;    
        end
        fracARTAge = (ageARTSubTots ./ (ageARTSubTots + ageHIVSubTots)); % fraction on ART by age
        agePopSubTots = ageARTSubTots + ageHIVSubTots; % total HIV-positives (on/off ART) by age
        minCoverLim = popCover{g}(ind) * minLim; % minimum ART coverage by age
        maxCoverLim = popCover{g}(ind) * maxLim; % maximum ART coverage by age
        popCoverInd = popCover{g}(ind); % desired population-level ART coverage

        % Calculate treat/artOut matrices to maintain ART coverage min/max by age
        [artOut , treat , maxAges , excMaxAges , minAges , excMinAges] = ...
            artMinMax(artOut , treat , minCoverLim , maxCoverLim , ...
            fracARTAge , ageVec , g , risk , 11 , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , 11 , ...
            agePopSubTots , dRange);
    end
end

% CD4 <= 350, from 2010 to 2013
if year >= 2010 && year < 2013
    if year >= 2010 && year < 2011
        yrs = 2007 : 1/stepsPerYear : 2011;
        ind = round(yrs , 4) == round(year , 4);
        popCover = {linspace(maxRateM(1) , maxRateM(2) , length(yrs)) , ...
            linspace(maxRateF(1) , maxRateF(2) , length(yrs))};
    elseif year >= 2011 && year < 2013
        yrs = 2011 : 1/stepsPerYear : 2016;
        ind = round(yrs , 4) == round(year , 4);
        popCover = {linspace(maxRateM(2) , maxRateM(3) , length(yrs)) , ...
            linspace(maxRateF(2) , maxRateF(3) , length(yrs))};
    end
    ageVec = [1 : age];
    dRange = [5 : 6];
    for g = 1 : gender   
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVSubTots = zeros(1 , age); % number HIV-positives by age
        for a = 11 : age 
            onArtAge = sumall(pop(hivInds(10 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            for d = 5 : 6
                for v = 1 : 5
                    totHivPosAge = totHivPosAge + sumall(pop(hivInds(d , v , g , a , : , :)));
                end
            end
            ageARTSubTots(1 , a) = onArtAge;
            ageHIVSubTots(1 , a) = totHivPosAge;  
        end 
        fracARTAge = (ageARTSubTots ./ (ageARTSubTots + ageHIVSubTots)); % fraction on ART by age
        agePopSubTots = ageARTSubTots + ageHIVSubTots; % total HIV-positives (on/off ART) by age
        
        minCoverLim = popCover{g}(ind) * minLim; % minimum ART coverage by age
        maxCoverLim = popCover{g}(ind) * maxLim; % maximum ART coverage by age
        popCoverInd = popCover{g}(ind); % desired population-level ART coverage

        % Calculate treat/artOut matrices to maintain ART coverage min/max by age
        [artOut , treat , maxAges , excMaxAges , minAges , excMinAges] = ...
            artMinMax(artOut , treat , minCoverLim , maxCoverLim , ...
            fracARTAge , ageVec , g , risk , 11 , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , 11 , ...
            agePopSubTots , dRange);
    end
end

% CD4 <= 500, from 2013 to 2016
if year >= 2013 && year < 2016
    yrs = 2011 : 1/stepsPerYear : 2016;
    ind = round(yrs , 4) == round(year , 4);
    popCover = {linspace(maxRateM(2) , maxRateM(3) , length(yrs)) , ...
        linspace(maxRateF(2) , maxRateF(3) , length(yrs))};
    ageVec = [1 : age];
    dRange = [4 : 6];
    for g = 1 : gender    
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVSubTots = zeros(1 , age); % number HIV-positives by age
        for a = 11 : age
            onArtAge = sumall(pop(hivInds(10 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            for d = 4 : 6
                for v = 1 : 5
                    totHivPosAge = totHivPosAge + sumall(pop(hivInds(d , v , g , a , : , :)));
                end
            end
            ageARTSubTots(1 , a) = onArtAge;
            ageHIVSubTots(1 , a) = totHivPosAge;       
        end
        fracARTAge = (ageARTSubTots ./ (ageARTSubTots + ageHIVSubTots)); % fraction on ART by age
        agePopSubTots = ageARTSubTots + ageHIVSubTots; % total HIV-positives (on/off ART) by age
        
        minCoverLim = popCover{g}(ind) * minLim; % minimum ART coverage by age
        maxCoverLim = popCover{g}(ind) * maxLim; % maximum ART coverage by age
        popCoverInd = popCover{g}(ind); % desired population-level ART coverage

        % Calculate treat/artOut matrices to maintain ART coverage min/max by age
        [artOut , treat , maxAges , excMaxAges , minAges , excMinAges] = ...
            artMinMax(artOut , treat , minCoverLim , maxCoverLim , ...
            fracARTAge , ageVec , g , risk , 11 , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , 11 , ...
            agePopSubTots , dRange);
    end
end

% Any CD4, after 2016
if year >= 2016
    if year >= 2016 && year < 2018
        yrs = 2016 : 1 / stepsPerYear : 2018;
        ind = round(yrs , 4) == round(year , 4);
        popCover = {linspace(maxRateM(3) , maxRateM(4) , length(yrs)) ,...
           linspace(maxRateF(3) , maxRateF(4) , length(yrs))};
        ageVec = [1 : age];
        dRange = [2 : 6];
    elseif year >= 2018
        yrs = 2018 : 1 / stepsPerYear : 2026;
        ind = round(yrs , 4) == round(year , 4);
        popCover = {linspace(maxRateM(4) , maxRateM(5) , length(yrs)) ,...
           linspace(maxRateF(4) , maxRateF(5) , length(yrs))};
        ageVec = [1 : age];
        dRange = [2 : 6];
    end
    for g = 1 : gender
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVSubTots = zeros(1 , age); % number HIV-positives by age
        for a = 11 : age
            onArtAge = sumall(pop(hivInds(10 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            for d = 2 : 6
                for v = 1 : 5
                    totHivPosAge = totHivPosAge + sumall(pop(hivInds(d , v , g , a , : , :)));
                end
            end
            ageARTSubTots(1 , a) = onArtAge;
            ageHIVSubTots(1 , a) = totHivPosAge;       
        end
        fracARTAge = (ageARTSubTots ./ (ageARTSubTots + ageHIVSubTots)); % fraction on ART by age
        agePopSubTots = ageARTSubTots + ageHIVSubTots; % total HIV-positives (on/off ART) by age
        
        if year < 2026
            minCoverLim = popCover{g}(ind) * minLim; % minimum ART coverage by age
            maxCoverLim = popCover{g}(ind) * maxLim; % maximum ART coverage by age
            popCoverInd = popCover{g}(ind); % desired population-level ART coverage
        elseif year >= 2026
            minCoverLim = popCover{g}(end) * minLim;
            maxCoverLim = popCover{g}(end) * maxLim;
            popCoverInd = popCover{g}(end);
        end     
        % Calculate treat/artOut matrices to maintain ART coverage min/max by age
        [artOut , treat , maxAges , excMaxAges , minAges , excMinAges] = ...
            artMinMax(artOut , treat , minCoverLim , maxCoverLim , ...
            fracARTAge , ageVec , g , risk , 11 , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , 11 , ...
            agePopSubTots , dRange);
    end  
end

%%
dPop = zeros(size(pop));
artTreat = zeros(disease , viral , gender , age , risk);
hivDeaths = zeros(gender , age , 1);

% Dropout from PrEP
prepOut = 0; % for now

for g = 1 : gender
    for a = 11 : age
        for r = 1 : risk
            
            hivPositiveArt = hivInds(10 , 6 , g , a , r , :); % allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r))
            
            for v = 1 : 5
                acuteInf = hivInds(2 , v , g , a , r , :); %allcomb(2 , v , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r));
                
                dPop(acuteInf) = dPop(acuteInf) ...
                    - (muHIV(a , 2) + kCD4(g , v , 1) + treat(2 , v , g , a , r)) ... % out: disease related mortality, ART coverage, CD4 progression. removed pie(2 , v , g , a , r)
                    .* pop(acuteInf) + artOut(g , a , r) * artDist(2 , v , g , a , r) ... % Distributed dropouts from ART
                    .* pop(hivPositiveArt);
                hivDeaths(g , a) = hivDeaths(g , a) + sumall(muHIV(a , 2) .* pop(acuteInf));
                
                artTreat(2 , v , g , a , r) = treat(2 , v , g , a , r) ... % keep track of distribution of people going on ART
                        .* sumall(pop(acuteInf)); % going on ART
                    
                % HIV-positive going on ART (d = 10)
                dPop(hivPositiveArt) = ...
                    dPop(hivPositiveArt)...
                    + treat(2 , v , g , a , r)... % rate of going on ART
                    .* pop(acuteInf); % going on ART    
                
                %                 artTreat(2 , v , g , a , r) = pie(2 , v , g , a , r) ... % keep track of distribution of people going on ART
                %                         .* sumall(pop(acuteInf)); % going on ART
                %
                % CD4 > 500 cells/uL (d = 3)
                % CD4 500-350 cells/uL (d = 4)
                % CD4 350-200 cells/uL (d = 5)
                % CD4 <200 cells/uL (d = 6)
                
                for d = 3 : 6
                    % get indices
                    cd4Curr = hivInds(d , v , g , a , r , :);
                    
                    cd4Prev = hivInds(d - 1 , v , g , a , r , :);
                    kCD4_next = 0;
                    if d ~= 6
                        kCD4_next = kCD4(g , v , d - 1); %  progression to next disease state (when d = 6 , kOn = 0 , else kOn = 1)
                    end
                    % calculate CD4 changes
                    dPop(cd4Curr) = ...
                        kCD4(g , v , d - 2) * pop(cd4Prev) ... % CD4 progression from previous disease state
                        + artOut(g , a , r) * artDist(d , v , g , a , r) ... % Distributed dropouts from ART
                        .* pop(hivPositiveArt)...
                        - (kCD4_next ... % progression to next disease state (when d = 6 , kOn = 0 , else kOn = 1)
                        + muHIV(a , d) ... % disease-related mortality
                        + treat(d , v , g , a , r))... % going on ART
                        .* pop(cd4Curr);
                    
                    hivDeaths(g , a) = hivDeaths(g , a) + sumall(muHIV(a , d) .* pop(cd4Curr));
                    
                    artTreat(d , v , g , a , r) = treat(d , v , g , a , r) ... % keep track of distribution of people going on ART
                        .* sumall(pop(cd4Curr)); % going on ART
                    
                    % HIV-positive going on ART (d = 10)
                    dPop(hivPositiveArt) = ...
                        dPop(hivPositiveArt)...
                        + treat(d , v , g , a , r)... % rate of going on ART
                        .* pop(cd4Curr); % going on ART
                end
            end
            
            % HIV Negative, uncircumcised (d = 1)
            hivNegative = hivInds(1 , 1 , g , a , r , :);  %allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r)
            
            hivNegativePrEP = hivInds(9 , 1 , g , a , r , :); % allcomb(9 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r))
            
            dPop(hivNegative) = prepOut... % rate of PrEP dropout
                .* (pop(hivNegativePrEP)) ... % PrEP dropouts from d = 9 (HIV-negative, uncircumcised and on PrEP)
                - treat(1 , 1 , g , a , r)... % rate of going on PrEP
                .* pop(hivNegative);
            
            % HIV-negative, circumcised, and no PrEP (d = 7)
            hivNegCirc = hivInds(7 , 1 , g , a , r , :); % allcomb(7 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r))
            
            hivNegCircPrep = hivInds(8 , 1 , g , a , r , :); % toInd(allcomb(8 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r));
            
            dPop(hivNegCirc) = prepOut ... % PrEP dropouts from d = 8
                * pop(hivNegCircPrep) - (treat(1 , 1 , g , a , r))...
                * pop(hivNegCirc);
            
            % HIV-negative, circumcised and on PrEP (d = 8)
            dPop(hivNegCircPrep) = treat(1 , 1 , g , a ,r)...
                * pop(hivNegCirc) - prepOut ...
                * pop(hivNegCircPrep); % prepOut to d = 1
            
            % HIV-negative, uncircumcised and on PrEP (d = 9)
            dPop(hivNegativePrEP) = treat(1 , 1 , g , a ,r)...
                * pop(hivNegative) - prepOut ...
                * pop(hivNegativePrEP); % prepOut to d = 1
            
            % ART dropout weighted by proportion of individuals with CD4 above/below 200 receiving ART
            % Dropout from ART (d = 10)
            dPop(hivPositiveArt) = ...
                dPop(hivPositiveArt)...
                - artOut(g , a , r) .* pop(hivPositiveArt); % artOut to d = 2:6 as determined by distribution matrix
        end
    end
end

%% Advance viral load
if size(pop , 1) ~= size(vlAdvancer , 2)
    pop = pop';
end
vlAdvanced = vlAdvancer * pop;

if size(vlAdvanced , 1) ~= size(dPop , 1)
    vlAdvanced = vlAdvanced';
end

dPop = dPop + vlAdvanced;

extraOuts{1} = hivDeaths;
extraOuts{2} = artTreat; %reshape(artTreat , [numel(artTreat) , 1]);
dPop = sparse(dPop);
