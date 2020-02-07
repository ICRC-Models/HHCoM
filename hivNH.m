% HIV Natural History and ART initiation
% Simulates progression through HIV states and takes treatments such as
% ART into account
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

function[dPop , extraOuts] = hivNH(t , pop , vlAdvancer , muHIV , artDist , ...
    kCD4 ,  maxRateM , maxRateF , minLim , maxLim , disease , viral , ...
    hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk , ...
    ageSexDebut , hivInds , stepsPerYear , year)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
artTreat = zeros(disease , viral , gender , age , risk); %zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk);
hivDeaths = zeros(gender , age , 1);

%% Calculate ART treatment coverage
artOut = zeros(gender , age ,risk);
artDist = reshape(artDist, [disease , viral , gender , age , risk]); % zeros(disease , viral , gender , age , risk);
treat = zeros(disease , viral , gender , age ,risk);

% CD4 <= 200, from 2004 to 2011
if year >= 2004 && year < 2011
    if year >= 2004 && year < 2007
        yrs = 2004 : 1 / stepsPerYear : 2007;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.0 , maxRateM(1) , length(yrs)) ,...
            linspace(0.0 , maxRateF(1) , length(yrs))};
    elseif year >= 2007 && year < 2011
        yrs = 2007 : 1 / stepsPerYear : 2011;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(maxRateM(1) , maxRateM(2) , length(yrs)) ,...
            linspace(maxRateF(1) , maxRateF(2) , length(yrs))};
    end
    ageVec = [1 : age];
    dRange = [7];
    for g = 1 : gender
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVSubTots = zeros(1 , age); % number HIV-positives by age
        for a = ageSexDebut : age
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            for v = 1 : 5
                totHivPosAge = totHivPosAge + sumall(pop(hivInds(7 , v , g , a , : , :)));
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
            fracARTAge , ageVec , g , risk , ageSexDebut , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , ageSexDebut , ...
            agePopSubTots , dRange);
    end
end

% CD4 <= 350, from 2011 to 2015
if year >= 2011 && year < 2015
    yrs = 2011 : 1/stepsPerYear : 2016;
    ind = round(yrs , 4) == round(year , 4);
    popCover = {linspace(maxRateM(2) , maxRateM(3) , length(yrs)) , ...
        linspace(maxRateF(2) , maxRateF(3) , length(yrs))};
    ageVec = [1 : age];
    dRange = [6 : 7];
    for g = 1 : gender   
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVSubTots = zeros(1 , age); % number HIV-positives by age
        for a = ageSexDebut : age 
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            for d = 6 : 7
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
            fracARTAge , ageVec , g , risk , ageSexDebut , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , ageSexDebut , ...
            agePopSubTots , dRange);
    end
end

% CD4 <= 500, from 2015 to 2016
if year >= 2015 && year < 2016
    yrs = 2011 : 1/stepsPerYear : 2016;
    ind = round(yrs , 4) == round(year , 4);
    popCover = {linspace(maxRateM(2) , maxRateM(3) , length(yrs)) , ...
        linspace(maxRateF(2) , maxRateF(3) , length(yrs))};
    ageVec = [1 : age];
    dRange = [5 : 7];
    for g = 1 : gender    
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVSubTots = zeros(1 , age); % number HIV-positives by age
        for a = ageSexDebut : age
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            for d = 5 : 7
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
            fracARTAge , ageVec , g , risk , ageSexDebut , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , ageSexDebut , ...
            agePopSubTots , dRange);
    end
end

% Any CD4, after 2016
if year >= 2016
    yrs = 2016 : 1 / stepsPerYear : 2030; % assuming 90-90-90 target reached by 2030
    ind = round(yrs , 4) == round(year , 4);
    popCover = {linspace(maxRateM(3) , maxRateM(4) , length(yrs)) ,...
       linspace(maxRateF(3) , maxRateF(4) , length(yrs))};
    ageVec = [1 : age];
    dRange = [3 : 7];
    for g = 1 : gender
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVSubTots = zeros(1 , age); % number HIV-positives by age
        for a = ageSexDebut : age
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            for d = 3 : 7
                for v = 1 : 5
                    totHivPosAge = totHivPosAge + sumall(pop(hivInds(d , v , g , a , : , :)));
                end
            end
            ageARTSubTots(1 , a) = onArtAge;
            ageHIVSubTots(1 , a) = totHivPosAge;       
        end
        fracARTAge = (ageARTSubTots ./ (ageARTSubTots + ageHIVSubTots)); % fraction on ART by age
        agePopSubTots = ageARTSubTots + ageHIVSubTots; % total HIV-positives (on/off ART) by age
        
        if year < 2030
            minCoverLim = popCover{g}(ind) * minLim; % minimum ART coverage by age
            maxCoverLim = popCover{g}(ind) * maxLim; % maximum ART coverage by age
            popCoverInd = popCover{g}(ind); % desired population-level ART coverage
        elseif year >= 2030
            minCoverLim = popCover{g}(end) * minLim;
            maxCoverLim = popCover{g}(end) * maxLim;
            popCoverInd = popCover{g}(end);
        end     
        % Calculate treat/artOut matrices to maintain ART coverage min/max by age
        [artOut , treat , maxAges , excMaxAges , minAges , excMinAges] = ...
            artMinMax(artOut , treat , minCoverLim , maxCoverLim , ...
            fracARTAge , ageVec , g , risk , ageSexDebut , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , ageSexDebut , ...
            agePopSubTots , dRange);
    end  
end

%% Apply CD4 count progression, HIV-associated mortality, ART treatment, and ART dropout
% for h = 1 : hpvVaxStates
%     for s = 1 : hpvNonVaxStates
%         for x = 1 : endpoints
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk

            hivPositiveArt = hivInds(8 , 6 , g , a , r , :);

            % CD4 progression for HIV-positives with acute infection
            for v = 1 : 5
                acuteInf = hivInds(3 , v , g , a , r , :);
                dPop(acuteInf) = dPop(acuteInf) ...
                    - (muHIV(a , 2) + kCD4(a , 1 , g) + treat(3 , v , g , a , r)) ... % out: disease related mortality, CD4 progression, ART coverage.
                    .* pop(acuteInf) + artOut(g , a , r) * artDist(3 , v , g , a , r) ... % Distributed dropouts from ART
                    .* pop(hivPositiveArt);

                % HIV-positive going on ART (d = 8)
                dPop(hivPositiveArt) = dPop(hivPositiveArt)...
                    + treat(3 , v , g , a , r)... % rate of going on ART
                    .* pop(acuteInf); % going on ART   

                hivDeaths(g , a) = hivDeaths(g , a) + sumall(muHIV(a , 2) .* pop(acuteInf));
                artTreat(3 , v , g , a , r) = treat(3 , v , g , a , r) .* sumall(pop(acuteInf)); % keep track of distribution of people going on ART

                % CD4 progression for HIV-positives advanced to decreased CD4 count
                for d = 4 : 7
                    cd4Curr = hivInds(d , v , g , a , r , :);
                    cd4Prev = hivInds(d - 1 , v , g , a , r , :);
                    kCD4_next = 0;
                    if d ~= 7
                        kCD4_next = kCD4(a , d - 2 , g); %  progression to next disease state
                    end
                    dPop(cd4Curr) = ...
                        kCD4(a , d - 3 , g) * pop(cd4Prev) ... % CD4 progression from previous disease state
                        + artOut(g , a , r) * artDist(d , v , g , a , r) ... % Distributed dropouts from ART
                        .* pop(hivPositiveArt)...
                        - (kCD4_next ... % progression to next disease state
                        + muHIV(a , d - 1) ... % disease-related mortality
                        + treat(d , v , g , a , r))... % going on ART
                        .* pop(cd4Curr);

                    % HIV-positive going on ART (d = 8)
                    dPop(hivPositiveArt) = dPop(hivPositiveArt)...
                        + treat(d , v , g , a , r)... % rate of going on ART
                        .* pop(cd4Curr); % going on ART

                    hivDeaths(g , a) = hivDeaths(g , a) + sumall(muHIV(a , d - 1) .* pop(cd4Curr));
                    artTreat(d , v , g , a , r) = treat(d , v , g , a , r) .* sumall(pop(cd4Curr)); % keep track of distribution of people going on ART                    
                end
            end

            % ART dropout weighted by proportion of individuals with CD4 above/below 200 receiving ART
            % Dropout from ART (d = 8)
            dPop(hivPositiveArt) = dPop(hivPositiveArt)...
                - artOut(g , a , r) .* pop(hivPositiveArt); % artOut to d = 3:7 as determined by distribution matrix
        end
    end
end
%         end
%     end
% end

%% Advance viral load
if size(pop , 1) ~= size(vlAdvancer , 2)
    pop = pop';
end
vlAdvanced = vlAdvancer * pop;

if size(vlAdvanced , 1) ~= size(dPop , 1)
    vlAdvanced = vlAdvanced';
end

dPop = dPop + vlAdvanced;

%% Save outputs and convert dPop to a column vector for output to ODE solver
extraOuts{1} = hivDeaths;
extraOuts{2} = reshape(artTreat , [numel(artTreat) , 1]); % artTreat

dPop = sparse(dPop);

