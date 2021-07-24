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

function[dPop , extraOuts] = hivNH(t , pop , vlAdvancer , muHIV , dMue , mue3 , mue4 , artDist , ...
    kCD4 ,  artYr_vec , artM_vec , artF_vec , minLim , maxLim , disease , viral , ...
    hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk , ...
    ageSexDebut , hivInds , stepsPerYear , year)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
artTreat = zeros(disease , viral , gender , age , risk); %zeros(disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk);
artDiscont = zeros(disease , viral , gender , age , risk);
hivDeaths = zeros(disease , gender , age);

%% Calculate background mortality rate for calculating HIV-associated mortality on ART
if (year < 2003)
    muART = zeros(age , gender);
elseif (year >= 2003) && (year < 2020)
    dt = (year - 2000) * stepsPerYear;
    mueYear = mue3 + dMue .* dt;
elseif (year >= 2020)
    mueYear = mue4;
end

%% Calculate ART treatment coverage
% Note: begin scale-up to the given coverage level in the prior year so that the desired coverage is reached BY that year
artOut = zeros(gender , age ,risk);
artDist = reshape(artDist, [disease , viral , gender , age , risk]); % zeros(disease , viral , gender , age , risk);
treat = zeros(disease , viral , gender , age ,risk);

% CD4 <= 200, from 2004 to 2011
if year >= 2003 && year < 2011
    % Calculate HIV-associated mortality on ART
    muART = 0.5 .* mueYear;
    % Calculate population-level ART coverage
    if year >= 2003 && year < 2004
        ind = (round(artYr_vec{1} , 4) == round(year , 4));
        popCover = {artM_vec{1} , artF_vec{1}};
    elseif year >= 2004 && year < 2005
        ind = (round(artYr_vec{2} , 4) == round(year , 4));
        popCover = {artM_vec{2} , artF_vec{2}};
    elseif year >= 2005 && year < 2006
        ind = (round(artYr_vec{3} , 4) == round(year , 4));
        popCover = {artM_vec{3} , artF_vec{3}};
    elseif year >= 2006 && year < 2007
        ind = (round(artYr_vec{4} , 4) == round(year , 4));
        popCover = {artM_vec{4} , artF_vec{4}};
    elseif year >= 2007 && year < 2008
        ind = (round(artYr_vec{5} , 4) == round(year , 4));
        popCover = {artM_vec{5} , artF_vec{5}};
    elseif year >= 2008 && year < 2009
        ind = (round(artYr_vec{6} , 4) == round(year , 4));
        popCover = {artM_vec{6} , artF_vec{6}};
    elseif year >= 2009 && year < 2010
        ind = (round(artYr_vec{7} , 4) == round(year , 4));
        popCover = {artM_vec{7} , artF_vec{7}};
    elseif year >= 2010 && year < 2011
        ind = (round(artYr_vec{8} , 4) == round(year , 4));
        popCover = {artM_vec{8} , artF_vec{8}};
    end
    ageVec = [1 : age];
    dRange = [7];
    for g = 1 : gender
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVeligSubTots = zeros(1 , age); % number HIV-positives with CD4 eligible for ART by age
        ageHIVallSubTots = zeros(1 , age); % number HIV-positives of all CD4 by age
        for a = ageSexDebut : age
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosEligAge = sumall(pop(hivInds(7 , 1 : 5 , g , a , : , :)));
            totHivPosAllAge = sumall(pop(hivInds(3 : 7 , 1 : 5 , g , a , : , :)));
            ageARTSubTots(1 , a) = onArtAge;
            ageHIVeligSubTots(1 , a) = totHivPosEligAge;
            ageHIVallSubTots(1 , a) = totHivPosAllAge;
        end
        fracARTAge = (ageARTSubTots ./ (ageARTSubTots + ageHIVallSubTots)); % fraction on ART by age
        agePopSubTots = ageARTSubTots + ageHIVallSubTots; % total HIV-positives (on/off ART) by age
        minCoverLim = popCover{g}(ind) * minLim; % minimum ART coverage by age
        maxCoverLim = popCover{g}(ind) * maxLim; % maximum ART coverage by age
        popCoverInd = popCover{g}(ind); % desired population-level ART coverage

        % Calculate treat/artOut matrices to maintain ART coverage min/max by age
        [artOut , treat , maxAges , excMaxAges , minAges , excMinAges] = ...
            artMinMax(artOut , treat , minCoverLim , maxCoverLim , ...
            fracARTAge , ageVec , ageHIVallSubTots , ageHIVeligSubTots , ...
            g , risk , ageSexDebut , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVallSubTots , ageHIVeligSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , ageSexDebut , ...
            agePopSubTots , dRange);
    end
end

% CD4 <= 350, from 2011 to 2015
if year >= 2011 && year < 2015
    % Calculate HIV-associated mortality on ART
    muART = 0.4 .* mueYear; %0.5
    % Calculate population-level ART coverage
   if year >= 2011 && year < 2012
        ind = (round(artYr_vec{9} , 4) == round(year , 4));
        popCover = {artM_vec{9} , artF_vec{9}};
    elseif year >= 2012 && year < 2014
        ind = (round(artYr_vec{10} , 4) == round(year , 4));
        popCover = {artM_vec{10} , artF_vec{10}}; 
%     elseif year >= 2013 && year < 2014
%         ind = (round(artYr_vec{10} , 4) == round(year , 4));
%         popCover = {artM_vec{10} , artF_vec{10}};
    elseif year >= 2014 && year < 2015
        ind = (round(artYr_vec{11} , 4) == round(year , 4));
        popCover = {artM_vec{11} , artF_vec{11}};
   end
    
    ageVec = [1 : age];
    dRange = [6 : 7];
    for g = 1 : gender   
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVeligSubTots = zeros(1 , age); % number HIV-positives with CD4 eligible for ART by age
        ageHIVallSubTots = zeros(1 , age); % number HIV-positives of all CD4 by age
        for a = ageSexDebut : age 
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosEligAge = 0;
            for d = 6 : 7
                totHivPosEligAge = totHivPosEligAge + sumall(pop(hivInds(d , 1 : 5 , g , a , : , :)));
            end
            totHivPosAllAge = sumall(pop(hivInds(3 : 7 , 1 : 5 , g , a , : , :)));
            ageARTSubTots(1 , a) = onArtAge;
            ageHIVeligSubTots(1 , a) = totHivPosEligAge;
            ageHIVallSubTots(1 , a) = totHivPosAllAge;
        end 
        fracARTAge = (ageARTSubTots ./ (ageARTSubTots + ageHIVallSubTots)); % fraction on ART by age
        agePopSubTots = ageARTSubTots + ageHIVallSubTots; % total HIV-positives (on/off ART) by age
        
        minCoverLim = popCover{g}(ind) * minLim; % minimum ART coverage by age
        maxCoverLim = popCover{g}(ind) * maxLim; % maximum ART coverage by age
        popCoverInd = popCover{g}(ind); % desired population-level ART coverage

        % Calculate treat/artOut matrices to maintain ART coverage min/max by age
        [artOut , treat , maxAges , excMaxAges , minAges , excMinAges] = ...
            artMinMax(artOut , treat , minCoverLim , maxCoverLim , ...
            fracARTAge , ageVec , ageHIVallSubTots , ageHIVeligSubTots , ...
            g , risk , ageSexDebut , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVallSubTots , ageHIVeligSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , ageSexDebut , ...
            agePopSubTots , dRange);
    end
end

% CD4 <= 500, from 2015 to 2016
if year >= 2015 && year < 2016
    % Calculate HIV-associated mortality on ART
    muART = 0.25 .* mueYear; %.25
    % Calculate population-level ART coverage
    ind = (round(artYr_vec{12} , 4) == round(year , 4));
    popCover = {artM_vec{12} , artF_vec{12}};
    ageVec = [1 : age];
    dRange = [5 : 7];
    for g = 1 : gender    
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVeligSubTots = zeros(1 , age); % number HIV-positives with CD4 eligible for ART by age
        ageHIVallSubTots = zeros(1 , age); % number HIV-positives of all CD4 by age
        for a = ageSexDebut : age
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosEligAge = 0;
            for d = 5 : 7
                totHivPosEligAge = totHivPosEligAge + sumall(pop(hivInds(d , 1 : 5 , g , a , : , :)));
            end
            totHivPosAllAge = sumall(pop(hivInds(3 : 7 , 1 : 5 , g , a , : , :)));
            ageARTSubTots(1 , a) = onArtAge;
            ageHIVeligSubTots(1 , a) = totHivPosEligAge;
            ageHIVallSubTots(1 , a) = totHivPosAllAge;
        end
        fracARTAge = (ageARTSubTots ./ (ageARTSubTots + ageHIVallSubTots)); % fraction on ART by age
        agePopSubTots = ageARTSubTots + ageHIVallSubTots; % total HIV-positives (on/off ART) by age
        
        minCoverLim = popCover{g}(ind) * minLim; % minimum ART coverage by age
        maxCoverLim = popCover{g}(ind) * maxLim; % maximum ART coverage by age
        popCoverInd = popCover{g}(ind); % desired population-level ART coverage

        % Calculate treat/artOut matrices to maintain ART coverage min/max by age
        [artOut , treat , maxAges , excMaxAges , minAges , excMinAges] = ...
            artMinMax(artOut , treat , minCoverLim , maxCoverLim , ...
            fracARTAge , ageVec , ageHIVallSubTots , ageHIVeligSubTots , ...
            g , risk , ageSexDebut , dRange);
        
        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVallSubTots , ageHIVeligSubTots , ageARTSubTots , maxAges , minAges , ...
            fracARTAge , minCoverLim , maxCoverLim , ageSexDebut , ...
            agePopSubTots , dRange);
    end
end

% Any CD4, after 2016
if year >= 2016
    % Calculate HIV-associated mortality on ART
    muART = 0.15 .* mueYear;
    % Calculate population-level ART coverage
   if year >= 2016 && year < 2017
        ind = (round(artYr_vec{13} , 4) == round(year , 4));
        popCover = {artM_vec{13} , artF_vec{13}};
    elseif year >= 2017 
        ind = (round(artYr_vec{14} , 4) == round(year , 4));
        popCover = {artM_vec{14} , artF_vec{14}}; 
   end
    ageVec = [1 : age];
    dRange = [3 : 7];
    for g = 1 : gender
        ageARTSubTots = zeros(1 , age); % number persons on ART by age
        ageHIVeligSubTots = zeros(1 , age); % number HIV-positives with CD4 eligible for ART by age
        ageHIVallSubTots = zeros(1 , age); % number HIV-positives of all CD4 by age
        for a = ageSexDebut : age
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosEligAge = 0;
            for d = 3 : 7
                totHivPosEligAge = totHivPosEligAge + sumall(pop(hivInds(d , 1 : 5 , g , a , : , :)));
            end
            totHivPosAllAge = sumall(pop(hivInds(3 : 7 , 1 : 5 , g , a , : , :)));
            ageARTSubTots(1 , a) = onArtAge;
            ageHIVeligSubTots(1 , a) = totHivPosEligAge;
            ageHIVallSubTots(1 , a) = totHivPosAllAge;
        end
        fracARTAge = (ageARTSubTots ./ (ageARTSubTots + ageHIVallSubTots)); % fraction on ART by age
        agePopSubTots = ageARTSubTots + ageHIVallSubTots; % total HIV-positives (on/off ART) by age
        
        if year < (2030 - 1)
            minCoverLim = popCover{g}(ind) * minLim; % minimum ART coverage by age
            maxCoverLim = popCover{g}(ind) * maxLim; % maximum ART coverage by age
            popCoverInd = popCover{g}(ind); % desired population-level ART coverage
        elseif year >= (2030 - 1)
            minCoverLim = popCover{g}(end) * minLim;
            maxCoverLim = popCover{g}(end) * maxLim;
            popCoverInd = popCover{g}(end);
        end     
       
        % Calculate treat/artOut matrices to maintain ART coverage min/max by age
        [artOut , treat , maxAges , excMaxAges , minAges , excMinAges] = ...
            artMinMax(artOut , treat , minCoverLim , maxCoverLim , ...
            fracARTAge , ageVec , ageHIVallSubTots , ageHIVeligSubTots , ...
            g , risk , ageSexDebut , dRange);

        % Calculate treat/artOut matrices to maintain population-level ART coverage
        [artOut , treat] = artPopCov(artOut , treat , excMaxAges , ...
            excMinAges , popCoverInd , g , risk , ...
            ageHIVallSubTots , ageHIVeligSubTots , ageARTSubTots , maxAges , minAges , ...
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

                artDiscontDist = (artDist(3 , v , g , a , r) / sumall(artDist(: , : , g , a , r)));
                aInd = a;
                while isnan(artDiscontDist) && (aInd > 1)
                    aInd = aInd - 1;
                    artDiscontDist = (artDist(3 , v , g , aInd , r) / sumall(artDist(: , : , g , aInd , r)));
                end
                artDiscontDist = max(0.0 , artDiscontDist);

                dPop(acuteInf) = dPop(acuteInf) ...
                    - (muHIV(a , 2) + kCD4(a , 1 , g) + treat(3 , v , g , a , r)) ... % out: disease related mortality, CD4 progression, ART coverage.
                    .* pop(acuteInf) + artOut(g , a , r) * artDiscontDist ... % Distributed dropouts from ART
                    .* pop(hivPositiveArt);

                % HIV-positive going on ART (d = 8)
                dPop(hivPositiveArt) = dPop(hivPositiveArt)...
                    + treat(3 , v , g , a , r)... % rate of going on ART
                    .* pop(acuteInf); % going on ART   

                hivDeaths(3 , g , a) = hivDeaths(3 , g , a) + sumall(muHIV(a , 2) .* pop(acuteInf));
                artTreat(3 , v , g , a , r) = treat(3 , v , g , a , r) .* sumall(pop(acuteInf)); % keep track of distribution of people going on ART
                artDiscont(3 , v , g , a , r) = artOut(g , a , r) * artDiscontDist .* sumall(pop(hivPositiveArt));

                % CD4 progression for HIV-positives advanced to decreased CD4 count
                for d = 4 : 7
                    cd4Curr = hivInds(d , v , g , a , r , :);
                    cd4Prev = hivInds(d - 1 , v , g , a , r , :);
                    kCD4_next = 0;
                    if d ~= 7
                        kCD4_next = kCD4(a , d - 2 , g); %  progression to next disease state
                    end

                    artDiscontDist = (artDist(d , v , g , a , r) / sumall(artDist(: , : , g , a , r)));
                    aInd = a;
                    while isnan(artDiscontDist) && (aInd > 1)
                        aInd = aInd - 1;
                        artDiscontDist = (artDist(d , v , g , aInd , r) / sumall(artDist(: , : , g , aInd , r)));
                    end
                    artDiscontDist = max(0.0 , artDiscontDist);
 
                    dPop(cd4Curr) = ...
                        kCD4(a , d - 3 , g) * pop(cd4Prev) ... % CD4 progression from previous disease state
                        + artOut(g , a , r) * artDiscontDist ... % Distributed dropouts from ART
                        .* pop(hivPositiveArt)...
                        - (kCD4_next ... % progression to next disease state
                        + muHIV(a , d - 1) ... % disease-related mortality
                        + treat(d , v , g , a , r))... % going on ART
                        .* pop(cd4Curr);

                    % HIV-positive going on ART (d = 8)
                    dPop(hivPositiveArt) = dPop(hivPositiveArt)...
                        + treat(d , v , g , a , r)... % rate of going on ART
                        .* pop(cd4Curr); % going on ART

                    hivDeaths(d , g , a) = hivDeaths(d , g , a) + sumall(muHIV(a , d - 1) .* pop(cd4Curr));
                    artTreat(d , v , g , a , r) = treat(d , v , g , a , r) .* sumall(pop(cd4Curr)); % keep track of distribution of people going on ART  
                    artDiscont(d , v , g , a , r) = artOut(g , a , r) * artDiscontDist .* sumall(pop(hivPositiveArt));
                end
            end

            % ART dropout weighted by proportion of individuals with CD4 above/below 200 receiving ART
            % Dropout from ART (d = 8)
            dPop(hivPositiveArt) = dPop(hivPositiveArt)...
                - (artOut(g , a , r) + muART(a , g)) .* pop(hivPositiveArt); % artOut to d = 3:7 as determined by distribution matrix
            hivDeaths(8 , g , a) = hivDeaths(8 , g , a) + sumall(muART(a , g) .* pop(hivPositiveArt));
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
extraOuts{2} = artTreat; %reshape(artTreat , [numel(artTreat) , 1]);
%extraOuts{3} = artDiscont;

dPop = sparse(dPop);

