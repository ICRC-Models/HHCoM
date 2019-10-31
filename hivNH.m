% HIV Progression
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

function[dPop , extraOuts] = hivNH(t , pop , vlAdvancer , artDist , muHIV , ...
    kCD4 ,  maxRateM1 , maxRateF1 , maxRateM2 , maxRateF2 , disease , viral , gender , age , risk , ...
    ageSexDebut , hivInds , stepsPerYear , year , sumall)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
artTreat = zeros(disease , viral , gender , age , risk);
hivDeaths = zeros(gender , age , 1);

%% Calculate ART treatment coverage
artOut = 0.0; %0.118; %0.0619; %0.03; % ART dropout rate
artDist = reshape(artDist, [disease , viral , gender , age , risk]);
treat = zeros(disease , viral , gender , age ,risk);

%CD4 > 200, from 2006 to 2013
if year >= 2006 && year < 2013
    yrs = 2006 : 1/stepsPerYear : 2013;
    ind = round(yrs , 4) == round(year , 4);
    for g = 1 : gender
        maxCover = {linspace(0 , maxRateM1 , length(yrs)) , ...
            linspace(0 , maxRateF1 , length(yrs))};
        onArt = sumall(pop(hivInds(8 , 6 , g , ageSexDebut : age , : , :)));
        aList = [];
        ageSubTots = [];
        totHivPosQualAges = 0; 
        toRedisAge = [];
        for a = ageSexDebut : age 
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            totHivPos = 0;
            for d = 3 : 6
                for v = 1 : 5
                    totHivPosAge = totHivPosAge + sumall(pop(hivInds(d , v , g , a , : , :)));
                    totHivPos = totHivPos + sumall(pop(hivInds(d , v , g , ageSexDebut : age , : , :)));
                end
            end
            fracARTAge = onArtAge / (onArtAge + totHivPosAge) * (1 - artOut);
            if fracARTAge >= 1.0
                toRedisAge = [toRedisAge totHivPosAge];
            elseif fracARTAge < 1.0
                aList = [aList a];
                ageSubTots = [ageSubTots totHivPosAge];
                totHivPosQualAges = totHivPosQualAges + totHivPosAge;
            end
        end 
        fracART = onArt / (onArt + totHivPos) * (1 - artOut);
        fracQualAges = ageSubTots ./ totHivPosQualAges;
        if year < 2013 && fracART < maxCover{g}(ind)
            cover = (maxCover{g}(ind) - fracART) ./ (1 - fracART);
            totToRedis = fracQualAges .* sumall(toRedisAge .* max(cover , 0));
            fracToRedis = totToRedis ./ ageSubTots;
            formatToRedis = ones(1 , 1 , 1 , length(aList) , 1);
            formatToRedis(:) = (max(cover , 0) + fracToRedis');
            treat(3 : 6 , 1 : 5 , g , aList , :) = bsxfun(@times , ...
                ones(4 , 5 , 1 , length(aList) , risk) , formatToRedis);
        elseif year >= 2013 && fracART < maxCover{g}(end)
            cover = (maxCover{g}(end) - fracART) ./ (1 - fracART);
            totToRedis = fracQualAges .* sumall(toRedisAge .* max(cover , 0));
            fracToRedis = totToRedis ./ ageSubTots;
            formatToRedis = ones(1 , 1 , 1 , length(aList) , 1);
            formatToRedis(:) = (max(cover , 0) + fracToRedis);
            treat(3 : 6 , 1 : 5 , g , aList , :) = bsxfun(@times , ...
                ones(4 , 5 , 1 , length(aList) , risk) , formatToRedis);
        end
    end
end

% CD4 <= 200, from 2004 to 2006
if year >= 2004 && year < 2013
    yrs = 2004 : 1 / stepsPerYear : 2006;
    ind = (round(yrs , 4) == round(year , 4));
    for g = 1 : gender
        maxCover = {linspace(0 , maxRateM1 , length(yrs)) ,...
            linspace(0 , maxRateF1 , length(yrs))};
        onArt = sumall(pop(hivInds(8 , 6 , g , ageSexDebut : age , : , :)));
        aList = [];
        ageSubTots = [];
        totBelow200QualAges = 0;
        toRedisAge = [];
        for a = ageSexDebut : age
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totBelow200Age = 0;
            totBelow200 = 0;
            for v = 1 : 5
                totBelow200Age = totBelow200Age + sumall(pop(hivInds(7 , v , g , a , : , :)));
                totBelow200 = totBelow200 + sumall(pop(hivInds(7 , v , g , ageSexDebut : age , : , :)));
            end
            fracARTAge = onArtAge / (onArtAge + totBelow200Age) * (1 - artOut);
            if fracARTAge >= 1.0
                toRedisAge = [toRedisAge totBelow200Age];
            elseif fracARTAge < 1.0
                aList = [aList a];
                ageSubTots = [ageSubTots totBelow200Age];
                totBelow200QualAges = totBelow200QualAges + totBelow200Age;
            end
        end
        fracART = onArt / (onArt + totBelow200) * (1 - artOut); 
        fracQualAges = ageSubTots ./ totBelow200QualAges;
        if year < 2006 && fracART < maxCover{g}(ind)
            cover = (maxCover{g}(ind) - fracART) ./ (1 - fracART);
            totToRedis = fracQualAges .* sumall(toRedisAge .* max(cover , 0));
            fracToRedis = totToRedis ./ ageSubTots;
            formatToRedis = ones(1 , 1 , 1 , length(aList) , 1);
            formatToRedis(:) = (max(cover , 0) + fracToRedis');
            treat(7 , 1 : 5 , g , aList , :) = bsxfun(@times , ...
                ones(1 , 5 , 1 , length(aList) , risk) , formatToRedis);
        elseif year >= 2006 && fracART < maxCover{g}(end)
            cover = (maxCover{g}(end) - fracART) ./ (1 - fracART);
            totToRedis = fracQualAges .* sumall(toRedisAge .* max(cover , 0));
            fracToRedis = totToRedis ./ ageSubTots;
            formatToRedis = ones(1 , 1 , 1 , length(aList) , 1);
            formatToRedis(:) = (max(cover , 0) + fracToRedis');
            treat(7 , 1 : 5 , g , aList , :) = bsxfun(@times , ...
                ones(1 , 5 , 1 , length(aList) , risk) , formatToRedis);
        end
    end
end

% CD4 >= 200, from 2013 to 2015
if year >= 2013 && year < 2015
    for g = 1 : gender  
        maxCover = {maxRateM1 , maxRateF1};
        onArt = sumall(pop(hivInds(8 , 6 , g , ageSexDebut : age , : , :)));
        aList = [];
        ageSubTots = [];
        totHivPosQualAges = 0;
        toRedisAge = [];
        for a = ageSexDebut : age
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            totHivPos = 0;
            for d = 3 : 7
                for v = 1 : 5
                    totHivPosAge = totHivPosAge + sumall(pop(hivInds(d , v , g , a , : , :)));
                    totHivPos = totHivPos + sumall(pop(hivInds(d , v , g , ageSexDebut : age , : , :)));
                end
            end
            fracARTAge = onArtAge / (onArtAge + totHivPosAge) * (1 - artOut);
            if fracARTAge >= 1.0
                toRedisAge = [toRedisAge totHivPosAge];
            elseif fracARTAge < 1.0
                aList = [aList a];
                ageSubTots = [ageSubTots totHivPosAge];
                totHivPosQualAges = totHivPosQualAges + totHivPosAge;
            end
        end
        fracART = onArt / (onArt + totHivPos) * (1 - artOut);
        fracQualAges = ageSubTots ./ totHivPosQualAges;
        if fracART < maxCover{g}
            cover = (maxCover{g} - fracART) ./ (1 - fracART);
            totToRedis = fracQualAges .* sumall(toRedisAge .* max(cover , 0));
            fracToRedis = totToRedis ./ ageSubTots;
            formatToRedis = ones(1 , 1 , 1 , length(aList) , 1);
            formatToRedis(:) = (max(cover , 0) + fracToRedis');
            treat(3 : 7 , 1 : 5 , g , aList , :) = bsxfun(@times , ...
                ones(5 , 5 , 1 , length(aList) , risk) , formatToRedis);
        end
    end
end

% CD4 >= 200, from 2015 to 2030
if year >= 2015
    yrs = 2015 : 1 / stepsPerYear : 2030; % assuming 90-90-90 target reached by 2030
    ind = round(yrs , 4) == round(year , 4);
    for g = 1 : gender
        maxCover = {linspace(maxRateM1 , maxRateM2 , length(yrs)) ,...
           linspace(maxRateF1 , maxRateF2 , length(yrs))};
        onArt = sumall(pop(hivInds(8 , 6 , g , ageSexDebut : age , : , :)));
        aList = [];
        ageSubTots = [];
        totHivPosQualAges = 0;
        toRedisAge = [];
        for a = ageSexDebut : age
            onArtAge = sumall(pop(hivInds(8 , 6 , g , a , : , :)));
            totHivPosAge = 0;
            totHivPos = 0;
            for d = 3 : 7
                for v = 1 : 5
                    totHivPosAge = totHivPosAge + sumall(pop(hivInds(d , v , g , a , : , :)));
                    totHivPos = totHivPos + sumall(pop(hivInds(d , v , g , ageSexDebut : age , : , :)));
                end
            end
            fracARTAge = onArtAge / (onArtAge + totHivPosAge) * (1 - artOut);
            if fracARTAge >= 1.0
                toRedisAge = [toRedisAge totHivPosAge];
            elseif fracARTAge < 1.0
                aList = [aList a];
                ageSubTots = [ageSubTots totHivPosAge];
                totHivPosQualAges = totHivPosQualAges + totHivPosAge;
            end
        end
        fracART = onArt / (onArt + totHivPos) * (1 - artOut);
        fracQualAges = ageSubTots ./ totHivPosQualAges;
        if year < 2030 && fracART < maxCover{g}(ind)
            cover = (maxCover{g}(ind) - fracART) ./ (1 - fracART);
            totToRedis = fracQualAges .* sumall(toRedisAge .* max(cover , 0));
            fracToRedis = totToRedis ./ ageSubTots;
            formatToRedis = ones(1 , 1 , 1 , length(aList) , 1);
            formatToRedis(:) = (max(cover , 0) + fracToRedis');
            treat(3 : 7 , 1 : 5 , g , aList , :) = bsxfun(@times , ...
                ones(5 , 5 , 1 , length(aList) , risk) , formatToRedis);
        elseif year >= 2030 && fracART < maxCover{g}(end)
            cover = (maxCover{g}(end)- fracART) ./ (1 - fracART);
            totToRedis = fracQualAges .* sumall(toRedisAge .* max(cover , 0));
            fracToRedis = totToRedis ./ ageSubTots;
            formatToRedis = ones(1 , 1 , 1 , length(aList) , 1);
            formatToRedis(:) = (max(cover , 0) + fracToRedis);
            treat(3 : 7 , 1 : 5 , g , aList , :) = bsxfun(@times , ...
                ones(5 , 5 , 1 , length(aList) , risk) , formatToRedis);
        end
    end
end

%% Apply CD4 count progression, HIV-associated mortality, ART treatment, and ART dropout
for g = 1 : gender
    for a = ageSexDebut : age
        for r = 1 : risk
            
            hivPositiveArt = hivInds(8 , 6 , g , a , r , :);
            
            % CD4 progression for HIV-positives with acute infection
            for v = 1 : 5
                acuteInf = hivInds(3 , v , g , a , r , :);
                dPop(acuteInf) = dPop(acuteInf) ...
                    - (muHIV(a , 2) + kCD4(g , v , 1) + treat(2 , v , g , a , r)) ... % out: disease related mortality, CD4 progression, ART coverage.
                    .* pop(acuteInf) + artOut * artDist(2 , v , g , a , r) ... % Distributed dropouts from ART
                    .* pop(hivPositiveArt);
                
                % HIV-positive going on ART (d = 8)
                dPop(hivPositiveArt) = dPop(hivPositiveArt)...
                    + treat(2 , v , g , a , r)... % rate of going on ART
                    .* pop(acuteInf); % going on ART   
                
                hivDeaths(g , a) = hivDeaths(g , a) + sumall(muHIV(a , 2) .* pop(acuteInf));
                artTreat(3 , v , g , a , r) = treat(2 , v , g , a , r) .* sumall(pop(acuteInf)); % keep track of distribution of people going on ART
                
                % CD4 progression for HIV-positives advanced to decreased CD4 count
                for d = 4 : 7
                    cd4Curr = hivInds(d , v , g , a , r , :);
                    cd4Prev = hivInds(d - 1 , v , g , a , r , :);
                    kCD4_next = 0;
                    if d ~= 7
                        kCD4_next = kCD4(g , v , d - 2); %  progression to next disease state
                    end
                    dPop(cd4Curr) = ...
                        kCD4(g , v , d - 3) * pop(cd4Prev) ... % CD4 progression from previous disease state
                        + artOut * artDist(d , v , g , a , r) ... % Distributed dropouts from ART
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
                - artOut .* pop(hivPositiveArt); % artOut to d = 3:7 as determined by distribution matrix
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

%% Save outputs and convert dPop to a column vector for output to ODE solver
extraOuts{1} = hivDeaths;
extraOuts{2} = artTreat; %reshape(artTreat , [numel(artTreat) , 1]);

dPop = sparse(dPop);

