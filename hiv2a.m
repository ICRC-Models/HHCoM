% HIV Progression
% Simulates progression through HIV states and takes treatments such as
% ART, PrEP, etc. into account
% Accepts:
% 1) Population matrix (pop)
% 2) Matrix describing the distribution of individuals (by disease status,
% viral load) that went on ART over the last 20 time steps as input (artDist).
% Returns:
% 1) dPop, a matrix of derivatives that describes the change in the
% population's subgroups due to HIV progression.
% 2) artTreat, a matrix describing the distribution of individuals who went
% on ART according to their disease and viral load status at the time they
% went on treatment.
function[dPop , extraOuts] = hiv2a(t , pop , vlAdvancer , artDist , muHIV , ...
    kCD4 ,  maxRateM1 , maxRateM2 , maxRateF1 , maxRateF2 , disease , ...
    viral , gender , age , risk , k , hivInds , ...
    stepsPerYear , year)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load constants and parameters
%%%
artOut = 0;%0.03; % 3% dropout
toInd = @(x) (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + (x(: , 6) - 1) * k(5) ...
    + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + (x(: , 3) - 1) * k(2) ...
    + (x(: , 2) - 1) * k(1) + x(: , 1);
sumall = @(x) sum(x(:));
artDist = reshape(artDist, [disease , viral , gender , age , risk]);
% disease related mortality
%
% pie [d x v x g x a x r]
% values obtained from k4_mat in MainMain.m of HIV model.
treat = zeros(disease , viral , gender , age ,risk);
treat(1 , : , : , 4 : end , :) = 0; % rate of going on PrEP
treat(2 , : , : , 4 : end , :) = 0; % ART coverage
treat(3 , : , : , 4 : end , :) = 0; % ART coverage
treat(4 , : , : , 4 : end , :) = 0; % ART coverage
treat(5 , : , : , 4 : end , :) = 0; % ART coverage
treat(6 , : , : , 4 : end , :) = 0; % ART coverage
kOn = ones(6 , 1);
kOn(6) = 0;

% if year < 2018
% maxRateM1 = 0.42;%1 - exp(-0.25);
% maxRateM2 = 0.42;%1 - exp(-0.35);
% maxRateF1 = 0.42;%1 - exp(-0.35);
% maxRateF2 = 0.42;%1 - exp(-0.5);
% end

% CD4 <= 200, from 2004 to 2011
if year >= 2003 && year < 2011
    if year >= 2003 && year < 2004 
        yrs = 2003 : 1/stepsPerYear : 2004;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.0 , 0.005731 , length(yrs)) ,...
            linspace(0.0 , 0.007745 , length(yrs))};
    elseif year >= 2004 && year < 2005
        yrs = 2004 : 1/stepsPerYear : 2005;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.005731 , 0.021779 , length(yrs)) ,...
            linspace(0.007745 , 0.029431 , length(yrs))};
    elseif year >= 2005 && year < 2006
        yrs = 2005 : 1/stepsPerYear : 2006;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.021779 ,0.047570 , length(yrs)) ,...
            linspace(0.029431 , 0.064284 , length(yrs))};
    elseif year >= 2006 && year < 2007
        yrs = 2006 : 1/stepsPerYear : 2007;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.047570 , 0.081958 , length(yrs)) ,...
            linspace(0.064284 , 0.110755 , length(yrs))};
    elseif year >= 2007 && year < 2008
        yrs = 2007 : 1/stepsPerYear : 2008;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.081958 , 0.115200 , length(yrs)) ,...
            linspace(0.110755 , 0.155676 , length(yrs))};
    elseif year >= 2008 && year < 2009
        yrs = 2008 : 1/stepsPerYear : 2009;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.115200 , 0.141565 , length(yrs)) ,...
            linspace(0.155676 , 0.191303 , length(yrs))};
    elseif year >= 2009 && year < 2010
        yrs = 2009 : 1/stepsPerYear : 2010;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.141565 , 0.175953 , length(yrs)) ,...
            linspace(0.191303 , 0.237774 , length(yrs))};
    elseif year >= 2010 && year < 2011
        yrs = 2010 : 1/stepsPerYear : 2016;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.175953 , 0.442133 , length(yrs)) ,...
            linspace(0.237774 , 0.597477 , length(yrs))};
    end
    for g = 1 : gender
        onArt = sumall(pop(hivInds(10 , 6 , g , 4 : age , 1 : risk , :)));
        totCD4elig = sumall(pop(hivInds(6 , 1 : 5 , g , 4 : age , 1 : risk , :)));
        totHivPos = sumall(pop(hivInds(2 : 6 , 1 : 5 , g , 4 : age , 1 : risk , :)));
        fracART = onArt * (1 - artOut) / (onArt + totHivPos);
        if fracART < popCover{g}(ind)
            cover = ((popCover{g}(ind) - fracART) ./ (1 - fracART)) .* ... % initiation fraction needed to meet population-level coverage
                (totHivPos ./ totCD4elig); % adjust coverage b/c only initiating persons with eligible CD4 = ((cover*HIVtot)/HIVelig)
            treat(6 , 1 : 5 , g , 4 : age , 1 : risk) = max(cover , 0);
        end   
    end
end

% CD4 <= 350, from 2011 to 2015
if year >= 2011 && year < 2015
    yrs = 2010 : 1/stepsPerYear : 2016;
    ind = (round(yrs , 4) == round(year , 4));
    popCover = {linspace(0.175953 , 0.442133 , length(yrs)) ,...
        linspace(0.237774 , 0.597477 , length(yrs))};
    for g = 1 : gender
        onArt = sumall(pop(hivInds(10 , 6 , g , 4 : age , 1 : risk , :)));
        totCD4elig = sumall(pop(hivInds(5 : 6 , 1 : 5 , g , 4 : age , 1 : risk , :)));
        totHivPos = sumall(pop(hivInds(2 : 6 , 1 : 5 , g , 4 : age , 1 : risk , :)));
        fracART = onArt * (1 - artOut) / (onArt + totHivPos);
        if fracART < popCover{g}(ind)
            cover = ((popCover{g}(ind) - fracART) ./ (1 - fracART)) .* ... % initiation fraction needed to meet population-level coverage
                (totHivPos ./ totCD4elig); % adjust coverage b/c only initiating persons with eligible CD4 = ((cover*HIVtot)/HIVelig)
            treat(5 : 6 , 1 : 5 , g , 4 : age , 1 : risk) = max(cover , 0);
        end   
    end
end

% CD4 <= 500, from 2015 to 2016
if year >= 2015 && year < 2016
    yrs = 2010 : 1/stepsPerYear : 2016;
    ind = (round(yrs , 4) == round(year , 4));
    popCover = {linspace(0.175953 , 0.442133 , length(yrs)) ,...
        linspace(0.237774 , 0.597477 , length(yrs))};
    for g = 1 : gender
        onArt = sumall(pop(hivInds(10 , 6 , g , 4 : age , 1 : risk , :)));
        totCD4elig = sumall(pop(hivInds(4 : 6 , 1 : 5 , g , 4 : age , 1 : risk , :)));
        totHivPos = sumall(pop(hivInds(2 : 6 , 1 : 5 , g , 4 : age , 1 : risk , :)));
        fracART = onArt * (1 - artOut) / (onArt + totHivPos);
        if fracART < popCover{g}(ind)
            cover = ((popCover{g}(ind) - fracART) ./ (1 - fracART)) .* ... % initiation fraction needed to meet population-level coverage
                (totHivPos ./ totCD4elig); % adjust coverage b/c only initiating persons with eligible CD4 = ((cover*HIVtot)/HIVelig)
            treat(4 : 6 , 1 : 5 , g , 4 : age , 1 : risk) = max(cover , 0);
        end         
    end
end
% Any CD4, after 2016
if year >= 2016
    if year >= 2016 && year < 2018
        yrs = 2016 : 1/stepsPerYear : 2018;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.442133 , 0.57 , length(yrs)) ,...
            linspace(0.597477 , 0.68 , length(yrs))};
    elseif year >= 2018
        yrs = 2018 : 1/stepsPerYear : 2025;
        ind = (round(yrs , 4) == round(year , 4));
        popCover = {linspace(0.57 , 0.65 , length(yrs)) ,...
            linspace(0.68 , 0.75 , length(yrs))};
    end
    for g = 1 : gender
        onArt = sumall(pop(hivInds(10 , 6 , g , 4 : age , 1 : risk , :)));
        totCD4elig = sumall(pop(hivInds(2 : 6 , 1 : 5 , g , 4 : age , 1 : risk , :)));
        totHivPos = sumall(pop(hivInds(2 : 6 , 1 : 5 , g , 4 : age , 1 : risk , :)));
        fracART = onArt * (1 - artOut) / (onArt + totHivPos);
        if year < 2025 && fracART < popCover{g}(ind)
            cover = ((popCover{g}(ind) - fracART) ./ (1 - fracART)) .* ... % initiation fraction needed to meet population-level coverage
                (totHivPos ./ totCD4elig); % adjust coverage b/c only initiating persons with eligible CD4 = ((cover*HIVtot)/HIVelig)
            treat(2 : 6 , 1 : 5 , g , 4 : age , 1 : risk) = max(cover , 0);
        elseif year >= 2025 && fracART < popCover{g}(end)
            cover = ((popCover{g}(end) - fracART) ./ (1 - fracART)) .* ... % initiation fraction needed to meet population-level coverage
                (totHivPos ./ totCD4elig); % adjust coverage b/c only initiating persons with eligible CD4 = ((cover*HIVtot)/HIVelig)
            treat(2 : 6 , 1 : 5 , g , 4 : age , 1 : risk) = max(cover , 0);  
        end  
    end
end

% see model notes for index values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dPop = zeros(size(pop));
artTreat = zeros(disease , viral , gender , age , risk);
hivDeaths = zeros(gender , age , 1);

% Dropout from PrEP
prepOut = 0; % for now

for g = 1 : gender
    for a = 4 : age
        for r = 1 : risk
            % HIV Negative, uncircumcised (d = 1)
            hivNegative = hivInds(1 , 1 , g , a , r , :);  %allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r)
            
            hivNegativePrEP = hivInds(9 , 1 , g , a , r , :); % allcomb(9 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r))
            
            dPop(hivNegative) = prepOut... % rate of PrEP dropout
                .* (pop(hivNegativePrEP)) ... % PrEP dropouts from d = 9 (HIV-negative, uncircumcised and on PrEP)
                - treat(1 , 1 , g , a , r)... % rate of going on PrEP
                .* pop(hivNegative);
            
            hivPositiveArt = hivInds(10 , 6 , g , a , r , :); % allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r))
            
            for v = 1 : 5
                acuteInf = hivInds(2 , v , g , a , r , :); %allcomb(2 , v , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r));
                
                dPop(acuteInf) = dPop(acuteInf) ...
                    - (muHIV(a , 2) + kCD4(g , v , 1) + treat(2 , v , g , a , r)) ... % out: disease related mortality, ART coverage, CD4 progression. removed pie(2 , v , g , a , r)
                    .* pop(acuteInf) + artOut * artDist(2 , v , g , a , r) ... % Distributed dropouts from ART
                    .* pop(hivPositiveArt);
                hivDeaths(g , a) = hivDeaths(g , a) + sumall(muHIV(a , 2) .* pop(acuteInf));
                
                artTreat(2 , v , g , a , r) = treat(2 , v , g , a , r) ... % keep track of distribution of people going on ART
                        .* sumall(pop(acuteInf)); % going on ART
                
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
                        + artOut * artDist(d , v , g , a , r) ... % Distributed dropouts from ART
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
                - artOut .* pop(hivPositiveArt) ... % artOut to d = 2:6 as determined by distribution matrix
                + treat(2 , v , g , a , r) * pop(acuteInf); 
        end
    end
end

% Advance viral load
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
