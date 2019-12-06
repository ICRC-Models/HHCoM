% Births and deaths module
% Simulates births and non-disease related deaths in a
% population.
% Accepts a population matrix as input and returns dPop, a matrix of
% derivatives that describes the change in the population's subgroups due
% to births and deaths.

% Aging module
% Ages population.
% Aged cohort is redistributed into new risk groups upon entry into new age
% group. Risk group status in previous age group does not affect risk group
% status in new age group. Risk redistribution is solely age dependent.
% Accepts a population vector as input and returns dPop, a vector of
% derivatives that describes the change in the population's subgroups due
% to aging.

function [dPop , extraOut] = bornAgeDieRisk(t , pop , year , ...
        gender , age , fivYrAgeGrpsOn , fertMat , fertMat2 , fertMat3 , hivFertPosBirth ,...
        hivFertNegBirth , hivFertPosBirth2 , hivFertNegBirth2 , hivFertPosBirth3 , ...
        hivFertNegBirth3 , deathMat , circMat , circMat2 , ...
        MTCTRate , circStartYear , ageInd , riskInd , riskDist , ...
        stepsPerYear , currYear , agesComb , noVaxScreen , noVaxXscreen , ...
        vaxScreen , vaxXscreen , hpvScreenStartYear)

%% Calculate MTCT rate
kHiv = MTCTRate(1); % year <= 2004
% linearly increase MTCT rate from 2004 to 2005, 2005 to 2008. Constant after 2008
if year > 2008
    kHiv = MTCTRate(3);
elseif year > 2005
    yrs = 2005 : 1 / stepsPerYear : 2008;
    mtctVec = linspace(MTCTRate(2) , MTCTRate(3) , length(yrs));
    ind = yrs == year;
    kHiv = mtctVec(ind);
elseif year > 2004
    yrs = 2004 : 1 / stepsPerYear : 2005;
    mtctVec  = linspace(MTCTRate(1) , MTCTRate(2) , length(yrs));
    ind = yrs == year;
    kHiv = mtctVec(ind);
end

%% Calculate births: HIV-negative and HIV-positive births
if year > 1990 && year <= 2010
    dt = (year - 1990) * stepsPerYear;
    dFertPos = (hivFertPosBirth2 - hivFertPosBirth) ...
        ./ ((2010 - 1990) * stepsPerYear);
    hivFertPosBirth = hivFertPosBirth + dFertPos .* dt;
    dFertNeg = (hivFertNegBirth2 - hivFertNegBirth) ...
        ./ ((2010 - 1990) * stepsPerYear);
    hivFertNegBirth = hivFertNegBirth + dFertNeg .* dt;
    dFertMat = (fertMat2 - fertMat) ...
        ./ ((2010 - 1990) * stepsPerYear);
    fertMat = fertMat + dFertMat .* dt;
elseif year > 2010 && year <=2020
    dt = (year - 2010) * stepsPerYear;
    dFertPos = (hivFertPosBirth3 - hivFertPosBirth2) ...
        ./ ((2020 - 2010) * stepsPerYear);
    hivFertPosBirth = hivFertPosBirth2 + dFertPos .* dt;
    dFertNeg = (hivFertNegBirth3 - hivFertNegBirth2) ...
        ./ ((2020 - 2010) * stepsPerYear);
    hivFertNegBirth = hivFertNegBirth2 + dFertNeg .* dt;
    dFertMat = (fertMat3 - fertMat2) ...
        ./ ((2020 - 2010) * stepsPerYear);
    fertMat = fertMat2 + dFertMat .* dt;
elseif year > 2020
    hivFertPosBirth = hivFertPosBirth3;
    hivFertNegBirth = hivFertNegBirth3;
    fertMat = fertMat3;
end

hivFertPosBirth = hivFertPosBirth .* kHiv;
hivFertNegBirth = hivFertNegBirth .* (1 - kHiv);

if size(pop , 1) ~= size(fertMat , 2)
    pop = pop';
end

births = fertMat * pop + hivFertNegBirth * pop;
hivBirths = hivFertPosBirth * pop;

%% Calculate deaths
deaths = deathMat * pop;

%% Calculate males receiving neonatal circumcision
circBirths = births * 0;
if (year > circStartYear) && (year <= currYear)
    circBirths = circMat * births;
elseif (year > currYear) && (year <= 2030)
    dt = (year - currYear) * stepsPerYear;
    dCircMat = (circMat2 - circMat) ...
        ./ ((2030 - currYear) * stepsPerYear);
    circMat = circMat + dCircMat .* dt;
    circBirths = circMat * births;
elseif year > 2030
    circBirths = circMat2 * births;
end

%% Aging and risk proportion redistribution
% Initialize dPop
dPop = zeros(size(pop));

% prospective population after accounting for births, deaths, and circumcision
prosPop = pop + circBirths + births + hivBirths + deaths;

for g = 1 : gender
    for a = 2 : age
        aPrev = ageInd(g , a - 1 , :);
        aCurr = ageInd(g , a , :);
        
        r1 = riskInd(g , a - 1 , 1 , :);
        r2 = riskInd(g , a - 1 , 2 , :);
        r3 = riskInd(g , a - 1 , 3 , :);
        r1To = riskInd(g , a , 1 , :);
        r2To = riskInd(g , a , 2 , :);
        r3To = riskInd(g , a , 3 , :);
        
        popR1Tot = sumall(pop(r1));
        popR2Tot = sumall(pop(r2));
        popR3Tot = sumall(pop(r3));
        
        % get prospective risk distribution if staying in same risk group when aging
        agedOut = (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* sumall(prosPop(aPrev)); % age 1/5th of previous age group
        agedProsp = agedOut + fivYrAgeGrpsOn * ((4/5) .* sumall(prosPop(aCurr))); % age 1/5th of previous age group into current age group
        riskTarget = agedProsp .* riskDist(a , : , g);
        riskNeed = riskTarget - ((4/5) .* [sumall(prosPop(r1To)) , sumall(prosPop(r2To)) , sumall(prosPop(r3To))]); % numbers needed to fill risk groups
        riskAvail = (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* [popR1Tot , popR2Tot , popR3Tot];
        riskDiff = riskNeed - riskAvail; % difference between numbers needed and available for each risk group
        riskFrac1 = 0;
        riskFrac2 = 0;
        riskFrac3 = 0;
        
        % find fraction of every compartment that must be moved to maintain risk group distribution
        if riskDiff(3) > 0 % if risk 3 deficient
            % start with moving from risk 2 to risk 3
            if riskAvail(2) > 0
                riskFrac2 = min(min(riskDiff(3) , riskAvail(2)) / popR2Tot , 1);
                dPop(r2To) = dPop(r2To) - pop(r2) .* riskFrac2;
                dPop(r3To) = dPop(r3To) + pop(r2) .* riskFrac2;
            end
            % if needed, move from risk 1 to risk 3
            if riskDiff(3) / riskAvail(2) > 1
                riskFrac1 = ...
                    min(min(riskAvail(1) , (riskDiff(3) - riskAvail(2))) / popR1Tot , 1);
                dPop(r1To) = dPop(r1To) - pop(r1) .* riskFrac1;
                dPop(r3To) = dPop(r3To) + pop(r1) .* riskFrac1;
            end
            riskAvail(1) = riskAvail(1) - sum(pop(r1) .* riskFrac1);
            riskAvail(2) = riskAvail(2) - sum(pop(r2) .* riskFrac2);
        end

        if riskDiff(2) > 0 % if risk 2 deficient
            % start with moving from risk 3 to risk 2
            if riskAvail(3) > 0
                riskFrac3 = min(min(riskDiff(2) , riskAvail(3)) / popR3Tot , 1);
                dPop(r3To) = dPop(r3To) - pop(r3) .* 0.99 .* riskFrac3;
                dPop(r2To) = dPop(r2To) + pop(r3) .* 0.99 .* riskFrac3;

            end
            % if needed, move from risk 1 to risk 2
            if riskDiff(2) / riskAvail(3) > 1
                riskFrac1 =...
                    min(min((riskDiff(2) - riskAvail(3)) , riskAvail(1)) / popR1Tot , 1);
                dPop(r1To) = dPop(r1To) - pop(r1) .* riskFrac1;
                dPop(r2To) = dPop(r2To) + pop(r1) .* riskFrac1;
            end
            riskAvail(1) = riskAvail(1) - sum(pop(r1) .* riskFrac1);
            riskAvail(3) = riskAvail(3) - sum(pop(r3) .* riskFrac3);
        end

        if riskDiff(1) > 0 % if risk 1 deficient
            % start with moving from risk 2 to risk 1
            if riskAvail(2) > 0
                riskFrac2 = min(min(riskDiff(1), riskAvail(2)) / popR2Tot , 1);
                dPop(r2To) = dPop(r2To) - pop(r2) .* riskFrac2;
                dPop(r1To) = dPop(r1To) + pop(r2) .* riskFrac2;
            end
            % if needed, move from risk 3 to risk 1
            if riskDiff(1) / riskAvail(2) > 1
                riskFrac3 = ...
                    min(min((riskDiff(1) - riskAvail(2)) , riskAvail(3)) / popR3Tot , 1);
                dPop(r3To) = dPop(r3To) - pop(r3) .* riskFrac3;
                dPop(r1To) = dPop(r1To) + pop(r3) .* riskFrac3;
            end
            riskAvail(2) = riskAvail(2) - sum(pop(r2) .* riskFrac2);
            riskAvail(3) = riskAvail(3) - sum(pop(r3) .* riskFrac3);
        end
        
        dPop(r1To) = dPop(r1To) + (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* pop(r1);
        dPop(r2To) = dPop(r2To) + (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* pop(r2);
        dPop(r3To) = dPop(r3To) + (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* pop(r3); 
        
        dPop(r1) = dPop(r1) - (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* pop(r1);
        dPop(r2) = dPop(r2) - (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* pop(r2);
        dPop(r3) = dPop(r3) - (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* pop(r3);

        % Remove screened status as people age out of screened age groups
        if (year >= hpvScreenStartYear)
            if any((a-1) == agesComb)
                aOut = find((a-1) == agesComb);
                noVaxScreend = dPop(noVaxScreen(:,aOut));
                dPop(noVaxScreen(:,aOut)) = dPop(noVaxScreen(:,aOut)) - noVaxScreend;
                dPop(noVaxXscreen(:,aOut)) = dPop(noVaxXscreen(:,aOut)) + noVaxScreend;

                vaxScreend = dPop(vaxScreen(:,aOut));
                dPop(vaxScreen(:,aOut)) = dPop(vaxScreen(:,aOut)) - vaxScreend;
                dPop(vaxXscreen(:,aOut)) = dPop(vaxXscreen(:,aOut)) + vaxScreend;
            end
        end
    end
    % age last age group
    dPop(r1To) = dPop(r1To) - (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* pop(r1To);
    dPop(r2To) = dPop(r2To) - (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* pop(r2To);
    dPop(r3To) = dPop(r3To) - (1.0/max(1 , (5*fivYrAgeGrpsOn))) .* pop(r3To);
end

% Account for births, deaths, circumcision, and aging
dPop = dPop + circBirths + births + hivBirths + deaths;

%% Save outputs and convert dPop to a column vector for output to ODE solver
extraOut{1} = abs(deaths);

dPop = sparse(dPop);

