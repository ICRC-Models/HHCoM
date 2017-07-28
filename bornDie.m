% Births, deaths, vaccination module
% Simulates births, non-disease related deaths, and vaccination in a
% population.
% More comments about births here....
% Accepts a population matrix as input and returns dPop, a matrix of
% derivatives that describes the change in the population's subgroups due
% to births, deaths, and vaccinations.
function[dPop , extraOut] = bornDie(t , pop , year , currStep , age , fertility , fertMat , hivFertMat ,...
    deathMat , circMat , vaxer , mtctVec , circStartYear , vaxStartYear , modelYr1 , stepsPerYear)

%fertility = zeros(age , 6);
year = floor(year);

if currStep > (2004 - modelYr1) * stepsPerYear
    if currStep <= (2005 - modelYr1) * stepsPerYear + 1
        kHiv = mtctVec(1);
    elseif currStep < (2008 - modelYr1) * stepsPerYear + 1
        i = currStep - (2005 - modelYr1) * stepsPerYear;
        kHiv = mtctVec(round(i));
    else
        kHiv = mtctVec(end);
    end
    for d = 2 : 6 % hiv infected
        for v = 1 : 5 % hiv infected
            for a = 1 : age
                hivInfected = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , 1 : periods , 2 , a , 1 : risk));
                hivFertMat(at(posMaleBirth , hivInfected)) = 0.5 * fertility(a , d) * kHiv;
                hivFertMat(at(posFemaleBirth , hivInfected)) = 0.5 * fertility(a , d) * kHiv;
                fertMat(at(negMaleBirth , hivInfected)) = 0.5 * fertility(a , d) * (1 - kHiv);
                fertMat(at(negFemaleBirth , hivInfected)) = 0.5 * fertility(a , d) * (1 - kHiv);
            end
        end
    end
end

if size(pop , 1) ~= size(fertMat , 2)
    pop = pop';
end

births = fertMat * pop;
hivBirths = hivFertMat * pop;
deaths = deathMat * pop;
vaxed = pop * 0;

if year >= vaxStartYear
    vaxed = vaxer * pop;
end

circBirths = births * 0;
if year > circStartYear
    circBirths = circMat * births;
end

extraOut{1} = abs(deaths);
dPop = circBirths + births + hivBirths + deaths + vaxed;