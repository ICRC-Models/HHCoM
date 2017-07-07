% Births, deaths, vaccination module
% Simulates births, non-disease related deaths, and vaccination in a
% population.
% More comments about births here....
% Accepts a population matrix as input and returns dPop, a matrix of
% derivatives that describes the change in the population's subgroups due
% to births, deaths, and vaccinations.

% Aging module
% Ages population.
% 1/5th of the previous age group progresses into the next age group each year.
% Aged cohort is redistributed into new risk groups upon entry into new age
% group. Risk group status in previous age group does not affect risk group
% status in new age group. Risk redistribution is solely age dependent.
% Accepts a population vector as input and returns dPop, a vector of
% derivatives that describes the change in the population's subgroups due
% to aging.

function [dPop , extraOut] = bornAgeDie(t , pop , ager , year , currStep , age , fertility , fertMat , hivFertPosBirth ,...
    hivFertNegBirth , deathMat , circMat , vaxer , MTCTRate , circStartYear , vaxStartYear , startYear , endYear, stepsPerYear)
%% births and deaths
%fertility = zeros(age , 6);
kHiv = MTCTRate(1); % year <= 2004
% linearly increase MTCT rate from 2004 to 2005, 2005 to 2008. Constant
% after 2008
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


hivFertPosBirth = hivFertPosBirth .* kHiv;
hivFertNegBirth = hivFertNegBirth .* (1 - kHiv);

if size(pop , 1) ~= size(fertMat , 2)
    pop = pop';
end

births = fertMat * pop + hivFertNegBirth * pop;
hivBirths = hivFertPosBirth * pop;
deaths = deathMat * pop;
vaxed = pop * 0;
vaxRate = linspace(0 , 0.3 , (endYear - startYear) * stepsPerYear);
yrs = vaxStartYear : 1 / stepsPerYear : endYear;
ind = yrs == year;
if year >= vaxStartYear
    vaxed = vaxer .* vaxRate(ind) * pop;
end

circBirths = births * 0;
if year > circStartYear
    circBirths = circMat * births;
end

%% aging

aged = ager * pop;

extraOut{1} = abs(deaths);
dPop = circBirths + births + hivBirths + deaths + vaxed + aged;