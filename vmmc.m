% Voluntary Male Medical Circumcision (VMMC)
% Scales up VMMC by age in HIV-negative males
% Accepts:
% 1) Population matrix (pop)
% Returns:
% 1) dPop, a matrix of derivatives that describes the change in the
% population's subgroups due to VMMC scale-up.

function[dPop , hivNegCirc] = vmmc(pop , circStartYear , circNatStartYear , ...
    vmmcYr_vec , vmmc_vec , circ_aVec , hivNegNonVMMCinds , hivNegVMMCinds , ...
    ageSexDebut , year)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
hivNegCirc = 0;

%% Calculate VMMC coverage by age
if year >= circStartYear && year < 2000
    ind = (round(vmmcYr_vec{1} , 4) == round(year , 4));
    periodInd = 1;
elseif year >= 2000 && year < 2008
    ind = (round(vmmcYr_vec{2} , 4) == round(year , 4));
    periodInd = 2;
elseif year >= 2008 && year < 2010
    ind = (round(vmmcYr_vec{3} , 4) == round(year , 4));
    periodInd = 3;
elseif year >= 2010 && year < 2012
    ind = (round(vmmcYr_vec{4} , 4) == round(year , 4));
    periodInd = 4;
elseif year >= 2012 && year < 2020
    ind = (round(vmmcYr_vec{5} , 4) == round(year , 4));
    periodInd = 5;
elseif year >= 2020 && year < 2030
    ind = (round(vmmcYr_vec{6} , 4) == round(year , 4));
    periodInd = 6;
elseif year >= 2030
    ind = length(vmmcYr_vec{6});
    periodInd = 6;
end

%% Apply VMMC
ageGroups = 2;
if year >= circNatStartYear
    ageGroups = length(circ_aVec);
end

for aInd = 1 : ageGroups
    a = circ_aVec{aInd};
    fracVMMC = sumall(pop(hivNegVMMCinds(a , :))) / ...
        (sumall(pop(hivNegNonVMMCinds(a , :))) + sumall(pop(hivNegVMMCinds(a , :))));
    if vmmc_vec{periodInd , aInd}(ind) - fracVMMC > 10 ^ -6 % proportion medically circumcised is below target level
        vmmcCover = max(0 , (vmmc_vec{periodInd , aInd}(ind) - fracVMMC) ./ (1 - fracVMMC)); % cirucmcise enough HIV-negative men in age group to reach target
        toCirc = vmmcCover .* pop(hivNegNonVMMCinds(a , :));
        dPop(hivNegNonVMMCinds(a , :)) = dPop(hivNegNonVMMCinds(a , :)) - toCirc;
        dPop(hivNegVMMCinds(a , :)) = dPop(hivNegVMMCinds(a , :)) + toCirc;
        hivNegCirc = hivNegCirc + sumall(toCirc);
    end
end

%% Save outputs and convert dPop to a column vector for output to ODE solver
dPop = sparse(dPop);

