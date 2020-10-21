% % % HIV testing campaign
% % % Calculates proportion of PLWHIV who are undiagnosed and diagnosed
% % % Tracks the number of persons tested
% % % Accepts:
% % % 1) Population matrix (pop)
% % % Returns:
% % %

function [nHivDiag , nHivUndiag] = hivTestPopStats(pop , nHivDiag , ...
    nHivUndiag , deaths , hivDeaths , ccDeath , aged1519 , hivInds , ...
    gar , disease , viral , gender , age , risk , hpvTypeGroups)

%% Calculate diagnosed and undiagnosed HIV
for g = 1 : gender
    nHivNeg = sumall(pop(hivInds(1 : 2 , 1 : viral , g , 1 : age , 1 : risk , :)));
    hivDiagProp = nHivDiag / (nHivDiag(1 , g) + nHivUndiag + nHivNeg);
    hivUndiagProp = nHivUndaig / (nHivDiag + nHivUndiag + nHivNeg);
    
    nHivDiag = nHivDiag + inc * propDiagOneYear ...
        + hivDiagProp * sumall(aged1519) ...
        - hivDiagProp * sumall(deaths(gar(g , 1 : age , 1 : risk , :))) ...
        - sumall(hivDeaths(1 : disease , g , 1 : age)) ...
        - (g - 1) * sumall(ccDeath(1 : disease , 1 : age , 1 : hpvTypeGroups));
    
    nHivUndiag = nHivUndiag + inc * propDiagOneYear ...
        + hivUndiagProp * sumall(aged1519) ...
        - hivUndiagProp * sumall(deaths(gar(g , 1 : age , 1 : risk , :)));
end
