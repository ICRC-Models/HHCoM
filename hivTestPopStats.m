% % % HIV testing campaign
% % % Calculates proportion of PLWHIV who are undiagnosed and diagnosed
% % % Tracks the number of persons tested
% % % Accepts:
% % % 1) Population matrix (pop)
% % % Returns:
% % %

function [nHivDiag , nHivUndiag] = hivTestPopStats(pop , propDiagOneYear , ...
    nHivDiag , nHivUndiag , newHiv , deaths , ...
    hivDeaths , ccDeath , aged1519 , aged7579 , hivInds , ...
    disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    gender , age , risk , hpvTypeGroups)

%% Calculate diagnosed and undiagnosed HIV
for g = 1 : gender
    % Account for new HIV cases
    nHivDiag(1 , g) = nHivDiag(1 , g) + propDiagOneYear * sumall(newHiv(1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , g , 4 : age , 1 : risk));
    nHivUndiag(1 , g) = nHivUndiag(1 , g) + (1 - propDiagOneYear) * sumall(newHiv(1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , g , 4 : age , 1 : risk));
    
    % Account for HIV and cervical cancer deaths; persons aging into 15-79 population
    hivDiagProp = nHivDiag(1 , g) / (nHivDiag(1 , g) + nHivUndiag(1 , g));
    
    nHivDiag(1 , g) = nHivDiag(1 , g) - hivDiagProp * sumall(deaths(1 , hivInds(3 : 8 , 1 : viral , g , 4 : age , 1 : risk , :))) ...
        - hivDiagProp * aged7579(1 , g) ...
        - sumall(hivDeaths(1 , 1 : disease , g , 4 : age)) ...
        - (g - 1) * sumall(ccDeath(1 , 3 : 8 , 4 : age , 1 : hpvTypeGroups));
    
    nHivUndiag(1 , g) = nHivUndiag(1 , g) + aged1519(1 , g) ...
        - (1 - hivDiagProp) * aged7579(1 , g) ...
        - (1 - hivDiagProp) * sumall(deaths(1 , hivInds(3 : 8 , 1 : viral , g , 4 : age , 1 : risk , :)));
end
