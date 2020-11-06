% HIV testing campaign
% Calculates proportion of PLWHIV who are undiagnosed and diagnosed
% Tracks the number of persons tested
% Accepts:
% 1) Population matrix (pop)
% Returns:
%

function [nTested , propHivDiagWCamp , nHivDiag , nHivUndiag] = hivTest(pop , ...
    hivTestCampCov , propHivDiagWCamp , nHivPos , nHivUndiag , nHivDiag , ...
    hivInds , viral , gender , age , risk)

%% Initialize output vectors
nTested = zeros(1 , gender);

%% Calculate diagnosed and undiagnosed HIV
for g = 1 : gender
    nHivNeg = sumall(pop(hivInds(1 : 2 , 1 : viral , g , 1 : age , 1 : risk , :)));
    
    nHivTestCampTargetPop = nHivNeg + nHivUndiag(1 , g);
    nTested(1 , g) = hivTestCampCov * nHivTestCampTargetPop;
    nHivDiag(1 , g) = nHivDiag(1 , g) + hivTestCampCov * nHivUndiag(1 , g);
    nHivUndiag(1 , g) = nHivUndiag(1 , g) - hivTestCampCov * nHivUndiag(1 , g);
    propHivDiagWCamp(1 , g) = nHivDiag(1 , g) / nHivPos;
end
