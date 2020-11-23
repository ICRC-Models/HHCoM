% HIV testing campaign
% Calculates proportion of PLWHIV who are undiagnosed and diagnosed
% Tracks the number of persons tested
% Accepts:
% 1) Population matrix (pop)
% Returns:
%

function [nTestedNeg , nTestedUndiag , propHivDiagWCamp , nHivDiag , nHivUndiag] = hivTest(pop , ...
    hivTestCampCov , nHivPos , nHivNeg , nHivUndiag , nHivDiag , ...
    hivInds , viral , gender , age , risk)

%% Initialize output vectors
nTested = zeros(1 , gender);

%% Calculate diagnosed and undiagnosed HIV
for g = 1 : gender
    nTestedUndiag(1 , g) = hivTestCampCov * nHivUndiag(1 , g);
    nTestedNeg(1 , g) = hivTestCampCov * nHivNeg(1 , g);
    nHivDiag(1 , g) = nHivDiag(1 , g) + hivTestCampCov * nHivUndiag(1 , g);
    nHivUndiag(1 , g) = nHivUndiag(1 , g) - hivTestCampCov * nHivUndiag(1 , g);
    propHivDiagWCamp(1 , g) = nHivDiag(1 , g) / (nHivDiag(1 , g) + nHivUndiag(1 , g));
end
 