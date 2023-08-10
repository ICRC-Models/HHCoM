% Note: likeFun only set up to use 5-year age groups

function negSumLogL = likeFun_kSympCalib(ccSymp , ccTreat , stageDist_1997_dObs , startYear , stepsPerYear)

mObs = [];
dMean = [];
dVar = [];

propTimeSpan = [((1997 - startYear) * stepsPerYear +1) : ((1998 -
startYear) * stepsPerYear +6)]; 
local = sum(ccSymp(propTimeSpan, 1, 1:end, 1:end), "all") + sum(ccTreat(propTimeSpan, 1, 1:end, 1:end), "all"); 
regional = sum(ccSymp(propTimeSpan, 2, 1:end, 1:end), "all") + sum(ccTreat(propTimeSpan, 2, 1:end, 1:end), "all"); 
distant = sum(ccSymp(propTimeSpan, 3, 1:end, 1:end), "all") + sum(ccTreat(propTimeSpan, 3, 1:end, 1:end), "all"); 
total = local+regional+distant; 
stageDist = [local/total; regional/total; distant/total]; 

mObs = [mObs ; stageDist]; 
dMean = [dMean ; stageDist_1997_dObs(: , 1)]; 
dVar = [dVar ; stageDist_1997_dObs(: , 2)]; 

%% Likelihood function
logL = -(0.5*log(2*pi)) - (0.5.*log(dVar)) - ((0.5.*(1./dVar)).*(mObs-dMean).^2);
negSumLogL = sum(logL); % CJB note: despite name, positive summed logL, used to be -sum(logL), the negative logL to be minimized
