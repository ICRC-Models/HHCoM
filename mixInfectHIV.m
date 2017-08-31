% Population Mixing and Infection Module
% Simulates mixing and HPV/HIV coinfection in a population.
% Models heterosexual interactions and corresponding transmissions.
% Accepts contact matrix [agecohort x sact x gender] and a population matrix
% as input and returns dPop, a matrix of derivatives that describes the
% change in the population's subgroups.
function [dPop , newInfs] = mixInfectHIV(t , pop , currStep , ...
    gar , hivSus , toHiv , mCurr , fCurr , ...
    mCurrArt , fCurrArt ,epsA_vec , epsR_vec , yr , modelYr1 , ...
    circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
    betaHIVF2M , betaHIVM2F , disease , viral , gender , age , risk , ...
    hpvStates , hpvTypes , k , periods , stepsPerYear , year)
sumall = @(x) sum(x(:));
%% mixInfect Constants
% load workspace and get constants

% process constants
dataYr1 = yr(1);
dataYrLast = yr(size(yr , 1));
% epsAge and epsRisk - extent of assortative mixing
now = currStep / stepsPerYear + modelYr1;
baseYrInd = max(find(now >= yr , 1, 'last') , 1); % get index of first year <= current year
baseYr = yr(baseYrInd);
if currStep < (dataYr1 - modelYr1) * stepsPerYear % assortativity in 1st year
    epsA = epsA_vec{1}(1);
    epsR = epsR_vec{1}(1);
elseif currStep < (dataYrLast - modelYr1) * stepsPerYear % assortativity between 1st and last year
    epsA = epsA_vec{baseYrInd}(currStep - (baseYr - modelYr1) * stepsPerYear + 1);
    epsR = epsR_vec{baseYrInd}(currStep - (baseYr - modelYr1) * stepsPerYear + 1);
else % assortativity in last year
    lastIndA = size(epsA_vec , 1);
    lastIndR = size(epsR_vec , 1);
    epsA = epsA_vec{lastIndA}(size(epsA_vec{lastIndA} , 2));
    epsR = epsR_vec{lastIndR}(size(epsR_vec{lastIndR} , 2));
end

% deltaR and deltaA - nature of assortative mixing (Kronecker delta)
% for all times
deltaR = eye(3 , 3);
% if currStep <= (2005 - modelYr1) * int
deltaAF = eye(16) .* 0.3 + diag(ones(15 , 1) .* 0.7 , 1);
deltaAM = eye(16) .* 0.3 + diag(ones(15 , 1) .* 0.7 , -1);

% after 2005
if currStep > (2005 - modelYr1) * stepsPerYear
    deltaAF = eye(16) .* 0.7 + diag(ones(15 , 1) .* 0.3 , 1);
    deltaAM = eye(16) .* 0.7 + diag(ones(15 , 1) .* 0.3 , -1);
    deltaR = eye(3);
    deltaAM(5 , 4) = 0.6;
    deltaAM(5 , 5) = 0.4;
end
deltaAF(4 , 4) = 0.8;
deltaAF(3 , 4) = 0;
deltaAF(4 , 5) = 0.2;
deltaAF(3 , 3) = 1;

deltaAM(4 , 4) = 0.8;
deltaAM(4 , 5) = 0.2;
deltaAM(4 , 3) = 0;

acts = actsPer; % acts per partnership, from loaded workspace [gender x risk]
% rate of partner change (contact)
% males
c(1 , : , :) = partnersM;
% females
c(2 , : , :) = partnersF;

% protection from circumcision
%circProtect
%condProtect
prepProtect = 0.75; % Ying et al. HIV Home HTC supplementary appendix
artProtect = 0.96; % Ying et al. HIV Home HTC supplementary appendix. change to 0.93?
%% MixInfect
% Mixing matrix
popSum = zeros(gender , age , risk);

for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            popSum(g , a , r) = sumall(pop(gar(g , a , r , :)));
        end
    end
end

ageNum = sum(c .* popSum , 3); % numerator for age portion, sum by r -> dim [g x a]
den = sum(ageNum , 2); % sum across a -> dim [g x 1]

% calculate fractions of population by age group and risk group
ageFraction = bsxfun(@rdivide , ageNum , den); % [g x a]
riskFraction = bsxfun(@rdivide , c .* popSum , ageNum); % [g x a x r]
% if denominator = 0 , set fraction to 0
ageFraction(isnan(ageFraction)) = 0; % 0 / 0 -> 0
riskFraction(isnan(riskFraction)) = 0; % 0 / 0 -> 0

% make age fraction x age fraction square matrix
ageFraction_M = bsxfun(@times , ageFraction(1 , : , :) , ones(age , age));
ageFraction_F = bsxfun(@times , ageFraction(2 , : , :) , ones(age , age));

% initialize risk fraction x risk fraction matrices
riskFraction_M = zeros(age, risk, risk); % [a x r x r]
riskFraction_F = riskFraction_M;
% prepare matrices containing risk info associated with age
for i = 1 : age % create square risk fraction x risk fraction matrices for each age group
    riskFraction_M(i , : , :) = bsxfun(@times , squeeze(riskFraction(1 , i , :)) , ones(risk , risk))'; % [r x r](age)
    riskFraction_F(i , : , :) = bsxfun(@times , squeeze(riskFraction(2 , i , :)) , ones(risk , risk))'; % [r x r](age)
end

% rho
rhoAgeF = epsA .* ageFraction_M + (1 - epsA) .* deltaAF; % [a x a]
rhoAgeM = epsA .* ageFraction_F + (1 - epsA) .* deltaAM; % [a x a]
rhoRiskM = zeros(age , risk , risk);
rhoRiskF = rhoRiskM;

for i = 1 : age
    rhoRiskF(i , : , :) = squeeze(epsR .* riskFraction_M(i , : , :))...
        + (1 - epsR) .* deltaR; % [a(i) x r x r] + [r x r] -> [a x r x r]
    rhoRiskM(i , : , :) = squeeze(epsR .* riskFraction_F(i , : , :))...
        + (1 - epsR) .* deltaR; % [a(i) x r x r] + [r x r] -> [a x r x r]
end


% Intialize rho matrices for males and females
rhoM = zeros(age, age, risk, risk);
rhoF = rhoM;
for i = 1 : age
    for ii = 1 : age
        for j = 1 : risk
            for jj = 1 : risk
                rhoM(i , ii , j , jj) = rhoAgeM(i , ii) * rhoRiskM(ii , j , jj);
                rhoF(i , ii , j , jj) = rhoAgeF(i , ii) * rhoRiskF(ii , j , jj);
            end
        end
    end
end
rho(1 , : , : , : , :) = rhoM;
rho(2 , : , : , : , :) = rhoF;

% calculate discrepancy between male and female reported contacts
mfRatio = zeros(age , age , risk , risk);
for aa = 1 : age
    for rr = 1 : risk
        if popSum(2 , aa , rr) ~= 0
            mfRatio(: , aa , : , rr) = ...
                max(popSum(1 , : , :) ./ popSum(2 , aa ,rr) , 0);
        end
    end
end

B = zeros(age , age , risk , risk);
cMale = partnersM; % [age x risk]
cFemale = partnersF; % [age x risk]
for i = 1 : age
    for ii = 1: age
        for j = 1 : risk
            for jj = 1: risk
                B(i, ii , j , jj) = cMale(i , j) * rhoM(i , ii , j , jj) ...
                    /(cFemale(ii , jj) * rhoF(ii , i , jj , j)) ...
                    * mfRatio(i , ii , j , jj);
            end
        end
    end
end
B(isnan(B)) = 0; % 0/0 , 1/0 -> 0
B(isinf(B)) = 0;
% adjust c to account for discrepancy in reporting
theta = 0.5; % assume adjusted contact rate is mostly compensating for female under-report
cAdjMale = zeros(age , age , risk , risk);
cAdjFemale = cAdjMale;
for i = 1 : age
    for j = 1 : risk
        cAdjMale(i , : , j , :) = cMale(i , j)...
            * B(i , : , j , :) .^ -(1 - theta);

        cAdjFemale(: , i , : , j) = cFemale(: , :)...
            .* squeeze(B(i , : , j , :)) .^ theta;
    end
end
cAdj(1 , : , : , : , :) = cAdjMale;
cAdj(2 , : , : , : , :) = cAdjFemale;
cAdj(isnan(cAdj)) = 0;
cAdj(isinf(cAdj)) = 0;
% HIV average betas
beta = zeros(gender , age , age , risk , risk);

% infection probability by viral load
for a = 1 : age
  for aa = 1 : age
    for r = 1 : risk
        for rr = 1 : risk
            if popSum(1 , a , r) ~= 0
                for v = 1 : 5 % viral load (up to vl = 6). Note: last index is (viral - 1) + 1. Done to align pop index with betaHIV index.
                    beta(1 , a , aa , r , rr) = beta(1 , a , aa , r , rr) - log(1 - betaHIVM2F(aa , rr , v)) ...
                      * sumall(pop(mCurr(a , r , v , :))) ./ popSum(1 , a , r);
                end
                beta(1 , a , aa , r , rr) = beta(1 , a , aa , r , rr) - log(1 - betaHIVM2F(aa , rr , 6)) ...
                * sumall(pop(mCurrArt(a , r , 1 , :)))  ./ popSum(1 , a , r);
            end
            if popSum(2 , a , r) ~= 0
                for v = 1 : 5 % viral load (up to vl = 6). Note: last index is (viral - 1) + 1. Done to align pop index with betaHIV index.
                    beta(2 , a , aa , r , rr) = beta(2 , a , aa , r , rr) - log(1 - betaHIVF2M(aa , rr , v))...
                      * sumall(pop(fCurr(a , r , v , :))) ./ popSum(2 , a , r);
                end
                beta(2 , a , aa , r , rr) = beta(2 , a , aa , r , rr) - log(1 -   betaHIVF2M(aa , rr , 6))...
                  * sumall(pop(fCurrArt(a , r , 1 , :))) ./ popSum(2 , a , r);
            end
        end
    end
  end
end

% lambda
lambda = zeros(gender, age , risk);
for g = 1 : gender
    gg = 2; % partner's gender
    if g == 2
        gg = 1;
    end
    for i = 1 : age
        for j = 1 : risk
            for ii = 1 : age
                for jj = 1 : risk
                    lambda(g , i , j) = ...
                        squeeze(lambda(g , i , j))...
                        + cAdj(g , i , ii , j , jj) * rho(g , i , ii , j , jj)...
                        * beta(gg , ii , i , jj , j);
                end
            end
        end
    end
end
if year >= 2002
    condUse = 0.5* 0.5;
end
cond = 1-(condProtect * condUse); % condom usage and condom protection rates
psi = ones(disease) .* cond;  % vector for protective factors. Scaled to reflect protection by contraception. Currently parameterized for HIV only.
psi(7) = 1 - circProtect .* cond;
psi(8) = (1 - circProtect) * (1 - prepProtect) .* cond;
dPop = zeros(size(pop));
newHiv = zeros(gender , age , risk); % incidence tally by gender
%% Infection
dVec = [1 , 7 : 9];
for a = 1 : age
    for r = 1 : risk % move age and risk iterators below fromState and toState to reduce unneccessary iterations
        if lambda(1 , a , r) > 10 ^ - 6 || lambda(2 , a , r) > 10 ^ -6 % only evaluate if lambda is non-zero
            for i = 1 : length(dVec)
                d = dVec(i);
                mSus = hivSus(d , 1 , a , r , :);
                fSus = hivSus(d , 2 , a , r , :);
                mTo = toHiv(1 , a , r , :); % set disease state = 2 (acute)
                fTo = toHiv(2 , a , r , :);

                mInfected = lambda(1 , a , r)...
                    .* psi(d) .* pop(mSus); % infected males
                fInfected = lambda(2 , a , r)...
                    .* psi(d) .* pop(fSus); % infected females

                % HIV incidence tracker
                newHiv(1 , a , r) = newHiv(1 , a , r) + sumall(mInfected);
                newHiv(2 , a , r) = newHiv(2 , a , r) + sumall(fInfected);

                dPop(mSus) = dPop(mSus) - mInfected; % efflux of infected males
                dPop(fSus) = dPop(fSus) - fInfected; % efflux of infected females

                dPop(mTo) = dPop(mTo) + mInfected; % influx of infected males
                dPop(fTo) = dPop(fTo) + fInfected; % influx of infected females
            end
        end
    end
end

newInfs{1} = newHiv;
% Convert to column vector for output to ODE solver
dPop = sparse(dPop);
