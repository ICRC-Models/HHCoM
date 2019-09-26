% Population Mixing and Infection Module
% Simulates mixing and HPV/HIV coinfection in a population.
% Models heterosexual interactions and corresponding transmissions.
% Accepts contact matrix [agecohort x sact x gender] and a population matrix
% as input and returns dPop, a matrix of derivatives that describes the
% change in the population's subgroups.
function [dPop , newInfs] = mixInfect(t , pop , ...
        gar , perPartnerHpv_vax , perPartnerHpv_nonV , maleActs , ...
        femaleActs , lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , ...
        epsR_vec , yr , circProtect , condProtect , condUse , ...
        partnersM , partnersF , hpv_hivMult , hpvSus , hpvImm , hpvVaxd , ...
        hpvVaxdScreen , hpvVaxd2 , hpvVaxd2Screen , hpvImmVaxd2 , hpvImmVaxd2Screen , ...
        hpvVaxd2NonV , hpvVaxd2NonVScreen , hpvImmVaxd2NonV , hpvImmVaxd2NonVScreen , ...
        hivSus , toHiv , mCurr , fCurr , mCurrArt , fCurrArt , betaHIVF2M , ...
        betaHIVM2F , disease , viral , gender , age , risk , hpvVaxStates , hpvNonVaxStates , ...
        vaxInds , nonVInds , vaxNonVInds , endpoints , intervens , startYear , stepsPerYear , year)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
newHpv = zeros(gender , disease , age , risk);
newImmHpv = newHpv;
newVaxHpv = newHpv;
newHiv = zeros(gender , age , risk);
    
%% Find epsAge and epsRisk according to the present year (extent of assortative mixing) 
% Random mixing (epsilon = 1), mixing proportional to relative sizes of all compartments
% Assortative mixing (epsilon = 0), mixing among groups similar by age/risk
dataYr1 = yr(1);
dataYrLast = yr(size(yr , 1));
baseYrInd = max(find(year >= yr , 1, 'last') , 1); % get index of first year <= current year
baseYr = yr(baseYrInd);
if year <= dataYr1 % assortativity in 1st year
    epsA = epsA_vec{1}(1);
    epsR = epsR_vec{1}(1);
elseif (year < dataYrLast) && (year > dataYr1) % assortativity between 1st and last year
    epsA = epsA_vec{baseYrInd}(round((year - baseYr) * stepsPerYear) + 1);
    epsR = epsR_vec{baseYrInd}(round((year - baseYr) * stepsPerYear) + 1);
elseif year >= dataYrLast % assortativity in last year and after
    lastIndA = size(epsA_vec , 1);
    lastIndR = size(epsR_vec , 1);
    epsA = epsA_vec{lastIndA}(size(epsA_vec{lastIndA} , 2));
    epsR = epsR_vec{lastIndR}(size(epsR_vec{lastIndR} , 2));
end

%% Assign deltaR and deltaA (nature of assortative mixing by age and gender; Kronecker delta)
deltaR = eye(3 , 3);
diagVec = [-2 , -1 , 1 , 2];
deltaAF = eye(80) .* (0.3*0.2);
deltaAM = eye(80) .* (0.3*0.2);
for i = 1 : length(diagVec)
    deltaAF = deltaAF + diag(ones(80-abs(diagVec(i)) , 1) .* (0.3*0.2) , diagVec(i));
    deltaAM = deltaAM + diag(ones(80-abs(diagVec(i)) , 1) .* (0.3*0.2) , diagVec(i));
end
diagVecF = [3 : 1 : 7];
diagVecM = [-3 : -1 : -7];
for j = 1 : length(diagVecF)
    deltaAF = deltaAF + diag(ones(80-abs(diagVecF(j)) , 1) .* (0.7*0.2) , diagVecF(j));
    deltaAM = deltaAM + diag(ones(80-abs(diagVecM(j)) , 1) .* (0.7*0.2) , diagVecM(j));
end

deltaAF(1:10,1:10) = eye(10) .* 0.2;
deltaAM(1:10,1:10) = eye(10) .* 0.2;
for i = 1 : length(diagVec)
    deltaAF(1:10,1:10) = deltaAF(1:10,1:10) + diag(ones(10-abs(diagVec(i)) , 1) .* 0.2 , diagVec(i));
    deltaAM(1:10,1:10) = deltaAM(1:10,1:10) + diag(ones(10-abs(diagVec(i)) , 1) .* 0.2 , diagVec(i));
end
deltaAF(11,9) = 0.0;
deltaAF(11:12,10) = 0.0;
deltaAF(3:12,11:19) = [zeros(8,9); ...
                       ones(1,8).*0.125 , 0.0; ...
                       ones(1,9).*(1/9)];
deltaAM(11,4) = 0.0;
deltaAM(11:12,5) = 0.0;
deltaAM(11:13,6) = 0.0;
deltaAM(11:14,7) = 0.0;
deltaAM(11:15,8) = 0.0;
deltaAM(11:16,9) = 0.0;
deltaAM(11:17,10) = 0.0;
deltaAM(8:17,11:19) = [zeros(3,9); ...
                       1/3 , 1/3 , 1/3 , zeros(1,6); ...
                       0.25 , 0.25 , 0.25 , 0.25 , zeros(1,5); ...
                       0.2 , 0.2 , 0.2 , 0.2 , 0.2 , zeros(1,4); ...
                       0.0 , 0.2 , 0.2 , 0.2 , 0.2 , 0.2 , zeros(1,3); ...
                       0.0 , 0.0 , 0.2 , 0.2 , 0.2 , 0.2 , 0.2 , 0.0 , 0.0; ...
                       ones(1,8).*(1/8) , 0.0; ...
                       ones(1,9).*(1/9)];

%% Calculate mixing matrix rho (pattern of sexual contact by gender, age, risk)
% partnership/ contact matrices
% males
c(1 , : , :) = partnersM;
% females
c(2 , : , :) = partnersF;

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
for i = 11 : age % create square risk fraction x risk fraction matrices for each age group
    riskFraction_M(i , : , :) = bsxfun(@times , squeeze(riskFraction(1 , i , :)) , ones(risk , risk))'; % [r x r](age)
    riskFraction_F(i , : , :) = bsxfun(@times , squeeze(riskFraction(2 , i , :)) , ones(risk , risk))'; % [r x r](age)
end

% rho
rhoAgeF = epsA .* ageFraction_M + (1 - epsA) .* deltaAF; % [a x a]
rhoAgeM = epsA .* ageFraction_F + (1 - epsA) .* deltaAM; % [a x a]
rhoRiskM = zeros(age , risk , risk);
rhoRiskF = rhoRiskM;

for i = 11 : age
    rhoRiskF(i , : , :) = squeeze(epsR .* riskFraction_M(i , : , :))...
        + (1 - epsR) .* deltaR; % [a(i) x r x r] + [r x r] -> [a x r x r]
    rhoRiskM(i , : , :) = squeeze(epsR .* riskFraction_F(i , : , :))...
        + (1 - epsR) .* deltaR; % [a(i) x r x r] + [r x r] -> [a x r x r]
end

% intialize rho matrices for males and females
rhoM = zeros(age, age, risk, risk);
rhoF = rhoM;
for i = 11 : age
    for ii = 11 : age
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

%% Adjust for discrepancies in contact reporting between male and female contacts
% calculate the discrepancy between the two populations
mfRatio = zeros(age , age , risk , risk);
for aa = 11 : age
    for rr = 1 : risk
        if popSum(2 , aa , rr) ~= 0
            mfRatio(11 : age , aa , : , rr) = ...
                max(popSum(1 , 11 : age , :) ./ popSum(2 , aa ,rr) , 0);
        end
    end
end

B = zeros(age , age , risk , risk);
cMale = partnersM; % [age x risk]
cFemale = partnersF; % [age x risk]
for i = 11 : age
    for ii = 11: age
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
theta = 0.5; % contact rate equally driven by contact rates reported by males and females
cAdjMale = zeros(age , age , risk , risk);
cAdjFemale = cAdjMale;
for i = 11 : age
    for j = 1 : risk
        cAdjMale(i , 11 : age , j , :) = cMale(i , j) ...
            .* rhoM(i , 11 : age , j , :)...
            .* B(i , 11 : age , j , :) .^ -(1 - theta) ...
            .* mfRatio(i , 11 : age , j , 1 : risk) .^ theta;
        
        cAdjFemale(11 : age , i , : , j) = cFemale(11 : age , :) ...
            .* squeeze(rhoF(11 : age , i , : , j)) ...
            .* squeeze(B(i , 11 : age , j , :)) .^ theta ...
            .* squeeze(mfRatio(i , 11 : age , j , 1 : risk)) .^ -(1 - theta);
    end
end
cAdj(1 , : , : , : , :) = cAdjMale;
cAdj(2 , : , : , : , :) = cAdjFemale;
cAdj(isnan(cAdj)) = 0;
cAdj(isinf(cAdj)) = 0;

%% HPV Infection
% calculate per-partnership probability of HPV transmission
sexPop = zeros(2 , 1);
sexPop(1) = sumall(popSum(1 , 11 : age , :)); % total sexually active males in population
sexPop(2) = sumall(popSum(2 , 11 : age , :)); % total sexually active females in population
beta = zeros(gender , age , risk , 3);
for g = 1 : gender
    for a = 11 : age
        for r = 1 : risk
            beta_hpvVax_M2F = 1 - (1 - perPartnerHpv_vax) ^ femaleActs(a , r); % per year per partner probability of HPV transmission
            beta_hpvVax_F2M = 1 - (1 - perPartnerHpv_vax) ^ maleActs(a , r);
            beta_hpvNonVax_M2F = 1 - (1 - perPartnerHpv_nonV) ^ femaleActs(a , r);
            beta_hpvNonVax_F2M = 1 - (1 - perPartnerHpv_nonV) ^ maleActs(a , r);
            % Per partnership beta
            beta_hpvVax(1 , a , r) = beta_hpvVax_M2F; % males to infect HPV-negative females [g x r]
            beta_hpvVax(2 , a , r) = beta_hpvVax_F2M;  % females to infect HPV-negative males [g x r]
            beta_hpvNonVax(1 , a , r) = beta_hpvNonVax_M2F; % males to infect HPV-negative females [g x r]
            beta_hpvNonVax(2 , a , r) = beta_hpvNonVax_F2M; % females to infect HPV-negative males [g x r]
            
            vax = vaxInds(g , a , r , :);
            nonV = nonVInds(g , a , r , :);
            vaxNonV = vaxNonVInds(g , a , r , :);
            vaxSum = sumall(pop(vax));
            nonVSum = sumall(pop(nonV));
            vaxNonVSum = sumall(pop(vaxNonV);
            p_vaxHPV = max(vaxSum / popSum(g , a , r) , 0); % proportion of sexually active with vax-type HPV. Max handles divide by 0 error.
            p_nonVHPV = max(nonVSum / popSum(g , a , r) , 0); % proportion of sexually active with non-vax-type HPV
            p_vaxNonVHPV = max(vaxNonVSum / popSum(g , a , r) , 0); % proportion of sexually active with vax-type and non-vax-type HPV
            % adjust betas for HPV transmission to account for proportion of population that is carrying HPV. 
            % Transmission probability throughout population is dependent on the "concentration" of HPV carriers in the population.
            beta(g , a , r , 1) = -log(1 - beta_hpvVax(g , a , r)) .* p_vaxHPV;
            beta(g , a , r , 2) = -log(1 - beta_hpvNonVax(g , a , r)) .* p_nonVHPV;
            beta(g , a , r , 3) = -log(1 - beta_hpvVax(g , a , r)*beta_hpvNonVax(g , a , r)) .* p_vaxNonVHPV;
        end
    end
end

% calculate lambda (HPV force of infection)
lambda = zeros(gender, age , risk , 3);
for g = 1 : gender
    gg = 2; % partner's gender
    if g == 2
        gg = 1;
    end
    for i = 11 : age
        for j = 1 : risk
            for ii = 11 : age
                for jj = 1 : risk
                    for z = 1 : 3
                        lambda(g , i , j , z) = ...
                            lambda(g , i , j , z)...
                            + cAdj(g , i , ii , j , jj)...
                            * beta(gg , ii , jj , z); % [3 x 1]
                    end
                end
            end
        end
    end
end

% calculate new HPV infections
for d = 1 : disease
    for a = 11 : age
        for r = 1 : risk % move age and risk iterators below fromState and toState to reduce unneccessary iterations
            if any(lambda(1 , a , r , 1 : hpvTypeGroups) > 10 ^ -6) || any(lambda(2 , a , r , 1 : hpvTypeGroups) > 10 ^ -6) ... % only evaluate if lambda is non-zero
                
                % Prepare indices
                % susceptible to vaccine-type HPV --> infected with vaccine-type HPV
                mhpvVaxSus = hpvVaxSus(d , 1 , a , r , :); % non-naturally immune
                fhpvVaxSus = hpvVaxSus(d , 2 , a , r , :);
                fhpvVaxImm = hpvVaxImm(d , 2 , a , r , :); % naturally immune, only females have natural immunity
                % vaccine-type HPV infection
                mhpvVaxInf = hpvVaxInf(d , 1 , a , r , :); % update to infected
                fhpvVaxInf = hpvVaxInf(d , 2 , a , r , :);
                
                % susceptible to non-vaccine-type HPV --> infected with non-vaccine-type HPV
                mhpvNonVaxSus = hpvNonVaxSus(d , 1 , a , r , :); % non-naturally immune
                fhpvNonVaxSus = hpvNonVaxSus(d , 2 , a , r , :);
                fhpvNonVaxImm = hpvNonVaxImm(d , 2 , a , r , :); % naturally immune, only females have natural immunity
                % non-vaccine-type HPV infection
                mhpvNonVaxInf = hpvNonVaxInf(d , 1 , a , r , :); % update to infected
                fhpvNonVaxInf = hpvNonVaxInf(d , 2 , a , r , :);
                
                % susceptible to vaccine-type and non-vaccine-type HPV --> infected with both HPV types
                mhpvVaxNonVaxSusSus = hpvVaxNonVaxSusSus(d , 1 , a , r , :);
                fhpvVaxNonVaxSusSus = hpvVaxNonVaxSusSus(d , 1 , a , r , :);
                mhpvVaxNonVaxSusImm = hpvVaxNonVaxSusImm(d , 1 , a , r , :);
                fhpvVaxNonVaxSusImm = hpvVaxNonVaxSusImm(d , 1 , a , r , :);
                mhpvVaxNonVaxImmSus = hpvVaxNonVaxImmSus(d , 1 , a , r , :);
                fhpvVaxNonVaxImmSus = hpvVaxNonVaxImmSus(d , 1 , a , r , :);
                mhpvVaxNonVaxImmImm = hpvVaxNonVaxImmImm(d , 1 , a , r , :);
                fhpvVaxNonVaxImmImm = hpvVaxNonVaxImmImm(d , 1 , a , r , :);
                % vaccine-type and non-vaccine type HPV infection
                mhpvVaxNonVaxInf = hpvVaxNonVaxInf(d , 1 , a , r , :);
                fhpvVaxNonVaxInf = hpvVaxNonVaxInf(d , 1 , a , r , :);
                
                
                % Set lambda multipliers based on CD4 count
                lambdaMultF = 1;
                lambdaMultM = 1;
                if (d > 3) && (d < 8) % CD4 > 500 -> CD4 < 200
                    lambdaMultF = hpv_hivMult(d - 3 , 1); %hTo - 1);
                    lambdaMultM = hpv_hivMult(d - 3 , 1); %hTo - 1);
                elseif d == 8
                    lambdaMultF = artHpvMult; 
                    lambdaMultM = artHpvMult;
                end

                %****************STOPPED EDITING HERE********************
                % Calculate infections
                % non-vaccinated susceptibles (screened or non-screened) (infection rate capped at 0.99) 
                % non-naturally immune
                mInfected = min(lambdaMultM * lambda(1 , a , r , toState)...
                    , 0.999) .* pop(mSus); % infected males
                fInfected = min(lambdaMultF * lambda(2 , a , r , toState)...
                    , 0.999) .* pop(fSus); % infected females
                % naturally immune
                %mInfImm = min(lambdaMultM * lambdaMultImm(a) * lambda(1 , a , r , toState) ...
                %    , 0.999) .* pop(mSusImm);
                fInfImm = min(lambdaMultF * lambdaMultImm(a) * lambda(2 , a , r , toState) ...
                    , 0.999) .* pop(fSusImm);

                % vaccinated susceptibles
                hVax = 1;
                if hTo == 3
                    hVax = 2;
                end
                vaxProtect = max(lambdaMultVax(a , hVax) , 0.0); % oHR has no vaccine protection

                % vaccinated with no infection history, non-naturally immune
                % non-screened
                mInfVax = min(lambdaMultM * vaxProtect * lambda(1 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(mSusVax);
                fInfVax = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(fSusVax);
                % screened
                mInfVaxScreen = min(lambdaMultM * vaxProtect * lambda(1 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(mSusVaxScreen);
                fInfVaxScreen = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(fSusVaxScreen);

                % vaccinated with vaccine-type infection history
                % non-screened
                mInfVax2 = min(lambdaMultM * vaxProtect * lambda(1 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(mSusVax2); % non-naturally immune
                fInfVax2 = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(fSusVax2);
                fInfImmVax2 = min(lambdaMultF * vaxProtect * lambdaMultImm(a) * lambda(2 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(fSusImmVax2); % naturally immune
                % screened
                mInfVax2Screen = min(lambdaMultM * vaxProtect * lambda(1 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(mSusVax2Screen); % non-naturally immune
                fInfVax2Screen = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(fSusVax2Screen);
                fInfImmVax2Screen = min(lambdaMultF * vaxProtect * lambdaMultImm(a) * lambda(2 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(fSusImmVax2Screen); % naturally immune

                % vaccinated with non-vaccine-type infection history
                % non-screened
                mInfVax2NonV = min(lambdaMultM * vaxProtect * lambda(1 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(mSusVax2NonV);
                fInfVax2NonV = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(fSusVax2NonV);
                fInfImmVax2NonV = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(fSusImmVax2NonV);
                % screened
                mInfVax2NonVScreen = min(lambdaMultM * vaxProtect * lambda(1 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(mSusVax2NonVScreen);
                fInfVax2NonVScreen = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(fSusVax2NonVScreen);
                fInfImmVax2NonVScreen = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                    , 0.999 * vaxProtect) .* pop(fSusImmVax2NonVScreen);


                % Incidence tracker
                % non-vaccinated new infections (screened or non-screened)
                % non-naturally immune
                newHpv(1 , d , a , r) = newHpv(1 , d , a , r) + sumall(mInfected);
                newHpv(2 , d , a , r) = newHpv(2 , d , a , r) + sumall(fInfected);
                % naturally immune
                %newImmHpv(1 , d , a , r) = newImmHpv(1 , d , a , r) + sumall(mInfImm);
                newImmHpv(2 , d , a , r) = newImmHpv(2 , d , a , r) + sumall(fInfImm);

                % vaccinated new infections
                newVaxHpv(1 , d , a , r) = newVaxHpv(1 , d , a , r) ...
                    + sumall(mInfVax) + sumall(mInfVaxScreen) + sumall(mInfVax2) + sumall(mInfVax2Screen) ...
                    + sumall(mInfVax2NonV) + sumall(mInfVax2NonVScreen);
                newVaxHpv(2 , d , a , r) = newVaxHpv(2 , d , a , r)...
                    + sumall(fInfVax) + sumall(fInfVaxScreen) + sumall(fInfVax2) + sumall(fInfVax2Screen) ...
                    + sumall(fInfImmVax2) + sumall(fInfImmVax2Screen) + sumall(fInfVax2NonV) + sumall(fInfVax2NonVScreen) ...
                    + sumall(fInfImmVax2NonV) + sumall(fInfImmVax2NonVScreen);


                % Adjust compartments
                % non-vaccinated (screened or non-screened)
                % non-naturally immune
                dPop(mSus) = dPop(mSus) - mInfected; % efflux of infected males
                dPop(fSus) = dPop(fSus) - fInfected; % efflux of infected females
                dPop(mTo) = dPop(mTo) + mInfected; % influx of infected males
                dPop(fTo) = dPop(fTo) + fInfected; % influx of infected females
                % naturally immune
                %dPop(mSusImm) = dPop(mSusImm) - mInfImm;
                dPop(fSusImm) = dPop(fSusImm) - fInfImm;
                %dPop(mToImm) = dPop(mToImm) + mInfImm;
                dPop(fToImm) = dPop(fToImm) + fInfImm;

                %vaccinated  
                dPop(mSusVax) = dPop(mSusVax) - mInfVax; % vaccinated w/ no infection history, non-naturally immune , non-screened
                dPop(fSusVax) = dPop(fSusVax) - fInfVax;

                dPop(mSusVaxScreen) = dPop(mSusVaxScreen) - mInfVaxScreen; % vaccinated w/ no infection history, non-naturally immune , screened
                dPop(fSusVaxScreen) = dPop(fSusVaxScreen) - fInfVaxScreen;

                dPop(mSusVax2) = dPop(mSusVax2) - mInfVax2; % vaccinated w/ vaccine-type infection history, non-naturally immune, non-screened
                dPop(fSusVax2) = dPop(fSusVax2) - fInfVax2;
                dPop(fSusImmVax2) = dPop(fSusImmVax2) - fInfImmVax2; % vaccinated w/ vaccine-type infection history, naturally immune, non-screened

                dPop(mSusVax2Screen) = dPop(mSusVax2Screen) - mInfVax2Screen; % vaccinated w/ vaccine-type infection history, non-naturally immune, screened
                dPop(fSusVax2Screen) = dPop(fSusVax2Screen) - fInfVax2Screen;
                dPop(fSusImmVax2Screen) = dPop(fSusImmVax2Screen) - fInfImmVax2Screen; % vaccinated w/ vaccine-type infection history, naturally immune, screened

                dPop(mSusVax2NonV) = dPop(mSusVax2NonV) - mInfVax2NonV; % vaccinated w/ non-vaccine-type infection history, non-naturally immune, non-screened
                dPop(fSusVax2NonV) = dPop(fSusVax2NonV) - fInfVax2NonV;
                dPop(fSusImmVax2NonV) = dPop(fSusImmVax2NonV) - fInfImmVax2NonV; % vaccinated w/ non-vaccine-type infection history, naturally immune, non-screened

                dPop(mSusVax2NonVScreen) = dPop(mSusVax2NonVScreen) - mInfVax2NonVScreen; % vaccinated w/ non-vaccine-type infection history, non-naturally immune, screened
                dPop(fSusVax2NonVScreen) = dPop(fSusVax2NonVScreen) - fInfVax2NonVScreen;
                dPop(fSusImmVax2NonVScreen) = dPop(fSusImmVax2NonVScreen) - fInfImmVax2NonVScreen; % vaccinated w/ non-vaccine-type infection history, naturally immune, screened

                if hTo == 2 % vaccinated and acquiring vaccine type infection
                    dPop(mToVax) = dPop(mToVax) + mInfVax + mInfVax2 + mInfVax2NonV;
                    dPop(fToVax) = dPop(fToVax) + fInfVax + fInfVax2 + fInfImmVax2 + fInfVax2NonV + fInfImmVax2NonV;

                    dPop(mToVaxScreen) = dPop(mToVaxScreen) + mInfVaxScreen + mInfVax2Screen + mInfVax2NonVScreen;
                    dPop(fToVaxScreen) = dPop(fToVaxScreen) + fInfVaxScreen + fInfVax2Screen + fInfImmVax2Screen + ...
                        fInfVax2NonVScreen + fInfImmVax2NonVScreen;

                elseif hTo == 3 % vaccinated and acquiring non-vaccine type infection
                    dPop(mToVaxNonV) = dPop(mToVaxNonV) + mInfVax + mInfVax2 + mInfVax2NonV;
                    dPop(fToVaxNonV) = dPop(fToVaxNonV) + fInfVax + fInfVax2 + fInfImmVax2 + fInfVax2NonV + fInfImmVax2NonV;

                    dPop(mToVaxNonVScreen) = dPop(mToVaxNonVScreen) + mInfVaxScreen + mInfVax2Screen + mInfVax2NonVScreen;
                    dPop(fToVaxNonVScreen) = dPop(fToVaxNonVScreen) + fInfVaxScreen + fInfVax2Screen + fInfImmVax2Screen + ...
                        fInfVax2NonVScreen + fInfImmVax2NonVScreen;
                end
            end
        end
    end
end

newInfs{1} = newHpv;
newInfs{2} = newImmHpv;
newInfs{3} = newVaxHpv;

%% Protection against HIV due to condom usage and circumcision
% find condom use according to the present year
condStart = 1995;
peakYear = 2000;
yrVec = condStart : 1 / stepsPerYear : peakYear;
condUseVec = linspace(0 , condUse , (peakYear - condStart) * stepsPerYear);
condUse = condUseVec(1); % year <= peakYear
if year < peakYear && year > condStart
    yrInd = year == yrVec;
    condUse = condUseVec(yrInd);
elseif year >= peakYear
    condUse = condUseVec(end);
end
cond = 1-(condProtect * condUse); % condom usage and condom protection rates

% calculate psi vector for protective factors. Scaled to reflect protection by contraception. Currently parameterized for HIV only.
psi = ones(disease) .* cond; % condom use only for all disease states
psi(2) = (1 - circProtect) .* cond; % condom use + circumcision protection for d=2

%% HIV Infection
% HIV average betas
beta = zeros(risk , age , gender , age , risk);

% infection probability by viral load
for a = 1 : age
%     for aa = 1 : age
        for r = 1 : risk
%             for rr = 1 : risk
                % force of infection: males to infect HIV-negative females 
                if popSum(1 , a , r) ~= 0
                    for v = 1 : 5 % viral load (up to vl = 6). Note: last index is (viral - 1) + 1. Done to align pop index with betaHIV index.
                        beta(: , : , 1 , a , r) = beta(: , : , 1 , a , r) - log(1 - betaHIVM2F(: , : , v)) ...
                            * sumall(pop(mCurr(a , r , v , :))) ./ popSum(1 , a , r);
                    end
                    beta(: , : , 1 , a , r) = beta(: , : , 1 , a , r) - log(1 - betaHIVM2F(: , : , 6)) ...
                        * sumall(pop(mCurrArt(a , r , 1 , :)))  ./ popSum(1 , a , r);
                end
                % force of infection: females to infect HIV-negative males 
                if popSum(2 , a , r) ~= 0
                    for v = 1 : 5 % viral load (up to vl = 6). Note: last index is (viral - 1) + 1. Done to align pop index with betaHIV index.
                        beta(: , : , 2 , a , r) = beta(: , : , 2 , a , r) - log(1 - betaHIVF2M(: , : , v))...
                            * sumall(pop(fCurr(a , r , v , :))) ./ popSum(2 , a , r);
                    end
                    beta(: , : , 2 , a , r) = beta(: , : , 2 , a , r) - log(1 - betaHIVF2M(: , : , 6))...
                        * sumall(pop(fCurrArt(a , r , 1 , :))) ./ popSum(2 , a , r);
                end
%             end
        end
%     end
end

% lambda
lambda = zeros(gender, age , risk);
for g = 1 : gender
    gg = 2; % partner's gender
    if g == 2
        gg = 1;
    end
    for i = 11 : age
        for j = 1 : risk
            for ii = 11 : age
                for jj = 1 : risk
                    lambda(g , i , j) = ...
                        lambda(g , i , j)...
                        + cAdj(g , i , ii , j , jj) ...
                        * beta(j , i , gg , ii , jj);
                end
            end
        end
    end
end
dVec = [1 , 7 : 9];
for a = 11 : age
    for r = 1 : risk % move age and risk iterators below fromState and toState to reduce unneccessary iterations
        if lambda(1 , a , r) > 10 ^ - 6 || lambda(2 , a , r) > 10 ^ -6 % only evaluate if lambda is non-zero
            for i = 1 : length(dVec)
                d = dVec(i);
                mSus = hivSus(d , 1 , a , r , :);
                fSus = hivSus(d , 2 , a , r , :);
                mTo = toHiv(1 , a , r , :); % set disease state = 2 (acute)
                fTo = toHiv(2 , a , r , :);

                mInfected = min(lambda(1 , a , r)...
                    .* psi(d) , 0.999) .* pop(mSus); % infected males
                fInfected = min(lambda(2 , a , r)...
                    .* psi(d) , 0.999) .* pop(fSus); % infected females

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

newInfs{4} = newHiv;
% Convert to column vector for output to ODE solver
dPop = sparse(dPop);

