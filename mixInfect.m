% Population Mixing and Infection Module
% Simulates mixing and HPV/HIV coinfection in a population.
% Models heterosexual interactions and corresponding transmissions.
% Accepts contact parameters and a population matrix
% as input and returns dPop, a matrix of derivatives that describes the
% change in the population's subgroups.
function [dPop , newInfs] = mixInfect(t , pop , ...
    stepsPerYear , year , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , gender , ...
    age , risk , fivYrAgeGrpsOn , hpvTypeGroups , ageSexDebut , gar , epsA_vec , epsR_vec , yr , ...
    partnersM , partnersF , partnersMmult, ...
    beta_hpvVax_mod , beta_hpvNonVax_mod , vaxInds , nonVInds , ...
    lambdaMultImm , lambdaMultVax , artHpvMult , hpv_hivMult , ...
    hpvVaxSus , hpvVaxImm , hpvVaxInf , hpvNonVaxSus , hpvNonVaxImm , hpvNonVaxInf , ...
    circProtect , condProtect , condUse , betaHIV_mod , ...
    d_partnersMmult,  ...
    hivSus , toHiv , hivCurr)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
newHpvVax = zeros(gender , disease , age , risk , intervens);
newImmHpvVax = newHpvVax;
newHpvNonVax = newHpvVax;
newImmHpvNonVax = newHpvVax;
newHiv = zeros(hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk);

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

if fivYrAgeGrpsOn
    %% Assign deltaR and deltaA (nature of assortative mixing by age and gender; Kronecker delta)
    deltaR = eye(3 , 3);
    deltaAF = eye(16) .* 0.3 + diag(ones(15 , 1) .* 0.7 , 1);
    deltaAM = eye(16) .* 0.3 + diag(ones(15 , 1) .* 0.7 , -1);
    deltaAF(1 , 1) = 0.0;
    deltaAF(1 , 2) = 0.0;
    deltaAF(2 , 2) = 0.0;
    deltaAF(2 , 3) = 0.0;
    deltaAM(1 , 1) = 0.0;
    deltaAM(2 , 1) = 0.0;
    deltaAM(2 , 2) = 0.0;
    deltaAM(3 , 2) = 0.0;
    
%     deltaR = eye(3 , 3);
%     deltaAF = eye(16) .* 0.3 + diag(ones(15 , 1) .* 0.7 , 1);
%     deltaAM = eye(16) .* 0.3 + diag(ones(15 , 1) .* 0.7 , -1);
%     deltaAF(4 , 4) = 1;
%     deltaAF(3 , 4) = 0;
%     deltaAF(4 , 5) = 0;
%     deltaAF(3 , 3) = 1;
%     deltaAM(4 , 4) = 1;
%     deltaAM(4 , 3) = 0;
%     deltaAM(3 , 2) = 0;
%     deltaAM(3 , 3) = 1;
else
    %% Assign deltaR and deltaA (nature of assortative mixing by age and gender; Kronecker delta)
    deltaR = eye(3 , 3);
    diagVec = [-2 , -1 , 1 , 2];
    deltaAF = eye(age) .* (0.3*0.2);
    deltaAM = eye(age) .* (0.3*0.2);
    for i = 1 : length(diagVec)
        deltaAF = deltaAF + diag(ones(age-abs(diagVec(i)) , 1) .* (0.3*0.2) , diagVec(i));
        deltaAM = deltaAM + diag(ones(age-abs(diagVec(i)) , 1) .* (0.3*0.2) , diagVec(i));
    end
    diagVecF = [3 : 1 : 7];
    diagVecM = [-3 : -1 : -7];
    for j = 1 : length(diagVecF)
        deltaAF = deltaAF + diag(ones(age-abs(diagVecF(j)) , 1) .* (0.7*0.2) , diagVecF(j));
        deltaAM = deltaAM + diag(ones(age-abs(diagVecM(j)) , 1) .* (0.7*0.2) , diagVecM(j));
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
end

%% Calculate mixing matrix rho (pattern of sexual contact by gender, age, risk)
% partnership/ contact matrices
% partnersMmult = [2 4];
if (year >= 1995) && (year < 1999)
    dt = (year - 1995) * stepsPerYear;
    partnersMmult(1) = partnersMmult(1) + d_partnersMmult(1) .* dt;
    partnersMmult(2) = partnersMmult(2) + d_partnersMmult(2) .* dt;
    partnersMmult(3) = partnersMmult(3) + d_partnersMmult(3) .* dt;
elseif year >= 1999
    partnersMmult(1) = 1.5;
    partnersMmult(2) = 2.5;
    partnersMmult(3) = 0.9;
end
partnersM(4:5, 1:3) = partnersM(4:5, 1:3) .* partnersMmult(1);
partnersF(4:5, 1:3) = partnersF(4:5, 1:3) .* partnersMmult(2);
partnersM(6:10, 1:3) = partnersM(6:10, 1:3) .* partnersMmult(3);
partnersF(6:10, 1:3) = partnersF(6:10, 1:3) .* partnersMmult(3);

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
for i = ageSexDebut : age % create square risk fraction x risk fraction matrices for each age group
    riskFraction_M(i , : , :) = bsxfun(@times , squeeze(riskFraction(1 , i , :)) , ones(risk , risk))'; % [r x r](age)
    riskFraction_F(i , : , :) = bsxfun(@times , squeeze(riskFraction(2 , i , :)) , ones(risk , risk))'; % [r x r](age)
end

% rho
rhoAgeF = epsA .* ageFraction_M + (1 - epsA) .* deltaAF; % [a x a]
rhoAgeM = epsA .* ageFraction_F + (1 - epsA) .* deltaAM; % [a x a]
rhoRiskM = zeros(age , risk , risk);
rhoRiskF = rhoRiskM;

for i = ageSexDebut : age
    rhoRiskF(i , : , :) = squeeze(epsR .* riskFraction_M(i , : , :))...
        + (1 - epsR) .* deltaR; % [a(i) x r x r] + [r x r] -> [a x r x r]
    rhoRiskM(i , : , :) = squeeze(epsR .* riskFraction_F(i , : , :))...
        + (1 - epsR) .* deltaR; % [a(i) x r x r] + [r x r] -> [a x r x r]
end

% intialize rho matrices for males and females
rhoM = zeros(age, age, risk, risk);
rhoF = rhoM;
for i = ageSexDebut : age
    for ii = ageSexDebut : age
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
for aa = ageSexDebut : age
    for rr = 1 : risk
        if popSum(2 , aa , rr) ~= 0
            mfRatio(ageSexDebut : age , aa , : , rr) = ...
                max(popSum(1 , ageSexDebut : age , :) ./ popSum(2 , aa ,rr) , 0);
        end
    end
end

B = zeros(age , age , risk , risk);
cMale = partnersM; % [age x risk]
cFemale = partnersF; % [age x risk]
for i = ageSexDebut : age
    for ii = ageSexDebut : age
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
for i = ageSexDebut : age
    for j = 1 : risk
        cAdjMale(i , ageSexDebut : age , j , :) = cMale(i , j) ...
            .* rhoM(i , ageSexDebut : age , j , :)...
            .* B(i , ageSexDebut : age , j , :) .^ -(1 - theta) ...
            .* mfRatio(i , ageSexDebut : age , j , 1 : risk) .^ theta;
        
        cAdjFemale(ageSexDebut : age , i , : , j) = cFemale(ageSexDebut : age , :) ...
            .* squeeze(rhoF(ageSexDebut : age , i , : , j)) ...
            .* squeeze(B(i , ageSexDebut : age , j , :)) .^ theta ...
            .* squeeze(mfRatio(i , ageSexDebut : age , j , 1 : risk)) .^ -(1 - theta);
    end
end
cAdj(1 , : , : , : , :) = cAdjMale;
cAdj(2 , : , : , : , :) = cAdjFemale;
cAdj(isnan(cAdj)) = 0;
cAdj(isinf(cAdj)) = 0;

%% Protection against HIV and HPV due to condom usage and circumcision
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

% calculate psi vectors for protective factors
% HIV
cond_hiv = 1-(condProtect(:,1) .* condUse); % condom usage and condom protection rates
psi_hiv = ones(gender,disease) .* cond_hiv; % condom use only for all disease states
psi_hiv(:,2) = (1 - circProtect(:,1)) .* cond_hiv; % condom use + circumcision protection for d=2
%HPV
cond_hpv = 1-(condProtect(:,2) * condUse); % condom usage and condom protection rates
psi_hpv = ones(gender,disease) .* cond_hpv;
psi_hpv(:,2) = (1 - circProtect(:,2)) .* cond_hpv;

%% HPV Infection
% calculate per-partnership probability of HPV transmission
beta = zeros(risk , age , gender , age , risk , hpvTypeGroups);

for g = 1 : gender
    for a = ageSexDebut : age
        for r = 1 : risk
            for v = 1 : viral
                for x = 1 : 3
                    vax = vaxInds(v , x , g , a , r , :);
                    vaxSum = sumall(pop(vax));
                    if vaxSum ~= 0.0
                        p_vaxHPV = max(vaxSum / popSum(g , a , r) , 0); % proportion of sexually active with vax-type HPV. Max handles divide by 0 error.
                        % adjust beta for HPV transmission to account for proportion of population that is carrying HPV. 
                        % Transmission probability throughout population is dependent on the "concentration" of HPV carriers in the population.
                        beta(: , : , g , a , r , 1) = beta(: , : , g , a , r , 1) - log(1 - beta_hpvVax_mod(: , : , v , x , g)) * p_vaxHPV;
                    end
                    nonV = nonVInds(v , x , g , a , r , :);
                    nonVSum = sumall(pop(nonV));
                    if nonVSum ~= 0.0
                        p_nonVHPV = max(nonVSum / popSum(g , a , r) , 0); % proportion of sexually active with non-vax-type HPV
                        % adjust beta for HPV transmission to account for proportion of population that is carrying HPV. 
                        % Transmission probability throughout population is dependent on the "concentration" of HPV carriers in the population.
                        beta(: , : , g , a , r , 2) = beta(: , : , g , a , r , 2) - log(1 - beta_hpvNonVax_mod(: , : , v , x , g)) * p_nonVHPV;
                    end
                end
            end
        end
    end
end

% calculate lambda (HPV force of infection)
lambda = zeros(gender, age , risk , hpvTypeGroups);
for g = 1 : gender
    gg = 2; % partner's gender
    if g == 2
        gg = 1;
    end
    for i = ageSexDebut : age
        for j = 1 : risk
            for ii = ageSexDebut : age
                for jj = 1 : risk
                    for z = 1 : hpvTypeGroups % 2 hpvTypeGroups: vaccine-type and non-vaccine-type
                        lambda(g , i , j , z) = ...
                            lambda(g , i , j , z)...
                            + cAdj(g , i , ii , j , jj)...
                            * beta(j , i , gg , ii , jj , z);
                    end
                end
            end
        end
    end
end

% calculate new HPV infections
for a = ageSexDebut : age
    for r = 1 : risk
        if any(lambda(1 , a , r , 1 : hpvTypeGroups) > 10 ^ -6) || any(lambda(2 , a , r , 1 : hpvTypeGroups) > 10 ^ -6) % only evaluate if lambda is non-zero
            for d = 1 : disease
                for p = 1 : intervens
                    % Prepare indices
                    % susceptible to vaccine-type HPV --> infected with vaccine-type HPV
                    mhpvVaxSus = hpvVaxSus(d , 1 , a , r , p , :); % non-naturally immune
                    fhpvVaxSus = hpvVaxSus(d , 2 , a , r , p , :);
                    fhpvVaxImm = hpvVaxImm(d , 2 , a , r , p , :); % naturally immune, only females have natural immunity
                    % vaccine-type HPV infection
                    mhpvVaxInf = hpvVaxInf(d , 1 , a , r , p , :); % update to infected
                    fhpvVaxInf = hpvVaxInf(d , 2 , a , r , p , :);

                    % susceptible to non-vaccine-type HPV --> infected with non-vaccine-type HPV
                    mhpvNonVaxSus = hpvNonVaxSus(d , 1 , a , r , p , :); % non-naturally immune
                    fhpvNonVaxSus = hpvNonVaxSus(d , 2 , a , r , p , :);
                    fhpvNonVaxImm = hpvNonVaxImm(d , 2 , a , r , p , :); % naturally immune, only females have natural immunity
                    % non-vaccine-type HPV infection
                    mhpvNonVaxInf = hpvNonVaxInf(d , 1 , a , r , p , :); % update to infected
                    fhpvNonVaxInf = hpvNonVaxInf(d , 2 , a , r , p , :);


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
                    % Set lambda multiplier for vaccination
                    if (p == 2) || (p == 4)
                        vaxProtect = lambdaMultVax(a , 1);
                    else
                        vaxProtect = 1.0; % no protection
                    end

                    % Calculate infections
                    % susceptible to vaccine-type HPV --> infected with vaccine-type HPV (infection rate capped at 0.99) 
                    mInfectedVax = min(lambdaMultM * vaxProtect * psi_hpv(1,d) * lambda(1 , a , r , 1)...
                        , 0.999 * vaxProtect) .* pop(mhpvVaxSus);
                    fInfectedVax = min(lambdaMultF * vaxProtect * psi_hpv(2,d) * lambda(2 , a , r , 1)...
                        , 0.999 * vaxProtect) .* pop(fhpvVaxSus);
                    fInfectedVaxImm = min(lambdaMultF * lambdaMultImm(a) * vaxProtect * psi_hpv(2,d) * lambda(2 , a , r , 1)...
                        , 0.999 * vaxProtect) .* pop(fhpvVaxImm);

                    % susceptible to non-vaccine-type HPV --> infected with non-vaccine-type HPV (infection rate capped at 0.99) 
                    mInfectedNonVax = min(lambdaMultM * psi_hpv(1,d) * lambda(1 , a , r , 2)...
                        , 0.999) .* pop(mhpvNonVaxSus);
                    fInfectedNonVax = min(lambdaMultF * psi_hpv(2,d) * lambda(2 , a , r , 2)...
                        , 0.999) .* pop(fhpvNonVaxSus);
                    fInfectedNonVaxImm = min(lambdaMultF * lambdaMultImm(a) * psi_hpv(2,d) * lambda(2 , a , r , 2)...
                        , 0.999) .* pop(fhpvNonVaxImm);
              

                    % Incidence tracker
                    % susceptible to vaccine-type HPV --> infected with vaccine-type HPV
                    % non-naturally immune
                    newHpvVax(1 , d , a , r , p) = newHpvVax(1 , d , a , r , p) + sumall(mInfectedVax);
                    newHpvVax(2 , d , a , r , p) = newHpvVax(2 , d , a , r , p) + sumall(fInfectedVax);
                    % naturally immune
                    newImmHpvVax(2 , d , a , r , p) = newImmHpvVax(2 , d , a , r , p) + sumall(fInfectedVaxImm);
                    
                    % susceptible to non-vaccine-type HPV --> infected with non-vaccine-type HPV
                    newHpvNonVax(1 , d , a , r , p) = newHpvNonVax(1 , d , a , r , p) + sumall(mInfectedNonVax);
                    newHpvNonVax(2 , d , a , r , p) = newHpvNonVax(2 , d , a , r , p) + sumall(fInfectedNonVax);
                    % naturally immune
                    newImmHpvNonVax(2 , d , a , r , p) = newImmHpvNonVax(2 , d , a , r , p) + sumall(fInfectedNonVaxImm);
                    
                    
                    % Adjust compartments
                    % susceptible to vaccine-type HPV --> infected with vaccine-type HPV
                    dPop(mhpvVaxSus) = dPop(mhpvVaxSus) - mInfectedVax;
                    dPop(fhpvVaxSus) = dPop(fhpvVaxSus) - fInfectedVax;
                    dPop(fhpvVaxImm) = dPop(fhpvVaxImm) - fInfectedVaxImm;
                    
                    dPop(mhpvVaxInf) = dPop(mhpvVaxInf) + mInfectedVax;
                    dPop(fhpvVaxInf) = dPop(fhpvVaxInf) + fInfectedVax + fInfectedVaxImm;
                    
                    % susceptible to non-vaccine-type HPV --> infected with non-vaccine-type HPV
                    dPop(mhpvNonVaxSus) = dPop(mhpvNonVaxSus) - mInfectedNonVax;
                    dPop(fhpvNonVaxSus) = dPop(fhpvNonVaxSus) - fInfectedNonVax;
                    dPop(fhpvNonVaxImm) = dPop(fhpvNonVaxImm) - fInfectedNonVaxImm;
                   
                    dPop(mhpvNonVaxInf) = dPop(mhpvNonVaxInf) + mInfectedNonVax;
                    dPop(fhpvNonVaxInf) = dPop(fhpvNonVaxInf) + fInfectedNonVax + fInfectedNonVaxImm;
                end
            end
        end
    end
end

newInfs{1} = newHpvVax;
newInfs{2} = newImmHpvVax;
newInfs{3} = newHpvNonVax;
newInfs{4} = newImmHpvNonVax;

%% HIV Infection
% HIV average betas
beta = zeros(risk , age , gender , age , risk);

% calculate HIV infection probability by viral load
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            if popSum(g , a , r) ~= 0.0
                for v = 1 : viral
                    for x = 1 : endpoints
                        beta(: , : , g , a , r) = beta(: , : , g , a , r) - log(1 - betaHIV_mod(: , : , v , x , g)) ...
                            * sumall(pop(hivCurr(v , x , g , a , r , :))) ./ popSum(g , a , r);
                    end
                end
            end
        end
    end
end

% calculate lambda (HIV force of infection)
lambda = zeros(gender, age , risk);
for g = 1 : gender
    gg = 2; % partner's gender
    if g == 2
        gg = 1;
    end
    for i = ageSexDebut : age
        for j = 1 : risk
            for ii = ageSexDebut : age
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

% calculate new HIV infections
for h = 1 : hpvVaxStates
    for s = 1 : hpvNonVaxStates
        for x = 1 : endpoints
            for a = ageSexDebut : age
                for r = 1 : risk
                    if lambda(1 , a , r) > 10 ^ - 6 || lambda(2 , a , r) > 10 ^ -6 % only evaluate if lambda is non-zero
                        for d = 1 : 2
                            % Prepare indices
                            mSus = hivSus(d , h , s , x , 1 , a , r , :);
                            fSus = hivSus(d , h , s , x , 2 , a , r , :);
                            mTo = toHiv(h , s , x , 1 , a , r , :); % set disease state = 3 (acute infection)
                            fTo = toHiv(h , s , x , 2 , a , r , :);

                            % Calculate infections
                            mInfected = min(lambda(1 , a , r)...
                                .* psi_hiv(1,d) , 0.999) .* pop(mSus); % infected males
                            fInfected = min(lambda(2 , a , r)...
                                .* psi_hiv(2,d) , 0.999) .* pop(fSus); % infected females

                            % HIV incidence tracker
                            newHiv(h , s , x , 1 , a , r) = newHiv(h , s , x , 1 , a , r) + sumall(mInfected);
                            newHiv(h , s , x , 2 , a , r) = newHiv(h , s , x , 2 , a , r) + sumall(fInfected);

                            % Adjust compartments
                            dPop(mSus) = dPop(mSus) - mInfected; % efflux of infected males
                            dPop(fSus) = dPop(fSus) - fInfected; % efflux of infected females
                            dPop(mTo) = dPop(mTo) + mInfected; % influx of infected males
                            dPop(fTo) = dPop(fTo) + fInfected; % influx of infected females
                        end
                    end
                end
            end
        end
    end
end

%% Save outputs and convert dPop to a column vector for output to ODE solver
newInfs{5} = newHiv;

dPop = sparse(dPop);

