% Population Mixing and Infection Module
% Simulates mixing and HPV/HIV coinfection in a population.
% Models heterosexual interactions and corresponding transmissions.
% Accepts contact matrix [agecohort x sact x gender] and a population matrix
% as input and returns dPop, a matrix of derivatives that describes the
% change in the population's subgroups.
function [dPop , newHpv , newImmHpv , newVaxHpv , newHiv] = mixInfect(pop , ...
        gar , perPartnerHpv , perPartnerHpv_lr , perPartnerHpv_nonV , maleActs , ...
        femaleActs , lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , ...
        epsR_vec , yr , circProtect , condProtect , condUse , actsPer , ...
        partnersM , partnersF , hpv_hivMult , hpvSus , hpvImm , hpvVaxd , ...
        hpvVaxdScreen , hpvVaxd2 , hpvImmVaxd2 , toHpv , toHpv_Imm , toHpv_Vax, ...
        toHpv_VaxScreen , toHpv_VaxNonV , toHpv_VaxNonVScreen , hivSus , toHiv , ...
        mCurr , fCurr , mCurrArt , fCurrArt , betaHIVF2M , betaHIVM2F , disease , ...
        viral , gender , age , risk , hpvStates , hpvTypes , hrInds , lrInds , ...
        hrlrInds , periods , startYear , stepsPerYear , year)
    
sumall = @(x) sum(x(:));
%% Process mixInfect constants

% epsAge and epsRisk - extent of assortative mixing
% % dataYr1 = yr(1);
% % dataYrLast = yr(size(yr , 1));
% % now = currStep / stepsPerYear + startYear;
% % baseYrInd = max(find(now >= yr , 1, 'last') , 1); % get index of first year <= current year
% % baseYr = yr(baseYrInd);
% % if currStep < (dataYr1 - startYear) * stepsPerYear % assortativity in 1st year
% %     epsA = epsA_vec{1}(1);
% %     epsR = epsR_vec{1}(1);
% % elseif currStep < (dataYrLast - startYear) * stepsPerYear % assortativity between 1st and last year
% %     epsA = epsA_vec{baseYrInd}(currStep - (baseYr - startYear) * stepsPerYear + 1);
% %     epsR = epsR_vec{baseYrInd}(currStep - (baseYr - startYear) * stepsPerYear + 1);
% % else % assortativity in last year
% %     lastIndA = size(epsA_vec , 1);
% %     lastIndR = size(epsR_vec , 1);
% %     epsA = epsA_vec{lastIndA}(size(epsA_vec{lastIndA} , 2));
% %     epsR = epsR_vec{lastIndR}(size(epsR_vec{lastIndR} , 2));
% % end
% Note: if reimplementing this, need to fix logic around currStep, see
% condom use or hpv screening vectors as an example
epsA = 0.3;
epsR = 0.3;

% deltaR and deltaA - nature of assortative mixing (Kronecker delta)
% for all times
deltaR = eye(3 , 3);
% if currStep <= (2005 - startYear) * int
% original
deltaAF = eye(16) .* 0.3 + diag(ones(15 , 1) .* 0.7 , 1);
deltaAM = eye(16) .* 0.3 + diag(ones(15 , 1) .* 0.7 , -1);
% after 2005
% if currStep > (2000 - startYear) * stepsPerYear
%     deltaAF = eye(16) .* 0.8 + diag(ones(15 , 1) .* 0.2 , 1);
%     deltaAM = eye(16) .* 0.8 + diag(ones(15 , 1) .* 0.2 , -1);
%     deltaR = eye(3);
% %     deltaAM(5 , 4) = 0.6;
% %     deltaAM(5 , 5) = 0.4;
% end
deltaAF(4 , 4) = 1;
deltaAF(3 , 4) = 0;
deltaAF(4 , 5) = 0;
deltaAF(3 , 3) = 1;
deltaAM(4 , 4) = 1;
deltaAM(4 , 3) = 0;
deltaAM(3 , 2) = 0;
deltaAM(3 , 3) = 1;

acts = actsPer; % acts per partnership, from loaded workspace [gender x risk] (not currently used)

% Rate of partner change (contact)
% males
c(1 , : , :) = partnersM;
% females
c(2 , : , :) = partnersF;

% Protection from circumcision
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
for i = 3 : age % create square risk fraction x risk fraction matrices for each age group
    riskFraction_M(i , : , :) = bsxfun(@times , squeeze(riskFraction(1 , i , :)) , ones(risk , risk))'; % [r x r](age)
    riskFraction_F(i , : , :) = bsxfun(@times , squeeze(riskFraction(2 , i , :)) , ones(risk , risk))'; % [r x r](age)
end

% rho
rhoAgeF = epsA .* ageFraction_M + (1 - epsA) .* deltaAF; % [a x a]
rhoAgeM = epsA .* ageFraction_F + (1 - epsA) .* deltaAM; % [a x a]
rhoRiskM = zeros(age , risk , risk);
rhoRiskF = rhoRiskM;

for i = 3 : age
    rhoRiskF(i , : , :) = squeeze(epsR .* riskFraction_M(i , : , :))...
        + (1 - epsR) .* deltaR; % [a(i) x r x r] + [r x r] -> [a x r x r]
    rhoRiskM(i , : , :) = squeeze(epsR .* riskFraction_F(i , : , :))...
        + (1 - epsR) .* deltaR; % [a(i) x r x r] + [r x r] -> [a x r x r]
end

% Intialize rho matrices for males and females
rhoM = zeros(age, age, risk, risk);
rhoF = rhoM;
for i = 3 : age
    for ii = 3 : age
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
for aa = 3 : age
    for rr = 1 : risk
        if popSum(2 , aa , rr) ~= 0
            mfRatio(3 : age , aa , : , rr) = ...
                max(popSum(1 , 3 : age , :) ./ popSum(2 , aa ,rr) , 0);
        end
    end
end

B = zeros(age , age , risk , risk);
cMale = partnersM; % [age x risk]
cFemale = partnersF; % [age x risk]
for i = 3 : age
    for ii = 3: age
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
theta = 0.5; 
cAdjMale = zeros(age , age , risk , risk);
cAdjFemale = cAdjMale;
for i = 3 : age
    for j = 1 : risk
        cAdjMale(i , 3 : age , j , :) = cMale(i , j) ...
            .* rhoM(i , 3 : age , j , :)...
            .* B(i , 3 : age , j , :) .^ -(1 - theta) ...
            .* mfRatio(i , 3 : age , j , 1 : risk) .^ theta;
        
        cAdjFemale(3 : age , i , : , j) = cFemale(3 : age , :) ...
            .* squeeze(rhoF(3 : age , i , : , j)) ...
            .* squeeze(B(i , 3 : age , j , :)) .^ theta ...
            .* squeeze(mfRatio(i , 3 : age , j , 1 : risk)) .^ -(1 - theta);
    end
end
cAdj(1 , : , : , : , :) = cAdjMale;
cAdj(2 , : , : , : , :) = cAdjFemale;
cAdj(isnan(cAdj)) = 0;
cAdj(isinf(cAdj)) = 0;

% Protection due to condom usage and circumcision
peakYear = 2000;
condStart = 1995;
yrVec = condStart : 1 / stepsPerYear : peakYear;
condUseVec = linspace(0 , 0.5 * 0.5 , (peakYear - condStart) * stepsPerYear);
condUse = condUseVec(1); % year >= peakYear
if year < peakYear && year > condStart
    yrInd = year == yrVec;
    condUse = condUseVec(yrInd);
elseif year >= peakYear
    condUse = condUseVec(end);
end
cond = 1-(condProtect * condUse); % condom usage and condom protection rates
psi = ones(disease) .* cond;  % vector for protective factors. Scaled to reflect protection by contraception. Currently parameterized for HIV only.
psi(7) = (1 - circProtect) .* cond;
psi(8) = (1 - circProtect) * (1 - prepProtect) .* cond; % no one on PrEP for now
dPop = zeros(size(pop));
newHiv = zeros(gender , age , risk); % incidence tally by gender

%% HPV Infection
% HPV parameters
% for 4V analysis
% beta_hrHPV_F2M = 1 - (1 - perPartnerHpv) ^ 12; % per year per partner probability
% beta_hrHPV_M2F = 1 - (1 - perPartnerHpv) ^ 12;
% 
% beta_lrHPV_F2M = 1 - (1 - perPartnerHpv_lr) ^ 12;
% beta_lrHPV_M2F = 1 - (1 - perPartnerHpv_lr) ^ 12;
% 
% beta_nonV_HPV_F2M = 1 - (1 - perPartnerHpv_nonV) ^ 12;
% beta_nonV_HPV_M2F = 1 - (1 - perPartnerHpv_nonV) ^ 12;

% Per partnership beta
% beta_hrHPV(2 , 1 : 3) = beta_hrHPV_F2M;  % HPV(-) males [g x r]
% beta_hrHPV(1 , 1 : 3) = beta_hrHPV_M2F; % HPV(-) females [g x r]
% beta_lrHPV(2 , 1 : 3) = beta_lrHPV_F2M; % HPV(-) males [g x r]
% beta_lrHPV(1 , 1 : 3) = beta_lrHPV_M2F; % HPV(-) females [g x r]
% beta_nonV_HPV(2 , 1 : 3) = beta_nonV_HPV_F2M; % HPV(-) males [g x r]
% beta_nonV_HPV(1 , 1 : 3) = beta_nonV_HPV_M2F; % HPV(-) females [g x r]
sexPop = zeros(2 , 1);
sexPop(1) = sumall(popSum(1 , 3 : age , :)); % total sexually active males in population
sexPop(2) = sumall(popSum(2 , 3 : age , :)); % total sexually active females in population
beta = zeros(gender , age , risk , 3);
for g = 1 : gender
    for a = 3 : age
        for r = 1 : risk
            beta_hrHPV_F2M = 1 - (1 - perPartnerHpv) ^ maleActs(a , r); % per year per partner probability
            beta_hrHPV_M2F = 1 - (1 - perPartnerHpv) ^ femaleActs(a , r);
            beta_lrHPV_F2M = 1 - (1 - perPartnerHpv_lr) ^ maleActs(a , r);
            beta_lrHPV_M2F = 1 - (1 - perPartnerHpv_lr) ^ femaleActs(a ,r);
            beta_nonV_HPV_F2M = 1 - (1 - perPartnerHpv_nonV) ^ maleActs(a , r);
            beta_nonV_HPV_M2F = 1 - (1 - perPartnerHpv_nonV) ^ femaleActs(a , r);
            % Per partnership beta
            beta_hrHPV(2 , a , r) = beta_hrHPV_F2M;  % HPV(-) males [g x r]
            beta_hrHPV(1 , a , r) = beta_hrHPV_M2F; % HPV(-) females [g x r]
            beta_lrHPV(2 , a , r) = beta_lrHPV_F2M; % HPV(-) males [g x r]
            beta_lrHPV(1 , a , r) = beta_lrHPV_M2F; % HPV(-) females [g x r]
            beta_nonV_HPV(2 , a , r) = beta_nonV_HPV_F2M; % HPV(-) males [g x r]
            beta_nonV_HPV(1 , a , r) = beta_nonV_HPV_M2F; % HPV(-) females [g x r]
            
            hrSum = 0;
            lrSum = 0;
            hrlrSum = 0;
            hr = hrInds(g , a , r , :);
            lr = lrInds(g , a , r , :);
            hrlr = hrlrInds(g , a , r , :);
            hrSum = hrSum + sumall(pop(hr));
            lrSum = lrSum + sumall(pop(lr));
            hrlrSum = hrlrSum + sumall(pop(hrlr));
            p_hrHPV = max(hrSum / popSum(g , a , r) , 0); % proportion of sexually active with 9vHPV. Max handles divide by 0 error.
            p_lrHPV = max(lrSum / popSum(g , a , r) , 0); % proportion of sexually active with 4vHPV
            p_hrlrHPV = max(hrlrSum / popSum(g , a , r) , 0);  % proportion of sexually active with non-vax type HPV
            p_hrHPV = p_hrHPV + p_hrlrHPV;  % total proportion with hr HPV
            p_lrHPV = p_lrHPV + p_hrlrHPV; % total proportion with lr HPV
            % adjust betas for HPV transmission to account for proportion of population
            % that is carrying HPV. Transmission probability throughout population
            % is dependent on the "concentration" of HPV carriers in the population.
            beta(g , a , r , 1) = -log(1 - beta_hrHPV(g , a , r)) .* p_hrHPV;
            beta(g , a , r , 2) = -log(1 - beta_lrHPV(g , a , r)) .* p_lrHPV;
            beta(g , a , r , 3) = -log(1 - beta_nonV_HPV(g , a , r)) .* p_hrlrHPV;
%             beta(g , a , r , 1) = -log(1 - beta_hrHPV(g , r)) .* p_hrHPV;
%             beta(g , a , r , 2) = -log(1 - beta_lrHPV(g , r)) .* p_lrHPV;
%             beta(g , a , r , 3) = -log(1 - beta_nonV_HPV(g , r)) .* p_hrlrHPV;
        end
    end
end

states = 3; % (3 HPV infection states)
% lambda
lambda = zeros(gender, age , risk , states);
for g = 1 : gender
    gg = 2; % partner's gender
    if g == 2
        gg = 1;
    end
    for i = 3 : age
        for j = 1 : risk
            for ii = 3 : age
                for jj = 1 : risk
                    for s = 1 : states
                        lambda(g , i , j , s) = ...
                            lambda(g , i , j , s)...
                            + cAdj(g , i , ii , j , jj)...
                            * beta(gg , ii , jj , s); % [3 x 1]
                    end
                end
            end
        end
    end
end
newHpv = zeros(gender , disease , age , risk);
newImmHpv = newHpv;
newVaxHpv = newHpv;
for d = 1 : disease
    for h = 1% : hpvTypes - 1 % coinfected compartments cannot acquire more infections
        for toState = 1 % only 1 HPV type % CJB note: used to be 1:states(3)
            hTo = toState + 1; 
            for a = 3 : age
                for r = 1 : risk % move age and risk iterators below fromState and toState to reduce unneccessary iterations
                    if lambda(1 , a , r , toState) > 10 ^ -6 || lambda(2 , a , r , toState) > 10 ^ -6 ... % only evaluate if lambda is non-zero
                            && hTo ~= h % and if not acquiring pre-existing infection
                        
                        % Prepare indices
                        % non-vaccinated, non-naturally immune, screened or non-screened
                        mSus = hpvSus(d , h , 1 , a , r , :);
                        fSus = hpvSus(d , h , 2 , a , r , :);
                        mTo = toHpv(d , hTo , 1 , a , r , :); % update HPV state to infected if just acquiring HPV
                        fTo = toHpv(d , hTo , 2 , a , r , :);

                        % non-vaccinated, naturally immune, screened or non-screened
                        %mSusImm = hpvImm(d , hTo , 1 , a , r , :); % only females have natural immunity
                        fSusImm = hpvImm(d , hTo , 2 , a , r , :);
                        %mToImm = toHpv_Imm(d , hTo , 1 , a , r , :);
                        fToImm = toHpv_Imm(d , hTo , 2 , a , r , :); % update HPV state to infected if just acquiring HPV
                        
                        % vaccinated w/ no infection history, non-naturally immune , non-screened
                        mSusVax = hpvVaxd(d , h , 1 , a , r , :); % vaccinated with no infection history
                        fSusVax = hpvVaxd(d , h , 2 , a , r , :);
                        mToVax = toHpv_Vax(d , hTo , 1 , a , r , :); % vaccinated and getting infected for the first time
                        fToVax = toHpv_Vax(d , hTo , 2 , a , r , :);
                        
                        % vaccinated w/ no infection history, non-naturally immune , screened
                        mSusVaxScreen = hpvVaxdScreen(d , h , 1 , a , r , :); % vaccinated + screened with no infection history
                        fSusVaxScreen = hpvVaxdScreen(d , h , 2 , a , r , :);
                        mToVaxScreen = toHpv_VaxScreen(d , hTo , 1 , a , r , :); % vaccinated + screened and getting infected for the first time
                        fToVaxScreen = toHpv_VaxScreen(d , hTo , 2 , a , r , :);
                        
                        % vaccinated w/ infection history, non-naturally immune, screened or non-screened
                        mSusVax2 = hpvVaxd2(d , h , 1 , a , r , :); % vaccinated with vaccine-type infection history
                        fSusVax2 = hpvVaxd2(d , h , 2 , a , r , :);
                        mToVax2 = hpvVaxd2(d , hTo , 1 , a , r , :); % vaccinated and getting infected by non-vaccine type
                        fToVax2 = hpvVaxd2(d , hTo , 2 , a , r , :);
                        
                        % vaccinated w/ infection history, naturally immune, screened or non-screened
                        fSusImmVax2 = hpvImmVaxd2(d , h , 2 , a , r , :);
                        
                        % for non-vaccine type infections (not currently used, instead include all hr HPV infections as "one type" and reduce vaccine efficacy to 90%)
                        mSusVax3 = toHpv_VaxNonV(d , h , 1 , a , r , :); % vaccinated with non-vaccine type infection history (not currently used)
                        fSusVax3 = toHpv_VaxNonV(d , h , 2 , a , r , :); % (not currently used)
                        mToVax3 = toHpv_VaxNonV(d , hTo , 1 , a , r , :); % vaccinated and getting infected by non-vaccine type, no vaccine-type infection history (not currently used)
                        fToVax3 = toHpv_VaxNonV(d , hTo , 2 , a , r , :); % (not currently used)
                        % Note: would need code to deal with toHpv_VaxNonVScreen if we add back in separate non-vaccine type infections
                        
                        
                        % Set lambda multipliers based on CD4 count
                        lambdaMultF = 1;
                        lambdaMultM = 1;
                        if d > 2 && d < 7% && toState < 3 % CD4 > 500 -> CD4 < 200
                            lambdaMultF = hpv_hivMult(d - 2 , hTo - 1);
                            lambdaMultM = hpv_hivMult(d - 2 , hTo - 1);
                        elseif d == 10
                            lambdaMultF = artHpvMult; 
                            lambdaMultM = artHpvMult;
                        end
 
                        
                        % Calculate infections
                        % normal susceptibles (non-vaccinated, non-naturally immune, screened or non-screened)
                        % infection rate capped at 0.99
                        mInfected = min(lambdaMultM * lambda(1 , a , r , toState)...
                            , 0.999) .* pop(mSus); % infected males
                        fInfected = min(lambdaMultF * lambda(2 , a , r , toState)...
                            , 0.999) .* pop(fSus); % infected females

                        % naturally immune (non-vaccinated, naturally immune, screened or non-screened)
                        %mInfImm = min(lambdaMultM * lambdaMultImm(a) * lambda(1 , a , r , toState) ...
                        %    , 0.999) .* pop(mSusImm);
                        fInfImm = min(lambdaMultF * lambdaMultImm(a) * lambda(2 , a , r , toState) ...
                            , 0.999) .* pop(fSusImm);
                        
                        % vaccinated susceptibles
                        hVax = 1;
                        if hTo == 3
                            hVax = 2;
                        end
                        vaxProtect = max(lambdaMultVax(a , hVax) , (hTo == 4) .* 1); % oHR has no vaccine protection
                        
                        mInfVax = min(lambdaMultM * vaxProtect * lambda(1 , a , r , toState) ...
                            , 0.999 * vaxProtect) .* pop(mSusVax); % vaccinated w/ no infection history, non-naturally immune , non-screened
                        fInfVax = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                            , 0.999 * vaxProtect) .* pop(fSusVax);
                        
                        mInfVaxScreen = min(lambdaMultM * vaxProtect * lambda(1 , a , r , toState) ...
                            , 0.999 * vaxProtect) .* pop(mSusVaxScreen); % vaccinated w/ no infection history, non-naturally immune , screened
                        fInfVaxScreen = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                            , 0.999 * vaxProtect) .* pop(fSusVaxScreen);
                                      
                        mInfVax2 = min(lambdaMultM * vaxProtect * lambda(1 , a , r , toState) ...
                            , 0.999 * vaxProtect) .* pop(mSusVax2); % vaccinated w/ infection history, non-naturally immune, screened or non-screened
                        fInfVax2 = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                            , 0.999 * vaxProtect) .* pop(fSusVax2);
                        
                        fInfImmVax2 = min(lambdaMultF * vaxProtect * lambdaMultImm(a) * lambda(2 , a , r , toState) ...
                            , 0.999 * vaxProtect) .* pop(fSusImmVax2); % vaccinated w/ infection history, naturally immune, screened or non-screened
                        
                        mInfVax3 = min(lambdaMultM * vaxProtect * lambda(1 , a , r , toState) ...
                            , 0.999 * vaxProtect) .* pop(mSusVax3);  % (currently not used = 0)
                        fInfVax3 = min(lambdaMultF * vaxProtect * lambda(2 , a , r , toState) ...
                            , 0.999 * vaxProtect) .* pop(fSusVax3); % (currently not used = 0)
                        
                        
                        % Incidence tracker
                        % normal susceptibles (non-vaccinated, non-naturally immune, screened or non-screened)
                        newHpv(1 , d , a , r) = newHpv(1 , d , a , r) + sumall(mInfected);
                        newHpv(2 , d , a , r) = newHpv(2 , d , a , r) + sumall(fInfected);

                        % naturally immune (non-vaccinated, naturally immune, screened or non-screened)
                        %newImmHpv(1 , d , a , r) = newImmHpv(1 , d , a , r) + sumall(mInfImm);
                        newImmHpv(2 , d , a , r) = newImmHpv(2 , d , a , r) + sumall(fInfImm);

                        % vaccinated
                        newVaxHpv(1 , d , a , r) = newVaxHpv(1 , d , a , r) ...
                            + sumall(mInfVax) + sumall(mInfVaxScreen) + sumall(mInfVax2) + sumall(mInfVax3);
                        newVaxHpv(2 , d , a , r) = newVaxHpv(2 , d , a , r)...
                            + sumall(fInfVax) + sumall(fInfVaxScreen) + sumall(fInfVax2) + sumall(fInfImmVax2) + sumall(fInfVax3);

                        
                        % Adjust compartments
                        % normal susceptibles (non-vaccinated, non-naturally immune, screened or non-screened)
                        dPop(mSus) = dPop(mSus) - mInfected; % efflux of infected males
                        dPop(fSus) = dPop(fSus) - fInfected; % efflux of infected females
                        dPop(mTo) = dPop(mTo) + mInfected; % influx of infected males
                        dPop(fTo) = dPop(fTo) + fInfected; % influx of infected females

                        % naturally immune (non-vaccinated, naturally immune, screened or non-screened)
                        %dPop(mSusImm) = dPop(mSusImm) - mInfImm;
                        dPop(fSusImm) = dPop(fSusImm) - fInfImm;
                        %dPop(mToImm) = dPop(mToImm) + mInfImm;
                        dPop(fToImm) = dPop(fToImm) + fInfImm;

                        %vaccinated
                        dPop(mSusVax) = dPop(mSusVax) - mInfVax; % vaccinated w/ no infection history, non-naturally immune , non-screened
                        dPop(fSusVax) = dPop(fSusVax) - fInfVax;
                        
                        dPop(mSusVaxScreen) = dPop(mSusVaxScreen) - mInfVaxScreen; % vaccinated w/ no infection history, non-naturally immune , screened
                        dPop(fSusVaxScreen) = dPop(fSusVaxScreen) - fInfVaxScreen;
                        
                        dPop(mSusVax2) = dPop(mSusVax2) - mInfVax2; % vaccinated w/ infection history, non-naturally immune, screened or non-screened
                        dPop(fSusVax2) = dPop(fSusVax2) - fInfVax2;
                        
                        dPop(fSusImmVax2) = dPop(fSusImmVax2) - fInfImmVax2; % vaccinated w/ infection history, naturally immune, screened or non-screened
                        
                        dPop(mSusVax3) = dPop(mSusVax3) - mInfVax3; % (currently not used = 0)
                        dPop(fSusVax3) = dPop(fSusVax3) - fInfVax3; % (currently not used = 0)
                        
                        if hTo ~= 4 % vaccinated and acquiring vaccine type
                            dPop(mToVax) = dPop(mToVax) + mInfVax + mInfVax3;
                            dPop(fToVax) = dPop(fToVax) + fInfVax + fInfVax3;
                            
                            dPop(mToVaxScreen) = dPop(mToVaxScreen) + mInfVaxScreen;
                            dPop(fToVaxScreen) = dPop(fToVaxScreen) + fInfVaxScreen;
                            
                            dPop(mToVax2) = dPop(mToVax2) + mInfVax2;
                            dPop(fToVax2) = dPop(fToVax2) + fInfVax2 + fInfImmVax2;
                        else % vaccinated and acquiring non-vaccine type (currently not used)
                            dPop(mToVax3) = dPop(mToVax3) + mInfVax + mInfVax3; % male: no previous vaccine-type infection history
                            dPop(fToVax3) = dPop(fToVax3) + fInfVax + fInfVax3; % female: no previous vaccine-type infection history
                            
                            dPop(mToVax) = dPop(mToVax) + mInfVax2; % male: previous vaccine-type infection history
                            dPop(fToVax) = dPop(fToVax) + fInfVax2; % female: previous vaccine-type infection history
                        end
                    end
                end
            end
        end
    end
end

newInfs{1} = newHpv;
newInfs{2} = newImmHpv;
newInfs{3} = newVaxHpv;

%% HIV Infection
% HIV average betas
beta = zeros(risk , age , gender , age , risk);

% infection probability by viral load
for a = 1 : age
%     for aa = 1 : age
        for r = 1 : risk
%             for rr = 1 : risk
                if popSum(1 , a , r) ~= 0
                    for v = 1 : 5 % viral load (up to vl = 6). Note: last index is (viral - 1) + 1. Done to align pop index with betaHIV index.
                        beta(: , : , 1 , a , r) = beta(: , : , 1 , a , r) - log(1 - betaHIVM2F(: , : , v)) ...
                            * sumall(pop(mCurr(a , r , v , :))) ./ popSum(1 , a , r);
                    end
                    beta(: , : , 1 , a , r) = beta(: , : , 1 , a , r) - log(1 - betaHIVM2F(: , : , 6)) ...
                        * sumall(pop(mCurrArt(a , r , 1 , :)))  ./ popSum(1 , a , r);
                end
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
    for i = 3 : age
        for j = 1 : risk
            for ii = 3 : age
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
for a = 3 : age
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