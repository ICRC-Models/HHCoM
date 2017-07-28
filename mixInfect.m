% Population Mixing and Infection Module
% Simulates mixing and HPV/HIV coinfection in a population.
% Models heterosexual interactions and corresponding transmissions.
% Accepts contact matrix [agecohort x sact x gender] and a population matrix
% as input and returns dPop, a matrix of derivatives that describes the
% change in the population's subgroups.
function [dPop , newInfs] = mixInfect(t , pop , betaOption , currStep , ...
    naive , coInf, gar , hivSus , hpvSus , toHiv , toHpv , mCurr , fCurr , ...
    mCurrArt , fCurrArt ,epsA_vec , epsR_vec , yr , modelYr1 , ...
    circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
    hpv_hivMult , betaHIVF2M , betaHIVM2F , beta_hrHPV_val , beta_lrHPV_val , disease , ...
    viral , gender , age , risk , hpvStates , hpvTypes , hrInds , lrInds ,...
    hrlrInds,  k , periods , stepsPerYear , year)
toInd = @(x) (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + (x(: , 6) - 1) * k(5) ...
    + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + (x(: , 3) - 1) * k(2) ...
    + (x(: , 2) - 1) * k(1) + x(: , 1);
sumall = @(x) sum(x(:));
% vlBeta = 1; % evaluate beta using VL or CD4
% if strcmp('cd4' , betaOption)
%     vlBeta = 0; % Will be used to determine beta calcuation method, i.e. by CD4 or VL
% end
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
deltaAF = eye(12) .* 0.3 + diag(ones(11 , 1) .* 0.7 , 1);
% deltaAF(12 , 11) = 0.7;
deltaAM = eye(12) .* 0.3 + diag(ones(11 , 1) .* 0.7 , -1);
% deltaAM(12 , 11) = 0.7;
deltaAM(1 , 2) = 0.7;
% deltaA = eye(12) .* 0.5 + diag(ones(11 , 1) .* 0.3 , 1) + diag(ones(11 , 1) .* 0.2 , -1); % orig 0.7, 0.3
% after 2005
if currStep > (2005 - modelYr1) * stepsPerYear
    deltaAF = eye(12) .* 0.7 + diag(ones(11 , 1) .* 0.3 , 1);
%     deltaAF(12 , 11) = 0.3;
    deltaAM = eye(12) .* 0.7 + diag(ones(11 , 1) .* 0.3 , -1);
%     deltaAM(1 , 2) = 0.3;
    deltaR = eye(3);
end

acts = actsPer; % acts per partnership, from loaded workspace [gender x risk]
% rate of partner change (contact)
% males
c(1 , : , :) = partnersM;
% females
c(2 , : , :) = partnersF;
% HIV parameters
% betaHIV_F = betaHIV;
% betaHIV_M = betaHIV;

% acts = [90 30 3; 70 20 3] .* 1.1; %Braunstein, Johnson; dimensions = gender x sact
% betaHIV = zeros(gender, risk, viral);
% per year, per partnership HIV beta
% betaHIVM = 1 - (bsxfun(@power, 1 - betaHIV_M , acts(1 , :)')); % HIV(-) males
% betaHIVF = 1 - (bsxfun(@power, 1 - betaHIV_F , acts(2 , :)')); % HIV(-) females

% HPV parameters
beta_hrHPV_F2M = 0.718;% 0.143;%beta_hrHPV_val; % from hpvData.m
beta_hrHPV_M2F = 0.718;%0.143;%beta_hrHPV_val;

beta_lrHPV_F2M = 0.718;%0.143;% beta_lrHPV_val; % from hpvData.m
beta_lrHPV_M2F = 0.718;%0.143;%beta_lrHPV_val;

% per year, per partnership betaHPV for high risk and low risk HPV. int = periods per
% year. Act values are per year.
% beta_hrHPV(1 , :) = 1 - (1 - beta_hrHPV_F2M) .^ (acts(1 , :)); % HPV(-) males [g x r]
% beta_hrHPV(2 , :) = 1 - (1 - beta_hrHPV_M2F) .^ (acts(2 , :)); % HPV(-) females [g x r]
% beta_lrHPV(1 , :) = 1 - (1 - beta_lrHPV_F2M) .^ (acts(1 , :)); % HPV(-) males [g x r]
% beta_lrHPV(2 , :) = 1 - (1 - beta_lrHPV_M2F) .^ (acts(2 , :)); % HPV(-) females [g x r]

% Change to per partnership beta
beta_hrHPV(1 , 1 : 3) = beta_hrHPV_F2M;  % HPV(-) males [g x r]
beta_hrHPV(2 , 1 : 3) = beta_hrHPV_M2F; % HPV(-) females [g x r]
beta_lrHPV(1 , 1 : 3) = beta_lrHPV_F2M; % HPV(-) males [g x r]
beta_lrHPV(2 , 1 : 3) = beta_lrHPV_M2F; % HPV(-) females [g x r]




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
%             ind = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r));
%             popSum(g , a , r) = sumall(pop(ind));
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
%                 B(i, ii , j , jj) = min(cMale(i , j) * rhoM(i , ii , j , jj) ...
%                     /(cFemale(ii , jj) * rhoF(ii , i , jj , j)) ...
%                     * mfRatio(i , ii , j , jj) , min(cMale(i , j) ...
%                     / cFemale(ii , jj) , max(cMale(i , j) , cFemale(ii, jj))));
            end
        end
    end
end
B(isnan(B)) = 0; % 0/0 , 1/0 -> 0
B(isinf(B)) = 0;
% adjust c to account for discrepancy in reporting
theta = 0.5; % assume adjusted contact rate is equally male and female driven
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
x = zeros(gender , age , risk);

% if vlBeta
% evaluate beta by viral load
for a = 1 : age
    for r = 1 : risk
        if popSum(1 , a , r) ~= 0
            for v = 1 : 5 % viral load (up to vl = 6). Note: last index is (viral - 1) + 1. Done to align pop v index with betaHIV v index.
                x(1 , a , r) = x(1 , a , r) + sumall(pop(mCurr(a , r , v , :))) * ...
                    betaHIVM2F(a , r , v) ./ popSum(1 , a , r);
            end
            x(1 , a , r) = x(1 , a , r) + sumall(pop(mCurrArt(a , r , 1 , :))) * ...
                betaHIVM2F(a , r , 6) ./ popSum(1 , a , r);
        end
        if popSum(2 , a , r) ~= 0
            for v = 1 : 5 % viral load (up to vl = 6). Note: last index is (viral - 1) + 1. Done to align pop v index with betaHIV v index.
                x(2 , a , r) = x(2 , a , r) + sumall(pop(fCurr(a , r , v , :))) * ...
                    betaHIVF2M(a , r , v) ./ popSum(2 , a , r);
            end
            x(2 , a , r) = x(2 , a , r) + sumall(pop(fCurrArt(a , r , 1 , :))) * ...
                betaHIVF2M(a , r , 6) ./ popSum(2 , a , r);
        end
    end
end
% end

% average betas for "high risk" and "low risk HPV"
y = repmat(beta_hrHPV , [1 , 1 , 12]); % replicate along age dimension -> [g x r x a]
z = repmat(beta_lrHPV , [1 , 1 , 12]);
z = permute(z , [1 3 2]); % [g x r x a] -> [g x a x r]
y = permute(y , [1 3 2]); % [g x r x a] -> [g x a x r]



sexPop = zeros(2 , 1);
sexPop(1) = sumall(popSum(1 , 3 : age , :)); % total sexually active males in population
sexPop(2) = sumall(popSum(2 , 3 : age , :)); % total sexually active females in population

for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            hrSum = 0;
            lrSum = 0;
            hrlrSum = 0;
            hr = hrInds(g , a , r , :);
            lr = lrInds(g , a , r , :);
            hrlr = hrlrInds(g , a , r , :);
            hrSum = hrSum + sumall(pop(hr));
            lrSum = lrSum + sumall(pop(lr));
            hrlrSum = hrlrSum + sumall(pop(hrlr));
            p_hrHPV = max(hrSum / popSum(g , a , r) , 0); % proportion of sexually active with hrHPV. Max handles divide by 0 error.
            p_lrHPV = max(lrSum / popSum(g , a , r) , 0); % proportion of sexually active with lrHPV
            p_hrlrHPV = max(hrlrSum / popSum(g , a , r) , 0);  % proportion of sexually active with hrHPV and lrHPV
            p_hrHPV = p_hrHPV + p_hrlrHPV;  % total proportion with hr HPV
            p_lrHPV = p_lrHPV + p_hrlrHPV; % total proportion with lr HPV
            % adjust betas for HPV transmission to account for proportion of population
            % that is carrying HPV. Transmission probability throughout population
            % is dependent on the "concentration" of HPV carriers in the population.
            y(g , a , r) = y(g , a , r) .* p_hrHPV;
            z(g , a , r) = z(g , a , r) .* p_lrHPV;
        end
    end
end

% coinfection transmission matrix
states = 8; % (3 single infection states, 3 dual infection states, 1 triple infection state, 1 completely susceptible state)
a = z .* y .* x; % probability of acquiring all infections in one time period
% beta dim [states x states x gender x age x risk]
xz = z .* x - a;
xy = y .* x - a;
yz = z .* y - a;
zero = zeros(gender , age , risk);
% Susceptible to all
beta(: , : , : , : , 1) = cat(4 , zero , zero , zero , zero , zero , zero , zero , zero);
%HIV positive
beta(: , : , : , : , 2) = cat(4 , y - y .* (x + z) + a ,  zero , zero , zero ,...
    zero , zero , zero , zero);
beta(: , : , : , : , 3) = cat(4 , z - z .* (x + y) + a ,  zero , zero , zero , ...
    zero , zero , zero , zero);
beta(: , : , : , : , 4) = cat(4 , yz , z , y , zero , zero , zero , zero , zero);  % probability of acquiring both y and z given current state
% HIV negative
beta(: , : , : , : , 5) = cat(4 , x - x .* (y + z) + a,  zero , zero , zero ,...
    zero , zero , zero , zero);
beta(: , : , : , : , 6) = cat(4 , xy , x , zero , zero , y , zero , zero , zero); % probability of acquiring both x and y given current state
beta(: , : , : , : , 7) = cat(4 , xz , zero , x , zero , z , zero , zero , zero); % probability of acquiring both x and z given current state
beta(: , : , : , : , 8) = cat(4 , a , x .* z , x .* y , x , z .* y , z , y , zero); % probability of acquiring x , y , and z given current state

% lambda
lambda = zeros(gender, age , risk , states , states);
for g = 1 : gender
    gg = 2; % partner's gender
    if g == 2
        gg = 1;
    end
    for i = 1 : age
        for j = 1 : risk
            for ii = 1 : age
                for jj = 1 : risk
                    lambda(g , i , j , : , :) = ...
                        squeeze(lambda(g , i , j , : , :))...
                        + cAdj(g , i , ii , j , jj) * rho(g , i , ii , j , jj)...
                        * squeeze(beta(gg , ii , jj , : , :)); % [2 x 1]
                    % when sum of lambda for all possible transitions > 1, adjust
                    % probabilities to reflect relative proportions of subpopulation
                    % that will move to different compartments or remain in same compartment.
%                     lambdaCurr = squeeze(lambda(g , i , j , : , :));
                end
%                 for s = 1 : states
%                     pTransTotal = sumall(lambdaCurr(s , :));
%                     if pTransTotal >= 1
%                         lambda(g , i , j , s , :) = ...
%                             lambda(g , i , j , s , :) ./ pTransTotal * 0.99;
%                     end
%                 end
            end
        end
    end
end
%disp(['Max lambda = ' , num2str(max(lambda(:)))])
cond = 1-(condProtect * condUse); % condom usage and condom protection rates
psi = ones(disease) .* cond;  % vector for protective factors. Scaled to reflect protection by contraception. Currently parameterized for HIV only.
psi(7) = 1 - circProtect .* cond;
psi(8) = (1 - circProtect) * (1 - prepProtect) .* cond;
%psi(10) = 1 - artProtect .* cond;
dPop = zeros(size(pop));
newHiv = zeros(gender , age , risk); % incidence tally by gender
newHpv = zeros(gender , age , risk);
% drawnow
% hist(lambda(lambda(:)>0))
% title('Lambda'); xlabel('Value'); ylabel('Frequency')
%% Infection
for d = 1 : disease
    hiv = d > 1 && d < 7 || d == 10;
    for h = 1 : hpvTypes
        hpv = h > 1;
        fromState = h + 4 * hiv;
        for toState = 2 + 3 * hiv : 8 % exclude state 1 (completely susceptible)
            for a = 1 : age
                for r = 1 : risk % move age and risk iterators below fromState and toState to reduce unneccessary iterations
                    if lambda(1 , a , r , fromState , toState)...
                            || lambda(2 , a , r , fromState , toState) % only evaluate if lambda is non-zero
                        hTo = mod(toState , 4); % index for HPV infection/coinfection
                        hTo(hTo == 0) = 4; % HPV coinfection
                        hivInfection = toState > 4 && not(hiv); % if HIV(-) experiences HIV infection
                        hpvInfection = not(hpv) && hTo ~= 1; % if no HPV and no hysterectomy and acquiring HPV
                        % intialize "from compartment(s)" vector(s) to
                        % completely naive compartments
                        mSus = naive(1 , a , r , :);
                        fSus = naive(2 , a , r , :);
                        % intialize "to compartment(s)" vector(s) to
                        % coinfected compartments
                        mTo = coInf(hTo , 1 , a , r , :);
                        fTo = coInf(hTo , 2 , a , r , :);
                        if not(hivInfection && hpvInfection) % no coinfection
                            if hivInfection
                                mSus = hivSus(d , h , 1 , a , r , :);
                                fSus = hivSus(d , h , 2 , a , r , :);
                                mTo = toHiv(hTo , 1 , a , r , :); % set disease state = 2 (acute)
                                fTo = toHiv(hTo , 2 , a , r , :); % and update viral load state to reflect HIV acquisition
                            else % hpvInfection
                                mSus = hpvSus(d , 1 , a , r , :);
                                fSus = hpvSus(d , 2 , a , r , :);
                                mTo = toHpv(d , hTo , 1 , a , r , :); % update HPV state to infected if just acquiring HPV
                                fTo = toHpv(d , hTo , 2 , a , r , :);
                            end
                        end
%                         assert(lambda(1 , a , r , fromState , toState) < 1 , ...
%                             ['lambda for 1 , ' , num2str(a) , ', ' num2str(r) ' >= 1'])
%                         assert(lambda(2 , a , r , fromState , toState) < 1 , ...
%                             ['lambda for 2 , ' , num2str(a) , ', ' num2str(r) ' >= 1'])
                        lambdaMult = 1;
                        if hpv && d > 2 && d < 7 % hpv infection & CD4 > 500 -> CD4 < 200
                            lambdaMult = hpv_hivMult(d - 2 , h - 1);
                        end
                        mInfected = lambdaMult * lambda(1 , a , r , fromState , toState)...
                            .* psi(d) .* pop(mSus); % infected males
                        fInfected = lambdaMult * lambda(2 , a , r , fromState , toState)...
                            .* psi(d) .* pop(fSus); % infected females
                        if hivInfection
                            newHiv(1 , a , r) = newHiv(1 , a , r) + sumall(mInfected);
                            newHiv(2 , a , r) = newHiv(2 , a , r) + sumall(fInfected);
                        else
                            newHpv(1 , a , r) = newHpv(1 , a , r) + sumall(mInfected);
                            newHpv(2 , a , r) = newHpv(2 , a , r) + sumall(fInfected);
                        end
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
newInfs{1} = newHiv;
newInfs{2} = newHpv;

% for testing. Remove later.
% if year == 2000
%     save('partnershipMats' , 'cAdj' , 'B');
% end
% Convert to column vector for output to ODE solver
dPop = sparse(dPop);
