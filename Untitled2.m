% ccDeath dims: [t , disease , viral , hpvTypes , age]

annualDiscRate = 0.03; % make input later

%% Calculate discount rate
discountRate = (1 + annualDiscRate) ^ (1 / stepsPerYear) - 1; % annual rate to rate corresponding to step size
discountVec = (1 + discountRate) .^ (-(0 : length(tVec) - 1));

%% Find years of life lost
YLL_normal = YLL_normal + deaths - lifeExp(a); % years of life lost to natural causes
YLL_hiv = YLL_hiv + hivDeaths - lifeExp(a); % years of life lost to HIV
YLL_CC = YLL_CC + sum(sum(sum(ccDeath(: , : , : , : , a) , 2) , 3) , 4) - lifeExp(a); % years of life lost to cervical cancer

YLL = YLL_normal + YLL_hiv + YLL_CC; % total years of life lost by age and time step
YLL_Disc = sum(sum(YLL .* discountVec)); % discounted years of life lost



