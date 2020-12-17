% The expected value of a beta distribution random variable mu(X) = a/(a+b).
% When a = b, the mean = 0.5, which means the expected value is in the middle of the distribution. 
% The variance: var(X) = ab/[((a+b)^2) * (a + b + 1)]

% Parameters for HIV transmission probability per sex act 
X = 0 : 0.001 : 0.25;
a_x = 1.0; 
b_x = 600;

distShape = betapdf(X, a_x, b_x);

figure
plot(X, distShape)

mu_x = a_x/(a_x+b_x)
var_x = (a_x*b_x)/(((a_x+b_x)^2) * (a_x + b_x + 1))

% Gray, 2001; Lancet. Unadjusted HIV transmission prob per act: 0·0011 (95%CI 0·0008–0·0015)

%% HPV transmission probability
Y = 0 : 0.001 : 0.25;
a_y = 3; 
b_y = 200;

distShape = betapdf(Y, a_y, b_y);

figure
plot(Y, distShape)

mu_y = a_y/(a_y+b_y)
var_y = (a_y*b_y)/(((a_y+b_y)^2) * (a_y + b_y + 1))

%% List of 100 random selected HIV, HPV, and HIV-HPV interaction paramer sets

file = [pwd , '\hiv_hpv_probs.xlsx']
probs = xlsread(file, 'B1:CW3')
paramDir = [pwd , '\Params\'];
% save(fullfile(paramDir ,'stochasticParamsets'), 'probs' );
dlmwrite([paramDir ,'stochasticParamsets.dat'], probs)