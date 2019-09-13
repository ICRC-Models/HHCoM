% Load parameters and other input data. Accepts steps per year as an
% argument
function[] = loadUp(stepsPerYear)
load('settings')

%% import and save population parameters
file = [pwd , '\Config\Population_data.xlsx'];
disp(['Loading data from ' , file , '...']);
disp('This may take a couple seconds...');
% demographic data
popInit = xlsread(file , 'Demographics' , 'B6:C21'); % [age x gender]; initial population size
riskDistM = xlsread(file , 'Demographics' , 'B54:D69'); % [age x risk]; male risk distribution
riskDistF = xlsread(file , 'Demographics' , 'E54:G69'); % [age x risk]; female risk distribution
mue = xlsread(file , 'Demographics' , 'B84:C99'); % [age x gender]; background mortality
fertility = xlsread(file , 'Demographics' , 'B112:G127'); % [age x disease]; fertility rate per year before 1995
fertility2 = xlsread(file , 'Demographics' , 'B133:G148'); % [age x disease]; fertility rate/2 for 2005 onwards
partnersM = xlsread(file , 'Demographics' , 'B160:D175'); % [age x risk]; male partnerships per year
partnersF = xlsread(file , 'Demographics' , 'E160:G175'); % [age x risk]; female partnerships per year
actsPer = xlsread(file , 'Demographics' , 'B238:D239'); % [gender x risk]; acts per partnership (not currently used)
epsA = xlsread(file , 'Demographics' , 'B185:B187'); % [year] <1985 , 1990 , >2000; force of infection mixing by age
epsR = xlsread(file , 'Demographics' , 'C185:C187'); % [year] <1985 , 1990 , >2000; force of infection mixing by sexual risk
yr = xlsread(file , 'Demographics' , 'A185:A187'); % years
savdir = [pwd , '\Params']; 
save(fullfile(savdir , 'popData') , 'popInit' , 'riskDistM' , 'riskDistF' , 'mue' , 'fertility' , 'partnersM' , ...
    'partnersF' , 'actsPer' , 'epsA' , 'epsR' , 'fertility2' , 'yr');

%% import and save HIV parameters
file = [pwd , '\Config\HIV_parameters.xlsx'];
disp(['Loading data from ' , file , '...']);
disp('This may take a couple seconds...');

% protection data
circ = xlsread(file , 'Protection' , 'B4 : C4'); % Circumcision rate [1 x year] <2009, 2012>
circProtect = xlsread(file , 'Protection' , 'B18'); % Protection conferred by circumcision
condProtect = xlsread(file , 'Protection' , 'B19'); % Protection conferred by condom use
condUse = xlsread(file , 'Protection' , 'B25'); % Average proportion of population using condoms (not currently used, reset in mixInfect.m)
% disease data
MTCTRate = xlsread(file , 'Disease Data' , 'B6:B8'); % <2004, 2005 , >2008; mother to child transmission rate
mtctVec = linspace(MTCTRate(1) , MTCTRate(end) , size(MTCTRate , 1) * 4);
muHIV = xlsread(file , 'Disease Data' , 'B20 : G35'); %[age x cd4], [12 x 6]; HIV mortality
kCD4(1 , : , :) = xlsread(file , 'Disease Data' , 'B44:E48'); % [gender x vl x cd4]; male cd4 progression by vl
kCD4(2 , : , :) = xlsread(file , 'Disease Data' , 'B52:E56'); % [gender x vl x cd4]; female cd4 progression by vl
% viral load progression by CD4 count
kVl(1 , : , :) = xlsread(file , 'Disease Data' , 'B65 : E69'); % [gender x cd4 x vl]; male vl progression by cd4
kVl(2 , : , :) = xlsread(file , 'Disease Data' , 'B73 : E77'); % [gender x cd4 x vl]; female vl progression by cd4

%% Actual data
% disp('Retrieving actual data')
% overallHivPrev_SA = xlsread('HIVParameters.xlsx' , 'Actual' , 'A22 : B32');
% overallHivPrev_KZN = xlsread('HIVParameters.xlsx' , 'Actual' , 'A33 : B43');
% overallHivPrev_KZN_2 = xlsread('HIVParameters.xlsx' , 'Actual' , 'E33 : H36');
% overallHivPrev_KZN_3 = xlsread('HIVParameters.xlsx' , 'Actual' , 'E22 : H27');
% overallHivPrev_KZN_3 = [overallHivPrev_KZN_3 ; xlsread('HIVParameters.xlsx' , 'Actual' , 'E33 : H36')];
% HivPrevWomen = xlsread('HIVParameters.xlsx' , 'Actual' , 'J5 : Q11'); %  age group (15 - 19 to 45 - 49) x year (2004-2011)
% HivPrevMen = xlsread('HIVParameters.xlsx' , 'Actual' , 'U5 : AB11');
% prevYears = xlsread('HIVParameters.xlsx' , 'Actual' , 'U4 : AB4');
% save('actual')
% disp('Actual data loaded')
% disp(' ')

%%
% Treatment intiation rates
toPrep = 0; % initially

% 2004 - 2006 begin ART
% below200Art_2004 = interp1([2004 , 2006] , [0 0.25] ,[2004 : 1 / stepsPerYear : 2006]); % ART coverage for persons with CD4 < 200
% above200Art_2004 = interp1([2004 , 2006] , [0 0.25] , [2004 : 1 / stepsPerYear : 2006]);
% pie4Vec_2004 = interp1([2004 , 2006] , [0 0.01] , [2004 : 1 / stepsPerYear : 2006]);
% % 2006 - 2013 begin ART
% below200Art_2006 = interp1([2006 , 2013] , [0.03 0.25] , [2006 : 1 / stepsPerYear : 2013]);
% pieVec5_2006 = interp1(linspace(2006 , 2013 , stepsPerYear) , [0.02 0.06]);
% pieVec4_2006 = interp1(linspace(2006 , 2013 , stepsPerYear) , [0.01 0.06]);
%
% if year >= 2014
%     pie(6 , : , : , 4 : end , :) = 0.06; % Home HTC study
%     pie(5 , : , : , 4 : end , :) = 0.02;
%     pie(4 , : , : , 4 : end , :) = 0.01;
% elseif year >= 2006
%     pie(6 , : , : , 4 : end, :) = 0.03; % ART coverage for persons with CD4 < 200
% end

% Treatment dropout rates
prepOut = xlsread(file , 'Disease Data' , 'B85');  %not currently used (reset in hiv2a)
artOut = xlsread(file , 'Disease Data' , 'C85');  %not currently used (reset in hiv2a)

save(fullfile(savdir , 'HIVParams')) % save HIV parameters to a workspace file

%% save general parameters
disease = 10;
viral = 6;
hpvTypes = 4;
hpvStates = 10;
periods = 6;
gender = 2;
age = 80;
risk = 3;

dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age , risk];

% index retrieval function
k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);

toInd = @(x) (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + (x(: , 6) - 1) * k(5) ...
    + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + (x(: , 3) - 1) * k(2) ...
    + (x(: , 2) - 1) * k(1) + x(: , 1);

sumall = @(x) sum(x(:));

savdir = [pwd , '\Params'];
save(fullfile(savdir ,'general'), 'disease' , 'viral' , 'hpvTypes' , 'hpvStates' , 'periods' ,...
    'gender' , 'age' , 'risk' , 'dim' , 'k' , 'sumall' , 'toInd' , 'circ' , ...
    'condUse' , 'kCD4' , 'kVl' , 'stepsPerYear'); % save general model parameters to a workspace file

%% save parameters for mixInfect
load('settings')
step = 1 / stepsPerYear;
epsA_vec = cell(size(yr , 1) - 1, 1); % save data over time interval in a cell array
epsR_vec = cell(size(yr , 1) - 1, 1);
for i = 1 : size(yr , 1) - 1          % interpolate epsA/epsR values at steps within period
    period = [yr(i) , yr(i + 1)];
    epsA_vec{i} = interp1(period , epsA(i : i + 1 , 1) , ...
        yr(i) : step : yr(i + 1));
    epsR_vec{i} = interp1(period , epsR(i : i + 1 , 1) , ...
        yr(i) : step : yr(i + 1));
end
modelYrLast = endYear;

save(fullfile(savdir , 'mixInfectParams')  , 'epsA_vec' , ...
    'epsR_vec' , 'yr' , 'modelYr1' , 'modelYrLast' , ...
    'circProtect' , 'condProtect')

analProp = [0 , 0;
    0.5111 , 0.4266;
    0.5111 , 0.4266]; % [risk x gender]; proportion practicing anal sex
analProp = analProp .* 0; % no anal transmission for now
analTrans = [138 ; 11] ./ 10 ^ 4; % [gender x 1]; anal transmission probability
vagTransM = 8 / 10 ^ 4 * ones(size(analProp , 1) , 1);
vagTransF = 4 / 10 ^ 4 * ones(size(analProp , 1) , 1);
transM = vagTransM .* (1 - analProp(: , 1)) + analTrans(1) * analProp(: , 1);
transF = vagTransF .* (1 - analProp(: , 2)) + analTrans(2) * analProp(: , 2);
betaHIV_F2M = bsxfun(@times , [7 1 5.8 6.9 11.9 0.04;
    7 1 5.8 6.9 11.9 0.04;
    7 1 5.8 6.9 11.9 0.04] , transF);
betaHIV_M2F = bsxfun(@ times , [7 1 5.8 6.9 11.9 0.04;
    7 1 5.8 6.9 11.9 0.04;
    7 1 5.8 6.9 11.9 0.04] , transM);

file = [pwd , '\Config\Population_data.xlsx'];
maleActs = xlsread(file , 'Demographics' , 'D199 : F214'); % [age x risk]; male acts
femaleActs = xlsread(file , 'Demographics' , 'D219 : F234'); % [age x risk]; female acts
betaHIVF2M = zeros(age , risk , viral);
betaHIVM2F = betaHIVF2M;
for a = 1 : age % calculate per-partnership probability of HIV transmission
    betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
    betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
end
save(fullfile(savdir , 'vlBeta') , 'betaHIVF2M' , 'betaHIVM2F' , 'betaHIV_F2M' , 'betaHIV_M2F' , ...
    'maleActs', 'femaleActs')

%% Load indices
% prompt = {'Calculate Indices? (1 = yes, 0 = no)'};
% tit = 'Calculate Indices?';
% input = inputdlg(prompt , tit);
% if str2double(input{1})
%     disp('Preparing indices...')
%     disp('This may take a while...')
%     syms d v h s p g a r
%     calcInds();
% else
%     disp('Skipped index calculation')
% end

%% import and save HPV parameters
disp(' ')
disp('Importing HPV data...')
file = [pwd , '\Config\HPV_parameters.xlsx'];
disp(['Loading data from ' , file , '...']);
disp('This may take a couple seconds...');

% demographic data
beta_hrHPV_val = xlsread(file , 'HPV' , 'B4'); % per-parnership probability of HPV transmission (common vaccine types 16/18/31/33) (not currently used, calibrated)
beta_lrHPV_val = xlsread(file , 'HPV' , 'B5'); % oncogenic non-vaccine type (not currently used, calibrated)
rHivHpv = xlsread(file , 'HPV' , 'B20 : D23'); % HIV multipliers (<= 1) for HPV regression rate (not currently used, calibrated)
hivCin2 = xlsread(file , 'HPV' , 'B34 : B36'); % HIV multipliers for precancer progression (not currently used, calibrated)
hivCin3 = hivCin2; % (not currently used, calibrated)

muCC = xlsread(file , 'Cervical Cancer' , 'B6 : D11'); % undetected CC mortality
muCC = min(muCC .* 12 , 0.99); % convert undetected CC mortality rate from monthly to yearly
muCC_det = xlsread(file , 'Cervical Cancer' , 'B20 : D25'); % detected CC mortality
muCC_det = min(muCC_det .* 12 , 0.99); % convert detected CC mortality rate from monthly to yearly
kRL = xlsread(file , 'Cervical Cancer' , 'B33'); % local -> regional  progression
kDR = xlsread(file , 'Cervical Cancer' , 'B34'); % regional -> distant progression
detCC = xlsread(file , 'Cervical Cancer' , 'B42 : B44'); % [region x 1] detection probability (not currently used)
kCC = xlsread(file , 'Cervical Cancer' , 'B54 : B62'); % [period] group size (not currently used)
hivCC = xlsread(file , 'Cervical Cancer' , 'B74 : B77'); % HIV multipliers (HR) for cervical cancer mortality rate (not currently used)
hivCC(end) = 1;

% Currently no screening in model
kPap = xlsread(file , 'Screening and Treatment' , 'B4'); % pap smear rate (not currently used)
hpvSens = xlsread(file , 'Screening and Treatment' , 'B17 : C17'); % [1 x Cin] 2 , 3; HPV primary screening sensitivity (not currently used)
cytoSens = xlsread(file , 'Screening and Treatment' , 'B18 : C18'); % [1 x Cin] 2 , 3; Cytology ASC-US or worse primary screening sensitivity (not currently used)
leep = 1- xlsread(file , 'Screening and Treatment' , 'A29'); % rate of recurrence/persistence with LEEP (not currently used)
screenFreq = xlsread(file , 'Screening and Treatment' , 'B39 : C40'); % [Test x HIV status] (2 x 2); screening frequency (not currently used)
screenCover = xlsread(file , 'Screening and Treatment' , 'B48'); % screening coverage (not currently used)
ageStart = xlsread(file , 'Screening and Treatment' , 'B56'); % screening age start (not currently used)
ageEnd = xlsread(file , 'Screening and Treatment' , 'B57'); % screening age end (not currently used)

% CIN transition data
%HPV 16/18/Vaccine type ohr
kCin1_Inf(: , 1) = xlsread(file , 'CIN Transition' , 'B5 : B20'); % HPV to CIN1
kCin2_Cin1(: , 1) = xlsread(file , 'CIN Transition' , 'C5 : C20'); % CIN1 to CIN2
kCin3_Cin2(: , 1) = xlsread(file , 'CIN Transition', 'D5 : D20'); %CIN2 to CIN3
kCC_Cin3(: , 1) = xlsread(file , 'CIN Transition', 'E5 : E20'); % CIN3 to unlocalized

rNormal_Inf(: , 1) = xlsread(file , 'CIN Transition' , 'F5 : F20'); % HPV to Well (natural immunity)
kInf_Cin1(: , 1) = xlsread(file , 'CIN Transition' , 'G5 : G20'); % CIN1 to HPV
kCin1_Cin2(: , 1) = xlsread(file , 'CIN Transition', 'H5 : H20'); % CIN2 to CIN1
kCin2_Cin3(: , 1) = xlsread(file , 'CIN Transition', 'I5 : I20'); % CIN3 to CIN2

% 9v type ohr
kCin1_Inf(: , 2) = xlsread(file , 'CIN Transition' , 'B25 : B40');
kCin2_Cin1(: , 2) = xlsread(file , 'CIN Transition' , 'C25 : C40');
kCin3_Cin2(: , 2) = xlsread(file , 'CIN Transition' , 'D25 : D40');
kCC_Cin3(: , 2) = xlsread(file , 'CIN Transition' , 'E25 : E40');

rNormal_Inf(: , 2) = xlsread(file , 'CIN Transition' , 'F25 : F40');
kInf_Cin1(: , 2) = xlsread(file , 'CIN Transition' , 'G25 : G40');
kCin1_Cin2(: , 2) = xlsread(file , 'CIN Transition', 'H25 : H40');
kCin2_Cin3(: , 2) = xlsread(file , 'CIN Transition', 'I25 : I40');

%Non-vaccine type ohr
kCin1_Inf(: , 3) = xlsread(file , 'CIN Transition' , 'B45 : B60');
kCin2_Cin1(: , 3) = xlsread(file , 'CIN Transition' , 'C45 : C60');
kCin3_Cin2(: , 3) = xlsread(file , 'CIN Transition', 'D45 : D60');
kCC_Cin3(: , 3) = xlsread(file , 'CIN Transition', 'E45 : E60');

rNormal_Inf(: , 3) = xlsread(file , 'CIN Transition' , 'F45 : F60');
kInf_Cin1(: , 3) = xlsread(file , 'CIN Transition' , 'G45 : G60');
kCin1_Cin2(: , 3) = xlsread(file , 'CIN Transition', 'H45 : H60');
kCin2_Cin3(: , 3) = xlsread(file , 'CIN Transition', 'I45 : I60');

hpv_hivMult = flipud(xlsread(file , 'HPV' , 'B46 : D49')); % [inc CD4 count x vaccine type]; HPV incidence multiplier by type and CD4 count

% Weight HPV transitions and HPV incidence multiplier according to type distribution
distWeight = [0.7 , 0.2 , 0.1];
kCin1_Inf = sum(bsxfun(@times , kCin1_Inf , distWeight) , 2);
kCin2_Cin1 = sum(bsxfun(@times , kCin2_Cin1 , distWeight) , 2);
kCin3_Cin2 = sum(bsxfun(@times , kCin3_Cin2 , distWeight) , 2);
kCC_Cin3 = sum(bsxfun(@times , kCC_Cin3 , distWeight) , 2);
rNormal_Inf = sum(bsxfun(@times , rNormal_Inf , distWeight) , 2);
kInf_Cin1 = sum(bsxfun(@times , kInf_Cin1 , distWeight) , 2);
kCin1_Cin2 = sum(bsxfun(@times , kCin1_Cin2 , distWeight) , 2);
kCin2_Cin3 = sum(bsxfun(@times , kCin2_Cin3 , distWeight) , 2);
hpv_hivMult = sum(bsxfun(@times , hpv_hivMult , distWeight) , 2);

hpv_hivClear = xlsread(file , 'CIN Transition' , 'B73 : B76'); % [dec CD4 count]; HPV clearance multipliers
c3c2Mults = xlsread(file , 'CIN Transition' , 'B85 : B88'); % [dec CD4 count]; CIN2 to CIN3 multipliers
c2c1Mults = xlsread(file , 'CIN Transition' , 'B97 : B100'); % [dec CD4 count]; CIN1 to CIN2 multipliers
save(fullfile(savdir , 'hpvData') , 'beta_hrHPV_val' , 'beta_lrHPV_val' , 'kCC' , ...
    'rHivHpv' , 'hivCin2' , 'hivCin3' , 'muCC' , 'muCC_det' , 'kRL' , 'kDR' , 'detCC' , 'hivCC' , ...
    'kPap' , 'hpvSens' , 'cytoSens' , 'leep' , 'screenFreq' , 'ageStart' , 'ageEnd',...
    'kInf_Cin1' , 'kCin1_Cin2' , 'kCin2_Cin3' , 'kCin2_Cin1', ...
    'kCin3_Cin2' , 'kCC_Cin3' , 'kCin1_Inf' , 'rNormal_Inf' , 'hpv_hivClear' , ...
    'c3c2Mults' , 'c2c1Mults' , 'hpv_hivMult' , 'screenCover')
disp('HPV data loaded.')
disp(' ')

%% import and save weights and costs
clear
disp('Retrieving weights and costs...')
file = [pwd , '\Config\Weights_costs.xlsx'];
savdir = [pwd , '\Params']; 
hivTreatCost = xlsread(file , 'Costs' , 'B4 : B7'); % [inc CD4 count]; average HIV hositalization costs (not currently used, hard-coded in vaxCEA analysis)
artTreatCost = xlsread(file , 'Costs' , 'B8'); % HIV hospitalization costs on ART  (not currently used, hard-coded in vaxCEA analysis)
kCCDet = xlsread(file , 'Costs' , 'B17:B19'); %(local, regional, distant); probability of symptom detection
kCCDet = min(kCCDet .* 12 , 0.99); % convert monthly to yearly rate
vaxPrice = xlsread(file , 'Costs' , 'B26'); % bivalent HPV vaccine costs per vaccinated girl  (not currently used, hard-coded into vaxCEA analysis)
ccCost = xlsread(file , 'Costs' , 'B34:B36'); % (local, regional, distant); CC costs (later re-defined in vaxCEA analysis)
hivTreatCost = flipud(hivTreatCost); % (not currently used)
save(fullfile(savdir , 'cost_weights'))
disp('Done')

%% import and save calibration data
% dimensions = [pos , N]
clear
file = [pwd , '\Config\Calibration_targets.xlsx'];
savdir = [pwd , '\Params']; 
cinPos2014_obs = xlsread(file , 'Calibration' , 'D2 : F11'); %CIN2/CIN3 Prevalence (HIV+) 2014, by age
cinNeg2014_obs = xlsread(file , 'Calibration' , 'D12 : F21'); %CIN2/CIN3 Prevalence (HIV-) 2014, by age

% hpv_hiv_2008_obs = xlsread(file , 'Calibration' , 'D42 : F51'); % HPV Prevalence in HIV+ Women (no CIN2/3) 2008, by age
% hpv_hivNeg_2008_obs = xlsread(file , 'Calibration' , 'D32 : F41'); % HPV Prevalence in HIV- Women (no CIN2/3) 2008, by age

hpv_hiv_obs = xlsread(file , 'Calibration' , 'D144 : F152'); % HPV Prevalence in HIV+ Women (All) 2014, by age
hpv_hivNeg_obs = xlsread(file , 'Calibration' , 'D153 : F161'); % HPV Prevalence in HIV- Women (All) 2014, by age

hpv_hivM2008_obs = xlsread(file , 'Calibration' , 'E52 : F55'); % HPV Prevalence in HIV+ Men, by age
hpv_hivMNeg2008_obs = xlsread(file , 'Calibration' , 'E56 : F59'); % HPV Prevalence in HIV- Men, by age

hivPrevM_obs = xlsread(file , 'Calibration' , 'D60 : F101'); % HIV Prevalence in Men 2003,2005,2006,2007,2008,2009, by age
hivPrevF_obs = xlsread(file , 'Calibration' , 'D102 : F143'); % HIV Prevalence in Women 2003,2005,2006,2007,2008,2009, by age
save(fullfile(savdir , 'calibData'))
clear

