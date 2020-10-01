% Load and save parameters, calibration data, indices, and matrices
function[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd, annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , partnersMmult, maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4,...
    mue , mue2 , mue3 , mue4 , epsA_vec , epsR_vec , yr , ...
    hivOn , betaHIV_mod , muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , waning , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear ,circNatStartYear , vaxStartYear , baseline , cisnet , who , whob , ...
    circProtect , condProtect , MTCTRate , hyst , OMEGA , ...
    ccInc2012_dObs , cc_dist_dObs , cin3_dist_dObs , cin1_dist_dObs , ...
    hpv_dist_dObs , cinPos2007_dObs , cin1_2010_dObs , cin2_2010_dObs, ...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hpv_all_dObs , hpv_hiv2009_dObs , ...
    hivPrevF_dObs , hivPrevM_dObs , hivPrevAll_dObs, popAgeDist_dObs , totPopSize_dObs , hivCurr , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ageInd , riskInd , ...
    hivNegNonVMMCinds , hivNegVMMCinds , vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3  , dFertMat3, d_partnersMmult, riskAdj, d_riskAdj, ...
    deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2(fivYrAgeGrpsOn , calibBool , pIdx , paramsSub , paramSet)

tic

paramDir = [pwd , '/Params/'];

%% Set and save general parameters

% Time
stepsPerYear = 6;
timeStep = 1 / stepsPerYear;
startYear = 1925;
currYear = 2020;
endYear = currYear; %2015; %currYear;
years = endYear - startYear;

% Compartments
disease = 8;
viral = 6;
hpvVaxStates = 7;
hpvNonVaxStates = 7;
endpoints = 4;
intervens = 4;
gender = 2;
age = 80 / max(1,fivYrAgeGrpsOn*5);
risk = 3;

hpvTypeGroups = 2;

dim = [disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , gender , age , risk];

% Index retrieval function
k = cumprod([disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , gender , age]);

toInd = @(x) (x(: , 9) - 1) * k(8) + (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + ...
    (x(: , 6) - 1) * k(5) + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + ...
    (x(: , 3) - 1) * k(2) + (x(: , 2) - 1) * k(1) + x(: , 1);

annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear));

%% Save demographic and behavioral parameters

% Import from Excel initial population size, male risk distribution, 

%background mortality, fertility, partnerships, and acts per partnership
file = [pwd , '/Config/Kenya_parameters_Feb20.xlsx'];
popInit = xlsread(file , 'Population' , 'B169:C184'); 
popInit = popInit .* 1000; 
riskDistM = xlsread(file , 'Sexual behavior' , 'F73:H88');
mue = zeros(age , gender);
mue(: , 1) = xlsread(file , 'Mortality' , 'C94:C109'); %1950
mue(: , 2) = xlsread(file , 'Mortality' , 'D94:D109');
mue2 = zeros(age , gender);
mue2(: , 1) = xlsread(file , 'Mortality' , 'G94:G109'); %1985
mue2(: , 2) = xlsread(file , 'Mortality' , 'H94:H109');
mue3 = zeros(age , gender);
mue3(: , 1) = xlsread(file , 'Mortality' , 'K94:K109'); % 2000
mue3(: , 2) = xlsread(file , 'Mortality' , 'L94:L109');
mue3(1, :)= mue3(1 , :) .* 2 ;
mue4 = zeros(age , gender);
mue4(: , 1) = xlsread(file , 'Mortality' , 'O94:O109'); % 2020
mue4(: , 2) = xlsread(file , 'Mortality' , 'P94:P109');
mue4(1, :)= mue4(1 , :) .* 2 ;
fertility = xlsread(file , 'Fertility' , 'D104:I119');
partnersM = xlsread(file , 'Sexual behavior' , 'O73:Q88');
partnersF = xlsread(file , 'Sexual behavior' , 'L73:N88');
maleActs = xlsread(file , 'Sexual behavior' , 'D168:F183');
femaleActs = xlsread(file , 'Sexual behavior' , 'D188:F203');
save(fullfile(paramDir ,'demoParamsFrmExcel'), 'popInit' , 'riskDistM' , ...
    'mue' , 'mue2' , 'mue3' , 'mue4' , 'fertility' , 'partnersM' ,  'partnersF' , 'maleActs' , 'femaleActs');

% Load pre-saved initial population size by age and gender, male risk distribution by age, 
% background mortality by age and gender, and fertility by age and gender
load([paramDir , 'demoParamsFrmExcel'] , 'popInit' , 'riskDistM' , 'mue' , ...
    'mue2' , 'mue3' , 'mue4' , 'fertility');

% Set female risk distribution
riskDistM(3, 1:3) = [0.971839 0.02640375 0.00175725]; % 95% of the original 15-19 dist
riskDistM(4, 1:3) = [0.887356 0.105615 0.007029]; % 300% of the original 15-19 dist
riskDistF = riskDistM;

% Calculate fertility2 and fertility3 matrices
if calibBool && any(36 == pIdx);
    idx = find(36 == pIdx);
    fertDeclineProp = paramSet(paramsSub{idx}.inds(:));
else
    fertDeclineProp = [0.7 ; 0.4; 0.6];
end
fertility2 = fertility .* fertDeclineProp(1,1);
fertility3 = fertility2 .* fertDeclineProp(2,1);
fertility4 = fertility3 .* fertDeclineProp(3,1);

% Male partners per year by age and risk group
if calibBool && any(1 == pIdx)
    idx = find(1 == pIdx);
    %partnersM = zeros(age , risk);
    %partnersMmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    load([paramDir , 'demoParamsFrmExcel'] , 'partnersM');
    partnersM(3 , 2: 3) = partnersM(3 , 2:3) + 1;
    partnersM(4:16, 2:3) = zeros;
    partnersM(4:6, 2) = paramSet(paramsSub{idx}.inds(1)) .* partnersM(4:6, 1);
    partnersM(7:9, 2) = paramSet(paramsSub{idx}.inds(2)) .* partnersM(7:9, 1);
    partnersM(10:age, 2) = paramSet(paramsSub{idx}.inds(3)) .* partnersM(10:age, 1);
    partnersM(4:6, 3) = paramSet(paramsSub{idx}.inds(4)) .* partnersM(4:6, 1);
    partnersM(7:10, 3) = paramSet(paramsSub{idx}.inds(5)) .* partnersM(7:10, 1);
    partnersM(11:age, 3) = paramSet(paramsSub{idx}.inds(6)) .* partnersM(11:age, 1);
%     partnersM(3 , 2) = paramSet(paramsSub{idx}.inds(7))*partnersM(3 , 3);
%     partnersM(3 , 1) = paramSet(paramsSub{idx}.inds(8))*partnersM(3 , 2);
%     for a = 4 : age
%         partnersM(a , 1:2) = (partnersM(a , 3)/partnersM(a-1 , 3)) .* partnersM(a-1 , 1:2);
%     end
else 
    load([paramDir , 'demoParamsFrmExcel'] , 'partnersM');
    
    partnersM(3 , 2:3) = partnersM(3, 2:3) + 1;
    partnersM(4 , 2:3) = partnersM(4, 2:3);
    partnersM(5:6 , 3) = partnersM(5:6 , 3) ;
    
end

% Female partners per year by age and risk group
if calibBool && any(2 == pIdx)
    idx = find(2 == pIdx);
    %partnersF = zeros(age , risk);
    %partnersFmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    load([paramDir , 'demoParamsFrmExcel'] , 'partnersF');
    partnersF(3, 2: 3) = partnersF(3 , 2:3) + 1;
    partnersF(4:16, 2:3) = zeros;
    partnersF(4:6, 2) = paramSet(paramsSub{idx}.inds(1)) .* partnersF(4:6, 1);
    partnersF(7:9, 2) = paramSet(paramsSub{idx}.inds(2)) .* partnersF(7:9, 1);
    partnersF(10:age, 2) = paramSet(paramsSub{idx}.inds(3)) .* partnersF(10:age, 1);
    partnersF(4:6, 3) = paramSet(paramsSub{idx}.inds(4)) .* partnersF(4:6, 1);
    partnersF(7:9, 3) = paramSet(paramsSub{idx}.inds(5)) .* partnersF(7:9, 1);
    partnersF(10:age, 3) = paramSet(paramsSub{idx}.inds(6)) .* partnersF(10:age, 1);
%     partnersF(3 , 2) = paramSet(paramsSub{idx}.inds(7))*partnersF(3 , 3);
%     partnersF(3 , 1) = paramSet(paramsSub{idx}.inds(8))*partnersF(3 , 2);
%     for a = 4 : age
%         partnersF(a , 1:2) = (partnersF(a , 3)/partnersF(a-1 , 3)) .* partnersF(a-1 , 1:2);
%     end
else
    load([paramDir , 'demoParamsFrmExcel'] , 'partnersF');
    partnersF(3 , 2: 3) = partnersF(3 , 2:3) + 1;
    partnersF(4 , 2: 3) = partnersF(4 , 2:3); % .* partnersMmult(2);
    partnersF(5:6 , 1: 3) = partnersF(5:6 , 1:3) ;
    partnersF(7:10 , 1: 3) = partnersF(7:10 , 1:3) ;
end    

% Male acts per partnership per year by age and risk group
if calibBool && any(8 == pIdx)
    idx = find(8 == pIdx);
    %maleActs = zeros(age , risk);
    load([paramDir , 'demoParamsFrmExcel'] , 'maleActs');
    %maleActsmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    maleActs(3:4, 1:risk) = maleActs(3:4, 1:risk) .* paramSet(paramsSub{idx}.inds(1)); 
    maleActs(5, 1:risk) = maleActs(5, 1:risk) .* paramSet(paramsSub{idx}.inds(2)) ;
    maleActs(6:7, 1:risk) = maleActs(6:7, 1:risk) .* paramSet(paramsSub{idx}.inds(3));
    maleActs(10:16, 1:risk) = maleActs(10:16, 1:risk).* paramSet(paramsSub{idx}.inds(4));
else
    load([paramDir , 'demoParamsFrmExcel'] , 'maleActs');
    maleActs(3:4, 1:risk) = maleActs(3:4, 1:risk) .* 5; 
    maleActs(5, 1:risk) = maleActs(5, 1:risk) .* 2.5 ;
    maleActs(6:7, 1:risk) = maleActs(6:7, 1:risk) .* 0.6;
    maleActs(10:16, 1:risk) = maleActs(10:16, 1:risk).* 1.5 ;
end

% Female acts per partnership per year by age and risk group
if calibBool && any(9 == pIdx)
    idx = find(9 == pIdx);
%     femaleActs = zeros(age , risk);
%     femaleActsmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    load([paramDir , 'demoParamsFrmExcel'] , 'femaleActs');
    femaleActs(3 : 4, 1:risk) = femaleActs(3 : 4, 1:risk) .* paramSet(paramsSub{idx}.inds(1));
    femaleActs(5, 1:risk) = femaleActs(5, 1:risk) .* paramSet(paramsSub{idx}.inds(2)) ;
    femaleActs(6:7, 1:risk) = femaleActs(6:7, 1:risk).*paramSet(paramsSub{idx}.inds(3)) ;
    femaleActs(10:16, 1:risk) = femaleActs(10:16, 1:risk).* paramSet(paramsSub{idx}.inds(4)) ;
else
    load([paramDir , 'demoParamsFrmExcel'] , 'femaleActs');
    femaleActs(3 : 4, 1:risk) = femaleActs(3 : 4, 1:risk) .* 6;
    femaleActs(5, 1:risk) = femaleActs(5, 1:risk) .* 2.6 ;
    femaleActs(6:7, 1:risk) = femaleActs(6:7, 1:risk).* 0.6 ;
    femaleActs(10:16, 1:risk) = femaleActs(10:16, 1:risk).* 1.5 ;
end

% Convert 5-year age groups to 1-year age groups
if ~fivYrAgeGrpsOn
    % Divide popInit age groups equally into five
    popInit_orig = popInit;
    [ageDim, valDim] = size(popInit_orig);
    popInit = zeros(ageDim*5 , valDim);
    for i = 1 : ageDim
        popInit(((i-1)*5+1) : i*5 , :) = ones(5 , valDim) .* (popInit_orig(i , :)./5);
    end

    % Replicate rates across single age groups for other variables
    vars5To1_nms = {'riskDistM' , 'riskDistF' , 'mue' , 'mue2' , 'mue3' , 'mue4' , ...
        'fertility' , 'fertility2' , 'fertility3' , 'fertility4' ,...
        'partnersM' , 'partnersF' , 'maleActs' , 'femaleActs'};
    vars5To1_vals = {riskDistM , riskDistF , mue , mue2 , mue3 , mue4 , ...
        fertility , fertility2 , fertility3 , fertility4 ,...
        partnersM , partnersF , maleActs , femaleActs};    
    for j = 1 : length(vars5To1_vals)
        valsA1 = age5To1(vars5To1_vals{j});
        assignin('base', vars5To1_nms{j} , valsA1);
    end
end

% Rename initial population size variables
mInit = popInit(: , 1); % initial male population size by age
fInit = popInit(: , 2); % initial female population size by age

% Rename risk distribution variable
riskDist = zeros(age , risk , 2);
riskDist(: , : , 1) = riskDistM;
riskDist(: , : , 2) = riskDistF;

ageSexDebut = (10/max(1 , fivYrAgeGrpsOn*5)+1);

% Mixing by age group
if calibBool && any(6 == pIdx);
    idx = find(6 == pIdx);
    %epsA = paramSet(paramsSub{idx}.inds(:));
    epsA = ones(3,1).*paramSet(paramsSub{idx}.inds(:));
else
    epsA = [0.3 ; 0.3 ; 0.3];
end
% Mixing by risk group
if calibBool && any(7 == pIdx);
    idx = find(7 == pIdx);
    %epsR = paramSet(paramsSub{idx}.inds(:));
    epsR = ones(3,1).*paramSet(paramsSub{idx}.inds(:));
else
    epsR = [0.3 ; 0.3 ; 0.3];
end
yr = [1985; 1990; 2000];
epsA_vec = cell(size(yr , 1) - 1, 1); % save data over time interval in a cell array
epsR_vec = cell(size(yr , 1) - 1, 1);
for i = 1 : size(yr , 1) - 1          % interpolate epsA/epsR values at steps within period
    period = [yr(i) , yr(i + 1)];
    epsA_vec{i} = interp1(period , epsA(i : i + 1 , 1) , ...
        yr(i) : timeStep : yr(i + 1));
    epsR_vec{i} = interp1(period , epsR(i : i + 1 , 1) , ...
        yr(i) : timeStep : yr(i + 1));
end

%% Save HIV natural history parameters
hivOn = 1; % bool to turn HIV on or off although model calibrated for HIV on

% Import from Excel HIV-associated death rate and CD4/VL transition rates
% file = [pwd , '/Config/HIV_parameters_Kenya.xlsx'];
% muHIV = xlsread(file , 'Disease Data' , 'B20:G35');
% kCD4male = xlsread(file , 'Disease Data' , 'C45:F60');
% kCD4female = xlsread(file , 'Disease Data' , 'C62:F77');
% kVlmale = xlsread(file , 'Disease Data' , 'C84:F99');
% kVlfemale = xlsread(file , 'Disease Data' , 'C101:F116');
% save(fullfile(paramDir ,'hivNHParamsFrmExcel'), 'muHIV' , 'kCD4male' , ...
%     'kCD4female' , 'kVlmale' , 'kVlfemale');

% Load pre-saved HIV-associated death rate and CD4/VL transition matrices by age and gender
load([paramDir , 'hivNHParamsFrmExcel'] , 'muHIV' , 'kCD4male' , ...
    'kCD4female' , 'kVlmale' , 'kVlfemale');

% Convert 5-year age groups to 1-year age groups
if ~fivYrAgeGrpsOn
    % Replicate rates across single age groups for other variables
    vars5To1_nms = {'muHIV' , 'kVlmale' , 'kVlfemale' , 'kCD4male' , 'kCD4female'};
    vars5To1_vals = {muHIV , kVlmale , kVlfemale , kCD4male , kCD4female};    
    for j = 1 : length(vars5To1_vals)
        valsA1 = age5To1(vars5To1_vals{j});
        assignin('base', vars5To1_nms{j} , valsA1);
    end
end

% HIV Vl state transition matrices
kVl = zeros(age , 4 , gender);
kVl(: , : , 1) = kVlmale;
kVl(: , : , 2) = kVlfemale;

% HIV CD4 state transition matrices
kCD4 = zeros(age , 4 , gender);
kCD4(: , : , 1) = kCD4male;
kCD4(: , : , 2) = kCD4female;

% Baseline HIV vaginal transmission rate
if calibBool && any(35 == pIdx);
    idx = find(35 == pIdx);
    baseVagTrans = paramSet(paramsSub{idx}.inds(:));
else
    baseVagTrans = [0.001]; %Thomson et al: For a 25-year-old woman not taking PrEP, 
                            %the HIV acquisition probability per condomless sex act with an HIV-infected 
                            %male partner with a viral load of 10 000 copies/mL was 0.0011 (95% CI 0.005–0.0019)
end

% HIV tranmission rate % make HIV M-> F trans the smae 
vagTransM = (baseVagTrans) * ones(risk , 1); % probability of transmission from male (insertive) to female (receptive) based on male's disease state; female acquisition 
vagTransF = baseVagTrans * ones(risk , 1); % probability of transmission from female (receptive) to male (insertive) based on female's disease state; male acquisition 
betaHIV_F2M = bsxfun(@times , [9.0 1.0 2.5 7.0 0.7 0.0; 9.0 1.0 2.5 7.0 0.7 0.0; 9.0 1.0 2.5 7.0 0.7 0.0] , vagTransF);
betaHIV_M2F = bsxfun(@times , [9.0 1.0 2.5 7.0 0.7 0.0; 9.0 1.0 2.5 7.0 0.7 0.0; 9.0 1.0 2.5 7.0 0.7 0.0] , vagTransM);
betaHIV_F2M_red = bsxfun(@times , [9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0; 9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0; 9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0] , vagTransF);
betaHIV_M2F_red = bsxfun(@times , [9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0; 9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0; 9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0] , vagTransM);
betaHIV = zeros(gender , age , risk , viral);
betaHIV_red = zeros(gender , age , risk , viral);
for a = 1 : age % calculate per-partnership probability of HIV transmission
    % force of infection: females to infect HIV-negative males, 
    % affected by betaHIV_F2M, probability of transmission from female (receptive) to male(insertive) based on female's disease state), and number of male acts
    betaHIV(2 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)'));
    betaHIV_red(2 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M_red , maleActs(a , :)'));
    % force of infection: males to infect HIV-negative females,
    % affected by betaHIV_M2F, probability of transmission from male (insertive) to female (receptive) based on male's disease state), and number of female acts
    betaHIV(1 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); 
    betaHIV_red(1 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F_red , femaleActs(a , :)')); 
end

betaHIV_mod = zeros(risk , age , viral , endpoints , gender);
for v = 1 : viral
    for x = 1 : endpoints
        for g = 1 : gender
            for a = 1 : age
                for r = 1 : risk
                    if (x > 1)
                        betaHIV_mod(r , a , v , x , g) = betaHIV_red(g , a , r , v);
                    else
                        betaHIV_mod(r , a , v , x , g) = betaHIV(g , a , r , v);
                    end
                end
            end
        end
    end
end

%% Import HPV/CIN/CC transition data from Excel
% file = [pwd , '/Config/HPV_parameters.xlsx'];
% 
% ageTrends = xlsread(file , 'CIN transitions by type' , 'C38 : J40');
% 
% kCin1_Inf_orig = zeros(1,2);
% kCin2_Cin1_orig = zeros(1,2);
% kCin3_Cin2_orig = zeros(1,2);
% kCC_Cin3_orig = zeros(1,2);
% rNormal_Inf_orig = zeros(1,2);
% kInf_Cin1_orig = zeros(1,2);
% kCin1_Cin2_orig = zeros(1,2);
% kCin2_Cin3_orig = zeros(1,2);
% % 9v HPV types
% kCin1_Inf_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'C25'); % HPV to CIN1, ages 10-24
% kCin2_Cin1_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'D25'); % CIN1 to CIN2
% kCin3_Cin2_orig(1, 1) = xlsread(file , 'CIN transitions by type', 'E25'); % CIN2 to CIN3
% kCC_Cin3_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'F25'); % CIN3 to unlocalized
% rNormal_Inf_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'G25'); % HPV to Well (natural immunity)
% kInf_Cin1_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'H25'); % CIN1 to HPV
% kCin1_Cin2_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'I25'); % CIN2 to CIN1
% kCin2_Cin3_orig(1 , 1) = xlsread(file , 'CIN transitions by type' , 'J25'); % CIN3 to CIN2
% % Non-9v HPV types
% kCin1_Inf_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'M25');  % HPV to CIN1, ages 10-24
% kCin2_Cin1_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'N25'); % CIN1 to CIN2
% kCin3_Cin2_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'O25'); %CIN2 to CIN3
% kCC_Cin3_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'P25'); % CIN3 to unlocalized
% rNormal_Inf_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'Q25'); % HPV to Well (natural immunity)
% kInf_Cin1_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'R25'); % CIN1 to HPV
% kCin1_Cin2_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'S25'); % CIN2 to CIN1
% kCin2_Cin3_orig(1 , 2) = xlsread(file , 'CIN transitions by type' , 'T25'); % CIN3 to CIN2

%% Set HPV/CIN/CC transitions by age group based on pre-saved or calibrated values

% Set transitions in first age group from pre-saved or calibrated values
load([paramDir , 'hpvCinCCTransFrmExcel'] , 'kCin1_Inf_orig' , 'kCin2_Cin1_orig' , ...
    'kCin3_Cin2_orig' , 'kCC_Cin3_orig' , 'rNormal_Inf_orig' , 'kInf_Cin1_orig' , ...
    'kCin1_Cin2_orig' , 'kCin2_Cin3_orig') % , 'ageTrends');

file = [pwd , '/Config/HPV_parameters.xlsx'];
ageTrends = xlsread(file , 'CIN transitions by type' , 'C45:J55');

kCin1_Inf = zeros(16,2);  %(:, 1) 9vHPV, (:, 2) non-9vHPV
kCin2_Cin1 = zeros(16,2);
kCin3_Cin2 = zeros(16,2);
kCC_Cin3 = zeros(16,2);
rNormal_Inf = zeros(16,2);
kInf_Cin1 = zeros(16,2);
kCin1_Cin2 = zeros(16,2);
kCin2_Cin3 = zeros(16,2);

% HPV to CIN1, ages 10-24
if calibBool && any(27 == pIdx)
    idx = find(27 == pIdx);
    kCin1_InfMult = paramSet(paramsSub{idx}.inds(:));
    kCin1_Inf(1 : 5 , 1) = kCin1_Inf_orig(1 , 1) * kCin1_InfMult(1);
    kCin1_Inf(1 : 5 , 2) = kCin1_Inf_orig(1 , 2) * kCin1_InfMult(2);
else
    kCin1_Inf(1 : 5 , 1) = kCin1_Inf_orig(1 , 1) * 2.2;
    kCin1_Inf(1 : 5 , 2) = kCin1_Inf_orig(1 , 2) * 3.2; % original multiplier = 1.6
end

% CIN1 to CIN2, ages 10-24
if calibBool && any(28 == pIdx)
    idx = find(28 == pIdx);
    kCin2_Cin1Mult = paramSet(paramsSub{idx}.inds(:));
    kCin2_Cin1(1 : 5 , 1) = kCin2_Cin1_orig(1 , 1) * kCin2_Cin1Mult(1);
    kCin2_Cin1(1 : 5 , 2) = kCin2_Cin1_orig(1 , 2) * kCin2_Cin1Mult(2);
else
    kCin2_Cin1(1 : 5 , 1) = kCin2_Cin1_orig(1 , 1)*2.8;
    kCin2_Cin1(1 : 5 , 2) = kCin2_Cin1_orig(1 , 2)*5; %original multiplier = 1.7
end

% CIN2 to CIN3, ages 10-24
if calibBool && any(29 == pIdx)
    idx = find(29 == pIdx);
    kCin3_Cin2Mult = paramSet(paramsSub{idx}.inds(:));
    kCin3_Cin2(1 : 5 , 1) = kCin3_Cin2_orig(1 , 1) * kCin3_Cin2Mult(1);
    kCin3_Cin2(1 : 5 , 2) = kCin3_Cin2_orig(1 , 2) * kCin3_Cin2Mult(2);
else
    kCin3_Cin2(1 : 5 , 1) = kCin3_Cin2_orig(1 , 1) * 1.3;
    kCin3_Cin2(1 : 5 , 2) = kCin3_Cin2_orig(1 , 2)* 3.3; %original multiplier 1.7
end

% CIN3 to unlocalized cancer, ages 10-24
if calibBool && any(30 == pIdx)
    idx = find(30 == pIdx);
    kCC_Cin3Mult = paramSet(paramsSub{idx}.inds(:));
    kCC_Cin3(1 : 5 , 1) = kCC_Cin3_orig(1 , 1) * kCC_Cin3Mult(1);
    kCC_Cin3(1 : 5 , 2) = kCC_Cin3_orig(1 , 2) * kCC_Cin3Mult(2);
else
    kCC_Cin3(1 : 5 , 1) = kCC_Cin3_orig(1 , 1) .* 2.1;
    kCC_Cin3(1 : 5 , 2) = kCC_Cin3_orig(1 , 2) .* 6.7;
end

% HPV to Well (or natural immunity), ages 10-24
if calibBool && any(31 == pIdx)
    idx = find(31 == pIdx);
    rNormal_InfMult = paramSet(paramsSub{idx}.inds(:));
    rNormal_Inf(1 : 5 , 1) = rNormal_Inf_orig(1 , 1) * rNormal_InfMult(1);
    rNormal_Inf(1 : 5 , 2) = rNormal_Inf_orig(1 , 2) * rNormal_InfMult(2);
else
    rNormal_Inf(1 : 5 , 1) = rNormal_Inf_orig(1 , 1) * 2.9;
    rNormal_Inf(1 : 5 , 2) = rNormal_Inf_orig(1 , 2) * 1.9;
end

% CIN1 to HPV, ages 10-24
if calibBool && any(32 == pIdx)
    idx = find(32 == pIdx);
    kInf_Cin1Mult = paramSet(paramsSub{idx}.inds(:));
    kInf_Cin1(1 : 5 , 1) = kInf_Cin1_orig(1 , 1) * kInf_Cin1Mult(1);
    kInf_Cin1(1 : 5 , 2) = kInf_Cin1_orig(1 , 2) * kInf_Cin1Mult(2);
else
    kInf_Cin1(1 : 5 , 1) = kInf_Cin1_orig(1 , 1) * 1.7;
    kInf_Cin1(1 : 5 , 2) = kInf_Cin1_orig(1 , 2) * .93;
end

% CIN2 to CIN1, ages 10-24
if calibBool && any(33 == pIdx)
    idx = find(33 == pIdx);
    kCin1_Cin2Mult = paramSet(paramsSub{idx}.inds(:));
    kCin1_Cin2(1 : 5 , 1) = kCin1_Cin2_orig(1 , 1) * kCin1_Cin2Mult(1);
    kCin1_Cin2(1 : 5 , 2) = kCin1_Cin2_orig(1 , 2) * kCin1_Cin2Mult(2);
else
    kCin1_Cin2(1 : 5 , 1) = kCin1_Cin2_orig(1 , 1) * 1.25;
    kCin1_Cin2(1 : 5 , 2) = kCin1_Cin2_orig(1 , 2) * 0.35;
end

% CIN3 to CIN2, ages 10-24
if calibBool && any(34 == pIdx)
    idx = find(34 == pIdx);
    kCin2_Cin3Mult = paramSet(paramsSub{idx}.inds(:));
    kCin2_Cin3(1 : 5 , 1) = kCin2_Cin3_orig(1 , 1) * kCin2_Cin3Mult(1);
    kCin2_Cin3(1 : 5 , 2) = kCin2_Cin3_orig(1 , 2) * kCin2_Cin3Mult(2);
else
    kCin2_Cin3(1 : 5 , 1) = kCin2_Cin3_orig(1 , 1)* 1 ;
    kCin2_Cin3(1 : 5 , 2) = kCin2_Cin3_orig(1 , 2) * 0.4;
end
% 
% % Apply age trends to 9v HPV transitions
% kCin1_Inf(6 : 10 , 1) = kCin1_Inf(1 , 1) * ageTrends(1,1); % ages 25-49
% kCin2_Cin1(6 : 10 , 1) = kCin2_Cin1(1 , 1) * ageTrends(1,2);
% kCin3_Cin2(6 : 10 , 1) = kCin3_Cin2(1, 1) * ageTrends(1,3);
% kCC_Cin3(6 : 10 , 1) = kCC_Cin3(1 , 1) * ageTrends(1,4);
% rNormal_Inf(6 : 10 , 1) = rNormal_Inf(1 , 1) * ageTrends(1,5) ; %include 20-24
% kInf_Cin1(6 : 10 , 1) = kInf_Cin1(1 , 1) * ageTrends(1,6);
% kCin1_Cin2(6 : 10 , 1) = kCin1_Cin2(1 , 1) * ageTrends(1,7);
% kCin2_Cin3(6 : 10 , 1) = kCin2_Cin3(1 , 1) * ageTrends(1,8);
% kCin1_Inf(11 : 14 , 1) = kCin1_Inf(1 , 1) * ageTrends(2,1); % ages 50-69
% kCin2_Cin1(11 : 14 , 1) = kCin2_Cin1(1 , 1) * ageTrends(2,2);
% kCin3_Cin2(11 : 14 , 1) = kCin3_Cin2(1 , 1) * ageTrends(2,3);
% kCC_Cin3(11 : 14 , 1) = kCC_Cin3(1 , 1) * ageTrends(2,4) * .8;
% rNormal_Inf(11 : 14 , 1) = rNormal_Inf(1 , 1) * ageTrends(2,5) ;
% kInf_Cin1(11 : 14 , 1) = kInf_Cin1(1 , 1) * ageTrends(2,6);
% kCin1_Cin2(11 : 14 , 1) = kCin1_Cin2(1 , 1) * ageTrends(2,7);
% kCin2_Cin3(11 : 14 , 1) = kCin2_Cin3(1 , 1) * ageTrends(2,8);
% kCin1_Inf(15 : 16 , 1) = kCin1_Inf(1 , 1) * ageTrends(3,1); % ages 70-79
% kCin2_Cin1(15 : 16 , 1) = kCin2_Cin1(1 , 1) * ageTrends(3,2);
% kCin3_Cin2(15 : 16 , 1) = kCin3_Cin2(1 , 1) * ageTrends(3,3);
% kCC_Cin3(15 : 16 , 1) = kCC_Cin3(1 , 1) * ageTrends(3,4) * .9 ;
% rNormal_Inf(15 : 16 , 1) = rNormal_Inf(1 , 1) * ageTrends(3,5) ;
% kInf_Cin1(15 : 16 , 1) = kInf_Cin1(1 , 1) * ageTrends(3,6);
% kCin1_Cin2(15 : 16 , 1) = kCin1_Cin2(1 , 1) * ageTrends(3,7);
% kCin2_Cin3(15 : 16 , 1) = kCin2_Cin3(1 , 1) * ageTrends(3,8);
% 
% % Apply age trends to non-9v HPV transitions
% kCin1_Inf(6 : 10 , 2) = kCin1_Inf(1 , 2) * ageTrends(1,1); % ages 25-49
% kCin2_Cin1(6 : 10 , 2) = kCin2_Cin1(1 , 2) * ageTrends(1,2);
% kCin3_Cin2(6 : 10 , 2) = kCin3_Cin2(1, 2) * ageTrends(1,3);
% kCC_Cin3(6 : 10 , 2) = kCC_Cin3(1 , 2) * ageTrends(1,4);
% rNormal_Inf(6 : 10 , 2) = rNormal_Inf(1 , 2) * ageTrends(1,5) ;
% kInf_Cin1(6 : 10 , 2) = kInf_Cin1(1 , 2) * ageTrends(1,6);
% kCin1_Cin2(6 : 10 , 2) = kCin1_Cin2(1 , 2) * ageTrends(1,7);
% kCin2_Cin3(6 : 10 , 2) = kCin2_Cin3(1 , 2) * ageTrends(1,8);
% kCin1_Inf(11 : 14 , 2) = kCin1_Inf(1 , 2) * ageTrends(2,1); % ages 50-69
% kCin2_Cin1(11 : 14 , 2) = kCin2_Cin1(1 , 2) * ageTrends(2,2);
% kCin3_Cin2(11 : 14 , 2) = kCin3_Cin2(1 , 2) * ageTrends(2,3);
% kCC_Cin3(11 : 14 , 2) = kCC_Cin3(1 , 2) * ageTrends(2,4) * .8 ;
% rNormal_Inf(11 : 14 , 2) = rNormal_Inf(1 , 2) * ageTrends(2,5);
% kInf_Cin1(11 : 14 , 2) = kInf_Cin1(1 , 2) * ageTrends(2,6);
% kCin1_Cin2(11 : 14 , 2) = kCin1_Cin2(1 , 2) * ageTrends(2,7);
% kCin2_Cin3(11 : 14 , 2) = kCin2_Cin3(1 , 2) * ageTrends(2,8);
% kCin1_Inf(15 : 16 , 2) = kCin1_Inf(1 , 2) * ageTrends(3,1); % ages 70-79
% kCin2_Cin1(15 : 16 , 2) = kCin2_Cin1(1 , 2) * ageTrends(3,2);
% kCin3_Cin2(15 : 16 , 2) = kCin3_Cin2(1 , 2) * ageTrends(3,3);
% kCC_Cin3(15 : 16 , 2) = kCC_Cin3(1 , 2) * ageTrends(3,4) * .9 ;
% rNormal_Inf(15 : 16 , 2) = rNormal_Inf(1 , 2) * ageTrends(3,5) ;
% kInf_Cin1(15 : 16 , 2) = kInf_Cin1(1 , 2) * ageTrends(3,6);
% kCin1_Cin2(15 : 16 , 2) = kCin1_Cin2(1 , 2) * ageTrends(3,7);
% kCin2_Cin3(15 : 16 , 2) = kCin2_Cin3(1 , 2) * ageTrends(3,8);

% HPV progression multipler 5-year age group
%%
for a = 1 : length(ageTrends)
    mult = 1;
    if a >= 7
        mult = 1.29 .* (1 - (1/a));
    end
    % Apply age trends to 9v HPV transitions
    kCin1_Inf(a + 5 , 1) = kCin1_Inf(1 , 1) * ageTrends(a,1) * mult; 
    kCin2_Cin1(a + 5 , 1) = kCin2_Cin1(1 , 1) * ageTrends(a,2) * mult;
    kCin3_Cin2(a + 5 , 1) = kCin3_Cin2(1, 1) * ageTrends(a,3) * mult;
    kCC_Cin3(a + 5 , 1) = kCC_Cin3(1 , 1) * ageTrends(a,4) * mult;
    rNormal_Inf(a + 5, 1) = rNormal_Inf(1 , 1) * ageTrends(a,5) * mult ; 
    kInf_Cin1(a + 5 , 1) = kInf_Cin1(1 , 1) * ageTrends(a,6) * mult;
    kCin1_Cin2(a + 5 , 1) = kCin1_Cin2(1 , 1) * ageTrends(a,7) * mult;
    kCin2_Cin3(a + 5 , 1) = kCin2_Cin3(1 , 1) * ageTrends(a,8) * mult;
    
    % Apply age trends to non-9v HPV transitions 
    kCin1_Inf(a + 5 , 2) = kCin1_Inf(1 , 2) * ageTrends(a,1) * mult;
    kCin2_Cin1(a + 5 , 2) = kCin2_Cin1(1 , 2) * ageTrends(a,2) * mult;
    kCin3_Cin2(a + 5 , 2) = kCin3_Cin2(1, 2) * ageTrends(a,3) * mult; 
    kCC_Cin3(a + 5, 2) = kCC_Cin3(1 , 2) * ageTrends(a,4) * mult;
    rNormal_Inf(a + 5 , 2) = rNormal_Inf(1 , 2) * ageTrends(a,5) * mult ;
    kInf_Cin1(a + 5 , 2) = kInf_Cin1(1 , 2) * ageTrends(a,6) * mult;
    kCin1_Cin2(a + 5 , 2) = kCin1_Cin2(1 , 2) * ageTrends(a,7) * mult;
    kCin2_Cin3(a + 5 , 2) = kCin2_Cin3(1 , 2) * ageTrends(a,8) * mult;
end

% % if calibBool && any(25 == pIdx) % CJB note: old code
% %     idx = find(25 == pIdx);
% %     kProgrsMult = paramSet(paramsSub{idx}.inds(:));
% % else
% %     kProgrsMult = 1.0;
% % end
% % 
% % if calibBool && any(26 == pIdx)
% %     idx = find(26 == pIdx);
% %     kRegrsMult = paramSet(paramsSub{idx}.inds(:));
% % else
% %     kRegrsMult = 1.0;
% % end
% % kCin1_Inf = [kCin1_Inf , kCin1_Inf .* kProgrsMult];
% % kCin2_Cin1 = [kCin2_Cin1 , kCin2_Cin1 .* kProgrsMult];
% % kCin3_Cin2 = [kCin3_Cin2 , kCin3_Cin2 .* kProgrsMult];
% % kCC_Cin3 = [kCC_Cin3 , kCC_Cin3 .* kProgrsMult];
% % rNormal_Inf = [rNormal_Inf , rNormal_Inf .* kRegrsMult];
% % kInf_Cin1 = [kInf_Cin1 , kInf_Cin1 .* kRegrsMult];
% % kCin1_Cin2 = [kCin1_Cin2 , kCin1_Cin2 .* kRegrsMult];
% % kCin2_Cin3 = [kCin2_Cin3 , kCin2_Cin3 .* kRegrsMult];

%% Save HPV natural history parameters
hpvOn = 1; % bool to turn HPV on or off although model not set up for HPV to be off

% Import from Excel CC-associated death rate
% file = [pwd , '/Config/HPV_parameters.xlsx'];
% muCC = xlsread(file , 'Cervical Cancer' , 'B6:D11');
% save(fullfile(paramDir ,'hpvNHParamsFrmExcel'), 'muCC');

% Load pre-saved CC-associated death rate
load([paramDir , 'hpvNHParamsFrmExcel'] , 'muCC' );

% Natural immunity multiplier
if calibBool && any(18 == pIdx)
    idx = find(18 == pIdx);
    lambdaMultImm = zeros(age , 1);
    lambdaMultImmmult = paramSet(paramsSub{idx}.inds(:));
    lambdaMultImm = lambdaMultImm .* lambdaMultImmmult;
else
    lambdaMultImm = zeros(16 , 1);
    lambdaMultImm(1 : 4) = 1 - 0.01;
    lambdaMultImm(5 : 8) = 1 - logspace(log10(0.33) , log10(0.405) , 4);
    lambdaMultImm(9 : 12) = 1 - logspace(log10(0.40) , log10(0.215) , 4);
    lambdaMultImm(13 : 16) = lambdaMultImm(12);
end

% Convert 5-year age groups to 1-year age groups
if ~fivYrAgeGrpsOn
    % Replicate rates across single age groups for other variables
    vars5To1_nms = {'kCin1_Inf' , 'kCin2_Cin1' , 'kCin3_Cin2' , 'kCC_Cin3' , ...
                 'rNormal_Inf' , 'kInf_Cin1' , 'kCin1_Cin2' , 'kCin2_Cin3' , 'lambdaMultImm'};
    vars5To1_vals = {kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , ...
                 rNormal_Inf , kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm};    
    for j = 1 : length(vars5To1_vals)
        valsA1 = age5To1(vars5To1_vals{j});
        assignin('base', vars5To1_nms{j} , valsA1);
    end
end

% HPV clearance multiplier for HIV-positive persons 
rImmuneHiv = [1.4167; 1.5682; 1.9722; 2.8333];

% HPV infection multiplier for HIV-positive persons
hpv_hivMult = [1.78; 1.99; 2.12; 2.32] .* 1.6;

% HPV clearance multiplier for HIV-positive persons
if calibBool && any(14 == pIdx)
    idx = find(14 == pIdx);
    hpv_hivClear = zeros(4 , 1);
    hpv_hivClear(1,1) = paramSet(paramsSub{idx}.inds(1));
    hpv_hivClear(2,1) = hpv_hivClear(1,1)*paramSet(paramsSub{idx}.inds(2));
    hpv_hivClear(3,1) = hpv_hivClear(2,1)*paramSet(paramsSub{idx}.inds(3));
    hpv_hivClear(4,1) = hpv_hivClear(3,1)*paramSet(paramsSub{idx}.inds(4));
else
    hpv_hivClear = [0.6; 0.55; 0.45; 0.30]; %original values [0.60; 0.55; 0.45; 0.30]
end

% CIN2 to CIN3 progression multiplier for HIV-positive women
if calibBool && any(15 == pIdx)
    idx = find(15 == pIdx);
    c3c2Mults = zeros(4 , 1);
    c3c2Mults(4,1) = paramSet(paramsSub{idx}.inds(3));
    c3c2Mults(3,1) = c3c2Mults(4,1)*paramSet(paramsSub{idx}.inds(2));
    c3c2Mults(2,1) = c3c2Mults(3,1)*paramSet(paramsSub{idx}.inds(1));
else
    c3c2Mults = [1.0; 1.1; 1.2; 1.3]; %original values [1.0; 1.8; 2.6; 5.5]
end

% CIN1 to CIN2 progression multiplier for HIV-positive women
if calibBool && any(16 == pIdx)
    idx = find(16 == pIdx);
    c2c1Mults = zeros(4 , 1);
    c2c1Mults(4,1) = paramSet(paramsSub{idx}.inds(3));
    c2c1Mults(3,1) = c3c2Mults(4,1)*paramSet(paramsSub{idx}.inds(2));
    c2c1Mults(2,1) = c3c2Mults(3,1)*paramSet(paramsSub{idx}.inds(1));
else
    c2c1Mults = [1.0; 1.4; 1.8; 2.2]; % orginal values [1.00; 1.75; 2.30; 2.66];
end

% HPV tranmission rates
if calibBool && any(10 == pIdx)
    idx = find(10 == pIdx);
    perPartnerHpv_vax = paramSet(paramsSub{idx}.inds(:));
else
    perPartnerHpv_vax = 0.0078;
end

if calibBool && any(11 == pIdx)
    idx = find(11 == pIdx);
    perPartnerHpv_nonV = paramSet(paramsSub{idx}.inds(:));
else
    perPartnerHpv_nonV = perPartnerHpv_vax .* 1.6;
end

% Decrease HPV transmission rate in women with cervical cancer as a proxy for decreased sexual activity
vagTrans_vax = ones(risk , 1) .* perPartnerHpv_vax; % [risk x 1]
vagTrans_nonV = ones(risk , 1) .* perPartnerHpv_nonV; % [risk x 1]
betaHPV_vax = bsxfun(@times , [1.0 0.5 0.1; 1.0 0.5 0.1; 1.0 0.5 0.1] , vagTrans_vax);
betaHPV_nonV = bsxfun(@times , [1.0 0.5 0.1; 1.0 0.5 0.1; 1.0 0.5 0.1] , vagTrans_nonV);
beta_hpvVax = zeros(gender , age , risk , 3); % age x risk x [normal transmission , reduced transmission (CC-regional / CC-distant) , reduced transmission (late-stage HIV)]
beta_hpvNonVax = zeros(gender , age , risk , 3);
for a = 1 : age % calculate per-partnership probability of HPV transmission
    % VACCINE-TYPE HPV
    % force of infection: females to infect HPV-negative males 
    % affected by betaHPV_F2M_vax, probability of transmission based on cervical cancer status/progression, and number of male acts
    beta_hpvVax(2 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHPV_vax , maleActs(a , :)'));
    % force of infection: males to infect HPV-negative females,
    % affected by betaHPV_M2F_vax, probability of transmission based on cervical cancer status/progression, and number of female acts
    beta_hpvVax(1 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHPV_vax , femaleActs(a , :)'));
    
    % NON-VACCINE-TYPE HPV
    % force of infection: females to infect HPV-negative males 
    % affected by betaHPV_F2M_nonV, probability of transmission based on cervical cancer status/progression, and number of male acts
    beta_hpvNonVax(2 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHPV_nonV , maleActs(a , :)'));
    % force of infection: males to infect HPV-negative females,
    % affected by betaHPV_M2F_nonV, probability of transmission based on cervical cancer status/progression, and number of female acts
    beta_hpvNonVax(1 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHPV_nonV , femaleActs(a , :)'));  
end

beta_hpvVax_mod = zeros(risk , age , viral , endpoints , gender);
beta_hpvNonVax_mod = zeros(risk , age , viral , endpoints , gender);
for v = 1 : viral
    for x = 1 : endpoints
        for g = 1 : gender
            for a = 1 : age
                for r = 1 : risk
                    if (v == 5)
                        beta_hpvVax_mod(r , a , v , x , g) = beta_hpvVax(g , a , r , 3);
                        beta_hpvNonVax_mod(r , a , v , x , g) = beta_hpvNonVax(g , a , r , 3);
                    elseif (x > 1)
                        beta_hpvVax_mod(r , a , v , x , g) = beta_hpvVax(g , a , r , 2);
                        beta_hpvNonVax_mod(r , a , v , x , g) = beta_hpvNonVax(g , a , r , 2);
                    else
                        beta_hpvVax_mod(r , a , v , x , g) = beta_hpvVax(g , a , r , 1);
                        beta_hpvNonVax_mod(r , a , v , x , g) = beta_hpvNonVax(g , a , r , 1);
                    end
                end
            end
        end
    end
end

% Cervical cancer progression
kRL = 0.02;
kDR = 0.025;

% Immunity
rImmune = 0.024; % Clearance rate from Immune to Normal; for HPV16, Johnson (2012)
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity

% Male HPV clearance
if calibBool && any(37 == pIdx)
    idx = find(37 == pIdx);
    maleHpvClearMult = paramSet(paramsSub{idx}.inds(:));
else
    maleHpvClearMult = 0.7;
end

% HIV and ART multipliers
if calibBool && any(21 == pIdx)
    idx = find(21 == pIdx);
    artHpvMult = paramSet(paramsSub{idx}.inds(:));
else
    artHpvMult = 1.0;
end

%% Save intervention parameters

%Import from Excel HIV intervention parameters
file = [pwd , '/Config/HIV_parameters_Kenya.xlsx'];
circProtect = xlsread(file , 'Protection' , 'B18');
condProtect = xlsread(file , 'Protection' , 'B19');
MTCTRate = xlsread(file , 'Disease Data' , 'B6:B8');
artVScov = xlsread(file , 'Protection' , 'A33:C46');    % [years , females , males] 
save(fullfile(paramDir ,'hivIntParamsFrmExcel'), 'circProtect' , ...
    'condProtect' , 'MTCTRate' , 'artVScov');

% Load pre-saved HIV intervention parameters
load([paramDir , 'hivIntParamsFrmExcel'] , 'circProtect' , ...
    'condProtect' , 'MTCTRate' , 'artVScov');

% Protection from circumcision and condoms
circProtect = [[circProtect; 0.3] , [0; 0]];  % HIV protection (changed from 30% to 45%) , HPV protection;  
condProtect = [ones(gender,1).*condProtect , [0; 0]];    % HIV protection , HPV protection

% Condom use
if calibBool && any(5 == pIdx);
    idx = find(5 == pIdx);
    condUse = paramSet(paramsSub{idx}.inds(:));
else
    if fivYrAgeGrpsOn
        condUse = 0.15; % * 0.5;
    else
        condUse = .05; %changed from 20%
    end
end

% Background hysterectomy ********NOT UPDATED!!!!!!!!!!!!!!!!!
hyst = 0; % bool to turn background hysterectomy on or off
OMEGA = zeros(age , 1); % hysterectomy rate

% ART coverage
artOutMult = 1.0; %0.95;
minLim = (0.70/0.81); % minimum ART coverage by age
maxLim = ((1-(0.78/0.81)) + 1); % maximum ART coverage by age, adjust to lower value to compensate for HIV-associated mortality
artYr = [(artVScov(:,1) - 1); (2030 - 1)]; % assuming 90-90-90 target reached by 2030
maxRateM = [artVScov(:,3) ./ 100 ; 0.55] .* artOutMult; % population-level ART coverage in males (72.9% if 90-90-90)
maxRateF = [artVScov(:,2) ./ 100 ; 0.65] .* artOutMult; % population-level ART coverage in females (72.9% if 90-90-90)
artYr_vec = cell(size(artYr , 1) - 1, 1); % save data over time interval in a cell array
artM_vec = cell(size(artYr , 1) - 1, 1);
artF_vec = cell(size(artYr , 1) - 1, 1);
for i = 1 : size(artYr , 1) - 1 % interpolate ART viral suppression coverages at steps within period
    period = [artYr(i) , artYr(i + 1)];
    artYr_vec{i} = interp1(period , artYr(i : i + 1 , 1) , ...
        artYr(i) : timeStep : artYr(i + 1));
    artM_vec{i} = interp1(period , maxRateM(i : i + 1 , 1) , ...
        artYr(i) : timeStep : artYr(i + 1));
    artF_vec{i} = interp1(period , maxRateF(i : i + 1 , 1) , ...
        artYr(i) : timeStep : artYr(i + 1));
end

% Intervention start years
hivStartYear = 1980;
circStartYear = 1980;
circNatStartYear = 2008;
vaxStartYear = 2014;
%%
% VMMC coverage
vmmcYr = [circStartYear; 2003; 2008; 2014; 2030];
circ_aVec = {4 , 5 , 6, [7:8] , [9:10], [11:age]}; % Ages: (15-19), (20-24), (25-29), (30-39), (40-49), (50+)
vmmcRate = [0.0 0.0 0.0 0.0 0 0; ... % 1980
            0.5 0.70 0.75 0.70 0.6 0.50; ... % 2003
            0.60 0.70 0.75 0.70 0.60 0.60; ... %2008
            0.70 0.75 0.80 0.75 0.70 0.60; ... %2014 
            0.70 0.75 0.80 0.75 0.70 0.60];   % 2030 ideal: 90% [year x age group]
vmmcYr_vec = cell(size(vmmcYr , 1) - 1 , 1); % save data over time interval in a cell array
vmmc_vec = cell(size(vmmcYr , 1) - 1 , length(circ_aVec));
for i = 1 : size(vmmcYr , 1) - 1 % interpolate VMMC coverages at steps within period
    period = [vmmcYr(i) , vmmcYr(i + 1)];
    xq = [vmmcYr(i) : timeStep : vmmcYr(i + 1)];
    vmmcYr_vec{i} = interp1(period , vmmcYr(i : i + 1 , 1) , xq);
    for aInd = 1 : length(circ_aVec)
        vmmc_vec{i,aInd} = interp1(period , vmmcRate(i : i + 1 , aInd) , xq);
    end
end
%%
% Vaccination
waning = 0;    % bool to turn waning on or off

% Screening timeframe
screenYrs = [2000; 2003; 2016; currYear; 2023; 2030; 2045];
hpvScreenStartYear = screenYrs(1);

% Screening test sensitivities
cytoSens = [0.0 , 0.57 , 0.57]; % pap smear
hpvSens = [0.0 , 0.881 , 0.881]; % careHPV
hpvSensWHO = [0.0 , 0.90 , 0.94]; % HPV test

% Baseline screening algorithm
baseline.screenCover = [0.0; 0.08; 0.16; 0.16; 0.16; 0.16; 0.16]; %[0.0; 0.08; 0.16; 0.16; 0.16; 0.16; 0.16];
%baseline.diseaseInds = [1 : disease];
baseline.screenAge = [35/max(1 , fivYrAgeGrpsOn*5)+1];
baseline.screenAgeMults = [1.0 / max(1 , fivYrAgeGrpsOn*5)];
baseline.testSens = cytoSens;
% cryoElig = [1.0 , 0.85 , 0.75 , 0.10 , 0.10 , 0.10];
baseline.colpoRetain = 0.361; % 0.721
baseline.cinTreatEff = [0.905 , 0.905 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766]; % cryotherapy/LEEP effectiveness by HIV status
baseline.cinTreatRetain = 0.25; % 0.51
baseline.cinTreatHpvPersist = 0.28; % HPV persistence with LEEP including treatment failure
baseline.cinTreatHpvPersistHivNeg = baseline.cinTreatHpvPersist - (1-baseline.cinTreatEff(1)); % 0.185; proportion of effectively treated HIV-negative women who have persistent HPV after LEEP
baseline.ccTreatRetain = 0.2; % 0.40
baseline.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    baseline.screenCover_vec{i} = interp1(period , baseline.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% CISNET screening algorithm
cisnet.screenCover = [0.0; 0.1; 0.30; 0.30; 0.30; 0.30; 0.30];
cisnet.screenAge = [35/max(1 , fivYrAgeGrpsOn*5)+1];
cisnet.screenAgeMults = [1.0 / max(1 , fivYrAgeGrpsOn*5)];
cisnet.testSens = cytoSens;
cisnet.colpoRetain = 0.361; % (compliance) * (CIN2+/CC correctly identified by same-day colposcopy)
cisnet.cinTreatEff = baseline.cinTreatEff;
cisnet.cinTreatRetain = 0.25;
cisnet.cinTreatHpvPersist = 0.28; % HPV persistence with cryotherapy including treatment failure
cisnet.cinTreatHpvPersistHivNeg = cisnet.cinTreatHpvPersist - (1-cisnet.cinTreatEff(1)); % proportion of effectively treated HIV-negative women who have persistent HPV after cryotherapy
cisnet.ccTreatRetain = 0.2;
cisnet.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    cisnet.screenCover_vec{i} = interp1(period , cisnet.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% WHO screening algorithm - version a
who.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.70; 0.90]; % CJB note: removed 90% screening compliance beginning in current year
who.testSens = hpvSensWHO;
who.colpoRetain = 1.0;
who.cinTreatEff = [1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0];
who.cinTreatRetain = 0.90; % treatment compliance
who.cinTreatHpvPersist = 0.0; % 100% treatment efficacy 
who.cinTreatHpvPersistHivNeg = who.cinTreatHpvPersist - (1-who.cinTreatEff(1)); % proportion of effectively treated HIV-negative women who have persistent HPV after treatment
who.ccTreatRetain = 0.90; % treatment compliance
who.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    who.screenCover_vec{i} = interp1(period , who.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% WHO screening algorithm - version b (to apply WHO screening parameters at different ages by HIV status)
whob.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.70; 0.90]; %CJB note: removed 90% screening compliance beginning in current year
whob.screenAge = [(35/max(1 , fivYrAgeGrpsOn*5)+1) , (45/max(1 , fivYrAgeGrpsOn*5)+1)];
whob.screenAgeMults = [(1.0 / max(1 , fivYrAgeGrpsOn*5)) , (1.0 / max(1 , fivYrAgeGrpsOn*5))];
whob.testSens = hpvSensWHO;
whob.colpoRetain = 1.0;
whob.cinTreatEff = [1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0];
whob.cinTreatRetain = 0.90; % treatment compliance
whob.cinTreatHpvPersist = 0.0; % 100% treatment efficacy
whob.cinTreatHpvPersistHivNeg = whob.cinTreatHpvPersist - (1-whob.cinTreatEff(1)); % proportion of effectively treated HIV-negative women who have persistent HPV after treatment
whob.ccTreatRetain = 0.90; % treatment compliance
whob.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    whob.screenCover_vec{i} = interp1(period , whob.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

%% Import and save calibration data
file = [pwd , '/Config/Calibration_targets_Kenya.xlsx'];

ccInc2012_dObs(: , 1) = xlsread(file , 'Calibration' , 'E10:E22'); % CC Incidence Rate 2012, by age
ccInc2012_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H10:I22');

cc_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'E2:E3'); % CC type distribution for 9v and non-9v HPV types
cc_dist_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H2:I3');

cin3_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'E8:E9'); % CIN3 type distribution for 9v and non-9v HPV types
cin3_dist_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H8:I9');

cin1_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'E6:E7'); % CIN1 type distribution for 9v and non-9v HPV types
cin1_dist_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H6:I7');

hpv_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'E4:E5'); % HPV type distribution for 9v and non-9v HPV types
hpv_dist_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H4:I5');

cinPos2007_dObs(: , 1) = xlsread(file , 'Calibration' , 'E23:E25'); % CIN2/CIN3 Prevalence (HIV+) 2007, by age
cinPos2007_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H23:I25');

cin1_2010_dObs(: , 1) = xlsread(file , 'Calibration' , 'E29:E30'); % CIN1 Prevalence 2010, by HIV status
cin1_2010_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H29:I30');

cin2_2010_dObs(: , 1) = xlsread(file , 'Calibration' , 'E31:E32'); % CIN2+ Prevalence 2010, by HIV status
cin2_2010_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H31:I32');

hpv_hiv_dObs(: , 1) = xlsread(file , 'Calibration' , 'E33:E36'); % HPV Prevalence in high risk HIV+ Women (high risk) 2006, by age
hpv_hiv_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H33:I36');
hpv_hivNeg_dObs(: , 1) = xlsread(file , 'Calibration' , 'E37:E40'); % HPV Prevalence in high risk HIV- Women (high risk) 2006, by age
hpv_hivNeg_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H37:I40');

hpv_all_dObs(: , 1) = xlsread(file , 'Calibration' , 'E41:E46'); % HPV Prevalence in all women, by age
hpv_all_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H41:I46');
hpv_hiv2009_dObs(: , 1) = xlsread(file , 'Calibration' , 'E47:E51'); % HPV Prevalence in HIV+ women, by age
hpv_hiv2009_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H47:I51');

hivPrevF_dObs(: , 1) = xlsread(file , 'Calibration' , 'E52:E72'); % HIV Prevalence in Women 2003, 2007, and 2009, by age
hivPrevF_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H52:I72');
hivPrevM_dObs(: , 1) = xlsread(file , 'Calibration' , 'E73:E96'); % HIV Prevalence in Men 2003, 2007, and 2009, by age
hivPrevM_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H73:I96');
hivPrevAll_dObs(: , 1) = xlsread(file , 'Calibration' , 'E97:E106'); % HIV Prevalence in all in 2012, by age
hivPrevAll_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H97:I106');
 
popAgeDist_dObs(: , 1) = xlsread(file , 'Calibration' , 'E107:E170');  % Population age distribution in 1979, 1989, 1990, and 2009
popAgeDist_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H107:I170'); 

totPopSize_dObs(: , 1) = xlsread(file , 'Calibration' , 'E171:E174'); % Total population size in 1979, 1989, 1990, and 2009
totPopSize_dObs(: , 2:3) = xlsread(file , 'Calibration' , 'H171:I174');
 
save(fullfile(paramDir , 'calibData'), 'ccInc2012_dObs' , 'cc_dist_dObs' , 'cin3_dist_dObs' , ...
    'cin1_dist_dObs' , 'hpv_dist_dObs' , 'cinPos2007_dObs' , 'cin1_2010_dObs' ,'cin2_2010_dObs', ...
    'hpv_hiv_dObs' , 'hpv_hivNeg_dObs' , 'hpv_all_dObs', 'hpv_hiv2009_dObs'  , ...
    'hivPrevF_dObs' , 'hivPrevM_dObs' , 'hivPrevAll_dObs', 'popAgeDist_dObs' , 'totPopSize_dObs')

load([paramDir , 'calibData'], 'ccInc2012_dObs' , 'cc_dist_dObs' , 'cin3_dist_dObs' , ...
    'cin1_dist_dObs' , 'hpv_dist_dObs' , 'cinPos2007_dObs' , 'cin1_2010_dObs' ,'cin2_2010_dObs', ...
    'hpv_hiv_dObs' , 'hpv_hivNeg_dObs' , 'hpv_all_dObs', 'hpv_hiv2009_dObs'  , ...
    'hivPrevF_dObs' , 'hivPrevM_dObs' , 'hivPrevAll_dObs', 'popAgeDist_dObs' , 'totPopSize_dObs');



%% Load indices *****************************************************************************************************************************************************************************
% disp('Preparing indices...')
% disp('This may take a while...')

%% mixInfect.m indices
gar = zeros(gender , age , risk , disease * viral * hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            gar(g , a , r , :) = sort(toInd(allcomb(1 : disease , 1 : viral ,...
                1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , r)));
        end
    end
end

vaxInds = zeros(viral , endpoints , gender , age , risk , disease * 5 * hpvNonVaxStates * intervens); % 5 HPV+ states
nonVInds = zeros(viral , endpoints , gender , age , risk , disease * hpvVaxStates * 5 * intervens);
for v = 1 : viral
    for x = 1 : endpoints
        for g = 1 : gender
            for a = 1 : age
                for r = 1 : risk
                    vaxInds(v , x , g , a , r , :) = ...
                        sort(toInd(allcomb(1 : disease , v , 2 : 6 , ...
                        1 : hpvNonVaxStates , x , 1 : intervens , g , a , r)));
                    nonVInds(v , x , g , a , r , :) = ...
                        sort(toInd(allcomb(1 : disease , v , 1 : hpvVaxStates , ...
                        2 : 6 , x , 1 : intervens , g , a , r)));
                end
            end
        end
    end
end

hpvVaxSus = zeros(disease, gender , age , risk , intervens , viral*hpvNonVaxStates*3);
hpvVaxImm = hpvVaxSus;
hpvVaxInf = hpvVaxSus;
hpvNonVaxSus = zeros(disease , gender , age , risk , intervens , viral*hpvVaxStates*3);
hpvNonVaxImm = hpvNonVaxSus;
hpvNonVaxInf = hpvNonVaxSus;
for d = 1 : disease
    for g = 1 : gender
        for a = 1 : age
            for r = 1 : risk
                for p = 1 : intervens
                    hpvVaxSus(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 : 3 , p , g , a , r)));
                    hpvVaxImm(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 : 3 , p , g , a , r)));
                    hpvVaxInf(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 2 , 1 : hpvNonVaxStates , 1 : 3 , p , g , a , r)));

                    hpvNonVaxSus(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 , 1 : 3 , p , g , a , r)));
                    hpvNonVaxImm(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 7 , 1 : 3 , p , g , a , r)));
                    hpvNonVaxInf(d , g , a , r , p , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 2 , 1 : 3 , p , g , a , r)));
                end
            end
        end
    end
end

hivCurr = zeros(viral , endpoints , gender , age , risk , 6 * hpvVaxStates * hpvNonVaxStates * intervens); % 6 HIV+ disease states
for v = 1 : viral
    for x = 1 : endpoints
        for g = 1 : gender
            for a = 1 : age
                for r = 1 : risk
                    hivCurr(v , x , g , a , r , :) = toInd(allcomb(3 : 8 , v , 1 : hpvVaxStates , ...
                        1 : hpvNonVaxStates , x , 1 : intervens , g , a , r));
                end
            end
        end
    end
end
 
hivSus = zeros(2 , hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk , intervens);
for d = 1 : 2
    for h = 1 : hpvVaxStates
        for s = 1 : hpvNonVaxStates
            for x = 1 : endpoints
                for g = 1 : gender
                    for a = 1 : age
                        for r = 1 : risk
                            hivSus(d , h , s , x , g , a , r , :) = ...
                                sort(toInd(allcomb(d , 1 , h , s , x , 1 : intervens , g , a , r)));
                        end
                    end
                end
            end
        end
    end
end

toHiv = zeros(hpvVaxStates , hpvNonVaxStates , endpoints , gender , age , risk , intervens);
for h = 1 : hpvVaxStates
    for s = 1 : hpvNonVaxStates
        for x = 1 : endpoints
            for g = 1 : gender
                for a = 1 : age
                    for r = 1 : risk
                        toHiv(h , s , x , g , a , r , :) = ...
                            sort(toInd(allcomb(3 , 1 , h , s , x , 1 : intervens , g , a , r)));
                    end
                end
            end
        end
    end
end
% disp('mixInfect indices loaded')

%% hivNH.m indices
hivInds = zeros(disease , viral , gender , age , risk , hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
for d = 1 : disease
    for v = 1 : viral
%         for h = 1 : hpvVaxStates
%             for s = 1 : hpvNonVaxStates
%                 for x = 1 : endpoints
                    for g = 1 : gender
                        for a = 1 : age
                            for r = 1 :risk
                                hivInds(d , v , g , a , r , :) = ...
                                    sort(toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , r)));
                            end
                        end
                    end
%                 end
%             end
%         end
    end
end
% disp('hivNH indices loaded')

%% hpvCCNH.m indices
% disp('Preparing indices for HPV modules...')

cin3hpvVaxIndsFrom = zeros(disease , hpvNonVaxStates , age , viral * intervens * risk * 3);
ccLochpvVaxIndsTo = cin3hpvVaxIndsFrom;
ccLochpvVaxIndsFrom = zeros(disease , hpvNonVaxStates , age , viral * intervens * risk);
ccReghpvVaxInds = ccLochpvVaxIndsFrom;
ccDisthpvVaxInds = ccLochpvVaxIndsFrom;
cin3hpvNonVaxIndsFrom = zeros(disease , hpvVaxStates , age , viral * intervens * risk * 3);
ccLochpvNonVaxIndsTo = cin3hpvNonVaxIndsFrom;
ccLochpvNonVaxIndsFrom = zeros(disease , hpvVaxStates , age , viral * intervens * risk);
ccReghpvNonVaxInds = ccLochpvNonVaxIndsFrom;
ccDisthpvNonVaxInds = ccLochpvNonVaxIndsFrom;
cin1hpvVaxInds = zeros(disease , age , viral * hpvNonVaxStates * intervens * risk * 3);
cin2hpvVaxInds = cin1hpvVaxInds;
cin3hpvVaxInds = cin1hpvVaxInds;
cin1hpvNonVaxInds = zeros(disease , age , viral * hpvVaxStates * intervens * risk * 3);
cin2hpvNonVaxInds = cin1hpvNonVaxInds;
cin3hpvNonVaxInds = cin1hpvNonVaxInds;
for d = 1 : disease
    for a = 1 : age
        for s = 1 : hpvNonVaxStates
            cin3hpvVaxIndsFrom(d , s , a , :) = sort(toInd(allcomb(d , 1 : viral , 5 , s , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
            ccLochpvVaxIndsTo(d , s , a , :) = sort(toInd(allcomb(d , 1 : viral , 6 , s , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
            ccLochpvVaxIndsFrom(d , s , a , :) = sort(toInd(allcomb(d , 1 : viral , 6 , s , 1 , 1 : intervens , 2 , a , 1 : risk)));
            ccReghpvVaxInds(d , s , a , :) = sort(toInd(allcomb(d , 1 : viral , 6 , s , 2 , 1 : intervens , 2 , a , 1 : risk)));
            ccDisthpvVaxInds(d , s , a , :) = sort(toInd(allcomb(d , 1 : viral , 6 , s , 3 , 1 : intervens , 2 , a , 1 : risk)));
        end
        for h = 1 : hpvVaxStates
            cin3hpvNonVaxIndsFrom(d , h , a , :) = sort(toInd(allcomb(d , 1 : viral , h , 5 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
            ccLochpvNonVaxIndsTo(d , h , a , :) = sort(toInd(allcomb(d , 1 : viral , h , 6 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
            ccLochpvNonVaxIndsFrom(d , h , a , :) = sort(toInd(allcomb(d , 1 : viral , h , 6 , 1 , 1 : intervens , 2 , a , 1 : risk)));
            ccReghpvNonVaxInds(d , h , a , :) = sort(toInd(allcomb(d , 1 : viral , h , 6 , 2 , 1 : intervens , 2 , a , 1 : risk)));
            ccDisthpvNonVaxInds(d , h , a , :) = sort(toInd(allcomb(d , 1 : viral , h , 6 , 3 , 1 : intervens , 2 , a , 1 : risk)));
        end
        cin1hpvVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 3 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
        cin2hpvVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 4 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
        cin3hpvVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 5 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
        cin1hpvNonVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 3 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
        cin2hpvNonVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 4 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
        cin3hpvNonVaxInds(d , a , :) = sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 5 , 1 : 3 , 1 : intervens , 2 , a , 1 : risk)));
    end
end

normalhpvVaxInds = zeros(disease , gender , age , viral * hpvNonVaxStates * intervens * risk * 3);
immunehpvVaxInds = normalhpvVaxInds;
infhpvVaxInds = normalhpvVaxInds;
normalhpvNonVaxInds = zeros(disease , gender , age , viral * hpvVaxStates * intervens * risk * 3);
immunehpvNonVaxInds = normalhpvNonVaxInds;
infhpvNonVaxInds = normalhpvNonVaxInds;
for g = 1 : gender
    for d = 1 : disease
        for a = 1 : age
                normalhpvVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
                normalhpvNonVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
                
                immunehpvVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 7 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
                immunehpvNonVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 7 , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
                
                infhpvVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 2 , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
                infhpvNonVaxInds(d , g , a , :) = ...
                    sort(toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 2 , 1 : 3 , 1 : intervens , g , a , 1 : risk)));
        end
    end      
end

% BACKGROUND HYSTERECTOMY NOT UPDATED!!!!!!!!!!!!!!!!!
% hystSusInds = zeros(disease , 9 , age , viral * hpvVaxStates * intervens * risk);
% hystInds = zeros(disease , age , viral * hpvVaxStates * intervens * risk);
% for a = 1 : age
%     for d = 1 : disease
%         for h = 1 : hpvVaxStates
%             hysthpvVaxSusInds(d , h , a , :) = toInd(allcomb(d , 1 : viral , h , ...
%                 1 : hpvNonVaxStates , 1 : intervens , 2 , a , 1 : risk));
%             hysthpvVaxInds(d , a , :) = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , ...
%                 8 , 4 , 1 : intervens , 2 , a , 1 : risk));
%         end
%         for s = 1 : hpvNonVaxStates
%             hysthpvNonVaxSusInds(d , s , a , :) = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , ...
%                 s , 1 : intervens , 2 , a , 1 : risk));
%             hysthpvNonVaxInds(d , a , :) = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , ...
%                 8 , 1 : intervens , 2 , a , 1 : risk));
%         end
%     end
% end

% disp('hpvCCNH indices loaded')

%% ageRisk.m indices
ageInd = zeros(gender , age , disease * viral * hpvVaxStates * hpvNonVaxStates * endpoints * intervens * risk);
riskInd = zeros(gender , age , risk , disease * viral * hpvVaxStates * hpvNonVaxStates * endpoints * intervens);
for g = 1 : gender
    for a = 1 : age
        ageInd(g , a , :) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , ...
            1 : hpvNonVaxStates, 1 : endpoints , 1 : intervens , g , a , 1 : risk));   
        for r = 1 : risk
            riskInd(g , a , r , :) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , ...
                1 : hpvNonVaxStates, 1 : endpoints , 1 : intervens , g , a , r));
        end
    end
end

%% vmmc.m indices
hivNegNonVMMCinds = zeros(age , hpvVaxStates * hpvNonVaxStates * intervens  * risk);
hivNegVMMCinds = hivNegNonVMMCinds;
for a = 1 : age
    hivNegNonVMMCinds(a , :) = toInd(allcomb(1 , 1 , 1 : hpvVaxStates , ...
        1 : hpvNonVaxStates , 1 , 1 : intervens , 1 , a , 1 : risk));
    hivNegVMMCinds(a , :) = toInd(allcomb(2 , 1 , 1 : hpvVaxStates , ...
        1 : hpvNonVaxStates , 1 , 1 : intervens , 1 , a , 1 : risk));
end


%% Make matrices ******************************************************************************************************************************************************************************
pop = spalloc(prod(dim) , 1 , prod(dim));

%% Viral load progression
% disp('Building viral load progression matrix')

xInds = [];
yInds = [];
vals = [];
for g = 1 : gender
    for a = 1 : age
        for d = 3 : 7
            % for v = 1; vlAcute
            vlAcute = toInd(allcomb(d , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , 1 : risk));
            xInds = [xInds; vlAcute];
            yInds = [yInds; vlAcute];
            vals = [vals; ones(length(vlAcute),1) .* ( -kVl(a , 1 , g) )];

            % for v = 2 : 4; Asymptomatic -> Pre-AIDS -> AIDS
            for v = 2 : 4
                vlCurr = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , 1 : risk));
                vlPrev = toInd(allcomb(d , (v - 1) , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , 1 : risk));
                xInds = [xInds; vlCurr; vlCurr]; % prev -> curr
                yInds = [yInds; vlPrev; vlCurr]; % curr -> (next)
                vals = [vals; ones(length(vlCurr),1) .* ( kVl(a , v-1 , g) ); ones(length(vlPrev),1) .* ( -kVl(a , v , g) )];
            end

            % for v = 5; Late-stage
            vl_50k = toInd(allcomb(d , 5 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , 1 : risk));
            vl_10kto50k = toInd(allcomb(d , 4 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , g , a , 1 : risk));
            xInds = [xInds; vl_50k];
            yInds = [yInds; vl_10kto50k];
            vals = [vals; ones(length(vl_50k),1) .* ( kVl(a , 4 , g) )];
        end
    end
end
vlAdvancer = sparse(xInds , yInds , vals , numel(pop) , numel(pop));
% disp('Finished building viral load progression matrix')

%% Fertility prior to 1990

% birth indices
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));
posMaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix for uninfected mothers
% disp('Building fertility matrix for uninfected mothers')
xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    hivUninf = toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    hivPosArt = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    xInds = [xInds; ones(length(hivUninf),1).*negMaleBirth; ones(length(hivUninf),1).*negFemaleBirth; ...
        ones(length(hivPosArt),1).*negMaleBirth; ones(length(hivPosArt),1).*negFemaleBirth];
    yInds = [yInds; hivUninf; hivUninf; hivPosArt; hivPosArt];
    vals = [vals; ones((length(hivUninf)*2+length(hivPosArt)*2),1) .* ( 0.5*fertility(a,1) )];
end
fertMat = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

% fertility matrix for infected mothers
% disp('Building fertility matrix for HIV-infected mothers')
xIndsPos = [];
yIndsPos = [];
valsPos = [];
xIndsNeg = [];
yIndsNeg = [];
valsNeg = [];
for d = 3 : 7 % hiv infected
    for v = 1 : viral % hiv infected
        for a = 1 : age
            hivInfected = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk));
            xIndsPos = [xIndsPos; ones(length(hivInfected),1).*posMaleBirth; ones(length(hivInfected),1).*posFemaleBirth];
            yIndsPos = [yIndsPos; hivInfected; hivInfected];
            valsPos = [valsPos; ones((length(hivInfected)*2),1) .* ( 0.5*fertility(a,d-1) )];   
            xIndsNeg = [xIndsNeg; ones(length(hivInfected),1).*negMaleBirth; ones(length(hivInfected),1).*negFemaleBirth];
            yIndsNeg = [yIndsNeg; hivInfected; hivInfected];
            valsNeg = [valsNeg; ones((length(hivInfected)*2),1) .* ( 0.5*fertility(a,d-1) )];   
        end
    end
end
hivFertPosBirth = sparse(xIndsPos , yIndsPos , valsPos , numel(pop) , numel(pop));
hivFertNegBirth = sparse(xIndsNeg , yIndsNeg , valsNeg , numel(pop) , numel(pop));

%% Fertility by 2010

% birth indices
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));
posMaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix for uninfected mothers
% disp('Building fertility2 matrix for uninfected mothers')
xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    hivUninf = toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    hivPosArt = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    xInds = [xInds; ones(length(hivUninf),1).*negMaleBirth; ones(length(hivUninf),1).*negFemaleBirth; ...
        ones(length(hivPosArt),1).*negMaleBirth; ones(length(hivPosArt),1).*negFemaleBirth];
    yInds = [yInds; hivUninf; hivUninf; hivPosArt; hivPosArt];
    vals = [vals; ones((length(hivUninf)*2+length(hivPosArt)*2),1) .* ( 0.5*fertility2(a,1) )];
end
fertMat2 = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

% fertility matrix for infected mothers
% disp('Building fertility2 matrix for HIV-infected mothers')
xIndsPos = [];
yIndsPos = [];
valsPos = [];
xIndsNeg = [];
yIndsNeg = [];
valsNeg = [];
for d = 3 : 7 % hiv infected
    for v = 1 : viral % hiv infected
        for a = 1 : age
            hivInfected = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk));
            xIndsPos = [xInds; ones(length(hivInfected),1).*posMaleBirth; ones(length(hivInfected),1).*posFemaleBirth];
            yIndsPos = [yInds; hivInfected; hivInfected];
            valsPos = [vals; ones((length(hivInfected)*2),1) .* ( 0.5*fertility2(a,d-1) )];   
            xIndsNeg = [xInds; ones(length(hivInfected),1).*negMaleBirth; ones(length(hivInfected),1).*negFemaleBirth];
            yIndsNeg = [yInds; hivInfected; hivInfected];
            valsNeg = [vals; ones((length(hivInfected)*2),1) .* ( 0.5*fertility2(a,d-1) )];   
        end
    end
end
hivFertPosBirth2 = sparse(xIndsPos , yIndsPos , valsPos , numel(pop) , numel(pop));
hivFertNegBirth2 = sparse(xIndsNeg , yIndsNeg , valsNeg , numel(pop) , numel(pop));

%% Fertility by 2020

% birth indices
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));
posMaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix for uninfected mothers
% disp('Building fertility3 matrix for uninfected mothers')
xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    hivUninf = toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    hivPosArt = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    xInds = [xInds; ones(length(hivUninf),1).*negMaleBirth; ones(length(hivUninf),1).*negFemaleBirth; ...
        ones(length(hivPosArt),1).*negMaleBirth; ones(length(hivPosArt),1).*negFemaleBirth];
    yInds = [yInds; hivUninf; hivUninf; hivPosArt; hivPosArt];
    vals = [vals; ones((length(hivUninf)*2+length(hivPosArt)*2),1) .* ( 0.5*fertility3(a,1) )];
end
fertMat3 = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

% fertility matrix for infected mothers
% disp('Building fertility3 matrix for HIV-infected mothers')
xIndsPos = [];
yIndsPos = [];
valsPos = [];
xIndsNeg = [];
yIndsNeg = [];
valsNeg = [];
for d = 3 : 7 % hiv infected
    for v = 1 : viral % hiv infected
        for a = 1 : age
            hivInfected = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk));
            xIndsPos = [xInds; ones(length(hivInfected),1).*posMaleBirth; ones(length(hivInfected),1).*posFemaleBirth];
            yIndsPos = [yInds; hivInfected; hivInfected];
            valsPos = [vals; ones((length(hivInfected)*2),1) .* ( 0.5*fertility3(a,d-1) )];   
            xIndsNeg = [xInds; ones(length(hivInfected),1).*negMaleBirth; ones(length(hivInfected),1).*negFemaleBirth];
            yIndsNeg = [yInds; hivInfected; hivInfected];
            valsNeg = [vals; ones((length(hivInfected)*2),1) .* ( 0.5*fertility3(a,d-1) )];   
        end
    end
end
hivFertPosBirth3 = sparse(xIndsPos , yIndsPos , valsPos , numel(pop) , numel(pop));
hivFertNegBirth3 = sparse(xIndsNeg , yIndsNeg , valsNeg , numel(pop) , numel(pop));

%% Fertility between 2020-2070

% birth indices
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));
posMaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix for uninfected mothers
% disp('Building fertility3 matrix for uninfected mothers')
xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    hivUninf = toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    hivPosArt = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
        1 : intervens , 2 , a , 1 : risk));
    xInds = [xInds; ones(length(hivUninf),1).*negMaleBirth; ones(length(hivUninf),1).*negFemaleBirth; ...
        ones(length(hivPosArt),1).*negMaleBirth; ones(length(hivPosArt),1).*negFemaleBirth];
    yInds = [yInds; hivUninf; hivUninf; hivPosArt; hivPosArt];
    vals = [vals; ones((length(hivUninf)*2+length(hivPosArt)*2),1) .* ( 0.5*fertility4(a,1) )];
%     hivUninf = toInd(allcomb(1 , 1 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , ...
%         1 : intervens , 2 , a , 1 : risk));
%     xInds = [xInds; ones(length(hivUninf),1).*negMaleBirth; ones(length(hivUninf),1).*negFemaleBirth];
%     yInds = [yInds; hivUninf; hivUninf];
%     vals = [vals; ones((length(hivUninf)*2),1) .* ( 0.5*fertility3(a,1) )];
end
fertMat4 = sparse(xInds , yInds , vals , numel(pop) , numel(pop));

% fertility matrix for infected mothers
% disp('Building fertility3 matrix for HIV-infected mothers')
xIndsPos = [];
yIndsPos = [];
valsPos = [];
xIndsNeg = [];
yIndsNeg = [];
valsNeg = [];
for d = 3 :7 % hiv infected
    for v = 1 : viral % hiv infected
        for a = 1 : age
            hivInfected = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : 3 , 1 : intervens , 2 , a , 1 : risk));
            xIndsPos = [xInds; ones(length(hivInfected),1).*posMaleBirth; ones(length(hivInfected),1).*posFemaleBirth];
            yIndsPos = [yInds; hivInfected; hivInfected];
            valsPos = [vals; ones((length(hivInfected)*2),1) .* ( 0.5*fertility4(a,d-1) )];   
            xIndsNeg = [xInds; ones(length(hivInfected),1).*negMaleBirth; ones(length(hivInfected),1).*negFemaleBirth];
            yIndsNeg = [yInds; hivInfected; hivInfected];
            valsNeg = [vals; ones((length(hivInfected)*2),1) .* ( 0.5*fertility4(a,d-1) )];   
        end
    end
end
hivFertPosBirth4 = sparse(xIndsPos , yIndsPos , valsPos , numel(pop) , numel(pop));
hivFertNegBirth4 = sparse(xIndsNeg , yIndsNeg , valsNeg , numel(pop) , numel(pop));

%% Fertility scale-down
dFertPos1 = (hivFertPosBirth2 - hivFertPosBirth) ./ ((1990 - 1970) * stepsPerYear);
dFertNeg1 = (hivFertNegBirth2 - hivFertNegBirth) ./ ((1990 - 1970) * stepsPerYear); 
dFertMat1 = (fertMat2 - fertMat) ./ ((1990 - 1970) * stepsPerYear);

dFertPos2 = (hivFertPosBirth3 - hivFertPosBirth2) ./ ((2020 - 1990) * stepsPerYear);
dFertNeg2 = (hivFertNegBirth3 - hivFertNegBirth2) ./ ((2020 - 1990) * stepsPerYear);
dFertMat2 = (fertMat3 - fertMat2) ./ ((2020 - 1990) * stepsPerYear);

dFertPos3 = (hivFertPosBirth4 - hivFertPosBirth3) ./ ((2070 - 2020) * stepsPerYear);
dFertNeg3 = (hivFertNegBirth4 - hivFertNegBirth3) ./ ((2070 - 2020) * stepsPerYear);
dFertMat3 = (fertMat4 - fertMat3) ./ ((2070 - 2020) * stepsPerYear);

%% partnersM multiplier 
partnersMmult = ones(1, 3); %[1.2 2.4 1.1];

if calibBool && any(40 == pIdx);
    idx = find(40 == pIdx);
    ptMult = paramSet(paramsSub{idx}.inds(:));
else
    ptMult = [1.2 2.4 1.1 1 1.4 1];
end

d_partnersMmult = ones(6, 60);
d_partnersMmult(1, 1:30) = linspace(1, ptMult(1), 5 * stepsPerYear); % multiplier for increasing pts in M aged 15 - 24
d_partnersMmult(2, 1:30) = linspace(1, ptMult(2), 5 * stepsPerYear);  % multiplier for increasing pts in F aged 15-24
d_partnersMmult(3, 1:30) = linspace(1, ptMult(3), 5 * stepsPerYear); % mulitplier for increasing pts in M and F aged 25 - 49
d_partnersMmult(4, 1:60) = linspace(ptMult(1), ptMult(4), 10 * stepsPerYear); % multiplier for decreasing pts in M aged 15 - 24
d_partnersMmult(5, 1:60) = linspace(ptMult(2), ptMult(5), 10 * stepsPerYear); % multiplier for decreasing pts in F aged 15 - 24
d_partnersMmult(6, 1:60) = linspace(ptMult(3), ptMult(6), 10 * stepsPerYear); % mulitplier for decreasing pts in M and F aged 25 - 49
%d_partnersMmult(5, 1:5) =-logspace(log10(1.2), log10(0.25), 5);

%% risk adjustment multiplier
riskAdj = 0.005;
d_riskAdj = (0 - riskAdj) ./ ((1994 - 1990) .* stepsPerYear);

%% Background death rate before 1950
% disp('Building death matrix')

xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    males = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    females = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    
    xInds = [xInds; males; females];
    yInds = [yInds; males; females];
    vals = [vals; ones(length(males),1) .* ( -mue(a,1) ); ones(length(females),1) .* ( -mue(a,2) )];
end
deathMat = sparse(xInds , yInds , vals , numel(pop) , numel(pop));
% disp('Death matrix complete')

%% Background death rate by 1985
% disp('Building second death matrix')

xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    males = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    females = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    
    xInds = [xInds; males; females];
    yInds = [yInds; males; females];
    vals = [vals; ones(length(males),1) .* ( -mue2(a,1) ); ones(length(females),1) .* ( -mue2(a,2) )];
end
deathMat2 = sparse(xInds , yInds , vals , numel(pop) , numel(pop));
% disp('Second death matrix complete')

%% Background death rate by 2000
% disp('Building third death matrix')

xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    males = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    females = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    
    xInds = [xInds; males; females];
    yInds = [yInds; males; females];
    vals = [vals; ones(length(males),1) .* ( -mue3(a,1) ); ones(length(females),1) .* ( -mue3(a,2) )];
end
deathMat3 = sparse(xInds , yInds , vals , numel(pop) , numel(pop));
% disp('Third death matrix complete')

%% Background death rate by 2020
% disp('Building fourth death matrix')

xInds = [];
yInds = [];
vals = [];
for a = 1 : age
    males = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    females = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , 1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    
    xInds = [xInds; males; females];
    yInds = [yInds; males; females];
    vals = [vals; ones(length(males),1) .* ( -mue4(a,1) ); ones(length(females),1) .* ( -mue4(a,2) )];
end
deathMat4 = sparse(xInds , yInds , vals , numel(pop) , numel(pop));
% disp('Fourth death matrix complete')

%% Mortality matrix scale-down
dDeathMat = (deathMat2 - deathMat) ./ ((1985 - 1950) * stepsPerYear);
dDeathMat2 = (deathMat3 - deathMat2) ./ ((2000 - 1985) * stepsPerYear);
dDeathMat3 = (deathMat4 - deathMat3) ./ ((2020 - 2000) * stepsPerYear);

%% Mortality rate scale-down during the years of ART for use in calculation of HIV-associated mortality on ART
dMue = (mue4 - mue3) ./ ((2020 - 2000) * stepsPerYear);


toc
