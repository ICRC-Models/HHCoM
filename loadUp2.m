% Load and save parameters, calibration data, indices, and matrices
function[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd , ...
    annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4 , ...
    mue , mue2 , mue3 , mue4 , epsA_vec , epsR_vec , ...
    yr , ...
    hivOn , betaHIV_mod , muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , c2c3Mults , c1c2Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , waning , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear , circNatStartYear , vaxStartYear , ...
    baseline , cisnet , who , whob , circProtect , condProtect , MTCTRate , ...
    hyst , OMEGA , ...
    ccInc2012_dObs , ccInc2018_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    cinPos2015_dObs , cinNeg2015_dObs , hpv_hiv_dObs , hpv_hivNeg_dObs , ...
    hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
    popAgeDist_dObs , totPopSize_dObs , ...
    hivCurr , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , ...
    hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , ageInd , riskInd , ...
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3 , dFertMat3 , deathMat , deathMat2 , deathMat3 , deathMat4 , ...
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
% background mortality, fertility, partnerships, and acts per partnership
% file = [pwd , '/Config/Population_data.xlsx'];
% popInit = xlsread(file , 'Demographics' , 'D6:E21');
% riskDistM = xlsread(file , 'Demographics' , 'B54:D69');
% mue = zeros(age , gender);
% mue(: , 1) = xlsread(file , 'Demographics' , 'AX84:AX99'); %'B84:B99'); % 1950
% mue(: , 2) = xlsread(file , 'Demographics' , 'BE84:BE99'); %'O84:O99');
% mue2 = zeros(age , gender);
% mue2(: , 1) = xlsread(file , 'Demographics' , 'AY84:AY99'); %'I84:I99'); % 1985
% mue2(: , 2) = xlsread(file , 'Demographics' , 'BF84:BF99'); %'V84:V99');
% mue3 = zeros(age , gender);
% mue3(: , 1) = xlsread(file , 'Demographics' , 'AZ84:AZ99'); %'J84:J99'); % 2000
% mue3(: , 2) = xlsread(file , 'Demographics' , 'BG84:BG99'); %'W84:W99');
% mue4 = zeros(age , gender);
% mue4(: , 1) = xlsread(file , 'Demographics' , 'BB84:BB99'); %'L84:L99'); % 2020
% mue4(: , 2) = xlsread(file , 'Demographics' , 'BI84:BI99'); %'Y84:Y99');
% fertility = xlsread(file , 'Demographics' , 'B115:G130');
% partnersM = xlsread(file , 'Demographics' , 'B163:D178');
% partnersF = xlsread(file , 'Demographics' , 'E163:G178');
% maleActs = xlsread(file , 'Demographics' , 'B202:D217');
% femaleActs = xlsread(file , 'Demographics' , 'B222:D237');
% save(fullfile(paramDir ,'demoParamsFrmExcel'), 'popInit' , 'riskDistM' , ...
%     'mue' , 'mue2' , 'mue3' , 'mue4' , 'fertility' , 'partnersM' , 'partnersF' , 'maleActs' , 'femaleActs');

% Load pre-saved initial population size by age and gender, male risk distribution by age, 
% background mortality by age and gender, and fertility by age and gender
load([paramDir , 'demoParamsFrmExcel'] , 'popInit' , 'riskDistM' , 'mue' , ...
    'mue2' , 'mue3' , 'mue4' , 'fertility');

% Set female risk distribution
riskDistF = riskDistM;

% Calculate fertility2 and fertility3 matrices
if calibBool && any(36 == pIdx);
    idx = find(36 == pIdx);
    fertDeclineProp = paramSet(paramsSub{idx}.inds(:));
else
    fertDeclineProp = [0.37 ; 0.75];
end
fertility2 = fertility .* fertDeclineProp(1,1);
fertility3 = fertility2 .* fertDeclineProp(2,1);
fertility4 = fertility3 .* 0.60;

% Male partners per year by age and risk group
if calibBool && any(1 == pIdx)
    idx = find(1 == pIdx);
    partnersM = zeros(age , risk);
    %partnersMmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    partnersM(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
    %partnersM(3:10, 1:risk) = [paramSet(rl).*paramSet(rm).*paramSet(rh) , paramSet(rm).*paramSet(rh) , paramSet(rh)];
    %partnersM(11:age , 1:risk) = ones(6,risk).*partnersM(10 , 1:risk);
    %partnersM(3:age , 1:risk) = partnersM(3:age , 1:risk) .* partnersMmult;
    %partnersM(10:age , 3) = ones(7 , 1);
    partnersM(4:6 , 3) = paramSet(paramsSub{idx}.inds(1:3));
    partnersM(3 , 3) = partnersM(4 , 3) .* 0.5;
    partnersM(7:9 , 3) = ones(3,1).*paramSet(paramsSub{idx}.inds(4));
    partnersM(10:age , 3) = ones(7,1).*paramSet(paramsSub{idx}.inds(5));
    partnersM(3 , 2) = paramSet(paramsSub{idx}.inds(6))*partnersM(3 , 3);
    partnersM(3 , 1) = paramSet(paramsSub{idx}.inds(7))*partnersM(3 , 2);
    for a = 4 : age
        partnersM(a , 1:2) = (partnersM(a , 3)/partnersM(a-1 , 3)) .* partnersM(a-1 , 1:2);
    end
else 
    load([paramDir , 'demoParamsFrmExcel'] , 'partnersM');
end

% Female partners per year by age and risk group
if calibBool && any(2 == pIdx)
    idx = find(2 == pIdx);
    partnersF = zeros(age , risk);
    %partnersFmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    partnersF(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
    %partnersF(3:10 , 1:risk) = [paramSet(rl).*paramSet(rm).*paramSet(rh) , paramSet(rm).*paramSet(rh) , paramSet(rh)];
    %partnersF(11:age , 1:risk) = ones(6,risk).*partnersF(10 , 1:risk);
    %partnersF(3:age , 1:risk) = partnersF(3:age , 1:risk) .* partnersFmult;
    %partnersF(10:age , 3) = ones(7 , 1);
    partnersF(4:6 , 3) = paramSet(paramsSub{idx}.inds(1:3));
    partnersF(3 , 3) = partnersF(4 , 3) .* 0.5;
    partnersF(7:9 , 3) = ones(3,1).*paramSet(paramsSub{idx}.inds(4));
    partnersF(10:age , 3) = ones(7,1).*paramSet(paramsSub{idx}.inds(5));
    partnersF(3 , 2) = paramSet(paramsSub{idx}.inds(6))*partnersF(3 , 3);
    partnersF(3 , 1) = paramSet(paramsSub{idx}.inds(7))*partnersF(3 , 2);
    for a = 4 : age
        partnersF(a , 1:2) = (partnersF(a , 3)/partnersF(a-1 , 3)) .* partnersF(a-1 , 1:2);
    end
else
    load([paramDir , 'demoParamsFrmExcel'] , 'partnersF');
end    

% Female acts per partnership per year by age and risk group
if calibBool && any(9 == pIdx)
    idx = find(9 == pIdx);
    femaleActs = zeros(age , risk);
    %femaleActsmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    femaleActs(1:2 , 1:risk) = zeros(2 , risk);
    %femaleActs(3:age , 1:risk) = [paramSet(rl) , paramSet(rm).*paramSet(rl) , paramSet(rh).*paramSet(rm).*paramSet(rl)];
    %femaleActs(3:age , 1:risk) = femaleActs(3:age , 1: risk) .* femaleActsmult;
    femaleActs(4:6 , 1) = paramSet(paramsSub{idx}.inds(1:3));
    femaleActs(3 , 1) = femaleActs(4 , 1) .* 0.5;
    femaleActs(7:9 , 1) = ones(3,1).*paramSet(paramsSub{idx}.inds(4));
    femaleActs(10:age , 1) = ones(7,1).*paramSet(paramsSub{idx}.inds(5));
    femaleActs(3 , 2) = 0.6 .* femaleActs(3 , 1);
    femaleActs(3 , 3) = 0.6 .* femaleActs(3 , 2);
    for a = 4 : age
        femaleActs(a , 2:3) = (femaleActs(a , 1)/femaleActs(a-1 , 1)) .* femaleActs(a-1 , 2:3);
    end
else
    load([paramDir , 'demoParamsFrmExcel'] , 'femaleActs');
end

% Male acts per partnership per year by age and risk group
if calibBool && any(8 == pIdx)
    idx = find(8 == pIdx);
    %maleActsmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    %maleActs(1:2 , 1:risk) = zeros(2 , risk);
    %maleActs(3:age , 1:risk) = [paramSet(rl) , paramSet(rm).*paramSet(rl) , paramSet(rh).*paramSet(rm).*paramSet(rl)];
    %maleActs(3:age , 1:risk) = maleActs(3:age , 1:risk) .* maleActsmult;
    maleActs = femaleActs;
    maleActs(3:5 , 1:risk) = maleActs(3:5 , 1:risk) .* paramSet(paramsSub{idx}.inds(1));
    maleActs(6:10 , 1:risk) = maleActs(6:10 , 1:risk) .* paramSet(paramsSub{idx}.inds(2));
else
    %load([paramDir , 'demoParamsFrmExcel'] , 'maleActs');
    maleActs(1:4 , 1:risk) = femaleActs(1:4 , 1:risk);
    for a = 5 : age
        maleActs(a , 1:risk) = femaleActs(a-1 , 1:risk);
    end
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
        'fertility' , 'fertility2' , 'fertility3' , 'fertility4' , ...
        'partnersM' , 'partnersF' , 'maleActs' , 'femaleActs'};
    vars5To1_vals = {riskDistM , riskDistF , mue , mue2 , mue3 , mue4 , ...
        fertility , fertility2 , fertility3 , fertility4 , ...
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
    epsA = [0.4921 ; 0.4921 ; 0.4921];
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
% file = [pwd , '/Config/HIV_parameters.xlsx'];
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
    baseVagTrans = [0.000575];
end

% HIV tranmission rate
vagTransM = (baseVagTrans*1.0) * ones(risk , 1); % probability of transmission from male (insertive) to female (receptive) based on male's disease state; female acquisition 
vagTransF = baseVagTrans * ones(risk , 1); % probability of transmission from female (receptive) to male (insertive) based on female's disease state; male acquisition 
betaHIV_F2M = bsxfun(@times , [9.0 1.0 2.5 7.0 0.7 0.0; 9.0 1.0 2.5 7.0 0.7 0.0; 9.0 1.0 2.5 7.0 0.7 0.0] , vagTransF);
betaHIV_M2F = bsxfun(@times , [9.0 1.0 2.5 7.0 0.7 0.0; 9.0 1.0 2.5 7.0 0.7 0.0; 9.0 1.0 2.5 7.0 0.7 0.0] , vagTransM);
betaHIV_F2M_red = bsxfun(@times , [9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0; 9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0; 9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0] , vagTransF);
betaHIV_M2F_red = bsxfun(@times , [9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0; 9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0; 9.0*0.5 1.0*0.5 2.5*0.5 7.0*0.5 0.7 0.0] , vagTransM);
betaHIV = zeros(gender , age , risk , viral);
betaHIV_red = zeros(gender , age , risk , viral);
for a = 1 : age % calculate per-partnership probability of HIV transmission
    % per-partnership beta: females to infect HIV-negative males, 
    % affected by betaHIV_F2M, probability of transmission from female (receptive) to male(insertive) based on female's disease state), and number of male acts
    betaHIV(2 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)'));
    betaHIV_red(2 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M_red , maleActs(a , :)'));
    % per-partnership beta: males to infect HIV-negative females,
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
    'kCin1_Cin2_orig' , 'kCin2_Cin3_orig' , 'ageTrends');

kCin1_Inf = zeros(16,2);
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
    kCin1_Inf(1 : 6 , 1) = kCin1_Inf_orig(1 , 1) * kCin1_InfMult(1);
    kCin1_Inf(1 : 6 , 2) = kCin1_Inf_orig(1 , 2) * kCin1_InfMult(2);
else
    kCin1_Inf(1 : 6 , 1) = kCin1_Inf_orig(1 , 1) * 0.8;
    kCin1_Inf(1 : 6 , 2) = kCin1_Inf_orig(1 , 2) * 2.0;
end

% CIN1 to CIN2, ages 10-24
if calibBool && any(28 == pIdx)
    idx = find(28 == pIdx);
    kCin2_Cin1Mult = paramSet(paramsSub{idx}.inds(:));
    kCin2_Cin1(1 : 6 , 1) = kCin2_Cin1_orig(1 , 1) * kCin2_Cin1Mult(1);
    kCin2_Cin1(1 : 6 , 2) = kCin2_Cin1_orig(1 , 2) * kCin2_Cin1Mult(2);
else
    kCin2_Cin1(1 : 6 , 1) = kCin2_Cin1_orig(1 , 1) * 1.4;
    kCin2_Cin1(1 : 6 , 2) = kCin2_Cin1_orig(1 , 2) * 2.0;
end

% CIN2 to CIN3, ages 10-24
if calibBool && any(29 == pIdx)
    idx = find(29 == pIdx);
    kCin3_Cin2Mult = paramSet(paramsSub{idx}.inds(:));
    kCin3_Cin2(1 : 6 , 1) = kCin3_Cin2_orig(1 , 1) * kCin3_Cin2Mult(1);
    kCin3_Cin2(1 : 6 , 2) = kCin3_Cin2_orig(1 , 2) * kCin3_Cin2Mult(2);
else
    kCin3_Cin2(1 : 6 , 1) = kCin3_Cin2_orig(1 , 1) * 1.3;
    kCin3_Cin2(1 : 6 , 2) = kCin3_Cin2_orig(1 , 2) * 2.0;
end

% CIN3 to unlocalized cancer, ages 10-24
if calibBool && any(30 == pIdx)
    idx = find(30 == pIdx);
    kCC_Cin3Mult = paramSet(paramsSub{idx}.inds(:));
    kCC_Cin3(1 : 6 , 1) = kCC_Cin3_orig(1 , 1) * kCC_Cin3Mult(1);
    kCC_Cin3(1 : 6 , 2) = kCC_Cin3_orig(1 , 2) * kCC_Cin3Mult(2);
else
    kCC_Cin3(1 : 6 , 1) = kCC_Cin3_orig(1 , 1) * 1.0;
    kCC_Cin3(1 : 6 , 2) = kCC_Cin3_orig(1 , 2) * 5.0;
end

% HPV to Well (or natural immunity), ages 10-24
if calibBool && any(31 == pIdx)
    idx = find(31 == pIdx);
    rNormal_InfMult = paramSet(paramsSub{idx}.inds(:));
    rNormal_Inf(1 : 6 , 1) = rNormal_Inf_orig(1 , 1) * rNormal_InfMult(1);
    rNormal_Inf(1 : 6 , 2) = rNormal_Inf_orig(1 , 2) * rNormal_InfMult(2);
else
    rNormal_Inf(1 : 6 , 1) = rNormal_Inf_orig(1 , 1) * 2.60;
    rNormal_Inf(1 : 6 , 2) = rNormal_Inf_orig(1 , 2) * 1.41;
end

% CIN1 to HPV, ages 10-24
if calibBool && any(32 == pIdx)
    idx = find(32 == pIdx);
    kInf_Cin1Mult = paramSet(paramsSub{idx}.inds(:));
    kInf_Cin1(1 : 6 , 1) = kInf_Cin1_orig(1 , 1) * kInf_Cin1Mult(1);
    kInf_Cin1(1 : 6 , 2) = kInf_Cin1_orig(1 , 2) * kInf_Cin1Mult(2);
else
    kInf_Cin1(1 : 6 , 1) = kInf_Cin1_orig(1 , 1) * 1.80;
    kInf_Cin1(1 : 6 , 2) = kInf_Cin1_orig(1 , 2) * 0.97;
end

% CIN2 to CIN1, ages 10-24
if calibBool && any(33 == pIdx)
    idx = find(33 == pIdx);
    kCin1_Cin2Mult = paramSet(paramsSub{idx}.inds(:));
    kCin1_Cin2(1 : 6 , 1) = kCin1_Cin2_orig(1 , 1) * kCin1_Cin2Mult(1);
    kCin1_Cin2(1 : 6 , 2) = kCin1_Cin2_orig(1 , 2) * kCin1_Cin2Mult(2);
else
    kCin1_Cin2(1 : 6 , 1) = kCin1_Cin2_orig(1 , 1) * 1.4;
    kCin1_Cin2(1 : 6 , 2) = kCin1_Cin2_orig(1 , 2) * 0.2;
end

% CIN3 to CIN2, ages 10-24
if calibBool && any(34 == pIdx)
    idx = find(34 == pIdx);
    kCin2_Cin3Mult = paramSet(paramsSub{idx}.inds(:));
    kCin2_Cin3(1 : 6 , 1) = kCin2_Cin3_orig(1 , 1) * kCin2_Cin3Mult(1);
    kCin2_Cin3(1 : 6 , 2) = kCin2_Cin3_orig(1 , 2) * kCin2_Cin3Mult(2);
else
    kCin2_Cin3(1 : 6 , 1) = kCin2_Cin3_orig(1 , 1) * 1.0;
    kCin2_Cin3(1 : 6 , 2) = kCin2_Cin3_orig(1 , 2) * 0.2;
end

% Apply age trends to 9v HPV transitions
kCin1_Inf(7 : 10 , 1) = kCin1_Inf(1 , 1) * (1.05 * linspace(1/4 , 1 , 4)); %ageTrends(1,1); % ages 25-49
kCin2_Cin1(7 : 10 , 1) = kCin2_Cin1(1 , 1) * (ageTrends(1,2) * linspace(1/4 , 1 , 4));
kCin3_Cin2(7 : 10 , 1) = kCin3_Cin2(1, 1) * (1.6 * linspace(1/4 , 1 , 4)); %ageTrends(1,3);
kCC_Cin3(7 : 10 , 1) = kCC_Cin3(1 , 1) * (2.0 * linspace(1/4 , 1 , 4)); %ageTrends(1,4);
rNormal_Inf(7 : 10 , 1) = rNormal_Inf(1 , 1) * linspace((1.0-1.0)*0.75+1.0 , 1.0 , 4); %ageTrends(1,5);
kInf_Cin1(7 : 10 , 1) = kInf_Cin1(1 , 1) * linspace((1.0-ageTrends(1,6))*0.75+ageTrends(1,6) , ageTrends(1,6) , 4);
kCin1_Cin2(7 : 10 , 1) = kCin1_Cin2(1 , 1) * linspace((1.0-ageTrends(1,7))*0.75+ageTrends(1,7) , ageTrends(1,7) , 4);
kCin2_Cin3(7 : 10 , 1) = kCin2_Cin3(1 , 1) * linspace((1.0-0.8)*0.75+0.8 , 0.8 , 4); %ageTrends(1,8);
kCin1_Inf(11 : 16 , 1) = kCin1_Inf(1 , 1) * linspace(1.05+((1.10-1.05)/6) , 1.10 , 6); %ageTrends(2,1); % ages 50-69
kCin2_Cin1(11 : 16 , 1) = kCin2_Cin1(1 , 1) * linspace(ageTrends(1,2)+((ageTrends(2,2)-ageTrends(1,2))/6) , ageTrends(2,2) , 6);
kCin3_Cin2(11 : 16 , 1) = kCin3_Cin2(1 , 1) * linspace(1.6+((2.0-1.6)/6) , 2.0 , 6); %ageTrends(2,3);
kCC_Cin3(11 : 16 , 1) = kCC_Cin3(1 , 1) * linspace(2.0+((3.0-2.0)/6) , 3.0 , 6); %ageTrends(2,4);
rNormal_Inf(11 : 16 , 1) = rNormal_Inf(1 , 1) * linspace((1.0-ageTrends(2,5))*(5/6)+ageTrends(2,5) , ageTrends(2,5) , 6);
kInf_Cin1(11 : 16 , 1) = kInf_Cin1(1 , 1) * linspace((ageTrends(1,6)-ageTrends(2,6))*(5/6)+ageTrends(2,6) , ageTrends(2,6) , 6);
kCin1_Cin2(11 : 16 , 1) = kCin1_Cin2(1 , 1) * linspace((ageTrends(1,7)-ageTrends(2,7))*(5/6)+ageTrends(2,7) , ageTrends(2,7) , 6);
kCin2_Cin3(11 : 16 , 1) = kCin2_Cin3(1 , 1) * linspace((0.8-0.6)*(5/6)+0.6 , 0.6 , 6); %ageTrends(2,8);
%kCin1_Inf(15 : 16 , 1) = kCin1_Inf(1 , 1) * 1.15; %ageTrends(3,1); % ages 70-79
%kCin2_Cin1(15 : 16 , 1) = kCin2_Cin1(1 , 1) * ageTrends(3,2);
%kCin3_Cin2(15 : 16 , 1) = kCin3_Cin2(1 , 1) * ageTrends(3,3);
%kCC_Cin3(15 : 16 , 1) = kCC_Cin3(1 , 1) * 4.0; %ageTrends(3,4);
%rNormal_Inf(15 : 16 , 1) = rNormal_Inf(1 , 1) * ageTrends(3,5);
%kInf_Cin1(15 : 16 , 1) = kInf_Cin1(1 , 1) * ageTrends(3,6);
%kCin1_Cin2(15 : 16 , 1) = kCin1_Cin2(1 , 1) * ageTrends(3,7);
%kCin2_Cin3(15 : 16 , 1) = kCin2_Cin3(1 , 1) * ageTrends(3,8);

% Apply age trends to non-9v HPV transitions
kCin1_Inf(7 : 10 , 2) = kCin1_Inf(1 , 2) * (1.05 * linspace(1/4 , 1 , 4)); %ageTrends(1,1); % ages 25-49
kCin2_Cin1(7 : 10 , 2) = kCin2_Cin1(1 , 2) * (ageTrends(1,2) * linspace(1/4 , 1 , 4));
kCin3_Cin2(7 : 10 , 2) = kCin3_Cin2(1, 2) * (1.6 * linspace(1/4 , 1 , 4)); %ageTrends(1,3);
kCC_Cin3(7 : 10 , 2) = kCC_Cin3(1 , 2) * (2.0 * linspace(1/4 , 1 , 4)); %ageTrends(1,4);
rNormal_Inf(7 : 10 , 2) = rNormal_Inf(1 , 2) * linspace((1.0-1.0)*0.75+1.0 , 1.0 , 4); %ageTrends(1,5);
kInf_Cin1(7 : 10 , 2) = kInf_Cin1(1 , 2) * linspace((1.0-ageTrends(1,6))*0.75+ageTrends(1,6) , ageTrends(1,6) , 4);
kCin1_Cin2(7 : 10 , 2) = kCin1_Cin2(1 , 2) * linspace((1.0-ageTrends(1,7))*0.75+ageTrends(1,7) , ageTrends(1,7) , 4);
kCin2_Cin3(7 : 10 , 2) = kCin2_Cin3(1 , 2) * linspace((1.0-0.8)*0.75+0.8 , 0.8 , 4); %ageTrends(1,8);
kCin1_Inf(11 : 16 , 2) = kCin1_Inf(1 , 2) * linspace(1.05+((1.10-1.05)/6) , 1.10 , 6); %ageTrends(2,1); % ages 50-69
kCin2_Cin1(11 : 16 , 2) = kCin2_Cin1(1 , 2) * linspace(ageTrends(1,2)+((ageTrends(2,2)-ageTrends(1,2))/6) , ageTrends(2,2) , 6);
kCin3_Cin2(11 : 16 , 2) = kCin3_Cin2(1 , 2) * linspace(1.6+((2.0-1.6)/6) , 2.0 , 6); %ageTrends(2,3);
kCC_Cin3(11 : 16 , 2) = kCC_Cin3(1 , 2) * linspace(2.0+((3.0-2.0)/6) , 3.0 , 6); %ageTrends(2,4);
rNormal_Inf(11 : 16 , 2) = rNormal_Inf(1 , 2) * linspace((1.0-ageTrends(2,5))*(5/6)+ageTrends(2,5) , ageTrends(2,5) , 6);
kInf_Cin1(11 : 16 , 2) = kInf_Cin1(1 , 2) * linspace((ageTrends(1,6)-ageTrends(2,6))*(5/6)+ageTrends(2,6) , ageTrends(2,6) , 6);
kCin1_Cin2(11 : 16 , 2) = kCin1_Cin2(1 , 2) * linspace((ageTrends(1,7)-ageTrends(2,7))*(5/6)+ageTrends(2,7) , ageTrends(2,7) , 6);
kCin2_Cin3(11 : 16 , 2) = kCin2_Cin3(1 , 2) * linspace((0.8-0.6)*(5/6)+0.6 , 0.6 , 6); %ageTrends(2,8);
%kCin1_Inf(15 : 16 , 2) = kCin1_Inf(1 , 2) * 1.15; %ageTrends(3,1); % ages 70-79
%kCin2_Cin1(15 : 16 , 2) = kCin2_Cin1(1 , 2) * ageTrends(3,2);
%kCin3_Cin2(15 : 16 , 2) = kCin3_Cin2(1 , 2) * ageTrends(3,3);
%kCC_Cin3(15 : 16 , 2) = kCC_Cin3(1 , 2) * 4.0; %ageTrends(3,4);
%rNormal_Inf(15 : 16 , 2) = rNormal_Inf(1 , 2) * ageTrends(3,5);
%kInf_Cin1(15 : 16 , 2) = kInf_Cin1(1 , 2) * ageTrends(3,6);
%kCin1_Cin2(15 : 16 , 2) = kCin1_Cin2(1 , 2) * ageTrends(3,7);
%kCin2_Cin3(15 : 16 , 2) = kCin2_Cin3(1 , 2) * ageTrends(3,8);

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
    lambdaMultImmmult = paramSet(paramsSub{idx}.inds(:));
    lambdaMultImm = ones(age , 1) .* lambdaMultImmmult;   
else
    lambdaMultImm = ones(age , 1) .* (0.57479*0.90);
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

% HPV immunity clearance multiplier for HIV-positive persons 
rImmuneHiv = [1.4167; 1.5682; 1.9722; 2.8333];

% HPV infection multiplier for HIV-positive persons
hpv_hivMult = [1.78; 1.99; 2.12; 2.32];

% HPV clearance multiplier for HIV-positive persons
if calibBool && any(14 == pIdx)
    idx = find(14 == pIdx);
    hpv_hivClear = zeros(4 , 1);
    hpv_hivClear(1,1) = paramSet(paramsSub{idx}.inds(1));
    hpv_hivClear(2,1) = hpv_hivClear(1,1)*paramSet(paramsSub{idx}.inds(2));
    hpv_hivClear(3,1) = hpv_hivClear(2,1)*paramSet(paramsSub{idx}.inds(3));
    hpv_hivClear(4,1) = hpv_hivClear(3,1)*paramSet(paramsSub{idx}.inds(4));
else
    hpv_hivClear = [0.60; 0.55; 0.45; 0.30];
end

% CIN2 to CIN3 progression multiplier for HIV-positive women
if calibBool && any(15 == pIdx)
    idx = find(15 == pIdx);
    c3c2multiplier = paramSet(paramsSub{idx}.inds(:));
    c3c2Mults = [1.0; 2.0; 2.7; 3.5];
    c3c2Mults(2 : end) = c3c2Mults(2 : end).*c3c2multiplier;
else
    c3c2Mults = [1.0; 2.0; 2.7; 3.5];
end

% CIN1 to CIN2 progression multiplier for HIV-positive women
if calibBool && any(16 == pIdx)
    idx = find(16 == pIdx);
    c2c1multiplier = paramSet(paramsSub{idx}.inds(:));
    c2c1Mults = [1.0; 1.9; 2.4; 2.7];
    c2c1Mults(2 : end) = c2c1Mults(2 : end).*c2c1multiplier;
else
    c2c1Mults = [1.0; 1.9; 2.4; 2.7];
end

% CIN3 to CIN2 regression multiplier for HIV-positive women
if calibBool && any(38 == pIdx)
    idx = find(38 == pIdx);
    c2c3multiplier = paramSet(paramsSub{idx}.inds(:));
    c2c3Mults = [0.60; 0.55; 0.45; 0.30].*c2c3multiplier;
else
    c2c3Mults = [0.60; 0.55; 0.45; 0.30];
end

% CIN2 to CIN1 regression multiplier for HIV-positive women
if calibBool && any(39 == pIdx)
    idx = find(39 == pIdx);
    c1c2multiplier = paramSet(paramsSub{idx}.inds(:));
    c1c2Mults = [0.60; 0.55; 0.45; 0.30].*c1c2multiplier;
else
    c1c2Mults = [0.60; 0.55; 0.45; 0.30];
end

% HPV tranmission rates
if calibBool && any(10 == pIdx)
    idx = find(10 == pIdx);
    perPartnerHpv_vax = paramSet(paramsSub{idx}.inds(:));
else
    perPartnerHpv_vax = 0.01;
end

if calibBool && any(11 == pIdx)
    idx = find(11 == pIdx);
    perPartnerHpv_nonV = paramSet(paramsSub{idx}.inds(:));
else
    perPartnerHpv_nonV = perPartnerHpv_vax;
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
    % per-partnership beta: females to infect HPV-negative males 
    % affected by betaHPV_F2M_vax, probability of transmission based on cervical cancer status/progression, and number of male acts
    beta_hpvVax(2 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHPV_vax , maleActs(a , :)'));
    % per-partnership beta: males to infect HPV-negative females,
    % affected by betaHPV_M2F_vax, probability of transmission based on cervical cancer status/progression, and number of female acts
    beta_hpvVax(1 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHPV_vax , femaleActs(a , :)'));
    
    % NON-VACCINE-TYPE HPV
    % per-partnership beta: females to infect HPV-negative males 
    % affected by betaHPV_F2M_nonV, probability of transmission based on cervical cancer status/progression, and number of male acts
    beta_hpvNonVax(2 , a , : , :) = 1 - (bsxfun(@power, 1 - betaHPV_nonV , maleActs(a , :)'));
    % per-partnership beta: males to infect HPV-negative females,
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
    maleHpvClearMult = 3.5;
end

% HIV and ART multipliers
if calibBool && any(21 == pIdx)
    idx = find(21 == pIdx);
    artHpvMult = paramSet(paramsSub{idx}.inds(:));
else
    artHpvMult = 1.0;
end

%% Save intervention parameters

% Import from Excel HIV intervention parameters
% file = [pwd , '/Config/HIV_parameters.xlsx'];
% circProtect = xlsread(file , 'Protection' , 'B18');
% condProtect = xlsread(file , 'Protection' , 'B19');
% MTCTRate = xlsread(file , 'Disease Data' , 'B6:B8');
% artVScov = xlsread(file , 'Protection' , 'A33:C41');    % [years , females , males] 
% save(fullfile(paramDir ,'hivIntParamsFrmExcel'), 'circProtect' , ...
%     'condProtect' , 'MTCTRate' , 'artVScov');

% Load pre-saved HIV intervention parameters
load([paramDir , 'hivIntParamsFrmExcel'] , 'circProtect' , ...
    'condProtect' , 'MTCTRate' , 'artVScov');

% Intervention start years
hivStartYear = 1980;
circStartYear = 1960; % Year VMMC begins in ages 15-19
circNatStartYear = 2010; % Year SA National VMMC program begins, older men circumcised
vaxStartYear = 2014;

% Protection from circumcision and condoms
circProtect = [[circProtect; 0.0] , [0.30; 0.0]];    % HIV protection , HPV protection
condProtect = [ones(gender,1).*condProtect , [0.46; 0.70]];    % HIV protection , HPV protection

% Condom use
if calibBool && any(5 == pIdx);
    idx = find(5 == pIdx);
    condUse = paramSet(paramsSub{idx}.inds(:));
else
    if fivYrAgeGrpsOn
        condUse = 0.5 * 0.5;
    else
        condUse = 0.20;
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
maxRateM = [artVScov(:,3) ; 0.729] .* artOutMult; % population-level ART coverage in males
maxRateF = [artVScov(:,2) ; 0.729] .* artOutMult; % population-level ART coverage in females
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

% VMMC coverage
vmmcYr = [circStartYear; 2000; 2008; 2010; 2012; 2017; 2030];
circ_aVec = {4 , 5 , [6:10] , [11:age]}; % Ages: (15-19), (20-24), (25-49), (50+)
vmmcRate = [0.04 0.06 0.0 0.0; ... % 1960
            0.10 0.13 0.0 0.0; ... % 2000
            0.114 0.161 0.0 0.0; ... % 2008
            0.143 0.201 0.14 0.12; ... % 2010
            0.172 0.242 0.191 0.143; ... % 2012
            0.459 0.42 0.318  0.204; ... % 2017 
            0.70  0.70 0.70   0.70];   % 2030 [year x age group]
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
baseline.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.48; 0.48];
baseline.diseaseInds = [1 : disease];
baseline.screenAge = [35/max(1 , fivYrAgeGrpsOn*5)+1];
baseline.screenAgeMults = [1.0]; % / max(1 , fivYrAgeGrpsOn*5)
baseline.testSens = cytoSens;
% cryoElig = [1.0 , 0.85 , 0.75 , 0.10 , 0.10 , 0.10];
baseline.colpoRetain = 0.72;
baseline.cinTreatEff = [0.905 , 0.905 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766]; % cryotherapy/LEEP effectiveness by HIV status
baseline.cinTreatRetain = 0.51;
baseline.cinTreatHpvPersist = 0.28; % HPV persistence with LEEP including treatment failure
baseline.cinTreatHpvPersistHivNeg = baseline.cinTreatHpvPersist - (1-baseline.cinTreatEff(1)); % 0.185; proportion of effectively treated HIV-negative women who have persistent HPV after LEEP
baseline.ccTreatRetain = 0.40;
baseline.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    baseline.screenCover_vec{i} = interp1(period , baseline.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% CISNET screening algorithm
cisnet.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.70; 0.90];
cisnet.screenAge = [(35/max(1 , fivYrAgeGrpsOn*5)+1) , (45/max(1 , fivYrAgeGrpsOn*5)+1)];
cisnet.screenAgeMults = [(1.0 / max(1 , fivYrAgeGrpsOn*5)) , (1.0 / max(1 , fivYrAgeGrpsOn*5))];
cisnet.testSens = hpvSens;
cisnet.colpoRetain = 0.81*0.85; % (compliance) * (CIN2+/CC correctly identified by same-day colposcopy)
cisnet.cinTreatEff = baseline.cinTreatEff;
cisnet.cinTreatRetain = 1.0;
cisnet.cinTreatHpvPersist = 0.48; % HPV persistence with cryotherapy including treatment failure
cisnet.cinTreatHpvPersistHivNeg = cisnet.cinTreatHpvPersist - (1-cisnet.cinTreatEff(1)); % proportion of effectively treated HIV-negative women who have persistent HPV after cryotherapy
cisnet.ccTreatRetain = 1.0;
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
% file = [pwd , '/Config/Calibration_targets.xlsx'];
% 
% ccInc2012_dObs(: , 1) = xlsread(file , 'Calibration' , 'D275 : D287'); % CC Incidence Rate 2012, by age
% ccInc2012_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H275 : I287');
% 
% ccInc2018_dObs(: , 1) = xlsread(file , 'Calibration' , 'D294 : D305'); % CC Incidence Rate 2018, by age
% ccInc2018_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H294 : I305');
% 
% cc_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'D183 : D184'); % CC type distribution for 9v and non-9v HPV types
% cc_dist_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H183 : I184');
% 
% cin3_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'D189 : D190'); % CIN3 type distribution for 9v and non-9v HPV types
% cin3_dist_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H189 : I190');
% 
% cin1_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'D187 : D188'); % CIN1 type distribution for 9v and non-9v HPV types
% cin1_dist_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H187 : I188');
% 
% hpv_dist_dObs(: , 1) = xlsread(file , 'Calibration' , 'D185 : D186'); % HPV type distribution for 9v and non-9v HPV types
% hpv_dist_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H185 : I186');
% 
% cinPos2002_dObs(: , 1) = xlsread(file , 'Calibration' , 'D2 : D11'); % CIN2/CIN3 Prevalence (HIV+) 2002, by age
% cinPos2002_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H2 : I11');
% cinNeg2002_dObs(: , 1) = xlsread(file , 'Calibration' , 'D12 : D21'); % CIN2/CIN3 Prevalence (HIV-) 2002, by age
% cinNeg2002_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H12 : I21');
% 
% cinPos2015_dObs(: , 1) = xlsread(file , 'Calibration' , 'D288 : D290'); % CIN1, CIN2, CIN3 Prevalence (HIV+) 2015, ages 30-65
% cinPos2015_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H288 : I290');
% cinNeg2015_dObs(: , 1) = xlsread(file , 'Calibration' , 'D291 : D293'); % CIN1, CIN2, CIN3 Prevalence (HIV-) 2015, ages 30-65
% cinNeg2015_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H291 : I293');
% 
% hpv_hiv_dObs(: , 1) = xlsread(file , 'Calibration' , 'D144 : D152'); % HPV Prevalence in HIV+ Women (All) 2014, by age
% hpv_hiv_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H144 : I152');
% hpv_hivNeg_dObs(: , 1) = xlsread(file , 'Calibration' , 'D153 : D161'); % HPV Prevalence in HIV- Women (All) 2014, by age
% hpv_hivNeg_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H153 : I161');
% 
% hpv_hivM2008_dObs(: , 1) = xlsread(file , 'Calibration' , 'D52 : D55'); % HPV Prevalence in HIV+ Men, by age
% hpv_hivM2008_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H52 : I55');
% hpv_hivMNeg2008_dObs(: , 1) = xlsread(file , 'Calibration' , 'D56 : D59'); % HPV Prevalence in HIV- Men, by age
% hpv_hivMNeg2008_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H56 : I59');
% 
% hivPrevM_dObs(: , 1) = xlsread(file , 'Calibration' , 'D60 : D101'); % HIV Prevalence in Men 2003,2005,2006,2007,2008,2009, by age
% hivPrevM_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H60 : I101');
% hivPrevM_dObs(: , 4 : 5) = xlsread(file , 'Calibration' , 'E60 : F101'); % raw data
% hivPrevF_dObs(: , 1) = xlsread(file , 'Calibration' , 'D102 : D143'); % HIV Prevalence in Women 2003,2005,2006,2007,2008,2009, by age
% hivPrevF_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H102 : I143');
% hivPrevF_dObs(: , 4 : 5) = xlsread(file , 'Calibration' , 'E102 : F143'); % raw data
% 
% popAgeDist_dObs(: , 1) = [xlsread(file , 'Calibration' , 'D191 : D206'); ... % Population age distribution in 1996, 2011, and 2019
%     xlsread(file , 'Calibration' , 'D223 : D238');
%     xlsread(file , 'Calibration' , 'D255 : D270')];
% popAgeDist_dObs(: , 2 : 3) = [xlsread(file , 'Calibration' , 'H191 : I206'); ...
%     xlsread(file , 'Calibration' , 'H223 : I238');
%     xlsread(file , 'Calibration' , 'H255 : I270')];
% 
% totPopSize_dObs(: , 1) = xlsread(file , 'Calibration' , 'D271 : D273'); % Total population size in 2001, 2011, and 2019
% totPopSize_dObs(: , 2 : 3) = xlsread(file , 'Calibration' , 'H271 : I273');
% 
% save(fullfile(paramDir , 'calibData'), 'ccInc2012_dObs' , 'ccInc2018_dObs' , 'cc_dist_dObs' , 'cin3_dist_dObs' , ...
%     'cin1_dist_dObs' , 'hpv_dist_dObs' , 'cinPos2002_dObs' , 'cinNeg2002_dObs' , ...
%     'cinPos2015_dObs' , 'cinNeg2015_dObs' , 'hpv_hiv_dObs' , 'hpv_hivNeg_dObs' , ...
%     'hpv_hivM2008_dObs' , 'hpv_hivMNeg2008_dObs' , 'hivPrevM_dObs' , 'hivPrevF_dObs' , ...
%     'popAgeDist_dObs' , 'totPopSize_dObs')

load([paramDir , 'calibData'], 'ccInc2012_dObs' , 'ccInc2018_dObs' , 'cc_dist_dObs' , 'cin3_dist_dObs' , ...
    'cin1_dist_dObs' , 'hpv_dist_dObs' , 'cinPos2002_dObs' , 'cinNeg2002_dObs' , ...
    'cinPos2015_dObs' , 'cinNeg2015_dObs' , 'hpv_hiv_dObs' , 'hpv_hivNeg_dObs' , ...
    'hpv_hivM2008_dObs' , 'hpv_hivMNeg2008_dObs' , 'hivPrevM_dObs' , 'hivPrevF_dObs' , ...
    'popAgeDist_dObs' , 'totPopSize_dObs');




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

%% Fertility by 2035

% birth indices
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));
posMaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(3 , 1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix for uninfected mothers
% disp('Building fertility4 matrix for uninfected mothers')
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
for d = 3 : 7 % hiv infected
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
dFertPos1 = (hivFertPosBirth2 - hivFertPosBirth) ./ ((2000 - 1960) * stepsPerYear);
dFertNeg1 = (hivFertNegBirth2 - hivFertNegBirth) ./ ((2000 - 1960) * stepsPerYear); 
dFertMat1 = (fertMat2 - fertMat) ./ ((2000 - 1960) * stepsPerYear);

dFertPos2 = (hivFertPosBirth3 - hivFertPosBirth2) ./ ((2020 - 2010) * stepsPerYear);
dFertNeg2 = (hivFertNegBirth3 - hivFertNegBirth2) ./ ((2020 - 2010) * stepsPerYear);
dFertMat2 = (fertMat3 - fertMat2) ./ ((2020 - 2010) * stepsPerYear);

dFertPos3 = (hivFertPosBirth4 - hivFertPosBirth3) ./ ((2035 - 2020) * stepsPerYear);
dFertNeg3 = (hivFertNegBirth4 - hivFertNegBirth3) ./ ((2035 - 2020) * stepsPerYear);
dFertMat3 = (fertMat4 - fertMat3) ./ ((2035 - 2020) * stepsPerYear);

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
