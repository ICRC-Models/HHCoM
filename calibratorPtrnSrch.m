function [negSumLogL] = calibratorPtrnSrch(paramSet)

%% Load parameters
paramDir = [pwd ,'/Params/'];
% load([paramDir,'settings'])
% load([paramDir,'popData'])
% load([paramDir,'HIVParams'])
% load([paramDir,'mixInfectParams'])
% load([paramDir,'vlBeta'])
% load([paramDir,'hpvData'])
% load([paramDir,'cost_weights'])
load([paramDir , 'calibratedParams'])
perPartnerHpv = 0.0045;
condUse = 0.5 * 0.5;
epsA = [0.3 ; 0.3 ; 0.3];
epsR = [0.3 ; 0.3 ; 0.3];
muHIV(11,2) = 0.02;
OMEGA = zeros(age , 1); % hysterectomy rate
paramDir = [pwd ,'/Params/'];
% Load parameters
load([paramDir,'general'])
load([paramDir,'calibData'])
% Load indices
load([paramDir,'mixInfectIndices'])
load([paramDir,'hivIndices'])
load([paramDir,'hpvIndices'])
load([paramDir,'ageRiskInds'])
% Load matrices
load([paramDir,'ager'])
load([paramDir,'vlAdvancer'])
load([paramDir,'fertMat'])
load([paramDir,'hivFertMats'])
load([paramDir,'fertMat2'])
load([paramDir,'hivFertMats2'])
load([paramDir,'deathMat'])
load([paramDir,'circMat'])
load([paramDir,'circMat2'])

% Load calibration parameter information
pIdx = load([paramDir,'pIdx_patternSrch_' , '22Aug19' , '_0.dat']);
[paramsAll] = genParamStruct();
paramsSub = cell(length(pIdx),1);
startIdx = 1;
for s = 1 : length(pIdx)
    paramsSub{s}.length = paramsAll{pIdx(s)}.length;
    paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1));
    startIdx = startIdx + paramsSub{s}.length;
end

% Model specs
hyst = 0;
hpvOn = 1;
hivOn = 1;

% Use newly calibrated parameters
% masterSetMatrix = load([paramDir , 'masterSets_calib_05Aug19_0.dat']); % load most recent parameter sample
% paramSet = masterSetMatrix(:,1424);
% pIdx = [1,2,5,6,7,8,9,10,15,16,17,19,22];
% % paramSet = [1.1884; 3.5575; 0.28806; 0.24017; 0.53502; 0.10428; ...
% %    0.21381; 0.43467; 12.663; 1.3583; 0.48721];    % 2179
% % paramSet = [4.7022; 0.47332; 0.81494; 0.70898; 0.76744; 0.38557; 0.26816; ...
% % 0.75175; 0.37118; 8.7002; 13.258; 0.35454; 0.57692]; % 1.7591e+05; 857
% % paramSet = [1.4894; 1.9707; 0.84234; 0.9383; 0.82112; 0.56634; 0.77342; ...
% %     0.53429; 0.25521; 3.4473; 6.7432; 0.7543; 0.67004];
% [paramsAll] = genParamStruct();
% paramsSub = cell(length(pIdx),1);
% startIdx = 1;
% for s = 1 : length(pIdx)
%     paramsSub{s}.length = paramsAll{pIdx(s)}.length;
%     paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1));
%     startIdx = startIdx + paramsSub{s}.length;
% end

% Time
c = fix(clock);
currYear = c(1); % get the current year
startYear = 1910; %1980;
endYear = 2015; % run to 2015 for calibration
years = endYear - startYear;
timeStep = 1 / stepsPerYear;

% Intervention start years
hivStartYear = 1980;
circStartYear = 1990;
vaxStartYear = 2014;

% Helper functions
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

%% Calibration parameters (FEED FROM LHS)
if any(1 == pIdx)
    idx = find(1 == pIdx);
    partnersMmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    partnersM(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
    %partnersM(3:10, 1:risk) = [paramSet(rl).*paramSet(rm).*paramSet(rh) , paramSet(rm).*paramSet(rh) , paramSet(rh)];
    %partnersM(11:age , 1:risk) = ones(6,risk).*partnersM(10 , 1:risk);
    partnersM(3:age , 1:risk) = partnersM(3:age , 1:risk) .* partnersMmult;
    % partnersM(10:age , 3) = ones(7 , 1);
end
if any(2 == pIdx)
    idx = find(2 == pIdx);
    partnersFmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    partnersF(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
    %partnersF(3:10 , 1:risk) = [paramSet(rl).*paramSet(rm).*paramSet(rh) , paramSet(rm).*paramSet(rh) , paramSet(rh)];
    %partnersF(11:age , 1:risk) = ones(6,risk).*partnersF(10 , 1:risk);
    partnersF(3:age , 1:risk) = partnersF(3:age , 1:risk) .* partnersFmult;
    % partnersF(10:age , 3) = ones(7 , 1);
end
% partnersM(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
% partnersF(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
% partnersM(3:age , 1:risk) = [paramSet(1:14) , paramSet(15:28) , paramSet(29:42)];
% partnersF(3:age , 1:risk) = [paramSet(43:56) , paramSet(57:70) , paramSet(71:84)];
% partnersM(10:age , 3) = ones(7 , 1);
% partnersF(10:age , 3) = ones(7 , 1);
if any(5 == pIdx);
    idx = find(5 == pIdx);
    condUse = paramSet(paramsSub{idx}.inds(:));
end
% condUse = paramSet(169);
if any(6 == pIdx);
    idx = find(6 == pIdx);
    %epsA = paramSet(paramsSub{idx}.inds(:));
    epsA = ones(3,1).*paramSet(paramsSub{idx}.inds(:));
end
if any(7 == pIdx);
    idx = find(7 == pIdx);
    %epsR = paramSet(paramsSub{idx}.inds(:));
    epsR = ones(3,1).*paramSet(paramsSub{idx}.inds(:));
end
% epsA = paramSet(170:172);
% epsR = paramSet(173:175);
if any(8 == pIdx)
    idx = find(8 == pIdx);
    maleActsmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    maleActs(1:2 , 1:risk) = zeros(2 , risk);
    %maleActs(3:age , 1:risk) = [paramSet(rl) , paramSet(rm).*paramSet(rl) , paramSet(rh).*paramSet(rm).*paramSet(rl)];
    maleActs(3:age , 1:risk) = maleActs(3:age , 1:risk) .* maleActsmult;
end
if any(9 == pIdx)
    idx = find(9 == pIdx);
    femaleActsmult = paramSet(paramsSub{idx}.inds(:));
    %rowL = paramsSub{idx}.length/3;
    %rl = paramsSub{idx}.inds(1:rowL);
    %rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    %rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    femaleActs(1:2 , 1:risk) = zeros(2 , risk);
    %femaleActs(3:age , 1: risk) = [paramSet(rl) , paramSet(rm).*paramSet(rl) , paramSet(rh).*paramSet(rm).*paramSet(rl)];
    femaleActs(3:age , 1: risk) = femaleActs(3:age , 1: risk) .* femaleActsmult;
end
% maleActs(1:2 , 1:risk) = zeros(2 , risk);
% femaleActs(1:2 , 1:risk) = zeros(2 , risk);
% maleActs(3:age , 1:risk) = [paramSet(176:189) , paramSet(190:203) , paramSet(204:217)];
% femaleActs(3:age , 1:risk) = [paramSet(218:231) , paramSet(232:245) , paramSet(246:259)];
if any(10 == pIdx)
    idx = find(10 == pIdx);
    perPartnerHpv = paramSet(paramsSub{idx}.inds(:)).*0.1;
end
% perPartnerHpv = paramSet(260);
% perPartnerHpv_lr = paramSet(261);
% perPartnerHpv_nonV = paramSet(262);
% hpv_hivMult = paramSet(263:266) .* 2.0;
% rNormal_Inf_orig = rNormal_Inf;
% rNormal_Inf = paramSet(267) .* rNormal_Inf_orig;
if any(15 == pIdx)
    idx = find(15 == pIdx);
    hpv_hivClear(1,1) = paramSet(paramsSub{idx}.inds(1));
    hpv_hivClear(2,1) = hpv_hivClear(1,1)*paramSet(paramsSub{idx}.inds(2));
    hpv_hivClear(3,1) = hpv_hivClear(2,1)*paramSet(paramsSub{idx}.inds(3));
    hpv_hivClear(4,1) = hpv_hivClear(3,1)*paramSet(paramsSub{idx}.inds(4));
end
% hpv_hivClear = paramSet(283:286) .* 0.5;
if any(16 == pIdx)
    idx = find(16 == pIdx);
    c3c2Mults(4,1) = paramSet(paramsSub{idx}.inds(3));
    c3c2Mults(3,1) = c3c2Mults(4,1)*paramSet(paramsSub{idx}.inds(2));
    c3c2Mults(2,1) = c3c2Mults(3,1)*paramSet(paramsSub{idx}.inds(1));
end
% c3c2Mults = paramSet(287:290) .* 2.0;
if any(17 == pIdx)
    idx = find(17 == pIdx);
    c2c1Mults(4,1) = paramSet(paramsSub{idx}.inds(3));
    c2c1Mults(3,1) = c3c2Mults(4,1)*paramSet(paramsSub{idx}.inds(2));
    c2c1Mults(2,1) = c3c2Mults(3,1)*paramSet(paramsSub{idx}.inds(1));
end
% c2c1Mults = paramSet(291:294) .* 2.0;
% kCCDet = paramSet(295:297);
if any(19 == pIdx)
    idx = find(19 == pIdx);
    lambdaMultImmmult = paramSet(paramsSub{idx}.inds(:));
    lambdaMultImm = lambdaMultImm .* lambdaMultImmmult;
end
% lambdaMultImm = paramSet(298:313);
% maxRateM_vec = paramSet(314:315);
% maxRateF_vec = paramSet(316:317);
if any(22 == pIdx)
    idx = find(22 == pIdx);
    artHpvMult = paramSet(paramsSub{idx}.inds(:));
end
% artHpvMult = paramSet(318) .* 2.0;
if any(25 == pIdx)
    idx = find(25 == pIdx);
    kCCcin3mult = paramSet(paramsSub{idx}.inds(:));
    kCC_Cin3 = kCC_Cin3 .* kCCcin3mult;
end
% kCC_Cin3
% kCD4(1,:,:) = [paramSet(319:323) , paramSet(324:328) , paramSet(329:333) , paramSet(334:338)];
% kCD4(2,:,:) = [paramSet(339:343) , paramSet(344:348) , paramSet(349:353) , paramSet(354:358)];
% kVl(1,:,:) = [paramSet(359:363) , paramSet(364:368) , paramSet(369:373) , paramSet(374:378)];
% kVl(2,:,:) = [paramSet(379:383) , paramSet(384:388) , paramSet(389:393) , paramSet(394:398)];

maxRateM_vec = [0.4 , 0.4];
maxRateF_vec = [0.55 , 0.55];
maxRateM1 = maxRateM_vec(1);
maxRateM2 = maxRateM_vec(2);
maxRateF1 = maxRateF_vec(1);
maxRateF2 = maxRateF_vec(2);

% % rImmuneHiv = hpv_hivClear;
% % lambdaMultImm = ones(age,1);
% % artHpvMult = 1.0;
% %
% % for a = 1 : age
% %     betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
% %     betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
% % end
% % betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
% % betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl
% 
epsA_vec = cell(size(yr , 1) - 1, 1); % save data over time interval in a cell array
epsR_vec = cell(size(yr , 1) - 1, 1);
for i = 1 : size(yr , 1) - 1          % interpolate epsA/epsR values at steps within period
    period = [yr(i) , yr(i + 1)];
    epsA_vec{i} = interp1(period , epsA(i : i + 1 , 1) , ...
        yr(i) : timeStep : yr(i + 1));
    epsR_vec{i} = interp1(period , epsR(i : i + 1 , 1) , ...
        yr(i) : timeStep : yr(i + 1));
end
% epsA_vec = epsA;
% epsR_vec = epsR;

%% Screening
screenYrs = [2000; 2003; 2016; currYear; 2023; 2030; 2045];
hpvScreenStartYear = screenYrs(1);
cytoSens = [0.0 , 0.0 , 0.57 , 0.57 , 0.57 , 0.57 , 0.57 , 0.0 , 0.0 , 0.0];

% Baseline screening algorithm
baseline.screenCover = [0.0; 0.18; 0.48; 0.48; 0.48; 0.48; 0.48];
baseline.screenAge = 8;
baseline.testSens = cytoSens;
% cryoElig = [1.0 , 0.85 , 0.75 , 0.10 , 0.10 , 0.10];
baseline.colpoRetain = 0.72;
baseline.cinTreatEff = [0.905 , 0.766 , 0.766 , 0.766 , 0.766 , 0.766 , 0.905 , 0.905 , 0.905 , 0.766]; % cryotherapy/LEEP effectiveness by HIV status
baseline.cinTreatRetain = 0.51;
baseline.cinTreatHpvPersist = 0.28; % HPV persistence with LEEP
baseline.ccTreatRetain = 0.40;

screenAlgorithm = 1;
screenAlgs{1} = baseline;
screenAlgs{1}.diseaseInds = [1 : disease];

screenAlgs{1}.screenCover_vec = cell(size(screenYrs , 1) - 1, 1); % save data over time interval in a cell array
for i = 1 : size(screenYrs , 1) - 1          % interpolate dnaTestCover values at steps within period
    period = [screenYrs(i) , screenYrs(i + 1)];
    screenAlgs{1}.screenCover_vec{i} = interp1(period , screenAlgs{1}.screenCover(i : i + 1 , 1) , ...
        screenYrs(i) : timeStep : screenYrs(i + 1));
end

% Create screening indices
screenAgeAll = zeros(disease , viral , hpvTypes , hpvStates , periods , length(screenAlgs{1}.screenAge) , risk);
screenAgeS = zeros(disease , viral , hpvTypes , hpvStates , 2 , length(screenAlgs{1}.screenAge) , risk);
noVaxNoScreen = zeros(disease , viral , hpvTypes , hpvStates , length(screenAlgs{1}.screenAge) , risk);
noVaxToScreen = noVaxNoScreen;
vaxNoScreen = noVaxNoScreen;
vaxToScreen = noVaxNoScreen;
noVaxToScreenTreatImm = zeros(disease , viral , length(screenAlgs{1}.screenAge) , risk);
vaxToScreenTreatImm = noVaxToScreenTreatImm;
noVaxToScreenTreatHpv = noVaxToScreenTreatImm;
vaxToScreenTreatHpv = noVaxToScreenTreatImm;
noVaxToScreenHyst = noVaxToScreenTreatImm;
vaxToScreenHyst = noVaxToScreenTreatImm;
noVaxScreen = zeros(disease*viral*hpvTypes*hpvStates*risk , length(screenAlgs{1}.screenAge));
noVaxXscreen = noVaxScreen;
vaxScreen = noVaxScreen;
vaxXscreen = noVaxScreen;
for aS = 1 : length(screenAlgs{1}.screenAge)
    a = screenAlgs{1}.screenAge(aS);
    
    for d = 1 : disease
        for v = 1 : viral
            for h = 1 : hpvTypes
                for s = 1 : hpvStates
                    for r = 1 : risk
                        screenAgeAll(d,v,h,s,:,aS,r) = toInd(allcomb(d , v , h , s , 1 : periods , 2 , a , r)); 
                        screenAgeS(d,v,h,s,:,aS,r) = toInd(allcomb(d , v , h , s , [4,6] , 2 , a , r));

                        noVaxNoScreen(d,v,h,s,aS,r) = sort(toInd(allcomb(d , v , h , s , 1 , 2 , a , r)));
                        noVaxToScreen(d,v,h,s,aS,r) = sort(toInd(allcomb(d , v , h , s , 6 , 2 , a , r)));
                        vaxNoScreen(d,v,h,s,aS,r) = sort(toInd(allcomb(d , v , h , s , 2 , 2 , a , r)));
                        vaxToScreen(d,v,h,s,aS,r) = sort(toInd(allcomb(d , v , h , s , 4 , 2 , a , r)));

                        noVaxToScreenTreatImm(d,v,aS,r) = toInd(allcomb(d , v , 1 , 10 , 6 , 2 , a , r));
                        vaxToScreenTreatImm(d,v,aS,r) = toInd(allcomb(d , v , 1 , 10 , 4 , 2 , a , r));
                        noVaxToScreenTreatHpv(d,v,aS,r) = toInd(allcomb(d , v , 2 , 1 , 6 , 2 , a , r));
                        vaxToScreenTreatHpv(d,v,aS,r) = toInd(allcomb(d , v , 2 , 1 , 4 , 2 , a , r));
                        noVaxToScreenHyst(d,v,aS,r) = toInd(allcomb(d , v , 1 , 8 , 6 , 2 , a , r));
                        vaxToScreenHyst(d,v,aS,r) = toInd(allcomb(d , v , 1 , 8 , 4 , 2 , a , r));
                    end
                end
            end
        end

    end

    % Create indices for removing screening status as people age out of screened age groups
    noVaxScreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 6 , ... 
        2 , a+1 , 1 : risk));
    noVaxXscreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 , ... 
        2 , a+1 , 1 : risk));
    vaxScreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 4 , ... 
        2 , a+1 , 1 : risk));
    vaxXscreen(:,aS) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 2 , ... 
        2 , a+1 , 1 : risk));
end

%% Vaccination
% actually bivalent vaccine, but to avoid adding additional compartments,
% we use nonavalent vaccine and then reduce coverage
vaxEff = 0.9;    

%Parameters for school-based vaccination regimen
vaxAge = 2;
vaxRate = 0.86*(0.7/0.9);    % (9 year-olds = 1/5th of age group) * (bivalent vaccine efficacy adjustment)
vaxG = 2;   % indices of genders to vaccinate (1 or 2 or 1,2)

% Parameters for waning
waning = 0;    % turn waning on or off

lambdaMultVaxMat = zeros(age , 1);   % age-based vector for modifying lambda based on vaccination status

% No waning
lambdaMultVaxMat(3 : age) = vaxEff;

% Waning
effPeriod = 20; % number of years that initial efficacy level is retained
wanePeriod = 20; % number of years over which initial efficacy level wanes
if waning 
    % Following a period (in years) where original efficacy is retained, 
    % specified by 'effPeriod' , linearly scale down vaccine efficacy 
    % to 0% over time period specificed by 'wanePeriod'
    % To make waning rate equal in all scenarios, the linear rate of 
    % waning is based on the least effective initial vaccine efficacy scenario.        
    kWane = vaxEff / round(wanePeriod / 5);     
    vaxInit = vaxEff;
    lambdaMultVaxMat(round(effPeriod / 5) + vaxAge - 1 : age) = ...
        max(0 , linspace(vaxInit , ...
        vaxInit - kWane * (1 + age - (round(wanePeriod / 5) + vaxAge)) ,...
        age - (round(wanePeriod / 5) + vaxAge) + 2)'); % ensures vaccine efficacy is >= 0
end
lambdaMultVax = 1 - lambdaMultVaxMat;

%% Initial population
mInit = popInit(: , 1);
fInit = popInit(: , 2);

% riskDistM(1:2 , 1:risk) = [1 , 0 , 0; 1 , 0 , 0];
% riskDistF(1:2 , 1:risk) = [1 , 0 , 0; 1 , 0 , 0];
% riskDistM(3:age , 1:risk) = [paramSet(85:98) , paramSet(99:112) , paramSet(113:126)];
% riskDistF(3:age , 1:risk) = [paramSet(127:140) , paramSet(141:154) , paramSet(155:168)];
% riskDist(: , : , 1) = riskDistM;
% riskDist(: , : , 2) = riskDistF;

MpopStruc = riskDistM;
FpopStruc = riskDistF;

mPop = zeros(age , risk);
fPop = mPop;

if startYear >= 1980;
    for i = 1 : age
        mPop(i , :) = MpopStruc(i, :).* mInit(i) ./ 1.25;
        fPop(i , :) = FpopStruc(i, :).* fInit(i) ./ 1.25;
    end
else
    for i = 1 : age
        mPop(i , :) = MpopStruc(i, :).* mInit(i) ./ (12*1.12);
        fPop(i , :) = FpopStruc(i, :).* fInit(i) ./ (12*1.12);
    end
end

dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age ,risk];
initPop = zeros(dim);
initPop(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop; % HIV-, acute infection, HPV Susceptible, no precancer, __, male
initPop(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop; % HIV-, acute infection, HPV Susceptible, no precancer, __, female
initPop_0 = initPop;

if (hivOn && (hivStartYear == startYear))
    initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = 0.005 / 2 .* ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3); % initial HIV infected male (age groups 4-6, med-high risk) (% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) .* (1 - 0.005 / 2); % moved to HIV infected
    initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) = 0.005 / 2 .*...
        initPop_0(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3); % initial HIV infected female (% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , 2 , 4 : 6 , 2 : 3) .* (1 - 0.005 / 2); % moved to HIV infected

    if hpvOn
        initPopHiv_0 = initPop;
        % HIV+ not infected by HPV
        % females
        initPop(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.3 .* ...
            initPopHiv_0(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
        % males
        initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.3 .* ...
            initPopHiv_0(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);

        % HIV+ infected by HPV
        % females
        initPop(3 , 2 , 2 , 1 , 1 , 2 , 4 : 6 , 1 : 3) = 0.7 .* ...
            initPopHiv_0(3 , 2 , 1 , 1 , 1 , 2 , 4 : 6 , 1 : 3);
        % males
        initPop(3 , 2 , 2 , 1 , 1 , 1 , 4 : 6 , 1 : 3) = 0.7 .* ...
            initPopHiv_0(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 1 : 3);
    end
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HIV infections.')

if (hpvOn && hivOn && (hivStartYear == startYear))
    infected = initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) * 0.20; % 20% initial HPV prevalence among age groups 4 - 9 (sexually active) (HIV-)
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - infected; % moved from HPV-

    % Omni-HPV type (transition rates weighted by estimated prevalence in population)
    initPop(1 , 1 , 2 , 1 , 1 , : , 4 : 9 , :) = infected; % moved to HPV+
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

if (hpvOn && ~hivOn) || (hpvOn && hivOn && (hivStartYear > startYear))
    infected = initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) * (0.2 * 0.9975); % initial HPV prevalence among age groups 4 - 9 (sexually active) (HIV-)
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - infected; % moved from HPV

    % Omni-HPV type (transition rates weighted by estimated prevalence in population)
    initPop(1 , 1 , 2 , 1 , 1 , : , 4 : 9 , :) = infected; % moved to HPV+
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

%% Simulation parameters
% lambdaMultVax = ones(age , 2);
fImm(1 : age) = 1; % all infected individuals who clear HPV get natural immunity

% Initialize vectors
s = 1 : timeStep : years + 1; % stepSize and steps calculated in loadUp.m
artDistMat = zeros(size(prod(dim) , 20)); % initialize artDistMat to track artDist over past 20 time steps
import java.util.LinkedList
artDistList = LinkedList();
popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8);
popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
newHiv = zeros(length(s) - 1 , gender , age , risk);
newHpv = zeros(length(s) - 1 , gender , disease , age , risk);
newImmHpv = newHpv;
newVaxHpv = newHpv;
newCC = zeros(length(s) - 1 , disease , hpvTypes , age);
ccDeath = newCC;
ccTreated = zeros(length(s) - 1 , disease , hpvTypes , age , 3); % 3 cancer stages: local, regional, distant
newScreen = zeros(length(s) - 1 , disease , viral , hpvTypes , hpvStates , risk , 2);
newTreatImm = newScreen;
newTreatHpv = newScreen;
newTreatHyst = newScreen;
vaxdSchool = zeros(length(s) - 1 , 1);
hivDeaths = zeros(length(s) - 1 , gender , age);
deaths = popVec;
artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
popVec(1 , :) = popIn;
tVec = linspace(startYear , endYear , size(popVec , 1));
k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);
artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0

%% Simulation
error = 0;
for i = 2 : length(s) - 1
    year = startYear + s(i) - 1;
%     currStep = round(s(i) * stepsPerYear);
    tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
    popIn = popVec(i - 1 , :);
    
    % Add HIV index cases at hivStartYear if not present in initial population
    if (hivOn && (hivStartYear > startYear) && (year == hivStartYear))
        % Initialize hiv cases in population
        popIn_init = popIn;
        
        % Create indices
        fromNonHivNonHpv = sort(toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1:gender , 4:6 , 1:risk))); 
        toHivNonHpv = sort(toInd(allcomb(3 , 2 , 1 , 1 , 1 , 1:gender , 4:6 , 1:risk)));
        fromNonHivHpv = sort(toInd(allcomb(1 , 1 , 2:4 , 1:hpvStates , 1 , 1:gender , 4:6 , 1:risk))); 
        toHivHpv = sort(toInd(allcomb(3 , 2 , 2:4 , 1:hpvStates , 1 , 1:gender , 4:6 , 1:risk)));

        % Distribute HIV infections (HPV-)        
        popIn(fromNonHivNonHpv) = (1 - 0.002) .* popIn_init(fromNonHivNonHpv);    % reduce non-HIV infected
        popIn(toHivNonHpv) = (0.002) .* popIn_init(fromNonHivNonHpv);    % increase HIV infected ( male/female, age groups 4-6, med-high risk) (% prevalence)

        % Distribute HIV infections (HPV+)
        popIn(fromNonHivHpv) = (1 - 0.001) .* popIn_init(fromNonHivHpv);    % reduce non-HIV infected
        popIn(toHivHpv) = (0.001) .* popIn_init(fromNonHivHpv);    % increase HIV infected ( male/female, age groups 4-6, med-high risk) (% prevalence)
    end

    if hpvOn
        % Progression/regression from initial HPV infection to
        % precancer stages and cervical cancer. Differential CC
        % detection by CC stage and HIV status/CD4 count.
        [~ , pop , newCC(i , : , : , :) , ccDeath(i , : , : , :) , ...
            ccTreated(i , : , : , : , :)] ...
            = ode4xtra(@(t , pop) ...
            hpvCCdet(t , pop , immuneInds , infInds , cin1Inds , ...
            cin2Inds , cin3Inds , normalInds , ccInds , ccRegInds , ccDistInds , ...
            ccTreatedInds , ccLocDetInds , ccRegDetInds , ccDistDetInds ,...
            kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
            kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf  ,...
            rNormal_Inf , hpv_hivClear , c3c2Mults , ...
            c2c1Mults , fImm , kRL , kDR , muCC , muCC_det , kCCDet , ...
            disease , age , hpvTypes , ...
            rImmuneHiv , hyst , hystInds , hystSusInds, OMEGA) , tspan , popIn);
        popIn = pop(end , :);
        if any(pop(end , :) <  0)
            disp('After hpv')
            error = 1;
            break
        end
        
        if (year >= hpvScreenStartYear)
            [dPop , newScreen(i , : , : , : , : , : , :) , ...
                newTreatImm(i , : , : , : , : , : , :) , ...
                newTreatHpv(i , : , : , : , : , : , :) , ...
                newTreatHyst(i , : , : , : , : , : , :)] ...
                = hpvScreen(popIn , disease , viral , hpvTypes , hpvStates , risk , ...
                screenYrs , screenAlgs , year , stepsPerYear , screenAgeAll , screenAgeS , ...
                noVaxNoScreen , noVaxToScreen , vaxNoScreen , vaxToScreen , ...
                noVaxToScreenTreatImm , vaxToScreenTreatImm , noVaxToScreenTreatHpv , ...
                vaxToScreenTreatHpv , noVaxToScreenHyst , vaxToScreenHyst , ...
                screenAlgorithm);
            pop(end , :) = pop(end , :) + dPop;
            popIn = pop(end , :); % for next module
            if any(pop(end , :) <  0)
                disp('After hpv screen')
                break
            end
        end
    end

    % HIV and HPV mixing and infection module. Protective effects of condom
    % coverage, circumcision, ART, PrEP (not currently used) are accounted for. 
    [~ , pop , newHpv(i , : , : , : , :) , newImmHpv(i , : , : , : , :) , ...
        newVaxHpv(i , : , : , : , :) , newHiv(i , : , : , :)] = ...
        ode4xtra(@(t , pop) mixInfect(t , pop , ...
        gar , perPartnerHpv , perPartnerHpv_lr , perPartnerHpv_nonV , maleActs , ...
        femaleActs , lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , ...
        epsR_vec , yr , circProtect , condProtect , condUse , actsPer , ...
        partnersM , partnersF , hpv_hivMult , hpvSus , hpvImm , hpvVaxd , ...
        hpvVaxdScreen , hpvVaxd2 , hpvImmVaxd2 , toHpv , toHpv_Imm , toHpv_Vax , ...
        toHpv_VaxScreen , toHpv_VaxNonV , toHpv_VaxNonVScreen , hivSus , toHiv , ...
        mCurr , fCurr , mCurrArt , fCurrArt , betaHIVF2M , betaHIVM2F , disease , ...
        viral , gender , age , risk , hpvStates , hpvTypes , hrInds , lrInds , ...
        hrlrInds , periods , startYear , stepsPerYear , year) , tspan , popIn);
    popIn = pop(end , :); % for next mixing and infection module
    if any(pop(end , :) < 0)
        disp('After mixInfect')
        error = 1;
        break
    end

    % HIV module, CD4 Progression, VL progression, ART initiation/dropout,
    % excess HIV mortality
    if (hivOn && (year >= hivStartYear))
        [~ , pop , hivDeaths(i , : , :) , artTreat] =...
            ode4xtra(@(t , pop) hiv2a(t , pop , vlAdvancer , artDist , muHIV , ...
            kCD4 ,  maxRateM1 , maxRateM2 , maxRateF1 , maxRateF2 , disease , ...
            viral , gender , age , risk , k , hivInds , ...
            stepsPerYear , year) , tspan , popIn);
        popIn = pop(end , :);
        artTreatTracker(i , : , : , : , :  ,:) = artTreat;
        artDistList.add(artTreat);
        if artDistList.size() >= stepsPerYear * 2
            artDistList.remove(); % remove CD4 and VL distribution info for people initiating ART more than 2 years ago
        end
        artDist = calcDist(artDistList , disease , viral , gender , age , ...
            risk);
        if any(pop(end , :) < 0)
            disp('After hiv')
            error = 1;
            break
        end
    end

    % Birth, aging, risk redistribution module
    [~ , pop , deaths(i , :) ] = ode4xtra(@(t , pop) ...
        bornAgeDieRisk(t , pop , year , ...
        gender , age , fertility , fertMat , fertMat2 , hivFertPosBirth ,...
        hivFertNegBirth , hivFertPosBirth2 , hivFertNegBirth2 , deathMat , circMat , circMat2 , ...
        MTCTRate , circStartYear , ageInd , riskInd , riskDist , ...
        stepsPerYear , currYear , screenAlgs{1}.screenAge , noVaxScreen , noVaxXscreen , ...
        vaxScreen , vaxXscreen , hpvScreenStartYear) , tspan , popIn);
    popIn = pop(end , :);
    if any(pop(end , :) < 0)
        disp('After bornAgeDieRisk')
        error = 1;
        break
    end
    
    if (year >= vaxStartYear)
        % HPV vaccination module- school-based vaccination regimen
        [dPop , vaxdSchool(i , :)] = hpvVaxSchool(popIn , k , ...
            disease , viral , risk , hpvTypes , hpvStates , ...
            periods , vaxG , vaxAge , vaxRate);
        pop(end , :) = pop(end , :) + dPop;
        if any(pop(end , :) < 0)
            disp('After hpvVaxSchool')
            break
        end
    end
    
    % add results to population vector
    popVec(i , :) = pop(end , :)';
end
popLast = popVec(end , :);
popVec = sparse(popVec); % compress population vectors

% allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
%        1 : periods , 2 , 4 : age , 1 : risk)); ...
%        toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%        1 : periods , 2 , 4 : age , 1 : risk))];
if error
    negSumLogL = -1;
    %ccInc = -1;
else
    negSumLogL = likeFunPtrnSrch(popVec , newCC , cinPos2014_obs , cinNeg2014_obs ,...
        hpv_hiv_2008_obs , hpv_hivNeg_2008_obs , hpv_hiv_obs , hpv_hivNeg_obs , ...
        hivPrevM_obs , hivPrevF_obs , hpv_hivM2008_obs , hpv_hivMNeg2008_obs , ...
        disease , viral , gender , age , risk , hpvTypes , hpvStates , ...
        periods , startYear , stepsPerYear);

    %ccInc = annlz(sum(sum(sum(newCC(end , ':' , : , 4 : age),2),3),4)) ./ ...
    %    (annlz(sum(popVec(end , allF) , 2) ./ stepsPerYear))* (10 ^ 5);
end

% negSumLogL
% pathModifier = 'toNow_081219calib_05Aug19_0_1424';
% savdir = [pwd , '\HHCoM_Results\'];
% save(fullfile(savdir , pathModifier) , 'tVec' ,  'popVec' , 'newHiv' ,...
%     'newImmHpv' , 'newVaxHpv' , 'newHpv' , 'hivDeaths' , 'deaths' , ...
%     'vaxdSchool' , 'newScreen' , 'newTreatImm' , 'newTreatHpv' , 'newTreatHyst' , ...
%     'newCC' , 'artTreatTracker' , 'artDist' , 'artDistList' , ... 
%     'startYear' , 'endYear' , 'popLast');
% showResults(pathModifier)

