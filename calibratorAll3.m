function [negSumLogL] = calibratorAll3(pIdx , paramsSub , paramSet)

%% Load parameters
paramDir = [pwd ,'/Params/'];
load([paramDir,'settings'])
load([paramDir,'popData'])
load([paramDir,'HIVParams'])
load([paramDir,'general'])
load([paramDir,'mixInfectParams'])
load([paramDir,'vlBeta'])
load([paramDir,'hpvData'])
load([paramDir,'cost_weights'])
load([paramDir,'calibData'])
load([paramDir,'mixInfectIndices'])
load([paramDir,'hivIndices'])
load([paramDir,'hpvIndices'])
load([paramDir,'hpvTreatIndices'])
load([paramDir,'ageRiskInds'])
load([paramDir,'vlAdvancer'])
load([paramDir,'fertMat'])
load([paramDir,'hivFertMats'])
load([paramDir,'fertMat2'])
load([paramDir,'hivFertMats2'])
load([paramDir,'circMat'])
load([paramDir,'circMat2'])
load([paramDir,'deathMat'])

% Model specs
hyst = 'off';
hpvOn = 1;
hivOn = 1;

% Use calibrated parameters
load([paramDir,'calibratedParams'])
perPartnerHpv = 0.0045;
%paramSetMatrix = load([paramDir,'params_calib_22Feb19.dat']);
%paramSet = paramSetMatrix(:,985);

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

% Helper functions
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

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
    infected = initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) * (0.05 * 0.9975); % initial HPV prevalence among age groups 4 - 9 (sexually active) (HIV-)
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
        initPop_0(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - infected; % moved from HPV

    % Omni-HPV type (transition rates weighted by estimated prevalence in population)
    initPop(1 , 1 , 2 , 1 , 1 , : , 4 : 9 , :) = infected; % moved to HPV+
end
assert(~any(initPop(:) < 0) , 'Some compartments negative after seeding HPV infections.')

%% Calibration parameters (FEED FROM LHS)
if any(1 == pIdx)
    idx = find(1 == pIdx);
    rowL = paramsSub{idx}.length/3;
    rl = paramsSub{idx}.inds(1:rowL);
    rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    partnersM(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
    partnersM(3:age , 1:risk) = [paramSet(rl) , paramSet(rm) , paramSet(rh)];
    partnersM(10:age , 3) = ones(7 , 1);
end
if any(2 == pIdx)
    idx = find(1 == pIdx);
    rowL = paramsSub{idx}.length/3;
    rl = paramsSub{idx}.inds(1:rowL);
    rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    partnersF(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
    partnersF(3:age , 1:risk) = [paramSet(rl) , paramSet(rm) , paramSet(rh)];
    partnersF(10:age , 3) = ones(7 , 1);
end
% partnersM(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
% partnersF(1:2 , 1:risk) = ones(2 , risk) .* 0.00001;
% partnersM(3:age , 1:risk) = [paramSet(1:14) , paramSet(15:28) , paramSet(29:42)];
% partnersF(3:age , 1:risk) = [paramSet(43:56) , paramSet(57:70) , paramSet(71:84)];
% partnersM(10:age , 3) = ones(7 , 1);
% partnersF(10:age , 3) = ones(7 , 1);
% condUse = paramSet(169);
if any(6 == pIdx);
    idx = find(6 == pIdx);
    epsA = paramSet(paramsSub{idx}.inds(:));
end
if any(7 == pIdx);
    idx = find(7 == pIdx);
    epsR = paramSet(paramsSub{idx}.inds(:));
end
% epsA = paramSet(170:172);
% epsR = paramSet(173:175);
if any(8 == pIdx)
    idx = find(8 == pIdx);
    rowL = paramsSub{idx}.length/3;
    rl = paramsSub{idx}.inds(1:rowL);
    rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    maleActs(1:2 , 1:risk) = zeros(2 , risk);
    maleActs(3:age , 1:risk) = [paramSet(rl) , paramSet(rm) , paramSet(rh)];
end
if any(9 == pIdx)
    idx = find(9 == pIdx);
    rowL = paramsSub{idx}.length/3;
    rl = paramsSub{idx}.inds(1:rowL);
    rm = paramsSub{idx}.inds(rowL+1 : rowL*2);
    rh = paramsSub{idx}.inds(rowL*2+1 : rowL*3);
    femaleActs(1:2 , 1:risk) = zeros(2 , risk);
    femaleActs(3:age , 1: risk) = [paramSet(rl) , paramSet(rm) , paramSet(rh)];
end
% maleActs(1:2 , 1:risk) = zeros(2 , risk);
% femaleActs(1:2 , 1:risk) = zeros(2 , risk);
% maleActs(3:age , 1:risk) = [paramSet(176:189) , paramSet(190:203) , paramSet(204:217)];
% femaleActs(3:age , 1:risk) = [paramSet(218:231) , paramSet(232:245) , paramSet(246:259)];
if any(10 == pIdx)
    idx = find(10 == pIdx);
    perPartnerHpv = paramSet(paramsSub{idx}.inds(:));
end
% perPartnerHpv = paramSet(260);
% perPartnerHpv_lr = paramSet(261);
% perPartnerHpv_nonV = paramSet(262);
% hpv_hivMult = paramSet(263:266) .* 2.0;
% rNormal_Inf_orig = rNormal_Inf;
% rNormal_Inf = paramSet(267) .* rNormal_Inf_orig;
% hpv_hivClear = paramSet(283:286) .* 0.5;
% c3c2Mults = paramSet(287:290) .* 2.0;
% c2c1Mults = paramSet(291:294) .* 2.0;
% kCCDet = paramSet(295:297);
% lambdaMultImm = paramSet(298:313);
% maxRateM_vec = paramSet(314:315);
% maxRateF_vec = paramSet(316:317);
% artHpvMult = paramSet(318) .* 2.0;
% kCD4(1,:,:) = [paramSet(319:323) , paramSet(324:328) , paramSet(329:333) , paramSet(334:338)];
% kCD4(2,:,:) = [paramSet(339:343) , paramSet(344:348) , paramSet(349:353) , paramSet(354:358)];
% kVl(1,:,:) = [paramSet(359:363) , paramSet(364:368) , paramSet(369:373) , paramSet(374:378)];
% kVl(2,:,:) = [paramSet(379:383) , paramSet(384:388) , paramSet(389:393) , paramSet(394:398)];
% 
% % maxRateM_vec = [0.4 , 0.4];
% % maxRateF_vec = [0.5 , 0.5];
% rImmuneHiv = hpv_hivClear;
% % lambdaMultImm = ones(age,1);
% % artHpvMult = 1.0;
% 
% for a = 1 : age
%     betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
%     betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
% end
% betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
% betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl

maxRateM1 = maxRateM_vec(1);
maxRateM2 = maxRateM_vec(2);
maxRateF1 = maxRateF_vec(1);
maxRateF2 = maxRateF_vec(2);

% epsA_vec = cell(size(yr , 1) - 1, 1); % save data over time interval in a cell array
% epsR_vec = cell(size(yr , 1) - 1, 1);
% for i = 1 : size(yr , 1) - 1          % interpolate epsA/epsR values at steps within period
%     period = [yr(i) , yr(i + 1)];
%     epsA_vec{i} = interp1(period , epsA(i : i + 1 , 1) , ...
%         yr(i) : timeStep : yr(i + 1));
%     epsR_vec{i} = interp1(period , epsR(i : i + 1 , 1) , ...
%         yr(i) : timeStep : yr(i + 1));
% end
epsA_vec = epsA;
epsR_vec = epsR;

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
vaxd = zeros(length(s) - 1 , 1);
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
    currStep = round(s(i) * stepsPerYear);
    tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
    popIn = popVec(i - 1 , :);
    
    % Add HIV index cases at hivStartYear if not present in initial population
    if (hivOn && (hivStartYear > startYear) && (year == hivStartYear))
        % Initialize hiv cases in population
        popIn_init = popIn;
        
        % Create indices
        fromNonHivNonHpv = sort(toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1:gender , 4:6 , 2:3))); 
        toHivNonHpv = sort(toInd(allcomb(3 , 2 , 1 , 1 , 1 , 1:gender , 4:6 , 2:3)));
        fromNonHivHpv = sort(toInd(allcomb(1 , 1 , 2:4 , 1:hpvStates , 1 , 1:gender , 4:6 , 1:3))); 
        toHivHpv = sort(toInd(allcomb(3 , 2 , 2:4 , 1:hpvStates , 1 , 1:gender , 4:6 , 1:3)));

        % Distribute HIV infections (HPV-)        
        popIn(fromNonHivNonHpv) = (1 - 0.007) .* popIn_init(fromNonHivNonHpv);    % reduce non-HIV infected
        popIn(toHivNonHpv) = (0.007) .* popIn_init(fromNonHivNonHpv);    % increase HIV infected ( male/female, age groups 4-6, med-high risk) (% prevalence)

        % Distribute HIV infections (HPV+)
        popIn(fromNonHivHpv) = (1 - 0.006) .* popIn_init(fromNonHivHpv);    % reduce non-HIV infected
        popIn(toHivHpv) = (0.006) .* popIn_init(fromNonHivHpv);    % increase HIV infected ( male/female, age groups 4-6, med-high risk) (% prevalence)
    end

    if hpvOn
        hystOption = 'on';
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
            disease , viral , age , hpvTypes , ...
            rImmuneHiv , vaccinated , hystOption) , tspan , popIn);
        popIn = pop(end , :);
        if any(pop(end , :) <  0)
            disp('After hpv')
            error = 1;
            break
        end

    end

    % HIV and HPV mixing and infection module. Protective effects of condom
    % coverage, circumcision, ART, PrEP (not currently used) are accounted for. 
    [~ , pop , newHpv(i , : , : , : , :) , newImmHpv(i , : , : , : , :) , ...
        newVaxHpv(i , : , : , : , :) , newHiv(i , : , : , :)] = ...
        ode4xtra(@(t , pop) mixInfect(t , pop , currStep , ...
        gar , perPartnerHpv , perPartnerHpv_lr , perPartnerHpv_nonV , maleActs , femaleActs , ...
        lambdaMultImm , lambdaMultVax , artHpvMult , epsA_vec , epsR_vec , yr , modelYr1 , ...
        circProtect , condProtect , condUse , actsPer , partnersM , partnersF , ...
        hpv_hivMult , hpvSus , hpvImm , toHpv_Imm , hpvVaxd , hpvVaxd2 , hpvImmVaxd2 ,  toHpv , toHpv_ImmVax , ...
        toHpv_ImmVaxNonV , hivSus , toHiv , mCurr , fCurr , mCurrArt , fCurrArt , ...
        betaHIVF2M , betaHIVM2F , disease , viral , gender , age , risk , hpvStates , hpvTypes , ...
        hrInds , lrInds , hrlrInds , periods , startYear , stepsPerYear , year) , tspan , popIn);
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
        artTreatTracker(i , : , : , : , :  ,:) = artTreat;
        artDistList.add(artTreat);
        if artDistList.size() >= stepsPerYear * 2
            artDistList.remove(); % remove CD4 and VL distribution info for people initiating ART more than 2 years ago
        end
        artDist = calcDist(artDistList , disease , viral , gender , age , ...
            risk);
        popIn = pop(end , :);
        if any(pop(end , :) < 0)
            disp('After hiv')
            error = 1;
            break
        end
    end

    % Birth, aging, risk redistribution module
    [~ , pop , deaths(i , :) ] = ode4xtra(@(t , pop) ...
        bornAgeDieRisk(t , pop , year , currStep ,...
        gender , age , risk , fertility , fertMat , fertMat2 , hivFertPosBirth ,...
        hivFertNegBirth , hivFertPosBirth2 , hivFertNegBirth2 , deathMat , circMat , circMat2 , ...
        MTCTRate , circStartYear , ageInd , riskInd , riskDist , startYear , ...
        endYear, stepsPerYear , currYear) , tspan , popIn);
    if any(pop(end , :) < 0)
        disp('After bornAgeDieRisk')
        error = 1;
        break
    end
    % add results to population vector
    popVec(i , :) = pop(end , :)';
end
popLast = popVec(end , :);
popVec = sparse(popVec); % compress population vectors

%allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
%        1 : periods , 2 , 4 : age , 1 : risk)); ...
%        toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%        1 : periods , 2 , 4 : age , 1 : risk))];
if error
    negSumLogL = -1;
    %ccInc = -1;
else
    negSumLogL = likeFun(popVec , newCC , cinPos2014_obs , cinNeg2014_obs ,...
        hpv_hiv_2008_obs , hpv_hivNeg_2008_obs , hpv_hiv_obs , hpv_hivNeg_obs , ...
        hivPrevM_obs , hivPrevF_obs , disease , viral , gender , age , risk , ...
        hpvTypes , hpvStates , periods , startYear , stepsPerYear);

    %ccInc = annlz(sum(sum(sum(newCC(end , ':' , : , 4 : age),2),3),4)) ./ ...
    %    (annlz(sum(popVec(end , allF) , 2) ./ stepsPerYear))* (10 ^ 5);
end

% savdir = [pwd , '\HHCoM_Results\'];%'H:\HHCoM_Results';
% save(fullfile(savdir , '02252019_paramSet985') , 'tVec' ,  'popVec' , 'newHiv' ,...
%     'newImmHpv' , 'newVaxHpv' , 'newHpv' , 'hivDeaths' , ...
%     'deaths' , 'newCC' , 'artTreatTracker' , 'artDist' , 'artDistList' , ... 
%     'startYear' , 'endYear' , 'popLast');
% disp(' ')
% disp('Simulation complete.')

