close all;clear all;clc
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'])
load([paramDir,'mixInfectIndices'])
load([paramDir,'vlAdvancer'])
load([paramDir,'fertMat'])
load([paramDir,'hivFertMats'])
load([paramDir,'fertMat2'])
load([paramDir,'hivFertMats2'])
load([paramDir,'deathMat'])
load([paramDir,'circMat'])
load([paramDir,'vaxer'])
load([paramDir,'mixInfectParams'])
load([paramDir,'popData'])
load([paramDir,'HIVParams'])
load([paramDir,'hivIndices'])
load([paramDir,'hpvIndices'])
load([paramDir,'ager'])
load([paramDir,'vlBeta'])
load([paramDir,'hpvTreatIndices'])
load([paramDir,'calibParams'])
load([paramDir,'vaxInds'])
% load([paramDir,'settings'])
load([paramDir,'hpvData'])
load([paramDir,'calibData'])
load([paramDir,'calibInitParams'])
load([paramDir , 'ageRiskInds'])
load([paramDir ,'cost_weights'])
import java.util.LinkedList
vaxerAger = ager;
vaxRate = 0;
startYear = 1980;
%% Initial population
load([paramDir , 'popData'])
load([paramDir,'calibInitParams'])
% load('initPop')
% simulation
mInit = popInit(: , 1);

fInit = popInit(: , 2);


% test!!!!
riskDistF = riskDistM;
partnersF = partnersM;

MpopStruc = riskDistM;
FpopStruc = riskDistF;

mPop = zeros(age , risk);
fPop = mPop;

for i = 1 : age
    mPop(i , :) = MpopStruc(i, :).* mInit(i) ./ 1.25;
    fPop(i , :) = FpopStruc(i, :).* fInit(i) ./ 1.25;
end

dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age ,risk];
initPop = zeros(dim);
initPop(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop;
initPop(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop;
initPop_0 = initPop;

% Intervention start years
circStartYear = 1990;
vaxStartYear = 2017;    
load([paramDir , 'HPV_calib16.dat']) 
distWeight = [0.7 , 0.2 , 0.1];
vaxMat = ager .* 0;
maxRateM_vec = [0.45 , 0.45];% maxRateM_arr{sim};
maxRateF_vec = [0.6 , 0.6];% maxRateF_arr{sim};
hpv_hivMult = sum(bsxfun(@times , hpv_hivMult , distWeight) , 2);
% for quickCalibrateModel

kCin1_Inf = HPV_calib16(1) .* kCin1_Inf;
rNormal_Inf = HPV_calib16(2) .* rNormal_Inf;
kCC_Cin3 = HPV_calib16(3) .* kCC_Cin3;
kCin3_Cin2 = HPV_calib16(4) .* kCin3_Cin2;

rImmuneHiv = HPV_calib16(5 : 8);
% c3c2Mults = HPV_calib16(9 : 12);
% c2c1Mults = HPV_calib16(13 : 16);
artHpvMult = HPV_calib16(9);
perPartnerHpv = HPV_calib16(10);
lambdaMultImm = HPV_calib16(11 : 26);

partnersM(4 , :) = partnersM(4 , :) .* [1.25 , 1.75 , 1.75];
partnersF(4 , :) = partnersF(4 , :) .* [1.25 , 1.75 , 1.75];
partnersM(5 , :) = partnersM(5 , :) .* [1.25 , 1.5 , 1.75];
partnersF(5 , :) = partnersF(5 , :) .* [1.25 , 1.5 , 1.75];

femaleActs(4 : 5 , :) = femaleActs(4 : 5 , :) .* 1.2 ;
femaleActs(6 : 10 , :) = femaleActs(6 : 10 , :) .* 0.9;
maleActs(4 : 5 , :) = maleActs(4 : 5 , :);

% for i = 0 : 2
%     maleActs(: , i + 1) = maleActs(: , i + 1) .* HPV_calib16(27 + i);
%     femaleActs(: , i + 1) = femaleActs(: , i + 1) .* HPV_calib16(30 + i);
% end
% 
% for i = 0 : 2
%     partnersM(: , i + 1) = partnersM(: , i + 1) .* HPV_calib16(33 + i);
%     partnersF(: , i + 1) = partnersF(: , i + 1) .* HPV_calib16(36 + i);
% end

for a = 1 : age
    betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
    betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
end
betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl

kCCDet = min(kCCDet .* 12 , 0.99); % convert monthly to yearly rate

riskDistF = riskDistM;
riskDist(: , : , 1) = riskDistM;
riskDist(: , : , 2) = riskDistF;
save([paramDir , 'calibratedParams'])