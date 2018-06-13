function quickCalibrateModel()
close all; clear all; clc
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'])
load([paramDir,'mixInfectIndices'])
load([paramDir,'vlAdvancer'])
load([paramDir,'fertMat'])
load([paramDir,'hivFertMats'])
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
load([paramDir,'settings'])
load([paramDir,'hpvData'])
load([paramDir,'calibData'])
% load([paramDir,'calibInitParams'])
import java.util.LinkedList

w = ones(4 , 1) ./ 4;
kCC_Cin3_Orig = kCC_Cin3;
kCin2_Cin3_Orig = kCin2_Cin3;
kCin2_Cin1_Orig = kCin2_Cin1;
kCin1_Cin2_Orig = kCin1_Cin2;
kCin3_Cin2_Orig = kCin3_Cin2;
rNormal_Inf_Orig = rNormal_Inf;

for i = 1 : 3
    rNormal_Inf(: , i) = conv(rNormal_Inf_Orig(: , i) , w , 'same');
    rNormal_Inf(end - 1 : end , i) = rNormal_Inf_Orig(end - 1 : end , i);
    kCC_Cin3(: , i) = conv(kCC_Cin3_Orig(: , i) , w , 'same');
    kCC_Cin3(end - 1 : end , i) = kCC_Cin3_Orig(end - 1 : end , i);
%     kCin2_Cin3(: , i) = conv(kCin2_Cin3_Orig(: , i) , w , 'same');
%     kCin2_Cin3(end - 1 : end , i) = kCin3_Cin2_Orig(end - 1 : end , i);
    kCin3_Cin2(: , i) = conv(kCin3_Cin2_Orig(: , i) , w , 'same');
    kCin3_Cin2(end - 1 : end , i) = kCin3_Cin2_Orig(end - 1 : end , i);
    kCin1_Cin2(: , i) = conv(kCin1_Cin2_Orig(: , i) , w , 'same');
    kCin1_Cin2(end - 1 : end , i) = kCin1_Cin2_Orig(end - 1 : end , i);
    kCin2_Cin1(: , i) = conv(kCin2_Cin1_Orig(: , i) , w , 'same');
    kCin2_Cin1(end - 1 : end , i) = kCin2_Cin1_Orig(end - 1 : end , i);
end
muCC = min(muCC .* 12 , 0.99); % convert cervical cancer mortality rate from monthly to yearly
%     fImm(4 : age) = 1; % RR(0.75; 0.5 , 0.92) fraction fully protected by immunity based on RR of natural immunity (Beachler, 2017)
perPartnerHpv = 0.08; % high risk HPV transmission risk per month
rImmuneHiv = 1 ./ hpv_hivClear;
lambdaMultImm(1 : 4) = 1 - 0.01;
lambdaMultImm(5 : 10) = 1 - logspace(log10(0.01) , log10(0.1) , 6);
lambdaMultImm(11 : age) = lambdaMultImm(10);
lambdaMultVax = ones(age , 2);

artHpvMult = 1;%hpv_hivMult(1 , 1);
perPartnerHpv = 0.015;
perPartnerHpv_lr = 0.08;
perPartnerHpv_nonV = 0.08;
rImmuneHiv = 2 ./ hpv_hivClear;

save([paramDir , 'calibInitParams']);
%% Continue from last calibration
load([paramDir,'calibInitParams'])
load([paramDir,'HPV_calib12.dat'])
rImmuneHiv = 1 ./ hpv_hivClear; 

lambdaMultImm = HPV_calib12(24 : 39);
hpv_hivClear = HPV_calib12(40 : 43);
hpvClearMult = HPV_calib12(44 : 47);
perPartnerHpv_lr = HPV_calib12(48);%0.1;
perPartnerHpv_nonV = HPV_calib12(49); %0.1;

% Weight HPV transitions according to type distribution

distWeight = [0.7 , 0.2 , 0.1];
kInf_Cin1 = sum(bsxfun(@times , kInf_Cin1 , distWeight) , 2);
kCin1_Cin2 = sum(bsxfun(@times , kCin1_Cin2 , distWeight) , 2);
kCin2_Cin3 = sum(bsxfun(@times , kCin2_Cin3 , distWeight) , 2);
kCin2_Cin1 = sum(bsxfun(@times , kCin2_Cin1 , distWeight) , 2);
kCin3_Cin2 = sum(bsxfun(@times , kCin3_Cin2 , distWeight) , 2);
kCC_Cin3 = sum(bsxfun(@times , kCC_Cin3 , distWeight) , 2);
kCin1_Inf = sum(bsxfun(@times , kCin1_Inf , distWeight) , 2);
rNormal_Inf = sum(bsxfun(@times , rNormal_Inf , distWeight) , 2);

vaxMat = ager .* 0;
maxRateM_vec = [0.60 , 0.60];% maxRateM_arr{sim};
maxRateF_vec = [0.70 , 0.70];% maxRateF_arr{sim};

maxRateM1 = 1 - exp(-maxRateM_vec(1));
maxRateM2 = 1 - exp(-maxRateM_vec(2));
maxRateF1 = 1 - exp(-maxRateF_vec(1));
maxRateF2 = 1 - exp(-maxRateF_vec(2));
load([paramDir,'fertMat'])
load([paramDir,'hivFertMats'])
load([paramDir,'fertMat2'])
load([paramDir,'hivFertMats2'])
lambdaMultVax = ones(age , 2);


partnersM(4 , :) = partnersM(4 , :) .* [1.25 , 1.75 , 1.75];
partnersF(4 , :) = partnersF(4 , :) .* [1.25 , 1.75 , 1.75];
partnersM(5 , :) = partnersM(5 , :) .* [1.25 , 1.5 , 1.75];
partnersF(5 , :) = partnersF(5 , :) .* [1.25 , 1.5 , 1.75];

femaleActs(4 : 5 , :) = femaleActs(4 : 5 , :) .* 1.2 ;
femaleActs(6 : 10 , :) = femaleActs(6 : 10 , :) .* 0.9;
maleActs(4 : 5 , :) = maleActs(4 : 5 , :);
for a = 1 : age
    betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
    betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
end
betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl

artHpvMult = 1;
hpv_hivMult = sum(bsxfun(@times , hpv_hivMult , distWeight) , 2);
%%
kCin1_Inf_Mult = 1;
rNormal_InfMult = 1;
kCC_Cin3_Mult = 1;
kCin3_Cin2_Mult = 1;
maleActsMult = [1 ; 1 ; 1];
femaleActsMult = [1 ; 1; 1];
partnersM_Mult = [1; 1; 1];
partnersF_Mult = [1; 1; 1];
hpvClearMult = 0.1 .* hpv_hivClear;

initParams = [kCin1_Inf_Mult;
    rNormal_InfMult;
    kCC_Cin3_Mult;
    rImmuneHiv; % HPV immunity clearance multiplier for HIV+
    c3c2Mults; % CIN2 to CIN3 progression multiplier for HIV+
    c2c1Mults; % CIN1 to CIN2 progression multiplier for HIV+
    artHpvMult; % HPV acquisition multiplier for HIV+
    perPartnerHpv;
    lambdaMultImm;
    hpv_hivClear;
    kCin3_Cin2_Mult;
    maleActsMult;
    femaleActsMult
    partnersM_Mult;
    partnersF_Mult];

lb = initParams .* 0.5;
ub = initParams .* 1.5; 

% max transition rate between CC to CIN is 0.5

ub(1) = min(0.3 / max(kCin1_Inf) , 1.5);
ub(2) = min(0.8 / max(kInf_Cin1) ,  1.5);
ub(3) = min(0.5 / max(kCC_Cin3) ,  1.5);
ub(4) = min(0.6 / max(kCin3_Cin2) , 1.5);

lb(5 : 6) = 1; % CIN progression multiplier for HIV+ >=1
lb(7) = 1; % HPV acquisition multiplier for HIV+ on ART >= 1
% ub(22) = max(hpv_hivMult(1 , :)); % max HPV acquisition multiplier for HIV+ on ART <= max acquisition multiplier for HIV+
ub(18) = 1.5; % max per partner transmission probability multiplier for HPV
ub(35 : 38) = 1; % HPV clearance multipliers <= 1 for HIV-positive

% A = zeros(length(initParams));
% for i = 0 : 2
%     A(i + 1 , 14 + i) = 1;
%     A(i + 1 , 15 + i) = -1;
% %     
% %     A(i + 4 , 18 + i) = 1;
% %     A(i + 4 , 19 + i) = -1;
%     
%     A(i + 7 , 10 + i) = 1;
%     A(i + 7 , 11 + i) = -1;
%     
%     A(i + 10 , 40 + i) = -1;
%     A(i + 10 , 41 + i) = 1;
% end
% 
% b = zeros(length(initParams) , 1);

options = psoptimset('UseParallel' , true , 'Cache' , 'on' ,...
    'CacheTol' , 0.1 , 'CompletePoll' , 'on' , 'TolMesh' , 0.1, ...
    'Display','iter','PlotFcn',@psplotbestf);
x = patternsearch(@calibrator, initParams , [] , [] , [] , [] , lb , ub , [] , options);
%%
file = 'HPV_calib13.dat';
paramDir = [pwd , '\Params\'];
csvwrite([paramDir, file] , x)
