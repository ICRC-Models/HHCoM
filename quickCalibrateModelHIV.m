function quickCalibrateModelHIV()
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

% Decreasing transition rates for young ages
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
muCC_det = min(muCC_det .* 12 , 0.99); % convert cervical cancer mortality rate from monthly to yearly
%     fImm(4 : age) = 1; % RR(0.75; 0.5 , 0.92) fraction fully protected by immunity based on RR of natural immunity (Beachler, 2017)
rImmuneHiv = 1 ./ hpv_hivClear;
rImmuneHiv = 2 ./ hpv_hivClear;
lambdaMultImm(1 : 4) = 1 - 0.01;
lambdaMultImm(5 : 10) = 1 - logspace(log10(0.01) , log10(0.1) , 6);
lambdaMultImm(11 : age) = lambdaMultImm(10);
lambdaMultVax = ones(age , 2);

artHpvMult = hpv_hivMult(1 , 1);
perPartnerHpv = 0.08; % high risk HPV transmission risk per month
perPartnerHpv_lr = 0.08;
perPartnerHpv_nonV = 0.08;

save([paramDir , 'calibInitParams']);
%% Continue from last calibration
load([paramDir,'calibInitParams'])
load([paramDir,'HPV_calib12.dat'])
for i = 1 : 3
    kCin1_Inf(: , i) = HPV_calib12(i) .* kCin1_Inf(: , i);
    rNormal_Inf(: , i) = HPV_calib12(3 + i) .* rNormal_Inf(: , i);
    kCC_Cin3(: , i) = HPV_calib12(6 + i) .* kCC_Cin3(: , i);
    kCin3_Cin2(: , i) = HPV_calib12(49 + i) .* kCin3_Cin2(: , i); 
end
% kCC_Cin3(: , 2) = kCC_Cin3(: , 3);
% kCin3_Cin2(: , 3) = 1.5 .* kCin3_Cin2(: , 3);
% kCC_Cin3(: , 2 : 3) = kCC_Cin3(: , 2 : 3) .* 1.25;
% rNormal_Inf(: , 2 : 3) = 0.5 .* rNormal_Inf(: , 2 : 3);

%commented out for testing
% for i = 0 : 2
%     maleActs(: , i + 1) = maleActs(: , i + 1) .* HPV_calib12(53 + i);
%     femaleActs(: , i + 1) = femaleActs(: , i + 1) .* HPV_calib12(56 + i);
% end

for a = 1 : age
    betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
    betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
end
betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl
rImmuneHiv = HPV_calib12(10 : 13);
% c3c2Mults = HPV_calib12(14 : 17);
% c2c1Mults = HPV_calib12(18 : 21);
perPartnerHpv= HPV_calib12(23);
lambdaMultImm = HPV_calib12(24 : 39);
hpv_hivClear = HPV_calib12(40 : 43);
hpvClearMult = HPV_calib12(44 : 47);
perPartnerHpv_lr = HPV_calib12(48);%0.1;
perPartnerHpv_nonV = HPV_calib12(49); %0.1;
% artHpvMult = HPV_calib12(22);

% Weight HPV transitions according to type distribution

distWeight = [0.6 , 0.3 , 0.1];
kInf_Cin1 = sum(bsxfun(@times , kInf_Cin1 , distWeight) , 2);
kCin1_Cin2 = sum(bsxfun(@times , kCin1_Cin2 , distWeight) , 2);
kCin2_Cin3 = sum(bsxfun(@times , kCin2_Cin3 , distWeight) , 2);
kCin2_Cin1 = sum(bsxfun(@times , kCin2_Cin1 , distWeight) , 2);
kCin3_Cin2 = sum(bsxfun(@times , kCin3_Cin2 , distWeight) , 2);
kCC_Cin3 = sum(bsxfun(@times , kCC_Cin3 , distWeight) , 2);
kCin1_Inf = sum(bsxfun(@times , kCin1_Inf , distWeight) , 2);
rNormal_Inf = sum(bsxfun(@times , rNormal_Inf , distWeight) , 2);

vaxMat = ager .* 0;
maxRateM_vec = [0.45 , 0.45];% maxRateM_arr{sim};
maxRateF_vec = [0.65 , 0.65];% maxRateF_arr{sim};

maxRateM1 = 1 - exp(-maxRateM_vec(1));
maxRateM2 = 1 - exp(-maxRateM_vec(2));
maxRateF1 = 1 - exp(-maxRateF_vec(1));
maxRateF2 = 1 - exp(-maxRateF_vec(2));
%%
partnersF = partnersM;

initParams = [partnersM(: , 1);
    partnersM(: , 2);
    partnersM(: , 3);
    partnersF(: , 1);
    partnersF(: , 2);
    partnersF(: , 3);
    maleActs(: , 1);
    maleActs(: , 2);
    maleActs(: , 3);
    femaleActs(: , 1);
    femaleActs(: , 2);
    femaleActs(: , 3)];

lb = initParams .* 0.5;
ub = initParams .* 1.5; 

% max transition rate between CC to CIN is 0.5

A = zeros(length(initParams));
for i = 1 : size(partnersM , 1)
    % partners
    % Fewer partners in low risk group than medium risk group (male)
    A(i , i) = 1;
    A(i , size(partnersM , 1) + i) = -1;
    % Fewer partners in medium risk group than high risk group (male)
    A(i + size(partnersM , 1) , size(partnersM , 1) + i) = 1;
    A(i + size(partnersM , 1) , 2 * size(partnersM , 1) + i - 1) = -1;
    
    % Fewer partners in low risk group than medium risk group (female)
    A(i + 2 * size(partnersM , 1) , 2 * size(partnersM , 1) + i) = 1;
    A(i + 2 * size(partnersM , 1) , 3 * size(partnersM , 1) + i - 1) = -1;
    % Fewer partners in medium risk group than high risk group (female)
    A(i + 3 * size(partnersM , 1) , 3 * size(partnersM , 1) + i) = 1;
    A(i + 3 * size(partnersM , 1) , 4 * size(partnersM , 1) + i - 1) = -1;
    
    % acts per
    
    % more acts in low risk than medium risk (male)
    A(i + 4 * size(partnersM , 1) , 4 * size(partnersM , 1) + i) = -1;
    A(i + 4 * size(partnersM , 1) , 5 * size(partnersM , 1) + i - 1) = 1;
    % more acts in medium risk than high risk (male)
    A(i + 5 * size(partnersM , 1) , 5 * size(partnersM , 1) + i) = -1;
    A(i + 5 * size(partnersM , 1) , 6 * size(partnersM , 1) + i - 1) = 1;
    % more acts in low risk than medium risk (female)
    A(i + 6 * size(partnersM , 1) , 6 * size(partnersM , 1) + i) = -1;
    A(i + 6 * size(partnersM , 1) , 7 * size(partnersM , 1) + i - 1) = 1;
    % more acts in medium risk than high risk (female)
    A(i + 7 * size(partnersM , 1) , 7 * size(partnersM , 1) + i) = -1;
    A(i + 7 * size(partnersM , 1) , 8 * size(partnersM , 1) + i - 1) = 1;
end

b = zeros(length(initParams) , 1);

options = psoptimset('UseParallel' , true , 'Cache' , 'on' ,...
    'CacheTol' , 0.1 , 'CompletePoll' , 'on' , 'TolMesh' , 0.1, ...
    'Display','iter','PlotFcn',@psplotbestf);
x = patternsearch(@calibratorHIV, initParams , [] , [] , [] , [] , lb , ub , [] , options);
%%
file = 'HIV_calib.dat';
paramDir = [pwd , '\Params\'];
csvwrite([paramDir, file] , x)
