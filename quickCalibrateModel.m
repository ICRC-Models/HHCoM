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
    rNormal_Inf(end - 1 : end , i) = conv(rNormal_Inf_Orig(end - 1 : end , i) , w , 'same');
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
% rNormal_Inf(: , 1) = rNormal_Inf(: , 1) .* 1.2;
% rNormal_Inf(: , 2) = rNormal_Inf(: , 2) .* 0.8;
% rNormal_Inf(: , 3) = rNormal_Inf(: , 3) .* 0.9;
rNormal_Inf = rNormal_Inf;
% c2c1Mults = 1.5 .* c2c1Mults;
% kCin1_Cin2(1 : end , :) = 1.2 * kCin1_Cin2(1 : end , :);
% kCin2_Cin1(6 : end , :) = 1.25 * kCin2_Cin1(6 : end , :); 
% kCin3_Cin2(10 : end , :) = 1.5 * kCin3_Cin2(10 : end , :);
% kCin2_Cin3 = 0.5 .* kCin2_Cin3;
% kCC_Cin3(7 : end , :) = 2 .* kCC_Cin3(7 : end , :);
muCC = min(muCC .* 12 , 0.99); % convert cervical cancer mortality rate from monthly to yearly
%     fImm(4 : age) = 1; % RR(0.75; 0.5 , 0.92) fraction fully protected by immunity based on RR of natural immunity (Beachler, 2017)
perPartnerHpv = 0.08; % high risk HPV transmission risk per month
rImmuneHiv = 1 ./ hpv_hivClear;
lambdaMultImm(1 : 4) = 1 - 0.01;
lambdaMultImm(5 : 10) = 1 - logspace(log10(0.01) , log10(0.1) , 6);
lambdaMultImm(11 : age) = lambdaMultImm(10);
lambdaMultVax = ones(age , 2);

artHpvMult = hpv_hivMult(1 , 1);
perPartnerHpv = 0.09;
perPartnerHpv_lr = 0.1;
perPartnerHpv_nonV = 0.1;
rImmuneHiv = 2 ./ hpv_hivClear;

save([paramDir , 'calibInitParams']);
%%
kCin1_Inf_Mults = [1 ; 1; 1];
kInf_Cin1_Mults = [1 ; 1 ; 1];
kCC_Cin3_Mults = [1 ; 1 ; 1];
kCin3_Cin2_Mults = [1 ; 1 ; 1];
hpvClearMult = 0.1 .* hpv_hivClear;
%% Continue from last calibration
load([paramDir,'calibInitParams'])
% load([paramDir,'HPV_calib5.dat'])
% 
% % for i = 1 : 3
% % %     kCin1_Inf_Mults(i) = HPV_calib5(i);
% %     kInf_Cin1_Mults(i) = HPV_calib5(3 + i);
% % %     kCC_Cin3_Mults(i) = HPV_calib5(6 + i);
% % end
% 
% rImmuneHiv = HPV_calib5(10 : 13);
% c3c2Mults = HPV_calib5(14 : 17);
% c2c1Mults = max(1 , HPV_calib5(18 : 21));
% artHpvMult = HPV_calib5(22);
% perPartnerHpv= HPV_calib5(23);
% lambdaMultImm = HPV_calib5(24 : 39);
% hpv_hivClear = HPV_calib5(40 : 43);
% hpvClearMult = HPV_calib5(44 : 47);
%%
initParams = [kCin1_Inf_Mults;
    kInf_Cin1_Mults;
    kCC_Cin3_Mults;
    rImmuneHiv; % HPV immunity clearance multiplier for HIV+
    c3c2Mults; % CIN2 to CIN3 progression multiplier for HIV+
    c2c1Mults; % CIN1 to CIN2 progression multiplier for HIV+
    artHpvMult; % HPV acquisition multiplier for HIV+
    perPartnerHpv;
    lambdaMultImm';
    hpv_hivClear;
    hpvClearMult;
    perPartnerHpv_lr;
    perPartnerHpv_nonV;
    kCin3_Cin2_Mults];

lb = initParams .* 0.5;
ub = initParams .* 2; 

% max transition rate between CC to CIN is 0.5
for i = 1 : 3
    ub(i) = min(0.3 / max(kCin1_Inf(: , i)) , 2);
    ub(3 + i) = min(0.8 / max(kInf_Cin1(: , i)) ,  2);
    ub(6 + i) = min(0.5 / max(kCC_Cin3(: , i)) ,  2);
    ub(49 + i) = min(0.6 / max(kCin3_Cin2(: , i)) , 2);
end
lb(14 : 21) = 1; % CIN progression multiplier for HIV+ >=1
lb(22) = 1; % HPV acquisition multiplier for HIV+ on ART >= 1
ub(22) = max(hpv_hivMult(1 , :)); % max HPV acquisition multiplier for HIV+ on ART <= max acquisition multiplier for HIV+
ub(23) = 0.1; % max per partner transmission probability for HPV
ub(40 : 47) = 1; % HPV clearance multipliers <= 1 for HIV-positive
ub(48 : 49) = 0.1; % max per partner transmission probability for HPV

A = zeros(length(initParams));
for i = 0 : 2
    A(i + 1 , 14 + i) = 1;
    A(i + 1 , 15 + i) = -1;
    
    A(i + 4 , 18 + i) = 1;
    A(i + 4 , 19 + i) = -1;
end

b = zeros(length(initParams) , 1);

options = psoptimset('UseParallel' , true , 'Cache' , 'on' ,...
    'CacheTol' , 0.1 , 'CompletePoll' , 'on' , 'TolMesh' , 0.1, ...
    'Display','iter','PlotFcn',@psplotbestf);
x = patternsearch(@calibrator, initParams , A , b , [] , [] , lb , ub , [] , options);
%%
file = 'HPV_calib6.dat';
paramDir = [pwd , '\Params\'];
csvwrite([paramDir, file] , x)
