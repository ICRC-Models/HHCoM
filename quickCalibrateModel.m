function quickCalibrateModel()
close all; clear all; clc
load('popData')
load('general');
load('settings');
load('mixInfectIndices')
load('vlAdvancer')
load('fertMat')
load('hivFertMats')
load('deathMat')
load('circMat')
load('vaxer')
load('mixInfectParams');
load('hpvData')
load('popData')
load('HIVParams')
load('hivIndices')
load('hpvIndices')
load('ager')
load('vlBeta')
load('hpvTreatIndices')
load('calibData')
load('calibParams')
load('vaxInds')

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
rNormal_Inf = rNormal_Inf .* 0.85;
% c2c1Mults = 1.5 .* c2c1Mults;
% kCin1_Cin2(1 : end , :) = 1.2 * kCin1_Cin2(1 : end , :);
% kCin2_Cin1(6 : end , :) = 1.25 * kCin2_Cin1(6 : end , :); 
% kCin3_Cin2(10 : end , :) = 1.5 * kCin3_Cin2(10 : end , :);
% kCin2_Cin3 = 0.5 .* kCin2_Cin3;
% kCC_Cin3(7 : end , :) = 2 .* kCC_Cin3(7 : end , :);
muCC = min(muCC .* 12 , 0.99); % convert cervical cancer mortality rate from yearly to monthly
%     fImm(4 : age) = 1; % RR(0.75; 0.5 , 0.92) fraction fully protected by immunity based on RR of natural immunity (Beachler, 2017)
perPartnerHpv = 0.1; % high risk HPV transmission risk per month
rImmuneHiv = 2 ./ hpv_hivClear;
lambdaMultImm(1 : 4) = 1 - 0.01;
lambdaMultImm(5 : 10) = 1 - logspace(log10(0.01) , log10(0.1) , 6);
lambdaMultImm(11 : age) = lambdaMultImm(10);
lambdaMultVax = ones(age , 2);

artHpvMult = hpv_hivMult(1);
perPartnerHpv = 0.08;
rImmuneHiv = 2 ./ hpv_hivClear;

save('calibInitParams');

kCin2_Cin3_Mults = [1 ; 1 ; 1];
kCin3_Cin2_Mults = [1 ; 1 ; 1];
kCC_Cin3_Mults = [1 ; 1 ; 1];



initParams = [kCin2_Cin3_Mults;
    kCin3_Cin2_Mults;
    kCC_Cin3_Mults
    rImmuneHiv;
    c3c2Mults;
    c2c1Mults;
    artHpvMult;
    perPartnerHpv;
    lambdaMultImm'];

lb = initParams .* 0.5;
ub = initParams .* 2; 

% max transition rate between CIN states is 0.5
for i = 1 : 3
    ub(i) = min(2 , 0.5 / max(kCin2_Cin3(: , i)));
    ub(3 + i) = min(0.5 / max(kCin3_Cin2(: , i)) ,  2);
    ub(6 + i) = min(0.5 / max(kCC_Cin3(: , i)) ,  2);
end
lb(22) = 1;
ub(22) = max(hpv_hivMult(1 , :));
ub(23) = 0.1;
options = psoptimset('UseParallel' , true , 'Cache' , 'on' ,...
    'CacheTol' , 0.1 , 'CompletePoll' , 'on' , 'TolMesh' , 0.1, ...
    'Display','iter','PlotFcn',@psplotbestf);
x = patternsearch(@calibrator, initParams , [] , [] , [] , [] , lb , ub , [] , options);
%%
file = 'HPV_calib3.dat';
csvwrite(file , x)
