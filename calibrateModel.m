function calibrateModel()
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

artHpvMult = hpv_hivMult * 0.5;
perPartnerHpv = 0.08;
rImmuneHiv = 1 ./ hpv_hivClear;

initParams = [kCin2_Cin3(: , 1);
kCin3_Cin2(: , 1);
kCC_Cin3(: , 1);
kCin2_Cin3(: , 2);
kCin3_Cin2(: , 2);
kCC_Cin3(: , 2);
kCin2_Cin3(: , 3);
kCin3_Cin2(: , 3);
kCC_Cin3(: , 3);
rImmuneHiv;
c3c2Mults;
c2c1Mults;
artHpvMult;
perPartnerHpv];

lb = initParams .* 0.4;
lb(end) = 0.01;
ub = initParams .* 3; 
ub(end) = 0.1;


options = optimoptions('patternsearch', 'UseParallel' , true , 'cache' , 'on');
x = patternsearch(@calibrator, initParams , [] , [] , [] , [] , lb , ub , [] , options);
%%
file = 'HPV_calib.dat';
csvwrite(file , x)
