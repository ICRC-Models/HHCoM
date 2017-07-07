close all; clear all; clc
disp('Calibrating to overall HIV prevalence from 1990 to 2012...')
% epsA_Initial = [0.7 , 0.3 , 0.1];
% epsR_Initial = [0.8 , 0.3 , 0.1];
% eps = [epsA_Initial ; epsR_Initial]';
% lb = 0.05 * ones(size(eps));
% ub = 0.99 * ones(size(eps));
gender = 2;
ages = 12;
risk = 3;
cMult = ones(5 , risk - 1 , gender);
lb = 0.05 * ones(size(cMult));
ub = 2.5 * ones(size(cMult));
% if isempty(gcp)
%     parpool local
% end
options = optimoptions('patternsearch', 'UseParallel' , true , 'cache' , 'on' , 'Display','iter','PlotFcn',@psplotbestf);
x = patternsearch(@prevCalib, cMult , [] , [] , [] , [] , lb , ub , [] , options);
%%
% disp('Optimized values for epsA:') 
% disp(num2str(x(: , 1)))
% disp(' ')
% disp('Optimized values for epsR:')
% disp(num2str(x(: , 2)))
% disp(' ')
% file = 'Prevalence calibration.dat';
disp('Optimized values for partnersM:') 
disp(num2str(x(: , : , 1)))
disp(' ')
disp('Optimized values for partnersF:')
disp(num2str(x(: , : , 2)))
disp(' ')
file = 'Prevalence calibration by partners.dat';
csvwrite(file , x);
disp(['Results saved to ' file]);