close all; clear all; clc
disp('Optimizing...')
disp(' ')
nvars = 54;
lb = ones(nvars , 1) .* 0.7;
ub = ones(nvars , 1) .* 3;
x0 = ones(nvars , 1);
startYear = 1980;
endYear = 2012;
options = optimset('PlotFcns' , @optimplotfval);
[x,fval,maxfval,exitflag,output]  = ...
    fminimax(@modelFit , x0  , [] , [] , [] , [] , lb , ub , [] , options);
disp(' ')
disp('Optimization complete')