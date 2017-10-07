close all; clear all; clc
%% parameters
endYear = 2100;
testParams = [0.9 , 0.5 , 0];
nTests = length(testParams);
calibPop = load('H:\HHCoM Results\results');
currPop = calibPop.popLast;
simVax(endYear , nTests , testParams , currPop)