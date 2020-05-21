function buildModel()
close all; clear all; clc
stepsPerYear = 6;
loadUp(stepsPerYear)
disp('This could take a while...')
disp('Data loaded')
makeMat
disp('Transition matrices constructed')
disp(' ')
disp('Running main.m')
main()