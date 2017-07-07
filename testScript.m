% Test Script
%% bornDie Test
close all; clear all; clc
risk = 3;
gender = 2;
age = 12;
disease = 10;
viral = 7;
hpvTypes = 5;
hpvStates = 7;
times = [];
zeroPop = zeros(disease , viral , hpvTypes , hpvStates , gender , age , risk);
dim = [disease , viral , hpvTypes , hpvStates , gender , age , risk];
tic
[dPop1 , ~ , ~] = bornDie(zeroPop);
times = [times toc];
disp('Test 1: Pop = 0')
passed = 0;
tests = 1;
if dPop1 ~= 0
    disp('Failed. Should have no individuals in population')
elseif dPop1(: , : , 5 , 6 , 2 , : , :) ~= 0 % check if vaxd
    disp('Failed. No individuals should be detected in vaccinated subpopulation with 0 population input.')
else
    disp('Passed')
    passed = passed + 1;
end

%% no adults
close all;
tests = tests + 1;
inPop2 = zeros(dim);
inPop2(: , : , : , : , : , 1 : 2 , :) = 1000;
tic
[dPop2 , bS2 , bI2] = bornDie(inPop2);
times = [times toc];
disp('Test 2: No reproducing individuals')
if bS2 ~= 0
    disp('Failed. bS != 0.')
elseif bI2 ~= 0
    disp('Failed. bI != 0.')
elseif dPop2(: , : , 5 , 6 , 2 , : , :) ~= 0 % check if vaxd
    disp('Failed. No individuals should be detected in vaccinated subpopulation. Eligible age not yet reached.')
else
    disp('Passed')
    passed = passed + 1;
end
%% reproducing adults
% females only
close all; 
tests = tests + 1;
inPop3 = zeros(dim);
inPop3(: , : , : , : , 2 , 3 : 9 , :) = 1000;
tic
[dPop3 , bS3 , bI3] = bornDie(inPop3);
times = [times toc];
disp('Test 3: Some reproducing females')
if bS3 <= 0
    disp('Failed. bS <= 0.')
elseif bI3 <= 0
    disp('Failed. bI <= 0.')
elseif dPop3(: , : , 5 , 6 , 2 , : , :) == 0 % check if vaxd
    disp('Failed. No individuals detected in vaccinated subpopulation.')
else
    disp('Passed')
    passed = passed + 1;
end

% females and males
tests = tests + 1;
inPop4 = zeros(dim);
inPop4(: , : , : , : , : , 3 : 9 , :) = 1000;
tic
[dPop4 , bS4 , bI4] = bornDie(inPop4);
times = [times toc];
disp('Test 4: Same number of reproducing females as in Test 3 with males added.')
if bS4 <= 0
    disp('Failed. bS <= 0.')
elseif bI4 <= 0
    disp('Failed. bI <= 0.')
elseif bI4 ~= bI3
    disp('Failed. bI with males and females != bI with females alone.')
elseif bS4 ~= bS3
    disp('Failed. bS with males and females != bS with females alone.')
elseif dPop4(: , : , 5 , 6 , 2 , : , :) == 0 % check if vaxd
    disp('Failed. No individuals detected in vaccinated subpopulation.')
else
    disp('Passed')
    passed = passed + 1;
end

% males only
close all;
tests = tests + 1;
inPop5 = zeros(dim);
inPop5(: , : , : , : , 1 , 3 : 9 , :) = 1000;
tic
[dPop5 , bS5 , bI5] = bornDie(inPop5);
times = [times toc];
disp('Test 5: Input only males of reproducing age')
if bS5 > 0 
    disp('Failed. bS > 0. Reproduction not dependent on males in model.')
elseif bI5 > 0
    disp('Failed. bI > 0. Reproduction not dependent on males in model.')
elseif dPop5(: , : , 5 , 6 , 2 , : , :) ~= 0 % check if vaxd
    disp('Failed. No males should be detected in vaccinated subpopulation. Model shoud only have female vaccination for now.')
else
    disp('Passed')
    passed = passed + 1;
end

%% Summary
tMean = mean(times);
t200 = 200 * tMean;
disp(['Passed ' , num2str(passed) , ' out of ' , num2str(tests) , ' tests.'])
disp(['Average runtime was ' , num2str(tMean) , ' s per iteration. (Max = ' , ...
    num2str(max(times)) , ' s, Min = ' , num2str(min(times)) , ' s)'])
disp(['Estimated runtime over 200 iterations = ' , num2str(t200) , ' s.'])


