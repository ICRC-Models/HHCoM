%% tests hiv
close all; clear all; clc
disp('hiv Tests')
disp(' ')
sumall = @(x) sum(x(:));
load('general')
times = [];
dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age , risk];
Mpop80 = [293191; %Nyanza pop scaled by Kenya
        221895;
        184334;
        149528;
        114647;
        89629;
        68154;
        54083;
        45403;
        37906;
        30905;
        25188].*0.88;
Msum80 = sum(Mpop80);

Fpop80 = [287674;
        218405;
        182359;
        151279;
        120244;
        95951;
        72704;
        57338;
        47338;
        38974;
        31547;
        25279].*0.88;
Fsum80 = sum(Fpop80);


MpopStruc = [.999 .0005 0.0005; %Extrapolated from RB
            .999 .0005 0.0005;
            .98 .015 .005;
            .80 .17 .03;
            .78 .20 .02;
            .65 .29 .06;
            .66 .28 .06;
            .68 .27 .05;
            .75 .20 .05;
            .78 .17 .05;
            .88 .08 .04;
            .96 .035 .005];

FpopStruc = [.998 .001 0.001; % --calibrated %Extrapolated from RB
            .998 .001 0.001;
            .975 .015 .01;
            .80 .17 .03;
            .62 .31 .05;
            .60 .35 .05;
            .65 .30 .05;
            .65 .30 .05;
            .78 .17 .05;
            .80 .16 .04;
            .85 .13 .02;
            .95 .045 .005];
%
mPop = zeros(age , risk);
fPop = mPop;
for i = 1 : age
    mPop(i , :) = MpopStruc(i, :).* Mpop80(i); %0.70
    fPop(i , :) = FpopStruc(i, :).* Fpop80(i);
end
%% Test 1: Zero population input
zeroPop = zeros(dim);
tic
artDist = ones(disease , viral , age , risk);
[dPop1 , ~] = hiv(zeroPop , artDist);
times = [times toc];
disp('Test 1: Pop = 0')
passed = 0;
tests = 1;
if dPop1 + zeroPop ~= zeroPop
    disp('Failed. Should have no individuals in population')
else
    disp('Passed')
    passed = passed + 1;
end
disp(' ')
%% Test 2: Healthy pop
pop2 = zeros(dim);
pop2(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop;
pop2(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop;
tic
[dPop2 , ~] = hiv(pop2 , artDist);
times = [times toc];
newPop2 = dPop2 + pop2;
disp('Test 2: Healthy population')
passed = 0;
tests = tests + 1;
if sum(dPop2(:)) ~= 0
    disp('Failed. No changes in population size when all individuals are healthy.')
elseif isequal(newPop2, pop2)
    disp('Failed. Some healthy individuals should go on PrEP.')
else
    disp('Passed')
    passed = passed + 1;
end
disp(' ')
%% Tests 3 - 72: CD4 progression when all individuals are placed in a single HIV disease stage compartment
i = 3;
for d = 1 : disease
    for v = 1 : viral
        disp(['Test ' , num2str(i) , ': CD4 progression when individuals are ' ...
            'placed in disease state ' , num2str(d) , ' and VL group ' , num2str(v)]) 
        i = i + 1;
        pop3 = zeros(dim);
        pop3(d , v , 1 , 1 , 1 , 1 , : , :) = mPop;
        pop3(d , v , 1 , 1 , 1 , 2 , : , :) = fPop;
        tic
        [dPop3 , ~] = hiv(pop3 , artDist);
        times = [times toc];
        newPop3 = dPop3 + pop3;
        tests = tests + 1;
        bool = 1;
        if any(newPop3(:) < -10 ^ -6) 
            disp('Failed. Should have no negative values in population.')
            bool = 0;
        end
        if d > 1 && d < 7 && sum(newPop3(:)) >= sum(pop3(:))
            disp('Failed. Some individuals should die due to disease associated mortality.')
            bool = 0;
        end
        if abs(sumall(pop3(d , : , 1 , 1 , 1 , : , : , :)) ...
                <= sumall(newPop3(d , : , 1 , 1 , 1 , : , : , :)))...
                && d > 1 && d < 7
            disp(['Failed. Should have progression in CD4 when' ...
            ' individuals are infected and not yet at max progression.'])
            bool = 0;
        elseif d > 1 && d < 7 && sumall(pop3(d , : , 1 , 1 , 1 , : , : , :)) ...
                < sumall(newPop3(d , : , 1 , 1 , 1 , : , : , :))
            disp(['Failed. New population should have fewer individuals in disease state ' ,...
                num2str(d) , ' after individuals in d = ' , num2str(d) , ...
                ' go through HIV module.'])
            bool = 0;
        elseif d > 1 && d < 6 && sumall(newPop3(d + 1 , : , 1 , 1 , 1 , : , : , :)) ...
                < sumall(pop3(d + 1, : , 1 , 1 , 1 , : , : , :))
            disp(['Failed. New population should have more individuals in disease state ' , ...
                num2str(d + 1), ' after individuals in d = ' , num2str(d) , ' go through ', ...
                'HIV module.'])
            bool = 0;
        elseif d > 1 && sumall(pop3(d - 1 , : , : , : , : , : , : , :)) ~= 0 
            disp('Failed. No individuals should be regressing into lower disease stage groups.')
            bool = 0;
        end
        if bool == 1
            disp('Passed')
            passed = passed + 1;
        end
        disp(' ')
    end
end
%% Summary
tMean = mean(times);
t240 = 240 * tMean;
disp(' ');
disp(['Passed ' , num2str(passed) , ' out of ' , num2str(tests) , ' tests.'])
disp(['Average runtime was ' , num2str(tMean) , ' s per iteration. (Max = ' , ...
    num2str(max(times)) , ' s, Min = ' , num2str(min(times)) , ' s, StDev = ',...
    num2str(std(times)), ' s)'])
disp(['Estimated runtime over 240 iterations = ' , num2str(t240) , ' s.'])
