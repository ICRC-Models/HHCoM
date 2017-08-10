%% tests vlAdv
close all; clear all; clc
disp('vlAdv Tests')
disp(' ')
sumall = @(x) sum(x(:));
risk = 3;
gender = 2;
age = 12;
disease = 10;
viral = 7;
hpvTypes = 5;
hpvStates = 10;
periods = 9;
times = [];
dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age , risk];
zeroPop = zeros(dim);
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
mPop = zeros(age , risk);
fPop = mPop;
for i = 1 : age
    mPop(i , :) = MpopStruc(i, :).* Mpop80(i); %0.70
    fPop(i , :) = FpopStruc(i, :).* Fpop80(i);
end
%% Test 1: Zero population input
tic
dPop1 = vlAdv(zeroPop);
times = [times toc];
disp('Test 1: Pop = 0')
passed = 0;
tests = 1;
newPop1 = dPop1 + zeroPop;
if newPop1 ~= zeroPop
    disp('Failed. Should have no individuals in population')
elseif any(newPop1(:)) < 0
    disp('Failed. Should have no negative values in population.')
else
    disp('Passed')
    passed = passed + 1;
end
disp(' ');

%% Test 2: Healthy pop
pop2 = zeros(dim);
pop2(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop;
pop2(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop;
tic
dPop2 = vlAdv(pop2);
times = [times toc];
newPop2 = dPop2 + pop2;
disp('Test 2: Healthy population')
tests = tests + 1;
if sum(dPop2(:)) ~= 0
    disp('Failed. No changes in population size in this module.')
elseif any(newPop2(:)) < 0
    disp('Failed. Should have no negative values in population.')
elseif ~isequal(newPop2, pop2)
    disp('Failed. Should have no viral load progression in population with no HIV.')
else
    disp('Passed')
    passed = passed + 1;
end
disp(' ')
%% Tests 3 - 72: VL progression when all individuals are placed in a single VL group
i = 3;
for d = 1 : disease
    for v = 1 : viral
        disp(['Test ' , num2str(i) , ': VL progression when individuals are ' ...
            'placed in disease state ' , num2str(d) , ' and VL group ' , num2str(v)]) 
        i = i + 1;
        pop3 = zeros(dim);
        pop3(d , v , 1 , 1 , 1 , 1 , : , :) = mPop;
        pop3(d , v , 1 , 1 , 1 , 2 , : , :) = fPop;
        tic
        dPop3 = vlAdv(pop3);
        times = [times toc];
        newPop3 = dPop3 + pop3;
        tests = tests + 1;
        bool = 1;
        if any(newPop3(:) < -10 ^ -6) 
            disp('Failed. Should have no negative values in population.')
            bool = 0;
        end
        if sum(dPop3(:)) ~= 0
            disp('Failed. Should have no changes in population size in this module.')
            bool = 0;
        end
        for d2 = 1 : disease
            if abs(sumall(pop3(d2 , : , : , : , : , : , : , :))...
                    - sumall(newPop3(d2 , : , : , : , : , : , : , :))) >  10 ^ -6
                disp('Failed. Disease state (CD4) should not change in this module.')
                bool = 0;
            end
        end
        if (d < 2 || d > 6) && abs(sumall(pop3(: , v , 1 , 1 , 1 , : , : , :)) ...
                - sumall(newPop3(: , v , 1 , 1 , 1 , : , : , :))) > 10^ -6
            disp(['Failed. Should have no progression in viral load when d = ' , num2str(d)])
            bool = 0;
        end
        if d >= 2 && d <= 6 && abs(sumall(pop3(: , v , 1 , 1 , 1 , : , : , :)) ...
                <= sumall(newPop3(: , v , 1 , 1 , 1 , : , : , :)))...
                && v < 6 && v > 1 ... % individuals do not progress in viral load beyond VL >50,000 copies/mL 
            disp(['Failed. Should have progression in viral load when' ...
            ' individuals are infected and are not yet at VL > 50,000 copies/mL.'])
            bool = 0;
        elseif abs(sumall(pop3(: , v , 1 , 1 , 1 , : , : , :)) ...
                - sumall(newPop3(: , v , 1 , 1 , 1 , : , : , :))) > 10 ^ -6 ...
                && v == 6
            disp(['Failed. Viral load levels of individuals with viral load beyond ',...
                'VL >50,000 copies/mL should not be progressing any further. They are already at ',...
                'the highest viral load level.'])
            bool = 0;
        elseif v < 7 && sumall(pop3(: , v + 1 , 1 , 1 , 1 , : , : , :)) ...
                > sumall(newPop3(: , v + 1 , 1 , 1 , 1 , : , : , :))
            disp(['Failed. Viral load group ' , num2str(v + 1) , ' should have more individuals than viral load group ', ...
                num2str(v) , ' after population goes through vlAdv.'])
            bool = 0;
        elseif v > 1 && sumall(pop3(: , v - 1 , 1 , 1 , 1 , : , : , :)) ~= 0 
            disp('Failed. No individuals should be regressing into lower viral load groups.')
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
