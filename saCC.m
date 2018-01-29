load('general')
% base CIN3 -> CC
o90 = load('H:\HHCoM_Results\Vax_0.9_wane_0.mat'); % 90 coverage
o70 = load('H:\HHCoM_Results\Vax_0.7_wane_0.mat'); % 70 coverage
o50 = load('H:\HHCoM_Results\Vax_0.5_wane_0.mat'); % 50 coverage
oNo = load('H:\HHCoM_Results\Vax_0_wane_0.mat'); % No vaccine

% 0.5x CIN3 -> CC
o90_0_5 = load('H:\HHCoM_Results\0.5_CIN3_CCVax_0.9_wane_0.mat'); %90 coverage
o70_0_5 = load('H:\HHCoM_Results\0.5_CIN3_CCVax_0.7_wane_0.mat'); %70 coverage
o50_0_5 = load('H:\HHCoM_Results\0.5_CIN3_CCVax_0.5_wane_0.mat'); %50 coverage
oNo_0_5 = load('H:\HHCoM_Results\0.5_CIN3_CCVax_0_wane_0.mat'); % No vaccine

%1.5x CIN3 -> CC
o90_1_5 = load('H:\HHCoM_Results\1.5_CIN3_CCVax_0.9_wane_0.mat'); %90 coverage
o70_1_5 = load('H:\HHCoM_Results\1.5_CIN3_CCVax_0.7_wane_0.mat'); %70 coverage
o50_1_5 = load('H:\HHCoM_Results\1.5_CIN3_CCVax_0.5_wane_0.mat'); %50 coverage
oNo_1_5 = load('H:\HHCoM_Results\1.5_CIN3_CCVax_0_wane_0.mat'); % No vaccine

% same time vector for all scenarios
tVec = o90.tVec;

% vector of weights for age standardization
wVec = zeros(age , 1);
wVec(5 : age) = [0.188 , 0.18 , 0.159 , 0.121 , 0.088 , 0.067 , 0.054 , ...
    0.046 , 0.038 , 0.029 , 0.017 , 0.013]; 

%% plot settings
set(0 , 'defaultlinelinewidth' , 2)
set(0,'defaultAxesFontSize',16)

%% populations
% base CIN3 -> CC
pop90Arr{1} = o90.popVec;
pop70Arr{1} = o70.popVec;
pop50Arr{1} = o50.popVec;
popNoArr{1} = oNo.popVec;

% 0.5x CIN3 -> CC
pop90Arr{2} = o90_0_5.popVec;
pop70Arr{2} = o70_0_5.popVec;
pop50Arr{2} = o50_0_5.popVec;
popNoArr{2} = oNo_0_5.popVec;

% 1.5x CIN3 -> CC
pop90Arr{3} = o90_1_5.popVec;
pop70Arr{3} = o70_1_5.popVec;
pop50Arr{3} = o50_1_5.popVec;
popNoArr{3} = oNo_1_5.popVec;

%% New cervical cancers
% base CIN3 -> CC
newCC_90Arr{1} = o90.newCC;
newCC_70Arr{1} = o70.newCC;
newCC_50Arr{1} = o50.newCC;
newCC_0Arr{1} = oNo.newCC;

% 0.5x CIN3 -> CC
newCC_90Arr{2} = o90_0_5.newCC;
newCC_70Arr{2} = o70_0_5.newCC;
newCC_50Arr{2} = o50_0_5.newCC;
newCC_0Arr{2} = oNo_0_5.newCC;

% 1.5x CIN3 -> CC
newCC_90Arr{3} = o90_1_5.newCC;
newCC_70Arr{3} = o70_1_5.newCC;
newCC_50Arr{3} = o50_1_5.newCC;
newCC_0Arr{3} = oNo_1_5.newCC;

%% SA of general cervical cancer incidence reduction (90/70/50 coverage, lifelong efficacy)

% General Incidence
% susceptible population
susGen = zeros(age , length(tVec) - 1);
pop90_susGen = {susGen , susGen , susGen};
pop70_susGen = pop90_susGen;
pop50_susGen = pop90_susGen;
popNo_susGen = pop90_susGen;
for i = 1 : size(pop90Arr , 2)
    pop90 = pop90Arr{i};
    pop70 = pop70Arr{i};
    pop50 = pop50Arr{i};
    popNo = popNoArr{i};
    for a = 1 : age
        ccSusGen = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , 2 , a , 1 : risk));
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , 2 , a , 1 : risk))];
        pop90_susGen{i}(a , :) = (sum(pop90(1 : end - 1 , ccSusGen) , 2) + sum(pop90(2 : end , ccSusGen) , 2)) * 0.5;
        pop70_susGen{i}(a , :) = (sum(pop70(1 : end - 1 , ccSusGen) , 2) + sum(pop70(2 : end , ccSusGen) , 2)) * 0.5;
        pop50_susGen{i}(a , :) = (sum(pop50(1 : end - 1 , ccSusGen) , 2) + sum(pop50(2 : end , ccSusGen) , 2)) * 0.5;
        popNo_susGen{i}(a , :) = (sum(popNo(1 : end - 1 , ccSusGen) , 2) + sum(popNo(2 : end , ccSusGen) , 2)) * 0.5;
    end
end

%cases in general population
ageCase = zeros(age , length(tVec));
genCC90Age = {ageCase , ageCase , ageCase};
genCC70Age = genCC90Age;
genCC50Age = genCC90Age;
genCCNoAge = genCC90Age;
for i = 1 : size(pop50_susGen , 2)
    genCC90Age{i} = squeeze(sum(sum(sum(newCC_90Arr{i} , 2),3),4))';
    genCC70Age{i} = squeeze(sum(sum(sum(newCC_70Arr{i} , 2),3),4))';
    genCC50Age{i}  = squeeze(sum(sum(sum(newCC_50Arr{i} , 2),3),4))';
    genCCNoAge{i} = squeeze(sum(sum(sum(newCC_0Arr{i} , 2),3),4))';
end
%%
% incidence
fac = 10 ^ 5;
i_gen90 = zeros(3 , length(tVec) - 1);
i_gen70 = i_gen90;
i_gen50 = i_gen90;
i_genNo = i_gen90;
genRelRed_90 = i_gen90;
genRelRed_70 = i_gen90;
genRelRed_50 = i_gen90;
figure()
for i = 1 : 3
    i_gen90(i , :) = sum(bsxfun(@times , squeeze(genCC90Age{i}(: , 2 : end)) ./ pop90_susGen{i} , wVec)) .* fac;
    i_gen70(i , :) = sum(bsxfun(@times , squeeze(genCC70Age{i}(: , 2 : end)) ./ pop70_susGen{i} , wVec)) .* fac;
    i_gen50(i , :) = sum(bsxfun(@times , squeeze(genCC50Age{i}(: , 2 : end)) ./ pop50_susGen{i} , wVec)) .* fac;
    i_genNo(i , :) = sum(bsxfun(@times , squeeze(genCCNoAge{i}(: , 2 : end)) ./ popNo_susGen{i} , wVec)) .* fac;
    % figure()
    % plot(tVec(2 : end) , i_gen90 , tVec(2 : end) , i_gen70 , tVec(2 : end) , i_gen50 , tVec(2 : end) , i_genNo)
    % title('General Cervical Cancer Incidence')
    % xlabel('Year'); ylabel('Incidence per 100,000')
    % legend('90% coverage' , '70% coverage' , '50% coverage' , 'No vaccination')

    % relative incidence reduction
%     figure()
    genRelRed_90(i , :) = (i_gen90(i , :) - i_genNo(i , :)) ./ i_genNo(i , :) * 100;
    genRelRed_70(i , :) = (i_gen70(i , :) - i_genNo(i , :)) ./ i_genNo(i , :) * 100;
    genRelRed_50(i , :) = (i_gen50(i , :) - i_genNo(i , :)) ./ i_genNo(i , :) * 100;
end
%%
% Sensitivity analysis plot of CC incidence reduction (90/70/50 coverage , 
% 0.5, 1, 1.5x CIN3 -> CC progression rate)
p1 = plot(tVec(2 : end) , genRelRed_90(1 , :) , 'b' , tVec(2 : end) , genRelRed_90(2 , :)...
   , 'b--' , tVec(2 : end) , genRelRed_90(3 , :) , 'b:');
hold on
p2 = plot(tVec(2 : end) , genRelRed_70(1 , :) , 'r' , tVec(2 : end) , genRelRed_70(2 , :)...
    , 'r--' , tVec(2 : end) , genRelRed_70(3 , :) , 'r:');
hold on
p3 = plot(tVec(2 : end) , genRelRed_50(1 , :) , 'm' , tVec(2 : end) , genRelRed_50(2 , :)...
    , 'm--' , tVec(2 : end) , genRelRed_50(3 , :) , 'm:');
hold off
legend('90% coverage' , '90% (0.5x CC)' , '90% (1.5x CC)',...
'70% coverage' , '70% (0.5x CC)' , '70% (1.5x CC)' , ...
'50% coverage' , '50% (0.5x CC)' , '50% (1.5x CC)')
title('General Cervical Cancer Incidence Reduction')
axis([tVec(1) , tVec(end) , -100 , 0])
xlabel('Year'); ylabel('Relative Difference (%)')

%% Export general incidence and reduction as csv

%% output reductions for general population
yr = 2088;
t = (2088 - 2018) * stepsPerYear;
Base = [genRelRed_90(1 , t) ; genRelRed_70(1 , t) ; genRelRed_50(1, t)];
LowCC = [genRelRed_90(2 , t) ; genRelRed_70(2 , t) ; genRelRed_50(2, t)];
HighCC = [genRelRed_90(3 , t) ; genRelRed_70(3 , t) ; genRelRed_50(3, t)];
Coverage = {'90%' ; '70%' ; '50%'};
T_Incidence = table(Coverage , Base , LowCC , HighCC);

% T = table(tVec(2 : end)' , i_gen90' , i_gen70' , i_gen50' , i_genNo' , genRelRed_90' , ...
%     genRelRed_70' , genRelRed_50');
% writetable(T , 'General_Incidence_stand.csv' , 'Delimiter' , ',')


% SA of general cervical cancer mortality reduction (90/70/50 coverage, lifelong efficacy