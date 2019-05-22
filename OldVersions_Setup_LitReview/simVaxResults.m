%%
o1 = load('H:\HHCoM Results\output1.mat'); % 90% coverage
o2 = load('H:\HHCoM Results\output2.mat'); % 50% coverage
o3 = load('H:\HHCoM Results\output3.mat'); % 0% coverage
tVec = o1.tVec;
%% new cervical cancers
newCC_90 = o1.newCC;
newCC_50 = o2.newCC;
newCC_0 = o3.newCC;

%% populations
pop90 = o1.popVec;
pop50 = o2.popVec;
pop0 = o3.popVec;
%% Cumulative Cancer Incidence

% cumulative cases
cum_newCC_90 = cumsum(newCC_90);
cum_newCC_50 = cumsum(newCC_50);
cum_newCC_0 = cumsum(newCC_0);

% cumulative cases averted 
averted90 = cum_newCC_0 - cum_newCC_90;
averted50 = cum_newCC_0 - cum_newCC_50;

% cumulative cases averted across population
genAvert90 = sum(sum(sum(sum(averted90,2),3),4),5);
genAvert50 = sum(sum(sum(sum(averted50,2),3),4),5);

figure()
plot(tVec , genAvert90 , tVec , genAvert50)
legend('90% coverage' , '50% coverage')
title('Cumulative Cases Averted Across General Population')
xlabel('Year'); ylabel('Cases')
%%
% General Incidence
% susceptible population
hpvSusGen = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 : periods , 2 , 4 : age , 1 : risk));
    toInd(allcomb(1 : disease , 1 : viral , 1 , 9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
pop90_susGen = sum(pop90(1 : end - 1 , hpvSusGen) , 2);
pop50_susGen = sum(pop50(1 : end - 1 , hpvSusGen) , 2);
pop0_susGen = sum(pop0(1 : end - 1 , hpvSusGen) , 2);

% cases in general population
genCC90 = sum(sum(sum(sum(newCC_90 , 2),3),4),5);
genCC50 = sum(sum(sum(sum(newCC_50 , 2),3),4),5);
genCC0  = sum(sum(sum(sum(newCC_0 , 2),3),4),5);
%%
newCC_90_15= newCC_90(: ,: ,:  , : , 4 : age);
newCC_50_15= newCC_90(: ,: , : , : , 4 : age);
newCC_0_15= newCC_90(: , : , : , : , 4 : age);

newCC_90_15all = sum(sum(sum(sum(newCC_90_15 , 2) , 3) , 4),5);
newCC_50_15all = sum(sum(sum(sum(newCC_50_15 , 2) , 3) , 4),5);
newCC_0_15all = sum(sum(sum(sum(newCC_0_15 , 2) , 3) , 4),5);
figure()
plot(tVec , newCC_90_15all , tVec , newCC_50_15all , tVec , newCC_0_15all)
legend('90% coverage' , '50% coverage' , 'No vaccination')

%%
% incidence
fac = 10 ^ 5;
i_gen90 = genCC90(2 : end) ./ pop90_susGen .* fac;
i_gen50 = genCC50(2 : end) ./ pop50_susGen .* fac;
i_gen0 = genCC0(2 : end) ./ pop0_susGen .* fac;
figure()
plot(tVec(2 : end) , i_gen90 , tVec(2 : end) , i_gen50 , tVec(2 : end) , i_gen0)
title('Cervical Cancer Incidence')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('90% coverage' , '50% coverage' , 'No vaccination')

% relative incidence reduction
figure()
genRelRed_90 = (i_gen90 - i_gen0) ./ i_gen0 * 100;
genRelRed_50 = (i_gen50 - i_gen0) ./ i_gen0 * 100;
plot(tVec(2 : end) , genRelRed_90 , tVec(2 : end) , genRelRed_50)
title('General Cervical Cancer Reduction')
xlabel('Year'); ylabel('Relative Difference (%)')
legend('90% coverage' , '50% coverage')
%%
plot(tVec(1 : end - 1) , pop90_susGen , tVec(1 : end -1) , pop50_susGen , 