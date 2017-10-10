%%
load('general')
o90 = load('H:\HHCoM_Results\Vax_0.9_wane_0.mat'); % 90% coverage
o70 = load('H:\HHCoM_Results\Vax_0.7_wane_0.mat'); % 70% coverage
o50 = load('H:\HHCoM_Results\Vax_0.5_wane_0.mat'); % 50% coverage
oNo = load('H:\HHCoM_Results\Vax_0_wane_0.mat'); % No vaccine
tVec = o90.tVec;
set(0 , 'defaultlinelinewidth' , 2)

%% CC associated deaths
% general
allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 2 , 4 : age , 1 : risk));
% All HIV-positive women (not on ART)
allHivF = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 2 , 4 : age , 1 : risk));
% All HIV-negative women
hivNeg = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    2 , 4 : age , 1 : risk));
inds = {':' , 2 : 6 , 1};
genArray = {allF , allHivF , hivNeg};
files = {'General_CCMortality_VaxCover' , 'HivAll_CCMortality_VaxCover' , 'HivNeg_CCMortality_VaxCover'};
plotTits = {'General Cervical Cancer' , 'HIV+ Cervical Cancer' , 'HIV- Cervical Cancer'};
for i = 1 : length(genArray)
    vNo_Mort = ...
        sum(sum(sum(sum(oNo.ccDeath(2 : end , inds{i} , : , : , :),2),3),4),5) ./ ...
        sum((oNo.popVec(1 : end - 1 , genArray{i}) ...
        + oNo.popVec(2 : end , genArray{i})) * 0.5 , 2) * fac;

    v90_Mort = ...
        sum(sum(sum(sum(o90.ccDeath(2 : end , inds{i} , : , : , :),2),3),4),5) ./ ...
        sum((o90.popVec(1 : end - 1 , genArray{i}) ...
        + o90.popVec(2 : end , genArray{i})) * 0.5 , 2) * fac;

    v70_Mort = ...
        sum(sum(sum(sum(o70.ccDeath(2 : end , inds{i} , : , : , :),2),3),4),5) ./ ...
        sum((o70.popVec(1 : end - 1 , genArray{i}) ...
        + o70.popVec(2 : end , genArray{i})) * 0.5 , 2) * fac;

    v50_Mort = ...
        sum(sum(sum(sum(o50.ccDeath(2 : end , inds{i} , : , : , :),2),3),4),5) ./ ...
        sum((o50.popVec(1 : end - 1 , genArray{i}) ...
        + o50.popVec(2 : end , genArray{i})) * 0.5 , 2) * fac;
    figure()
    plot(tVec(2 : end) , vNo_Mort , tVec(2 : end) , v90_Mort , ...
        tVec(2 : end) , v70_Mort , tVec(2 : end) , v50_Mort)
    title([plotTits{i} , ' Mortality'])
    xlabel('Year'); ylabel('Mortality per 100,000')
    legend('No vaccination' , '90% coverage' , '70% coverage' , ...
        '50% coverage')
    % Reduction
    v90_Red = (v90_Mort - vNo_Mort) ./ vNo_Mort * 100;
    v70_Red = (v70_Mort - vNo_Mort) ./ vNo_Mort * 100;
    v50_Red = (v50_Mort - vNo_Mort) ./ vNo_Mort * 100;

    figure()
    plot(tVec(2 : end) , v90_Red , tVec(2 : end) , v70_Red , ...
        tVec(2 : end) , v50_Red)
    title([plotTits{i} , ' Mortality Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('90% coverage' , '70% coverage' , '50% coverage')
    axis([tVec(2) tVec(end) -100 0])

    T = table(tVec(2 : end)' , v90_Mort , v70_Mort , v50_Mort , ...
        v90_Red , v70_Red , v50_Red);
    writetable(T , [files{i} , '.csv'] , 'Delimiter' , ',')
end

%% By CD4
dVec = [2 : 6 , 10];
tits = {'Acute' , 'CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , 'CD4 < 200' , ...
    'ART'};
filenames = {'AcuteMort' , 'CD4_500Mort' , 'CD4_500_350Mort' , 'CD4_350_200Mort' , ...
    'CD4_200Mort' , 'ARTMort'};
for d = 1 : length(dVec)
    hiv_ccSus = [toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , 4 : age , 1 : risk)) ; toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];

    vNo_wNoMort = ...
        sum(sum(sum(sum(oNo.ccDeath(2 : end , dVec(d) , : , : , :),2),3),4),5) ./ ...
        sum((oNo.popVec(1 : end - 1 , hiv_ccSus) ...
        + oNo.popVec(2 : end , hiv_ccSus)) * 0.5 , 2) * fac;

    v90_Mort = ...
        sum(sum(sum(sum(o90.ccDeath(2 : end , dVec(d) , : , : , :),2),3),4),5) ./ ...
        sum((o90.popVec(1 : end - 1 , hiv_ccSus) ...
        + o90.popVec(2 : end , hiv_ccSus)) * 0.5 , 2) * fac;

    v70_Mort = ...
        sum(sum(sum(sum(o70.ccDeath(2 : end , dVec(d) , : , : , :),2),3),4),5) ./ ...
        sum((o70.popVec(1 : end - 1 , hiv_ccSus) ...
        + o70.popVec(2 : end , hiv_ccSus)) * 0.5 , 2) * fac;

    v50_Mort = ...
        sum(sum(sum(sum(o50.ccDeath(2 : end , dVec(d) , : , : , :),2),3),4),5) ./ ...
        sum((o50.popVec(1 : end - 1 , hiv_ccSus) ...
        + o50.popVec(2 : end , hiv_ccSus)) * 0.5 , 2) * fac;

    figure(104)
    subplot(3 , 2 , d)
    plot(tVec(2 : end) , vNo_wNoMort , tVec(2 : end) , v90_Mort , ...
        tVec(2 : end) , v70_Mort , tVec(2 : end) , v50_Mort)
    title([tits{i} , ' Mortality'])
    xlabel('Year'); ylabel('Mortality per 100,000')
    legend('No vaccination' , '90% coverage' , '70% coverage' ,...
        '50% coverage' , 'Location' , 'northeastoutside')
    % Reduction
    v90_Red = (v90_Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v70_Red = (v70_Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v50_Red = (v50_Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;

    figure(105)
    subplot(3 , 2 ,d)
    plot(tVec(2 : end) , v90_Red , tVec(2 : end) , v70_Red , ...
        tVec(2 : end) , v50_Red)
    title([tits{d} , ' Mortality Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('90% coverage' , '70% coverage' , '50% coverage', ...
        'Location' , 'northeastoutside')
    axis([tVec(2) tVec(end) -100 0])

    T = table(tVec(2 : end)' , v90_Mort , v70_Mort , v50_Mort , ...
        v90_Red , v70_Red , v50_Red);
    writetable(T , ['VaxCover_' , filenames{d} , '.csv'] , 'Delimiter' , ',')
end
%% new cervical cancers
newCC_90 = o90.newCC;
newCC_70 = o70.newCC;
newCC_50 = o50.newCC;
newCC_0 = oNo.newCC;

%% populations
pop90 = o90.popVec;
pop70 = o70.popVec;
pop50 = o50.popVec;
popNo = oNo.popVec;
%%
% General Incidence
% susceptible population
ccSusGen = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , 2 , 4 : age , 1 : risk));
    toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 8 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
pop90_susGen = sum(pop90(1 : end - 1 , ccSusGen) , 2);
pop70_susGen = sum(pop70(1 : end - 1 , ccSusGen) , 2);
pop50_susGen = sum(pop50(1 : end - 1 , ccSusGen) , 2);
popNo_susGen = sum(popNo(1 : end - 1 , ccSusGen) , 2);

% cases in general population
genCC90 = sum(sum(sum(sum(newCC_90 , 2),3),4),5);
genCC70 = sum(sum(sum(sum(newCC_70 , 2),3),4),5);
genCC50  = sum(sum(sum(sum(newCC_50 , 2),3),4),5);
genCCNo = sum(sum(sum(sum(newCC_0 , 2),3),4),5);
%%
% incidence
fac = 10 ^ 5;
i_gen90 = genCC90(2 : end) ./ pop90_susGen .* fac;
i_gen70 = genCC70(2 : end) ./ pop70_susGen .* fac;
i_gen50 = genCC50(2 : end) ./ pop50_susGen .* fac;
i_genNo = genCCNo(2 : end) ./ popNo_susGen .* fac;
figure()
plot(tVec(2 : end) , i_gen90 , tVec(2 : end) , i_gen70 , tVec(2 : end) , i_gen50 , tVec(2 : end) , i_genNo)
title('General Cervical Cancer Incidence')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('90% coverage' , '70% coverage' , '50% coverage' , 'No vaccination')

% relative incidence reduction
figure()
genRelRed_90 = (i_gen90 - i_genNo) ./ i_genNo * 100;
genRelRed_70 = (i_gen70 - i_genNo) ./ i_genNo * 100;
genRelRed_50 = (i_gen50 - i_genNo) ./ i_genNo * 100;
plot(tVec(2 : end) , genRelRed_90 , tVec(2 : end) , genRelRed_70 , tVec(2 : end) , genRelRed_50)
title('General Cervical Cancer Reduction')
axis([tVec(1) , tVec(end) , -100 , 0])
xlabel('Year'); ylabel('Relative Difference (%)')
legend('90% coverage' , '70% coverage' , '50% coverage')
%% Export general incidence and reduction as csv
T = table(tVec(2 : end)' , i_gen90 , i_gen70 , i_gen50 , i_genNo , genRelRed_90 , ...
    genRelRed_70 , genRelRed_50);
writetable(T , 'General_Incidence.csv' , 'Delimiter' , ',')

%% HIV specific incidence
%% By CD4
% incidence
fac = 10 ^ 5;
figure()
tits = {'Acute' , 'CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , 'CD4 < 200' , ...
    'ART'};
filenames = {'Acute' , 'CD4_500' , 'CD4_500_350' , 'CD4_350_200' , 'CD4_200' , ...
    'ART'};
hiv_vec = [2 : 6 , 10];
for i = 1 : length(hiv_vec)
    d = hiv_vec(i);
    hivCC90 = sum(sum(sum(newCC_90(: , d , : , : , :),3),4),5);
    hivCC70 = sum(sum(sum(newCC_70(: , d , : , : , :),3),4),5);
    hivCC50 = sum(sum(sum(newCC_50(: , d , : , : , :),3),4),5);
    hivCCNo = sum(sum(sum(newCC_0(: , d , : , : , :),3),4),5);

    hivSus = [toInd(allcomb(d , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , 4 : age , 1 : risk)) ; toInd(allcomb(d , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
    pop90_susHiv = sum(pop90(1 : end - 1 , hivSus) , 2);
    pop70_susHiv = sum(pop70(1 : end - 1 , hivSus) , 2);
    pop50_susHiv = sum(pop50(1 : end - 1 , hivSus) , 2);
    popNo_susHiv = sum(popNo(1 : end - 1 , hivSus) , 2);


    h_gen90 = hivCC90(2 : end) ./ pop90_susHiv .* fac;
    h_gen70 = hivCC70(2 : end) ./ pop70_susHiv .* fac;
    h_gen50 = hivCC50(2 : end) ./ pop50_susHiv .* fac;
    h_genNo = hivCCNo(2 : end) ./ popNo_susHiv .* fac;
    subplot(3 , 2 , i)
    plot(tVec(2 : end) , h_gen90 , tVec(2 : end) , h_gen70 , tVec(2 : end) , ...
        h_gen50 , tVec(2 : end) , h_genNo)
    title(['Cervical Cancer Incidence ', tits{i}])
    xlabel('Year'); ylabel('Incidence per 100,000')

    % Export HIV-positive incidence as csv
    T = table(tVec(2 : end)' , h_gen90 , h_gen70 , h_gen50 , h_genNo);
    writetable(T , [filenames{i} , '_Incidence.csv'] , 'Delimiter' , ',')
end
legend('90% coverage' , '70% coverage' , '50% coverage' , 'No vaccination')

% relative incidence reduction
figure()
for i = 1 : length(hiv_vec)
    d = hiv_vec(i);
    hivCC90 = sum(sum(sum(newCC_90(: , d , : , : , :),3),4),5);
    hivCC70 = sum(sum(sum(newCC_70(: , d , : , : , :),3),4),5);
    hivCC50 = sum(sum(sum(newCC_50(: , d , : , : , :),3),4),5);
    hivCCNo = sum(sum(sum(newCC_0(: , d , : , : , :),3),4),5);

    hivSus = [toInd(allcomb(d , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , 4 : age , 1 : risk)) ; toInd(allcomb(d , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
    pop90_susHiv = sum(pop90(1 : end - 1 , hivSus) , 2);
    pop70_susHiv = sum(pop70(1 : end - 1 , hivSus) , 2);
    pop50_susHiv = sum(pop50(1 : end - 1 , hivSus) , 2);
    popNo_susHiv = sum(popNo(1 : end - 1 , hivSus) , 2);

    h_gen90 = hivCC90(2 : end) ./ pop90_susHiv .* fac;
    h_gen70 = hivCC70(2 : end) ./ pop70_susHiv .* fac;
    h_gen50 = hivCC50(2 : end) ./ pop50_susHiv .* fac;
    h_genNo = hivCCNo(2 : end) ./ popNo_susHiv .* fac;

    hivRelRed_90 = (h_gen90 - h_genNo) ./ h_genNo * 100;
    hivRelRed_70 = (h_gen70 - h_genNo) ./ h_genNo * 100;
    hivRelRed_50 = (h_gen50 - h_genNo) ./ h_genNo * 100;

    subplot(3 , 2 , i)
    plot(tVec(2 : end) , hivRelRed_90 , tVec(2 : end) , hivRelRed_70 , tVec(2 : end) , hivRelRed_50)
    title(['Cervical Cancer Reduction ', tits{i}])
    xlabel('Year'); ylabel('Relative Difference (%)')
    axis([tVec(1) , tVec(end) , -100 , 0])
    % Export HIV-positive reduction as csv
    T = table(tVec(2 : end)' , hivRelRed_90 , hivRelRed_70 , hivRelRed_50);
    writetable(T , [filenames{i} , '_Incidence.csv'] , 'Delimiter' , ',')
end
legend('90% coverage' , '70% coverage' , '50% coverage')

%% Aggregate (without ART)

hivAllCC90 = sum(sum(sum(sum(newCC_90(: , 2 : 6 , : , : , :),2),3),4),5);
hivAllCC70 = sum(sum(sum(sum(newCC_70(: , 2 : 6 , : , : , :),2),3),4),5);
hivAllCC50 = sum(sum(sum(sum(newCC_50(: , 2 : 6 , : , : , :),2),3),4),5);
hivAllCCNo = sum(sum(sum(sum(newCC_0(: , 2 : 6 , : , : , :),2),3),4),5);

hivAllSus = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
    2 , 4 : age , 1 : risk)); toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes ,...
    8 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
pop90_susHiv = sum(pop90(1 : end - 1 , hivAllSus) , 2);
pop70_susHiv = sum(pop70(1 : end - 1 , hivAllSus) , 2);
pop50_susHiv = sum(pop50(1 : end - 1 , hivAllSus) , 2);
popNo_susHiv = sum(popNo(1 : end - 1 , hivAllSus) , 2);

hivAllInc90 = hivAllCC90(2 : end) ./ pop90_susHiv * fac;
hivAllInc70 = hivAllCC70(2 : end) ./ pop70_susHiv * fac;
hivAllInc50 = hivAllCC50(2 : end) ./ pop50_susHiv * fac;
hivAllIncNo = hivAllCCNo(2 : end) ./ popNo_susHiv * fac;
figure()
plot(tVec(2 : end) , hivAllInc90 , tVec(2 :end) , hivAllInc70 , ...
    tVec(2 : end) , hivAllInc50 , tVec(2 : end) , hivAllIncNo)
title('Cervical Cancer Incidence Among All HIV+')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('90% coverage' , '70% coverage' , '50% coverage' , 'No vaccination')

hivAllRed_90 = (hivAllInc90 - hivAllIncNo) ./ hivAllIncNo * 100;
hivAllRed_70 = (hivAllInc70 - hivAllIncNo) ./ hivAllIncNo * 100;
hivAllRed_50 = (hivAllInc50 - hivAllIncNo) ./ hivAllIncNo * 100;
figure()
plot(tVec(2 : end) , hivAllRed_90 , tVec(2 :end) , hivAllRed_70 , ...
    tVec(2 : end) , hivAllRed_50)
title('Cervical Cancer Reduction Among All HIV+')
xlabel('Year'); ylabel('Reduction (%)')
legend('90% coverage' , '70% coverage' , '50% coverage')
axis([tVec(1) , tVec(end) , -100 , 0])
% Export HIV+ aggregate cervical cancer incidence and reduction as csv
T = table(tVec(2 : end)' , hivAllInc90 , hivAllInc70 , hivAllInc50 , ...
    hivAllIncNo , hivAllRed_90 , hivAllRed_70 , hivAllRed_50);
    writetable(T , 'AllHiv_Incidence.csv' , 'Delimiter' , ',')
%% HIV-
% incidence
figure()
hivNegCC90 = sum(sum(sum(newCC_90(: , 1 , : , : , :),3),4),5);
hivNegCC70 = sum(sum(sum(newCC_70(: , 1 , : , : , :),3),4),5);
hivNegCC50 = sum(sum(sum(newCC_50(: , 1 , : , : , :),3),4),5);
hivNegCCNo = sum(sum(sum(newCC_0(: , 1 , : , : , :),3),4),5);

hivNegSus = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
    2 , 4 : age , 1 : risk)) ; toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , ...
    9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
pop90_susHivNeg = sum(pop90(1 : end - 1 , hivNegSus) , 2);
pop70_susHivNeg = sum(pop70(1 : end - 1 , hivNegSus) , 2);
pop50_susHivNeg = sum(pop50(1 : end - 1 , hivNegSus) , 2);
popNo_susHivNeg = sum(popNo(1 : end - 1 , hivNegSus) , 2);


h_neg90 = hivNegCC90(2 : end) ./ pop90_susHivNeg .* fac;
h_neg70 = hivNegCC70(2 : end) ./ pop70_susHivNeg .* fac;
h_neg50 = hivNegCC50(2 : end) ./ pop50_susHivNeg .* fac;
h_negNo = hivNegCCNo(2 : end) ./ popNo_susHivNeg .* fac;
plot(tVec(2 : end) , h_neg90 , tVec(2 : end) , h_neg70 , tVec(2 : end) , ...
    h_neg50 , tVec(2 : end) , h_negNo)
title('Cervical Cancer Incidence in HIV-')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('90% coverage' , '70% coverage' , '50% coverage' , 'No vaccination')
% Export HIV-negative incidence as csv
T = table(tVec(2 : end)' , h_neg90 , h_neg70 , h_neg50 , h_negNo);
writetable(T , 'HIVNeg_Incidence.csv' , 'Delimiter' , ',')

% relative incidence reduction
figure()
hivSus = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
    2 , 4 : age , 1 : risk)) ; toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , ...
    9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
pop90_susHivNeg = sum(pop90(1 : end - 1 , hivSus) , 2);
pop70_susHivNeg = sum(pop70(1 : end - 1 , hivSus) , 2);
pop50_susHivNeg = sum(pop50(1 : end - 1 , hivSus) , 2);
popNo_susHivNeg = sum(popNo(1 : end - 1 , hivSus) , 2);

h_gen90 = hivNegCC90(2 : end) ./ pop90_susHivNeg .* fac;
h_gen70 = hivNegCC70(2 : end) ./ pop70_susHivNeg .* fac;
h_gen50 = hivNegCC50(2 : end) ./ pop50_susHivNeg .* fac;
h_genNo = hivNegCCNo(2 : end) ./ popNo_susHivNeg .* fac;

hivNegRelRed_90 = (h_neg90 - h_negNo) ./ h_negNo * 100;
hivNegRelRed_70 = (h_neg70 - h_negNo) ./ h_negNo * 100;
hivNegRelRed_50 = (h_neg50 - h_negNo) ./ h_negNo * 100;

plot(tVec(2 : end) , hivNegRelRed_90 , tVec(2 : end) , hivNegRelRed_70 , ...
    tVec(2 : end) , hivNegRelRed_50)
title('Cervical Cancer Reduction among HIV-')
xlabel('Year'); ylabel('Relative Difference (%)')
axis([tVec(1) , tVec(end) , -100 , 0])
legend('90% coverage' , '70% coverage' , '50% coverage')

% Export HIV-negative reduction as csv
T = table(tVec(2 : end)' , hivNegRelRed_90 , hivNegRelRed_70 , hivNegRelRed_50);
writetable(T , 'HIVNeg_Reduction.csv' , 'Delimiter' , ',')

%% Waning
v90_w20 = load('H:\HHCoM Results\Vax_0.9_wane_20.mat');
v90_w15 = load('U:\HHCoM Results\Vax_0.9_wane_15.mat');
v90_w10 = load('U:\HHCoM Results\Vax_0.9_wane_10.mat');
v90_w0 = load('U:\HHCoM Results\Vax_0.9_wane_0.mat');
v0_w0 = load('U:\HHCoM Results\Vax_0_wane_0.mat');
%% Deaths
% general
allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 2 , 4 : age , 1 : risk));
% All HIV-positive women (not on ART)
allHivF = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 2 , 4 : age , 1 : risk));
% All HIV-negative women
hivNeg = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    2 , 4 : age , 1 : risk));
inds = {':' , 2 : 6 , 1};
genArray = {allF , allHivF , hivNeg};
files = {'General_CCMortality' , 'allHiv_CCMortality' , 'hivNegMortality'};
plotTits = {'General Cervical Cancer' , 'HIV+ Cervical Cancer' , 'HIV- Cervical Cancer'};
for i = 1 : length(genArray)
    vNo_wNoMort = ...
        sum(sum(sum(sum(v0_w0.ccDeath(2 : end , inds{i} , : , : , :),2),3),4),5) ./ ...
        sum((v0_w0.popVec(1 : end - 1 , genArray{i}) ...
        + v0_w0.popVec(2 : end , genArray{i})) * 0.5 , 2) * fac;

    v90_w20Mort = ...
        sum(sum(sum(sum(v90_w20.ccDeath(2 : end , inds{i} , : , : , :),2),3),4),5) ./ ...
        sum((v90_w20.popVec(1 : end - 1 , genArray{i}) ...
        + v90_w20.popVec(2 : end , genArray{i})) * 0.5 , 2) * fac;

    v90_w15Mort = ...
        sum(sum(sum(sum(v90_w15.ccDeath(2 : end , inds{i} , : , : , :),2),3),4),5) ./ ...
        sum((v90_w15.popVec(1 : end - 1 , genArray{i}) ...
        + v90_w15.popVec(2 : end , genArray{i})) * 0.5 , 2) * fac;

    v90_w10Mort = ...
        sum(sum(sum(sum(v90_w10.ccDeath(2 : end , inds{i} , : , : , :),2),3),4),5) ./ ...
        sum((v90_w10.popVec(1 : end - 1 , genArray{i}) ...
        + v90_w10.popVec(2 : end , genArray{i})) * 0.5 , 2) * fac;

    v90_wNoMort = ...
        sum(sum(sum(sum(v90_w0.ccDeath(2 : end , inds{i} , : , : , :),2),3),4),5) ./ ...
        sum((v90_w0.popVec(1 : end - 1 , genArray{i}) ...
        + v90_w0.popVec(2 : end , genArray{i})) * 0.5 , 2) * fac;
    figure()
    plot(tVec(2 : end) , vNo_wNoMort , tVec(2 : end) , v90_wNoMort , ...
        tVec(2 : end) , v90_w20Mort , tVec(2 : end) , v90_w15Mort , ...
        tVec(2 : end) , v90_w10Mort)
    title([plotTits{i} , ' Mortality'])
    xlabel('Year'); ylabel('Mortality per 100,000')
    legend('No vaccination' , 'No Waning' , '20 years' , '15 years' , '10 years')
    % Reduction
    v90_wNoRed = (v90_wNoMort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w20Red = (v90_w20Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w15Red = (v90_w15Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w10Red = (v90_w10Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;

    figure()
    plot(tVec(2 : end) , v90_wNoRed , tVec(2 : end) , v90_w20Red , ...
        tVec(2 : end) , v90_w15Red , tVec(2 : end) , v90_w10Red)
    title([plotTits{i} , ' Mortality Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('No Waning' , '20 years' , '15 years' , '10 years')
    axis([tVec(2) tVec(end) -100 0])

    T = table(tVec(2 : end)' , v90_w20Mort , v90_w15Mort , v90_w10Mort , ...
        v90_wNoMort , v90_wNoRed , v90_w20Red , v90_w15Red , v90_w10Red);
    writetable(T , ['waning_', files{i} , '.csv'] , 'Delimiter' , ',')
end
%%
% By CD4
dVec = [2 : 6 , 10];
tits = {'Acute' , 'CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , 'CD4 < 200' , ...
    'ART'};
filenames = {'AcuteMort' , 'CD4_500Mort' , 'CD4_500_350Mort' , 'CD4_350_200Mort' , ...
    'CD4_200Mort' , 'ARTMort'};
for d = 1 : length(dVec)
    hiv_ccSus = [toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , 4 : age , 1 : risk)) ; toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];

    vNo_wNoMort = ...
        sum(sum(sum(sum(v0_w0.ccDeath(2 : end , dVec(d) , : , : , :),2),3),4),5) ./ ...
        sum((v0_w0.popVec(1 : end - 1 , hiv_ccSus) ...
        + v0_w0.popVec(2 : end , hiv_ccSus)) * 0.5 , 2) * fac;

    v90_w20Mort = ...
        sum(sum(sum(sum(v90_w20.ccDeath(2 : end , dVec(d) , : , : , :),2),3),4),5) ./ ...
        sum((v90_w20.popVec(1 : end - 1 , hiv_ccSus) ...
        + v90_w20.popVec(2 : end , hiv_ccSus)) * 0.5 , 2) * fac;

    v90_w15Mort = ...
        sum(sum(sum(sum(v90_w15.ccDeath(2 : end , dVec(d) , : , : , :),2),3),4),5) ./ ...
        sum((v90_w15.popVec(1 : end - 1 , hiv_ccSus) ...
        + v90_w15.popVec(2 : end , hiv_ccSus)) * 0.5 , 2) * fac;

    v90_w10Mort = ...
        sum(sum(sum(sum(v90_w10.ccDeath(2 : end , dVec(d) , : , : , :),2),3),4),5) ./ ...
        sum((v90_w10.popVec(1 : end - 1 , hiv_ccSus) ...
        + v90_w10.popVec(2 : end , hiv_ccSus)) * 0.5 , 2) * fac;

    v90_wNoMort = ...
        sum(sum(sum(sum(v90_w0.ccDeath(2 : end , dVec(d) , : , : , :),2),3),4),5) ./ ...
        sum((v90_w0.popVec(1 : end - 1 , hiv_ccSus) ...
        + v90_w0.popVec(2 : end , hiv_ccSus)) * 0.5 , 2) * fac;
    figure(102)
    subplot(3 , 2 , d)
    plot(tVec(2 : end) , vNo_wNoMort , tVec(2 : end) , v90_wNoMort , ...
        tVec(2 : end) , v90_w20Mort , tVec(2 : end) , v90_w15Mort , ...
        tVec(2 : end) , v90_w10Mort)
    title([tits{i} , ' Mortality'])
    xlabel('Year'); ylabel('Mortality per 100,000')
    legend('No vaccination' , 'No Waning' , '20 years' , '15 years' , '10 years' ,...
        'Location' , 'northeastoutside')
    % Reduction
    v90_wNoRed = (v90_wNoMort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w20Red = (v90_w20Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w15Red = (v90_w15Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w10Red = (v90_w10Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;

    figure(103)
    subplot(3 , 2 ,d)
    plot(tVec(2 : end) , v90_wNoRed , tVec(2 : end) , v90_w20Red , ...
        tVec(2 : end) , v90_w15Red , tVec(2 : end) , v90_w10Red)
    title([tits{d} , ' Mortality Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('No Waning' , '20 years' , '15 years' , '10 years' , ...
        'Location' , 'northeastoutside')
    axis([tVec(2) tVec(end) -100 0])

    T = table(tVec(2 : end)' , v90_w20Mort , v90_w15Mort , v90_w10Mort , ...
        v90_wNoMort , v90_wNoRed , v90_w20Red , v90_w15Red , v90_w10Red);
    writetable(T , [filenames{d} , '.csv'] , 'Delimiter' , ',')
end


%%
% General susceptibles
ccSus = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
    2 , 4 : age , 1 : risk)) ; toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
    9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
v90_w0_inc  = sum(sum(sum(sum(v90_w0.newCC(2 : end , : , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w0.popVec(1 : end - 1 , ccSus) , 2) * fac;
v90_w20_inc  = sum(sum(sum(sum(v90_w20.newCC(2 : end , : , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w20.popVec(1 : end - 1 , ccSus) , 2) * fac;
v90_w15_inc  = sum(sum(sum(sum(v90_w15.newCC(2 : end , : , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w15.popVec(1 : end - 1 , ccSus) , 2) * fac;
v90_w10_inc  = sum(sum(sum(sum(v90_w10.newCC(2 : end , : , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w10.popVec(1 : end - 1 , ccSus) , 2) * fac;

figure()
plot(tVec(2 : end) , i_genNo , tVec(2 : end) , v90_w0_inc , tVec(2 : end) , v90_w20_inc , ...
    tVec(2 : end) , v90_w15_inc , tVec(2 : end) , v90_w10_inc)
title('Vaccine Waning Period and General Cervical Cancer Incidence')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('No vaccination' , 'No waning' , '20 years' , '15 years' , '10 years')

% Relative reduction
wNo = (v90_w0_inc - i_genNo) ./ i_genNo * 100;
w20 = (v90_w20_inc - i_genNo) ./ i_genNo * 100;
w15 = (v90_w15_inc - i_genNo) ./ i_genNo * 100;
w10 = (v90_w10_inc - i_genNo) ./ i_genNo * 100;
figure()
plot(tVec(2 : end) , wNo , tVec(2 : end) , w20 , ...
    tVec(2 : end) , w15 , tVec(2 : end) , w10)
axis([tVec(2) tVec(end) -100 0])
title('Vaccine Waning Period and Reduction of General Cervical Cancer Incidence')
xlabel('Year'); ylabel('Reduction (%)')
legend('No waning' , '20 years' , '15 years' , '10 years')
% Export general incidence/reduction as csv
T = table(tVec(2 : end)' , v90_w0_inc , v90_w20_inc , v90_w15_inc , v90_w10_inc , ...
    wNo , w20 , w15 , w10);
writetable(T , 'GenInc_Waning.csv' , 'Delimiter' , ',')
%% HIV-negative
% incidence
ccNegSus = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
    2 , 4 : age , 1 : risk)) ; toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , ...
    9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
v90_w0_incNeg  = sum(sum(sum(sum(v90_w0.newCC(2 : end , 1 , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w0.popVec(1 : end - 1 , ccNegSus) , 2) * fac;
v90_w20_incNeg  = sum(sum(sum(sum(v90_w20.newCC(2 : end , 1 , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w20.popVec(1 : end - 1 , ccNegSus) , 2) * fac;
v90_w15_incNeg  = sum(sum(sum(sum(v90_w15.newCC(2 : end , 1 , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w15.popVec(1 : end - 1 , ccNegSus) , 2) * fac;
v90_w10_incNeg  = sum(sum(sum(sum(v90_w10.newCC(2 : end , 1 , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w10.popVec(1 : end - 1 , ccNegSus) , 2) * fac;

figure()
plot(tVec(2 : end) , h_genNo , tVec(2 : end) , v90_w0_incNeg , tVec(2 : end) ,...
    v90_w20_incNeg , tVec(2 : end) , v90_w15_incNeg , tVec(2 : end) , v90_w10_incNeg)
title('Vaccine Waning Period and HIV- Cervical Cancer Incidence')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('No vaccination' , 'No waning' , '20 years' , '15 years' , '10 years')

% Relative reduction
wNoNeg = (v90_w0_incNeg - h_genNo) ./ h_genNo * 100;
w20Neg = (v90_w20_incNeg - h_genNo) ./ h_genNo * 100;
w15Neg = (v90_w15_incNeg - h_genNo) ./ h_genNo * 100;
w10Neg = (v90_w10_incNeg - h_genNo) ./ h_genNo * 100;
figure()
plot(tVec(2 : end) , wNoNeg , tVec(2 : end) , w20Neg , ...
    tVec(2 : end) , w15Neg , tVec(2 : end) , w10Neg)
axis([tVec(2) tVec(end) -100 0])
title('Vaccine Waning Period and Reduction of HIV- Cervical Cancer Incidence')
xlabel('Year'); ylabel('Reduction (%)')
legend('No waning' , '20 years' , '15 years' , '10 years')
% Export general incidence/reduction as csv
T = table(tVec(2 : end)' , v90_w0_incNeg , v90_w20_incNeg , v90_w15_incNeg , ...
    v90_w10_incNeg , wNoNeg , w20Neg , w15Neg , w10Neg);
writetable(T , 'HIVNeg_Waning.csv' , 'Delimiter' , ',')
%% HIV-positive
%% Aggregated
% Incidence
hiv_ccSus = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
    2 , 4 : age , 1 : risk)) ; toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , ...
    9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
v90_w0_incHiv  = sum(sum(sum(sum(v90_w0.newCC(2 : end , 2 : 6 , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w0.popVec(1 : end - 1 , hiv_ccSus) , 2) * fac;
v90_w20_incHiv  = sum(sum(sum(sum(v90_w20.newCC(2 : end , 2 : 6 , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w20.popVec(1 : end - 1 , hiv_ccSus) , 2) * fac;
v90_w15_incHiv  = sum(sum(sum(sum(v90_w15.newCC(2 : end , 2 : 6 , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w15.popVec(1 : end - 1 , hiv_ccSus) , 2) * fac;
v90_w10_incHiv  = sum(sum(sum(sum(v90_w10.newCC(2 : end , 2 : 6 , : , : , :)...
    ,2),3),4),5) ./ sum(v90_w10.popVec(1 : end - 1 , hiv_ccSus) , 2) * fac;
figure()
plot(tVec(2 : end) , hivAllIncNo , tVec(2 : end) , v90_w0_incHiv , tVec(2 : end) , ...
    v90_w20_incHiv , tVec(2 : end) , v90_w15_incHiv , tVec(2 : end) , v90_w10_incHiv)
title('Vaccine Waning Period and General Cervical Cancer Incidence in HIV+')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('No vaccination' , 'No waning' , '20 years' , '15 years' , '10 years' , ...
    'Location' , 'northeastoutside')

% Relative reduction
wNoRed = (v90_w0_incHiv - hivAllIncNo) ./ hivAllIncNo * 100;
w20Red = (v90_w20_incHiv - hivAllIncNo) ./ hivAllIncNo * 100;
w15Red = (v90_w15_incHiv - hivAllIncNo) ./ hivAllIncNo * 100;
w10Red = (v90_w10_incHiv - hivAllIncNo) ./ hivAllIncNo * 100;
figure()
plot(tVec(2 : end) , wNoRed , tVec(2 : end) , w20Red , tVec(2 : end) , ...
    w15Red , tVec(2 : end) , w10Red)
title('Vaccine Waning Period and General Cervical Cancer Reduction in HIV+')
xlabel('Year'); ylabel('Reduction(%)')
axis([tVec(2) tVec(end) -100 0 ]);
legend('No vaccination' , 'No waning' , '20 years' , '15 years' , '10 years' , ...
    'Location' , 'northeastoutside')
% Export general incidence/reduction as csv
T = table(tVec(2 : end)' ,hivAllIncNo , v90_w0_incHiv , v90_w20_incHiv , ...
    v90_w15_incHiv , v90_w10_incHiv , wNoRed , w20Red , w15Red , w10Red);
writetable(T , 'AllHivInc_Waning.csv' , 'Delimiter' , ',')

%% By CD4
dVec = [2 : 6 , 10];
tits = {'Acute' , 'CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , 'CD4 < 200' , ...
    'ART'};
filenames = {'Acute' , 'CD4_500' , 'CD4_500_350' , 'CD4_350_200' , 'CD4_200' , ...
    'ART'};
for d = 1 : length(dVec)
    % Incidence
    hiv_ccSus = [toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , 4 : age , 1 : risk)) ; toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
    v90_wNo_incCD4  = sum(sum(sum(v90_w0.newCC(2 : end , dVec(d) , : , : , :)...
        ,3),4),5) ./ sum(v90_w0.popVec(1 : end - 1 , hiv_ccSus) , 2) * fac;
    v90_w20_incCD4  = sum(sum(sum(v90_w20.newCC(2 : end , dVec(d) , : , : , :)...
        ,3),4),5) ./ sum(v90_w20.popVec(1 : end - 1 , hiv_ccSus) , 2) * fac;
    v90_w15_incCD4  = sum(sum(sum(v90_w15.newCC(2 : end , dVec(d) , : , : , :)...
        ,3),4),5) ./ sum(v90_w15.popVec(1 : end - 1 , hiv_ccSus) , 2) * fac;
    v90_w10_incCD4  = sum(sum(sum(v90_w10.newCC(2 : end , dVec(d) , : , : , :)...
        ,3),4),5) ./ sum(v90_w10.popVec(1 : end - 1 , hiv_ccSus) , 2) * fac;
    % base (no vaccine)
    hivCCNo = sum(sum(sum(newCC_0(: , dVec(d) , : , : , :),3),4),5);
    popNo_susHiv = sum(popNo(1 : end - 1 , hiv_ccSus) , 2);
    h_genNo = hivCCNo(2 : end) ./ popNo_susHiv .* fac;

    figure(100)
    subplot(3 , 2 , d)
    plot(tVec(2 : end) , h_genNo , tVec(2 : end) , v90_wNo_incCD4 , tVec(2 : end) , ...
        v90_w20_incCD4 , tVec(2 : end) , v90_w15_incCD4 , tVec(2 : end) , v90_w10_incCD4)
    title(['Vaccine Waning Period and CC Incidence in ' , tits{d}])
    xlabel('Year'); ylabel('Incidence per 100,000')
    legend('No vaccination' , 'No waning' , '20 years' , '15 years' , '10 years' , ...
        'Location' , 'northeastoutside')

    % Relative reduction
    wNoRedCD4 = (v90_wNo_incCD4 - h_genNo) ./ h_genNo * 100;
    w20RedCD4 = (v90_w20_incCD4 - h_genNo) ./ h_genNo * 100;
    w15RedCD4 = (v90_w15_incCD4 - h_genNo) ./ h_genNo * 100;
    w10RedCD4 = (v90_w10_incCD4 - h_genNo) ./ h_genNo * 100;
    figure(101)
    subplot(3 , 2 , d)
    plot(tVec(2 : end) , wNoRedCD4 , tVec(2 : end) , w20RedCD4 , tVec(2 : end) , ...
        w15RedCD4 , tVec(2 : end) , w10RedCD4)
    axis([tVec(2) tVec(end) -100 0 ]);
    title(['Vaccine Waning Period and CC Reduction in ' tits{d}])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('No waning' , '20 years' , '15 years' , '10 years' , ...
        'Location' , 'northeastoutside')
    % Export general incidence/reduction as csv
    T = table(tVec(2 : end)' , v90_wNo_incCD4 , v90_w20_incCD4 , ...
        v90_w15_incCD4 , v90_w10_incCD4 , wNoRedCD4 , w20RedCD4 , w15RedCD4 , ...
        w10RedCD4);
    writetable(T , [filenames{d} , '_Waning.csv'] , 'Delimiter' , ',')
end
