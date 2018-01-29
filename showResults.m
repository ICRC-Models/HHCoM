function[] = showResults()
%%
load('actual')
load('calibData')
% load('C:\Users\nicktzr\Google Drive\ICRC\CISNET\Results\to2017')
load('H:\HHCoM_Results\to2017')
%% Plot Settings

% colors = [241, 90, 90;
%           240, 196, 25;
%           78, 186, 111;
%           45, 149, 191;
%           149, 91, 165]/255;
% 
% set(groot, 'DefaultAxesColor', [10, 10, 10]/255);
% set(groot, 'DefaultFigureColor', [10, 10, 10]/255);
% set(groot, 'DefaultFigureInvertHardcopy', 'off');
% set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')
% set(groot, 'DefaultAxesColorOrder', colors);
% set(groot, 'DefaultLineLineWidth', 3);
% set(groot, 'DefaultTextColor', [1, 1, 1]);
% set(groot, 'DefaultAxesXColor', [1, 1, 1]);
% set(groot, 'DefaultAxesYColor', [1, 1, 1]);
% set(groot , 'DefaultAxesZColor' , [1 , 1 ,1]);
% set(0,'defaultAxesFontSize',14)
% ax = gca;
% ax.XGrid = 'on';
% ax.XMinorGrid = 'on';
% ax.YGrid = 'on';
% ax.YMinorGrid = 'on';
% ax.GridColor = [1, 1, 1];
% ax.GridAlpha = 0.4;
reset(0)
set(0 , 'defaultlinelinewidth' , 2)
%% Plot results
% gropInds();
% load('groupedInds');u
% Total HIV positive
hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk));
hivPop = sum(popVec(: , hivInds) , 2);
artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk));
art = sum(popVec(: , artInds) , 2);
% Compared to Africa Center data
overallHivPrev_KZN_AC(1 , :) = 1990 : 2009;
overallHivPrev_KZN_AC(2 , :) = [0.464072571
    0.985438052
    1.506803533
    5.576509907
    8.126044402
    13.04177608
    12.54905705
    16.61876343
    19.50632609
    19.52064932
    22.2391979
    20.22439723
    22.09787539
    22.78825495
    25.16877536
    26.19622822
    25.36548102
    27.2380043
    27.42134161
    28.44974934];
figure()
popTot = popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    1 : 2 , 4 : 10, 1 : risk)));
plot(tVec , (hivPop + art) ./ sum(popTot , 2) * 100 , overallHivPrev_KZN_AC(1 , :) , overallHivPrev_KZN_AC(2 , :) , '*')
xlabel('Year'); ylabel('Proportion of Population (%)'); title('HIV Prevalence (Ages 15-49)')
legend('Model' , 'KZN Actual (Africa Center Data)')
%% HIV status by age
% hivAge = zeros(age , length(tVec));
ageGroup = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
figure()
for g = 1 : gender
    for a = 1 : age
        hivAgeInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a , 1 : risk));
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a , 1 : risk));
        %         hivAge(a , :) = sum(popVec(: , hivAgeInds) , 2);
        hivAgeRel = bsxfun(@rdivide , sum(popVec(: , hivAgeInds) , 2)' , sum(popVec(: , ageInds) , 2)') * 100;
        subplot(4 , 4 , a)
        hold on
        plot(tVec , hivAgeRel);
        xlabel('Year'); ylabel('% HIV'); title([' Age group ' , ageGroup{a} , ' HIV Prevalence'])
    end
end
legend('Male' , 'Female')


hivM = hivPrevM_obs(: , 2) ./ hivPrevM_obs(: , 3) .* 100;
hivF = hivPrevF_obs(: , 2) ./hivPrevF_obs(: , 3) .* 100;
prevYears = unique(hivPrevF_obs(: , 1));
gen = {'Male' , 'Female'};
for g = 1 : gender
    hivPrevs = hivM;
    if g == 2
        hivPrevs = hivF;
    end
    figure()
    for a = 4 : 10
        hivAgeInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a , 1 : risk));
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a , 1 : risk));
        hivAge(a , :) = sum(popVec(: , hivAgeInds) , 2);
        hivAgeRel = bsxfun(@rdivide , hivAge(a , :)' , sum(popVec(: , ageInds) , 2)) * 100;
        subplot(4 , 2 , a - 3)
        plot(tVec , hivAgeRel , prevYears , hivPrevs((a - 3) : 7 : end) ,'o');
        xlabel('Year'); ylabel('% HIV'); title([gen{g} , ' Age group ' , ageGroup{a} , ' HIV Prevalence'])
    end
    legend('Model' , 'Africa Center Data')
end

%% HIV mortality by age
figure()
bar(tVec , hivDeaths , 'stacked')
xlabel('Year'); ylabel('Deaths'); title('HIV-associated Deaths')
legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthWest')
figure()
bar(tVec , bsxfun(@rdivide , hivDeaths , sum(popVec , 2)) , 'stacked')
xlabel('Year'); ylabel('Deaths relative to population'); title('HIV-associated Deaths')
legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthWest')
%% Gender
mInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 1 , 1 : age , 1 : risk));
fInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 2 , 1 : age , 1 : risk));
plot(tVec , sum(popVec(: , mInds) , 2))
hold on
plot(tVec , sum(popVec(: , fInds) , 2))
legend('Males' , 'Females')
xlabel('Year'); ylabel('Population')
%% Total on ART
figure()
artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 1 : gender , 4 : 10 , 1 : risk));
artPop = sum(popVec(: , artInds) , 2);
hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : gender , 4 : 10 , 1 : risk));
hivPop = sum(popVec(: , hivInds) , 2);
hiv_art = [100 * hivPop ./ sum(popTot , 2), 100 * artPop ./ sum(popTot , 2)];
area(tVec , hiv_art); %art ./ sum(popVec , 2) , tVec , hiv ./ sum(popVec , 2))
xlabel('Year')
ylabel('Proportion of Population (%)')
title('Relative HIV Prevalence')
legend('Untreated', 'On ART' , 'Location' , 'NorthWest')
%%
figure()
artActual = [0	0	0	0	1	2	3	6, ...
    9	14	19	27	34	40	45	48];
yrsArtActual = [2000	2001	2002	2003	2004	2005 ...
    2006	2007	2008	2009	2010	2011	2012	2013 ...
    2014	2015];

plot(tVec , 100 * artPop ./ (hivPop + art) , yrsArtActual , artActual , '*')
xlabel('Year')
ylabel('Proportion of HIV Population')
title('Proportion on ART')
legend('Model' , 'Observed')
%% Overall HPV
genders = {'Male' , 'Female'};
figure()
for g = 1 : 2
    % General
    hpvInds = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : 7, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hpvPop = sum(popVec(: , hpvInds) , 2);
    popTot = popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : 10 , 1 : risk)));
    % HIV+
    hpvHivInds = toInd(allcomb(2 : 6 , 1 : viral , 2 : 4 , 1 : 7, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hpvHivPop = sum(popVec(: , hpvHivInds) , 2);
    popHivTot = popVec(: , toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : 10 , 1 : risk)));
    %ART
    hpvArtInds = toInd(allcomb(10, 6 , 2 : 4 , 1 : 7, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hpvArtPop = sum(popVec(: , hpvArtInds) , 2);
    popArtTot = popVec(: , toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : 10 , 1 : risk)));
    %HIV-
    hpvHivNegInds = toInd(allcomb(1 , 1 , 2 : 4 , 1 : 7, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hpvHivNegPop = sum(popVec(: , hpvHivNegInds) , 2);
    popHivNegTot = popVec(: , toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : 10 , 1 : risk)));

    subplot(2 , 1 , g)
    plot(tVec , 100 * hpvPop ./ sum(popTot , 2))
    hold on
    plot(tVec , 100 * hpvHivNegPop ./ sum(popHivNegTot , 2))
    plot(tVec , 100 * hpvHivPop ./ sum(popHivTot , 2))
    plot(tVec , 100 * hpvArtPop ./ sum(popArtTot , 2))
    axis([tVec(1) tVec(end) 0 100])
    % axis([-inf inf 0 min(max(hpvPop./sum(popTot , 2) * 100) + 10 , 100)])
    xlabel('Year'); ylabel('Prevalence (%)'); title([genders{g} , ' HPV Prevalence (15-49)'])
    legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Location' , 'NorthEastOutside')
end

%% General HPV Prevalence by Age in 2017
ageGroup = {'15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
hpv2017 = zeros(age - 4 + 1 , 1);
hpvHIV2017 = hpv2017;
hpvNeg2017 = hpv2017;
hpvImm2017 = hpv2017;

for a = 4 : age
    hpvInds = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : 2, ...
        1 : periods , 2 , a , 1 : risk));
    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));
    hpv2017(a - 3) = sum(popVec((2017 - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((2017 - startYear) * stepsPerYear , ageInds)) * 100;

    % immune
    hpvImmInds = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 10, ...
        1 : periods , 2 , a , 1 : risk));
    hpvImm2017(a - 3) = sum(popVec((2017 - startYear) * stepsPerYear , hpvImmInds))...
        ./ sum(popVec((2017 - startYear) * stepsPerYear , ageInds)) * 100;
end
figure()
plot(1 : length(hpv2017) , hpv2017 , 'o-')
set(gca , 'xtick' , 1 : length(hpv2017) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('Age Specific hrHPV Prevalence in 2017')
figure()
plot(1 : length(hpvImm2017) , hpvImm2017 , 'o-')
set(gca , 'xtick' , 1 : length(hpvImm2017) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('hrHPV Immunity Prevalence in 2017')

for a = 4 : age
    hrInds = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : 2, ...
        1 : periods , 2 , a , 1 : risk));

    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));
    hpv2017(a - 3) = (sum(popVec((2017 - startYear) * stepsPerYear , hrInds)))...
        ./ sum(popVec((2017 - startYear) * stepsPerYear , ageInds)) * 100;
    
    % HIV+
    hrHIVInds = [toInd(allcomb(2 : 6 , 1 : viral , 2 : 4 , 1 : 2, ...
        1 : periods , 2 , a , 1 : risk));toInd(allcomb(10 , 6 , 2 : 4 , 1 : 2, ...
        1 : periods , 2 , a , 1 : risk));];

    ageHIVInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk)); toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , 2 , a , 1 : risk))];
    hpvHIV2017(a - 3) = (sum(popVec((2017 - startYear) * stepsPerYear , hrHIVInds)))...
        ./ sum(popVec((2017 - startYear) * stepsPerYear , ageHIVInds)) * 100;
    
    % HIV-
    hrNegInds = toInd(allcomb(1 , 1 , 2 : 4 , 1 : 2, ...
        1 : periods , 2 , a , 1 : risk));

    ageNegInds = toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));
    hpvNeg2017(a - 3) = (sum(popVec((2017 - startYear) * stepsPerYear , hrNegInds)))...
        ./ sum(popVec((2017 - startYear) * stepsPerYear , ageNegInds)) * 100;

    % immune
    hrImmInds = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 10 , ...
        1 : periods , 2 , a , 1 : risk));

    hpvImm2017(a - 3) = (sum(popVec((2017 - startYear) * stepsPerYear , hrImmInds)))...
        ./ sum(popVec((2017 - startYear) * stepsPerYear , ageInds)) * 100;
end

% McDonald 2014
hpvHivObs(: , 1) = [0.75
0.61
0.60
0.55
0.46
0.42
0.43
0.54
0.35];

hpvHivObs(: , 2) = [0.63
0.54
0.54
0.47
0.42
0.34
0.32
0.35
0.16];

hpvHivObs(: ,3) = [0.87
0.67
0.66
0.62
0.51
0.50
0.55
0.72
0.53];

hpvNegObs(: , 1) = [0.60
0.38
0.24
0.20
0.19
0.18
0.13
0.17
0.15];

hpvNegObs(: , 2) = [0.53
0.34
0.21
0.17
0.18
0.16
0.11
0.14
0.12];

hpvNegObs(: , 3) = [0.67
0.41
0.27
0.23
0.21
0.20
0.15
0.19
0.18];

hpvHivObs = hpvHivObs * 100;
hpvNegObs = hpvNegObs * 100;
figure()
plot(1 : length(hpv2017) , hpv2017 , 'co-')
hold on
plot(1 : length(hpvHIV2017) , hpvHIV2017 , 'bo-');
hold on
plot(1 : length(hpvNeg2017) , hpvNeg2017 , 'o-')
set(gca , 'xtickLabel' , ageGroup);

% general
% yPosError = abs(hrHpvObs(: , 3) - hrHpvObs(: , 1));
% yNegError = abs(hrHpvObs(: , 2) - hrHpvObs(: , 1));
% errorbar(1 : length(hrHpvObs) , hrHpvObs(: , 1) , yNegError , yPosError , 'rs')
% HIV+
yPosError = abs(hpvHivObs(: , 3) - hpvHivObs(: , 1));
yNegError = abs(hpvHivObs(: , 2) - hpvHivObs(: , 1));
errorbar(1 : length(hpvHivObs) , hpvHivObs(: , 1) , yNegError , yPosError , 'bs')
%HIV-
hold on
yPosError = abs(hpvNegObs(: , 3) - hpvNegObs(: , 1));
yNegError = abs(hpvNegObs(: , 2) - hpvNegObs(: , 1));
errorbar(1 : length(hpvNegObs) , hpvNegObs(: , 1) , yNegError , yPosError , 'rs')

set(gca , 'xtick' , 1 : length(hpvNegObs) , 'xtickLabel' , ageGroup);
legend('General' , 'HIV+' , 'HIV-' , 'McDonald 2014 - HIV+' , 'McDonald 2014 - HIV-')
xlabel('Age Group'); ylabel('Prevalence (%)')
title('Age Specific hrHPV Prevalence in 2017')


figure()
plot(1 : length(hpvImm2017) , hpvImm2017 , 'o-')
set(gca , 'xtick' , 1 : length(hpvNegObs) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('hrHPV Immunity Prevalence in 2017')

%% HIV-HPV Prevalence by Age in 2017
yearPrev = 2017;
hpvHivPosInd = [toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvTypes , 1 : 4 , ...
    1 : periods , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(10, 6 ,...
    2 : hpvTypes , 1 : 4 , 1 : periods , 2 , 4 : 13 , 1 : risk))];
hpvHivNegInd = toInd(allcomb(1 , 1 , 2 : hpvTypes , 1 : 4 , 1 : periods , 2, ...
    4 : 13 , 1 : risk));
genHpvInd = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 4 , ...
    1 : periods , 2 , 4 : 13 , 1 : risk));

popHivInd = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 2 , 4 : 13 , 1 : risk)); toInd(allcomb(10, 6 ,...
    1 : hpvTypes , 1 : hpvStates , 1 : periods , 2 , 4 : 13 , 1 : risk))];
popNegInd = toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 2 , 4 : 13 , 1 : risk));
popGenInd = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 2 , 4 : 13 , 1 : risk));

popGen = sum(popVec((yearPrev - startYear) * stepsPerYear , popGenInd));
popHiv = sum(popVec((yearPrev - startYear) * stepsPerYear , popHivInd));
popNeg = sum(popVec((yearPrev - startYear) * stepsPerYear , popNegInd));


hpvHivPos = sum(popVec((yearPrev - startYear) * stepsPerYear , hpvHivPosInd))...
    ./ popHiv * 100;
hpvHivNeg = sum(popVec((yearPrev - startYear) * stepsPerYear , hpvHivNegInd))...
    ./ popNeg * 100;
genHpv = sum(popVec((yearPrev - startYear) * stepsPerYear , genHpvInd))...
    ./ popGen * 100;

%Mbulawa KZN
kzn_HpvData = [0.39 , 0.35 , 0.44; % gen hpv
    0.52, 0.46 , 0.58; % hpvHivPos
    0.23 , 0.17 , 0.28] .* 100; %hpvHivNeg

upper = abs(kzn_HpvData(: , 2) - kzn_HpvData(: , 1));
lower = abs(kzn_HpvData(: , 3) - kzn_HpvData(: , 1));
hpvGroup = {'General Population' , 'HIV+' , 'HIV-'};
figure()
bar([genHpv , hpvHivPos , hpvHivNeg])
hold on
errorbar(1 : 3, kzn_HpvData(: , 1) , lower , upper , 'rs') 
set(gca , 'xtick' , 1 : 3 , 'xtickLabel' , hpvGroup);
legend('Model' , 'Mbulawa')
title(['HPV Prevalence by HIV Status in ' , num2str(yearPrev)])
ylabel('Prevalence (%)')

%%
figure()
popGen = sum(popVec(: , popGenInd) , 2);
popHiv = sum(popVec(: , popHivInd) , 2);
popNeg = sum(popVec(: , popNegInd) , 2);
hpvHivPos = sum(popVec(: , hpvHivPosInd) , 2) ./ popHiv * 100;
hpvHivNeg = sum(popVec(:, hpvHivNegInd) , 2) ./ popNeg * 100;
genHpv = sum(popVec(: , genHpvInd) , 2) ./ popGen * 100; 
plot(tVec , genHpv)
hold on
plot(tVec , hpvHivPos)
hold on
plot(tVec , hpvHivNeg)
hold on
% plot([tVec(1) , tVec(end)] , [kzn_HpvData(1 , 1) , kzn_HpvData(1 , 1)] , 'r-'); 
plot(tVec , ones(length(tVec) , 1) .* kzn_HpvData(1 , 1) , 'r-')
hold on
plot(tVec , ones(length(tVec) , 1) .* kzn_HpvData(2 , 1) , 'b-')
hold on
plot(tVec , ones(length(tVec) , 1) .* kzn_HpvData(3 , 1) , 'c-')
xlabel('Year'); ylabel('HPV Prevalence'); title('HPV Prevalence by HIV Status')
legend('General Population' , 'HIV+' , 'HIV-' , 'Mbulawa (General)' , ...
    'Mbulawa (HIV+)' , 'Mbulawa (HIV-)' , 'location' , 'NorthEastOutside')
%% CIN2/3 prevalence by HIV status
cinPos2017 = zeros(hpvTypes - 1 , age - 4 + 1 , 1);
cinNeg2017 = cinPos2017;
ageGroup = {'15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
for a = 4 : age
    for h = 2 : hpvTypes
    % HIV+
    cinInds = toInd(allcomb(2 : 6 , 1 : viral , h , 3 : 4, ...
        1 : periods , 2 , a , 1 : risk));
    ageInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));
    cinPos2017(h - 1 , a - 3) = (sum(popVec((2017 - startYear) * stepsPerYear , cinInds)))...
        ./ sum(popVec((2017 - startYear) * stepsPerYear , ageInds)) * 100;

    cinNegInds = toInd(allcomb(1, 1 : viral , h , 3 : 4, ...
        1 : periods , 2 , a , 1 : risk));
    cinNeg2Inds = toInd(allcomb(7 : disease , 1 : viral , h , 3 : 4, ...
        1 : periods , 2 , a , 1 : risk));
    ageNegInds = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));
    ageNeg2Inds = toInd(allcomb(7 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));
    cinNeg2017(h - 1 , a - 3) = (sum(popVec((2017 - startYear) * stepsPerYear , cinNegInds))...
        + sum(popVec((2017 - startYear) * stepsPerYear , cinNeg2Inds)))...
        ./ (sum(popVec((2017 - startYear) * stepsPerYear , ageNegInds)) +...
        sum(popVec((2017 - startYear) * stepsPerYear , ageNeg2Inds))) * 100;
    end
end

% McDonald 2014
cinPosAct(: , 1) = [0.125
0.054
0.128
0.154
0.081
0.054
0.079
0.071
0.077
0.077]; % mean

cinPosAct(: , 2) = [0.03
0.02
0.09
0.10
0.05
0.02
0.02
0.00
0.00
0.00
]; % lb

cinPosAct(: , 3) = [0.22
0.08
0.17
0.21
0.11
0.09
0.14
0.17
0.22
0.22]; % ub

figure()
subplot(2 , 1 , 1)
cinPosAct = cinPosAct .* 100; % convert to %
yPosError = abs(cinPosAct(: , 3) - cinPosAct(: , 1));
yNegError = abs(cinPosAct(: , 2) - cinPosAct(: , 1));
plot(1 : length(cinPos2017) , cinPos2017 ,'o-')
hold on
errorbar(1 : length(cinPosAct) , cinPosAct(: , 1) , yNegError , yPosError , 'rs')
legend('16/18' , 'Non-4v HR' , 'oHR' , 'McDonald 2014')
set(gca , 'xtick' , 1 : length(cinPosAct) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('Age Specific CIN2/3 Prevalence Among HIV+ in 2017')


cinNegAct(: , 1) = [0.016
0.027
0.021
0.036
0.029
0.031
0.031
0.021
0.014
0.014]; % mean

cinNegAct(: , 2) = [0.00
0.02
0.01
0.02
0.02
0.02
0.02
0.01
0.00
0.00]; % lb

cinNegAct(: , 3) = [0.03
0.04
0.03
0.05
0.04
0.04
0.04
0.03
0.03
0.03]; % ub


subplot(2 , 1 , 2)
cinNegAct = cinNegAct .* 100; % convert to %
plot(1 : length(cinNeg2017) , cinNeg2017 , 'o-')
hold on
yPosError = abs(cinNegAct(: , 3) - cinNegAct(: , 1));
yNegError = abs(cinNegAct(: , 2) - cinNegAct(: , 1));
errorbar(1 : length(cinNegAct) , cinNegAct(: , 1) , yNegError , yPosError , 'rs')
legend('16/18' , 'Non-4V HR' , 'oHR' , 'McDonald 2014')
set(gca , 'xtick' , 1 : length(cinNegAct) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('Age Specific CIN2/3 Prevalence Among HIV- in 2017')

% CIN Prevalence among HIV+ Women

cinHiv = [0.46	0.18 0.09] .* 100; % Observed, Johannesburg, Cynthia Firnhaber
cinHiv2017 = zeros(3 , 1);
for a = 4 : age
    for cin = 2 : 4
        % HIV+
        cinInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , cin - 1, ...
            1 : periods , 2 , a , 1 : risk));
        ageInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            2 , a , 1 : risk));
        cinArtInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , cin - 1, ...
            1 : periods , 2 , a , 1 : risk));
        ageArtInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , cin - 1 , 1 : periods , ...
            2 , a , 1 : risk));
        cinHiv2017(cin - 1) = (sum(popVec((2017 - startYear) * stepsPerYear , cinInds)) + sum(popVec((2017 - startYear) * stepsPerYear , cinArtInds)))...
            ./ (sum(popVec((2017 - startYear) * stepsPerYear , ageInds)) + sum(popVec((2017 - startYear) * stepsPerYear , ageArtInds)))* 100;
    end
end

figure()
% cinGroup = {'CIN 1' , 'CIN 2' , 'CIN 3'};
plot(1 : length(cinHiv2017) , cinHiv2017 , 'o-' , 1 : length(cinHiv) , cinHiv , 'rs');
ylabel('Prevalence (%)')
% set(gca , 'xtickLabel' , cinGroup);
xlabel('CIN Stage')
legend('Model' , 'Firnhaber')
title('CIN Prevalence in HIV Positive Women')

%% CIN2/3 2008
% cinPos2008 = zeros(age - 4 + 1 , 1);
% cinNeg2008 = cinPos2008;
% ageGroup = {'15 - 19' , '20 -24' , '25 - 29' ,...
%     '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
%     '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
% for a = 4 : age
%     % HIV+
%     cinInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 3 : 4, ...
%         1 : periods , 2 , a , 1 : risk));
%     ageInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , a , 1 : risk));
%     cinPos2008(a - 3) = (sum(popVec((2008 - startYear) * stepsPerYear , cinInds)))...
%         ./ sum(popVec((2008 - startYear) * stepsPerYear , ageInds)) * 100;
% 
%     cinNegInds = toInd(allcomb(1, 1 : viral , 1 : hpvTypes , 3 : 4, ...
%         1 : periods , 2 , a , 1 : risk));
%     cinNeg2Inds = toInd(allcomb(7 : disease , 1 : viral , 1 : hpvTypes , 3 : 4, ...
%         1 : periods , 2 , a , 1 : risk));
%     ageNegInds = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , a , 1 : risk));
%     ageNeg2Inds = toInd(allcomb(7 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         2 , a , 1 : risk));
%     cinNeg2008(a - 3) = (sum(popVec((2008 - startYear) * stepsPerYear , cinNegInds))...
%         + sum(popVec((2008 - startYear) * stepsPerYear , cinNeg2Inds)))...
%         ./ (sum(popVec((2008 - startYear) * stepsPerYear , ageNegInds)) +...
%         sum(popVec((2008 - startYear) * stepsPerYear , ageNeg2Inds))) * 100;
% end
% 
% % Allan 2008
% cinPosAct(: , 1) = [0.125 % mean
%     0.054
%     0.128
%     0.154
%     0.081
%     0.054
%     0.079
%     0.071
%     0.077];
% 
% cinPosAct(: , 2) = [0.03 % lb
%     0.02
%     0.09
%     0.10
%     0.05
%     0.02
%     0.02
%     0.00
%     0.00];
% 
% cinPosAct(: , 3) = [0.22 % ub
%     0.08
%     0.17
%     0.21
%     0.11
%     0.09
%     0.14
%     0.17
%     0.22];
% 
% figure()
% subplot(2 , 1 , 1)
% cinPosAct = cinPosAct .* 100; % convert to %
% yPosError = abs(cinPosAct(: , 3) - cinPosAct(: , 1));
% yNegError = abs(cinPosAct(: , 2) - cinPosAct(: , 1));
% plot(1 : length(cinPos2008) , cinPos2008 ,'o-')
% hold on
% errorbar(1 : length(cinPosAct) , cinPosAct(: , 1) , yNegError , yPosError , 'rs')
% legend('Model' , 'Allan 2008')
% set(gca , 'xtick' , 1 : length(cinPosAct) , 'xtickLabel' , ageGroup);
% xlabel('Age Group'); ylabel('Prevalence (%)')
% title('Age Specific CIN2/3 Prevalence Among HIV+ in 2008')
% 
% 
% cinNegAct(: , 1) = [0.016
%     0.027
%     0.021
%     0.036
%     0.029
%     0.031
%     0.031
%     0.021
%     0.014];
% 
% cinNegAct(: , 2) = [0.00
%     0.02
%     0.01
%     0.02
%     0.02
%     0.02
%     0.02
%     0.01
%     0.00];
% 
% cinNegAct(: , 3) = [0.03
%     0.04
%     0.03
%     0.05
%     0.04
%     0.04
%     0.04
%     0.03
%     0.03];
% 
% 
% subplot(2 , 1 , 2)
% cinNegAct = cinNegAct .* 100; % convert to %
% plot(1 : length(cinNeg2008) , cinNeg2008 , 'o-')
% hold on
% yPosError = abs(cinNegAct(: , 3) - cinNegAct(: , 1));
% yNegError = abs(cinNegAct(: , 2) - cinNegAct(: , 1));
% errorbar(1 : length(cinNegAct) , cinNegAct(: , 1) , yNegError , yPosError , 'rs')
% legend('Model' , 'Allan 2008')
% set(gca , 'xtick' , 1 : length(cinNegAct) , 'xtickLabel' , ageGroup);
% xlabel('Age Group'); ylabel('Prevalence (%)')
% title('Age Specific CIN2/3 Prevalence Among HIV- in 2008')
% 
% % CIN Prevalence among HIV+ Women
% 
% cinHiv = [0.46	0.18 0.09] .* 100; % Observed, Johannesburg, Cynthia Firnhaber
% cinHiv2008 = zeros(3 , 1);
% for a = 4 : age
%     for cin = 2 : 4
%         % HIV+
%         cinInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , cin - 1, ...
%             1 : periods , 2 , a , 1 : risk));
%         ageInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%             2 , a , 1 : risk));
%         cinArtInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , cin - 1, ...
%             1 : periods , 2 , a , 1 : risk));
%         ageArtInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , cin - 1 , 1 : periods , ...
%             2 , a , 1 : risk));
%         cinHiv2008(cin - 1) = (sum(popVec((2008 - startYear) * stepsPerYear , cinInds)) + sum(popVec((2008 - startYear) * stepsPerYear , cinArtInds)))...
%             ./ (sum(popVec((2008 - startYear) * stepsPerYear , ageInds)) + sum(popVec((2008 - startYear) * stepsPerYear , ageArtInds)))* 100;
%     end
% end
% 
% figure()
% % cinGroup = {'CIN 1' , 'CIN 2' , 'CIN 3'};
% plot(1 : length(cinHiv2008) , cinHiv2008 , 'o-' , 1 : length(cinHiv) , cinHiv , 'rs');
% ylabel('Prevalence (%)')
% % set(gca , 'xtickLabel' , cinGroup);
% xlabel('CIN Stage')
% legend('Model' , 'Firnhaber')
% title('CIN Prevalence in HIV Positive Women')
%% HPV status by age
% hivAge = zeros(age , length(tVec));
ageGroup = {'15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
figure()
for g = 1 : gender
    for a = 4 : age
        hpvAgeInds = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : 7 , 1 : periods , ...
            g , a , 1 : risk));
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a , 1 : risk));
        hpvAgeRel = bsxfun(@rdivide , sum(popVec(: , hpvAgeInds) , 2)' , sum(popVec(: , ageInds) , 2)') * 100;
        subplot(5 , 3 , a - 3)
        hold on
        plot(tVec , hpvAgeRel);
        xlabel('Year'); ylabel('% HPV'); title([' Age group ' , ageGroup{a - 3} , ' HPV Prevalence'])
    end
end
legend('Male' , 'Female')
%% Cervical cancer
%% prevalence
ccAgeRel = zeros(age , 1);
ccAgeNegRel = ccAgeRel;
ccAgePosRel = ccAgeRel;
ccNegPosArt = zeros(age , 2);
ccArtRel = ccAgeRel;
for a = 1 : age
    % Total population
    ccInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 5 : 7 , 1 : periods , ...
        2 , a  , 1 : risk));
    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            2 , a , 1 : risk));
    ccAgeRel(a) = sum(popVec((2017 - startYear) * stepsPerYear , ccInds) , 2) ...
        / sum(popVec((2017 - startYear) * stepsPerYear , ageInds) , 2) * 100;

    % HIV Negative
    ccHivNegInds = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 5 : 7 , 1 : periods , ...
        2 , a  , 1 : risk));
    ageNegInds = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            2 , a , 1 : risk));
    ccAgeNegRel(a) = sum(popVec((2017 - startYear) * stepsPerYear , ccHivNegInds) , 2) ...
        / (sum(popVec((2017 - startYear) * stepsPerYear , ageNegInds) , 2)) * 100;

    % HIV Positive
    ccHivPosInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 5 : 7 , 1 : periods , ...
        2 , a  , 1 : risk));
    agePosInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));
    ccAgePosRel(a) = sum(popVec((2017 - startYear) * stepsPerYear , ccHivPosInds) , 2) ...
        / (sum(popVec((2017 - startYear) * stepsPerYear , agePosInds) , 2)) * 100;

    % On ART
    ccArtInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 5 : 7 , 1 : periods , ...
        2 , a  , 1 : risk));
    ageArtInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));
    ccArtRel(a) = sum(popVec((2017 - startYear) * stepsPerYear , ccArtInds) , 2) ...
        / (sum(popVec((2017 - startYear) * stepsPerYear , ageArtInds) , 2)) * 100;

    % Each group relative to total female population
    % HIV-
    ccNegPosArt(a , 1) = sum(popVec((2017 - startYear) * stepsPerYear , ccHivNegInds) , 2) ...
        / sum(popVec((2017 - startYear) * stepsPerYear , ageInds) , 2) * 100;
    % HIV+
    ccNegPosArt(a , 2) = sum(popVec((2017 - startYear) * stepsPerYear , ccHivPosInds) , 2) ...
        / sum(popVec((2017 - startYear) * stepsPerYear , ageInds) , 2) * 100;
    % ART
    ccNegPosArt(a , 3) = sum(popVec((2017 - startYear) * stepsPerYear , ccArtInds) , 2) ...
        / sum(popVec((2017 - startYear) * stepsPerYear , ageInds) , 2) * 100;

end

figure()
plot(1 : length(ccAgeRel) , ccAgeRel , '-o' , 1 : length(ccAgeNegRel) , ...
    ccAgeNegRel , '-o' , 1 : length(ccAgePosRel) , ccAgePosRel , '-o' , ...
    1 : length(ccArtRel) , ccArtRel , '-o');
xlabel('Age Group'); ylabel('Prevalence (%)')
set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
title('Cervical Cancer Prevalence in 2017')
legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Location' , 'NorthWest')

figure()
bar(1 : length(ccAgeRel) , ccNegPosArt , 'stacked')
hold on
plot(1 : length(ccAgeRel) , ccAgeRel , '-o')
xlabel('Age Group'); ylabel('Prevalence (%)')
set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
title('Cervical Cancer Prevalence in 2017')
legend('HIV-' , 'HIV+' , 'ART' , 'Location' , 'NorthWest')
%% Incidence
ccIncYears = [2017 , 2003];
ccAgeRel = zeros(age , length(ccIncYears));
ccAgeNegRel = ccAgeRel;
ccAgePosRel = zeros(age , 4 , length(ccIncYears));
ccNegPosArt = zeros(age , 3 , length(ccIncYears));
ccNegPosArtTot = ccNegPosArt;
ccArtRel = ccAgeRel;
fScale = 10^5;
ageGroup = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 - 34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
for a = 1 : age
    % Total population
    ageInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 8 : 10 , 1 : periods , ...
        2 , a , 1 : risk))];
    ccAgeRel(a , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , 1 : disease , 1 : viral , 1 : hpvTypes , a) , 2) , 3) , 4) ...
        ./ ((sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , ageInds) , 2) ...
        + sum(popVec(((ccIncYears - startYear)) * stepsPerYear , ageInds) , 2)) ./ 2) * fScale;

    % HIV Negative
    ageNegInds = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)); toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 8 : 10 , 1 : periods , ...
        2 , a , 1 : risk))];
    ccAgeNegRel(a , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , 1 , 1 : viral , 1 : hpvTypes , a) , 2) , 3) , 4) ...
        ./ ((sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , ageNegInds) , 2) ...
        + sum(popVec(((ccIncYears - startYear)) * stepsPerYear , ageNegInds) , 2)) ./ 2) * fScale;
    
    
    % Acute and CD4 > 500
    agePosInds = [toInd(allcomb(2 : 3 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)); toInd(allcomb(2 : 3 , 1 : viral , 1 : hpvTypes , 8 : 10 , 1 : periods , ...
        2 , a , 1 : risk))];
    ccAgePosRel(a , 1 , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , 2 : 3 , 1 : viral , 1 : hpvTypes , a), 2) , 3) , 4) ...
        ./ ((sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , agePosInds) , 2) ...
        + sum(popVec(((ccIncYears - startYear)) * stepsPerYear , agePosInds) , 2)) ./ 2) * fScale;
    for d = 4 : 6
    % HIV Positive CD4 500-350 -> CD4 < 200
    agePosInds = [toInd(allcomb(d , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)); toInd(allcomb(d , 1 : viral , 1 : hpvTypes , 8 : 10 , 1 : periods , ...
        2 , a , 1 : risk))];
    ccAgePosRel(a , d - 2 , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , d , 1 : viral , 1 : hpvTypes , a), 2) , 3) , 4) ...
        ./ ((sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , agePosInds) , 2) ...
        + sum(popVec(((ccIncYears - startYear)) * stepsPerYear , agePosInds) , 2)) ./ 2) * fScale;
    end

    % On ART
    ageArtInds = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)); toInd(allcomb(10 , 6 , 1 : hpvTypes , 8 : 10 , 1 : periods , ...
        2 , a , 1 : risk))];
    ccArtRel(a , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , 10 , 6 , 1 : hpvTypes , a) , 2) , 3) , 4) ...
        ./ ((sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , ageArtInds) , 2) ...
        + sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , ageArtInds) , 2)) ./ 2) * fScale;

    % Proportion of cervical cancers by HIV/ART status and age
    % Total by age
    ageTotal = sum(sum(sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , ...
        toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk))), 2 ) , 3), 4) + sum(sum(sum(popVec((ccIncYears - startYear - 1) * stepsPerYear , ...
        toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 8 : 10 , 1 : periods , ...
        2 , a , 1 : risk))), 2) , 3) , 4);
        
    % HIV-
    ccNegPosArt(a , 1 , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , 1 , 1 : viral , 1 : hpvTypes , a), 2) , 3) , 4) ...
        ./ ageTotal;
    % HIV+
    ccNegPosArt(a , 2 , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , 2 : 6 , 1 : viral , 1 : hpvTypes , a) , 2) , 3) , 4) ...
        ./ ageTotal;
    % ART
    ccNegPosArt(a , 3 , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , 10 , 6 , 1 : hpvTypes , a), 2) , 3) , 4) ...
        ./ ageTotal;
    
    
     ageAllPosInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)); toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 8 : 10 , 1 : periods , ...
        2 , a , 1 : risk))];    
    % HIV-
    ccNegPosArtTot(a , 1 , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , 1 , 1 : viral , 1 : hpvTypes , a), 2) , 3) , 4) ...
        ./ ((sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , ageNegInds) , 2) ...
        + sum(popVec(((ccIncYears - startYear)) * stepsPerYear , ageNegInds) , 2)) ./ 2) * fScale;
    % HIV+
    ccNegPosArtTot(a , 2 , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , 2 : 6 , 1 : viral , 1 : hpvTypes , a) , 2) , 3) , 4) ...
        ./ ((sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , ageAllPosInds) , 2) ...
        + sum(popVec(((ccIncYears - startYear)) * stepsPerYear , ageAllPosInds) , 2)) ./ 2) * fScale;
    % ART
    ccNegPosArtTot(a , 3 , :) = sum(sum(sum(newCC((ccIncYears - startYear) * stepsPerYear , 10 , 6 , 1 : hpvTypes , a), 2) , 3) , 4) ...
        ./ ((sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , ageArtInds) , 2) ...
        + sum(popVec(((ccIncYears - startYear) - 1) * stepsPerYear , ageArtInds) , 2)) ./ 2) * fScale;;
end




globocan =[0
2.646467154
8.848389036
45.1937379
53.40682334
63.4
68.3
70.7
73
77.4
82.7
88.6
95.2];

globocan_ub = [0.00
2.65
8.85
45.19
53.41
67.05
80.83
78.97
128.87
105.27
118.70
111.81
95.20];

globocan_lb = [0.00
0.00
0.41
9.97
12.61
25.00
45.69
36.31
50.55
57.08
62.69
42.43
52.01];

ccNegPosArt = ccNegPosArt .* fScale;

for y = 1 : length(ccIncYears)
    ccIncYear = ccIncYears(y);
    figure()
    plot(1 : size(ccAgeRel , 1) , ccAgeRel(: , y) , '-o' , 1 : length(ccAgeNegRel(: , y)) , ...
        ccAgeNegRel(: , y) , '-o' , 1 : size(ccAgePosRel , 1) , ccAgePosRel(: , : , y) , '-o' , ...
        1 : size(ccArtRel , 1) , ccArtRel(: , y) , '-o');
    xlabel('Age Group'); ylabel('Incidence per 100,000')
    set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
    title(['Cervical Cancer Incidence in ' num2str(ccIncYear)])
%     legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Location' , 'NorthWest');
%     legend('General' , 'HIV-' , ' Acute and CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , ...
%         'CD4 < 200' , 'ART' , 'Location' , 'NorthWest')
    hold on
    % globocan data
    plot(4 : age , globocan , '-' , 4 : age , globocan_ub , 'r--' , 4 : age , globocan_lb , 'r--')
%     legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Globocan' , 'Upper Bound' , 'Lower Bound' , ...
%         'Location' , 'NorthEastOutside')
    legend('General' , 'HIV-' , ' Acute and CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , ...
        'CD4 < 200' , 'ART' , 'Globocan' , 'Upper Bound' , 'Lower Bound' , ...
        'Location' , 'NorthEastOutside')

    figure()
    bar(1 : length(ccNegPosArt(: , : , y)) , ccNegPosArt(: , : , y), 'stacked')
    xlabel('Age Group'); ylabel('Incidence per 100,000')
    set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
    title(['Cervical Cancer Incidence Distribution in ' , num2str(ccIncYear)])
    legend('HIV-' , 'HIV+' , 'ART')
    
    figure()
     plot(1 : size(ccAgeRel , 1) , ccAgeRel(: , y) , '-o' , 1 : size(ccNegPosArtTot, 1) , ccNegPosArtTot(: , 1 , y) , '-o' , 1 : size(ccNegPosArtTot, 1) , ...
        ccNegPosArtTot(: , 2 , y) , '-o' , 1 : size(ccNegPosArtTot, 1) , ccNegPosArtTot(: , 3 , y) , '-o');
    hold on
    plot(4 : age , globocan , '-' , 4 : age , globocan_ub , 'k--' , 4 : age , globocan_lb , 'k--')
    title(['Cervical Cancer Incidence Distribution in ' , num2str(ccIncYear)])
    legend('General' , 'HIV-' , 'HIV+' , 'ART' , 'Globocan' , 'Upper Bound' , 'Lower Bound')
    xlabel('Age Group'); ylabel('Incidence per 100,000')
    set(gca , 'xtick' , 1 : length(ccAgeRel) , 'xtickLabel' , ageGroup);
    title(['Cervical Cancer Incidence in ' num2str(ccIncYear)])
end
%% Cervical cancer incidence type distribution
newCCTotal = sum(sum(sum(sum(newCC(: , : , : , : , :) , 2) , 3) , 4), 5);
newCCType = zeros(size(newCC , 1) , 3);
for h = 2 : hpvTypes
    newCCType(: , h - 1) = sum(sum(sum(newCC(: , : , : , h  , :) , 2) , 3) , 5) ...
        ./ newCCTotal;
end
figure(); area(tVec , newCCType)
legend('HPV 16/18' , 'Non-4v HPV' , 'Non-Vaccine HPV')
title('Cervical Cancer Incidence Type Distribution')
%% cervical cancer incidence by age over time
ageTotal = zeros(length(tVec) , age);
for a = 1 : age
    ageTotal(: , a) = sum(popVec(1 : size(popVec , 1) , ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
                2 , a , 1 : risk))), 2) + sum(popVec(1 : size(popVec , 1) , ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 8 : 10 , 1 : periods , ...
                2 , a , 1 : risk))) , 2);
end

newCCAge = (sum(sum(sum(newCC(1 : size(newCC , 1) , 1 : disease , 1 : viral ,...
    1 : hpvTypes , 1 : age) , 2) , 3) , 4));

ccIncEvo = bsxfun(@rdivide , permute(newCCAge , [1 , 5 , 4 , 3 , 2]) * fScale , ageTotal);
figure()
mesh(1 : age , tVec , ccIncEvo)
set(gca , 'yLim' , [tVec(1) tVec(end)]);
set(gca , 'xtick' , 1 : age , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Year'); zlabel('Incidence per 100,000')
title('Cervical Cancer Incidence')
%%
hold on
[px , py] = meshgrid(1 : age , 2017 ...
    .* ones(age , 1));
pz = bsxfun(@times , ones(size(px , 1) , size(py , 1)) , linspace(0 , max(ccIncEvo(:)) * 1.2 , size(px , 1)));
m = surf(px , py , pz' , 'edgecolor' , 'r');
set(m , 'facecolor' , 'r')
alpha(0.4)
%% cervical cancer incidence by age over time in HIV-
ageTotal = zeros(length(tVec) , age);
for a = 1 : age
    ageTotal(: , a) = sum(popVec(1 : size(popVec , 1) , ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
                2 , a , 1 : risk))), 2) + sum(popVec(1 : size(popVec , 1) , ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 8 : 10 , 1 : periods , ...
                2 , a , 1 : risk))) , 2);
end

newCCAge = (sum(sum(sum(newCC(1 : size(newCC , 1) , 1 , 1 ,...
    1 : hpvTypes , 1 : age) , 2) , 3) , 4));

ccIncEvo = bsxfun(@rdivide , permute(newCCAge , [1 , 5 , 4 , 3 , 2]) * fScale , ageTotal);
figure()
mesh(1 : age , tVec , ccIncEvo)
set(gca , 'yLim' , [tVec(1) tVec(end)]);
set(gca , 'xtick' , 1 : age , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Year'); zlabel('Incidence per 100,000')
title('HIV- Cervical Cancer Incidence')
%%
hold on
[px , py] = meshgrid(1 : age , 2017 ...
    .* ones(age , 1));
pz = bsxfun(@times , ones(size(px , 1) , size(py , 1)) , linspace(0 , max(ccIncEvo(:)) * 1.2 , size(px , 1)));
m = surf(px , py , pz' , 'edgecolor' , 'r');
set(m , 'facecolor' , 'r')
alpha(0.4)
% hold on
% [px , py] = meshgrid(1 : age , 2003 ...
%     .* ones(age , 1));
% pz = bsxfun(@times , ones(size(px , 1) , size(py , 1)) , linspace(0 , max(ccIncEvo(:)) * 1.2 , size(px , 1)));
% m = surf(px , py , pz' , 'edgecolor' , 'c');
% set(m , 'facecolor' , 'c')
% alpha(0.4)
%% Set up recording parameters (optional), and record
% OptionZ.FrameRate=15;OptionZ.Duration=10;OptionZ.Periodic=true;
% CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], 'CCInc',OptionZ)
%% Total population
%%
% Male (0 - 14)
m0_14all = zeros(length(tVec) , 1);
% Female(0 - 14)
f0_14all = zeros(length(tVec) , 1);
% Male (15 -64)
m15_64all = zeros(length(tVec) , 1);
% Female (15- 64)
f15_64all = zeros(length(tVec) , 1);
actualYrs = [1985 , 1996 , 2001 , 2011];
%% Compare simulation population size to actual population size
% extract results to plot against actual data
disp('Retrieving relevant data points to plot')
disp(' ')
for i = 1 : length(tVec)
    m0_14Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
        1 : hpvStates , 1 : periods , 1 , 1 : 3 , 1 : risk));
    m0_14all(i) = sumall(popVec(i , m0_14Inds));

    f0_14Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
        1 : hpvStates , 1 : periods , 2 , 1 : 3 , 1 : risk));
    f0_14all(i) = sumall(popVec(i , f0_14Inds));

    m15_64Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
        1 : hpvStates , 1 :  periods , 1 , 4 : 13 , 1 : risk));
    m15_64all(i) = sumall(popVec(i , m15_64Inds));

    f15_64Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
        1 : hpvStates , 1 : periods , 2 , 4 : 13 , 1 : risk));
    f15_64all(i) = sumall(popVec(i , f15_64Inds));
end

disp('Plotting model results against actual data points')
disp(' ')
%% Total population
% actual values
% males 0 - 14
actual(1 , :) = [1107298
    1526597
    1669704
    1658047];

% females 0 -14
actual(2 , :) = [1108238
    1537145
    1675104
    1621472];

% males 15 - 64
actual(3 , :) = [1549446
    2292625
    2659850
    3046456];

% females 15 - 64
actual(4 , :) = [1655333
    2711917
    3133521
    3433274];
actualYrs = [1985 , 1996 , 2001 , 2011];
popTotal = zeros(1 , size(popVec , 1));
for i = 1 : size(popVec , 1)
    popTotal(i) = sumall(popVec(i , :));
end
figure()
plot(tVec , popTotal , actualYrs , sum(actual , 1) , 'o')
xlabel('Year') ; ylabel('Population'); title('Population Size')
legend('Projected' , 'Actual')
figure()
plot(tVec , m0_14all , actualYrs , actual(1 , :) , 'o')
title('Males (0 - 14)')
xlabel('Year'); ylabel('Population Size')
legend('Projected' , 'Actual')
figure()
plot(tVec , f0_14all , actualYrs , actual(2 , :) , 'o')
title('Females (0 - 14)')
xlabel('Year'); ylabel('Population Size')
legend('Projected' , 'Actual')
figure()
plot(tVec , m15_64all , actualYrs , actual(3 , :) , 'o')
title('Males (15 - 64)')
xlabel('Year'); ylabel('Population Size')
legend('Projected' , 'Actual')
figure()
plot(tVec , f15_64all , actualYrs , actual(4 , :) , 'o')
title('Females (15 - 64)')
xlabel('Year'); ylabel('Population Size')
legend('Projected' , 'Actual')
%%
% Population by age group over time
popByAge = zeros(length(tVec) , age);
for i = 1 : length(tVec)
    for a = 1 : age
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
            1 : periods , 1 : gender , a , 1 : risk));
        popByAge(i , a) = sum(popVec(i , ageInds) , 2);
    end
end
%%
figure()
area(tVec , bsxfun(@rdivide , popByAge , sum(popVec , 2)))
xlabel('Year'); ylabel('Relative Population Size'); title('Population')
legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')

%%
figure()
ages = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
surf(1 : age , tVec , popByAge);
set(gca , 'xtickLabel' , ages);
set(gca , 'xLim' , [1 age]);
set(gca , 'yLim' , [tVec(1) tVec(end)]);
xlabel('Ages'); ylabel('Year'); zlabel('Population Size'); title('Population by Age')
%% Disease incidence
figure()
for g = 1 : 2
    for a = 1 : age
        subplot(4 , 4 , a)
        hivSusInds = toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , a , 1 : risk));
        hivSus = (sum(popVec(1 : end - 1 , hivSusInds) , 2) + (sum(popVec(2 : end , hivSusInds) , 2) ...
            - sum(newHiv(1 : end - 1 , g , a , :) , 4)));
        plot(tVec(1 : end - 1) , bsxfun(@rdivide , sum(newHiv(1 : end - 1 , g , a , :) , 4) , hivSus) * 100)
        hold on
        xlabel('Year'); ylabel('Rate'); title([' Age group ' , ages{a} , ' HIV Incidence per 100'])
    end
end
legend('Male' , 'Female')

figure()
for g = 1 : 2
    hpvSusInds = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 , ...
        1 : periods , g , 1 : age , 1 : risk));
    hpvSus = (sum(popVec(1 : end - 1 , hpvSusInds) , 2) + (sum(popVec(2 : end , hpvSusInds) , 2) ...
        - sum(newHpv(1 : end - 1 , g) , 2)));
    plot(tVec(1 : end - 1) , bsxfun(@rdivide , sum(newHpv(1 : end - 1 , :) , 2) , hpvSus) * 100)
    axis([startYear , endYear , 0 , 100])    
    hold on
end
legend('Male' , 'Female')
xlabel('Year'); ylabel('Rate Per 100'); title('HPV Incidence')

figure()
for g = 1 : 2
    for a = 1 : age
        subplot(4 , 4 , a)
        hpvSusInds = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 : hpvStates , ...
            1 : periods , g , a , 1 : risk));
        hpvSus = (sum(popVec(1 : end - 1 , hpvSusInds) , 2));
        plot(tVec(1 : end - 1) , bsxfun(@rdivide , sum(sum(newHpv(2 : end , g , 1 : disease , a , :), 5) , 3) + sum(sum(newImmHpv(2 : end , g , 1 : disease , a , :), 5 ) , 3) , hpvSus) * 100)
        hold on
        xlabel('Year'); ylabel('Rate Per 100'); title([' Age group ' , ages{a} , ' HPV Incidence'])
        axis([startYear , endYear , 0 , 100])  
    end
end
legend('Male' , 'Female')

% HIV+ HPV incidence
figure()
for g = 1 : 2
    for a = 1 : age
        subplot(4 , 4 , a)
        hpvHivSusInds = toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvStates , ...
            1 : periods , g , a , 1 : risk));
        hpvHivSus = (sum(popVec(1 : end - 1 , hpvHivSusInds) , 2));
        plot(tVec(1 : end - 1) , bsxfun(@rdivide , sum(sum(newHpv(2 : end , g , 2 : 6 , a , :) , 5 ), 3) + sum(sum(newImmHpv(2 : end , g , 2 : 6 , a , :) , 5 ) , 3) , hpvHivSus) * 100)
        hold on
        xlabel('Year'); ylabel('Rate Per 100'); title([' Age group ' , ages{a} , ' HIV+ HPV Incidence'])
        axis([startYear , endYear , 0 , 100])
    end
end
legend('Male' , 'Female')
%%
% HIV+ vs HIV- HPV incidence
figure()
for g = 2
    for a = 1 : age
        subplot(4 , 4 , a)
        hpvSusInds = toInd(allcomb(1 , 1 , 1 , 1 : hpvStates , ...
            1 : periods , g , a , 1 : risk));
        hpvSus = (sum(popVec(1 : end - 1 , hpvSusInds) , 2));
        hpvHivSusInds = toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvStates , ...
            1 : periods , g , a , 1 : risk));
        hpvHivSus = (sum(popVec(1 : end - 1 , hpvHivSusInds) , 2));
        hpvArtSusInds = toInd(allcomb(10 , 6 , 1 , 1 , ...
            1 : periods , g , a , 1 : risk));
        hpvArtSus = (sum(popVec(1 : end - 1 , hpvArtSusInds) , 2));
        plot(tVec(1 : end - 1) , bsxfun(@rdivide , sum(sum(newHpv(2 : end , g , 2 : 6 , a , :) , 5 ), 3) + sum(sum(newImmHpv(2 : end , g , 2 : 6 , a , :) , 5 ) , 3) , hpvHivSus) * 100)
        hold on
        plot(tVec(1 : end - 1) , bsxfun(@rdivide , sum(sum(newHpv(2 : end , g , 1 , a , :), 5) , 3) + sum(sum(newImmHpv(2 : end , g , 1 , a , :), 5 ) , 3) , hpvSus) * 100)
        hold on
        plot(tVec(1 : end - 1) , bsxfun(@rdivide , sum(sum(newHpv(2 : end , g , 10 , a , :), 5) , 3) + sum(sum(newImmHpv(2 : end , g , 10 , a , :), 5 ) , 3) , hpvArtSus) * 100)
        axis([startYear , endYear , 0 , 100])
        xlabel('Year'); ylabel('Rate Per 100'); title([ages{a} , ' Female HPV Incidence'])
    end
end
legend('HIV+' , 'HIV-' , 'ART')
%%
% HIV+ HPV incidence by risk
% for r = 1 : risk
%     figure()
%     for g = 2
%         for a = 1 : age
%             subplot(4 , 4 , a)
%             hpvSusInds = toInd(allcomb(1 , 1 : viral , 1 , 1 : hpvStates , ...
%                 1 : periods , g , a , r));
%             hpvSus = (sum(popVec(1 : end - 1 , hpvSusInds) , 2));
%             hpvHivSusInds = toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvStates , ...
%                 1 : periods , g , a , r));
%             hpvHivSus = (sum(popVec(1 : end - 1 , hpvHivSusInds) , 2));
%             hpvArtSusInds = toInd(allcomb(10 , 1 : viral , 1 , 1 : hpvStates , ...
%                 1 : periods , g , a , r));
%             hpvArtSus = (sum(popVec(1 : end - 1 , hpvArtSusInds) , 2));
%             plot(tVec(1 : end - 1) , bsxfun(@rdivide , sum(sum(newHpv(2 : end , g , 2 : 6 , a , r) , 5 ), 3) + sum(sum(newImmHpv(2 : end , g , 2 : 6 , a , r) , 5 ) , 3) , hpvHivSus) * 100)
%             hold on
%             plot(tVec(1 : end - 1) , bsxfun(@rdivide , sum(sum(newHpv(2 : end , g , 1 : disease , a , r), 5) , 3) + sum(sum(newImmHpv(2 : end , g , 1 : disease , a , r), 5 ) , 3) , hpvSus) * 100)
%             hold on
%             plot(tVec(1 : end - 1) , bsxfun(@rdivide , sum(sum(newHpv(2 : end , g , 10 , a , r), 5) , 3) + sum(sum(newImmHpv(2 : end , g , 10 , a , r), 5 ) , 3) , hpvArtSus) * 100)
%             axis([startYear , endYear , 0 , 100])
%             xlabel('Year'); ylabel('Rate Per 100'); title([ageGroup{a} , ' Female HPV Incidence (Risk ' , num2str(r) ')'])
%         end
%     end
%     legend('HIV+' , 'HIV-' , 'ART')
%     hold off
% end


%% VL breakdown
hivInf = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 1 : gender , 1 : age , 1 : risk));
totalHiv = sum(popVec(: , hivInf) , 2);
viralPop = zeros(length(tVec) , 6);
for v = 1 : viral
    viralGroup = toInd(allcomb(2 : 6 , v , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , 1 : risk));
    viralPop(: , v) = sum(popVec(: , viralGroup) , 2);
end
figure()
area(tVec , bsxfun(@rdivide , viralPop , totalHiv));
legend('Acute Infection' , 'VL < 1000' , 'VL 1,000 - 10,000' , 'VL 10,000 - 50,000' ,...
    'VL > 50,000' , 'Location' , 'NorthEastOutside')
xlabel('Year')
ylabel('Proportion of HIV infected')
figure()
bar(tVec , bsxfun(@rdivide , viralPop , totalHiv) , 'stacked')
legend('Acute Infection' , 'VL < 1000' , 'VL 1,000 - 10,000' , 'VL 10,000 - 50,000' ,...
    'VL > 50,000' , 'Location' , 'NorthEastOutside')
xlabel('Year')
ylabel('Proportion of HIV infected')
%%
bar(tVec , viralPop , 'stacked')
legend('Acute Infection' , 'VL < 1000' , 'VL 1,000 - 10,000' , 'VL 10,000 - 50,000' ,...
    'VL > 50,000' , 'Location' , 'NorthEastOutside')
xlabel('Year')
ylabel('HIV Infected')
%% CD4 breakdown
cd4Pop = zeros(length(tVec) , 5);
for d = 2 : 6
    cd4Group = toInd(allcomb(d , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , 1 : risk));
    cd4Pop(: , d - 1) = sum(popVec(: , cd4Group) , 2);
end
figure()
area(tVec , bsxfun(@rdivide , cd4Pop , totalHiv));
legend('Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
    'CD4 <= 200 cells/uL' , 'Location' , 'NorthEastOutside')
xlabel('Year')
ylabel('Proportion of HIV infected')
figure()
bar(tVec , bsxfun(@rdivide , viralPop , totalHiv) , 'stacked')
legend('Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
    'CD4 <= 200 cells/uL' , 'Location' , 'NorthEastOutside')
xlabel('Year')
ylabel('Proportion of HIV infected')
%%
figure()
bar(tVec , cd4Pop , 'stacked')
legend('Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
    'CD4 <= 200 cells/uL' , 'Location' , 'NorthEastOutside')
xlabel('Year')
ylabel('HIV Infected')
%% HIV by age group
hivAge = zeros(length(tVec) , 12);
for a = 1 : age
    hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , a , 1 : risk));
    hivArt = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , a , 1 : risk));
    hivAge(: , a) = sum(popVec(: , hivPos) , 2) + sum(popVec(: , hivArt) , 2);
end
hivPosAllInd = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
hivArtAllInd = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
hivPosAll = sum(popVec(: , hivPosAllInd) , 2 ) + sum(popVec(: , hivArtAllInd),2);
figure()
subplot(1 , 2 , 1)
area(tVec , bsxfun(@rdivide , hivAge , hivPosAll));
title('HIV Status by Age Group')
xlabel('Year')
ylabel('Proportion of HIV Positive')
legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')

% HIV by risk group
hivRisk = zeros(length(tVec) , risk);
for r = 1 : risk
    hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , r));
    hivArt = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , r));
    hivRisk(: , r) = sum(popVec(: , hivPos) , 2) + sum(popVec(: , hivArt) , 2);
end
subplot(1 , 2 , 2)
area(tVec , bsxfun(@rdivide , hivRisk , hivPosAll));
title('HIV Status by Risk Group')
xlabel('Year')
ylabel('Proportion of HIV Positive')
legend('Low' , 'Medium' , 'High' , 'Location' , 'NorthEastOutside')

%%
figure()
subplot(2,1,1)
area(tVec , bsxfun(@rdivide , hivAge , sum(popVec , 2)));
title('HIV Prevalence by Age Group')
xlabel('Year')
ylabel('Proportion of HIV Positive')
legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside')

subplot(2 , 1 , 2)
area(tVec , bsxfun(@rdivide , hivRisk , sum(popVec , 2)));
title('HIV Prevalence by Risk Group')
xlabel('Year')
ylabel('Proportion of HIV Positive')
legend('Low' , 'Medium' , 'High' , 'Location' , 'NorthEastOutside')


%% HIV by age and risk
hivAgeRisk = zeros(length(tVec) , age , risk);
for a = 1 : age
    for r = 1 : risk
        hivPos = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 1 : gender , a , r));
        hivArt = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 1 : gender , a , r));
        hivAgeRisk(: , a , r) = sum(popVec(: , hivPos) , 2) + sum(popVec(: , hivArt) , 2);
    end
end



%% ART treatment tracker
cd4ARTFrac = zeros(length(tVec) , 5);
for i = 1 : length(tVec)
    currTot = sumall(artTreatTracker(i , 2 : 6 , 1 : 5 , 1 : gender , 1 : age , 1 :risk));
    for d = 2 : 6
        curr = sumall(artTreatTracker(i , d , 1 : 5 , 1 : gender , 1 : age , 1 : risk));
        cd4ARTFrac(i , d - 1) = curr / currTot;
    end
end

figure()
area(tVec , cd4ARTFrac)
legend('Acute Infection' , 'CD4 > 500 cells/uL' , 'CD4 500 - 350 cells/uL' , 'CD4 350-200 cells/uL' ,...
    'CD4 <= 200 cells/uL' , 'Location' , 'NorthEastOutside')
xlabel('Year')
ylabel('Initiated ART')
% %% CD4 vs VL
% cd4Vl = zeros(length(tVec) , 6 , viral);
% dHivPos = [2 : 6 , 10];
% f = zeros(length(tVec) , 1);
% set(gca, 'nextplot','replacechildren', 'Visible','off');
% % mov(1:f) = struct('cdata',[], 'colormap',[]);
% % set(0,'DefaultFigureVisible','off')
% vid = VideoWriter([date '_cd4_vs_vl_heatmap']);
% cd4 = {'Acute infection' , 'CD4 > 500' , 'CD4 500-350' , ...
%     'CD4 350-200' , 'CD4 <= 200' , 'HIV-positive on ART'};
% vl = {'Acute infection' , 'VL < 1000' , 'VL 1000-10,000', 'VL 10,000-50,000' , ...
%     'VL > 50,000' , 'HIV-positive on ART'};
% vlHeatMap = flip(vl);
% for i = 1 : length(tVec)
%     for j = 1 : length(dHivPos)
%         d = dHivPos(j);
%         for v = 1: viral
%             curr = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates ,...
%                 1 : periods , 1 : gender , 1 : age , 1 : risk));
%             cd4Vl(i , j , v) = sumall(popVec(i , curr));
%         end
%     end
% end
% %%
% reset(gca)
% figure()
% open(vid);
% x = 1 : 6;
% y = 1 : 6;
% % [x , y] = meshgrid(x, y);
% % xq = 1 : 0.1 : 6;
% % yq = 1 : 0.1 : 6;
% % [xq , yq] = meshgrid(xq , yq);
% % for i = 1 : length(tVec)
% %     zq = griddata(x , y , squeeze(cd4Vl(i , : , :)) , xq , yq , 'nearest');
% %     mesh(xq , yq , zq)
% %     set(gca , 'xtick' , 1 : 6);
% %     set(gca, 'XTick', 1 : 6, 'XTickLabel', cd4)
% %
% %     set(gca , 'ytick' , 1 : 6);
% %     set(gca, 'YTick', 1 : 6, 'YTickLabel', vl)
% %
% %     set(gca , 'yTickLabelRotation' , 15)
% % %     set(gca , 'zLim' , [0 max(max(cd4Vl(i , : , :))) + 100]);
% %     set(gca , 'zLim' , [0 max(cd4Vl(:))]);
% %     title(['Year : ' , num2str(floor(tVec(i)))])
% %     writeVideo(vid , getframe(gcf));
% % end
%
% for i = 1 : length(tVec)
%     %     %     s = stem3(1 : 6 , 1 : 6 , squeeze(cd4Vl(i , : , :)) , '--');
%     %     %     s.Color = 'red';
%     %     %     s.MarkerFaceColor = 'red';
%     %     %     s.MarkerSize = 8;
%     %     %     set(gca , 'xtick' , 6 : 1);
%     %     %     set(gca , 'ytick' , 1 : 6);
%     %     %     set(gca , 'yTickLabelRotation' , 15)
%     %     %     set(gca , 'zLim' , [1 max(cd4Vl(:))]);
%     colormap('hot')
%     grid on
%     caxis([1 max(cd4Vl(:))])
%     colorbar
%     imagesc(flipud(squeeze(cd4Vl(i , : , :))))
%     set(gca, 'XTick', 1 : 6, 'XTickLabel', cd4)
%     set(gca , 'XTickLabelRotation' , 90)
%     set(gca, 'YTick', 1 : 6, 'YTickLabel', vlHeatMap)
%     title(['Year : ' , num2str(floor(tVec(i)))])
%     writeVideo(vid , getframe(gcf));
% end
% close(gcf);
% close(vid);
% winopen([date '_cd4_vs_vl_heatmap.avi'])
% set(0,'DefaultFigureVisible','on')
% reset(gca)
%%
% for f = 1 : numel(figHandles)
%     %   saveas(figHandles(f),sprintf('figure_%d.jpg',f))
%     baseFileName = sprintf('figure_%d.jpg',f);
%     % Specify some particular, specific folder:
%     fullFileName = fullfile('C:\Users\nicktzr\Google Drive\ICRC\CISNET\Model\Figures', baseFileName);
%     figure(f); % Activate the figure again.
%     export_fig(fullFileName); % Using export_fig instead of saveas.
% end
%%
% resultOut()
