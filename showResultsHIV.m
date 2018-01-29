function[] = showResultsHIV()
%%
load('actual')
load('calibData')
% load('C:\Users\nicktzr\Google Drive\ICRC\CISNET\Results\to2017')
load('H:\HHCoM_Results\hiv_to2017')
%% Plot Settings

colors = [241, 90, 90;
          240, 196, 25;
          78, 186, 111;
          45, 149, 191;
          149, 91, 165]/255;

set(groot, 'DefaultAxesColor', [10, 10, 10]/255);
set(groot, 'DefaultFigureColor', [10, 10, 10]/255);
set(groot, 'DefaultFigureInvertHardcopy', 'off');
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on')
set(groot, 'DefaultAxesColorOrder', colors);
set(groot, 'DefaultLineLineWidth', 3);
set(groot, 'DefaultTextColor', [1, 1, 1]);
set(groot, 'DefaultAxesXColor', [1, 1, 1]);
set(groot, 'DefaultAxesYColor', [1, 1, 1]);
set(groot , 'DefaultAxesZColor' , [1 , 1 ,1]);
set(0,'defaultAxesFontSize',14)
ax = gca;
ax.XGrid = 'on';
ax.XMinorGrid = 'on';
ax.YGrid = 'on';
ax.YMinorGrid = 'on';
ax.GridColor = [1, 1, 1];
ax.GridAlpha = 0.4;
% reset(0)
% set(0 , 'defaultlinelinewidth' , 2)
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
for g = 1 : gender
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : 10 , 1 : risk));
    artPop = sum(popVec(: , artInds) , 2);
    art = sum(popVec(: , artInds) , 2);
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hivPop = sum(popVec(: , hivInds) , 2);
    plot(tVec , 100 * artPop ./ (hivPop + art)) 
    hold on
end
hold on
plot(yrsArtActual , artActual , '*')
xlabel('Year')
ylabel('Proportion of HIV Population (%)')
title('Proportion on ART')
legend('Male' , 'Female' , 'General (Observed)')

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
