c = fix(clock);
currYear = c(1); % get the current year
yearNow = round((currYear - startYear) * stepsPerYear);
% Total HIV positive
hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : gender , 4 : 10 , 1 : risk));
hivPopBase = sum(artBase.popVec(: , hivInds) , 2);
hivPopMed = sum(artMed.popVec(: , hivInds) , 2);
hivPopHigh = sum(artHigh.popVec(: , hivInds) , 2);
hivPopMH = sum(artMH.popVec(: , hivInds) , 2);

artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : gender , 4 : 10 , 1 : risk));
artPopBase = sum(artBase.popVec(: , artInds) , 2);
artPopMed = sum(artMed.popVec(: , artInds) , 2);
artPopHigh = sum(artHigh.popVec(: , artInds) , 2);
artPopMH = sum(artMH.popVec(: , artInds) , 2);

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
totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    1 : gender , 4 : 10, 1 : risk));
basePopTot = artBase.popVec(: , totInds);
medPopTot = artMed.popVec(: , totInds);
highPopTot = artHigh.popVec(: , totInds);
mhPopTot = artMH.popVec(: , totInds);

basePrev = 100 * (hivPopBase + artPopBase) ./ sum(basePopTot , 2);
medPrev = 100 * (hivPopMed + artPopMed) ./ sum(medPopTot , 2);
highPrev = 100 * (hivPopHigh + artPopHigh) ./ sum(highPopTot , 2);
mhPrev = 100 * (hivPopMH + artPopMH) ./ sum(mhPopTot , 2);

plot(tVec , medPrev  , tVec , mhPrev , tVec , highPrev , tVec , basePrev)
hold on
plot(overallHivPrev_KZN_AC(1 , :) , overallHivPrev_KZN_AC(2 , :) , '*')
xlabel('Year'); ylabel('Proportion of Population (%)'); title('HIV Prevalence (Ages 15-49)')
xlim([tVec(1) , tVec(end)])
legend('65% Both' , 'M:65%, F:80%' , '80% Both' , 'M:45%, F:65%' , 'KZN Actual (Africa Center Data)')

%% Change in HIV incidence
susInds = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    1 : gender , 4 : 10, 1 : risk));

annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) ...
    / stepsPerYear)); % calculates total for year

% population at risk (average within year)
basePopSus = annlz(sum(artBase.popVec(1 : end , susInds) , 2)) ./ stepsPerYear;
medPopSus = annlz(sum(artMed.popVec(1 : end , susInds) , 2)) ./ stepsPerYear;
highPopSus = annlz(sum(artHigh.popVec(1 : end , susInds) , 2)) ./ stepsPerYear;
mhPopSus = annlz(sum(artMH.popVec(1 : end , susInds) , 2)) ./ stepsPerYear;

fac = 10 ^ 2;

baseInc = annlz(sum(sum(sum(artBase.newHiv(1 : end , 1 : gender , 4 : 10 , 1 : risk), 2), 3 ), 4))...
    ./ basePopSus * fac;
medInc = annlz(sum(sum(sum(artMed.newHiv(1 : end , 1 : gender , 4 : 10 , 1 : risk), 2), 3 ), 4))...
    ./ medPopSus * fac;
highInc = annlz(sum(sum(sum(artHigh.newHiv(1 : end , 1 : gender , 4 : 10 , 1 : risk), 2), 3 ), 4))...
    ./ highPopSus * fac;
mhInc = annlz(sum(sum(sum(artMH.newHiv(1 : end , 1 : gender , 4 : 10 , 1 : risk), 2), 3 ), 4))...
    ./ mhPopSus * fac;

figure()
plot(tVec(1 : stepsPerYear : end) , medInc , ...
    tVec(1 : stepsPerYear : end) , mhInc , ...
    tVec(1 : stepsPerYear : end) , highInc , ...
    tVec(1 : stepsPerYear : end) , baseInc);
xlim([tVec(yearNow) - 1 , tVec(end)])
xlabel('Year'); ylabel('Incidence per 100');
title('HIV Incidence')
legend('65% Both' , 'M:65%, F:80%' , '80% Both' , 'M:45%, F:65%')

deltaMed = (medInc - baseInc) ./ baseInc * 100;
deltaHigh = (highInc - baseInc) ./ baseInc * 100;
deltaMH = (mhInc - baseInc) ./ baseInc * 100;

figure()
plot(tVec(1 : stepsPerYear : end) , deltaMed ,  ...
    tVec(1 : stepsPerYear : end) , deltaMH , ...
    tVec(1 : stepsPerYear : end) , deltaHigh);
xlim([tVec(yearNow) - 1 , tVec(end)])
xlabel('Year'); ylabel('Change (%)');
title('Change in Incidence')
legend('65% Both' , 'M:65%, F:80%' , '80% Both')

% Change in HIV-related mortality
fac = 10 ^ 2;

baseMort = annlz(sum(sum(artBase.hivDeaths(1 : end , 1 : gender , 1 : age), 2), 3)) ...
    ./ basePopSus * fac;
medMort = annlz(sum(sum(artMed.hivDeaths(1 : end  , 1 : gender , 1 : age) , 2), 3)) ...
    ./ medPopSus * fac;
highMort = annlz(sum(sum(artHigh.hivDeaths(1 : end , 1 : gender , 1 : age) , 2), 3)) ...
    ./ highPopSus * fac;
mhMort = annlz(sum(sum(artMH.hivDeaths(1 : end , 1 : gender , 1 : age) , 2), 3))...
    ./ mhPopSus * fac;

figure()
plot(tVec(1 : stepsPerYear : end) , medMort , ...
    tVec(1 : stepsPerYear : end) , mhMort , ...
    tVec(1 : stepsPerYear : end) , highMort , ...
    tVec(1 : stepsPerYear : end) , baseMort);
xlim([tVec(yearNow) - 1 , tVec(end)])
xlabel('Year'); ylabel('Mortality per 100');
title(['HIV Mortality'])
legend('65% Both' , 'M:65%, F:80%' , '80% Both' , 'M:45%, F:65%')

deltaMortMed = (medMort - baseMort) ./ baseMort * 100;
deltaMortHigh = (highMort - baseMort) ./ baseMort * 100;
deltaMortMH = (mhMort - baseMort) ./ baseMort * 100;

figure()
plot(tVec(1 : stepsPerYear : end) , deltaMortMed ,...
    tVec(1 : stepsPerYear : end) , deltaMortMH , ...
    tVec(1 : stepsPerYear : end) , deltaMortHigh);
xlim([tVec(yearNow) - 1 , tVec(end)])
xlabel('Year'); ylabel('Change (%)');
title(['Change in Mortality'])
legend('65% Both' , 'M:65%, F:80%' , '80% Both')

%% summary table
yr_2030 = (2030 - startYear);% * stepsPerYear;
yr_2040 = (2040 - startYear);% * stepsPerYear;
yr_2050 = (2050 - startYear);% * stepsPerYear - 1;

% Change in mortality
deltaM_Med = [deltaMortMed(yr_2030) ; deltaMortMed(yr_2040) ; ...
    deltaMortMed(yr_2050)];
deltaM_High = [deltaMortHigh(yr_2030) ; deltaMortHigh(yr_2040) ; ...
    deltaMortHigh(yr_2050)];
deltaM_MH = [deltaMortMH(yr_2030) ; deltaMortMH(yr_2040) ; ...
    deltaMortMH(yr_2050)];

% Change in incidence
deltaIMed = [deltaMed(yr_2030) ; deltaMed(yr_2040) ; ...
    deltaMed(yr_2050)];
deltaIHigh = [deltaHigh(yr_2030) ; deltaHigh(yr_2040) ; ...
    deltaHigh(yr_2050)];
deltaIMH = [deltaMH(yr_2030) ; deltaMH(yr_2040) ; ...
    deltaMH(yr_2050)];

baseP = [basePrev(yr_2030) ; basePrev(yr_2040) ; basePrev(yr_2050)];
medP = [medPrev(yr_2030) ; medPrev(yr_2040) ; medPrev(yr_2050)];
highP = [highPrev(yr_2030) ; highPrev(yr_2040) ; highPrev(yr_2050)];
mhP = [mhPrev(yr_2030) ; mhPrev(yr_2040) ; mhPrev(yr_2050)];

hivSumm = table({'2030' ; '2040' ; '2050'} , baseP , medP, mhP , highP ,...
    deltaIMed , deltaIMH , deltaIHigh , deltaM_Med , deltaM_MH , deltaM_High);
hivSumm.Properties.VariableNames{1} = 'Year';
disp('Impact Summary for General Population')
disp(hivSumm)
filename = ['General_artCompare.xlsx'];
writetable(hivSumm , filename)
disp(['Results for general population written to ' , filename])