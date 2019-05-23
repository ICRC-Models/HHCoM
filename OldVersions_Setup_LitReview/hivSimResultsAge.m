artBase = load('H:/HHCoM_Results/M_0.48_F_0.68.mat');
artAge = load('H:\HHCoM_Results\ART_gender_age.mat');

% tVec ,  popVec , newHiv ,...
%             newImage2Pv , newVaxHpv , newHpv , deaths , hivDeaths , ...
%             ccDeath , newCC , artTreatTracker , currYear , endYear , popLast

tVec = artBase.tVec;
gendLabel = {'Males' , 'Females'};
c = fix(clock);
currYear = c(1); % get the current year
yearNow = round((currYear - startYear) * stepsPerYear);
%% ART uptake
% base
figure()
artActual = [0	0	0	0	1	2	3	6, ...
    9	14	19	27	34	40	45	48];
yrsArtActual = [2000	2001	2002	2003	2004	2005 ...
    2006	2007	2008	2009	2010	2011	2012	2013 ...
    2014	2015];
for g = 1 : gender
    artInds1 = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 1 : 6 , 1 : risk));
    artInds2 = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 7 : age , 1 : risk));
    artPop1 = sum(artBase.popVec(: , artInds1) , 2);
    artPop2 = sum(artBase.popVec(: , artInds2) , 2);
    hivInds1 = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 1 : 6 , 1 : risk));
    hivInds2 = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 7 : age , 1 : risk));
    hivPop1 = sum(artBase.popVec(: , hivInds1) , 2);
    hivPop2 = sum(artBase.popVec(: , hivInds2) , 2);
    plot(tVec , 100 * artPop1 ./ (hivPop1 + artPop1) , ...
        tVec , 100 * artPop2 ./ (hivPop2 + artPop2))
    hold on
end
hold on
plot(yrsArtActual , artActual , '*')

xlabel('Year')
ylabel('Proportion of HIV Population (%)')
title('Proportion on ART (Base)')
legend('Male (<30)' , 'Male (>30)' , 'Female (<30)' , 'Female(>30)' , ...
    'General (Observed)')

figure()
artActual = [0	0	0	0	1	2	3	6, ...
    9	14	19	27	34	40	45	48];
yrsArtActual = [2000	2001	2002	2003	2004	2005 ...
    2006	2007	2008	2009	2010	2011	2012	2013 ...
    2014	2015];
for g = 1 : gender
    artInds1 = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 1 : 6 , 1 : risk));
    artInds2 = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 7 : age , 1 : risk));
    artPop1 = sum(artAge.popVec(: , artInds1) , 2);
    artPop2 = sum(artAge.popVec(: , artInds2) , 2);
    hivInds1 = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 1 : 6 , 1 : risk));
    hivInds2 = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 7 : age , 1 : risk));
    hivPop1 = sum(artAge.popVec(: , hivInds1) , 2);
    hivPop2 = sum(artAge.popVec(: , hivInds2) , 2);
    plot(tVec , 100 * artPop1 ./ (hivPop1 + artPop1) , ...
        tVec , 100 * artPop2 ./ (hivPop2 + artPop2))
    hold on
end
hold on
plot(yrsArtActual , artActual , '*')

xlabel('Year')
ylabel('Proportion of HIV Population (%)')
title('Proportion on ART (Differential Uptake by Age and Gender)')
legend('Male (<30)' , 'Male (>30)' , 'Female (<30)' , 'Female(>30)' , ...
    'General (Observed)')

%% Plot results
for g = 1 : gender
    % Total HIV positive
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hivInds1 = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 1 : 6 , 1 : risk));
    hivInds2 = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 7 : age , 1 : risk));
    hivPopBase = sum(artBase.popVec(: , hivInds) , 2);
    hivPopAge = sum(artAge.popVec(: , hivInds) , 2);
    hivPopAge1 = sum(artAge.popVec(: , hivInds1) , 2);
    hivPopAge2 = sum(artAge.popVec(: , hivInds2) , 2);
    
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    artInds1 = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 1 : 6 , 1 : risk));
    artInds2 = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 7 : age , 1 : risk));
    
    artPopBase = sum(artBase.popVec(: , artInds) , 2);
    artPopAge = sum(artAge.popVec(: , artInds) , 2);
    artPopAge1 = sum(artAge.popVec(: , artInds1) , 2);
    artPopAge2 = sum(artAge.popVec(: , artInds2) , 2);
    
    
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
        g , 4 : 10 , 1 : risk));
    totInds1 = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 1 : 6 , 1 : risk));
    totInds2 = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 7 : age , 1 : risk));
    
    basePopTot = artBase.popVec(: , totInds);
    agePopTot = artAge.popVec(: , totInds);
    agePopTot1 = artAge.popVec(: , totInds1);
    agePopTot2 = artAge.popVec(: , totInds2);

    basePrev = 100 * (hivPopBase + artPopBase) ./ sum(basePopTot , 2);
    agePrev = 100 * (hivPopAge + artPopAge) ./ sum(agePopTot , 2);
    agePrev1 = 100 * (hivPopAge1 + artPopAge1) ./ sum(agePopTot1 , 2);
    agePrev2 = 100 * (hivPopAge2 + artPopAge2) ./ sum(agePopTot2 , 2);
    
    plot(tVec , agePrev , tVec , agePrev1 , tVec , agePrev2 , tVec , basePrev)
%     hold on
    %     plot(overallHivPrev_KZN_AC(1 , :) , overallHivPrev_KZN_AC(2 , :) , '*')
    xlabel('Year'); ylabel('Proportion of Population (%)'); title(['HIV Prevalence (Ages 15-49) in ' , gendLabel{g}])
    xlim([tVec(1) , tVec(end)])
    legend('Age-varying uptake' , '<30' , '>=30' , 'Base');% , 'KZN Actual (Africa Center Data)')
    
    %% Change in HIV incidence
    susInds = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : 10 , 1 : risk));
    susInds1 = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 1 : 6 , 1 : risk));
    susInds2 = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 7 : age , 1 : risk));
    
    % population at risk (use midpoint)
    basePopSus = (sum(artBase.popVec(1 : end - 1 , susInds) , 2) + ...
        sum(artBase.popVec(2 : end , susInds) , 2)) ./ 2;
    agePopSus =(sum(artAge.popVec(1 : end - 1 , susInds) , 2) + ...
        sum(artAge.popVec(2 : end , susInds) , 2)) ./ 2;
    agePopSus1 =(sum(artAge.popVec(1 : end - 1 , susInds1) , 2) + ...
        sum(artAge.popVec(2 : end , susInds1) , 2)) ./ 2;
    agePopSus2 =(sum(artAge.popVec(1 : end - 1 , susInds2) , 2) + ...
        sum(artAge.popVec(2 : end , susInds2) , 2)) ./ 2;
    
    fac = 10 ^ 5;
    
    baseInc = sum(sum(sum(artBase.newHiv(2 : end , g , 4 : 10 , 1 : risk), 2), 3 ), 4)...
        ./ basePopSus * fac;
    ageInc =sum(sum(sum(artAge.newHiv(2 : end , g , 4 : 10 , 1 : risk), 2), 3 ), 4)...
        ./ agePopSus * fac;
    ageInc1 =sum(sum(sum(artAge.newHiv(2 : end , g , 1 : 6 , 1 : risk), 2), 3 ), 4)...
        ./ agePopSus1 * fac;
    ageInc2 =sum(sum(sum(artAge.newHiv(2 : end , g , 7 : age , 1 : risk), 2), 3 ), 4)...
        ./ agePopSus2 * fac;
    
    figure()
    plot(tVec(2 : end) , ageInc , tVec(2 : end) , ageInc1 , ...
        tVec(2 : end) , ageInc2 , tVec(2 : end) , baseInc);
    xlim([tVec(yearNow) , tVec(end)])
    xlabel('Year'); ylabel('Incidence per 100,000');
    title(['HIV Incidence in ' , gendLabel{g}])
    legend('General age varying uptake' , '<30' , ' >=30' , 'Base')
    
    deltaAge = (ageInc - baseInc) ./ baseInc * 100;
    deltaAge1 = (ageInc1 - baseInc) ./ baseInc * 100;
    deltaAge2 = (ageInc2 - baseInc) ./ baseInc * 100;
    
    figure()
    plot(tVec(2 : end) , deltaAge , tVec(2 : end) , deltaAge1 , ...
        tVec(2 : end) , deltaAge2);
    xlim([tVec(yearNow) , tVec(end)])
    xlabel('Year'); ylabel('Change (%)');
    title(['Change in Incidence in ' , gendLabel{g}])
    legend('General age varying uptake' , '<30' , ' >=30')
    % Change in HIV-related mortality
    fac = 10 ^ 5;
    
    baseMort = sum(sum(artBase.hivDeaths(2 : end , g , 1 : age), 2), 3) ./ basePopSus * fac;
    ageMort = sum(sum(artAge.hivDeaths(2 : end  , g , 1 : age) , 2), 3) ./ agePopSus * fac;
    ageMort1 = sum(sum(artAge.hivDeaths(2 : end  , g , 1 : 6) , 2), 3) ./ agePopSus1 * fac;
    ageMort2 = sum(sum(artAge.hivDeaths(2 : end  , g , 7 : age) , 2), 3) ./ agePopSus2 * fac;
    
    figure()
    plot(tVec(2 : end) , ageMort , tVec(2 : end) , ageMort1 , tVec(2 : end) , ageMort2);
    xlim([tVec(yearNow) , tVec(end)])
    xlabel('Year'); ylabel('Mortality per 100,000');
    title(['HIV Mortality in ' , gendLabel{g}])
    legend('General age varying uptake' , '<30' , ' >=30')

    deltaMortAge = (ageMort - baseMort) ./ baseMort * 100;
    deltaMortAge1 = (ageMort1 - baseMort) ./ baseMort * 100;
    deltaMortAge2 = (ageMort2 - baseMort) ./ baseMort * 100;
    
    figure()
    plot(tVec(2 : end) , deltaMortAge , tVec(2 : end) , deltaMortAge1 , ...
        tVec(2 : end) , deltaMortAge2);
    xlim([tVec(yearNow) , tVec(end)])
    xlabel('Year'); ylabel('Change (%)');
    title(['Change in Mortality in ' , gendLabel{g}])
    legend('General age varying uptake' , '<30' , ' >=30')

    %% summary table
    yr_2030 = (2030 - startYear) * stepsPerYear;
    yr_2040 = (2040 - startYear) * stepsPerYear;
    yr_2050 = (2050 - startYear) * stepsPerYear - 1;
    
    % Change in mortality
    deltaM_Age = [deltaMortAge(yr_2030) ; deltaMortAge(yr_2040) ; ...
        deltaMortAge(yr_2050)];
    deltaM_Age1 = [deltaMortAge1(yr_2030) ; deltaMortAge1(yr_2040) ; ...
        deltaMortAge1(yr_2050)];
    deltaM_Age2 = [deltaMortAge2(yr_2030) ; deltaMortAge2(yr_2040) ; ...
        deltaMortAge2(yr_2050)];
    
    % Change in incidence 
    deltaI_Age = [deltaAge(yr_2030) ; deltaAge(yr_2040) ; ...
        deltaAge(yr_2050)];
    deltaI_Age1 = [deltaAge1(yr_2030) ; deltaAge1(yr_2040) ; ...
        deltaAge1(yr_2050)];
    deltaI_Age2 = [deltaAge2(yr_2030) ; deltaAge2(yr_2040) ; ...
        deltaAge2(yr_2050)];
    
    baseP = [basePrev(yr_2030) ; basePrev(yr_2040) ; basePrev(yr_2050)];
    ageP = [agePrev(yr_2030) ; agePrev(yr_2040) ; agePrev(yr_2050)];
    age1P = [agePrev1(yr_2030) ; agePrev1(yr_2040) ; agePrev1(yr_2050)];
    age2P = [agePrev2(yr_2030) ; agePrev2(yr_2040) ; agePrev2(yr_2050)];
    
    hivSumm = table({'2030' ; '2040' ; '2050'} , baseP , ageP, age2P , age1P ,...
        deltaI_Age , deltaI_Age2 , deltaI_Age1 , deltaM_Age , deltaM_Age2 , deltaM_Age1);
    hivSumm.Properties.VariableNames{1} = 'Year';
    disp(['Impact summary for ' , gendLabel{g}])
    disp(hivSumm)
    filename = [gendLabel{g} , '_artAgeGendCompare.xlsx'];
    writetable(hivSumm , filename)
    disp(['Results for ' , gendLabel{g} , ' written to ' , filename])
end

%hivSimResultsGenAge()
