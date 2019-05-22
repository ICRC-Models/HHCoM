artAge = load('H:\HHCoM_Results\ART_gender_age.mat');
artMed = load('H:/HHCoM_Results/M_0.7_F_0.68.mat');
artHigh = load('H:/HHCoM_Results/M_0.85_F_0.83.mat');
artMH = load('H:/HHCoM_Results/M_0.7_F_0.83.mat');

% tVec ,  popVec , newHiv ,...
%             newImmHpv , newVaxHpv , newHpv , deaths , hivDeaths , ...
%             ccDeath , newCC , artTreatTracker , currYear , endYear , popLast

tVec = artAge.tVec;
gendLabel = {'Males' , 'Females'};
%% Plot results
for g = 1 : gender
    % Total HIV positive
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hivPopAge = sum(artAge.popVec(: , hivInds) , 2);
    hivPopMed = sum(artMed.popVec(: , hivInds) , 2);
    hivPopHigh = sum(artHigh.popVec(: , hivInds) , 2);
    hivPopMH = sum(artMH.popVec(: , hivInds) , 2);
    
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    artPopAge = sum(artAge.popVec(: , artInds) , 2);
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
        g , 4 : 10, 1 : risk));
    agePopTot = artAge.popVec(: , totInds);
    medPopTot = artMed.popVec(: , totInds);
    highPopTot = artHigh.popVec(: , totInds);
    mhPopTot = artMH.popVec(: , totInds);
    
    agePrev = 100 * (hivPopAge + artPopAge) ./ sum(agePopTot , 2);
    medPrev = 100 * (hivPopMed + artPopMed) ./ sum(medPopTot , 2);
    highPrev = 100 * (hivPopHigh + artPopHigh) ./ sum(highPopTot , 2);
    mhPrev = 100 * (hivPopMH + artPopMH) ./ sum(mhPopTot , 2);
    
    plot(tVec , medPrev  , tVec , mhPrev , tVec , highPrev , tVec , agePrev)
    hold on
    %     plot(overallHivPrev_KZN_AC(1 , :) , overallHivPrev_KZN_AC(2 , :) , '*')
    xlabel('Year'); ylabel('Proportion of Population (%)'); title(['HIV Prevalence (Ages 15-49) in ' , gendLabel{g}])
    xlim([tVec(1) , tVec(end)])
    legend('65% Both' , 'M:65%, F:80%' , '80% Both' , 'Gender & Age Dependent');% , 'KZN Actual (Africa Center Data)')
    
    %% Change in HIV incidence
    susInds = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : 10, 1 : risk));
    
    % population at risk (use midpoint)
    agePopSus = (sum(artAge.popVec(1 : end - 1 , susInds) , 2) + ...
        sum(artAge.popVec(2 : end , susInds) , 2)) ./ 2;
    medPopSus =(sum(artMed.popVec(1 : end - 1 , susInds) , 2) + ...
        sum(artMed.popVec(2 : end , susInds) , 2)) ./ 2;
    highPopSus = (sum(artHigh.popVec(1 : end - 1 , susInds) , 2) + ...
        sum(artHigh.popVec(2 : end , susInds) , 2)) ./ 2;
    mhPopSus = (sum(artMH.popVec(1 : end - 1 , susInds) , 2) + ...
        sum(artHigh.popVec(2 : end , susInds) , 2)) ./ 2;
    
    fac = 10 ^ 5;
    
    ageInc = sum(sum(sum(artAge.newHiv(2 : end , g , 4 : 10 , 1 : risk), 2), 3 ), 4)...
        ./ agePopSus * fac;
    medInc =sum(sum(sum(artMed.newHiv(2 : end , g , 4 : 10 , 1 : risk), 2), 3 ), 4)...
        ./ medPopSus * fac;
    highInc = sum(sum(sum(artHigh.newHiv(2 : end , g , 4 : 10 , 1 : risk), 2), 3 ), 4)...
        ./ highPopSus * fac;
    mhInc = sum(sum(sum(artMH.newHiv(2 : end , g , 4 : 10 , 1 : risk), 2), 3 ), 4)...
        ./ highPopSus * fac;
    
    figure()
    plot(tVec(2 : end) , medInc , tVec(2 : end ) , mhInc , ...
        tVec( 2 : end) , highInc , tVec(2 : end) , ageInc);
    xlim([tVec(1) , tVec(end)])
    xlabel('Year'); ylabel('Incidence per 100,000');
    title(['HIV Incidence in ' , gendLabel{g}])
    legend('65% Both' , 'M:65%, F:80%' , '80% Both' , 'Gender & Age Dependent')
    
    deltaMed = (medInc - ageInc) ./ ageInc * 100;
    deltaHigh = (highInc - ageInc) ./ ageInc * 100;
    deltaMH = (mhInc - ageInc) ./ ageInc * 100;
    
    figure()
    plot(tVec(2 : end) , deltaMed , tVec(2 : end) , deltaMH , ...
        tVec(2 : end) , deltaHigh);
    xlim([tVec(1) , tVec(end)])
    xlabel('Year'); ylabel('Change (%)');
    title(['Change in Incidence in ' , gendLabel{g}])
    legend('65% Both' , 'M:65%, F:80%' , '80% Both')
    
    % Change in HIV-related mortality
    fac = 10 ^ 5;
    
    ageMort = sum(sum(artAge.hivDeaths(2 : end , g , 1 : age), 2), 3) ./ agePopSus * fac;
    medMort = sum(sum(artMed.hivDeaths(2 : end  , g , 1 : age) , 2), 3) ./ medPopSus * fac;
    highMort = sum(sum(artHigh.hivDeaths(2 : end , g , 1 : age) , 2), 3) ./ highPopSus * fac;
    mhMort = sum(sum(artMH.hivDeaths(2 : end , g , 1 : age) , 2), 3) ./ mhPopSus * fac;
    
    figure()
    plot(tVec(2 : end) , medMort , tVec(2 : end) , mhMort , ...
        tVec( 2 : end) , highMort , tVec(2 : end) , ageMort);
    xlim([tVec(1) , tVec(end)])
    xlabel('Year'); ylabel('Mortality per 100,000');
    title(['HIV Mortality in ' , gendLabel{g}])
    legend('65% Both' , 'M:65%, F:80%' , '80% Both' , 'Gender & Age Dependent')
    
    deltaMortMed = (medMort - ageMort) ./ ageMort * 100;
    deltaMortHigh = (highMort - ageMort) ./ ageMort * 100;
    deltaMortMH = (mhMort - ageMort) ./ ageMort * 100;
    
    figure()
    plot(tVec(2 : end) , deltaMortMed , tVec(2 : end) , deltaMortMH , ...
        tVec(2 : end) , deltaMortHigh);
    xlim([tVec(1) , tVec(end)])
    xlabel('Year'); ylabel('Change (%)');
    title(['Change in Mortality in ' , gendLabel{g}])
    legend('65% Both' , 'M:65%, F:80%' , '80% Both')
    
    %% summary table
    yr_2030 = (2030 - startYear) * stepsPerYear;
    yr_2040 = (2040 - startYear) * stepsPerYear;
    yr_2050 = (2050 - startYear) * stepsPerYear - 1;
    
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
    
    baseP = [agePrev(yr_2030) ; agePrev(yr_2040) ; agePrev(yr_2050)];
    medP = [medPrev(yr_2030) ; medPrev(yr_2040) ; medPrev(yr_2050)];
    highP = [highPrev(yr_2030) ; highPrev(yr_2040) ; highPrev(yr_2050)];
    mhP = [mhPrev(yr_2030) ; mhPrev(yr_2040) ; mhPrev(yr_2050)];
    
%     hivSumm = table({'2030' ; '2040' ; '2050'} , baseP , medP, mhP , highP ,...
%         deltaIMed , deltaIMH , deltaIHigh , deltaM_Med , deltaM_MH , deltaM_High);
%     hivSumm.Properties.VariableNames{1} = 'Year';
%     disp(['Impact summary for ' , gendLabel{g}])
%     disp(hivSumm)
%     filename = [gendLabel{g} , '_artCompare.xlsx'];
%     writetable(hivSumm , filename)
%     disp(['Results for ' , gendLabel{g} , ' written to ' , filename])
end

% hivSimResultsGen()
