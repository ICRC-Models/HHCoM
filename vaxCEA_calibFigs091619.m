function vaxCEA_calibFigs091619()

%% LOAD PARAMETERS
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'],'stepsPerYear','circ','condUse','disease','viral',...
    'hpvTypes','hpvStates','periods','gender','age','risk','dim','k','toInd','sumall')
paramDir = [pwd , '\Params\'];
load([paramDir,'calibData'])

% Helper functions
sumall = @(x) sum(x(:));
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

% Time
c = fix(clock); % get time
currYear = c(1); % get the current year from time

%% LOAD SAVED RESULTS
dirName_calibCurr1 = 'toNow_090319calib_22Aug19_baseline';
dirName_calibCurr2 = 'toNow_090319calib_22Aug19_6_1';
dirName_calibCurr3 = 'toNow_090319calib_22Aug19_6_2';
dirName_calibCurr4 = 'toNow_090319calib_22Aug19_6_3';
dirName_calibCurr5 = 'toNow_090319calib_22Aug19_6_4';
dirName_calibCurr6 = 'toNow_090319calib_22Aug19_6_5';
dirName_calibCurr7 = 'toNow_090319calib_22Aug19_6_6';
dirName_calibCurr8 = 'toNow_090319calib_22Aug19_6_7';
dirName_calibCurr9 = 'toNow_090319calib_22Aug19_6_8';
dirName_calibCurr10 = 'toNow_090919calib_22Aug19_10_3607'; 
dirName_calibCurr11 = 'toNow_090919calib_22Aug19_11_4628';
 
currVec = {dirName_calibCurr1 , dirName_calibCurr2 , dirName_calibCurr3 , ...
    dirName_calibCurr4 , dirName_calibCurr5 , dirName_calibCurr6 , dirName_calibCurr7, ...
    dirName_calibCurr8, dirName_calibCurr9 , dirName_calibCurr10, dirName_calibCurr11};

for j = 1 : 1 %length(currVec)
    % Load results
    currModifier = currVec{j};
    load([pwd , '\HHCoM_Results\' , currModifier]); % ***SET ME***: name for historical run file
    
    % Plot settings
    reset(0)
    set(0 , 'defaultlinelinewidth' , 2)

    %% HIV status by age
%     ageGroup = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
%         '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
%         '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
%     linColor = {'k' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , ...
%         '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , ...
%         '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]'};
% 
%     % 2010-2016 AC data
%     hivPrevF_val = [9.29	9.02	10.45	9.33	10.37	11.00	9.35
%     31.41	31.68	30.64	33.95	34.56	34.12	33.42
%     53.27	51.72	50.80	51.33	51.94	53.98	52.41
%     59.18	61.35	58.66	64.90	62.57	64.71	63.09
%     53.97	54.08	58.77	65.12	65.28	64.66	66.95
%     42.69	43.27	45.29	49.16	54.25	56.37	61.28
%     32.34	34.30	39.18	41.47	48.21	49.57	50.23
%     ];
% 
%     hivPrevM_val = [1.60	1.85	2.75	3.46	2.87	3.95	4.50
%     9.56	8.02	9.87	9.65	11.86	7.19	8.02
%     28.99	21.92	24.88	29.84	35.40	27.65	27.31
%     46.47	44.51	39.49	47.22	46.35	41.64	42.08
%     52.03	44.30	49.61	63.33	51.41	52.05	51.35
%     41.73	41.53	51.55	51.64	59.40	52.69	51.18
%     36.64	37.12	33.01	40.00	40.54	44.52	52.17
%     ];
% 
%     hivM = hivPrevM_obs(: , 2) ./ hivPrevM_obs(: , 3) .* 100;
%     hivF = hivPrevF_obs(: , 2) ./hivPrevF_obs(: , 3) .* 100;
%     prevYears = unique(hivPrevF_obs(: , 1));
% %     prevYears2 = [2010 : 2016];
%     gen = {'Male' , 'Female'};
%     for g = 1 : gender
%         hivPrevs = hivM;
% %         hivPrevs2 = hivPrevM_val;
%         if g == 2
%             hivPrevs = hivF;
% %             hivPrevs2 = hivPrevF_val;
%         end
%         for a = 4 : 15
%             hivAgeInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%                 g , a , 1 : risk)); toInd(allcomb(10, 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%                 g , a , 1 : risk))];
%             ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%                 g , a , 1 : risk));
%             hivAge(: , a) = sum(popVec(: , hivAgeInds) , 2);
%             hivAgeRel = bsxfun(@rdivide , hivAge(: , a) , sum(popVec(: , ageInds) , 2)) * 100;
%             subplot(4 , 3 , a - 3)
%             if a <= 10
%                 plot(tVec , hivAgeRel , 'Color' , linColor{j})
%                 hold all;
%                 plot(prevYears , hivPrevs((a - 3) : 7 : end) ,'ro'); % , ...
%                  %   prevYears2 , hivPrevs2((a - 3) : 7 : end) , 'bo');
%                 hold all;
%             else
%                 plot(tVec , hivAgeRel , 'Color' , linColor{j});
%                 hold all;
%             end
%             xlabel('Year'); ylabel('Prevalence (%)'); title([gen{g} , 's ages ' , ageGroup{a}]) % , ' HIV Prevalence'])
%             xlim([1980 2015])
%         end
%         legend('Model' , 'Africa Center Data (Calibration)'); % , 'Africa Center Data (Validation)')
%     end
%     hold all;

     %% General HPV Prevalence by Age in 2002
%     ageGroup = {'15 - 19' , '20 -24' , '25 - 29' ,...
%         '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
%         '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
%     linColor = {'k' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , ...
%         '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , ...
%         '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]'};
%     hpv2017 = zeros(age - 4 + 1 , 1);
%     hpvHIV2017 = hpv2017;
%     hpvNeg2017 = hpv2017;
% 
%     for a = 4 : age
%         hrInds = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : 4, ...
%             1 : periods , 2 , a , 1 : risk));
% 
%         ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%             2 , a , 1 : risk));
%         hpv2017(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , hrInds)))...
%             ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
% 
%         % HIV+
%         hrHIVInds = [toInd(allcomb(2 : 6 , 1 : viral , 2 : 4 , 1 : 4, ...
%             1 : periods , 2 , a , 1 : risk));toInd(allcomb(10 , 6 , 2 : 4 , 1 : 4, ...
%             1 : periods , 2 , a , 1 : risk));];
% 
%         ageHIVInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%             2 , a , 1 : risk)); toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
%             1 : periods , 2 , a , 1 : risk))];
%         hpvHIV2017(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , hrHIVInds)))...
%             ./ sum(popVec((2002 - startYear) * stepsPerYear , ageHIVInds)) * 100;
% 
%         % HIV-
%         hrNegInds = toInd(allcomb(1 , 1 , 2 : 4 , 1 : 4, ...
%             1 : periods , 2 , a , 1 : risk));
% 
%         ageNegInds = toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%             2 , a , 1 : risk));
%         hpvNeg2017(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , hrNegInds)))...
%             ./ sum(popVec((2002 - startYear) * stepsPerYear , ageNegInds)) * 100;
%     end
% 
%     % McDonald 2014
%     hpvHivObs(: , 1) = [0.75
%     0.61
%     0.60
%     0.55
%     0.46
%     0.42
%     0.43
%     0.54
%     0.35];
% 
%     hpvHivObs(: , 2) = [0.63
%     0.54
%     0.54
%     0.47
%     0.42
%     0.34
%     0.32
%     0.35
%     0.16];
% 
%     hpvHivObs(: ,3) = [0.87
%     0.67
%     0.66
%     0.62
%     0.51
%     0.50
%     0.55
%     0.72
%     0.53];
% 
%     hpvNegObs(: , 1) = [0.60
%     0.38
%     0.24
%     0.20
%     0.19
%     0.18
%     0.13
%     0.17
%     0.15];
% 
%     hpvNegObs(: , 2) = [0.53
%     0.34
%     0.21
%     0.17
%     0.18
%     0.16
%     0.11
%     0.14
%     0.12];
% 
%     hpvNegObs(: , 3) = [0.67
%     0.41
%     0.27
%     0.23
%     0.21
%     0.20
%     0.15
%     0.19
%     0.18];
% 
%     hpvHivObs = hpvHivObs * 100;
%     hpvNegObs = hpvNegObs * 100;
%     plot(1 : length(hpv2017) , hpv2017 ,'o-' , 'Color' , linColor{j})
%     hold on
% %     plot(1 : length(hpvHIV2017) , hpvHIV2017 , 'bo-');
% %     hold on
% %     plot(1 : length(hpvNeg2017) , hpvNeg2017 , 'r-')
% %     set(gca , 'xtickLabel' , ageGroup);
% 
%     % general
%     % yPosError = abs(hrHpvObs(: , 3) - hrHpvObs(: , 1));
%     % yNegError = abs(hrHpvObs(: , 2) - hrHpvObs(: , 1));
%     % errorbar(1 : length(hrHpvObs) , hrHpvObs(: , 1) , yNegError , yPosError , 'rs')
%     % HIV+
%     yPosError = abs(hpvHivObs(: , 3) - hpvHivObs(: , 1));
%     yNegError = abs(hpvHivObs(: , 2) - hpvHivObs(: , 1));
%     errorbar(1 : length(hpvHivObs) , hpvHivObs(: , 1) , yNegError , yPosError , 'bs')
%     %HIV-
%     hold on
%     yPosError = abs(hpvNegObs(: , 3) - hpvNegObs(: , 1));
%     yNegError = abs(hpvNegObs(: , 2) - hpvNegObs(: , 1));
%     errorbar(1 : length(hpvNegObs) , hpvNegObs(: , 1) , yNegError , yPosError , 'rs')
% 
%     set(gca , 'xtick' , 1 : length(hpvNegObs) , 'xtickLabel' , ageGroup);
%     legend('General' , 'HIV+' , 'HIV-' , 'McDonald 2014 - HIV+' , 'McDonald 2014 - HIV-')
%     xlabel('Age Group'); ylabel('Prevalence (%)')
%     title('Age Specific hrHPV Prevalence in 2017')
%     hold all;


%% CIN2/3 Prevalence All HR HPV combined
cinPos2017 = zeros(age - 4 + 1 , 1);
cinNeg2017 = cinPos2017;
ageGroup = {'15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};
linColor = {'k' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , ...
        '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , ...
        '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]' , '[0.5, 0.5, 0.5]'};
for a = 4 : age
    % HIV+
    cinInds = [toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvTypes , 3 : 4, ...
        1 : periods , 2 , a , 1 : risk));toInd(allcomb(10 , 6 , 2 : hpvTypes , 3 : 4, ...
        1 : periods , 2 , a , 1 : risk))];
    ageInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk))];
    cinPos2017(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , cinInds)))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
    
    cinNegInds = [toInd(allcomb(1, 1 : viral , 2 : hpvTypes , 3 : 4, ...
        1 : periods , 2 , a , 1 : risk));...
        toInd(allcomb(7 : disease , 1 : 5 , 2 : hpvTypes , 3 : 4, ...
        1 : periods , 2 , a , 1 : risk))];
    ageNegInds = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk));...
        toInd(allcomb(7 : disease , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        2 , a , 1 : risk))];
    cinNeg2017(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , cinNegInds)))...
        ./ (sum(popVec((2002 - startYear) * stepsPerYear , ageNegInds))) * 100;
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

subplot(2 , 1 , 1)
cinPosAct = cinPosAct .* 100; % convert to %
yPosError = abs(cinPosAct(: , 3) - cinPosAct(: , 1));
yNegError = abs(cinPosAct(: , 2) - cinPosAct(: , 1));
plot(1 : length(cinPos2017) , cinPos2017 ,'o-' , 'Color' , linColor{j})
hold on
errorbar(1 : length(cinPosAct) , cinPosAct(: , 1) , yNegError , yPosError , 'rs')
legend('HR HPV CIN 2/3' , 'McDonald 2014')
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('Age Specific CIN 2/3 Prevalence Among HIV+ in 2017')
ylim([0 25])

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
plot(1 : length(cinNeg2017) , cinNeg2017 , 'o-' , 'Color' , linColor{j})
hold on
yPosError = abs(cinNegAct(: , 3) - cinNegAct(: , 1));
yNegError = abs(cinNegAct(: , 2) - cinNegAct(: , 1));
errorbar(1 : length(cinNegAct) , cinNegAct(: , 1) , yNegError , yPosError , 'rs')
legend('HR HPV CIN 2/3' , 'McDonald 2014')
set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Prevalence (%)')
title('Age Specific CIN 2/3 Prevalence Among HIV- in 2017')
ylim([0 25])
    
end

