all90 = load('H:\HHCoM_Results\all90.mat');
art70_90 = load('H:\HHCoM_Results\art70_90.mat');
art70_95 = load('H:\HHCoM_Results\art70_95.mat');
load('general')
load('settings')
c = fix(clock);
currYear = c(1); % get the current year
yearNow = round((currYear - startYear) * stepsPerYear);

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
%% Age-aggregated Disease Incidence
wVec = zeros(age , 1);
wVec(5 : age) = [0.188 , 0.18 , 0.159 , 0.121 , 0.088 , 0.067 , 0.054 , ...
    0.046 , 0.038 , 0.029 , 0.017 , 0.013]; 
figure()
newHiv_Arr = {all90.newHiv , art70_90.newHiv , art70_95.newHiv};
popVec_Arr = {all90.popVec , art70_90.popVec , art70_95.popVec};
incMat = zeros(age , size(all90.popVec , 1) / stepsPerYear);

inc = {incMat , incMat , incMat , incMat};
tVec = all90.tVec;
for i = 1 : length(newHiv_Arr)
    for a = 4 : age
        newHiv = sum(sum(sum(newHiv_Arr{i}(1 : end , 1 : gender , a , :)...
            ,2),3),4);
        popVec = popVec_Arr{i};
        hivSusInds = toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 1 : gender , a , 1 : risk));
        hivSus = sum(popVec(1 : end , hivSusInds) , 2);
        hivSus_Year = sum(reshape(hivSus , stepsPerYear , size(hivSus , 1) ...
            / stepsPerYear)) ./ stepsPerYear; % average susceptible population size per year
        newHiv_Year = sum(reshape(newHiv , stepsPerYear , size(newHiv , 1) ...
            /stepsPerYear)); % total new HIV infections per year
        inc{i}(a , :) = newHiv_Year ./ hivSus_Year .* 100;
    end
end

for i = 1 : length(newHiv_Arr)
    incAS = sum(bsxfun(@times , inc{i} , wVec));
    plot(tVec(1 : stepsPerYear : end) , incAS)
    xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
    hold on
end
xlabel('Year'); ylabel('Incidence per 100'); title('HIV Incidence')
% Women aged 16-29:  35% on treatment (30% suppressed)
% Women aged 30+: 60% on treatment (55% suppressed)
% Men aged 16-29: 25% on treatment (20% suppressed)
% Men aged 30+: 50% on treatment (40% suppressed).
legend('Base (90% all)' , '70/90' , '70/95', ...
    'Location' , 'northeastoutside')
%% Difference in Age-standardized Disease Incidence
% relative to base case, i.e. 90% coverage for all
figure()
for i = 2 : length(newHiv_Arr)
    incAS = sum(bsxfun(@times , inc{i} , wVec)) ...
        - sum(bsxfun(@times , inc{1} , wVec));
    plot(tVec(1 : stepsPerYear : end) , incAS)
    xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
    hold on
end
xlabel('Year'); ylabel('Difference in Incidence per 100'); 
title('Difference in HIV Incidence (vs Base Case)')
% Women aged 16-29:  35% on treatment (30% suppressed)
% Women aged 30+: 60% on treatment (55% suppressed)
% Men aged 16-29: 25% on treatment (20% suppressed)
% Men aged 30+: 50% on treatment (40% suppressed).
legend('70/90' , '70/95', ...
    'Location' , 'northeastoutside')

% percentage form
figure()
for i = 2 : length(newHiv_Arr)
    incAS = (sum(bsxfun(@times , inc{i} , wVec)) ...
        - sum(bsxfun(@times , inc{1} , wVec))) ...
        ./ sum(bsxfun(@times , inc{1} , wVec)) * 100;
    plot(tVec(1 : stepsPerYear : end) , incAS)
    xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
    hold on
end
xlabel('Year'); ylabel('Difference in Incidence (%)'); 
title('Difference in HIV Incidence (vs Base Case)')
% Women aged 16-29:  35% on treatment (30% suppressed)
% Women aged 30+: 60% on treatment (55% suppressed)
% Men aged 16-29: 25% on treatment (20% suppressed)
% Men aged 30+: 50% on treatment (40% suppressed).
legend('70/90' , '70/95', ...
    'Location' , 'northeastoutside')

%% Non Age-standardized Disease Incidence
figure()
inc = {incMat , incMat , incMat , incMat};
tVec = all90.tVec;
for i = 1 : length(newHiv_Arr)
    newHiv = sum(sum(sum(newHiv_Arr{i}(1 : end , 1 : gender , 4 : age , :)...
        ,2),3),4);
    popVec = popVec_Arr{i};
    hivSusInds = toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 4 : age , 1 : risk));
    hivSus = sum(popVec(1 : end , hivSusInds) , 2);
    hivSus_Year = sum(reshape(hivSus , stepsPerYear , size(hivSus , 1) ...
        / stepsPerYear)) ./ stepsPerYear; % average susceptible population size per year
    newHiv_Year = sum(reshape(newHiv , stepsPerYear , size(newHiv , 1) ...
        /stepsPerYear)); % total new HIV infections per year
    inc{i} = newHiv_Year ./ hivSus_Year .* 100;
end

for i = 1 : length(newHiv_Arr)
    incAS = sum(bsxfun(@times , inc{i} , wVec));
    plot(tVec(1 : stepsPerYear : end) , incAS)
    xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
    hold on
end
xlabel('Year'); ylabel('Non-Age-Standardized Incidence per 100'); title('HIV Incidence')
% Women aged 16-29:  35% on treatment (30% suppressed)
% Women aged 30+: 60% on treatment (55% suppressed)
% Men aged 16-29: 25% on treatment (20% suppressed)
% Men aged 30+: 50% on treatment (40% suppressed).
legend('Base (90% all)' , '70/90' , '70/95', ...
    'Location' , 'northeastoutside')

%% ART uptake
ageArr = {4 : 6 , 7 : age};
figure()
for i = 1 : length(ageArr)
    hivInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , 1 : 2 , ageArr{i} , 1 : risk)); toInd(allcomb(10 , 6 , ...
        1 : hpvTypes , 1 : hpvStates, 1 : periods , 1 : 2 , ageArr{i} , 1 : risk))];
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , 1 : 2 , ageArr{i} , 1 : risk));
    plot(tVec , sum(all90.popVec(: , artInds) , 2) ./ sum(all90.popVec(: , hivInds) , 2) * 100 ,...
        tVec , sum(art70_90.popVec(: , artInds) , 2) ./ sum(art70_90.popVec(: , hivInds) , 2) * 100 , ...
        tVec , sum(art70_95.popVec(: , artInds) , 2) ./ sum(art70_95.popVec(: , hivInds) , 2) * 100)
    hold on
end
hold off
legend('90/90 (<30)' , '70/90 (<30)' , '70/95 (<30)', ...
    '90/90 (30+)' , '70/90 (30+)' , '70/95 (30+)', ...
    'Location' , 'northeastoutside')
title('Overall ART Coverage')
xlabel('Year'); ylabel('Coverage (%)')

%% HIV prevalence
hivInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk)); toInd(allcomb(10 , 6 , ...
    1 : hpvTypes , 1 : hpvStates, 1 : periods , 1 : 2 , 4 : 10 , 1 : risk))];
allInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk))];
artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk));
figure()
plot(tVec , sum(all90.popVec(: , hivInds) , 2) ./ sum(all90.popVec(: , allInds) , 2) * 100 ,...
    tVec , sum(art70_90.popVec(: , hivInds) , 2) ./ sum(art70_90.popVec(: , allInds) , 2) * 100 , ...
    tVec , sum(art70_95.popVec(: , hivInds) , 2) ./ sum(art70_95.popVec(: , allInds) , 2) * 100)
legend('Base (90% all)' , '70/90' , '70/95', ...
    'Location' , 'northeastoutside')
title('HIV Prevalence')
xlabel('Year'); ylabel('Prevalence (%)')

%% Circumcised
circInds = [toInd(allcomb(7, 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 , 4 : 10 , 1 : risk))];
allMInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 , 4 : 10 , 1 : risk))];
figure()
plot(tVec , sum(all90.popVec(: , circInds) , 2) ./ sum(all90.popVec(: , allMInds) , 2) * 100 ,...
    tVec , sum(art70_90.popVec(: , circInds) , 2) ./ sum(art70_90.popVec(: , allMInds) , 2) * 100 , ...
    tVec , sum(art70_95.popVec(: , circInds) , 2) ./ sum(art70_95.popVec(: , allMInds) , 2) * 100)
legend('Base (90% all)' , '70/90' , '70/95', ...
    'Location' , 'northeastoutside')
title('Circumcision')
xlabel('Year'); ylabel('Coverage (%)')