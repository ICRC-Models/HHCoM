noDap = load('H:\HHCoM_Results\noDap.mat');
dap10 = load('H:\HHCoM_Results\dap10.mat');
dap20 = load('H:\HHCoM_Results\dap20.mat');
paramDir = [pwd , '\Params\'];
load([paramDir , 'general'])
load([paramDir , 'settings'])
c = fix(clock);
currYear = c(1); % get the current year
yearNow = round((currYear - startYear) * stepsPerYear);

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
set(0, 'DefaultAxesFontSize',16)
%% Age-standardized Disease Incidence
wVec = zeros(age , 1);
wVec(5 : age) = [0.188 , 0.18 , 0.159 , 0.121 , 0.088 , 0.067 , 0.054 , ...
    0.046 , 0.038 , 0.029 , 0.017 , 0.013]; 
figure()
% tVec = noDap.tVec;
tVec = linspace(startYear , endYear , size(noDap.popVec , 1));
newHiv_Arr = {noDap.newHiv , dap10.newHiv , dap20.newHiv};
popVec_Arr = {noDap.popVec , dap10.popVec , dap20.popVec};
incMat = zeros(age , length(tVec) / stepsPerYear);

inc = {incMat , incMat , incMat , incMat};
for i = 1 : length(newHiv_Arr)
    for a = 4 : age
        newHiv = sum(sum(sum(newHiv_Arr{i}(1 : end , 1 : gender , a , :)...
            ,2),3),4);
        popVec = popVec_Arr{i};
        hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , a , 1 : risk)); ...
        toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , a , 1 : risk))];
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
legend('Base' , '10% Dapivirine coverage' , ...
    '20% Dapivirine coverage' , ... 
    'Location' , 'northeastoutside')
%% Difference in Age-standardized Disease Incidence
% relative to base case, i.e. 70% coverage for all
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
legend('10% Dapivirine coverage' , ...
    '20% Dapivirine coverage' , ... 
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
legend('10% Dapivirine coverage' , ...
    '20% Dapivirine coverage' , ... 
    'Location' , 'northeastoutside')


%% Non Age-standardized Disease Incidence
inc = {incMat , incMat , incMat , incMat};
newHiv_Year = zeros(length(newHiv_Arr) , length(tVec) / stepsPerYear);
for i = 1 : length(newHiv_Arr)
    newHiv = sum(sum(sum(newHiv_Arr{i}(1 : end , 1 : gender , 4 : age , :)...
        ,2),3),4);
    popVec = popVec_Arr{i};
    hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 4 : age , 1 : risk)); ...
        toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 4 : age , 1 : risk))];
    hivSus = sum(popVec(1 : end , hivSusInds) , 2);
    hivSus_Year = sum(reshape(hivSus , stepsPerYear , size(hivSus , 1) ...
        / stepsPerYear)) ./ stepsPerYear; % average susceptible population size per year
    newHiv_Year(i , :) = sum(reshape(newHiv , stepsPerYear , size(newHiv , 1) ...
        /stepsPerYear)); % total new HIV infections per year
    inc{i} = newHiv_Year(i , :) ./ hivSus_Year .* 100;
end

figure()
for i = 1 : length(newHiv_Arr)
    plot(tVec(1 : stepsPerYear : end) , inc{i})
    xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
    hold on
end
xlabel('Year'); ylabel('Non-Age-Standardized Incidence per 100'); title('HIV Incidence')
legend('Base' , '10% Dapivirine coverage' , ...
    '20% Dapivirine coverage' , ... 
    'Location' , 'northeastoutside')


figure()
for i = 1 : length(newHiv_Arr)
    plot(tVec(1 : stepsPerYear : end) , newHiv_Year(i , :))
    xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
    hold on
end
xlabel('Year'); ylabel('New HIV Cases'); title('New HIV Cases Per Year')
legend('No Dapivirine' , '10% Dapivirine coverage' , ...
    '20% Dapivirine coverage' , ... 
    'Location' , 'northeastoutside')



%% Overall ART uptake
figure()

hivInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : age , 1 : risk)); toInd(allcomb(10 , 6 , ...
    1 : hpvTypes , 1 : hpvStates, 1 : periods , 1 : 2 , 4 : age , 1 : risk))];
artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : age , 1 : risk));
plot(tVec , sum(noDap.popVec(: , artInds) , 2) ./ sum(noDap.popVec(: , hivInds) , 2) * 100 ,...
    tVec , sum(dap10.popVec(: , artInds) , 2) ./ sum(dap10.popVec(: , hivInds) , 2) * 100 ,...
    tVec , sum(dap20.popVec(: , artInds) , 2) ./ sum(dap20.popVec(: , hivInds) , 2) * 100)

legend('No Dapivirine' , '10% Dapivirine' , '20% Dapivirine' , ...
    'Location' , 'northeastoutside')
title('Overall ART Coverage')
xlabel('Year'); ylabel('Coverage (%)')

%%
figure()
plot(tVec , sum(noDap.popVec(: , artInds) , 2) ./ sum(noDap.popVec(: , hivInds) , 2) * 100) 
title('Overall ART Coverage')
xlabel('Year'); ylabel('Coverage (%)')

%% ART uptake by age
% ageArr = {4 : 6 , 7 : age};
% figure()
% for i = 1 : length(ageArr)
%     hivInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
%         1 : periods , 1 : 2 , ageArr{i} , 1 : risk)); toInd(allcomb(10 , 6 , ...
%         1 : hpvTypes , 1 : hpvStates, 1 : periods , 1 : 2 , ageArr{i} , 1 : risk))];
%     artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
%         1 : periods , 1 : 2 , ageArr{i} , 1 : risk));
%     plot(tVec , sum(noDap.popVec(: , artInds) , 2) ./ sum(noDap.popVec(: , hivInds) , 2) * 100 ,...
%         tVec , sum(dap10.popVec(: , artInds) , 2) ./ sum(dap10.popVec(: , hivInds) , 2) * 100 ,...
%         tVec , sum(dap20.popVec(: , artInds) , 2) ./ sum(dap20.popVec(: , hivInds) , 2) * 100)
%     hold on
% end
% hold off
% legend('No Dapivirine (<30)' , '10% Dapivirine (<30)' , '20% Dapivirine (<30)' , ...
%     'No Dapivirine (30+)' , '10% Dapivirine (30+)' , '20% Dapivirine (30+)' ,  ...
%     'Location' , 'northeastoutside')
% title('Overall ART Coverage')
% xlabel('Year'); ylabel('Coverage (%)')

%% HIV prevalence
hivInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk)); toInd(allcomb(10 , 6 , ...
    1 : hpvTypes , 1 : hpvStates, 1 : periods , 1 : 2 , 4 : 10 , 1 : risk))];
allInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk))];
artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : 10 , 1 : risk));
figure()
plot(tVec , sum(noDap.popVec(: , hivInds) , 2) ./ sum(noDap.popVec(: , allInds) , 2) * 100 ,...
    tVec , sum(dap10.popVec(: , hivInds) , 2) ./ sum(dap10.popVec(: , allInds) , 2) * 100 ,...
    tVec , sum(dap20.popVec(: , hivInds) , 2) ./ sum(dap20.popVec(: , allInds) , 2) * 100)
legend('Base' , '10% Dapivirine coverage' , ...
    '20% Dapivirine coverage' , ... 
    'Location' , 'northeastoutside')
title('HIV Prevalence')
xlabel('Year'); ylabel('Prevalence (%)')

%% Circumcised
% circInds = [toInd(allcomb(7, 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
%     1 : periods , 1 , 4 : 10 , 1 : risk))];
% allMInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
%     1 : periods , 1 , 4 : 10 , 1 : risk))];
% figure()
% plot(tVec , sum(all70.popVec(: , circInds) , 2) ./ sum(all70.popVec(: , allMInds) , 2) * 100 ,...
%     tVec , sum(all90.popVec(: , circInds) , 2) ./ sum(all90.popVec(: , allMInds) , 2) * 100 ,...
%     tVec , sum(art70_90.popVec(: , circInds) , 2) ./ sum(art70_90.popVec(: , allMInds) , 2) * 100 , ...
%     tVec , sum(art70_95.popVec(: , circInds) , 2) ./ sum(art70_95.popVec(: , allMInds) , 2) * 100)
% legend('70/70' , '90/90' , '70/90' , '70/95', ...
%     'Location' , 'northeastoutside')
% title('Circumcision')
% xlabel('Year'); ylabel('Coverage (%)')

%% Dapivirine
dapInds = [toInd(allcomb(7 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 2 , 4 : 10 , 1 : risk))];
allFInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 2 , 4 : 10 , 1 : risk))];
figure()
plot(tVec , sum(noDap.popVec(: , dapInds) , 2) ./ sum(noDap.popVec(: , allFInds) , 2) * 100 ,...
    tVec , sum(dap10.popVec(: , dapInds) , 2) ./ sum(dap10.popVec(: , allFInds) , 2) * 100 ,...
    tVec , sum(dap20.popVec(: , dapInds) , 2) ./ sum(dap20.popVec(: , allFInds) , 2) * 100)
legend('No Dapivirine' , '10% Dapivirine coverage' , ...
    '20% Dapivirine coverage' , ... 
    'Location' , 'northeastoutside')
title('Dapivirine Coverage')
xlabel('Year'); ylabel('Coverage (%)')