all90 = load('H:\HHCoM_Results\all90.mat');
art70_90 = load('H:\HHCoM_Results\art70_90.mat');
art70_95 = load('H:\HHCoM_Results\art70_95.mat');
all70 = load('H:\HHCoM_Results\all70.mat');
paramDir = [pwd , '\Params\'];
load([paramDir,'general'])
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
set(0 , 'DefaultAxesFontSize' , 16)
newHiv_Arr = {all70.newHiv , art70_90.newHiv , art70_95.newHiv , all90.newHiv};
popVec_Arr = {all70.popVec , art70_90.popVec , art70_95.popVec , all90.popVec};
incMat = zeros(age , size(all90.popVec , 1) / stepsPerYear);

%% Age-standardized Disease Incidence
% wVec = zeros(age , 1);
% wVec(5 : age) = [0.188 , 0.18 , 0.159 , 0.121 , 0.088 , 0.067 , 0.054 , ...
%     0.046 , 0.038 , 0.029 , 0.017 , 0.013]; 
% figure()
% inc = {incMat , incMat , incMat , incMat};
% tVec = all90.tVec;
% tVec = linspace(startYear , endYear , size(all90.popVec , 1));
% for i = 1 : length(newHiv_Arr)
%     for a = 4 : age
%         newHiv = sum(sum(sum(newHiv_Arr{i}(1 : end , 1 : gender , a , :)...
%             ,2),3),4);
%         popVec = popVec_Arr{i};
%         hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
%             1 : periods , 1 : gender , a , 1 : risk)); ...
%             toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
%             1 : periods , 1 : gender , a , 1 : risk))];
%         hivSus = sum(popVec(1 : end , hivSusInds) , 2);
%         hivSus_Year = sum(reshape(hivSus , stepsPerYear , size(hivSus , 1) ...
%             / stepsPerYear)) ./ stepsPerYear; % average susceptible population size per year
%         newHiv_Year = sum(reshape(newHiv , stepsPerYear , size(newHiv , 1) ...
%             /stepsPerYear)); % total new HIV infections per year
%         inc{i}(a , :) = newHiv_Year ./ hivSus_Year .* 100;
%     end
% end
% 
% for i = 1 : length(newHiv_Arr)
%     incAS = sum(bsxfun(@times , inc{i} , wVec));
%     plot(tVec(1 : stepsPerYear : end) , incAS)
%     xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
%     hold on
% end
% xlabel('Year'); ylabel('Age-Standardized Incidence per 100'); title('HIV Incidence')
% % Women aged 16-29:  35% on treatment (30% suppressed)
% % Women aged 30+: 60% on treatment (55% suppressed)
% % Men aged 16-29: 25% on treatment (20% suppressed)
% % Men aged 30+: 50% on treatment (40% suppressed).
% legend('70% ART coverage in all ages' , '70% in 15-29 year olds; 90% in 30+' , ...
%     '70% in 15-29 year olds; 95% in 30+', '90% ART coverage in all ages' , ... 
%     'Location' , 'northeastoutside')
% %% Difference in Age-standardized Disease Incidence
% % relative to base case, i.e. 70% coverage for all
% figure()
% for i = 2 : length(newHiv_Arr)
%     incAS = sum(bsxfun(@times , inc{i} , wVec)) ...
%         - sum(bsxfun(@times , inc{1} , wVec));
%     plot(tVec(1 : stepsPerYear : end) , incAS)
%     xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
%     hold on
% end
% xlabel('Year'); ylabel('Difference in Incidence per 100'); 
% title('Difference in HIV Incidence (vs Base Case)')
% % Women aged 16-29:  35% on treatment (30% suppressed)
% % Women aged 30+: 60% on treatment (55% suppressed)
% % Men aged 16-29: 25% on treatment (20% suppressed)
% % Men aged 30+: 50% on treatment (40% suppressed).
% legend('90/90' , '70/90' , '70/95', ...
%     'Location' , 'northeastoutside')
% 
% % percentage form
% figure()
% for i = 2 : length(newHiv_Arr)
%     incAS = (sum(bsxfun(@times , inc{i} , wVec)) ...
%         - sum(bsxfun(@times , inc{1} , wVec))) ...
%         ./ sum(bsxfun(@times , inc{1} , wVec)) * 100;
%     plot(tVec(1 : stepsPerYear : end) , incAS)
%     xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
%     hold on
% end
% xlabel('Year'); ylabel('Difference in Incidence (%)'); 
% title('Difference in Age-Standardized HIV Incidence (vs Base Case)')
% % Women aged 16-29:  35% on treatment (30% suppressed)
% % Women aged 30+: 60% on treatment (55% suppressed)
% % Men aged 16-29: 25% on treatment (20% suppressed)
% % Men aged 30+: 50% on treatment (40% suppressed).
% legend('70% in 15-29 year olds; 90% in 30+' , ...
%     '70% in 15-29 year olds; 95% in 30+', '90% ART coverage in all ages' , ...
%     'Location' , 'northeastoutside')

%% Non Age-standardized Disease Incidence
inc = {incMat , incMat , incMat , incMat};
tVec = all90.tVec;
newHiv_Year = zeros(length(newHiv_Arr) , length(tVec) / stepsPerYear);

for i = 1 : length(newHiv_Arr)
    newHiv = sum(sum(sum(newHiv_Arr{i}(1 : end , 1 : gender , 4 : age , :)...
        ,2),3),4);
    popVec = popVec_Arr{i};
    hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 4 : age , 1 : risk)); ...
        toInd(allcomb(7 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
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
xlabel('Year'); ylabel('Incidence per 100'); title('HIV Incidence')
% Women aged 16-29:  35% on treatment (30% suppressed)
% Women aged 30+: 60% on treatment (55% suppressed)
% Men aged 16-29: 25% on treatment (20% suppressed)
% Men aged 30+: 50% on treatment (40% suppressed).
legend('70% ART coverage in all ages' , '70% in 15-29 year olds; 90% in 30+' , ...
    '70% in 15-29 year olds; 95% in 30+', '90% ART coverage in all ages' , ... 
    'Location' , 'northeastoutside')

figure()
for i = 1 : length(newHiv_Arr)
    plot(tVec(1 : stepsPerYear : end) , newHiv_Year(i , :))
    xlim([tVec(yearNow - stepsPerYear) , tVec(end)])
    hold on
end
xlabel('Year'); ylabel('New HIV Cases'); title('New HIV Cases Per Year')
% Women aged 16-29:  35% on treatment (30% suppressed)
% Women aged 30+: 60% on treatment (55% suppressed)
% Men aged 16-29: 25% on treatment (20% suppressed)
% Men aged 30+: 50% on treatment (40% suppressed).
legend('70% ART coverage in all ages' , '70% in 15-29 year olds; 90% in 30+' , ...
    '70% in 15-29 year olds; 95% in 30+', '90% ART coverage in all ages' , ... 
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
    plot(tVec , sum(all70.popVec(: , artInds) , 2) ./ sum(all70.popVec(: , hivInds) , 2) * 100 ,...
        tVec , sum(all90.popVec(: , artInds) , 2) ./ sum(all90.popVec(: , hivInds) , 2) * 100 ,...
        tVec , sum(art70_90.popVec(: , artInds) , 2) ./ sum(art70_90.popVec(: , hivInds) , 2) * 100 , ...
        tVec , sum(art70_95.popVec(: , artInds) , 2) ./ sum(art70_95.popVec(: , hivInds) , 2) * 100)
    hold on
end
hold off
legend('70/70 (<30)' , '90/90 (<30)' , '70/90 (<30)' , '70/95 (<30)', ...
    '70/70 (30+)' , '90/90 (30+)' , '70/90 (30+)' , '70/95 (30+)', ...
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
plot(tVec , sum(all70.popVec(: , hivInds) , 2) ./ sum(all70.popVec(: , allInds) , 2) * 100 ,...
    tVec , sum(all90.popVec(: , hivInds) , 2) ./ sum(all90.popVec(: , allInds) , 2) * 100 ,...
    tVec , sum(art70_90.popVec(: , hivInds) , 2) ./ sum(art70_90.popVec(: , allInds) , 2) * 100 , ...
    tVec , sum(art70_95.popVec(: , hivInds) , 2) ./ sum(art70_95.popVec(: , allInds) , 2) * 100)
legend('70/70' , '90/90' , '70/90' , '70/95', ...
    'Location' , 'northeastoutside')
title('HIV Prevalence')
xlabel('Year'); ylabel('Prevalence (%)')

%% Circumcised
circInds = [toInd(allcomb(7, 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 , 4 : 10 , 1 : risk))];
allMInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 , 4 : 10 , 1 : risk))];
figure()
plot(tVec , sum(all70.popVec(: , circInds) , 2) ./ sum(all70.popVec(: , allMInds) , 2) * 100 ,...
    tVec , sum(all90.popVec(: , circInds) , 2) ./ sum(all90.popVec(: , allMInds) , 2) * 100 ,...
    tVec , sum(art70_90.popVec(: , circInds) , 2) ./ sum(art70_90.popVec(: , allMInds) , 2) * 100 , ...
    tVec , sum(art70_95.popVec(: , circInds) , 2) ./ sum(art70_95.popVec(: , allMInds) , 2) * 100)
legend('70/70' , '90/90' , '70/90' , '70/95', ...
    'Location' , 'northeastoutside')
title('Circumcision')
xlabel('Year'); ylabel('Coverage (%)')

%% Age distribution

allInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : 2 , 4 : age , 1 : risk))];
ages = {4 : 6 , 7 : age}; % 15-29, 30+
all70Prop = zeros(size(ages , 2) , length(all90.tVec));
all90Prop = all70Prop;
prop70_90 = all70Prop;
prop70_95 = all70Prop;

for i = 1 : length(ages)
    a = ages{i};
    ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , 1 : 2 , a , 1 : risk));
    all70Prop(i , :) = sum(all70.popVec(: , ageInds) , 2) ./ sum(all70.popVec(: , allInds) , 2) * 100;
    all90Prop(i , :) =  sum(all90.popVec(: , ageInds) , 2) ./ sum(all90.popVec(: , allInds) , 2) * 100;
    prop70_90(i , :) = sum(art70_90.popVec(: , ageInds) , 2) ./ sum(art70_90.popVec(: , allInds) , 2) * 100;
    prop70_95(i , :) = sum(art70_95.popVec(: , ageInds) , 2) ./ sum(art70_95.popVec(: , allInds) , 2) * 100;
end

ageDist70 = mean(all70Prop , 2);
stDev70 = std(all70Prop , 0  , 2);
ageDistAll90 = mean(all90Prop , 2);
stDev90 = std(all90Prop , 0 , 2);
ageDist70_90 = mean(prop70_90 , 2);
stDev70_90 = std(prop70_90 , 0 , 2);
ageDist70_95 = mean(prop70_95 , 2);
stDev70_95 = std(prop70_95 , 0 , 2);

ageDist = table({'15-29' ; '30+'} , ageDist70 , stDev70 , ageDistAll90 ,...
    stDev90 , ageDist70_90 , stDev70_90 , ageDist70_95 , stDev70_95);
ageDist.Properties.VariableNames{1} = 'AgeGroup';
file = 'ageDistArt.xlsx';
writetable(ageDist , file)
disp(ageDist)

%%
figure()
plot(tVec , all90Prop)    
xlabel('Year'); ylabel('Age Distribution (%)')
title('Age Distribution Under 90% ART for All Ages Scenario')
legend('15-29' , '30+')