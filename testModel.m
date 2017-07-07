% Test script. Allows for module by module testing.
close all; clear all; clc
profile clear;
% [~ , startYear , endYear , stepsPerYr , IntSteps] = Menu();
% disp('Done with input')
%%
% calcInds(stepsPerYear); % uncomment if recalculating indices
stepsPerYear = 4;
% Get parameter values and load model
prompt = {'Start year' , 'End year' , ...
    'Model HIV? (1 = yes , 0 = no)' , ...
    'Model HPV? (1 = yes , 0 = no)' , 'Load new data? (1 = yes , 0 = no)'};
tit = 'Model Inputs';
inputs = inputdlg(prompt , tit);
loadNew = str2double(inputs{5});
if isempty(inputs)
    return
end
disp('Initializing. Standby...')
disp(' ');
% Choose whether to model HPV
hpvOn = str2double(inputs{4});
if hpvOn
    disp('HPV module activated.')
end
% choose whether to model hysterectomy
hyst = 'off';
% choose whether to model HIV
hivOn = str2double(inputs{3});
if hivOn
    disp('HIV module activated')
end
disp(' ')
startYear = str2double(inputs{1});
endYear = str2double(inputs{2});
if loadNew
    loadUp();
else
    disp('Skipped parameter load up.')
end
% Load parameters and constants for main
load('general')
years = endYear - startYear;
%% Initial population 
load('popData')
mInit = popInit(: , 1);
MsumInit = sum(mInit);

fInit = popInit(: , 2);
FsumInit = sum(fInit);

MpopStruc = riskDistM;
FpopStruc = riskDistF;

mPop = zeros(age , risk);
fPop = mPop;

for i = 1 : age
    mPop(i , :) = MpopStruc(i, :).* mInit(i);
    fPop(i , :) = FpopStruc(i, :).* fInit(i);
end

dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age ,risk];
initPop = zeros(dim);
initPop(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop;
initPop(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop;

if hivOn
    initPop(6 , 4 , 1 , 1 , 1 , 1 , 6 , 2) = 10000; % initial HIV infected male
    initPop(6 , 4 , 1 , 1 , 1 , 2 , 5 , 2) = 10000; % initial HIV infected female
end

if hpvOn
    % initPop(1 , 1 , 2 : 4 , 1 , 1 , : , 4 : 9 , :) = 2; % initial HPV hr and lr infecteds (test)
    initPop(1 , 1 , 2 : 4 , 1 , 1 , : , 4 : 9 , :) = sum(initPop(:)) * 0.10 / 108;
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
        initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - sum(initPop(:)) * 0.10 / 36;
end
%% Simulation
profile on
disp(' ')
% Initialize vectors
timeStep = 1 / stepsPerYear;
s = 1 : timeStep : years + 1; % stepSize and steps calculated in loadUp.m
% artDistMat = zeros(size(prod(dim) , 20)); % initialize artDistMat to track artDist over past 20 time steps
% performance tracking
runtimes = zeros(size(s , 2) - 2 , 1);
import java.util.LinkedList
artDistList = LinkedList();
% popVec = zeros(disease , viral, hpvTypes , hpvStates , periods, ...
%    gender , age , risk , length(s));
popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8); 
popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
popVec(1 , :) = popIn;
artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0
disp(['Simulating period from ' num2str(startYear) ' to ' num2str(endYear)])
disp(' ')
disp('Simulation running (test mode)...')
disp(' ')
% Equations for each module in simulation
% Evaluate diff eqs over one time interval. Use output from preceding
% differential equations as input to successive differential equations.
prompt = {'Test mixInfect (1 = yes, 0 = no)' , 'Test vlAdv' , 'Test hiv' , ...
    'Test hpv' , 'Test cinAdv' , 'Test hpvTreat' , 'Test agePop' , 'Test bornDie'};
tit = 'Test options (1 = yes , 0 = no)';
inputs = inputdlg(prompt , tit);
if isempty(inputs)
    return
end
mixInfectTest = str2double(inputs{1});
vlAdvTest = str2double(inputs{2});
hivTest = str2double(inputs{3});
hpvTest = str2double(inputs{4});
cinAdvTest = str2double(inputs{5});
hpvTreatTest = str2double(inputs{6});
agePopTest = str2double(inputs{7});
bornDieTest = str2double(inputs{8});
pop = popIn';
newHiv = zeros(length(s) - 1 , gender);
newHpv = newHiv;
hivDeaths = zeros(length(s) - 1 , age);
for i = 2 : length(s) - 1
    tic
    year = floor(startYear + s(i) - 1);
    currStep = round(s(i) * stepsPerYear);
    disp(['current step = ' num2str(startYear + s(i) - 1) ' ('...
        num2str(length(s) - i) ' time steps remaining until year ' ...
        num2str(endYear) ')'])
    tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
    if mixInfectTest
        [~ , pop , newHiv(i , :) , newHpv(i , :)] = ode4xtra(@(t , pop) mixInfect(t , pop , 'vl' , ...
        currStep) , tspan , popIn);
    end
    if hivOn
        if vlAdvTest
            [~ , pop] = ode4x(@vlAdv , tspan , pop(end , :));
        end
        if hivTest
            [~ , pop , hivDeaths(i , :) , artTreat] = ode4xtra(@(t , pop) hiv(t , pop , artDist , year) ,...
                tspan , pop(end , :));
        end
        popCopy = pop;
        
        if size(artDistList) >= 20
            artDistList.remove(); % remove earlier artDist matrix more than "20 time steps old"
        else
            artDistList.add(artTreat);
        end
        artDist = calcDist(artDistList);
    end
    if hpvOn
        if hpvTest
            [~ , pop] = ode4x(@(t , pop) hpv(t , pop , hyst) , ...
                tspan , pop(end , :));
        end
        if cinAdvTest
            [~ , pop] = ode4x(@cinAdv , tspan , pop(end , :));
        end
        if hpvTreatTest
            [~ , pop] = ode4x(@(t , pop) hpvTreat(t , pop , hyst , year)...
                , tspan , pop(end , :)); 
        end
    end
    if bornDieTest
%         origPop = pop(end , 1) + pop(end , 7201);
        [~ , pop] = ode4x(@(t , pop) bornDie(t , pop , year , currStep) , tspan , pop(end , :));
%         deltaPop = (pop(end , 1) + pop(end , 7201)) - origPop;
%         disp(num2str(deltaPop))
    end
    if agePopTest
%         origPop = pop(end , :);
        [~ , pop] = ode4x(@agePop , tspan , pop(end , :));
%         deltaPop = sumall(pop(end , :)) - sumall(origPop(:));
%         disp([num2str(deltaPop), ' after age module'])
    end
    % add results to population vector
    popVec(i , :) = pop(end , :);
    popIn = popVec(i , :);
    runtimes(i - 1) = toc;
end
disp(['Reached year ' num2str(endYear)])
popVec = sparse(popVec); % compress population vectors
disp(' ')
disp('Simulation complete.')
profile viewer
tVec = linspace(startYear , endYear , size(popVec , 1));
%% 
figure()
plot(1 : size(runtimes , 1) , runtimes)
xlabel('Step'); ylabel('Time(s)')
title('Runtimes')
%% 
avgRuntime = mean(runtimes); % seconds
stdRuntime = std(runtimes); % seconds
disp(['Total runtime: ' , num2str(sum(runtimes) / 3600) , ' hrs' , ' (' , num2str(sum(runtimes) / 60) , ' mins)']); 
disp(['Average runtime per step: ' , num2str(avgRuntime / 60) , ' mins (' , num2str(avgRuntime) , ' secs)']);
disp(['Standard deviation: ' , num2str(stdRuntime / 60) , ' mins (' , num2str(stdRuntime) , ' secs)']);
figure()
h = histogram(runtimes);
title('Runtimes')
ylabel('Frequency')
xlabel('Times (s)')
% %% load target values
% load('targetVals')
%% Plot results
load('groupedInds');
% Total HIV positive
hiv = simPlot(popVec , tVec , hivInds , 'Year' , 'Population' , 'HIV Positive Individuals');
%% CD4 < 500
cd4_500 = simPlot(popVec , tVec , cd4_500Inds , 'Year' , 'Population' , 'CD4 < 500');
%% HIV status by age
% for a = 1 : age
%     hivAgeInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         1 : gender , a , 1 : risk));
%     simPlot(popVec , tVec , hivAgeInds , 'Year' , 'Population' , ['Age group ' , num2str(a) , ' with HIV']);
%     artAgeInds = toInd(allcomb(10 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
%         1 : gender , a , 1 : risk));
%     simPlot(popVec , tVec , artAgeInds , 'Year' , 'Population' , ['Age group ' , num2str(a) , ' on ART']);
% end
%% All ages
% relative
% for a = 1 : age
%     ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%         1 : periods , 1 : gender , a , 1 : risk));
%     simPlot(popVec , tVec , ageInds , 'Year' , 'Population Size' , ['Age group ' , num2str(a)]);
% end
% % absolute
% for a = 1 : age
%     ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
%         1 : periods , 1 : gender , a , 1 : risk));
%     simPlotAbs(popVec , tVec , ageInds , 'Year' , 'Population Size' , ['Age group ' , num2str(a)]);
% end
%% Total on ART
art = simPlot(popVec , tVec , artInds, 'Year' , 'Population' , 'Individuals on ART');
%%
figure()
plot(tVec , art ./ sum(popVec , 2)' , tVec , hiv ./ sum(popVec , 2)')
xlabel('Year')
ylabel('Proportion of Population')
legend('On ART', 'HIV+')
%% HPV
hpvTot = simPlot(popVec , tVec , hpvInds , 'Year' , 'Individuals' , 'HPV Prevalence');

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
        1 : hpvStates , 1 :  periods , 1 , 4 : age , 1 : risk));
    m15_64all(i) = sumall(popVec(i , m15_64Inds));
    
    f15_64Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
        1 : hpvStates , 1 : periods , 2 , 4 : age , 1 : risk));
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
    3046456] * 0.9;

% females 15 - 64
actual(4 , :) = [1655333
    2711917
    3133521
    3433274] * 0.9;
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

%% Disease incidence
figure()
%plot(tVec , newHiv)
area(tVec , bsxfun(@rdivide , newHiv , sum(popVec , 2)) * 100)
xlabel('Year'); ylabel('New cases (%)'); title('Relative HIV Incidence')
legend('Male' , 'Female')
figure()
area(tVec , bsxfun(@rdivide , newHpv , sum(popVec , 2)) * 100)
xlabel('Year'); ylabel('New cases (%)'); title('Relative HPV Incidence')
legend('Male' , 'Female')
