function vaxCEA_CISNETvaxCompare(pathModifier)
% modifying function vaxCEA() to produce outputs for CISNET vaxination
% scenario comparison model runs

waning = 0;    % turn waning on or off

%% Load parameters
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'],'stepsPerYear','circ','condUse','disease','viral',...
    'hpvTypes','hpvStates','periods','gender','age','risk','dim','k','toInd','sumall','modelYr1')
% Load results
nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
curr = load([pwd , '\HHCoM_Results\toNow.mat']); % Population up to current year

% Helper functions
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

% Time
c = fix(clock); % get time
currYear = c(1); % get the current year from time

vaxResult = cell(nSims , 1);
resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxSimResult'];
if waning
    resultFileName = [pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , 'vaxWaneSimResult'];
end
parfor n = 1 : nSims
    % load results from vaccine run into cell array
    vaxResult{n} = load([resultFileName , num2str(n), '.mat']);
    % concatenate vectors/matrices of population up to current year to population
    % matrices for years past current year
    vaxResult{n}.popVec = [curr.popVec(1 : end  , :) ; vaxResult{n}.popVec];
    vaxResult{n}.newHpv= [curr.newHpv(1 : end , : , : , : , :) ; vaxResult{n}.newHpv];
    vaxResult{n}.newImmHpv= [curr.newImmHpv(1 : end , : , : , : , :) ; vaxResult{n}.newImmHpv];
%     vaxResult{n}.ccDeath = [curr.ccDeath(1 : end - 1 , : , : , :) ; vaxResult{n}.ccDeath];
    vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :) ; vaxResult{n}.newCC];
    vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , :) ; vaxResult{n}.newHiv];
    vaxResult{n}.tVec = [curr.tVec(1 : end) , vaxResult{n}.tVec];
%     vaxResult{n}.ccTreated = [curr.ccTreated(1 : end - 1) , vaxResult{n}.ccTreated];
end

% Find no vaccine scenario
noVaxInd = -1;
for n = 1 : nSims
    if vaxResult{n}.vaxEff == 0
        noVaxInd = n;
    end
end
noV = vaxResult{noVaxInd};
tVec = noV.tVec;
tVecYr = tVec(1 : stepsPerYear : end);

%% Population Size
figure()
for n = 1 : length(vaxResult)
    plot(tVec , sum(vaxResult{n}.popVec , 2) , 'DisplayName' , ...
        ['Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
        'Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
    legend('-DynamicLegend')
    hold on
end
title('Population Size')
xlabel('Year'); ylabel('Individuals')
hold off

%% Population Size by Gender
figure()
mInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 1 , 1 : age , 1 : risk));
fInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 2 , 1 : age , 1 : risk));
for n = 1 : length(vaxResult)
    plot(tVec , sum(vaxResult{n}.popVec(: , mInds) , 2) , 'DisplayName' , ...
        ['Male, Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
        'Male, Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
    legend('-DynamicLegend')
    hold on
    plot(tVec , sum(vaxResult{n}.popVec(: , fInds) , 2) , 'DisplayName' , ...
        ['Female, Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100)) '% ,', ...
        'Female, Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '%'])
    legend('-DynamicLegend')
    hold on
end
%legend('Males' , 'Females')
title('Population Size')
xlabel('Year'); ylabel('Individuals')
hold off

%% Population Size by Gender and Age
for n = 1 : length(vaxResult)
    for g = 1 : gender
        figure()
        for a = 1 : age
            ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
                1 : periods , g , a , 1 : risk));
            
            plot(tVec , sum(vaxResult{n}.popVec(: , ageInds) , 2))
            hold on
        end
        xlabel('Year'); ylabel('Individuals'); title(['Population Size, ',...
            num2str(g), ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100))...
            '% ,', ' Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) '% ,']);
        legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
            '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
            '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside');
   end
end
hold off

%% Number HPV Infected by Gender and Age
for n = 1 : length(vaxResult)
    for g = 1 : gender
        figure()
        for a = 1 : age
            hpvInfInds = [toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 , ...
                1 : periods , g , a , 1 : risk))];
            
            plot(tVec , sum(vaxResult{n}.popVec(: , hpvInfInds) , 2))
            hold on
        end
        xlabel('Year'); ylabel('Individuals'); title(['HPV Infected, ',...
            num2str(g), ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100))...
            '% ,', ' Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) '% ,']);
        legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
            '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
            '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside');
   end
end
hold off

%% Number Vaccinated Women by Age

for n = 1 : length(vaxResult)
    figure()
    for a = 1 : age
        vaxInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 , ...
            1 : periods , 2 , a , 1 : risk))];

        plot(tVec , sum(vaxResult{n}.popVec(: , vaxInds) , 2))
        hold on
    end
    xlabel('Year'); ylabel('Individuals'); title(['Vaccinated Women, ', ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100))...
        '% ,', ' Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) '% ,']);
    legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
        '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
        '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside');
end
hold off

%% Number HPV Infected Vaccinated Women by Age

for n = 1 : length(vaxResult)
    figure()
    for a = 1 : age
        vaxHpvInfInds = [toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 9 , ...
            1 : periods , 2 , a , 1 : risk))];

        plot(tVec , sum(vaxResult{n}.popVec(: , vaxHpvInfInds) , 2))
        hold on
    end
    xlabel('Year'); ylabel('Individuals'); title(['HPV Infected Vaccinated Women, ', ': Efficacy: ' , num2str(round(vaxResult{n}.vaxEff * 100))...
        '% ,', ' Coverage: ' , num2str(round(vaxResult{n}.vaxRate * 100)) '% ,']);
    legend('0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
        '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' ,...
        '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79' , 'Location' , 'NorthEastOutside');
end
hold off
