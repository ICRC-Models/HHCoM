function vaxCEA_vclir061019 %(pathModifier)

close all; clear all;

waning = 0;    % turn waning on or off
pathModifier = '060719_VCLIR_SCE8';

%% Load parameters
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'],'stepsPerYear','circ','condUse','disease','viral',...
    'hpvTypes','hpvStates','periods','gender','age','risk','dim','k','toInd','sumall','modelYr1')

sumall = @(x) sum(x(:));

% Load results
nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
curr = load([pwd , '\HHCoM_Results\toNow_060719_VCLIR_SCE8']); % Population up to current year

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
    vaxResult{n}.newVaxHpv= [curr.newVaxHpv(1 : end , : , : , : , :) ; vaxResult{n}.newVaxHpv];
    %vaxResult{n}.ccDeath = [curr.ccDeath(1 : end - 1 , : , : , :) ; vaxResult{n}.ccDeath];
    vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :) ; vaxResult{n}.newCC];
    vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , :) ; vaxResult{n}.newHiv];
    vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :) ; vaxResult{n}.artTreatTracker];
    vaxResult{n}.tVec = [curr.tVec(1 : end) , vaxResult{n}.tVec];
%     vaxResult{n}.ccTreated = [curr.ccTreated(1 : end - 1) , vaxResult{n}.ccTreated];
end

% Find no vaccine scenario
% noVaxInd = -1;
% for n = 1 : nSims
%     if vaxResult{n}.vaxEff == 0
%         noVaxInd = n;
%     end
% end
noVaxInd = nSims;

noV = vaxResult{noVaxInd};
tVec = noV.tVec;
tVecYr = tVec(1 : stepsPerYear : end);

%% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)

%% Population Size
figure()
fac = 10 ^ 5;
inds = {[1:age],1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
xlC = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q'];
noV_F = [];
incV_F = [];
noV_M = [];
incV_M = [];

for aV = 1 : length(inds)
    a = inds{aV};
    totalPop = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , 1 : risk));
    females = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 2 , a , 1 : risk));
    males = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 , a , 1 : risk));

    noV_F = [noV_F , (sum(noV.popVec(1:stepsPerYear:end , females) , 2) ./ sum(noV.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
    incV_F = [incV_F , (sum(vaxResult{1}.popVec(1:stepsPerYear:end , females) , 2) ./ ...
        sum(vaxResult{1}.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];

    noV_M = [noV_M , (sum(noV.popVec(1:stepsPerYear:end , males) , 2) ./ sum(noV.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
    incV_M = [incV_M , (sum(vaxResult{1}.popVec(1:stepsPerYear:end , males) , 2) ./ ...
        sum(vaxResult{1}.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
end

plot(tVec(1:stepsPerYear:end) , noV_F(:,1))
hold all;
plot(tVec(1:stepsPerYear:end) , incV_F(:,1) , '--')
hold all;
plot(tVec(1:stepsPerYear:end) , noV_M(:,1))
hold all;
plot(tVec(1:stepsPerYear:end) , incV_M(:,1) , '--')
title('Total number of men/women (per 100,000 total population)')
xlabel('Year'); ylabel('Individuals')
xlim([1910 2120]);
legend('Women: noV' , 'Women: 50%' , 'Men: noV' , 'Men: 50%');

% Save results
fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\' , pathModifier , 'ByAge_PS_all' , '.xlsx'];
%%xlRange = [xlC(aV) , '1' , ':' , xlC(aV) , num2str(length(noV_F))]; %xlswrite(filename,A,sheet,xlRange)
%     if aV ~= 1 %exist(fname , 'file') == 2
%         M1 = xlsread(fname , 'PSnoVF');
%         M1 = catpad(2 , M1 , noV_F);
%         xlswrite(fname , M1 , 'PSnoVF')
%         M2 = xlsread(fname , 'PSincVF');
%         M2 = catpad(2 , M2 , incV_F);
%         xlswrite(fname , M2 , 'PSincVF')
%         M3 = xlsread(fname , 'PSnoVM');
%         M3 = catpad(2 , M3 , noV_M);
%         xlswrite(fname , M3 , 'PSnoVM')
%         M4 = xlsread(fname , 'PSincVM');
%         M4 = catpad(2 , M4 , incV_M);
%         xlswrite(fname , M4 , 'PSincVM')
%     elseif aV == 1
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV_F] , 'noVF')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , incV_F] , 'incVF')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV_M] , 'noVM')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , incV_M] , 'incVM')
%     end

%% HPV prevalence
figure()
fac = 10 ^ 5;
inds = {[1:age],1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
noV_FHPV = [];
incV_FHPV = [];
noV_MHPV = [];
incV_MHPV = [];

for aV = 1 : length(inds)
    a = inds{aV};
    totalPop = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , 1 : risk));
    femalesHPV_hpv = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 , ...
        [1,6] , 2 , a , 1 : risk));
    femalesHPV_cin1 = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 2 , ...
        [1,6] , 2 , a , 1 : risk));
    femalesHPV_cin23 = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , [3:4] , ...
        [1,6] , 2 , a , 1 : risk));
    femalesHPV_cc = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , [5:7] , ...
        [1,6] , 2 , a , 1 : risk));
    malesHPV_hpv = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 , ...
        [1,6] , 1 , a , 1 : risk));
    malesHPV_cin1 = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 2 , ...
        [1,6] , 1 , a , 1 : risk));
    malesHPV_cin23 = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , [3:4] , ...
        [1,6] , 1 , a , 1 : risk));
    malesHPV_cc = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , [5:7] , ...
        [1,6] , 1 , a , 1 : risk));

    noV_FHPV = [noV_FHPV , ((sum(noV.popVec(1:stepsPerYear:end , femalesHPV_hpv) , 2).*0.14 + sum(noV.popVec(1:stepsPerYear:end , femalesHPV_cin1) , 2).*0.16 + ...
        sum(noV.popVec(1:stepsPerYear:end , femalesHPV_cin23) , 2).*0.41 + sum(noV.popVec(1:stepsPerYear:end , femalesHPV_cc) , 2).*0.50) ...
        ./ sum(noV.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
    %noV_FHPV = round(noV_FHPV,20);
    incV_FHPV = [incV_FHPV , ((sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesHPV_hpv) , 2).*0.14 + sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesHPV_cin1) , 2).*0.16 + ...
        sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesHPV_cin23) , 2).*0.41 + sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesHPV_cc) , 2).*0.50) ...
        ./ sum(vaxResult{1}.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
    %incV_FHPV = round(incV_FHPV,20);

    noV_MHPV = [noV_MHPV , ((sum(noV.popVec(1:stepsPerYear:end , malesHPV_hpv) , 2).*0.14 + sum(noV.popVec(1:stepsPerYear:end , malesHPV_cin1) , 2).*0.16 + ...
        sum(noV.popVec(1:stepsPerYear:end , malesHPV_cin23) , 2).*0.41 + sum(noV.popVec(1:stepsPerYear:end , malesHPV_cc) , 2).*0.50) ...
        ./ sum(noV.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
    %noV_MHPV = round(noV_MHPV,20);
    incV_MHPV = [incV_MHPV , ((sum(vaxResult{1}.popVec(1:stepsPerYear:end , malesHPV_hpv) , 2).*0.14 + sum(vaxResult{1}.popVec(1:stepsPerYear:end , malesHPV_cin1) , 2).*0.16 + ...
        sum(vaxResult{1}.popVec(1:stepsPerYear:end , malesHPV_cin23) , 2).*0.41 + sum(vaxResult{1}.popVec(1:stepsPerYear:end , malesHPV_cc) , 2).*0.50) ...
        ./ sum(vaxResult{1}.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
    %incV_MHPV = round(incV_MHPV,20);
end

plot(tVec(1:stepsPerYear:end)' , noV_FHPV(:,1))
hold all;
plot(tVec(1:stepsPerYear:end)' , incV_FHPV(:,1) , '--')
hold all;
plot(tVec(1:stepsPerYear:end)' , noV_MHPV(:,1))
hold all;
plot(tVec(1:stepsPerYear:end)' , incV_MHPV(:,1) , '--')
title('Total number of HPV-infected men/women (per 100,000 total population)')
xlabel('Year'); ylabel('Individuals')
xlim([1910 2120]);
legend('Women: noV' , 'Women: 50%' , 'Men: noV' , 'Men: 50%');
    
% Save results
fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\' , pathModifier , 'ByAge_HPV_all' , '.xlsx'];
%     if aV > 1 %exist(fname , 'file') == 2
%         M1 = xlsread(fname , 'HPVnoVFHPV');
%         M1 = catpad(2 , M1 , noV_FHPV);
%         xlswrite(fname , M1 , 'HPVnoVFHPV')
%         M2 = xlsread(fname , 'HPVincVFHPV');
%         M2 = catpad(2 , M2 , incV_FHPV);
%         xlswrite(fname , M2 , 'HPVincVFHPV')
%         M3 = xlsread(fname , 'HPVnoVMHPV');
%         M3 = catpad(2 , M3 , noV_MHPV);
%         xlswrite(fname , M3 , 'HPVnoVMHPV')
%         M4 = xlsread(fname , 'HPVincVMHPV');
%         M4 = catpad(2 , M4 , incV_MHPV);
%         xlswrite(fname , M4 , 'HPVincVMHPV')
%     elseif aV == 1
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV_FHPV] , 'noVFHPV')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , incV_FHPV] , 'incVFHPV')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV_MHPV] , 'noVMHPV')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , incV_MHPV] , 'incVMHPV')
%     end

%% Vaccinated women
figure()
fac = 10 ^ 5;
inds = {[1:age],1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
noV_Fvax = [];
incV_Fvax = [];
noV_FvaxHPV = [];
incV_FvaxHPV = [];
noV_FnoVaxHPV = [];
incV_FnoVaxHPV = [];
noV_FnoVax = [];
incV_FnoVax = [];

for aV = 1 : length(inds)
    a = inds{aV};
    totalPop = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , 1 : risk));
    femalesVax = [toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , ...
        [1,6] , 2 , a , 1 : risk)); toInd(allcomb(1: disease , 1 : viral , ...
        1: hpvTypes , 1 : hpvStates , [2,4] , 2 , a , 1 : risk))];
    femalesVaxHPV = toInd(allcomb(1: disease , 1 : viral , 2 : hpvTypes , ...
        1 : 7 , [2,4] , 2 , a , 1 : risk));
    femalesNoVaxHPV_hpv = toInd(allcomb(1: disease , 1 : viral , 2 : hpvTypes , ...
        1 , [1,6] , 2 , a , 1 : risk));
    femalesNoVaxHPV_cin1 = toInd(allcomb(1: disease , 1 : viral , 2 : hpvTypes , ...
        2 , [1,6] , 2 , a , 1 : risk));
    femalesNoVaxHPV_cin23 = toInd(allcomb(1: disease , 1 : viral , 2 : hpvTypes , ...
        [3:4] , [1,6] , 2 , a , 1 : risk));
    femalesNoVaxHPV_cc = toInd(allcomb(1: disease , 1 : viral , 2 : hpvTypes , ...
        [5:7] , [1,6] , 2 , a , 1 : risk));
    femalesNoVax = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , [1:8,10] , ...
        [1,6] , 2 , a , 1 : risk));
    
    noV_Fvax = [noV_Fvax , (sum(noV.popVec(1:stepsPerYear:end , femalesVax) , 2) ./ sum(noV.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
    incV_Fvax = [incV_Fvax , (sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesVax) , 2) ./ ...
        sum(vaxResult{1}.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];

    noV_FvaxHPV = [noV_FvaxHPV , (sum(noV.popVec(1:stepsPerYear:end , femalesVaxHPV) , 2) ./ sum(noV.popVec(1:stepsPerYear:end , totalPop) , 2)) .* (fac * 0.0)];
    incV_FvaxHPV = [incV_FvaxHPV , (sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesVaxHPV) , 2) ./ ...
        sum(vaxResult{1}.popVec(1:stepsPerYear:end , totalPop) , 2)) .* (fac * 0.0)];

    noV_FnoVaxHPV = [noV_FnoVaxHPV , ((sum(noV.popVec(1:stepsPerYear:end , femalesNoVaxHPV_hpv) , 2).*0.14 + sum(noV.popVec(1:stepsPerYear:end , femalesNoVaxHPV_cin1) , 2).*0.16 + ...
        sum(noV.popVec(1:stepsPerYear:end , femalesNoVaxHPV_cin23) , 2).*0.41 + sum(noV.popVec(1:stepsPerYear:end , femalesNoVaxHPV_cc) , 2).*0.50) ...
        ./ sum(noV.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
    %noV_FnoVaxHPV = round(noV_FnoVaxHPV,20);
    incV_FnoVaxHPV = [incV_FnoVaxHPV , ((sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesNoVaxHPV_hpv) , 2).*0.14 + sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesNoVaxHPV_cin1) , 2).*0.16 + ...
        sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesNoVaxHPV_cin23) , 2).*0.41 + sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesNoVaxHPV_cc) , 2).*0.50) ...
        ./ sum(vaxResult{1}.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
    %incV_FnoVaxHPV = round(incV_FnoVaxHPV,20);
    
    noV_FnoVax = [noV_FnoVax , (sum(noV.popVec(1:stepsPerYear:end , femalesNoVax) , 2) ./ sum(noV.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
    incV_FnoVax = [incV_FnoVax , (sum(vaxResult{1}.popVec(1:stepsPerYear:end , femalesNoVax) , 2) ./ ...
        sum(vaxResult{1}.popVec(1:stepsPerYear:end , totalPop) , 2)) .* fac];
end

plot(tVec(1:stepsPerYear:end) , noV_Fvax(:,1))
hold all;
plot(tVec(1:stepsPerYear:end) , incV_Fvax(:,1) , '--')
hold all;
plot(tVec(1:stepsPerYear:end) , noV_FvaxHPV(:,1))
hold all;
plot(tVec(1:stepsPerYear:end) , incV_FvaxHPV(:,1) , '--')
hold all;
plot(tVec(1:stepsPerYear:end) , noV_FnoVaxHPV(:,1))
hold all;
plot(tVec(1:stepsPerYear:end) , incV_FnoVaxHPV(:,1) , '--')
hold all;
plot(tVec(1:stepsPerYear:end) , noV_FnoVax(:,1))
hold all;
plot(tVec(1:stepsPerYear:end) , incV_FnoVax(:,1) , '--')
title('Total number of vaccinated/vax+HPV/unvax+HPV women (per 100,000 total population)')
xlabel('Year'); ylabel('Individuals')
xlim([1910 2120]);
legend('Vax: noV' , 'Vax: 50%' , 'Vax+HPV: noV' , 'Vax+HPV: 50%' , 'NoVax+HPV: noV' , ...
    'NoVax+HPV: 50%' , 'NoVax: noV' , 'NoVax: 50%');

% Save results
fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\' , pathModifier , 'ByAge_Vax_all' , '.xlsx'];
%     if exist(fname , 'file') == 2
%         M1 = xlsread(fname , 'VnoVFvax');
%         M1 = catpad(2 , M1 , noV_Fvax);
%         xlswrite(fname , M1 , 'VnoVFvax')
%         M2 = xlsread(fname , 'VincVFvax');
%         M2 = catpad(2 , M2 , incV_Fvax);
%         xlswrite(fname , M2 , 'VincVFvax')
%         M3 = xlsread(fname , 'VnoVFvaxHPV');
%         M3 = catpad(2 , M3 , noV_FvaxHPV);
%         xlswrite(fname , M3 , 'VnoVFvaxHPV')
%         M4 = xlsread(fname , 'VincVFvaxHPV');
%         M4 = catpad(2 , M4 , incV_FvaxHPV);
%         xlswrite(fname , M4 , 'VincVFvaxHPV')
%     else
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV_Fvax] , 'VnoVFvax')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , incV_Fvax] , 'VincVFvax')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV_FvaxHPV] , 'VnoVFvaxHPV')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , incV_FvaxHPV] , 'VincVFvaxHPV')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV_FnoVaxHPV] , 'VnoVFnoVaxHPV')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , incV_FnoVaxHPV] , 'VincVFnoVaxHPV')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , noV_FnoVax] , 'VnoVFnoVax')
xlswrite(fname , [tVec(1 : stepsPerYear : end)' , incV_FnoVax] , 'VincVFnoVax')
%     end

