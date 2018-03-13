function simVaxResultOut_AgeStand()
%%
load('general')
o90 = load('H:\HHCoM_Results\Vax_0.9_wane_0.mat'); % 90% coverage
o70 = load('H:\HHCoM_Results\Vax_0.7_wane_0.mat'); % 70% coverage
o50 = load('H:\HHCoM_Results\Vax_0.5_wane_0.mat'); % 50% coverage
oNo = load('H:\HHCoM_Results\Vax_0_wane_0.mat'); % No vaccine
tVec = o90.tVec;
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 
midMat = zeros(stepsPerYear , size(o90.popVec , 1) / stepsPerYear);
midMat(1 , :) = 1;
midMat(end , :) = 1;
midAnn = @(x) sum(midMat .* reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) / 2;
%% Plot Settings

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
%%
ageGroups = age - 4 + 1;
wVec = zeros(age , 1);
wVec(5 : age) = [0.188 , 0.18 , 0.159 , 0.121 , 0.088 , 0.067 , 0.054 , ...
    0.046 , 0.038 , 0.029 , 0.017 , 0.013]; 
% if ageStand == 0
%     wVec = ones(age , 1) ./ age;
% end
%% HPV Incidence
inds = {':' , 2 : 6 , 1 , 10};
files = {'General_Hpv_VaxCover' , 'HivAll_Hpv_VaxCover' , 'HivNeg_Hpv_VaxCover' , 'ART_HPV_VaxCover'};
plotTits = {'General HPV' , 'HIV-Positive + HPV' , 'HIV-Negative + HPV' , ...
    'HIV-Positive on ART'};
fac = 10 ^ 5;
vNo_HpvAge = zeros(ageGroups , length(tVec) / stepsPerYear);
v90_HpvAge = vNo_HpvAge;
v70_HpvAge = vNo_HpvAge;
v50_HpvAge = vNo_HpvAge;
fac = 10 ^ 5;
for i = 1 : length(inds)
    for a = 1 : age
        % general
        allF = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        % All HIV-positive women (not on ART)
        allHivF = toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        % All HIV-negative women
        hivNeg = toInd(allcomb(1 , 1 : viral , 1 , 1 : hpvStates , 1 : periods , ...
            2 , a , 1 : risk));
        % Women on ART
        artF = toInd(allcomb(10 , 6 , 1 , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        
        genArray = {allF , allHivF , hivNeg , artF};
        
        vNo_HpvAge(a , :) = ...
            annlz(sum(sum(sum(sum(oNo.newHpv(: , 2 , inds{i} , a , :),2),3),4),5)) ./ ...
            (annlz(sum(oNo.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac;
       
        v90_HpvAge(a , :) = ...
            annlz(sum(sum(sum(sum(o90.newHpv(: , 2 , inds{i} , a , :),2),3),4),5)) ./ ...
           	(annlz(sum(o90.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
        v70_HpvAge(a , :) = ...
            annlz(sum(sum(sum(sum(o70.newHpv(: , 2 , inds{i} , a , :),2),3),4),5)) ./ ...
            (annlz(sum(o70.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
        v50_HpvAge(a , :) = ...
            annlz(sum(sum(sum(sum(o50.newHpv(: , 2 , inds{i} , a , :),2),3),4),5)) ./ ...
            (annlz(sum(o50.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
    end
    
    % Perform age standardization
    vNo_Hpv = sum(bsxfun(@times , vNo_HpvAge , wVec));
    v90_Hpv = sum(bsxfun(@times , v90_HpvAge , wVec));
    v70_Hpv = sum(bsxfun(@times , v70_HpvAge , wVec));
    v50_Hpv = sum(bsxfun(@times , v50_HpvAge , wVec));
    
    figure()
    plot(tVec(1 : stepsPerYear : end) , vNo_Hpv , tVec(1 : stepsPerYear : end) , v90_Hpv , ...
        tVec(1 : stepsPerYear : end) , v70_Hpv , tVec(1 : stepsPerYear : end) , v50_Hpv)
    title([plotTits{i} , ' Incidence'])
    xlabel('Year'); ylabel('Incidence per 100,000')
    legend('No vaccination' , '90% coverage' , '70% coverage' , ...
        '50% coverage')
    % Reduction
    v90_Red = (v90_Hpv - vNo_Hpv) ./ vNo_Hpv * 100;
    v70_Red = (v70_Hpv - vNo_Hpv) ./ vNo_Hpv * 100;
    v50_Red = (v50_Hpv - vNo_Hpv) ./ vNo_Hpv * 100;
    
    figure()
    plot(tVec(1 : stepsPerYear : end) , v90_Red , tVec(1 : stepsPerYear : end) , v70_Red , ...
        tVec(1 : stepsPerYear : end) , v50_Red)
    title([plotTits{i} , ' Incidence Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('90% coverage' , '70% coverage' , '50% coverage')
    axis([tVec(2) tVec(end) -100 0])
    
    T = table(tVec(1 : stepsPerYear : end)' , v90_Hpv' , v70_Hpv' , v50_Hpv' , ...
        v90_Red' , v70_Red' , v50_Red');
    writetable(T , [files{i} , '_stand.csv'] , 'Delimiter' , ',')
end
%% HPV Prevalence
inds = {':' , 2 : 6 , 1 , 10};
files = {'General_HpvPrev_VaxCover' , 'HivAll_HpvPrev_VaxCover' , 'HivNegPrev_Hpv_VaxCover' , ...
    'ART_HpvPrev_VaxCover'};
plotTits = {'General HPV' , 'HIV-Positive + HPV' , 'HIV-Negative + HPV' , ...
    'HIV-Positive on ART'};
fac = 10 ^ 5;
vNo_HpvPrevAge = zeros(age , length(tVec));
v90_HpvPrevAge = vNo_HpvPrevAge;
v70_HpvPrevAge = vNo_HpvPrevAge;
v50_HpvPrevAge = vNo_HpvPrevAge;
fac = 10 ^ 5;
% ageTotal = zeros(length(tVec) , age);
% for a = 1 : age
%     ageTotal(: , a) = sum(popVec(1 : size(popVec , 1) , ...
%         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
%         2 , a , 1 : risk))), 2) + sum(popVec(1 : size(popVec , 1) , ...
%         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 8 : 10 , 1 : periods , ...
%         2 , a , 1 : risk))) , 2);
% end
for i = 1
    for a = 3 : age
        % all
        % general
        allFTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        % All HIV-positive women (not on ART)
        allHivFTot = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        % All HIV-negative women
        hivNegTot = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            2 , a , 1 : risk));
        % ART
        artTot = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            2 , a , 1 : risk));
        
        % HPV positive
        % general with HPV
        allF = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        % HIV-positive women (not on ART)with HPV
        allHivF = toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        % HIV-negative women with HPV
        hivNeg = toInd(allcomb(1 , 1 : viral , 2 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            2 , a , 1 : risk));
         % ART and HPV
        art = toInd(allcomb(10 , 6 , 2 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            2 , a , 1 : risk));
        
        
        
        infArray = {allF , allHivF , hivNeg , art};
        totArray = {allFTot, allHivFTot , hivNegTot , artTot};
        
        vNo_HpvPrevAge(a , :) = ...
            sum(sum(sum(sum(sum(sum(sum(sum(oNo.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
            sum(sum(sum(sum(sum(sum(sum(sum(oNo.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9);
        
        v90_HpvPrevAge(a , :) = ...
            sum(sum(sum(sum(sum(sum(sum(sum(o90.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
            sum(sum(sum(sum(sum(sum(sum(sum(o90.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9);
        
        v70_HpvPrevAge(a , :) = ...
            sum(sum(sum(sum(sum(sum(sum(sum(o70.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
            sum(sum(sum(sum(sum(sum(sum(sum(o70.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9);
        
        v50_HpvPrevAge(a , :) = ...
            sum(sum(sum(sum(sum(sum(sum(sum(o50.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
            sum(sum(sum(sum(sum(sum(sum(sum(o50.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9);
        
        
    end
    T_HpvPrevAge = table(tVec' ,  vNo_HpvPrevAge' , v90_HpvPrevAge' , ...
        v70_HpvPrevAge' , v50_HpvPrevAge');
    writetable(T_HpvPrevAge , [files{i} , '_stand.csv'] , 'Delimiter' , ',')
end
%%
figure()
ageGroup = {'10 - 14' , '15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 - 34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};

m1 = mesh(1 : age , tVec , vNo_HpvPrevAge' * 100);
set(m1 , 'edgecolor' , 'r')
alpha(0.1)
hold on
m2 = mesh(1 : age , tVec , v90_HpvPrevAge' * 100);
set(m2 , 'edgecolor' , 'b')
alpha(0.1)
set(gca , 'yLim' , [tVec(1) tVec(end)]);
set(gca , 'xtick' , 3 : age , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Year'); zlabel('Prevalence (%)')
title('HPV Prevalence')
legend('No Vaccination' , '90% coverage')

% Vax start year
% hold on
% [px , py] = meshgrid(1 : age , 2018 ...
%     .* ones(age , 1));
% pz = bsxfun(@times , ones(size(px , 1) , size(py , 1)) , linspace(0 , 100 , size(px , 1)));
% m = surf(px , py , pz' , 'edgecolor' , 'r');
% set(m , 'facecolor' , 'r')
% legend('No Vaccination' , '90% Coverage' , 'Vaccination Start')
% alpha(0.4)
%%
figure()
v90_HpvPrevRed = (v90_HpvPrevAge - vNo_HpvPrevAge) ./ vNo_HpvPrevAge .* 100;
mesh(1 : age , tVec , v90_HpvPrevRed')
set(gca , 'yLim' , [tVec(1) tVec(end)]);
set(gca , 'xtick' , 3 : age , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Year'); zlabel('Prevalence Reduction (%)')
title('HPV Prevalence Reduction (90% coverage)')
% Vax start year
% hold on
% [px , py] = meshgrid(1 : age , 2018 ...
%     .* ones(age , 1));
% pz = bsxfun(@times , ones(size(px , 1) , size(py , 1)) , linspace(0 , 100 , size(px , 1)));
% m = surf(px , py , pz' , 'edgecolor' , 'r');
% set(m , 'facecolor' , 'r')
% alpha(0.4)
yr = 2088; % 70 years post-vaccination start
t = (yr - 2018) * stepsPerYear;
v90_HpvPrevRed_2088 = sum(wVec(3 : end) .* v90_HpvPrevRed(3 : end , t));
%% CC associated deaths
inds = {':' , 2 : 6 , 1};
files = {'General_CCMortality_VaxCover' , 'HivAll_CCMortality_VaxCover' , 'HivNeg_CCMortality_VaxCover'};
plotTits = {'General Cervical Cancer' , 'HIV+ Cervical Cancer' , 'HIV- Cervical Cancer'};
fac = 10 ^ 5;
vNo_MortAge = zeros(ageGroups , length(tVec) / stepsPerYear);
v90_MortAge = vNo_MortAge;
v70_MortAge = vNo_MortAge;
v50_MortAge = vNo_MortAge;

for i = 1 : length(inds)
    for a = 1 : age
        % general
        allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        % All HIV-positive women (not on ART)
        allHivF = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        % All HIV-negative women
        hivNeg = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            2 , a , 1 : risk));
        
        genArray = {allF , allHivF , hivNeg};
        
        vNo_MortAge(a , :) = ...
            annlz(sum(sum(sum(oNo.ccDeath(: , inds{i} , : , : , a),2),3),4)) ./ ...
            (annlz(sum(oNo.popVec(: , genArray{i}) , 2)) ./ stepsPerYear) * fac;
        
        v90_MortAge(a , :) = ...
            annlz(sum(sum(sum(o90.ccDeath(: , inds{i} , : , : , a),2),3),4)) ./ ...
            (annlz(sum(o90.popVec(: , genArray{i}) , 2)) ./ stepsPerYear) * fac;
        
        v70_MortAge(a , :) = ...
            annlz(sum(sum(sum(o70.ccDeath(: , inds{i} , : , : , a),2),3),4)) ./ ...
            (annlz(sum(o70.popVec(: , genArray{i}) , 2)) ./ stepsPerYear) * fac;
        
        v50_MortAge(a , :) = ...
            annlz(sum(sum(sum(o50.ccDeath(:, inds{i} , : , : , a),2),3),4)) ./ ...
            (annlz(sum(o50.popVec(:, genArray{i}) , 2)) ./ stepsPerYear) * fac;
    end
    
    % Perform age standardization
    vNo_Mort = sum(bsxfun(@times , vNo_MortAge , wVec));
    v90_Mort = sum(bsxfun(@times , v90_MortAge , wVec));
    v70_Mort = sum(bsxfun(@times , v70_MortAge , wVec));
    v50_Mort = sum(bsxfun(@times , v50_MortAge , wVec));
    
    figure()
    plot(tVec(1 : stepsPerYear : end - 1) , vNo_Mort , tVec(1 : stepsPerYear : end - 1) , v90_Mort , ...
        tVec(1 : stepsPerYear : end - 1) , v70_Mort , tVec(1 : stepsPerYear : end - 1) , v50_Mort)
    title([plotTits{i} , ' Mortality'])
    xlabel('Year'); ylabel('Mortality per 100,000')
    legend('No vaccination' , '90% coverage' , '70% coverage' , ...
        '50% coverage')
    % Reduction
    v90_Red = (v90_Mort - vNo_Mort) ./ vNo_Mort * 100;
    v70_Red = (v70_Mort - vNo_Mort) ./ vNo_Mort * 100;
    v50_Red = (v50_Mort - vNo_Mort) ./ vNo_Mort * 100;
    
    figure()
    plot(tVec(1 : stepsPerYear : end - 1) , v90_Red , tVec(1 : stepsPerYear : end - 1) , v70_Red , ...
        tVec(1 : stepsPerYear : end - 1) , v50_Red)
    title([plotTits{i} , ' Mortality Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('90% coverage' , '70% coverage' , '50% coverage')
    axis([tVec(2) tVec(end) -100 0])
    
    T = table(tVec(1 : stepsPerYear : end - 1)' , v90_Mort' , v70_Mort' , v50_Mort' , ...
        v90_Red' , v70_Red' , v50_Red');
    writetable(T , [files{i} , '_stand.csv'] , 'Delimiter' , ',')
end
%% By CD4
dVec = [2 : 6 , 10];
tits = {'Acute' , 'CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , 'CD4 < 200' , ...
    'ART'};
filenames = {'AcuteMort' , 'CD4_500Mort' , 'CD4_500_350Mort' , 'CD4_350_200Mort' , ...
    'CD4_200Mort' , 'ARTMort'};

vNo_wNoMortAge = zeros(age , length(tVec) / stepsPerYear);

for d = 1 : length(dVec)
    for a = 1 : age
        hiv_ccSusAge = [toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , 1 : 8 , 1 : periods , ...
            2 , 4 : age , 1 : risk)) ; toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , ...
            9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];

        vNo_wNoMortAge(a , :) = ...
            annlz(sum(sum(sum(sum(oNo.ccDeath(:, dVec(d) , : , : , :),2),3),4),5)) ./ ...
            (annlz(sum(oNo.popVec(: , hiv_ccSusAge) , 2)) ./ stepsPerYear) * fac;

        v90_MortAge(a , :) = ...
            annlz(sum(sum(sum(sum(o90.ccDeath(: , dVec(d) , : , : , :),2),3),4),5)) ./ ...
            (annlz(sum(o90.popVec(: , hiv_ccSusAge) , 2)) ./ stepsPerYear) * fac;

        v70_MortAge(a , :) = ...
            annlz(sum(sum(sum(sum(o70.ccDeath(: , dVec(d) , : , : , :),2),3),4),5)) ./ ...
            (annlz(sum(o70.popVec(: , hiv_ccSusAge) , 2)) ./ stepsPerYear) * fac;

        v50_MortAge(a , :) = ...
            annlz(sum(sum(sum(sum(o50.ccDeath(: , dVec(d) , : , : , :),2),3),4),5)) ./ ...
            (annlz(sum(o50.popVec(: , hiv_ccSusAge) , 2)) ./ stepsPerYear) * fac;
    end
    
%     hiv_ccSus = sum(bsxfun(@times , hiv_ccSusAge , wVec));
    vNo_wNoMort = sum(bsxfun(@times , vNo_wNoMortAge , wVec));
    v90_Mort = sum(bsxfun(@times , v90_MortAge , wVec));
    v70_Mort = sum(bsxfun(@times , v70_MortAge , wVec));
    v50_Mort = sum(bsxfun(@times , v50_MortAge , wVec));
    
    figure(104)
    subplot(3 , 2 , d)
    plot(tVec(1 : stepsPerYear : end) , vNo_wNoMort , tVec(1 : stepsPerYear : end) , v90_Mort , ...
        tVec(1 : stepsPerYear : end) , v70_Mort , tVec(1 : stepsPerYear : end) , v50_Mort)
    title([tits{d} , ' Mortality'])
    xlabel('Year'); ylabel('Mortality per 100,000')
    legend('No vaccination' , '90% coverage' , '70% coverage' ,...
        '50% coverage' , 'Location' , 'northeastoutside')
    % Reduction
    v90_Red = (v90_Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v70_Red = (v70_Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v50_Red = (v50_Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;

    figure(105)
    subplot(3 , 2 ,d)
    plot(tVec(1 : stepsPerYear : end) , v90_Red , tVec(1 : stepsPerYear : end) , v70_Red , ...
        tVec(1 : stepsPerYear : end) , v50_Red)
    title([tits{d} , ' Mortality Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('90% coverage' , '70% coverage' , '50% coverage', ...
        'Location' , 'northeastoutside')
    axis([tVec(2) tVec(end) -100 0])

    T = table(tVec(1 : stepsPerYear : end)' , v90_Mort' , v70_Mort' , v50_Mort' , ...
        v90_Red' , v70_Red' , v50_Red');
    writetable(T , ['VaxCover_' , filenames{d} , '_stand.csv'] , 'Delimiter' , ',')
end
%% new cervical cancers
newCC_90 = o90.newCC;
newCC_70 = o70.newCC;
newCC_50 = o50.newCC;
newCC_0 = oNo.newCC;

%% populations
pop90 = o90.popVec;
pop70 = o70.popVec;
pop50 = o50.popVec;
popNo = oNo.popVec;
%%
% General Incidence
% susceptible population
pop90_susGen = zeros(age , length(tVec) / stepsPerYear);
pop70_susGen = pop90_susGen;
pop50_susGen = pop90_susGen;
popNo_susGen = pop90_susGen;
for a = 1 : age
    ccSusGen = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , 2 , a , 1 : risk));
        toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , 2 , a , 1 : risk))];
    pop90_susGen(a , :) = annlz(sum(pop90(: , ccSusGen) , 2)) ./ stepsPerYear;
    pop70_susGen(a , :) = annlz(sum(pop70(: , ccSusGen) , 2)) ./ stepsPerYear;
    pop50_susGen(a , :) = annlz(sum(pop50(: , ccSusGen) , 2)) ./ stepsPerYear;
    popNo_susGen(a , :) = annlz(sum(popNo(: , ccSusGen) , 2)) ./ stepsPerYear;
end
% cases in general population
genCC90Age = zeros(age , length(tVec) / stepsPerYear);
genCC70Age = genCC90Age;
genCC50Age = genCC90Age;
genCCNoAge = genCC90Age;
for a = 1 : age
    genCC90Age(a , :) = annlz(squeeze(sum(sum(sum(newCC_90(: , : ,: , a) , 2),3),4)));
    genCC70Age(a , :) = annlz(squeeze(sum(sum(sum(newCC_70(: , : ,: , a) , 2),3),4)));
    genCC50Age(a , :)  = annlz(squeeze(sum(sum(sum(newCC_50(: , : ,: , a) , 2),3),4)));
    genCCNoAge(a , :) = annlz(squeeze(sum(sum(sum(newCC_0(: , : ,: , a) , 2),3),4)));
end
%%
% incidence
fac = 10 ^ 5;
i_gen90 = sum(bsxfun(@times , genCC90Age ./ pop90_susGen , wVec)) .* fac;
i_gen70 = sum(bsxfun(@times , genCC70Age ./ pop70_susGen , wVec)) .* fac;
i_gen50 = sum(bsxfun(@times , genCC50Age ./ pop50_susGen , wVec)) .* fac;
i_genNo = sum(bsxfun(@times , genCCNoAge ./ popNo_susGen , wVec)) .* fac;
figure()
plot(tVec(1 : stepsPerYear : end) , i_gen90 , tVec(1 : stepsPerYear : end) , ...
    i_gen70 , tVec(1 : stepsPerYear : end) , i_gen50 , tVec(1 : stepsPerYear : end) , i_genNo)
title('General Cervical Cancer Incidence')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('90% coverage' , '70% coverage' , '50% coverage' , 'No vaccination')

% relative incidence reduction
figure()
genRelRed_90 = (i_gen90 - i_genNo) ./ i_genNo * 100;
genRelRed_70 = (i_gen70 - i_genNo) ./ i_genNo * 100;
genRelRed_50 = (i_gen50 - i_genNo) ./ i_genNo * 100;
plot(tVec(1 : stepsPerYear : end) , genRelRed_90 , tVec(1 : stepsPerYear : end) , genRelRed_70 , tVec(1 : stepsPerYear : end) , genRelRed_50)
title('General Cervical Cancer Reduction')
axis([tVec(1) , tVec(end) , -100 , 0])
xlabel('Year'); ylabel('Relative Difference (%)')
legend('90% coverage' , '70% coverage' , '50% coverage')
%% Export general incidence and reduction as csv
T = table(tVec(1 : stepsPerYear : end)' , i_gen90' , i_gen70' , i_gen50' , i_genNo' , genRelRed_90' , ...
    genRelRed_70' , genRelRed_50');
writetable(T , 'General_Incidence_stand.csv' , 'Delimiter' , ',')

%% HIV specific incidence
%% By CD4
% incidence
fac = 10 ^ 5;
figure()
tits = {'Acute' , 'CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , 'CD4 < 200' , ...
    'ART'};
filenames = {'Acute' , 'CD4_500' , 'CD4_500_350' , 'CD4_350_200' , 'CD4_200' , ...
    'ART'};
hiv_vec = [2 : 6 , 10];
hivCC90 = zeros(age , length(tVec) / stepsPerYear);
hivCC70 = hivCC90;
hivCC50 = hivCC90;
hivCCNo = hivCC90;
pop90_susHiv = zeros(age , length(tVec) / stepsPerYear);
pop70_susHiv = pop90_susHiv;
pop50_susHiv = pop90_susHiv;
popNo_susHiv = pop90_susHiv;
for i = 1 : length(hiv_vec)
    for a = 1 : age
        d = hiv_vec(i);
        hivCC90(a , :) = annlz(sum(sum(newCC_90(: , d , : , a),3),4))';
        hivCC70(a , :) = annlz(sum(sum(newCC_70(: , d , : , a),3),4))';
        hivCC50(a , :) = annlz(sum(sum(newCC_50(: , d , : , a),3),4))';
        hivCCNo(a , :) = annlz(sum(sum(newCC_0(: , d , : , a),3),4))';
        
        hivSus = [toInd(allcomb(d , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , a , 1 : risk)) ; toInd(allcomb(d , 1 : viral , 1 : hpvTypes , ...
            9 : 10 , 1 : periods , 2 , a , 1 : risk))];
        pop90_susHiv(a , :) = annlz(sum(pop90(: , hivSus) , 2)) ./ stepsPerYear;
        pop70_susHiv(a , :) = annlz(sum(pop70(: , hivSus) , 2)) ./ stepsPerYear;
        pop50_susHiv(a , :) = annlz(sum(pop50(: , hivSus) , 2)) ./ stepsPerYear;
        popNo_susHiv(a , :) = annlz(sum(popNo(: , hivSus) , 2)) ./ stepsPerYear;
    end
    
    % Perform age standardization
    h_gen90 = sum(bsxfun(@times , hivCC90 ./ pop90_susHiv , wVec)) .* fac;
    h_gen70 = sum(bsxfun(@times , hivCC70 ./ pop70_susHiv , wVec)) .* fac;
    h_gen50 = sum(bsxfun(@times , hivCC50 ./ pop50_susHiv , wVec)) .* fac;
    h_genNo = sum(bsxfun(@times , hivCCNo ./ popNo_susHiv , wVec)) .* fac;
    
    subplot(3 , 2 , i)
    plot(tVec(1 : stepsPerYear : end) , h_gen90 , tVec(1 : stepsPerYear : end) , h_gen70 , tVec(1 : stepsPerYear : end) , ...
        h_gen50 , tVec(1 : stepsPerYear : end) , h_genNo)
    title(['Cervical Cancer Incidence ', tits{i}])
    xlabel('Year'); ylabel('Incidence per 100,000')

    % Export HIV-positive incidence as csv
    T = table(tVec(1 : stepsPerYear : end)' , h_gen90' , h_gen70' , h_gen50' , h_genNo');
    writetable(T , [filenames{i} , '_Incidence_stand.csv'] , 'Delimiter' , ',')
end
legend('90% coverage' , '70% coverage' , '50% coverage' , 'No vaccination')

% relative incidence reduction
figure()
for i = 1 : length(hiv_vec)
    for a = 1 : age
        d = hiv_vec(i);
        hivCC90(a , :) = annlz(sum(sum(newCC_90(: , d , : , a),3),4))';
        hivCC70(a , :) = annlz(sum(sum(newCC_70(: , d , : , a),3),4))';
        hivCC50(a , :) = annlz(sum(sum(newCC_50(: , d , : , a),3),4))';
        hivCCNo(a , :) = annlz(sum(sum(newCC_0(: , d , : , a),3),4))';
        
        hivSus = [toInd(allcomb(d , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , a , 1 : risk)) ; toInd(allcomb(d , 1 : viral , 1 : hpvTypes , ...
            9 : 10 , 1 : periods , 2 , a , 1 : risk))];
        pop90_susHiv(a , :) = annlz(sum(pop90(: , hivSus) , 2)) ./ stepsPerYear;
        pop70_susHiv(a , :) = annlz(sum(pop70(: , hivSus) , 2)) ./ stepsPerYear;
        pop50_susHiv(a , :) = annlz(sum(pop50(: , hivSus) , 2)) ./ stepsPerYear;
        popNo_susHiv(a , :) = annlz(sum(popNo(: , hivSus) , 2)) ./ stepsPerYear;
    end
    h_gen90 = sum(bsxfun(@times , hivCC90 ./ pop90_susHiv , wVec)).* fac;
    h_gen70 = sum(bsxfun(@times , hivCC70 ./ pop70_susHiv , wVec)) .* fac;
    h_gen50 = sum(bsxfun(@times , hivCC50 ./ pop50_susHiv , wVec)) .* fac;
    h_genNo = sum(bsxfun(@times , hivCCNo ./ popNo_susHiv , wVec)) .* fac;

    hivRelRed_90 = (h_gen90 - h_genNo) ./ h_genNo * 100;
    hivRelRed_70 = (h_gen70 - h_genNo) ./ h_genNo * 100;
    hivRelRed_50 = (h_gen50 - h_genNo) ./ h_genNo * 100;

    subplot(3 , 2 , i)
    plot(tVec(1 : stepsPerYear : end) , hivRelRed_90 , tVec(1 : stepsPerYear : end) , hivRelRed_70 , tVec(1 : stepsPerYear : end) , hivRelRed_50)
    title(['Cervical Cancer Reduction ', tits{i}])
    xlabel('Year'); ylabel('Relative Difference (%)')
    axis([tVec(1) , tVec(end) , -100 , 0])
    % Export HIV-positive reduction as csv
    T = table(tVec(1 : stepsPerYear : end)' , hivRelRed_90' , hivRelRed_70' , hivRelRed_50');
    writetable(T , [filenames{i} , '_Incidence_stand.csv'] , 'Delimiter' , ',')
end
legend('90% coverage' , '70% coverage' , '50% coverage')

%% Acute and CD4 > 500

fac = 10 ^ 5;
figure()
tit = 'Acute and CD4 > 500';
filename = 'Acute_CD4_500';
hivCC90 = zeros(age , length(tVec) / stepsPerYear);
hivCC70 = hivCC90;
hivCC50 = hivCC90;
hivCCNo = hivCC90;
vec = [2 : 3];

for a = 1 : age
    hivCC90(a , :) = annlz(sum(sum(sum(newCC_90(: , vec , : , a),2),3),4));
    hivCC70(a , :) = annlz(sum(sum(sum(newCC_70(: , vec , : , a),2),3),4));
    hivCC50(a , :) = annlz(sum(sum(sum(newCC_50(: , vec , : , a),2),3),4));
    hivCCNo(a , :) = annlz(sum(sum(sum(newCC_0(: , vec , : , a),2),3),4));
    
    hivSus = [toInd(allcomb(d , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)) ; toInd(allcomb(vec , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , a , 1 : risk))];
    pop90_susHiv(a , :) = annlz(sum(pop90(: , hivSus) , 2)) ./ stepsPerYear;
    pop70_susHiv(a , :) = annlz(sum(pop70(: , hivSus) , 2)) ./ stepsPerYear;
    pop50_susHiv(a , :) = annlz(sum(pop50(: , hivSus) , 2)) ./ stepsPerYear;
    popNo_susHiv(a , :) = annlz(sum(popNo(: , hivSus) , 2)) ./ stepsPerYear;
end

h_gen90 = sum(bsxfun(@times , hivCC90 ./ pop90_susHiv , wVec)) .* fac;
h_gen70 = sum(bsxfun(@times , hivCC70 ./ pop70_susHiv , wVec)) .* fac;
h_gen50 = sum(bsxfun(@times , hivCC50 ./ pop50_susHiv , wVec)) .* fac;
h_genNo = sum(bsxfun(@times , hivCCNo ./ popNo_susHiv , wVec)) .* fac;

plot(tVec(1 : stepsPerYear : end) , h_gen90 , tVec(1 : stepsPerYear : end) , h_gen70 , tVec(1 : stepsPerYear : end) , ...
    h_gen50 , tVec(1 : stepsPerYear : end) , h_genNo)
title(['Cervical Cancer Incidence ', tit])
xlabel('Year'); ylabel('Incidence per 100,000')

% Export HIV-positive incidence as csv
T = table(tVec(1 : stepsPerYear : end)' , h_gen90' , h_gen70' , h_gen50' , h_genNo');
writetable(T , [filename , '_Incidence_stand.csv'] , 'Delimiter' , ',')

legend('90% coverage' , '70% coverage' , '50% coverage' , 'No vaccination')

% relative incidence reduction
figure()
for a = 1 : age
    hivCC90(a , :) = annlz(sum(sum(sum(newCC_90(: , vec , : , a),2),3),4));
    hivCC70(a , :) = annlz(sum(sum(sum(newCC_70(: , vec , : , a),2),3),4));
    hivCC50(a , :) = annlz(sum(sum(sum(newCC_50(: , vec , : , a),2),3),4));
    hivCCNo(a , :) = annlz(sum(sum(sum(newCC_0(: , vec , : , a),2),3),4));

    hivSus = [toInd(allcomb(vec , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)) ; toInd(allcomb(vec , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , a , 1 : risk))];
    pop90_susHiv(a , :) = annlz(sum(pop90(: , hivSus) , 2));
    pop70_susHiv(a , :) = annlz(sum(pop70(: , hivSus) , 2));
    pop50_susHiv(a , :) = annlz(sum(pop50(: , hivSus) , 2));
    popNo_susHiv(a , :) = annlz(sum(popNo(: , hivSus) , 2));
end

h_gen90 = sum(bsxfun(@times , hivCC90 ./ pop90_susHiv , wVec)) .* fac;
h_gen70 = sum(bsxfun(@times , hivCC70 ./ pop70_susHiv , wVec)) .* fac;
h_gen50 = sum(bsxfun(@times , hivCC50 ./ pop50_susHiv , wVec)) .* fac;
h_genNo = sum(bsxfun(@times , hivCCNo ./ popNo_susHiv , wVec)) .* fac;

hivRelRed_90 = (h_gen90 - h_genNo) ./ h_genNo * 100;
hivRelRed_70 = (h_gen70 - h_genNo) ./ h_genNo * 100;
hivRelRed_50 = (h_gen50 - h_genNo) ./ h_genNo * 100;

plot(tVec(1 : stepsPerYear : end) , hivRelRed_90 , tVec(1 : stepsPerYear : end) , hivRelRed_70 , tVec(1 : stepsPerYear : end) , hivRelRed_50)
title(['Cervical Cancer Reduction ', tit])
xlabel('Year'); ylabel('Relative Difference (%)')
axis([tVec(1) , tVec(end) , -100 , 0])
% Export HIV-positive reduction as csv
T = table(tVec(1 : stepsPerYear : end)' , hivRelRed_90' , hivRelRed_70' , hivRelRed_50');
writetable(T , [filename , '_Incidence_stand.csv'] , 'Delimiter' , ',')

legend('90% coverage' , '70% coverage' , '50% coverage')

%% Aggregate (without ART)
hivAllCC90 = zeros(age , length(tVec) / stepsPerYear);
hivAllCC70 = hivAllCC90;
hivAllCC50 = hivAllCC90;
hivAllCCNo = hivAllCC90;
for a = 1 : age
    hivAllCC90(a , :) = annlz(sum(sum(sum(newCC_90(: , 2 : 6 , : , a),2),3),4));
    hivAllCC70(a , :) = annlz(sum(sum(sum(newCC_70(: , 2 : 6 , : , a),2),3),4));
    hivAllCC50(a , :) = annlz(sum(sum(sum(newCC_50(: , 2 : 6 , : , a),2),3),4));
    hivAllCCNo(a , :) = annlz(sum(sum(sum(newCC_0(: , 2 : 6 , : , a),2),3),4));
    
    hivAllSus = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)); toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes ,...
        8 : 10 , 1 : periods , 2 , a , 1 : risk))];
    pop90_susHiv(a , :) = annlz(sum(pop90(: , hivAllSus) , 2)) ./ stepsPerYear;
    pop70_susHiv(a , :) = annlz(sum(pop70(: , hivAllSus) , 2)) ./ stepsPerYear;
    pop50_susHiv(a , :) = annlz(sum(pop50(: , hivAllSus) , 2)) ./ stepsPerYear;
    popNo_susHiv(a , :) = annlz(sum(popNo(: , hivAllSus) , 2)) ./ stepsPerYear;
end
hivAllInc90 = sum(bsxfun(@times , hivAllCC90 ./ pop90_susHiv , wVec))* fac;
hivAllInc70 = sum(bsxfun(@times , hivAllCC70 ./ pop70_susHiv , wVec)) * fac;
hivAllInc50 = sum(bsxfun(@times , hivAllCC50 ./ pop50_susHiv , wVec)) * fac;
hivAllIncNo = sum(bsxfun(@times , hivAllCCNo ./ popNo_susHiv , wVec)) * fac;
figure()

plot(tVec(1 : stepsPerYear : end) , hivAllInc90 , tVec(1 : stepsPerYear : end) , hivAllInc70 , ...
    tVec(1 : stepsPerYear : end) , hivAllInc50 , tVec(1 : stepsPerYear : end) , hivAllIncNo)
title('Cervical Cancer Incidence Among All HIV+')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('90% coverage' , '70% coverage' , '50% coverage' , 'No vaccination')

hivAllRed_90 = (hivAllInc90 - hivAllIncNo) ./ hivAllIncNo * 100;
hivAllRed_70 = (hivAllInc70 - hivAllIncNo) ./ hivAllIncNo * 100;
hivAllRed_50 = (hivAllInc50 - hivAllIncNo) ./ hivAllIncNo * 100;
figure()
plot(tVec(1 : stepsPerYear : end) , hivAllRed_90 , tVec(1 : stepsPerYear : end) , hivAllRed_70 , ...
    tVec(1 : stepsPerYear : end) , hivAllRed_50)
title('Cervical Cancer Reduction Among All HIV+')
xlabel('Year'); ylabel('Reduction (%)')
legend('90% coverage' , '70% coverage' , '50% coverage')
axis([tVec(1) , tVec(end) , -100 , 0])
% Export HIV+ aggregate cervical cancer incidence and reduction as csv
T = table(tVec(1 : stepsPerYear : end)' , hivAllInc90' , hivAllInc70' , hivAllInc50' , ...
    hivAllIncNo' , hivAllRed_90' , hivAllRed_70' , hivAllRed_50');
    writetable(T , 'AllHiv_Incidence.csv' , 'Delimiter' , ',')
%% HIV-
% incidence
figure()
hivNegCC90 = zeros(age , length(tVec) / stepsPerYear);
hivNegCC70 = hivNegCC90;
hivNegCC50 = hivNegCC90;
hivNegCCNo = hivNegCC90;
pop90_susHivNeg = zeros(age , length(tVec) / stepsPerYear);
pop70_susHivNeg = pop90_susHivNeg;
pop50_susHivNeg = pop90_susHivNeg;
popNo_susHivNeg = pop90_susHivNeg;
for a = 1 : age
    hivNegCC90(a , :) = annlz(sum(sum(newCC_90(: , 1 , : , a),3),4));
    hivNegCC70(a , :) = annlz(sum(sum(newCC_70(: , 1 , : , a),3),4));
    hivNegCC50(a , :) = annlz(sum(sum(newCC_50(: , 1 , : , a),3),4));
    hivNegCCNo(a , :) = annlz(sum(sum(newCC_0(: , 1 , : , a),3),4));

    hivNegSus = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)) ; toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , a , 1 : risk))];
    pop90_susHivNeg(a , :) = annlz(sum(pop90(: , hivNegSus) , 2)) ./ stepsPerYear; 
    pop70_susHivNeg(a , :) = annlz(sum(pop70(: , hivNegSus) , 2)) ./ stepsPerYear;
    pop50_susHivNeg(a , :) = annlz(sum(pop50(: , hivNegSus) , 2)) ./ stepsPerYear;
    popNo_susHivNeg(a , :) = annlz(sum(popNo(: , hivNegSus) , 2)) ./ stepsPerYear;
end

h_neg90 = sum(bsxfun(@times , hivNegCC90 ./ pop90_susHivNeg , wVec)) .* fac;
h_neg70 = sum(bsxfun(@times , hivNegCC70 ./ pop70_susHivNeg , wVec)) .* fac;
h_neg50 = sum(bsxfun(@times , hivNegCC50 ./ pop50_susHivNeg , wVec)) .* fac;
h_negNo = sum(bsxfun(@times , hivNegCCNo ./ popNo_susHivNeg , wVec)) .* fac;
plot(tVec(1 : stepsPerYear : end) , h_neg90 , tVec(1 : stepsPerYear : end) , h_neg70 , tVec(1 : stepsPerYear : end) , ...
    h_neg50 , tVec(1 : stepsPerYear : end) , h_negNo)
title('Cervical Cancer Incidence in HIV-')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('90% coverage' , '70% coverage' , '50% coverage' , 'No vaccination')
% Export HIV-negative incidence as csv
T = table(tVec(1 : stepsPerYear : end)' , h_neg90' , h_neg70' , h_neg50' , h_negNo');
writetable(T , 'HIVNeg_Incidence_stand.csv' , 'Delimiter' , ',')

% relative incidence reduction
figure()
hivNegRelRed_90 = (h_neg90 - h_negNo) ./ h_negNo * 100;
hivNegRelRed_70 = (h_neg70 - h_negNo) ./ h_negNo * 100;
hivNegRelRed_50 = (h_neg50 - h_negNo) ./ h_negNo * 100;

plot(tVec(1 : stepsPerYear : end) , hivNegRelRed_90 , tVec(1 : stepsPerYear : end) , hivNegRelRed_70 , ...
    tVec(1 : stepsPerYear : end) , hivNegRelRed_50)
title('Cervical Cancer Reduction among HIV-')
xlabel('Year'); ylabel('Relative Difference (%)')
axis([tVec(1) , tVec(end) , -100 , 0])
legend('90% coverage' , '70% coverage' , '50% coverage')

% Export HIV-negative reduction as csv
T = table(tVec(1 : stepsPerYear : end)' , hivNegRelRed_90' , hivNegRelRed_70' , hivNegRelRed_50');
writetable(T , 'HIVNeg_Reduction_stand.csv' , 'Delimiter' , ',')

%% Waning
v90_w20 = load('H:\HHCoM_Results\Vax_0.9_wane_20.mat');
v90_w15 = load('H:\HHCoM_Results\Vax_0.9_wane_15.mat');
v90_w10 = load('H:\HHCoM_Results\Vax_0.9_wane_10.mat');
v90_w0 = load('H:\HHCoM_Results\Vax_0.9_wane_0.mat');
v0_w0 = load('H:\HHCoM_Results\Vax_0_wane_0.mat');

%% Deaths

inds = {':' , 2 : 6 , 1};
files = {'General_CCMortality' , 'allHiv_CCMortality' , 'hivNegMortality'};
plotTits = {'General Cervical Cancer' , 'HIV+ Cervical Cancer' , 'HIV- Cervical Cancer'};
v90_w20MortAge = vNo_wNoMortAge;
v90_w15MortAge = vNo_wNoMortAge;
v90_w10MortAge = vNo_wNoMortAge;
v90_wNoMortAge = vNo_wNoMortAge;
for i = 1 : length(inds)
    for a = 1 : age
        % general
        allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        % All HIV-positive women (not on ART)
        allHivF = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , 2 , a , 1 : risk));
        % All HIV-negative women
        hivNeg = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            2 , a , 1 : risk));
        
        genArray = {allF , allHivF , hivNeg};
               
        vNo_wNoMortAge(a , :) = ...
            annlz(sum(sum(sum(v0_w0.ccDeath(: , inds{i} , : , a),2),3),4)) ./ ...
            (annlz(sum(v0_w0.popVec(: , genArray{i}) , 2)) / stepsPerYear) * fac;
        
        v90_w20MortAge(a , :) = ...
            annlz(sum(sum(sum(v90_w20.ccDeath(: , inds{i} , : , a),2),3),4)) ./ ...
            (annlz(sum(v90_w20.popVec(: , genArray{i}) , 2)) / stepsPerYear) * fac;
        
        v90_w15MortAge(a , :) = ...
            annlz(sum(sum(sum(v90_w15.ccDeath(:, inds{i} , : , a),2),3),4)) ./ ...
            (annlz(sum(v90_w15.popVec(:, genArray{i}) , 2)) / stepsPerYear) * fac;
        
        v90_w10MortAge(a , :) = ...
            annlz(sum(sum(sum(v90_w10.ccDeath(: , inds{i} , : , a),2),3),4)) ./ ...
            (annlz(sum(v90_w10.popVec(: , genArray{i}) , 2)) / stepsPerYear) * fac;
        
        v90_wNoMortAge(a , :) = ...
            annlz(sum(sum(sum(v90_w0.ccDeath(: , inds{i} , : , a),2),3),4)) ./ ...
            (annlz(sum(v90_w0.popVec(: , genArray{i}) , 2)) / stepsPerYear) * fac;
    end
    
    vNo_wNoMort = sum(bsxfun(@times , vNo_wNoMortAge , wVec));
    v90_w20Mort = sum(bsxfun(@times , v90_w20MortAge , wVec));
    v90_w15Mort = sum(bsxfun(@times , v90_w15MortAge , wVec));
    v90_w10Mort = sum(bsxfun(@times , v90_w10MortAge , wVec));
    v90_wNoMort = sum(bsxfun(@times , v90_wNoMortAge , wVec));
    
    figure()
    plot(tVec(1 : stepsPerYear : end) , vNo_wNoMort , tVec(1 : stepsPerYear : end) , v90_wNoMort , ...
        tVec(1 : stepsPerYear : end) , v90_w20Mort , tVec(1 : stepsPerYear : end) , v90_w15Mort , ...
        tVec(1 : stepsPerYear : end) , v90_w10Mort)
    title([plotTits{i} , ' Mortality'])
    xlabel('Year'); ylabel('Mortality per 100,000')
    legend('No vaccination' , 'No Waning' , '20 years' , '15 years' , '10 years')
    % Reduction
    v90_wNoRed = (v90_wNoMort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w20Red = (v90_w20Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w15Red = (v90_w15Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w10Red = (v90_w10Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    
    figure()
    plot(tVec(1 : stepsPerYear : end) , v90_wNoRed , tVec(1 : stepsPerYear : end) , v90_w20Red , ...
        tVec(1 : stepsPerYear : end) , v90_w15Red , tVec(1 : stepsPerYear : end) , v90_w10Red)
    title([plotTits{i} , ' Mortality Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('No Waning' , '20 years' , '15 years' , '10 years')
    axis([tVec(2) tVec(end) -100 0])
    
    T = table(tVec(1 : stepsPerYear : end)' , v90_w20Mort' , v90_w15Mort' , v90_w10Mort' , ...
        v90_wNoMort' , v90_wNoRed' , v90_w20Red' , v90_w15Red' , v90_w10Red');
    writetable(T , ['waning_', files{i} , '_stand.csv'] , 'Delimiter' , ',')
end
%%
% By CD4

dVec = [2 : 6 , 10];
tits = {'Acute' , 'CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , 'CD4 < 200' , ...
    'ART'};
filenames = {'AcuteMort' , 'CD4_500Mort' , 'CD4_500_350Mort' , 'CD4_350_200Mort' , ...
    'CD4_200Mort' , 'ARTMort'};

for d = 1 : length(dVec)
    for a = 1 : age
    hiv_ccSus = [toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , 1 : 8 , 1 : periods , ...
        2 , a , 1 : risk)) ; toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , a , 1 : risk))];

    vNo_wNoMortAge(a , :) = ...
        annlz(sum(sum(sum(v0_w0.ccDeath(: , dVec(d) , : , a),2),3),4)) ./ ...
        (annlz(sum(v0_w0.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;

    v90_w20MortAge(a , :) = ...
        annlz(sum(sum(sum(v90_w20.ccDeath(: , dVec(d) , : , a),2),3),4)) ./ ...
        (annlz(sum(v90_w20.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;

    v90_w15MortAge(a , :) = ...
        annlz(sum(sum(sum(v90_w15.ccDeath(: , dVec(d) , : , a),2),3),4)) ./ ...
        (annlz(sum(v90_w15.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;

    v90_w10MortAge(a , :) = ...
        annlz(sum(sum(sum(v90_w10.ccDeath(:, dVec(d) , : , a),2),3),4)) ./ ...
        (annlz(sum(v90_w10.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;

    v90_wNoMortAge(a , :) = ...
        annlz(sum(sum(sum(v90_w0.ccDeath(: , dVec(d) , : , a),2),3),4)) ./ ...
        (annlz(sum(v90_w0.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    end
    vNo_wNoMort = sum(bsxfun(@times , vNo_wNoMortAge , wVec));
    v90_w20Mort = sum(bsxfun(@times , v90_w20MortAge , wVec));
    v90_w15Mort = sum(bsxfun(@times , v90_w15MortAge , wVec));
    v90_w10Mort = sum(bsxfun(@times , v90_w10MortAge , wVec));
    v90_wNoMort = sum(bsxfun(@times , v90_wNoMortAge , wVec));
    figure(102)
    subplot(3 , 2 , d)
    plot(tVec(1 : stepsPerYear : end) , vNo_wNoMort , tVec(1 : stepsPerYear : end) , v90_wNoMort , ...
        tVec(1 : stepsPerYear : end) , v90_w20Mort , tVec(1 : stepsPerYear : end) , v90_w15Mort , ...
        tVec(1 : stepsPerYear : end) , v90_w10Mort)
    title([tits{d} , ' Mortality'])
    xlabel('Year'); ylabel('Mortality per 100,000')
    legend('No vaccination' , 'No Waning' , '20 years' , '15 years' , '10 years' ,...
        'Location' , 'northeastoutside')
    % Reduction
    v90_wNoRed = (v90_wNoMort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w20Red = (v90_w20Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w15Red = (v90_w15Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
    v90_w10Red = (v90_w10Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;

    figure(103)
    subplot(3 , 2 ,d)
    plot(tVec(1 : stepsPerYear : end) , v90_wNoRed , tVec(1 : stepsPerYear : end) , v90_w20Red , ...
        tVec(1 : stepsPerYear : end) , v90_w15Red , tVec(1 : stepsPerYear : end) , v90_w10Red)
    title([tits{d} , ' Mortality Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('No Waning' , '20 years' , '15 years' , '10 years' , ...
        'Location' , 'northeastoutside')
    axis([tVec(2) tVec(end) -100 0])

    T = table(tVec(1 : stepsPerYear : end)' , v90_w20Mort' , v90_w15Mort' , v90_w10Mort' , ...
        v90_wNoMort' , v90_wNoRed' , v90_w20Red' , v90_w15Red' , v90_w10Red');
    writetable(T , [filenames{d} , '_stand.csv'] , 'Delimiter' , ',')
end
%% Acute and CD4 > 500
vec = [2 : 3];
tit = 'Acute and CD4 > 500';
filename = 'Acute_CD4_500Mort';
for a = 1 : age
    hiv_ccSus = [toInd(allcomb(vec , 1 : viral , 1 : hpvTypes , 1 : 8 , 1 : periods , ...
        2 , a , 1 : risk)) ; toInd(allcomb(vec , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , a , 1 : risk))];
    
    vNo_wNoMortAge(a , :) = ...
        annlz(sum(sum(sum(v0_w0.ccDeath(: , vec , : , a),2),3),4)) ./ ...
        (annlz(sum(v0_w0.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    
    v90_w20MortAge(a , :) = ...
        annlz(sum(sum(sum(v90_w20.ccDeath(: , vec , : , a),2),3),4)) ./ ...
        (annlz(sum(v90_w20.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    
    v90_w15MortAge(a , :) = ...
        annlz(sum(sum(sum(v90_w15.ccDeath(: , vec , : , a),2),3),4)) ./ ...
        (annlz(sum(v90_w15.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    
    v90_w10MortAge(a , :) = ...
        annlz(sum(sum(sum(v90_w10.ccDeath(: , vec , : , a),2),3),4)) ./ ...
        (annlz(sum(v90_w10.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    
    v90_wNoMortAge(a , :) = ...
        annlz(sum(sum(sum(v90_w0.ccDeath(: , vec , : , a),2),3),4)) ./ ...
        (annlz(sum(v90_w0.popVec(:, hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
end

vNo_wNoMort = sum(bsxfun(@times , vNo_wNoMortAge , wVec));
v90_w20Mort = sum(bsxfun(@times , v90_w20MortAge , wVec));
v90_w15Mort = sum(bsxfun(@times , v90_w15MortAge , wVec));
v90_w10Mort = sum(bsxfun(@times , v90_w10MortAge , wVec));
v90_wNoMort = sum(bsxfun(@times , v90_wNoMortAge , wVec));

figure(109)
plot(tVec(1 : stepsPerYear : end) , vNo_wNoMort , tVec(1 : stepsPerYear : end) , v90_wNoMort , ...
    tVec(1 : stepsPerYear : end) , v90_w20Mort , tVec(1 : stepsPerYear : end) , v90_w15Mort , ...
    tVec(1 : stepsPerYear : end) , v90_w10Mort)
title([tit , ' Mortality'])
xlabel('Year'); ylabel('Mortality per 100,000')
legend('No vaccination' , 'No Waning' , '20 years' , '15 years' , '10 years' ,...
    'Location' , 'northeastoutside')
% Reduction
v90_wNoRed = (v90_wNoMort - vNo_wNoMort) ./ vNo_wNoMort * 100;
v90_w20Red = (v90_w20Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
v90_w15Red = (v90_w15Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;
v90_w10Red = (v90_w10Mort - vNo_wNoMort) ./ vNo_wNoMort * 100;

figure(110)
plot(tVec(1 : stepsPerYear : end) , v90_wNoRed , tVec(1 : stepsPerYear : end) , v90_w20Red , ...
    tVec(1 : stepsPerYear : end) , v90_w15Red , tVec(1 : stepsPerYear : end) , v90_w10Red)
title([tit , ' Mortality Reduction'])
xlabel('Year'); ylabel('Reduction (%)')
legend('No Waning' , '20 years' , '15 years' , '10 years' , ...
    'Location' , 'northeastoutside')
axis([tVec(2) tVec(end) -100 0])

T = table(tVec(1 : stepsPerYear : end)' , v90_w20Mort' , v90_w15Mort' , v90_w10Mort' , ...
    v90_wNoMort' , v90_wNoRed' , v90_w20Red' , v90_w15Red' , v90_w10Red');
writetable(T , [filename , '_stand.csv'] , 'Delimiter' , ',')

%%
% General susceptibles
v90_w0_incAge = zeros(age , length(tVec) / stepsPerYear);
v90_w20_incAge = v90_w0_incAge;
v90_w15_incAge = v90_w0_incAge;
v90_w10_incAge = v90_w0_incAge;

for a = 1 : age
    ccSus = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)) ; toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , a , 1 : risk))];
    v90_w0_incAge(a , :)  = annlz(sum(sum(sum(sum(v90_w0.newCC(: , : , : , a)...
        ,2),3),4),5)) ./ (annlz(sum(v90_w0.popVec(: , ccSus) , 2)) ./ stepsPerYear) * fac;
    v90_w20_incAge(a , :)  = annlz(sum(sum(sum(sum(v90_w20.newCC(: , : , : , a)...
        ,2),3),4),5)) ./ (annlz(sum(v90_w20.popVec(: , ccSus) , 2)) ./ stepsPerYear) * fac;
    v90_w15_incAge(a , :)  = annlz(sum(sum(sum(sum(v90_w15.newCC(: , : , : , a)...
        ,2),3),4),5)) ./ (annlz(sum(v90_w15.popVec(: , ccSus) , 2)) ./ stepsPerYear) * fac;
    v90_w10_incAge(a , :)  = annlz(sum(sum(sum(sum(v90_w10.newCC(: , : , : , a)...
        ,2),3),4),5)) ./ (annlz(sum(v90_w10.popVec(: , ccSus) , 2)) ./ stepsPerYear) * fac;
end

v90_w0_inc = sum(bsxfun(@times , v90_w0_incAge , wVec));
v90_w20_inc = sum(bsxfun(@times , v90_w20_incAge , wVec));
v90_w15_inc = sum(bsxfun(@times , v90_w15_incAge , wVec));
v90_w10_inc = sum(bsxfun(@times , v90_w10_incAge , wVec));

figure()
plot(tVec(1 : stepsPerYear : end) , i_genNo , tVec(1 : stepsPerYear : end) , v90_w0_inc , tVec(1 : stepsPerYear : end) , v90_w20_inc , ...
    tVec(1 : stepsPerYear : end) , v90_w15_inc , tVec(1 : stepsPerYear : end) , v90_w10_inc)
title('Vaccine Waning Period and General Cervical Cancer Incidence')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('No vaccination' , 'No waning' , '20 years' , '15 years' , '10 years')

% Relative reduction
wNo = (v90_w0_inc - i_genNo) ./ i_genNo * 100;
w20 = (v90_w20_inc - i_genNo) ./ i_genNo * 100;
w15 = (v90_w15_inc - i_genNo) ./ i_genNo * 100;
w10 = (v90_w10_inc - i_genNo) ./ i_genNo * 100;
figure()
plot(tVec(1 : stepsPerYear : end) , wNo , tVec(1 : stepsPerYear : end) , w20 , ...
    tVec(1 : stepsPerYear : end) , w15 , tVec(1 : stepsPerYear : end) , w10)
axis([tVec(2) tVec(end) -100 0])
title('Vaccine Waning Period and Reduction of General Cervical Cancer Incidence')
xlabel('Year'); ylabel('Reduction (%)')
legend('No waning' , '20 years' , '15 years' , '10 years')
% Export general incidence/reduction as csv
T = table(tVec(1 : stepsPerYear : end)' , v90_w0_inc' , v90_w20_inc' , v90_w15_inc' , v90_w10_inc' , ...
    wNo' , w20' , w15' , w10');
writetable(T , 'GenInc_Waning_stand.csv' , 'Delimiter' , ',')
%% HIV-negative
% incidence

v90_w0_incNegAge = zeros(age , length(tVec) / stepsPerYear);
v90_w20_incNegAge = v90_w0_incNegAge;
v90_w15_incNegAge = v90_w0_incNegAge;
v90_w10_incNegAge = v90_w0_incNegAge;

for a = 1 : age
    ccNegSus = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)) ; toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , a , 1 : risk))];
    v90_w0_incNegAge(a , :)  = annlz(sum(sum(sum(v90_w0.newCC(: , 1 , : , a)...
        ,2),3),4)) ./ (annlz(sum(v90_w0.popVec(: , ccNegSus) , 2)) ./ stepsPerYear) * fac;
    v90_w20_incNegAge(a , :)  = annlz(sum(sum(sum(v90_w20.newCC(: , 1 , : , a)...
        ,2),3),4)) ./ (annlz(sum(v90_w20.popVec(: , ccNegSus) , 2)) ./ stepsPerYear) * fac;
    v90_w15_incNegAge(a , :)  = annlz(sum(sum(sum(v90_w15.newCC(: , 1 , : , a)...
        ,2),3),4)) ./ (annlz(sum(v90_w15.popVec(: , ccNegSus) , 2)) ./ stepsPerYear) * fac;
    v90_w10_incNegAge(a , :)  = annlz(sum(sum(sum(v90_w10.newCC(: , 1 , : , a)...
        ,2),3),4)) ./ (annlz(sum(v90_w10.popVec(: , ccNegSus) , 2)) ./ stepsPerYear) * fac;
end

v90_w0_incNeg = sum(bsxfun(@times , v90_w0_incNegAge , wVec));
v90_w20_incNeg = sum(bsxfun(@times , v90_w20_incNegAge , wVec));
v90_w15_incNeg = sum(bsxfun(@times , v90_w15_incNegAge , wVec));
v90_w10_incNeg = sum(bsxfun(@times , v90_w10_incNegAge , wVec));

figure()
plot(tVec(1 : stepsPerYear : end) , h_negNo , tVec(1 : stepsPerYear : end) , v90_w0_incNeg , tVec(1 : stepsPerYear : end) ,...
    v90_w20_incNeg , tVec(1 : stepsPerYear : end) , v90_w15_incNeg , tVec(1 : stepsPerYear : end) , v90_w10_incNeg)
title('Vaccine Waning Period and HIV- Cervical Cancer Incidence')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('No vaccination' , 'No waning' , '20 years' , '15 years' , '10 years')

% Relative reduction
wNoNeg = (v90_w0_incNeg - h_negNo) ./ h_negNo * 100;
w20Neg = (v90_w20_incNeg - h_negNo) ./ h_negNo * 100;
w15Neg = (v90_w15_incNeg - h_negNo) ./ h_negNo * 100;
w10Neg = (v90_w10_incNeg - h_negNo) ./ h_negNo * 100;
figure()
plot(tVec(1 : stepsPerYear : end) , wNoNeg , tVec(1 : stepsPerYear : end) , w20Neg , ...
    tVec(1 : stepsPerYear : end) , w15Neg , tVec(1 : stepsPerYear : end) , w10Neg)
axis([tVec(2) tVec(end) -100 0])
title('Vaccine Waning Period and Reduction of HIV- Cervical Cancer Incidence')
xlabel('Year'); ylabel('Reduction (%)')
legend('No waning' , '20 years' , '15 years' , '10 years')
% Export general incidence/reduction as csv
T = table(tVec(1 : stepsPerYear : end)' , v90_w0_incNeg' , v90_w20_incNeg' , v90_w15_incNeg' , ...
    v90_w10_incNeg' , wNoNeg' , w20Neg' , w15Neg' , w10Neg');
writetable(T , 'HIVNeg_Waning_stand.csv' , 'Delimiter' , ',')
%% HIV-positive
%% Aggregated
v90_w0_incHiv = zeros(age , length(tVec) / stepsPerYear);
v90_w20_incHiv = v90_w0_incHiv;
v90_w15_incHiv = v90_w0_incHiv;
v90_w10_incHiv = v90_w0_incHiv;
for a = 1 : age
    % Incidence
    hiv_ccSus = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , a , 1 : risk)) ; toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , a , 1 : risk))];
    v90_w0_incHiv(a , :)  = annlz(sum(sum(sum(v90_w0.newCC(: , 2 : 6 , : , a)...
        ,2),3),4)) ./ (annlz(sum(v90_w0.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    v90_w20_incHiv(a , :)  = annlz(sum(sum(sum(v90_w20.newCC(: , 2 : 6 , : , a)...
        ,2),3),4)) ./ (annlz(sum(v90_w20.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    v90_w15_incHiv(a , :)  = annlz(sum(sum(sum(v90_w15.newCC(: , 2 : 6 , : , a)...
        ,2),3),4)) ./ (annlz(sum(v90_w15.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear)* fac;
    v90_w10_incHiv(a , :)  = annlz(sum(sum(sum(v90_w10.newCC(: , 2 : 6 , : , a)...
        ,2),3),4)) ./ (annlz(sum(v90_w10.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear)* fac;
end

v90_w0_incHiv = sum(bsxfun(@times , v90_w0_incHiv , wVec));
v90_w20_incHiv = sum(bsxfun(@times , v90_w20_incHiv , wVec));
v90_w15_incHiv = sum(bsxfun(@times , v90_w15_incHiv , wVec));
v90_w10_incHiv = sum(bsxfun(@times , v90_w10_incHiv , wVec));

figure()
plot(tVec(1 : stepsPerYear : end) , hivAllIncNo , tVec(1 : stepsPerYear : end) , v90_w0_incHiv , tVec(1 : stepsPerYear : end) , ...
    v90_w20_incHiv , tVec(1 : stepsPerYear : end) , v90_w15_incHiv , tVec(1 : stepsPerYear : end) , v90_w10_incHiv)
title('Vaccine Waning Period and General Cervical Cancer Incidence in HIV+ (90% coverage)')
xlabel('Year'); ylabel('Incidence per 100,000')
legend('No vaccination' , 'No waning' , '20 years' , '15 years' , '10 years' , ...
    'Location' , 'northeastoutside')

% Relative reduction
wNoRed = (v90_w0_incHiv - hivAllIncNo) ./ hivAllIncNo * 100;
w20Red = (v90_w20_incHiv - hivAllIncNo) ./ hivAllIncNo * 100;
w15Red = (v90_w15_incHiv - hivAllIncNo) ./ hivAllIncNo * 100;
w10Red = (v90_w10_incHiv - hivAllIncNo) ./ hivAllIncNo * 100;
figure()
plot(tVec(1 : stepsPerYear : end) , wNoRed , tVec(1 : stepsPerYear : end) , w20Red , tVec(1 : stepsPerYear : end) , ...
    w15Red , tVec(1 : stepsPerYear : end) , w10Red)
title('Vaccine Waning Period and General Cervical Cancer Reduction in HIV+')
xlabel('Year'); ylabel('Reduction(%)')
axis([tVec(2) tVec(end) -100 0 ]);
legend('No waning' , '20 years' , '15 years' , '10 years' , ...
    'Location' , 'northeastoutside')
% Export general incidence/reduction as csv
T = table(tVec(1 : stepsPerYear : end)' ,hivAllIncNo' , v90_w0_incHiv' , v90_w20_incHiv' , ...
    v90_w15_incHiv' , v90_w10_incHiv' , wNoRed' , w20Red' , w15Red' , w10Red');
writetable(T , 'AllHivInc_Waning_stand.csv' , 'Delimiter' , ',')

%% By CD4
dVec = [2 : 6 , 10];
tits = {'Acute' , 'CD4 > 500' , 'CD4 500-350' , 'CD4 350-200' , 'CD4 < 200' , ...
    'ART'};
filenames = {'Acute' , 'CD4_500' , 'CD4_500_350' , 'CD4_350_200' , 'CD4_200' , ...
    'ART'};
v90_wNo_incCD4 = zeros(age , length(tVec) / stepsPerYear);
v90_w20_incCD4 = v90_wNo_incCD4;
v90_w15_incCD4 = v90_wNo_incCD4;
v90_w10_incCD4 = v90_wNo_incCD4;
hivCCNo = v90_w10_incCD4;

for d = 1 : length(dVec)
    for a = 1 : age
        % Incidence
        hiv_ccSus = [toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , a , 1 : risk)) ; toInd(allcomb(dVec(d) , 1 : viral , 1 : hpvTypes , ...
            9 : 10 , 1 : periods , 2 , a , 1 : risk))];
        v90_wNo_incCD4(a , :)  = annlz(sum(sum(v90_w0.newCC(: , dVec(d) , : , a)...
            ,3),4)) ./ (annlz(sum(v90_w0.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
        v90_w20_incCD4(a , :)  = annlz(sum(sum(v90_w20.newCC(: , dVec(d) , : , a)...
            ,3),4)) ./ (annlz(sum(v90_w20.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
        v90_w15_incCD4(a , :)  = annlz(sum(sum(v90_w15.newCC(: , dVec(d) , : , a)...
            ,3),4)) ./ (annlz(sum(v90_w15.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
        v90_w10_incCD4(a , :)  = annlz(sum(sum(v90_w10.newCC(: , dVec(d) , : , a)...
            ,3),4)) ./ (annlz(sum(v90_w10.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
        hivCCNo(a , :) = annlz(sum(sum(newCC_0(: , dVec(d) , : , a),3),4)) ...
            ./ (annlz(sum(popNo(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
        
    end
    
    v90_wNo_incCD4 = sum(bsxfun(@times , v90_wNo_incCD4 , wVec));
    v90_w20_incCD4 = sum(bsxfun(@times , v90_w20_incCD4 , wVec));
    v90_w15_incCD4 = sum(bsxfun(@times , v90_w15_incCD4 , wVec));
    v90_w10_incCD4 = sum(bsxfun(@times , v90_w10_incCD4 , wVec));
    
    % base (no vaccine)
    h_genNo = sum(bsxfun(@times , hivCCNo , wVec));
    %     hivCCNo = sum(sum(sum(newCC_0(: , dVec(d) , : , : , :),3),4),5);
    %     popNo_susHiv = sum(popNo(1 : end - 1 , hiv_ccSus) , 2);
    %     h_genNo = hivCCNo(2 : end) ./ popNo_susHiv .* fac;
    
    figure(100)
    subplot(3 , 2 , d)
    plot(tVec(1 : stepsPerYear : end) , h_genNo , tVec(1 : stepsPerYear : end) , v90_wNo_incCD4 , tVec(1 : stepsPerYear : end) , ...
        v90_w20_incCD4 , tVec(1 : stepsPerYear : end) , v90_w15_incCD4 , tVec(1 : stepsPerYear : end) , v90_w10_incCD4)
    title(['Vaccine Waning Period and CC Incidence in ' , tits{d}])
    xlabel('Year'); ylabel('Incidence per 100,000')
    legend('No vaccination' , 'No waning' , '20 years' , '15 years' , '10 years' , ...
        'Location' , 'northeastoutside')
    
    % Relative reduction
    wNoRedCD4 = (v90_wNo_incCD4 - h_genNo) ./ h_genNo * 100;
    w20RedCD4 = (v90_w20_incCD4 - h_genNo) ./ h_genNo * 100;
    w15RedCD4 = (v90_w15_incCD4 - h_genNo) ./ h_genNo * 100;
    w10RedCD4 = (v90_w10_incCD4 - h_genNo) ./ h_genNo * 100;
    figure(101)
    subplot(3 , 2 , d)
    plot(tVec(1 : stepsPerYear : end) , wNoRedCD4 , tVec(1 : stepsPerYear : end) , w20RedCD4 , tVec(1 : stepsPerYear : end) , ...
        w15RedCD4 , tVec(1 : stepsPerYear : end) , w10RedCD4)
    axis([tVec(2) tVec(end) -100 0 ]);
    title(['Vaccine Waning Period and CC Reduction in ' tits{d}])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('No waning' , '20 years' , '15 years' , '10 years' , ...
        'Location' , 'northeastoutside')
    % Export general incidence/reduction as csv
    T = table(tVec(1 : stepsPerYear : end)' , v90_wNo_incCD4' , v90_w20_incCD4' , ...
        v90_w15_incCD4' , v90_w10_incCD4' , wNoRedCD4' , w20RedCD4' , w15RedCD4' , ...
        w10RedCD4');
    writetable(T , [filenames{d} , '_Waning_stand.csv'] , 'Delimiter' , ',')
end

%% Acute and CD4 > 500
vec = [2 : 3];
tit = 'Acute and  CD4 > 500';
filename= 'Acute_500';
% Incidence
for a = 1 : age
    hiv_ccSus = [toInd(allcomb(vec , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
        2 , 4 : age , 1 : risk)) ; toInd(allcomb(vec , 1 : viral , 1 : hpvTypes , ...
        9 : 10 , 1 : periods , 2 , 4 : age , 1 : risk))];
    v90_wNo_incCD4(a , :)  = annlz(sum(sum(sum(sum(v90_w0.newCC(: , vec , : , :)...
        ,2),3),4),5)) ./ (annlz(sum(v90_w0.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    v90_w20_incCD4(a , :)  = annlz(sum(sum(sum(sum(v90_w20.newCC(: , vec , : , :)...
        ,2),3),4),5)) ./ (annlz(sum(v90_w20.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    v90_w15_incCD4(a , :)  = annlz(sum(sum(sum(sum(v90_w15.newCC(: , vec , : , :)...
        ,2),3),4),5)) ./ (annlz(sum(v90_w15.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    v90_w10_incCD4(a , :)  = annlz(sum(sum(sum(sum(v90_w10.newCC(: , vec , : , :)...
        ,2),3),4),5)) ./ (annlz(sum(v90_w10.popVec(: , hiv_ccSus) , 2)) ./ stepsPerYear) * fac;
    % base (no vaccine)
    hivCCNo = annlz(sum(sum(sum(sum(newCC_0(: , vec , : , :), 2), 3),4),5));
    popNo_susHiv = annlz(sum(popNo(: , hiv_ccSus) , 2)) ./ stepsPerYear;
    h_genNo = hivCCNo ./ popNo_susHiv .* fac;
end

v90_wNo_incCD4 = sum(bsxfun(@times , v90_wNo_incCD4, wVec));
v90_w20_incCD4 = sum(bsxfun(@times , v90_w20_incCD4, wVec));
v90_w15_incCD4 = sum(bsxfun(@times , v90_w15_incCD4, wVec));
v90_w10_incCD4 = sum(bsxfun(@times , v90_w10_incCD4, wVec));

figure(111)
plot(tVec(1 : stepsPerYear : end) , h_genNo , tVec(1 : stepsPerYear : end) , v90_wNo_incCD4 , tVec(1 : stepsPerYear : end) , ...
    v90_w20_incCD4 , tVec(1 : stepsPerYear : end) , v90_w15_incCD4 , tVec(1 : stepsPerYear : end) , v90_w10_incCD4)
title(['Vaccine Waning Period and CC Incidence in ' , tit])
xlabel('Year'); ylabel('Incidence per 100,000')
legend('No vaccination' , 'No waning' , '20 years' , '15 years' , '10 years' , ...
    'Location' , 'northeastoutside')

% Relative reduction
wNoRedCD4 = (v90_wNo_incCD4 - h_genNo) ./ h_genNo * 100;
w20RedCD4 = (v90_w20_incCD4 - h_genNo) ./ h_genNo * 100;
w15RedCD4 = (v90_w15_incCD4 - h_genNo) ./ h_genNo * 100;
w10RedCD4 = (v90_w10_incCD4 - h_genNo) ./ h_genNo * 100;
figure(101)
plot(tVec(1 : stepsPerYear : end) , wNoRedCD4 , tVec(1 : stepsPerYear : end) , w20RedCD4 , tVec(1 : stepsPerYear : end) , ...
    w15RedCD4 , tVec(1 : stepsPerYear : end) , w10RedCD4)
axis([tVec(2) tVec(end) -100 0 ]);
title(['Vaccine Waning Period and CC Reduction in ' tit])
xlabel('Year'); ylabel('Reduction (%)')
legend('No waning' , '20 years' , '15 years' , '10 years' , ...
    'Location' , 'northeastoutside')
% Export general incidence/reduction as csv
T = table(tVec(1 : stepsPerYear : end)' , v90_wNo_incCD4' , v90_w20_incCD4' , ...
    v90_w15_incCD4' , v90_w10_incCD4' , wNoRedCD4 , w20RedCD4 , w15RedCD4 , ...
    w10RedCD4);
writetable(T , [filename , '_Waning_stand.csv'] , 'Delimiter' , ',')

