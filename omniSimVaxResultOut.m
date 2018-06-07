function omniSimVaxResultOut()
%%
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'])
noV = load('H:\HHCoM_Results\VaxCover_0_Eff_0.9.mat');
c90_2vFull = load('H:\HHCoM_Results\VaxCover_0.9_Eff_0.7.mat'); 
c70_2vFull = load('H:\HHCoM_Results\VaxCover_0.7_Eff_0.7.mat'); 
c70_9vFull = load('H:\HHCoM_Results\VaxCover_0.7_Eff_0.9.mat'); 
c70_2vPartial = load('H:\HHCoM_Results\VaxCover_0.7_Eff_0.63.mat');
c90_9vFull = load('H:\HHCoM_Results\VaxCover_0.9_Eff_0.9.mat'); 
c90_2vPartial = load('H:\HHCoM_Results\VaxCover_0.9_Eff_0.63.mat');
tVec = c90_2vFull.tVec;
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 
% midMat = zeros(stepsPerYear , size(o90.popVec , 1) / stepsPerYear);
% midMat(1 , :) = 1;
% midMat(end , :) = 1;
% midAnn = @(x) sum(midMat .* reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) / 2;
c = fix(clock);
currYear = c(1); % get the current year
%% Plot Settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)
%% HPV prevalence
ageGroups = age - 4 + 1;

inds = {':' , 2 : 6 , 1 , 10};
files = {'General_HpvPrev_VaxCover' , 'HivAll_HpvPrev_VaxCover' , 'HivNegPrev_Hpv_VaxCover' , ...
    'ART_HpvPrev_VaxCover'};
plotTits = {'General HPV' , 'HIV-Positive + HPV' , 'HIV-Negative + HPV' , ...
    'HIV-Positive on ART'};
vNo_HpvPrevAge = zeros(age , length(noV.tVec));
c70_2vPartial_HpvPrevAge = vNo_HpvPrevAge;
c90_2vFull_HpvPrevAge = vNo_HpvPrevAge;
c90_9vFull_HpvPrevAge = vNo_HpvPrevAge;
gen = {'Male' , 'Female'};
for g = 1 : 2
    for i = 1 : length(files)
        for a = 3 : age
            % all
            % general
            allFTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
                1 : periods , g , a , 1 : risk));
            % All HIV-positive women
            allHivFTot = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
                1 : periods , g , a , 1 : risk));...
                toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
                g , a , 1 : risk))];
            % All HIV-negative women
            hivNegTot = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
                g , a , 1 : risk));
            % ART
            artTot = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
                g , a , 1 : risk));
            
            % HPV positive
            % general with HPV
            allF = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 9 , ...
                1 : periods , g , a , 1 : risk));
            % HIV-positive women with HPV
            allHivF = [toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvTypes , 1 : 9 , ...
                1 : periods , g , a , 1 : risk)); ...
                toInd(allcomb(10 , 6 , 2 : hpvTypes , 1 : 9 , 1 : periods , ...
                g , a , 1 : risk))];
            % HIV-negative women with HPV
            hivNeg = toInd(allcomb(1 , 1 : viral , 2 : hpvTypes , 1 : 9 , 1 : periods , ...
                g , a , 1 : risk));
            % ART and HPV
            art = toInd(allcomb(10 , 6 , 2 : hpvTypes , 1 : 9 , 1 : periods , ...
                g , a , 1 : risk));
            
            
            
            infArray = {allF , allHivF , hivNeg , art};
            totArray = {allFTot, allHivFTot , hivNegTot , artTot};
            
            vNo_HpvPrevAge(a , :) = ...
                sum(sum(sum(sum(sum(sum(sum(sum(noV.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
                sum(sum(sum(sum(sum(sum(sum(sum(noV.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9);
            
            c70_2vPartial_HpvPrevAge(a , :) = ...
                sum(sum(sum(sum(sum(sum(sum(sum(c70_2vPartial.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
                sum(sum(sum(sum(sum(sum(sum(sum(c70_2vPartial.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9);
            
            c90_2vFull_HpvPrevAge(a , :) = ...
                sum(sum(sum(sum(sum(sum(sum(sum(c90_2vFull.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
                sum(sum(sum(sum(sum(sum(sum(sum(c90_2vFull.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9);
            
            c90_9vFull_HpvPrevAge(a , :) = ...
                sum(sum(sum(sum(sum(sum(sum(sum(c90_9vFull.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
                sum(sum(sum(sum(sum(sum(sum(sum(c90_9vFull.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9);
            
            
        end
        T_HpvPrevAge = table(tVec' ,  vNo_HpvPrevAge' , c70_2vPartial_HpvPrevAge' , ...
            c90_2vFull_HpvPrevAge' , c90_9vFull_HpvPrevAge');
        writetable(T_HpvPrevAge , [gen{g} , files{i} , '_stand.csv'] , 'Delimiter' , ',')
    end
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
m2 = mesh(1 : age , tVec , c70_2vPartial_HpvPrevAge' * 100);
set(m2 , 'edgecolor' , 'b')
alpha(0.1)
set(gca , 'yLim' , [tVec(1) tVec(end)]);
set(gca , 'xtick' , 3 : age , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Year'); zlabel('Prevalence (%)')
title('HPV Prevalence in General Female Population')
legend('No Vaccination' , '70% coverage (Bivalent, 90% effective)')

%% HPV prevalence reduction by age
hpvPrevRed_c70_2vPartial = ...
    (c70_2vPartial_HpvPrevAge - vNo_HpvPrevAge) ./ vNo_HpvPrevAge;
figure()
mesh(1 : age , tVec , hpvPrevRed_c70_2vPartial' * 100)
set(gca , 'yLim' , [tVec(1) tVec(end)]);
set(gca , 'xtick' , 3 : age , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Year'); zlabel('Prevalence Reduction (%)')
title('HPV Prevalence Reduction')

%% Overall HPV prevalence and reduction
files = {'General_HpvPrev' , 'HivAll_HpvPrev' , 'HivNegPrev_Hpv' , ...
    'ART_HpvPrev'};
gen = {'Male' , 'Female'};
for g = 1 : 2
    for i = 1 : length(files)
            % all
            % general
            allFTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
                1 : periods , g , 4 : age , 1 : risk));
            % All HIV-positive women
            allHivFTot = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
                1 : periods , 2 , 4 : age , 1 : risk)); ...
                toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
                g , 4 : age , 1 : risk))];
            % All HIV-negative women
            hivNegTot = toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
                g , 4 : age , 1 : risk));
            % ART
            artTot = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
                g , 4 : age , 1 : risk));

            % HPV positive
            % general with HPV
            allF = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 9 , ...
                1 : periods , g , 4 : age , 1 : risk));
            % HIV-positive women with HPV
            allHivF = [toInd(allcomb(2 : 6 , 1 : viral , 2 : hpvTypes , 1 : 9 , ...
                1 : periods , g , 4 : age , 1 : risk)); ...
                toInd(allcomb(10 , 6 , 2 : hpvTypes , 1 : 9 , 1 : periods , ...
                g , 4 : age , 1 : risk))];
            % HIV-negative women with HPV
            hivNeg = toInd(allcomb(1 , 1 : viral , 2 : hpvTypes , 1 : 9 , 1 : periods , ...
                g , 4 : age , 1 : risk));
             % ART and HPV
            art = toInd(allcomb(10 , 6 , 2 : hpvTypes , 1 : 9 , 1 : periods , ...
                g , 4 : age , 1 : risk));



            infArray = {allF , allHivF , hivNeg , art};
            totArray = {allFTot, allHivFTot , hivNegTot , artTot};

            vNo_HpvPrev = ...
                sum(sum(sum(sum(sum(sum(sum(sum(noV.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
                sum(sum(sum(sum(sum(sum(sum(sum(noV.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9) * 100;

            c70_2vPartial_HpvPrev =...
                sum(sum(sum(sum(sum(sum(sum(sum(c70_2vPartial.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
                sum(sum(sum(sum(sum(sum(sum(sum(c70_2vPartial.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9) * 100;

            c90_2vFull_HpvPrev = ...
                sum(sum(sum(sum(sum(sum(sum(sum(c90_2vFull.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
                sum(sum(sum(sum(sum(sum(sum(sum(c90_2vFull.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9) * 100;

            c90_9vFull_HpvPrev= ...
                sum(sum(sum(sum(sum(sum(sum(sum(c90_9vFull.popVec(: , infArray{i}) ,2),3),4),5),6),7),8),9) ./ ...
                sum(sum(sum(sum(sum(sum(sum(sum(c90_9vFull.popVec(: , totArray{i}) ,2),3),4),5),6),7),8),9) * 100;

            c70_2vPartial_HpvPrevRed = (c70_2vPartial_HpvPrev - vNo_HpvPrev) ./ vNo_HpvPrev .* 100;

            c90_2vFull_HpvPrevRed = (c90_2vFull_HpvPrev - vNo_HpvPrev) ./ vNo_HpvPrev .* 100;

            c90_9vFull_HpvPrevRed = (c90_9vFull_HpvPrev - vNo_HpvPrev) ./ vNo_HpvPrev .* 100;

        T_HpvPrevAge = table(tVec' , c70_2vPartial_HpvPrevRed , ...
            c90_2vFull_HpvPrevRed , c90_9vFull_HpvPrevRed);
        writetable(T_HpvPrevAge , [files{i} , gen{g}, '_red.csv'] , 'Delimiter' , ',')
        
        
        noVaxing = vNo_HpvPrev(1 : stepsPerYear : end);
        vax70_2vPartial = c70_2vPartial_HpvPrev(1 : stepsPerYear :end);
        vax90_9vFull = c90_9vFull_HpvPrev(1 : stepsPerYear :end);
        vax90_2vFull =  c90_2vFull_HpvPrev(1 : stepsPerYear :end);
        
        T_HpvPrev = table(tVec(1 : stepsPerYear :end)' , noVaxing , vax70_2vPartial , ...
            vax90_2vFull , vax90_9vFull);
        writetable(T_HpvPrev , [files{i} , gen{g}, '.csv'] , 'Delimiter' , ',')

        figure()
        plot(tVec , vNo_HpvPrev , tVec , c70_2vPartial_HpvPrev , tVec , c90_2vFull_HpvPrev , ...
            tVec , c90_9vFull_HpvPrev)
        xlabel('Year'); ylabel('Prevalence (%)')
        title(['HPV Prevalence Among ' , gen{g}])
        legend('No Coverage' , '70% Coverage (Partial 2v)' , ...
            '90% coverage (Full 2v)' , '90% coverage (Full 9v)')

        figure()
        plot(tVec , c70_2vPartial_HpvPrevRed , tVec , c90_2vFull_HpvPrevRed , ...
            tVec , c90_9vFull_HpvPrevRed)
        xlabel('Year'); ylabel('Prevalence Reduction (%)')
        title(['HPV Prevalence Reduction Among ' , gen{g}])
        legend('70% Coverage (Partial 2v)' , ...
            '90% coverage (Full 2v)' , '90% coverage (Full 9v)')
    end
end
%%
figure()
ageGroup = {'10 - 14' , '15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 - 34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};

m1 = mesh(1 : age , tVec , c90_2vFull_HpvPrevAge' * 100);
set(m1 , 'edgecolor' , 'r')
alpha(0.1)
hold on
m2 = mesh(1 : age , tVec , c90_9vFull_HpvPrevAge' * 100);
set(m2 , 'edgecolor' , 'b')
alpha(0.1)
set(gca , 'yLim' , [tVec(1) tVec(end)]);
set(gca , 'xtick' , 3 : age , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Year'); zlabel('Prevalence (%)')
title('HPV Prevalence')
legend('90% coverage (2v, 90% effective)' , '90% coverage (9v, 90% effective)')

%% HPV Incidence
inds = {':' , 2 : 6 , 1 , 10};
files = {'General_Hpv_Inc' , 'HivAll_Hpv_Inc' , 'HivNeg_Hpv_Inc' , 'ART_HPV_Inc'};
plotTits = {'General HPV' , 'HIV-Positive + HPV' , 'HIV-Negative + HPV' , ...
    'HIV-Positive on ART'};
fac = 10 ^ 2;
noV_Hpv = zeros(1 , length(tVec) / stepsPerYear);
% c70_2vPartial_Inc = noV_Hpv;
% c90_9vFullInc = noV_HpvAge;
% c90_2vFullInc = noV_HpvAge;
% v90_2vFullInc = noV_HpvAge;
gen = {'Males' , 'Females'};
for g = 1 : gender
    for i = 1 : length(inds)
        % general
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 : hpvStates , ...
            1 : periods ,g , 4 : 10 , 1 : risk));...
            toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 9 : 10 , ...
            1 : periods ,g , 4 : 10 , 1 : risk))];
        % All HIV-positive  (not on ART)
        allHivF = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 : hpvStates , ...
            1 : periods , g , 4 : 10 , 1 : risk));
            toInd(allcomb(2 : 6 , 1 : viral , 2 : 4 , 9 : 10 , ...
            1 : periods , g , 4 : 10 , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 , 1 : hpvStates , ...
            1 : periods , g, 4 : 10 , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 2 : 4 , 9 : 10 , ...
            1 : periods , g, 4 : 10 , 1 : risk))];
        % All HIV-negative
        hivNeg = [toInd(allcomb(1 , 1 : viral , 1 , 1 : hpvStates , 1 : periods , ...
            g , 4 : 10 , 1 : risk)); ...
            toInd(allcomb(1 , 1 : viral , 2 : 4 , 9 : 10 , 1 : periods , ...
            g , 4 : 10 , 1 : risk))];
        %  on ART
        artF = [toInd(allcomb(10 , 6 , 1 , 1 : hpvStates , ...
            1 : periods , g, 4 : 10 , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 2 : 4 , 9 : 10 , ...
            1 : periods , g, 4 : 10 , 1 : risk))];

        genArray = {allF , allHivF , hivNeg , artF};

        noV_Hpv = ...
            annlz(sum(sum(sum(sum(noV.newHpv(: , g , inds{i} , 4 : 10 , :),2),3),4),5)) ./ ...
            (annlz(sum(noV.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac;

        c70_2vPartial_Inc = ...
            annlz(sum(sum(sum(sum(c70_2vPartial.newHpv(: , g , inds{i} , 4 : 10 , :),2),3),4),5)) ./ ...
            (annlz(sum(c70_2vPartial.popVec(: , genArray{i}) , 2)) ./ stepsPerYear) * fac;

        c90_9vFullInc = ...
            annlz(sum(sum(sum(sum(c90_9vFull.newHpv(: , g , inds{i} , 4 : 10 , :),2),3),4),5)) ./ ...
            (annlz(sum(c90_9vFull.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;

        c90_2vFullInc = ...
            annlz(sum(sum(sum(sum(c90_2vFull.newHpv(: , g , inds{i} , 4 : 10 , :),2),3),4),5)) ./ ...
            (annlz(sum(c90_2vFull.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;

        figure()
        plot(tVec(1 : stepsPerYear : end) , noV_Hpv , tVec(1 : stepsPerYear : end) , c70_2vPartial_Inc , ...
            tVec(1 : stepsPerYear : end) , c90_9vFullInc , tVec(1 : stepsPerYear : end) , c90_2vFullInc)
        title([plotTits{i} , ' Incidence'])
        xlabel('Year'); ylabel('Incidence per 100')
        legend('No vaccination' , '70% Coverage (Partial 2v)' , ...
            '90% coverage (Full 9v)' ,...
            '90% coverage (Full 2v)')
        % Reduction
        c90_2vFull_Red = (c90_2vFullInc - noV_Hpv) ./ noV_Hpv * 100;
        c90_9vFull_Red = (c90_9vFullInc - noV_Hpv) ./ noV_Hpv * 100;
        c70_2vPartial_Red = (c70_2vPartial_Inc - noV_Hpv) ./ noV_Hpv * 100;

        figure()
        plot(tVec(1 : stepsPerYear : end) , c90_2vFull_Red , ...
            tVec(1 : stepsPerYear : end) , c90_9vFull_Red , ...
            tVec(1 : stepsPerYear : end) , c70_2vPartial_Red)
        title([plotTits{i} , ' Incidence Reduction'])
        xlabel('Year'); ylabel('Reduction (%)')
        legend('90% coverage (Full 2v)' , '90% coverage (Full 9v)'...
            , '70% Coverage (Partial 2v)')
        axis([tVec(2) tVec(end) -100 0])

        T = table(tVec(1 : stepsPerYear : end)' , noV_Hpv', c70_2vPartial_Inc' , ...
            c90_9vFullInc' , c90_2vFullInc' , ...
            c90_2vFull_Red' , c90_9vFull_Red' , c70_2vPartial_Red');
        writetable(T , [gen{g} , files{i} , '.csv'] , 'Delimiter' , ',')
    end
end

%% CC Incidence
inds = {':' , 2 : 6 , 1 , 10};
files = {'CC_General_Hpv_VaxCover' , 'CC_HivAll_Hpv_VaxCover' , ...
    'CC_HivNeg_Hpv_VaxCover' , 'CC_ART_HPV_VaxCover'};
plotTits = {'General' , 'HIV-Positive' , 'HIV-Negative' , ...
    'HIV-Positive on ART'};
fac = 10 ^ 5;
noV_Hpv = zeros(1 , length(tVec) / stepsPerYear);
% c70_2vPartial_Inc = noV_Hpv;
% c90_9vFullInc = noV_HpvAge;
% c90_2vFullInc = noV_HpvAge;
% v90_2vFullInc = noV_HpvAge;
for i = 1 : length(inds)
        % general
        allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 4 : age , 1 : risk))];
        % All HIV-positive women 
        allHivF = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 4 : age , 1 : risk));...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 4 : age , 1 : risk))];
        % All HIV-negative women
        hivNeg = [toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 1 : 4 , 1 : periods , ...
            2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(1 , 1 : viral , 1 : hpvTypes , 9 : 10 , 1 : periods , ...
            2 , 4 : age , 1 : risk))];
        % Women on ART
        artF = [toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : 4 , ...
            1 : periods , 2 , 4 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 9 : 10 , ...
            1 : periods , 2 , 4 : age , 1 : risk))];
        
        genArray = {allF , allHivF , hivNeg , artF};
        
        noV_Hpv = ...
            annlz(sum(sum(sum(noV.newCC(: , inds{i} , : , 4 : 10),2),3),4)) ./ ...
            (annlz(sum(noV.popVec(: , genArray{i}) , 2) ./ stepsPerYear))* fac;
       
        c70_2vPartial_Inc = ...
            annlz(sum(sum(sum(c70_2vPartial.newCC(: , inds{i} , : , 4 : 10),2),3),4)) ./ ...
           	(annlz(sum(c70_2vPartial.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
        c90_9vFullInc = ...
            annlz(sum(sum(sum(c90_9vFull.newCC(: , inds{i} , : , 4 : 10),2),3),4)) ./ ...
            (annlz(sum(c90_9vFull.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
        c90_2vFullInc = ...
            annlz(sum(sum(sum(c90_2vFull.newCC(: , inds{i} , : , 4 : 10),2),3),4)) ./ ...
            (annlz(sum(c90_2vFull.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac;
        
    figure()
    plot(tVec(1 : stepsPerYear : end) , noV_Hpv , tVec(1 : stepsPerYear : end) , c70_2vPartial_Inc , ...
        tVec(1 : stepsPerYear : end) , c90_9vFullInc , tVec(1 : stepsPerYear : end) , c90_2vFullInc)
    title([plotTits{i} , ' Cervical Cancer Incidence'])
    xlabel('Year'); ylabel('Incidence per 100,000')
    legend('No vaccination' , '70% Coverage (Partial 2v)' , ...
        '90% coverage (Full 9v)' ,...
        '90% coverage (Full 2v)')
    % Reduction
    c90_2vFull_Red = (c90_2vFullInc - noV_Hpv) ./ noV_Hpv * 100;
    c90_9vFull_Red = (c90_9vFullInc - noV_Hpv) ./ noV_Hpv * 100;
    c70_2vPartial_Red = (c70_2vPartial_Inc - noV_Hpv) ./ noV_Hpv * 100;
    
    figure()
    plot(tVec(1 : stepsPerYear : end) , c90_2vFull_Red , ...
        tVec(1 : stepsPerYear : end) , c90_9vFull_Red , ...
        tVec(1 : stepsPerYear : end) , c70_2vPartial_Red)
    title([plotTits{i} , ' Cervical Cancer Incidence Reduction'])
    xlabel('Year'); ylabel('Reduction (%)')
    legend('90% coverage (Full 2v)' , '90% coverage (Full 9v)'...
        , '70% Coverage (Partial 2v)')
    axis([tVec(2) tVec(end) -100 0])
    
    T = table(tVec(1 : stepsPerYear : end)' , c70_2vPartial_Inc' , ...
        c90_9vFullInc' , c90_2vFullInc' , ...
        c90_2vFull_Red' , c90_9vFull_Red' , c70_2vPartial_Red');
    writetable(T , [files{i} , '.csv'] , 'Delimiter' , ',')
end


%%
figure()
artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 1 : gender , 4 : 10 , 1 : risk));
artPop = sum(noV.popVec(: , artInds) , 2);
hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : gender , 4 : 10 , 1 : risk));
popTot = noV.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
    1 : hpvStates , 1 : periods , 1 : 2 , 4 : 10, 1 : risk)));
hivPop = sum(noV.popVec(: , hivInds) , 2);
hiv_art = [100 * hivPop ./ sum(popTot , 2), 100 * artPop ./ sum(popTot , 2)];
area(tVec , hiv_art); %art ./ sum(popVec , 2) , tVec , hiv ./ sum(popVec , 2))
xlabel('Year')
ylabel('Proportion of Population (%)')
title('Relative HIV Prevalence')
legend('Untreated', 'On ART' , 'Location' , 'NorthWest')
figure()
plot(tVec , 100 * artPop ./ (hivPop + artPop))
xlabel('Year')
ylabel('Proportion of HIV Population')
title('Proportion on ART')
legend('Model' , 'Observed')

%% HPV Precancer/CC in ART

hpvInfInds = toInd(allcomb(10 , 6 , 2 , 1 , ...
    1 : periods , 2 , 4 : 10 , 1 : risk));
cin1Inds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 2 , ...
    1 : periods , 2 , 4 : 10 , 1 : risk));
cin2Inds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 3 , ...
    1 : periods , 2 , 4 : 10 , 1 : risk));
cin3Inds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 4 , ...
    1 : periods , 2 , 4 : 10 , 1 : risk));
ccInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 5 : 7 , ...
    1 : periods , 2 , 4 : 10 , 1 : risk));
hpvNegInds = [toInd(allcomb(10 , 6 , 1 , 1 , ...
    1 : periods , 2 , 4 : 10 , 1 : risk));
    toInd(allcomb(10 , 6 , 1 : hpvTypes , 10 , ...
    1 : periods , 2 , 4 : 10 , 1 : risk))];
vaxInds = [toInd(allcomb(10 , 6 , 1 , 9 , ...
    1 : periods , 2 , 4 : 10 , 1 : risk));
    toInd(allcomb(10 , 6 , 2 , 1 , ...
    2 : periods , 2 , 4 : 10 , 1 : risk))];

cin1_art = sum(noV.popVec(: , cin1Inds) , 2);
cin2_art = sum(noV.popVec(: , cin2Inds) , 2);
cin3_art = sum(noV.popVec(: , cin3Inds) , 2);
cc_art = sum(noV.popVec(: , ccInds) , 2);
hpvNeg = sum(noV.popVec(: , hpvNegInds) , 2);
vaxd = sum(noV.popVec(: , vaxInds) , 2);
hpvInfd = sum(noV.popVec(: , hpvInfInds) , 2);
artIndsF = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 2 , 4 : 10 , 1 : risk));

figure()
artPopF = sum(noV.popVec(: , artIndsF) , 2);
art_hpv = [cin1_art ./ artPopF , cin2_art ./ artPopF , ...
    cin3_art ./ artPopF , cc_art ./ artPopF , hpvInfd ./ artPopF , ...
    hpvNeg ./ artPopF , vaxd ./ artPopF];
area(tVec , art_hpv);
xlabel('Year')
ylabel('Female ART Population Size')
title('Relative Prevalence of Pre-cancer/CC Among Women on ART')
legend('CIN1' , 'CIN2' , 'CIN3' , 'CC' , 'HPV Infected' , 'HPV Negative' , ...
    'Vaccinated' , 'Location' , 'NorthWest')

%% Vaccine Prevalence by Gender
gen = {'Male' , 'Female'};
for g = 2
    vaxInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        2 : periods , g , 4 : 10 , 1 : risk));...
        toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , ...
        1 , g , 4 : 10 , 1 : risk))];
    vaxPop = sum(c70_2vPartial.popVec(: , vaxInds) , 2);
    popTot = c70_2vPartial.popVec(: , toInd(allcomb(1 : disease ,...
        1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : 10, 1 : risk)));
    vaxProp = vaxPop ./ sum(popTot , 2) * 100;
    filename = '2v_Vaxd.xlsx';
%     xlswrite(filename, [tVec(1 : stepsPerYear : end)' , vaxProp(1 : stepsPerYear : end)] , sheet)
    figure()
    plot(tVec , vaxProp)
    title([gen{g} , ' Vaccine Prevalence'])
    xlabel('Year'); ylabel('Proportion Vaccinated')
end

%% Vaccine Prevalence by Gender
gen = {'Male' , 'Female'};
for g = 2
    figure()
    for a = 3 : age
        vaxInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            2 : periods , g , a , 1 : risk));...
            toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , ...
            1 : periods , g , a , 1 : risk))];
        vaxPop = sum(c70_2vPartial.popVec(: , vaxInds) , 2);
        popTot = c70_2vPartial.popVec(: , toInd(allcomb(1 : disease ,...
            1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a, 1 : risk)));
        vaxProp = vaxPop ./ sum(popTot , 2) * 100;
        filename = '2v_Vaxd.xlsx';
        %     xlswrite(filename, [tVec(1 : stepsPerYear : end)' , vaxProp(1 : stepsPerYear : end)] , sheet)
        plot(tVec , vaxProp)
        title([gen{g} , ' Vaccine Prevalence'])
        xlabel('Year'); ylabel('Proportion Vaccinated')
        hold on
    end
end
legend('3' ,'4' ,'5', '6', '7' , '8' , '9' ,'10' ,'11' , '12' , '13' ,...
    '14', '15' ,'16')
%% Infected Vaccinated Prevalence by Gender
gen = {'Male' , 'Female'};
for g = 2
    figure()
    for a = 3 : age
        vaxInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            2 : periods , g , a , 1 : risk))];
        vaxPop = sum(c70_2vPartial.popVec(: , vaxInds) , 2);
        popTot = c70_2vPartial.popVec(: , toInd(allcomb(1 : disease ,...
            1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a, 1 : risk)));
        vaxProp = vaxPop ./ sum(popTot , 2) * 100;
        filename = '2v_Vaxd.xlsx';
        %     xlswrite(filename, [tVec(1 : stepsPerYear : end)' , vaxProp(1 : stepsPerYear : end)] , sheet)
        plot(tVec , vaxProp)
        title([gen{g} , ' Infected and Vaccinated Prevalence'])
        xlabel('Year'); ylabel('Proportion Infected and Vaccinated')
        hold on
    end
end
legend('3' ,'4' ,'5', '6', '7' , '8' , '9' ,'10' ,'11' , '12' , '13' ,...
    '14', '15' ,'16')
