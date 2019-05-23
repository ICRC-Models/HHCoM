function vaxCompare()
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
%% Vaccine
%% HIV incidence
yrs = noV.tVec(end) - noV.tVec(1);
hivInc = zeros(yrs , gender);
for g = 1 : 2
    hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : 10 , 1 : risk)); ...
        toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk))];
    hivSus = annlz(sum(c70_2vPartial.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
    hivInc(: , g) = (annlz(sum(sum(c70_2vPartial.newHiv(: , g , 4 : 10 , :) ...
        , 3) , 4)) ./ hivSus * 100)';
end
filename = 'HIV Incidence Per 100 Vax.xlsx';
xlswrite(filename, [(tVec(1 : stepsPerYear: end))' , hivInc])
%% HIV Prevalence by Gender
gen = {'Male' , 'Female'};
for g = 1 : gender
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hivPop = sum(c70_2vPartial.popVec(: , hivInds) , 2);
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    popTot = c70_2vPartial.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : 10, 1 : risk)));
    art = sum(c70_2vPartial.popVec(: , artInds) , 2);
    prev = (hivPop + art) ./ sum(popTot , 2) * 100;
    sheet = gen{g};
    filename = 'Gender Specific HIV Prevalence Vax.xlsx';
    xlswrite(filename, [tVec(1 : stepsPerYear : end)' , prev(1 : stepsPerYear : end)] , sheet)
end

%% ART coverage by gender
artCover = zeros(yrs , gender);
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; 

for g = 1 : gender
    art = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 4 : 10 , 1 : risk));
    allHiv = [toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 4 : 10 , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 4 : 10 , 1 : risk))];
    artCover(: , g) = annAvg(sum(c70_2vPartial.popVec(: , art) , 2) ...
        ./ sum(c70_2vPartial.popVec(: , allHiv) , 2));
end
        
filename = 'ART Coverage Vax.xlsx';
xlswrite(filename, [tVec(1 : stepsPerYear : end)' , artCover])



%% No vaccine
%% HIV incidence
yrs = tVec(end) - tVec(1);
hivInc = zeros(yrs , gender);
for g = 1 : 2
    hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : 10 , 1 : risk)); ...
        toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk))];
    hivSus = annlz(sum(noV.popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
    hivInc(: , g) = (annlz(sum(sum(noV.newHiv(: , g , 4 : 10 , :) ...
        , 3) , 4)) ./ hivSus * 100)';
end
filename = 'HIV Incidence Per 100 No Vax.xlsx';
xlswrite(filename, [(tVec(1 : stepsPerYear: end))' , hivInc])
%% HIV Prevalence by Gender
gen = {'Male' , 'Female'};
for g = 1 : gender
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hivPop = sum(noV.popVec(: , hivInds) , 2);
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    popTot = noV.popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : 10, 1 : risk)));
    art = sum(noV.popVec(: , artInds) , 2);
    prev = (hivPop + art) ./ sum(popTot , 2) * 100;
    sheet = gen{g};
    filename = 'Gender Specific HIV Prevalence No Vax.xlsx';
    xlswrite(filename, [tVec(1 : stepsPerYear : end)' , prev(1 : stepsPerYear : end)] , sheet)
end

%% ART coverage by gender
artCover = zeros(yrs , gender);
for g = 1 : gender
    art = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 4 : 10 , 1 : risk));
    allHiv = [toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 4 : 10 , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 4 : 10 , 1 : risk))];
    artCover(: , g) = annAvg(sum(noV.popVec(: , art) , 2) ...
        ./ sum(noV.popVec(: , allHiv) , 2));
end
        
filename = 'ART Coverage No Vax.xlsx';
xlswrite(filename, [tVec(1 : stepsPerYear : end)' , artCover])


%% HIV Prevalence by HPV status and Gender VAX
gen = {'Male' , 'Female'};
hpvStat = {'HPV Negative' , 'HPV Positive'};
hArr = {1 , 2 : 4};
for h = 1 : length(hArr)
    for g = 1 : gender
        hivInds = toInd(allcomb(2 : 6 , 1 : viral , hArr{h} , 1 : hpvStates, ...
            1 : periods , g , 4 : 10, 1 : risk));
        hivPop = sum(c70_2vPartial.popVec(: , hivInds) , 2);
        artInds = toInd(allcomb(10 , 6 , hArr{h} , 1 : hpvStates, ...
            1 : periods , g , 4 : 10, 1 : risk));
        popTot = c70_2vPartial.popVec(: , toInd(allcomb(1 : disease , 1 : viral , hArr{h} , 1 : hpvStates , 1 : periods , ...
            g , 4 : 10, 1 : risk)));
        art = sum(c70_2vPartial.popVec(: , artInds) , 2);
        prev = (hivPop + art) ./ sum(popTot , 2) * 100;
        tit = hpvStat{1};
        if h == 2
            tit = hpvStat{2};
        end
        sheet = [gen{g} , tit];
        filename = 'Vax_Gender Specific HIV Prevalence by HPV Status.xlsx';
        xlswrite(filename, [tVec(1 : stepsPerYear : end)' , prev(1 : stepsPerYear : end)] , sheet)
    end
end
%% hiv Prevalence by HPV status and Gender NON VAX
gen = {'Male' , 'Female'};
hpvStat = {'HPV Negative' , 'HPV Positive'};
hArr = {1 , 2 : 4};
for h = 1 : length(hArr)
    for g = 1 : gender
        hivInds = toInd(allcomb(2 : 6 , 1 : viral , hArr{h} , 1 : hpvStates, ...
            1 : periods , g , 4 : 10, 1 : risk));
        hivPop = sum(noV.popVec(: , hivInds) , 2);
        artInds = toInd(allcomb(10 , 6 , hArr{h} , 1 : hpvStates, ...
            1 : periods , g , 4 : 10, 1 : risk));
        popTot = noV.popVec(: , toInd(allcomb(1 : disease , 1 : viral , hArr{h} , 1 : hpvStates , 1 : periods , ...
            g , 4 : 10, 1 : risk)));
        art = sum(noV.popVec(: , artInds) , 2);
        prev = (hivPop + art) ./ sum(popTot , 2) * 100;
        tit = hpvStat{1};
        if h == 2
            tit = hpvStat{2};
        end
        sheet = [gen{g} , tit];
        filename = 'NoVax_Gender Specific HIV Prevalence by HPV Status.xlsx';
        xlswrite(filename, [tVec(1 : stepsPerYear : end)' , prev(1 : stepsPerYear : end)] , sheet)
    end
end
%% HPV Incidence by HIV status


