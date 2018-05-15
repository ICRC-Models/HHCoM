function[] = outputResults()
%%
paramDir = [pwd , '\Params\'];
savedir = [pwd , '\Results\'];
load([paramDir,'actual'])
load([paramDir,'calibData'])
% load('C:\Users\nicktzr\Google Drive\ICRC\CISNET\Results\to2017')
load('H:\HHCoM_Results\toNow')
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; 
yrs = endYear - startYear;
%% HIV Prevalence by Gender
gen = {'Male' , 'Female'};
for g = 1 : gender
    hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    hivPop = sum(popVec(: , hivInds) , 2);
    artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : 10 , 1 : risk));
    popTot = popVec(: , toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : 10, 1 : risk)));
    art = sum(popVec(: , artInds) , 2);
    prev = (hivPop + art) ./ sum(popTot , 2) * 100;
    sheet = gen{g};
    filename = 'Gender Specific HIV Prevalence.xlsx';
    xlswrite(filename, [tVec(1 : stepsPerYear : end)' , prev(1 : stepsPerYear : end)] , sheet)
end
%% HIV Prevalence by age
% hivAge = zeros(age , length(tVec));
ageGroup = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};

gen = {'Male' , 'Female'};
for g = 1 : gender
    for a = 4 : 10
        hivAgeInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a , 1 : risk)); toInd(allcomb(10, 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a , 1 : risk))];
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a , 1 : risk));
        hivAge(a , :) = sum(popVec(: , hivAgeInds) , 2);
        hivAgeRel = bsxfun(@rdivide , hivAge(a , :)' , sum(popVec(: , ageInds) , 2)) * 100;
        sheet = [gen{g} , ' ' , ageGroup{a}];
        filename = 'Age Specific HIV Prevalence.xlsx';
        xlswrite(filename, [tVec(1 : stepsPerYear : end)' , hivAgeRel(1 : stepsPerYear : end)] , sheet)
    end
end

%% CD4 distribution
for d = 2 : 6
    cd4Group = toInd(allcomb(d , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 1 : gender , 1 : age , 1 : risk));
    cd4Pop(: , d - 1) = sum(popVec(: , cd4Group) , 2);
end
hivInf = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
    1 : periods , 1 : gender , 1 : age , 1 : risk));
totalHiv = sum(popVec(: , hivInf) , 2);
cd4Dist = bsxfun(@rdivide , cd4Pop , totalHiv);
cd4DistAll_Ann = zeros(size(cd4Dist , 1) / stepsPerYear , size(cd4Dist , 2));
for i = 1 : size(cd4Dist , 2)
    cd4DistAll_Ann(: , i) = annAvg(cd4Dist(: , i));
end
cd4Dist_Ann(: , 1) = sum(cd4DistAll_Ann(: , 1 : 3) , 2); % 350+
cd4Dist_Ann(: , 2) = sum(cd4DistAll_Ann(: , 4) , 2); % 200-350
cd4Dist_Ann(: , 3) = sum(cd4DistAll_Ann(: , 5) , 2); % < 200
filename = 'CD4 Distribution.xlsx';
xlswrite(filename, [tVec(1 : stepsPerYear : end)' , cd4Dist_Ann])

%% ART coverage by gender
artCover = zeros(yrs , gender);
for g = 1 : gender
    art = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 4 : age , 1 : risk));
    allHiv = [toInd(allcomb(2 : 6 , 1 : 5 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 4 : age , 1 : risk)); ...
            toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 4 : age , 1 : risk))];
    artCover(: , g) = annAvg(sum(popVec(: , art) , 2) ./ sum(popVec(: , allHiv) , 2));
end
        
filename = 'ART Coverage.xlsx';
xlswrite(filename, [(startYear : endYear - 1)' , artCover])
    
%% HIV Incidence
hivInc = zeros(yrs , gender);
for g = 1 : 2
    hivSusInds = [toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , g , 4 : age , 1 : risk)); ...
        toInd(allcomb(7 : 9 , 1 , 1 : hpvTypes , 1 : hpvStates, ...
        1 : periods , g , 4 : age , 1 : risk))];
    hivSus = annlz(sum(popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
    hivInc(: , g) = annlz(sum(sum(newHiv(: , g , 4 : age , :) ...
        , 3) , 4)) ./ hivSus * 100;
end
filename = 'HIV Incidence Per 100.xlsx';
xlswrite(filename, [(startYear : endYear - 1)' , hivInc])

%% HPV Incidence
hpvInc = zeros(yrs , gender);
for g = 1 : 2
    hpvSusInds = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , ...
        1 : periods , g , 4 : age , 1 : risk))];
    hpvSus = annlz(sum(popVec(: , hpvSusInds) , 2)) ./ stepsPerYear;
    hpvInc(: , g) = annlz(sum(sum(sum(newHpv(: , g , : , 4 : age , :) ...
        , 3) , 4), 5)) ./ hpvSus * 100;
end
filename = 'HPV Incidence Per 100.xlsx';
xlswrite(filename, [(startYear : endYear - 1)' , hpvInc])

%% HPV Prevalence by Gender
gen = {'Male' , 'Female'};
for g = 1 : gender
    hpvInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 9, ...
        1 : periods , g , 4 : age , 1 : risk));
    popTot = popVec(: , toInd(allcomb(1 : disease , 1 : viral , ...
        1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
        g , 4 : age, 1 : risk)));
    hpvPrev = sum(popVec(: , hpvInds) , 2) ./ sum(popTot , 2) * 100;
    sheet = gen{g};
    filename = 'Gender Specific HPV Prevalence.xlsx';
    xlswrite(filename, [tVec(1 : stepsPerYear : end)' , hpvPrev(1 : stepsPerYear : end)] , sheet)
end

%% HPV By Age
gen = {'Male' , 'Female'};
hpvAgeRel = zeros(age , length(tVec));
for g = 1 : gender
    for a = 4 : 10
        hpvAgeInds = [toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a , 1 : risk))];
        ageInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , a , 1 : risk));
        hpvAge(a , :) = sum(popVec(: , hivAgeInds) , 2);
        hpvAgeRel(a , :) = bsxfun(@rdivide , hpvAge(a , :)' , sum(popVec(: , ageInds) , 2)) * 100;
        sheet = [gen{g}];
        filename = 'Age Specific HPV Prevalence.xlsx';
        xlswrite(filename, [tVec(1 : stepsPerYear : end)' , hpvAgeRel(: , 1 : stepsPerYear : end)'] , sheet)
    end
end

%% HIV Prevalence by HPV status and Gender
gen = {'Male' , 'Female'};
hpvStat = {'HPV Negative' , 'HPV Positive'};
hArr = {1 , 2 : 4};
for h = 1 : length(hArr)
    for g = 1 : gender
        hivInds = toInd(allcomb(2 : 6 , 1 : viral , hArr{h} , 1 : hpvStates, ...
            1 : periods , g , 4 : age , 1 : risk));
        hivPop = sum(popVec(: , hivInds) , 2);
        artInds = toInd(allcomb(10 , 6 , hArr{h} , 1 : hpvStates, ...
            1 : periods , g , 4 : age , 1 : risk));
        popTot = popVec(: , toInd(allcomb(1 : disease , 1 : viral , hArr{h} , 1 : hpvStates , 1 : periods , ...
            g , 4 : age, 1 : risk)));
        art = sum(popVec(: , artInds) , 2);
        prev = (hivPop + art) ./ sum(popTot , 2) * 100;
        tit = hpvStat{1};
        if h == 2
            tit = hpvStat{2};
        end
        sheet = [gen{g} , tit];
        filename = 'Gender Specific HIV Prevalence by HPV Status.xlsx';
        xlswrite(filename, [tVec(1 : stepsPerYear : end)' , prev(1 : stepsPerYear : end)] , sheet)
    end
end

%% HIV Prevalence by age, gender, HPV status
% hivAge = zeros(age , length(tVec));
ageGroup = {'0 - 4' , '5 - 9' , '10 - 14' , '15 - 19' , '20 -24' , '25 - 29' ,...
    '30 -34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};

gen = {'Male' , 'Female'};
hivAgeRel = zeros(age , length(tVec));
for g = 1 : gender
    for a = 4 : 10
        hivAgeInds = [toInd(allcomb(2 : 6 , 1 : viral , 1 , 1 , 1 : periods , ...
            g , a , 1 : risk)); toInd(allcomb(2 : 6 , 1 : viral , 2 : 4 , 10 , 1 : periods , ...
            g , a , 1 : risk)); toInd(allcomb(10, 6 , 1 , 1 , 1 : periods , ...
            g , a , 1 : risk)); toInd(allcomb(10, 6 , 2 : 4 , 10 , 1 : periods , ...
            g , a , 1 : risk));];
        ageInds = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 : periods , ...
            g , a , 1 : risk));toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 10 , 1 : periods , ...
            g , a , 1 : risk))];
        hivAge(a , :) = sum(popVec(: , hivAgeInds) , 2);
        hivAgeRel(a , :) = bsxfun(@rdivide , hivAge(a , :)' , sum(popVec(: , ageInds) , 2)) * 100;
        filename = 'Age Specific HIV Prevalence in HPV Negative.xlsx';
        sheet = gen{g};
        xlswrite(filename, [tVec(1 : stepsPerYear : end)' , hivAgeRel(: , 1 : stepsPerYear : end)'] , sheet)
    end
end

gen = {'Male' , 'Female'};
hivAgeRel = zeros(age , length(tVec));
for g = 1 : gender
    for a = 4 : 10
        hivAgeInds = [toInd(allcomb(2 : 6 , 1 : viral , 2 : 4 , 1 : 7 , 1 : periods , ...
            g , a , 1 : risk)); toInd(allcomb(10, 6 , 2 : 4 , 10 , 1 : periods , ...
            g , a , 1 : risk));];
        ageInds = [toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : 7 , 1 : periods , ...
            g , a , 1 : risk))];
        hivAge(a , :) = sum(popVec(: , hivAgeInds) , 2);
        hivAgeRel(a , :) = bsxfun(@rdivide , hivAge(a , :)' , sum(popVec(: , ageInds) , 2)) * 100;
        filename = 'Age Specific HIV Prevalence in HPV Positive.xlsx';
        sheet = gen{g};
        xlswrite(filename, [tVec(1 : stepsPerYear : end)' , hivAgeRel(: , 1 : stepsPerYear : end)'] , sheet)
    end
end


%% HIV Incidence by Gender and HPV Status
% hivInc = zeros(yrs , gender);
% for h = 1 : length(hArr)
%     for g = 1 : 2
%         hivSusInds = [toInd(allcomb(1 , 1 , hArr{h} , 1 : hpvStates , ...
%             1 : periods , g , 4 : 10 , 1 : risk)); ...
%             toInd(allcomb(7 : 9 , 1 , hArr{h} , 1 : hpvStates, ...
%             1 : periods , g , 4 : 10 , 1 : risk))];
%         hivSus = annlz(sum(popVec(: , hivSusInds) , 2)) ./ stepsPerYear;
%         hivInc(: , g) = annlz(sum(sum(newHiv(: , g , 4 : 10 , :) ...
%             , 3) , 4)) ./ hivSus * 100;
%         tit = hpvStat{1};
%         if h == 2
%             tit = hpvStat{2};
%         end
%     end
%     filename = ['HIV Incidence Per 100 in ' , tit , '.xlsx'];
%     xlswrite(filename, [(startYear : endYear - 1)' , hivInc])
% end

%% HPV Prevalence by Gender and HIV status
gen = {'Male' , 'Female'};
hivStat = {'HIV-Negative' , 'HIV Positive'};
dArr = {[1 , 7 : 9] , [2 : 6 , 10]};
for d = 1 : length(dArr)
    for g = 1 : gender
        hpvInds = toInd(allcomb(dArr{d} , 1 : viral , 2 : hpvTypes , 1 : 7, ...
            1 : periods , g , 4 : age , 1 : risk));
        popTot = popVec(: , toInd(allcomb(dArr{d} , 1 : viral , ...
            1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
            g , 4 : age, 1 : risk)));
        hpvPrev = sum(popVec(: , hpvInds) , 2) ./ sum(popTot , 2) * 100;
        tit = hivStat{1};
        if d == 2
            tit = hivStat{2};
        end
        sheet = [gen{g} , tit];
        filename = 'Gender Specific HPV Prevalence by HIV.xlsx';
        xlswrite(filename, [tVec(1 : stepsPerYear : end)' , hpvPrev(1 : stepsPerYear : end)] , sheet)
    end
end

%% HPV By Age, Gender, and HIV status
gen = {'Male' , 'Female'};
hivStat = {'HIV-Negative' , 'HIV Positive'};
dArr = {[1 , 7 : 9] , [2 : 6 , 10]};
hpvAgeRel = zeros(age , length(tVec));
for d = 1 : length(dArr)
    for g = 1 : gender
        for a = 4 : 10
            hpvAgeInds = [toInd(allcomb(dArr{d} , 1 : viral , 2 : hpvTypes , 1 : 7 , 1 : periods , ...
                g , a , 1 : risk))];
            ageInds = toInd(allcomb(dArr{d} , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
                g , a , 1 : risk));
            hpvAge(a , :) = sum(popVec(: , hpvAgeInds) , 2);
            hpvAgeRel(a  , :) = bsxfun(@rdivide , hpvAge(a , :)' , sum(popVec(: , ageInds) , 2)) * 100;
            tit = hivStat{1};
            if d == 2
                tit = hivStat{2};
            end    
        end
        sheet = [gen{g} , ' ' , tit];
        filename = 'Age Specific HPV Prevalence by HIV.xlsx';
        xlswrite(filename, [tVec(1 : stepsPerYear : end)' , hpvAgeRel(: , 1 : stepsPerYear : end)'] , sheet)
    end
end

%% HPV incidence by Gender and HIV Status
hpvInc = zeros(yrs , gender);
hivStat = {'HIV-Negative' , 'HIV Positive'};
dArr = {[1 , 7 : 9] , [2 : 6 , 10]};
for d = 1 : length(dArr)
    for g = 1 : 2
        hpvSusInds = [toInd(allcomb(dArr{d} , 1 : viral , 1 , 1 , ...
            1 : periods , g , 4 : age , 1 : risk))];
        hpvSus = annlz(sum(popVec(: , hpvSusInds) , 2)) ./ stepsPerYear;
        hpvInc(: , g) = annlz(sum(sum(sum(newHpv(: , g , dArr{d} , 4 : age , :) ...
            , 3) , 4), 5)) ./ hpvSus * 100;
        tit = hivStat{1};
        if d == 2
            tit = hivStat{2};
        end
    end
    filename = ['HPV Incidence Per 100 in ' , tit , '.xlsx'];
    xlswrite(filename, [(startYear : endYear - 1)' , hpvInc])
end
