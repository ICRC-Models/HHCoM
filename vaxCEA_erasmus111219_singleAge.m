function vaxCEA_erasmus11219_singleAge %(pathModifier)

pathModifier = '021420_singleAge_noBaseScreen_noBaseVax_2018_artLimsPop_origHivInit_Erasmus';

waning = 0;    % turn waning on or off

%% Load parameters
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'],'stepsPerYear','disease','viral','hpvTypes','hpvStates','periods',...
    'gender','age','risk','k','toInd','sumall')

sumall = @(x) sum(x(:));

% Load results
nSims = size(dir([pwd , '\HHCoM_Results\Vaccine' , pathModifier, '\' , '*.mat']) , 1);
curr = load([pwd , '\HHCoM_Results\toNow_021420_singleAge_noBaseScreen_noBaseVax_2018_artLimsPop_origHivInit']); % Population up to current year

% Helper functions
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % sums 1 year worth of values
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) ./ stepsPerYear; % finds average value of a quantity within a given year

% Time
c = fix(clock); % get time
currYear = 2018; %c(1); % get the current year from time

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
    vaxResult{n}.popVec = [curr.popVec(1 : end  , :) ; vaxResult{n}.popVec(2 : end , :)];
    vaxResult{n}.newHpv= [curr.newHpv(1 : end , : , : , : , :) ; vaxResult{n}.newHpv(2 : end , : , : , : ,:)];
    vaxResult{n}.newImmHpv= [curr.newImmHpv(1 : end , : , : , : , :) ; vaxResult{n}.newImmHpv(2 : end , : , : , : , :)];
    vaxResult{n}.newVaxHpv= [curr.newVaxHpv(1 : end , : , : , : , :) ; vaxResult{n}.newVaxHpv(2 : end , : , : , : , :)];
    vaxResult{n}.newCC = [curr.newCC(1 : end , : , : , :) ; vaxResult{n}.newCC(2 : end , : , : ,:)];
    vaxResult{n}.newHiv = [curr.newHiv(1 : end , : , : , :) ; vaxResult{n}.newHiv(2 : end , : , : ,:)];
    vaxResult{n}.artTreatTracker = [curr.artTreatTracker(1 : end , :  , : , : , : , :) ; vaxResult{n}.artTreatTracker(2 : end , : , : , : , : , :)];
    %temp_tVec = linspace(1910 , 2018 , 649);
    %vaxResult{n}.tVec = [temp_tVec(1 : end) , vaxResult{n}.tVec(2 : end)];
    vaxResult{n}.tVec = [curr.tVec(1 : end) , vaxResult{n}.tVec(2 : end)];
end

noVaxInd = nSims;

noV = vaxResult{noVaxInd};
tVec = noV.tVec;

%% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 2)

%% ART coverage
name = {'male-','female-'};
inds = {1 , 2};

for n = 1 : nSims
    for gInd = 1 : 2
        g = inds{gInd};
        hivInds = toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 16 : age , 1 : risk));
        artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
            1 : periods , g , 16 : age , 1 : risk)); 

        hivPop = sum(vaxResult{n}.popVec(: , hivInds) , 2);
        artPop = sum(vaxResult{n}.popVec(: , artInds) , 2);
        artPrev = (100 * (artPop ./ hivPop));

        fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\ART_' , ...
             'Cov- ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
        sname = [name{gInd} , 'ARTcoverage'];
        xlswrite(fname , [tVec(1 : stepsPerYear : end)' , artPrev(1 : stepsPerYear : end)] , sname);
    end
end

%% HIV prevalance by gender and HPV status
name = {'HIV_male-','HIV_female-','All_pop-'};
inds = {1 , 2 , [1:2]};

for n = 1 : nSims
    for gInd = 1 : 3
        g = inds{gInd};
        hivInds = toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 16 : 50 , 1 : risk));
        allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
            1 : periods , g , 16 : 50 , 1 : risk)); 

        hivInds_hpvPos = toInd(allcomb([2:6,10] , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
            1 : periods , g , 16 : 50 , 1 : risk));
        hpvPosInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
            1 : periods , g , 16 : 50 , 1 : risk));

        hivInds_hpvNeg = [toInd(allcomb([2:6,10] , 1 : viral , 1 , [1,8:9] , ...
            1 : periods , g , 16 : 50 , 1 : risk)); ...
            toInd(allcomb([2:6,10] , 1 : viral , 2 : hpvTypes , 10 , ...
            1 : periods , g , 16 : 50 , 1 : risk))];
        hpvNegInds = [toInd(allcomb(1 : disease , 1 : viral , 1 , [1,8:9] , ...
            1 : periods , g , 16 : 50 , 1 : risk));
            toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 10 , ...
            1 : periods , g , 16 : 50 , 1 : risk))];

        hivPop = sum(vaxResult{n}.popVec(: , hivInds) , 2);
        allPop = sum(vaxResult{n}.popVec(: , allInds) , 2);

        hivPop_hpvPos = sum(vaxResult{n}.popVec(: , hivInds_hpvPos) , 2);
        hivPop_hpvNeg = sum(vaxResult{n}.popVec(: , hivInds_hpvNeg) , 2);
        hpvPosPop = sum(vaxResult{n}.popVec(: , hpvPosInds) , 2);
        hpvNegPop = sum(vaxResult{n}.popVec(: , hpvNegInds) , 2);

        hivPrev = (100 * (hivPop ./ allPop));
        hivPrev_hpvPos = (100 * (hivPop_hpvPos ./ hpvPosPop));
        hivPrev_hpvNeg = (100 * (hivPop_hpvNeg ./ hpvNegPop));

        fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\HIV_' , ...
             'Cov- ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
        sname = [name{gInd} , 'HIVprev_hpvStatus'];
        xlswrite(fname , [tVec(1 : stepsPerYear : end)' , hivPrev(1 : stepsPerYear : end) , hivPrev_hpvPos(1 : stepsPerYear : end) , hivPrev_hpvNeg(1 : stepsPerYear : end)] , sname);
    end
end

%% HIV prevalance by gender, HPV status, and age
name = {'HIV_male-','HIV_female-','All_pop-'};
inds = {1,2,[1:2]};
inds2 = {11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,[51:age]};
aVec_all = [];
aVec_hpvP = [];
aVec_hpvN = [];

for n = 1 : nSims
    for gInd = 1 : 3
        g = inds{gInd};
        aVec_all = [];
        aVec_hpvP = [];
        aVec_hpvN = [];
        for aInd = 1 : 8
            a = inds2{aInd};
            hivInds = toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
                1 : periods , g , a , 1 : risk));
            allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
                1 : periods , g , a , 1 : risk)); 

            hivInds_hpvPos = toInd(allcomb([2:6,10] , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
                1 : periods , g , a , 1 : risk));
            hpvPosInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
                1 : periods , g , a , 1 : risk));

            hivInds_hpvNeg = [toInd(allcomb([2:6,10] , 1 : viral , 1 , [1,8:9] , ...
                1 : periods , g , a , 1 : risk)); ...
                toInd(allcomb([2:6,10] , 1 : viral , 2 , 10 , ...
                1 : periods , g , a , 1 : risk))];
            hpvNegInds = [toInd(allcomb(1 : disease , 1 : viral , 1 , [1,8:9] , ...
                1 : periods , g , a , 1 : risk));
                toInd(allcomb(1 : disease , 1 : viral , 2 , 10 , ...
                1 : periods , g , a , 1 : risk))];

            hivPop = sum(vaxResult{n}.popVec(: , hivInds) , 2);
            allPop = sum(vaxResult{n}.popVec(: , allInds) , 2);

            hivPop_hpvPos = sum(vaxResult{n}.popVec(: , hivInds_hpvPos) , 2);
            hivPop_hpvNeg = sum(vaxResult{n}.popVec(: , hivInds_hpvNeg) , 2);
            hpvPosPop = sum(vaxResult{n}.popVec(: , hpvPosInds) , 2);
            hpvNegPop = sum(vaxResult{n}.popVec(: , hpvNegInds) , 2);

            hivPrev = (100 * (hivPop ./ allPop));
            hivPrev_hpvPos = (100 * (hivPop_hpvPos ./ hpvPosPop));
            hivPrev_hpvNeg = (100 * (hivPop_hpvNeg ./ hpvNegPop));
            
            aVec_all = [aVec_all, hivPrev];
            aVec_hpvP = [aVec_hpvP, hivPrev_hpvPos];
            aVec_hpvN = [aVec_hpvN, hivPrev_hpvNeg];
        end
        
        fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\HIV_' , ...
             'Cov- ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
        sname = [name{gInd} , 'byAge_all'];
        xlswrite(fname , [tVec(1 : stepsPerYear : end)' , aVec_all(1 : stepsPerYear : end, :)] , sname);
        sname = [name{gInd} , 'byAge_HPV+'];
        xlswrite(fname , [tVec(1 : stepsPerYear : end)' , aVec_hpvP(1 : stepsPerYear : end, :)] , sname);
        sname = [name{gInd} , 'byAge_HPV-'];
        xlswrite(fname , [tVec(1 : stepsPerYear : end)' , aVec_hpvN(1 : stepsPerYear : end, :)] , sname);
    
    end
end

%% HIV incidence
name = {'male-','female-','allPop-'};
iName = {'allH'}; % , 'hpvPos' , 'hpvNeg'};
indsG = {1 , 2 , [1:2]};
%inds = {':' , [2:6,10] , [1,7:9]};
fac = 1000;

for n = 1 : nSims
   for gInd = 1 : 3
        g = indsG{gInd};
        for i = 1 % : length(inds)
            % General
            allF = toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
                1 : periods , g , 16:50 , 1 : risk));
            % HPV-positive women
            %hpvPos = toInd(allcomb([1,7:9] , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
            %    1 : periods , g , 16:50 , 1 : risk));
            % HPV-negative women
            %hpvNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 , [1,8:9] , 1 : periods , ...
            %    g , 16:50, 1 : risk)); ...
            %    toInd(allcomb([1,7:9] , 1 : viral , 2 : hpvTypes , 10 , 1 : periods , ...
            %    g , 16:50 , 1 : risk))];
            genArray = {allF}; % , hpvPos , hpvNeg};

            hivIncRef = ...
                (annlz(sum(sum(sum(vaxResult{n}.newHiv(: , g , 16:50 , :),2),3),4)) ./ ...
                (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
            vaxResult{n}.hivIncRef = hivIncRef;
            vaxResult{n}.hivInc = vaxResult{n}.hivIncRef;
            
            fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\HIVincidence_' , ...
                'Cov- ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
            sname = [name{gInd} , 'HIVinc_' , iName{i}];
            xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.hivInc'] , sname);
    
        end 
   end
end

%% HPV prevalance by gender and HIV status
name = {'male-','female-','allPop-'};
inds = {1 , 2 , [1:2]};

for n = 1 : nSims
    for gInd = 1 : 3
        g = inds{gInd};
        hpvInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
            1 : periods , g , 16 : 50 , 1 : risk));
        allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
            1 : periods , g , 16 : 50 , 1 : risk)); 

        hpvInds_hivPos = toInd(allcomb([2:6,10] , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
            1 : periods , g , 16 : 50 , 1 : risk));
        hivPosInds = toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 16 : 50 , 1 : risk));

        hpvInds_hivNeg = toInd(allcomb([1,7:9] , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
            1 : periods , g , 16 : 50 , 1 : risk));
        hivNegInds = toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            1 : periods , g , 16 : 50 , 1 : risk));

        hpvPop = sum(vaxResult{n}.popVec(: , hpvInds) , 2);
        allPop = sum(vaxResult{n}.popVec(: , allInds) , 2);

        hpvPop_hivPos = sum(vaxResult{n}.popVec(: , hpvInds_hivPos) , 2);
        hpvPop_hivNeg = sum(vaxResult{n}.popVec(: , hpvInds_hivNeg) , 2);
        hivPosPop = sum(vaxResult{n}.popVec(: , hivPosInds) , 2);
        hivNegPop = sum(vaxResult{n}.popVec(: , hivNegInds) , 2);

        hpvPrev = (100 * (hpvPop ./ allPop));
        hpvPrev_hivPos = (100 * (hpvPop_hivPos ./ hivPosPop));
        hpvPrev_hivNeg = (100 * (hpvPop_hivNeg ./ hivNegPop));

        fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\HPV_' , ...
             'Cov- ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
        sname = [name{gInd} , 'HPVprev_hivStatus'];
        xlswrite(fname , [tVec(1 : stepsPerYear : end)' , hpvPrev(1 : stepsPerYear : end) , hpvPrev_hivPos(1 : stepsPerYear : end) , hpvPrev_hivNeg(1 : stepsPerYear : end)] , sname);
    end
end

%% HPV prevalance by gender, HIV status, and age
name = {'male-','female-','allPop-'};
inds = {1 , 2 , [1:2]};
inds2 = {11:15,16:20,21:25,26:30,31:35,36:40,41:45,46:50,[51:age]};
aVec_all = [];
aVec_hivP = [];
aVec_hivN = [];

for n = 1 : nSims
    for gInd = 1 : 3
        g = inds{gInd};
        aVec_all = [];
        aVec_hivP = [];
        aVec_hivN = [];
        for aInd = 1 : 8
            a = inds2{aInd};
            
            hpvInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
                1 : periods , g , a , 1 : risk));
            allInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
                1 : periods , g , a , 1 : risk)); 

            hpvInds_hivPos = toInd(allcomb([2:6,10] , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
                1 : periods , g , a , 1 : risk));
            hivPosInds = toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
                1 : periods , g , a , 1 : risk));

            hpvInds_hivNeg = toInd(allcomb([1,7:9] , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
                1 : periods , g , a , 1 : risk));
            hivNegInds = toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
                1 : periods , g , a , 1 : risk));
            
            hpvPop = sum(vaxResult{n}.popVec(: , hpvInds) , 2);
            allPop = sum(vaxResult{n}.popVec(: , allInds) , 2);
            
            hpvPop_hivPos = sum(vaxResult{n}.popVec(: , hpvInds_hivPos) , 2);
            hpvPop_hivNeg = sum(vaxResult{n}.popVec(: , hpvInds_hivNeg) , 2);
            hivPosPop = sum(vaxResult{n}.popVec(: , hivPosInds) , 2);
            hivNegPop = sum(vaxResult{n}.popVec(: , hivNegInds) , 2);

            hpvPrev = (100 * (hpvPop ./ allPop));
            hpvPrev_hivPos = (100 * (hpvPop_hivPos ./ hivPosPop));
            hpvPrev_hivNeg = (100 * (hpvPop_hivNeg ./ hivNegPop));
          
            aVec_all = [aVec_all, hpvPrev];
            aVec_hivP = [aVec_hivP, hpvPrev_hivPos];
            aVec_hivN = [aVec_hivN, hpvPrev_hivNeg];
        end
        
        fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\HPV_' , ...
             'Cov- ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
        sname = [name{gInd} , 'byAge_all'];
        xlswrite(fname , [tVec(1 : stepsPerYear : end)' , aVec_all(1 : stepsPerYear : end, :)] , sname);
        sname = [name{gInd} , 'byAge_HIV+'];
        xlswrite(fname , [tVec(1 : stepsPerYear : end)' , aVec_hivP(1 : stepsPerYear : end, :)] , sname);
        sname = [name{gInd} , 'byAge_HIV-'];
        xlswrite(fname , [tVec(1 : stepsPerYear : end)' , aVec_hivN(1 : stepsPerYear : end, :)] , sname);
    
    end
end

%% HPV incidence
name = {'HPV_male-','HPV_female-','All_pop-'};
iName = {'allH' , 'hivPos' , 'hivNeg'};
indsG = {1 , 2 , [1:2]};
inds = {':' , [2:6,10] , [1,7:9]};
fac = 1000;

for n = 1 : nSims
   for gInd = 1 : 3
        g = indsG{gInd};
        for i = 1 : length(inds)
            % General
            allH = [toInd(allcomb(1 : disease , 1 : viral , 1 , [1,9] , ...
                1 : periods , g , 16:age , 1 : risk)); ...
                toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 10 , ...
                1 : periods , g , 16:age , 1 : risk))];
            % HIV-positive
            hivPos = [toInd(allcomb([2:6,10] , 1 : viral , 1 , [1,9] , ...
                1 : periods , g , 16:age , 1 : risk)); ...
                toInd(allcomb([2:6,10] , 1 : viral , 1 : hpvTypes , 10 , ...
                1 : periods , g , 16:age , 1 : risk))];
            % All HIV-negative
            hivNeg = [toInd(allcomb([1,7:9] , 1 : viral , 1 , [1,9] , 1 : periods , ...
                g , 16:age, 1 : risk)); ...
                toInd(allcomb([1,7:9] , 1 : viral , 1 : hpvTypes , 10 , 1 : periods , ...
                g , 16:age , 1 : risk))];
            genArray = {allH , hivPos , hivNeg};

            hpvIncRef = ...
                ((annlz(sum(sum(sum(sum(vaxResult{n}.newHpv(: , g , inds{i} , 16:age , :),2),3),4),5))+annlz(sum(sum(sum(sum(vaxResult{n}.newImmHpv(: , g , inds{i} , 16:age , :),2),3),4),5))+annlz(sum(sum(sum(sum(vaxResult{n}.newVaxHpv(: , g , inds{i} , 16:age , :),2),3),4),5))) ./ ...
                (annlz(sum(vaxResult{n}.popVec(: , genArray{i}) , 2) ./ stepsPerYear)) * fac);
            vaxResult{n}.hpvIncRef = hpvIncRef;
            vaxResult{n}.hpvInc = vaxResult{n}.hpvIncRef;
            
            fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\HPVincidence_' , ...
                'Cov- ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
            sname = [name{gInd} , 'HPVinc_' , iName{i}];
            xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.hpvInc'] , sname);
    
        end 
   end
end

%% HPV incidence- vax, non-vax
name = 'HPV_female-';
indsG = 2;
fac = 1000;

for n = 1 : nSims
   for gInd = 1
        g = indsG;
        for i = 1 : length(inds)
            % General
            vax = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , ...
                2 : 5 , g , 16:age , 1 : risk)); ...
                toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 10 , ...
                2 : 5 , g , 16:age , 1 : risk)); ...
                toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , ...
                [1,6] , g , 16:age , 1 : risk))];
            nVax = [toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , ...
                [1,6] , g , 16:age , 1 : risk)); ...
                toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 10 , ...
                [1,6] , g , 16:age , 1 : risk))];

            hpvIncRef_vax = ...
                (annlz(sum(sum(sum(sum(vaxResult{n}.newVaxHpv(: , g , : , 16:age , :),2),3),4),5)) ./ ...
                (annlz(sum(vaxResult{n}.popVec(: , vax) , 2) ./ stepsPerYear)) * fac);
            vaxResult{n}.hpvIncRef_vax = hpvIncRef_vax;
            vaxResult{n}.hpvInc_vax = vaxResult{n}.hpvIncRef_vax;
            
            hpvIncRef_nVax = ...
                ((annlz(sum(sum(sum(sum(vaxResult{n}.newHpv(: , g , : , 16:age , :),2),3),4),5))+annlz(sum(sum(sum(sum(vaxResult{n}.newImmHpv(: , g , : , 16:age , :),2),3),4),5))) ./ ...
                (annlz(sum(vaxResult{n}.popVec(: , nVax) , 2) ./ stepsPerYear)) * fac);
            vaxResult{n}.hpvIncRef_nVax = hpvIncRef_nVax;
            vaxResult{n}.hpvInc_nVax = vaxResult{n}.hpvIncRef_nVax;
            
            fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\HPVincidenceVax_' , ...
                'Cov- ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
            sname = [name , 'HPVinc_vax-nVax'];
            xlswrite(fname , [tVec(1 : stepsPerYear : end)' , vaxResult{n}.hpvInc_vax' , vaxResult{n}.hpvInc_nVax'] , sname);
    
        end 
   end
end

%% HPV prevalence- vax, non-vax
name = {'female-'};
inds = {2};

for n = 1 : nSims
    for gInd = 1
        g = inds{gInd};
        hpvVaxInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
            2 : 5 , g , 16 : 50 , 1 : risk));
        allVaxInds = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            2 : 5 , g , 16 : 50 , 1 : risk)); ...
            toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , ...
            [1,6] , g , 16 : 50 , 1 : risk))];
            
        hpvNvaxInds = toInd(allcomb(1 : disease , 1 : viral , 2 : hpvTypes , 1 : 7 , ...
            [1,6] , g , 16 : 50 , 1 : risk));
        allNvaxInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
            [1,6] , g , 16 : 50 , 1 : risk));

        hpvPop = sum(vaxResult{n}.popVec(: , hpvVaxInds) , 2);
        allPop = sum(vaxResult{n}.popVec(: , allVaxInds) , 2);

        hpvNvaxPop = sum(vaxResult{n}.popVec(: , hpvNvaxInds) , 2);
        allNvaxPop = sum(vaxResult{n}.popVec(: , allNvaxInds) , 2);

        hpvPrev = (100 * (hpvPop ./ allPop));
        hpvNvaxPrev = (100 * (hpvNvaxPop ./ allNvaxPop));

        fname = [pwd , '\HHCoM_Results\Vaccine' , pathModifier , '\HPVvax_' , ...
             'Cov- ' , num2str(round(vaxResult{n}.vaxRate * 100)) , '.xlsx'];
        sname = [name{gInd} , 'HPVprev_vax-noVax'];
        xlswrite(fname , [tVec(1 : stepsPerYear : end)' , hpvPrev(1 : stepsPerYear : end) , hpvNvaxPrev(1 : stepsPerYear : end)] , sname);
    end
end

