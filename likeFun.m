% Note: likeFun only set up to use 5-year age groups

function negSumLogL = likeFun(popVec , newCC , cinPos2007_dObs , cin1_2010_dObs ,...
    cin2_2010_dObs, hpv_hiv_dObs , hpv_hivNeg_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
    hivPrevAll_dObs, hpv_all_dObs , hpv_hiv2009_dObs ,  cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , ccInc2012_dObs ,popAgeDist_dObs , totPopSize_dObs , ...
    toInd , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , ...
    age , gender , risk , startYear , stepsPerYear , annlz)

mObs = [];
dMean = [];
dVar = [];

%% CIN2/3 prevalence among HIV+ women aged 20-50 in 2007
cinPos2007 = zeros(3 , 1);
yr = 2007;
ageVec = {[5:6],[7:8],[9:10]};
for aV = 1 : length(ageVec)
    a = ageVec{aV};
    % HIV-positive (on and not on ART)
    cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cinPos2007(aV) = (sum(popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
        ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
end

mObs = [mObs ; cinPos2007];
dMean = [dMean ; cinPos2007_dObs(: , 2)];
dVar = [dVar ; cinPos2007_dObs(: , 3)];

%% CIN1 prevalence among women aged 20-50 in 2010 by HIV status 
cin1_2010 = zeros(2 , 1);
yr = 2010;
dVec = {[1:2],[3:8]};
for dV = 1 : length(dVec)
    d = dVec{dV};
    
    cinInds = unique([toInd(allcomb(d , 1 : viral , 3 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 5 : 10 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
        [1 : 5 , 7] , 3 , 1 , 1 : intervens , 2 , 5 : 10 , 1 : risk))]);
    ageInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 5 : 10 , 1 : risk));
    cin1_2010(dV) = (sum(popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
        ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
end

mObs = [mObs ; cin1_2010];
dMean = [dMean ; cin1_2010_dObs(: , 2)];
dVar = [dVar ; cin1_2010_dObs(: , 3)];

%% CIN2/3 prevalence among women aged 20-50 in 2010 by HIV status 

yr = 2010;
dVec = {[1:2],[3:8]};
for dV = 1 : length(dVec)
    d = dVec{dV};
    
    cinInds = unique([toInd(allcomb(d , 1 : viral , 4:5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 5 : 10 , 1 : risk)); toInd(allcomb(d , 1 : viral , ...
        [1 : 5 , 7] , 4:5 , 1 , 1 : intervens , 2 , 5 : 10 , 1 : risk))]);
    ageInds = toInd(allcomb(d , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , 5 : 10 , 1 : risk));
    cin2_2010(dV) = (sum(popVec((yr - startYear) * stepsPerYear +1 , cinInds)))...
        ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
end

mObs = [mObs ; cin2_2010'];
dMean = [dMean ; cin2_2010_dObs(: , 2)];
dVar = [dVar ; cin2_2010_dObs(: , 3)];

%% HPV Prevalence in high risk women
%hpv_hiv = zeros(4 , 1);
yr = 2006;
ageVec = {[4:5],[6],[7:8],[9:10]};
for aV = 1 : length(ageVec)
    a = ageVec{aV};
    %HIV-positive
    hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 1 , a , 3)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 1 , a , 3))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 , a , 3));
    hpv_hiv(aV) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear+1 , ageInds));
    %HIV-negative
    hpvInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 1 , a , 3)); toInd(allcomb(1:2 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5, 1 , 1 : intervens , 1 , a , 3))]);
    ageInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 , a , 3));
    hpv_hivNeg(aV) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear+1 , ageInds));
end

mObs = [mObs ; hpv_hiv'; hpv_hivNeg' ];
dMean = [dMean ; hpv_hiv_dObs(: , 2); hpv_hivNeg_dObs(:, 2)];
dVar = [dVar ; hpv_hiv_dObs(: , 3); hpv_hivNeg_dObs(:, 3)];

%% HIV prevalence by age and sex 
hivYearVec = unique(hivPrevM_dObs(: ,1));
hivAgeM = zeros(8 , length(hivYearVec));
hivAgeF = zeros(7 , length(hivYearVec));

for t = 1 : length(hivYearVec)
    for a = 4 : 10
        hivFInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        artFInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        totFInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hivAgeF(a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivFInds)) ...
            + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artFInds))) ...
            / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totFInds));
    end
end
for t = 1 : length(hivYearVec)
    for a = 4 : 11
        hivMInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        artMInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        totMInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        hivAgeM(a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivMInds)) ...
            + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artMInds))) ...
            / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totMInds));
    end
end
mObs = [mObs ; hivAgeM(:) ; hivAgeF(:)];
dMean = [dMean ; hivPrevM_dObs(: , 2) ; hivPrevF_dObs(: , 2)];
dVar =  [dVar ;  hivPrevM_dObs(: , 3) ; hivPrevF_dObs(: , 3)];

%% HIV Prevalence in men and women 2012
%hivPrevM = zeros(4 , 1);
hivYearVec = unique(hivPrevAll_dObs(: ,1));
for a = 4 : 13
        hivInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        artInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        totInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        hivPrevAll(a - 3) =  (sum(popVec((hivYearVec - startYear) * stepsPerYear +1 , hivInds)) ...
            + sum(popVec((hivYearVec - startYear) * stepsPerYear +1 , artInds))) ...
            / sum(popVec((hivYearVec - startYear) * stepsPerYear +1 , totInds));
end

mObs = [mObs ; hivPrevAll'];
dMean = [dMean ; hivPrevAll_dObs(: , 2)];
dVar = [dVar ; hivPrevAll_dObs(: , 3)];

%% HPV prevalence in all women (not including CIN 2+) by age
%hpv_all = zeros(6 , 1);

for a = 6 : 11 % 15-19 -> 55-65
    hpvInds = unique([toInd(allcomb(1 : disease , 1 : viral , 2 : 3 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 2 : 3 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : disease, 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpv_all(a - 5) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
end

mObs = [mObs ; hpv_all'];
dMean = [dMean ; hpv_all_dObs(: , 2)];
dVar = [dVar ; hpv_all_dObs(: , 3)];

%% HPV prevalence in HIV+ women (including CIN)
yr = 2009;

ageVec = {[4:6], 7, 8, 9, [10:11]};
for aV = 1 : length(ageVec)
    a = ageVec{aV}; % 15-19 ->  55-65
    hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpv_hiv2009(aV) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
end

mObs = [mObs ; hpv_hiv2009'];
dMean = [dMean ; hpv_hiv2009_dObs(: , 2)];
dVar = [dVar ; hpv_hiv2009_dObs(: , 3)];

%% Cervical cancer type distribution
yr = 2011;
ccInds_vax = toInd(allcomb(1 : disease , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
ccInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , 6 , ...
    1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
ccInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
        1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

for i = 1 : length(yr)
    cc_vax = sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , ccInds_vax) , 2)...
        ./ sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);
    cc_nonVax = sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , ccInds_nonVax) , 2)...
        ./ sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);

    mObs = [mObs; cc_vax; cc_nonVax];
    dMean = [dMean ; cc_dist_dObs(: , 2)];
    dVar =  [dVar ; cc_dist_dObs(: , 3)];
end

%% CIN3 type distribution
yr = 2011;
cin3Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 5 , [1 : 5 , 7] , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cin3Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 4 , 7] , 5 , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cin3Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

for i = 1 : length(yr)
    cin3_vax = sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , cin3Inds_vax) , 2)...
        ./ sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);
    cin3_nonVax = sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , cin3Inds_nonVax) , 2)...
        ./ sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);

    mObs = [mObs; cin3_vax; cin3_nonVax];
    dMean = [dMean ; cin3_dist_dObs(: , 2)];
    dVar =  [dVar ; cin3_dist_dObs(: , 3)];
end

%% CIN1 type distribution
yr = 2011;
cin1Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cin1Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 2 , 7] , 3 , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cin1Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
   
for i = 1 : length(yr)   
    cin1_vax = sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , cin1Inds_vax) , 2)...
        ./ sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , cin1Inds_tot) , 2);
    cin1_nonVax = sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , cin1Inds_nonVax) , 2)...
        ./ sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , cin1Inds_tot) , 2);

    mObs = [mObs; cin1_vax; cin1_nonVax];
    dMean = [dMean ; cin1_dist_dObs(: , 2)];
    dVar =  [dVar ; cin1_dist_dObs(: , 3)];
end

%% HPV type distribution
yr = 2011;
hpvInds_vax = toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
hpvInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 2 , 7] , 3 , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
hpvInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    
for i = 1 : length(yr)     
    hpv_vax = sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , hpvInds_vax) , 2)...
        ./ sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , hpvInds_tot) , 2);
    hpv_nonVax = sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , hpvInds_nonVax) , 2)...
        ./ sum(popVec(((yr(i) - startYear) * stepsPerYear +1) , hpvInds_tot) , 2);

    mObs = [mObs; hpv_vax; hpv_nonVax];
    dMean = [dMean ; hpv_dist_dObs(: , 2)];
    dVar =  [dVar ; hpv_dist_dObs(: , 3)];
end

%% Cervical cancer incidence in 2012 by age
incTimeSpan = [((2012 - startYear) * stepsPerYear +1) : ((2012 - startYear) * stepsPerYear +6)];
fac = 10 ^ 5;

ccInc2012 = zeros(13 , 1);
for a = 4 : age
    allF = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    % Calculate incidence
    ccInc2012(a - 3 , 1) = ...
        (annlz(sum(sum(sum(newCC(incTimeSpan , : , a , :),2),3),4)) ./ ...
        (annlz(sum(popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);
end

mObs = [mObs ; ccInc2012];
dMean = [dMean ; ccInc2012_dObs(: , 2)];
dVar =  [dVar ; ccInc2012_dObs(: , 3)];

%% Population age distribution
popYearVec = unique(popAgeDist_dObs(: ,1));
popProp = zeros(age, length(popYearVec));

for t = 1 : length(popYearVec)
    for a = 1 : age
        popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
        popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
        popProp(a,t) = sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popAge),2) ./ ...
            sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTot),2);
    end
end

mObs = [mObs ; popProp(:)];
dMean = [dMean ; popAgeDist_dObs(: , 2)];
dVar =  [dVar ;  popAgeDist_dObs(: , 3)];

%% Total population size
popYearVec = unique(totPopSize_dObs(: ,1));
popSize = zeros(1 , length(popYearVec));

for t = 1 : length(popYearVec)
    popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
    popSize(t) = sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTot),2);
end

mObs = [mObs ; popSize(:)];
dMean = [dMean ; totPopSize_dObs(: , 2)];
dVar =  [dVar ;  totPopSize_dObs(: , 3)];

%% Likelihood function
logL = -(0.5*log(2*pi)) - (0.5.*log(dVar)) - ((0.5.*(1./dVar)).*(mObs-dMean).^2);
negSumLogL = sum(logL); % CJB note: despite name, positive summed logL, used to be -sum(logL), the negative logL to be minimized
