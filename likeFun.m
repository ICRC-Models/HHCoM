% Note: likeFun only set up to use 5-year age groups

function negSumLogL = likeFun(popVec , newCC , cinPos2002_dObs , cinNeg2002_dObs ,...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
    hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , ccInc2011_dObs , cc_dist_dObs , ...
    cin3_dist_dObs , cin1_dist_dObs , hpv_dist_dObs , popAgeDist_dObs , toInd , ...
    disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , ...
    age , gender , risk , startYear , stepsPerYear , annlz)

%% CIN2/3 prevalence by HIV status
cinPos2002 = zeros(10 , 1);
cinNeg2002 = cinPos2002;

for a = 4 : 13 % 15-19 -> 60-64
    % HIV-positive (on and not on ART)
    cinInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cinPos2002(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cinInds)))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear +1 , ageInds));
    % HIV-negative
    cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cinNeg2002(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear +1 , cinNegInds)))...
        ./ (sum(popVec((2002 - startYear) * stepsPerYear +1 , ageNegInds)));
end

mObs = [cinPos2002; cinNeg2002];
dMean = [cinPos2002_dObs(: , 2) ; cinNeg2002_dObs(: , 2)];
dVar = [cinPos2002_dObs(: , 3) ; cinNeg2002_dObs(: , 3)];

%% HR HPV Prevalence in HIV+ men
hpv_hivM2008 = zeros(4 , 1);

ageVec = {[4:5],[6:7],[8:9],[10:13]};
for aV = 1 : length(ageVec)
    a = ageVec{aV};
    hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 , [1 : 2 , 7] , ...
        1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    hpv_hivM2008(aV) = sum(popVec((2008 - startYear) * stepsPerYear +1 , hpvInds))...
        ./ sum(popVec((2008 - startYear) * stepsPerYear+1 , ageInds));
end

mObs = [mObs ; hpv_hivM2008];
dMean = [dMean ; hpv_hivM2008_dObs(: , 2)];
dVar = [dVar ; hpv_hivM2008_dObs(: , 3)];

%% HR HPV Prevalence in HIV- men
hpv_hivMNeg2008 = zeros(4 , 1);

ageVec = {[4:5],[6:7],[8:9],[10:13]};
for aV = 1 : length(ageVec)
    a = ageVec{aV};
    hpvInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 , [1 : 2 , 7] , ...
        1 , 1 : intervens , 1 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 2 , 7] , 2 , 1 , 1 : intervens , 1 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
    hpv_hivMNeg2008(aV) = sum(popVec((2008 - startYear) * stepsPerYear +1 , hpvInds))...
        ./ sum(popVec((2008 - startYear) * stepsPerYear +1 , ageInds));
end

mObs = [mObs ; hpv_hivMNeg2008];
dMean = [dMean ; hpv_hivMNeg2008_dObs(: , 2)];
dVar = [dVar ; hpv_hivMNeg2008_dObs(: , 3)];

%% HPV prevalence in HIV- women (including CIN)
yr = 2002;
hpv_hivNeg = zeros(9 , 1);

for a = 4 : 12 % 15-19 ->  55-65
    hpvInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpv_hivNeg(a - 3) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
end

mObs = [mObs ; hpv_hivNeg];
dMean = [dMean ; hpv_hivNeg_dObs(: , 2)];
dVar = [dVar ; hpv_hivNeg_dObs(: , 3)];

%% HPV prevalence in HIV+ women (including CIN)
hpv_hiv = zeros(9 , 1);

for a = 4 : 12 % 15-19 -> 55-65
    hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpv_hiv(a - 3) = sum(popVec((yr - startYear) * stepsPerYear +1 , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear +1 , ageInds));
end

mObs = [mObs ; hpv_hiv];
dMean = [dMean ; hpv_hiv_dObs(: , 2)];
dVar = [dVar ; hpv_hiv_dObs(: , 3)];

%% Cervical cancer type distribution
typeDistYearVec = [2011 , 2012 , 2013 , 2014 , 2015];
ccInds_vax = toInd(allcomb(1 : disease , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
    1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
ccInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 5 , 7] , 6 , ...
    1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk));
ccInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 6 , 1 : hpvNonVaxStates , ...
        1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 5 , 7] , 6 , 1 : 3 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

for i = 1 : length(typeDistYearVec)
    cc_vax = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_vax) , 2)...
        ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);
    cc_nonVax = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_nonVax) , 2)...
        ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , ccInds_tot) , 2);

    mObs = [mObs; cc_vax; cc_nonVax];
    dMean = [dMean ; cc_dist_dObs(: , 2)];
    dVar =  [dVar ; cc_dist_dObs(: , 3)];
end

%% CIN3 type distribution
typeDistYearVec = [2011 , 2012 , 2013 , 2014 , 2015];
cin3Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 5 , [1 : 5 , 7] , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cin3Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 4 , 7] , 5 , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cin3Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 4 , 7] , 5 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);

for i = 1 : length(typeDistYearVec)
    cin3_vax = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_vax) , 2)...
        ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);
    cin3_nonVax = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_nonVax) , 2)...
        ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin3Inds_tot) , 2);

    mObs = [mObs; cin3_vax; cin3_nonVax];
    dMean = [dMean ; cin3_dist_dObs(: , 2)];
    dVar =  [dVar ; cin3_dist_dObs(: , 3)];
end

%% CIN1 type distribution
typeDistYearVec = [2011 , 2012 , 2013 , 2014 , 2015];
cin1Inds_vax = toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cin1Inds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 2 , 7] , 3 , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
cin1Inds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
   
for i = 1 : length(typeDistYearVec)   
    cin1_vax = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_vax) , 2)...
        ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_tot) , 2);
    cin1_nonVax = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_nonVax) , 2)...
        ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , cin1Inds_tot) , 2);

    mObs = [mObs; cin1_vax; cin1_nonVax];
    dMean = [dMean ; cin1_dist_dObs(: , 2)];
    dVar =  [dVar ; cin1_dist_dObs(: , 3)];
end

%% HPV type distribution
typeDistYearVec = [2011 , 2012 , 2013 , 2014 , 2015];
hpvInds_vax = toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
hpvInds_nonVax = toInd(allcomb(1 : disease , 1 : viral , [1 : 2 , 7] , 3 , ...
    1 , 1 : intervens , 2 , 1 : age , 1 : risk));
hpvInds_tot = unique([toInd(allcomb(1 : disease , 1 : viral , 3 , [1 : 3 , 7] , ...
        1 , 1 : intervens , 2 , 1 : age , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , ...
        [1 : 2 , 7] , 3 , 1 , 1 : intervens , 2 , 1 : age , 1 : risk))]);
    
for i = 1 : length(typeDistYearVec)     
    hpv_vax = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_vax) , 2)...
        ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_tot) , 2);
    hpv_nonVax = sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_nonVax) , 2)...
        ./ sum(popVec(((typeDistYearVec(i) - startYear) * stepsPerYear +1) , hpvInds_tot) , 2);

    mObs = [mObs; hpv_vax; hpv_nonVax];
    dMean = [dMean ; hpv_dist_dObs(: , 2)];
    dVar =  [dVar ; hpv_dist_dObs(: , 3)];
end

%% Cervical cancer incidence in 2011 applied to 2009 --> CJB note: need to change year and add more years to this, not functional
% incTimeSpan = [((2009 - startYear) * stepsPerYear +1) : ((2009 - startYear) * stepsPerYear +6)];
% fac = 10 ^ 5;
% worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
%     247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
%     23477 9261 2155];
% 
% ccIncRef = 0.0;
% for aInd = 1:age+4
%     if aInd >= age
%         a = age;
%     end
%     allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
%         1 : periods , 2 , a , 1 : risk)); ...
%         toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
%         1 : periods , 2 , a , 1 : risk))];
%     % Calculate incidence
%     if aInd <= age    
%         ccIncRefA = ...
%             (annlz(sum(sum(sum(newCC(incTimeSpan , : , : , a),2),3),4)) ./ ...
%             (annlz(sum(popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac) ...
%             .* (worldStandard_WP2015(aInd));
%         if (i == 4) && (a < 3) && (max(annlz(sum(sum(sum(newCC(incTimeSpan , : , : , a),2),3),4))) == 0.0)
%             ccIncRefA = 0.0;
%         end
%     elseif aInd > age
%         ccIncRefA = ...
%             (annlz(sum(sum(sum(newCC(incTimeSpan , : , : , a),2),3),4)) ./ ...
%             (annlz(sum(popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);
%         ccIncRefA = [(ones(1,aInd-a).*ccIncRefA(1,1)) , ccIncRefA(1,1:end-(aInd-a))];
%         ccIncRefA = ccIncRefA .* worldStandard_WP2015(aInd);
%     end
%     ccIncRef = ccIncRef + ccIncRefA;
% end
% ccInc = ccIncRef ./ (sum(worldStandard_WP2015(1:age+4)));
% 
% mObs = [mObs ; ccInc];
% dMean = [dMean ; 43.0];
% dVar =  [dVar ; ___];

%% Cervical cancer incidence in 2011 by age
incTimeSpan = [((2011 - startYear) * stepsPerYear +1) : ((2011 - startYear) * stepsPerYear +6)];
fac = 10 ^ 5;

ccInc2011 = zeros(13 , 1);
for a = 4 : age
    allF = toInd(allcomb(1 : disease , 1 : viral , [1:5,7] , [1:5,7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk));
    % Calculate incidence
    ccInc2011(a - 3 , 1) = ...
        (annlz(sum(sum(sum(newCC(incTimeSpan , : , a , :),2),3),4)) ./ ...
        (annlz(sum(popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);
end

mObs = [mObs ; ccInc2011];
dMean = [dMean ; ccInc2011_dObs(: , 2)];
dVar =  [dVar ; ccInc2011_dObs(: , 3)];

%% HIV prevalence
hivYearVec = unique(hivPrevM_dObs(: ,1));
hivAgeM = zeros(7 , length(hivYearVec));
hivAgeF = hivAgeM;
for t = 1 : length(hivYearVec)
    for a = 4 : 10
        hivMInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        artMInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        totMInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 , a , 1 : risk));
        hivFInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        artFInds = toInd(allcomb(8 , 6 , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        totFInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
        hivAgeM(a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivMInds)) ...
            + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artMInds))) ...
            / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totMInds));
        hivAgeF(a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , hivFInds)) ...
            + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , artFInds))) ...
            / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear +1 , totFInds));
    end
end
mObs = [mObs ; hivAgeM(:) ; hivAgeF(:)];
dMean = [dMean ; hivPrevM_dObs(: , 2) ; hivPrevF_dObs(: , 2)];
dVar =  [dVar ;  hivPrevM_dObs(: , 3) ; hivPrevF_dObs(: , 3)];

%% Population age distribution
popYearVec = unique(popAgeDist_dObs(: ,1));
popProp2019 = zeros(age, length(popYearVec));
for t = 1 : length(popYearVec)
    for a = 1 : age
        popAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , a , 1 : risk));
        popTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
            1 : endpoints , 1 : intervens , 1 : gender , 1 : age , 1 : risk));
        popProp2019(a,t) = sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popAge),2) ./ ...
            sum(popVec(((popYearVec(t) - startYear) * stepsPerYear +1) , popTot),2);
    end
end
mObs = [mObs ; popProp2019(:)];
dMean = [dMean ; popAgeDist_dObs(: , 2)];
dVar =  [dVar ;  popAgeDist_dObs(: , 3)];

%% Likelihood function
logL = -(0.5*log(2*pi)) - (0.5.*log(dVar)) - ((0.5.*(1./dVar)).*(mObs-dMean).^2);
negSumLogL = sum(logL); % CJB note: despite name, positive summed logL, used to be -sum(logL), the negative logL to be minimized
