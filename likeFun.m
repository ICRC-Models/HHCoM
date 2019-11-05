% Note: likeFun only set up to use 5-year age groups

function negSumLogL = likeFun(popVec , newCC , cinPos2002_dObs , cinNeg2002_dObs ,...
    hpv_hiv_dObs , hpv_hivNeg_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
    hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , ccInc2011_dObs , toInd , ...
    disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , ...
    age , risk , startYear , stepsPerYear)

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

%% HPV Prevalence in HIV+ Women in 2008 (no CIN2/3)
% % hpv_hiv_2008 = zeros(10 , 1);
% % 
% % for a = 4 : 13 % 15-19 -> 60-64
% %     hpvInds = unique([toInd(allcomb(3 : 7 , 1 : viral , 2 , [1 : 2 , 7] , ...
% %         1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 7 , 1 : viral , ...
% %         [1 : 2 , 7] , 2 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
% %     ageInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %         1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
% %     hpv_hiv_2008(a - 3) = sum(popVec((2008 - startYear) * stepsPerYear +1 , hpvInds))...
% %         ./ sum(popVec((2008 - startYear) * stepsPerYear +1 , ageInds));
% % end
% % 
% % mObs = [mObs ; hpv_hiv_2008];
% % dMean = [dMean ; hpv_hiv_2008_dObs(: , 2)];
% % dVar = [dVar ; hpv_hiv_2008_dObs(: , 3)];

%% HPV Prevalence in HIV- Women in 2008 (no CIN2/3)   ****CJB Note: should ART really be included here??????
% % hpv_hivNeg_2008 = zeros(10 , 1);
% % 
% % for a = 4 : 13 % 15-19 -> 60-64
% %     hpvInds = unique([toInd(allcomb([1 : 2 , 8] , 1 : viral , 2 , [1 : 2 , 7] , ...
% %         1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb([1 : 2 , 8] , ...
% %         1 : viral , [1 : 2 , 7] , 2 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
% %     ageInds = toInd(allcomb([1 : 2 , 8] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %         1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
% %     hpv_hivNeg_2008(a - 3) = hpv_hivNeg_2008(a - 3) + sum(popVec((2008 - startYear) * stepsPerYear +1 , hpvInds))...
% %         ./ sum(popVec((2008 - startYear) * stepsPerYear +1 , ageInds));
% % end
% % 
% % mObs = [mObs ; hpv_hivNeg_2008];
% % dMean = [dMean ; hpv_hivNeg_2008_dObs(: , 2)];
% % dVar = [dVar ; hpv_hivNeg_2008_dObs(: , 3)];

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

%% Cervical cancer incidence type distribution --> CJB note: don't use, need to figure out how to represent this as a normal distribution
% newCCTotal = sum(sum(sum(newCC(: , : , : , :) , 2) , 3) , 4);
% newCCType = zeros(size(newCC , 1) , 2);
% for z = 1 : hpvTypeGroups
%     newCCType(: , z) = sum(sum(newCC(: , : , :  , z) , 2) , 3) ./ newCCTotal;
% end
% mObs = [mObs; mean(newCCType(: , 1)); mean(newCCType(: , 2))];
% nPos = [nPos ; 90 ; 10];
% N =  [N ; 100 ; 100];

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

ccInc2011 = zeros(14 , 1)';
for a = 3 : age
    allF = [toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 4 , ...
        1 : periods , 2 , a , 1 : risk)); ...
        toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 : 10 , ...
        1 : periods , 2 , a , 1 : risk))];
    % Calculate incidence
  
    ccInc2011(a - 2 , 1) = ...
        (annlz(sum(sum(sum(newCC(incTimeSpan , : , : , a),2),3),4)) ./ ...
        (annlz(sum(popVec(incTimeSpan , allF) , 2) ./ stepsPerYear)) * fac);
end

mObs = [mObs ; ccInc2011];
dMean = [dMean ; ccInc2011_dObs(: , 2)];
dVar =  [dVar ; ccInc2011_dObs(: , 3)];

%% HIV
hivYearVec = unique(hivPrevM_obs(: ,1));
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

%% Likelihood function
logL = -(0.5*log(2*pi)) - (0.5.*log(dVar)) - ((0.5.*(1./dVar))*(mObs-dMean).^2);
negSumLogL = sum(logL); % CJB note: despite name, positive summed logL, used to be -sum(logL), the negative logL to be minimized
