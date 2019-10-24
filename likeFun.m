function negSumLogL = likeFun(popVec , newCC , cinPos2002_obs , cinNeg2002_obs ,...
    hpv_hiv_2008_obs , hpv_hivNeg_2008_obs , hpv_hiv_obs , hpv_hivNeg_obs , ...
    hivPrevM_obs , hivPrevF_obs , hpv_hivM2008_obs , hpv_hivMNeg2008_obs , toInd , ...
    disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , intervens , ...
    gender , age , risk , startYear , stepsPerYear)

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
    cinPos2002(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , cinInds)))...
        ./ sum(popVec((2002 - startYear) * stepsPerYear , ageInds)) * 100;
    % HIV-negative
    cinNegInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 4 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 4 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageNegInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    cinNeg2002(a - 3) = (sum(popVec((2002 - startYear) * stepsPerYear , cinNegInds)))...
        ./ (sum(popVec((2002 - startYear) * stepsPerYear , ageNegInds))) * 100;
end

pPos = [cinPos2002; cinNeg2002];
N = [cinPos2002_obs(: , 3) ; cinNeg2002_obs(: , 3)];
nPos = [cinPos2002_obs(: , 2) ; cinNeg2002_obs(: , 2)];

%% HPV Prevalence in HIV+ Women in 2008 (no CIN2/3)
% % hpv_hiv_2008 = zeros(10 , 1);
% % 
% % for a = 4 : 13 % 15-19 -> 60-64
% %     hpvInds = unique([toInd(allcomb(3 : 7 , 1 : viral , 2 , [1 : 2 , 7] , ...
% %         1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 7 , 1 : viral , ...
% %         [1 : 2 , 7] , 2 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
% %     ageInds = toInd(allcomb(3 : 7 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %         1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
% %     hpv_hiv_2008(a - 3) = sum(popVec((2008 - startYear) * stepsPerYear , hpvInds))...
% %         ./ sum(popVec((2008 - startYear) * stepsPerYear , ageInds)) * 100;
% % end
% % 
% % pPos = [pPos ; hpv_hiv_2008];
% % N = [N ; hpv_hiv_2008_obs(: , 3)];
% % nPos = [nPos ; hpv_hiv_2008_obs(: , 2)];

%% HPV Prevalence in HIV- Women in 2008 (no CIN2/3)   ****CJB Note: should ART really be included here??????
% % hpv_hivNeg_2008 = zeros(10 , 1);
% % 
% % for a = 4 : 13 % 15-19 -> 60-64
% %     hpvInds = unique([toInd(allcomb([1 : 2 , 8] , 1 : viral , 2 , [1 : 2 , 7] , ...
% %         1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb([1 : 2 , 8] , ...
% %         1 : viral , [1 : 2 , 7] , 2 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
% %     ageInds = toInd(allcomb([1 : 2 , 8] , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
% %         1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
% %     hpv_hivNeg_2008(a - 3) = hpv_hivNeg_2008(a - 3) + sum(popVec((2008 - startYear) * stepsPerYear , hpvInds))...
% %         ./ sum(popVec((2008 - startYear) * stepsPerYear , ageInds)) * 100;
% % end
% % 
% % pPos = [pPos ; hpv_hivNeg_2008];
% % N = [N ; hpv_hivNeg_2008_obs(: , 3)];
% % nPos = [nPos ; hpv_hivNeg_2008_obs(: , 2)];

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
    hpv_hivM2008(aV) = sum(popVec((2008 - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((2008 - startYear) * stepsPerYear , ageInds)) * 100;
end

pPos = [pPos ; hpv_hivM2008];
N = [N ; hpv_hivM2008_obs(: , 2)];
nPos = [nPos ; hpv_hivM2008_obs(: , 1)];

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
    hpv_hivMNeg2008(aV) = sum(popVec((2008 - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((2008 - startYear) * stepsPerYear , ageInds)) * 100;
end

pPos = [pPos ; hpv_hivMNeg2008];
N = [N ; hpv_hivMNeg2008_obs(: , 2)];
nPos = [nPos ; hpv_hivMNeg2008_obs(: , 1)];

%% HPV prevalence in HIV- women (including CIN)
yr = 2002;
hpv_hivNeg = zeros(9 , 1);

for a = 4 : 12 % 15-19 ->  55-65
    hpvInds = unique([toInd(allcomb(1 : 2 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(1 : 2 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(1 : 2 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpv_hivNeg(a - 3) = sum(popVec((yr - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;
end

pPos = [pPos ; hpv_hivNeg];
N = [N ; hpv_hivNeg_obs(: , 3)];
nPos = [nPos ; hpv_hivNeg_obs(: , 2)];

%% HPV prevalence in HIV+ women (including CIN)
hpv_hiv = zeros(9 , 1);

for a = 4 : 12 % 15-19 -> 55-65
    hpvInds = unique([toInd(allcomb(3 : 8 , 1 : viral , 2 : 5 , [1 : 5 , 7] , ...
        1 , 1 : intervens , 2 , a , 1 : risk)); toInd(allcomb(3 : 8 , 1 : viral , ...
        [1 : 5 , 7] , 2 : 5 , 1 , 1 : intervens , 2 , a , 1 : risk))]);
    ageInds = toInd(allcomb(3 : 8 , 1 : viral , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
        1 : endpoints , 1 : intervens , 2 , a , 1 : risk));
    hpv_hiv(a - 3) = sum(popVec((yr - startYear) * stepsPerYear , hpvInds))...
        ./ sum(popVec((yr - startYear) * stepsPerYear , ageInds)) * 100;
end

pPos = [pPos ; hpv_hiv];
N = [N ; hpv_hiv_obs(: , 3)];
nPos = [nPos ; hpv_hiv_obs(: , 2)];

%% Cervical cancer incidence type distribution    ****CJB note: not updated
% % newCCTotal = sum(sum(sum(newCC(: , : , : , :) , 2) , 3) , 4);
% % newCCType = zeros(size(newCC , 1) , 3);
% % for h = 2 : hpvVaxStates
% %     newCCType(: , h - 1) = sum(sum(newCC(: , : , h  , :) , 2) , 4) ...
% %         ./ newCCTotal;
% % end
% % pPos = [pPos; mean(newCCType(: , 1)); mean(newCCType(: , 2));  mean(newCCType(: , 3))];
% % nPos = [nPos ; 70 ; 20 ; 10];
% % N =  [N ; 100 ; 100 ; 100];

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
        hivAgeM(a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear , hivMInds)) ...
            + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear , artMInds))) ...
            / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear , totMInds)) * 100;
        hivAgeF(a - 3 , t) =  (sum(popVec((hivYearVec(t) - startYear) * stepsPerYear , hivFInds)) ...
            + sum(popVec((hivYearVec(t) - startYear) * stepsPerYear , artFInds))) ...
            / sum(popVec((hivYearVec(t) - startYear) * stepsPerYear , totFInds)) * 100;
    end
end
pPos = [pPos; hivAgeM(:) ; hivAgeF(:)];
nPos = [nPos ; hivPrevM_obs(: , 2) ; hivPrevF_obs(: , 2)];
N =  [N ;  hivPrevM_obs(: , 3) ; hivPrevF_obs(: , 3)];

%% Likelihood function
pPos = pPos ./ 100; % scale percent probabilities to decimals
logL = nPos .* log(pPos) + (N - nPos) .* log(1 - pPos); % log likelihoods for binomial events
negSumLogL = sum(logL); % -sum(logL); % negative logL to be minimized --> CJB: positive logL
