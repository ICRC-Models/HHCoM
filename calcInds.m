% Calculates indices
function calcInds(stepsPerYear)
% Compartment parameters
disease = 10;
viral = 6;
hpvTypes = 4;
hpvStates = 10;
periods = 3;
gender = 2;
age = 16;
risk = 3;

% index retrieval function
k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);

toInd = @(x) (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + (x(: , 6) - 1) * k(5) ...
    + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + (x(: , 3) - 1) * k(2) ...
    + (x(: , 2) - 1) * k(1) + x(: , 1);

sumall = @(x) sum(x(:));
% Indices
%% Load indices
disp('Preparing indices...')
disp('This may take a while...')
%% mixInfect.m indices
mCurr = zeros(age , risk , viral - 1 , 5 * hpvStates * hpvTypes * periods); % 5 HIV+ disease states
fCurr = zeros(age , risk , viral - 1 , 5 * hpvStates * hpvTypes * periods);
mCurrArt = zeros(age , risk , 1 , 1 * hpvStates * hpvTypes * periods); % 1 HIV+ ART disease state
fCurrArt = zeros(age , risk , 1 , 1 * hpvStates * hpvTypes * periods);
for a = 1 : age
    for r = 1 : risk
        for v = 1 : 5 % viral load (up to vl = 6). Note: last index is (viral - 1) + 1. Done to align pop v index with betaHIV v index.
            mCurr(a , r , v , :) = toInd(allcomb(2 : 6 , v , 1 : hpvTypes , ...
                1 : hpvStates , 1 : periods , 1 , a , r)); %mCurrInd(v , a , r));
            fCurr(a , r , v , :) = toInd(allcomb(2 : 6 , v , 1 : hpvTypes , ...
                1 : hpvStates , 1 : periods , 2 , a , r)); %fCurrInd(v , a , r));
        end
        mCurrArt(a , r , 1 , :) = toInd(allcomb(10 , 6 , 1 : hpvTypes , ...
            1 : hpvStates , 1 : periods , 1 , a , r)); %mCurrInd(v , a , r));
        fCurrArt(a , r , 1 , :) = toInd(allcomb(10 , 6 , 1 : hpvTypes , ...
            1 : hpvStates , 1 : periods , 2 , a , r)); %fCurrInd(v , a , r));
    end
end

gar = zeros(gender , age , risk , disease * viral * hpvTypes * hpvStates * periods);

for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            gar(g , a , r , :) = sort(toInd(allcomb(1 : disease , 1 : viral ,...
                1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r)));
        end
    end
end



naive = zeros(gender , age , risk , periods);

for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            naive(g , a , r , :) = sort(toInd(allcomb(1 , 1 , 1 , 1 , 1 : periods , g , a , r)));
        end
    end
end

coInf = zeros(hpvTypes , gender , age , risk , periods);

for h = 1 : hpvTypes
    for g = 1 : gender
        for a = 1 : age
            for r = 1 : risk
                coInf(h , g , a , r , :) = sort(toInd(allcomb(2 , 2 , h , 2 , 1 : periods , g , a , r)));
            end
        end
    end
end

hivSus = zeros(disease , gender , age , risk , hpvStates * hpvTypes * periods);

for d = 1 : disease
    for g = 1 : gender
        for a = 1 : age
            for r = 1 : risk
                hivSus(d , g , a , r , :) =...
                    sort(toInd(allcomb(d , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r)));
            end
        end
    end
end

hpvSus = zeros(disease, hpvTypes , gender , age , risk , viral * periods * 7);
hpvImm = zeros(disease , gender , age , risk , viral * periods);
hpvVax = hpvImm;
for d = 1 : disease
    for g = 1 : gender
        for a = 1 : age
            for r = 1 : risk
                for h = 1 : hpvTypes
                    hpvSus(d , h , g , a , r , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , h , 1 : 7 , 1 : periods , g , a , r)));
                    hpvImm(d , h , g , a , r , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , h , 10 , 1 : periods , g , a , r)));
                    hpvVaxd(d , h , g , a , r , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , h , 9 , 1 : periods , g , a , r)));
                end
            end
        end
    end
end

toHiv = zeros(gender , age , risk , hpvTypes * hpvStates * periods);
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            toHiv(g , a , r , :) = ...
                sort(toInd(allcomb(2 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r)));
        end
    end
end

toHpv = zeros(disease , hpvTypes , gender , age , risk , viral * periods * 7);

for d = 1 : disease
    for g = 1 : gender
        for a = 1 : age
            for h = 1 : hpvTypes
                for r = 1 : risk
                    toHpv(d , h , g , a , r , :) = ...
                        sort(toInd(allcomb(d , 1 : viral , h , 1 : 7 , 1 : periods , g , a , r)));
                end
            end
        end
    end
end

hrInds = zeros(gender , age , risk , disease * viral * 9 * periods);
lrInds = hrInds;
hrlrInds = hrInds;

for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            hrInds(g , a , r , :) = ...
                sort(toInd(allcomb(1 : disease , 1 : viral , 2 , ...
                1 : 9 , 1 : periods , g , a , r)));
            lrInds(g , a , r , :) = ...
                sort(toInd(allcomb(1 : disease , 1 : viral , 3 , ...
                1 : 9 , 1 : periods , g , a , r)));
            hrlrInds(g , a , r , :) = ...
                sort(toInd(allcomb(1 : disease , 1 : viral , 4 , ...
                1 : 9 , 1 : periods , g , a , r)));
        end
    end
end

save('mixInfectIndices' , 'naive' , 'coInf' , 'hivSus' , 'hpvSus' , 'toHiv' , ...
    'toHpv' , 'mCurr' , 'fCurr' , 'mCurrArt' , 'fCurrArt' , 'gar' , 'hrInds' , ...
    'lrInds' , 'hrlrInds' , 'hpvImm' , 'hpvVaxd')
disp('mixInfect indices loaded')
%% hiv.m , treatDist.m indices
% hivInds(d , v , g , a , r) = makeVec('d' , 'v' , 1 : hpvTypes , 1 : hpvStates , ...
%     1 : periods , 'g' , 'a' , 'r');
% hivInds = matlabFunction(hivInds);

hivInds = zeros(disease , viral , gender , age , risk , hpvTypes * hpvStates * periods);
for d = 1 : disease
    for v = 1 : viral
        for g = 1 : gender
            for a = 1 : age
                for r = 1 :risk
                    hivInds(d , v , g , a , r , :) = ...
                        sort(toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , ...
                        1 : periods , g , a , r)));
                end
            end
        end
    end
end

save('hivIndices' , 'hivInds')
disp('treatDist indices loaded')
%% vlAdv.m indices
% vlInds(d , v , g) = makeVec(d , v , 1 : hpvTypes , 1 : hpvStates , ...
%     1 : periods , g , 1 : age , 1 : risk);
% vlInds = matlabFunction(vlInds , 'Optimize' , false);
% save('vlIndices' , 'vlInds')
%% bornDie.m indices
% hivInd(d , v , g , a) = makeVec('d' , 'v' , 1 : hpvTypes , 1 : hpvStates , ...
%     1 : periods , 'g' , 'a' , 1 : risk);
% hivInd = matlabFunction(hivInd);
% birthsInds(d , g , r) = makeVec('d' , 1 , 1 , 1 , 1 , 'g' , 1 , 'r');
% birthsInds = matlabFunction(birthsInds);
% % deathInds(g , a) = allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
% %     1 : periods , g , a , 1 : risk);
% % deathInds = matlabFunction(deathInds , 'Optimize' , false);
% vaxInds(h , s , g , a) = makeVec(1 : disease , 1 : viral , 'h' , 's' , ...
%     1 : periods , 'g' , 'a' , 1 : risk);
% vaxInds = matlabFunction(vaxInds);
% save('bornDieIndices' , 'hivInd' , 'birthsInds' , 'vaxInds')
% disp('bornDie indices loaded')

%% hpv.m indices
load('hpvData')
disp('Preparing indices for HPV modules...')
disp('This might take a while...')


ccInds = zeros(disease , viral , hpvTypes , age , risk);
ccRegInds = ccInds;
ccDistInds = ccInds;
cin1Inds = ccInds;
cin2Inds = ccInds;
cin3Inds = ccInds;

for d = 1 : disease
    for v = 1 : viral
        for h = 2 : hpvTypes
            for a = 1 : age
                ccInds(d , v , h , a , :) = sort(toInd(allcomb(d , v , h , 5 , 1 , 2 , a , 1 : risk)));
                ccRegInds(d , v , h , a, :) = sort(toInd(allcomb(d , v , h , 6 , 1 , 2 , a , 1 : risk)));
                ccDistInds(d , v , h , a , :) = sort(toInd(allcomb(d , v , h , 7 , 1 , 2 , a , 1 : risk)));
                cin1Inds(d , v , h , a , :) = sort(toInd(allcomb(d , v , h , 2 , 1 , 2 , a , 1 : risk)));
                cin2Inds(d , v , h , a , :) = sort(toInd(allcomb(d , v , h , 3 , 1 , 2 , a , 1 : risk)));
                cin3Inds(d , v , h , a , :) = sort(toInd(allcomb(d , v , h , 4 , 1 , 2 , a , 1 : risk)));
            end
        end
    end
end

screen35PlusInds = zeros(disease , viral , hpvTypes , (age - 8 + 1) * risk);
screen25_35Inds = zeros(disease , viral , hpvTypes , 2 * risk);

for d = 1 : disease
    for v = 1 : viral
        for h = 2 : hpvTypes
            screen35PlusInds(d , v , h , :) = sort(toInd(allcomb(d , v , h , 4 , 1 , 2 , 8 : age , 1 : risk)));
            screen25_35Inds(d , v , h , :) = sort(toInd(allcomb(d , v , h , 4 , 1 , 2 , 6 : 7 , 1 : risk)));
        end
    end
end


ccRInds = zeros(disease , viral , hpvTypes , hpvStates , periods , age * risk);
cc2SusInds = zeros(disease , viral , age * risk);


for d = 1 : disease
    for v = 1 : viral
        cc2SusInds(d , v , :) = sort(toInd(allcomb(d , v , 1 , 1 , 1 , 2 , 1 : age , 1 : risk)));
        for h = 2 : hpvTypes
            for s = 1 : hpvStates
                for p = 1 : periods
                    ccRInds(d , v , h , s , p , :) = ...
                        sort(toInd(allcomb(d , v , h , s , p , 2 , 1 : age , 1 : risk)));
                end
            end
        end
    end
end

normalInds = zeros(disease , viral , gender , age , risk);
infInds = zeros(disease , viral , hpvTypes , gender , age , risk);
immuneInds = infInds;
for g = 1 : gender
    for d = 1 : disease
        for v = 1 : viral
            for a = 1 : age
                normalInds(d , v , g , a , :) = ...
                    sort(toInd(allcomb(d , v , 1 , 1 , 1 , g , a , 1 : risk)));
                for h = 2 : hpvTypes
                    immuneInds(d , v , h , g , a , :) = ...
                        sort(toInd(allcomb(d , v , h , 10 , 1 , g , a , 1 : risk)));
                    infInds(d , v , h , g , a , :) = ...
                        sort(toInd(allcomb(d , v , h , 1 , 1 , g , a , 1 : risk)));
                end
            end
        end
    end
end

save('hpvIndices' , 'infInds' , 'cin1Inds' , 'cin2Inds' , 'cin3Inds' , 'normalInds' , ...
    'ccRInds' , 'screen35PlusInds' , 'screen25_35Inds' , 'ccInds' , 'ccRegInds' , 'ccDistInds' ,'immuneInds')
disp('hpv indices loaded')
%% hpvTreat.m indices
ccRInds = zeros(disease , viral , hpvTypes , hpvStates , periods , age , risk);
ccSusInds = zeros(disease , viral , hpvStates , age , risk);

for d = 1 : disease
    for v = 1 : viral
        for h = 2 : hpvTypes
            for s = 5 : 7 % 5 - 7 cervical cancer
                for a = 1 : age
                    for p = 1 : periods
                        ccRInds(d , v , h , s , p , a , :) = toInd(allcomb(d , v , h , s , p , 2 , a , 1 : risk));
                    end
                end

            end
        end
    end
end

for d = 1 : disease
    for v = 1 : viral
        for h = 1 : hpvTypes
            for a = 1 : age
                ccSusInds(d , v , h , a , :) = toInd(allcomb(d , v , h , 1 , 1 , 2 , a , 1 : risk));
            end
        end
    end
end
getHystPopInds = zeros(age , disease * viral * hpvTypes * 7 * risk * periods);
hystPopInds = zeros(age , disease * viral * hpvTypes * periods * risk);
for a = 1 : age
    getHystPopInds(a , :) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : 7 , ...
        1 : periods , 2 , a , 1 : risk));
    hystPopInds(a , :) = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 8 , ...
        1 : periods , 2 , a , 1 : risk));
end

save('hpvTreatIndices' , 'ccRInds' , 'ccSusInds' , 'getHystPopInds' , 'hystPopInds')
disp('hpvTreatIndices loaded')
%% cinAdv.m indices
inf1 = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : hpvStates ,...
    1 , 1 : gender , 1 : age , 1 : risk));
%kAdv = 1 / (stepsPerYear .* pCinSize); % pCinSize is a vector detailing the interval sizes of each CIN period group
kCC = 1 / (stepsPerYear .* pCCSize);% pCCSize is a vector detailing the interval sizes of each CC period group
local = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 5 , 1 : periods ,...
        1 : gender , 1 : age , 1 : risk));
regional = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 6 , 1 : periods ,...
        1 : gender , 1 : age , 1 : risk));
distant = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 7 , 1 : periods ,...
        1 : gender , 1 : age , 1 : risk));
save('cinAdvData' , 'kCC')
save('cinAdvIndices' , 'inf1' , 'local' , 'regional' , 'distant')
disp('cinAdv indices loaded')
%% agePop.m indices
% Pre-calculate agePop indices
% genPrevInds = zeros(gender , age - 1 , risk , 7 * viral * hpvTypes * hpvStates * periods);
% prepPrevInds = zeros(gender , age - 1 , risk , 2 * viral * hpvTypes * hpvStates * periods);
% artPrevInds = zeros(gender , age - 1 , risk , 1 * viral * hpvTypes * hpvStates * periods);
%
% genNextInds = zeros(gender , age - 1 , risk , 7 * viral * hpvTypes * hpvStates * periods);
% prepNextInds = zeros(gender , age - 1 , risk , 2 * viral * hpvTypes * hpvStates * periods);
% artNextInds = zeros(gender , age - 1 , risk , 1 * viral * hpvTypes * hpvStates * periods);
%
% genLastInds = zeros(gender , risk , 7 * viral * hpvTypes * hpvStates * periods);
% prepLastInds = zeros(gender , risk , 2 * viral * hpvTypes * hpvStates * periods);
% artLastInds = zeros(gender , risk , 1 * viral * hpvTypes * hpvStates * periods);
%
% for g = 1 : gender
%     for a = 2 : age
%         for r = 1 : risk
%             %indices for previous generation
%             genPrevInds (g , a - 1 , r , :) =...
%                 toInd(allcomb(1 : 7 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a - 1 , r));
%             prepPrevInds(g , a - 1 , r , :) = ...
%                 toInd(allcomb(8 : 9 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a - 1 , r));
%             artPrevInds(g , a - 1 , r , :) = ...
%                 toInd(allcomb(10 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a - 1 , r));
%             for rr = 1 : risk
%                 % indices for subsequent generation
%                 genNextInds(g , a - 1 , rr , :) = ...
%                     toInd(allcomb(1 : 7 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , rr));
%                 prepNextInds(g , a - 1 , rr , :) = ...
%                     toInd(allcomb(8 : 9 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , rr));
%                 artNextInds(g , a - 1 , rr , :) = ...
%                     toInd(allcomb(10 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , rr));
%             end
%         end
%     end
% end
%
% for g = 1 : gender
%     for r = 1 : risk
%         genLastInds(g , r , :) = ...
%             toInd(allcomb(1 : 7 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , age , r));
%         prepLastInds(g , r , :) = ...
%             toInd(allcomb(8 : 9 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , age , r));
%         artLastInds(g , r , :) = ...
%             toInd(allcomb(10 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , age , r));
%     end
% end
% genPrevInds = sort(genPrevInds);
% prepPrevInds = sort(prepPrevInds);
% artPrevInds = sort(artPrevInds);
% genNextInds = sort(genNextInds);
% prepNextInds = sort(prepNextInds);
% artNextInds = sort(artNextInds);
% genLastInds = sort(genLastInds);
% prepLastInds = sort(prepLastInds);
% artLastInds = sort(artLastInds);
%
% save('ageIndices' , 'genPrevInds' , 'prepPrevInds' , 'artPrevInds' , 'genNextInds' , 'prepNextInds' ,...
%     'artNextInds' , 'genLastInds' , 'prepLastInds' , 'artLastInds');

%% calcDist.m indices
% N/A
disp('Done')
disp('All indices loaded.')
disp(' ')
