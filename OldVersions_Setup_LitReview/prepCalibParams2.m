close all;clear all;clc
paramDir = [pwd ,'\Params\'];

load ([paramDir,'calibratedParams'])
load([paramDir,'popData'])
load([paramDir,'HIVParams'])
load([paramDir,'general'])
load([paramDir,'mixInfectParams'])
load([paramDir,'vlBeta'])
load([paramDir,'hpvData'])
load([paramDir ,'cost_weights'])
load([paramDir,'calibData'])
load([paramDir,'calibInitParams\'],'muCC','muCC_det','hivCC','leep','kCin2_Cin1','kCin1_Cin2','kCin2_Cin3')

%% Initial population
mInit = popInit(: , 1);
fInit = popInit(: , 2);

% test!!!!
riskDistF = riskDistM;
partnersF = partnersM;

MpopStruc = riskDistM;
FpopStruc = riskDistF;

mPop = zeros(age , risk);
fPop = mPop;

for i = 1 : age
    mPop(i , :) = MpopStruc(i, :).* mInit(i) ./ 1.25;
    fPop(i , :) = FpopStruc(i, :).* fInit(i) ./ 1.25;
end

dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age ,risk];
initPop = zeros(dim);
initPop(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop;
initPop(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop;
initPop_0 = initPop;

hpv_hivMult = sum(bsxfun(@times , hpv_hivMult , distWeight) , 2);

partnersM(4 , :) = partnersM(4 , :) .* [1.25 , 1.75 , 1.75];
partnersF(4 , :) = partnersF(4 , :) .* [1.25 , 1.75 , 1.75];
partnersM(5 , :) = partnersM(5 , :) .* [1.25 , 1.5 , 1.75];
partnersF(5 , :) = partnersF(5 , :) .* [1.25 , 1.5 , 1.75];

femaleActs(4 : 5 , :) = femaleActs(4 : 5 , :) .* 1.2 ;
femaleActs(6 : 10 , :) = femaleActs(6 : 10 , :) .* 0.9;
maleActs(4 : 5 , :) = maleActs(4 : 5 , :);

for a = 1 : age
    betaHIVF2M(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_F2M , maleActs(a , :)')); % HIV(-) males
    betaHIVM2F(a , : , :) = 1 - (bsxfun(@power, 1 - betaHIV_M2F , femaleActs(a , :)')); % HIV(-) females
end
betaHIVM2F = permute(betaHIVM2F , [2 1 3]); % risk, age, vl
betaHIVF2M = permute(betaHIVF2M , [2 1 3]); % risk, age, vl

kCCDet = min(kCCDet .* 12 , 0.99); % convert monthly to yearly rate

riskDistF = riskDistM;
riskDist(: , : , 1) = riskDistM;
riskDist(: , : , 2) = riskDistF;

save([paramDir , 'calibratedParams'])