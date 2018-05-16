function vaxCEA()

%% load results
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'])
noV = load('H:\HHCoM_Results\CEA_VaxCover_0_Eff_0.9.mat');
c90_2vFull = load('H:\HHCoM_Results\CEA_VaxCover_0.9_Eff_0.7.mat'); 
c70_2vFull = load('H:\HHCoM_Results\CEA_VaxCover_0.7_Eff_0.7.mat'); 
c70_9vFull = load('H:\HHCoM_Results\CEA_VaxCover_0.7_Eff_0.9.mat'); 
c70_2vPartial = load('H:\HHCoM_Results\CEA_VaxCover_0.7_Eff_0.63.mat');
c90_9vFull = load('H:\HHCoM_Results\CEA_VaxCover_0.9_Eff_0.9.mat'); 
c90_2vPartial = load('H:\HHCoM_Results\CEA_VaxCover_0.9_Eff_0.63.mat');
tVec = c90_2vFull.tVec;
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); 
% midMat = zeros(stepsPerYear , size(o90.popVec , 1) / stepsPerYear);
% midMat(1 , :) = 1;
% midMat(end , :) = 1;
% midAnn = @(x) sum(midMat .* reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)) / 2;
c = fix(clock);
currYear = c(1); % get the current year

%% Calculate life years saved

% No intervention case
basePopTot = sum(noV.popVec , 2);

c90_2vFull.lys = sum(c90_2vFull.popVec , 2) - basePopTot;
c70_2vFull.lys = sum(c70_2vFull.popVec , 2) - basePopTot;
c70_9vFull.lys = sum(c70_9vFull.popVec , 2) - basePopTot;
c70_2vPartial.lys = sum(c70_2vPartial.popVec , 2) - basePopTot;
c90_9vFull.lys = sum(c90_9vFull.popVec , 2) - basePopTot;
c90_2vPartial.lys = sum(c90_2vPartial.popVec , 2) - basePopTot;

%% Calculate annual cost of vaccination
cost2v = 27; % cost of 2 doses of bivalent vaccine
c90_2vFull.vaxCost = annlz(c90_2vFull.vaxd) * cost2v;
% c70_2vFull.vaxCost = c70_2vFull.vaxd;
% c70_9vFull.vaxCost = sc70_9vFull.vaxd;
% c70_2vPartial.vaxCost = c70_2vPartial.vaxd;
% c90_9vFull.vaxCost = c90_9vFull.vaxd;
% c90_2vPartial.vaxCost = c90_2vPartial.vaxd;

% NPV of vaccination cost
discountRate = 0.03; % discount rate of 3% per annum
c90_2vFull.vaxCostNPV = pvvar(c90_2vFull.vaxCost , discountRate);
%% Find price threshold



