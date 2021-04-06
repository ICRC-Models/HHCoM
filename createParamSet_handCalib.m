function [] = createParamSet_handCalib(tstep_abc , date_abc)
t_curr = tstep_abc;
date = date_abc;

%% Save parameter index values
pIdx = [5, 10];    % save indices into paramsAll cell array for parameters you want to include in calibration
file = ['pIdx_calib_' , date , '_' , num2str(t_curr) , '.dat'];
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , pIdx)

%% Create a matrix of trial parameter sets (size: [# parameters x # trials])
paramSetNew = [0.2, 0.3, 0.4; ... % condUse (index in paramsAll = 5)
               0.005, 0.01, 0.015]; % perPartnerHpv_vax (index in paramsAll = 10)

%% Save matrix of trial parameter sets
file = ['paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat'];
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , paramSetNew)

