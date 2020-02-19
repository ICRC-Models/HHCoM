% Calculates the negative log likelihood of a sample of potential parameter sets
% Accepts:
% 1) paramSetIdx (index in larger sample matrix to start selection of sub-set)
% Saves:
% 1) File: negSumLogL_calib_[date].dat (negative log likelihood for each parameter set in sub-set)

function lhsPatternSrch_prt2(SetIdx , t_curr , dateIn)

% dateIn = '7Oct19';
% 
% tstep = 0;
% t_curr = tstep;

%% Load parameters
paramDir = [pwd ,'\Params\'];
load([paramDir,'settings']);
load([paramDir,'general']);
paramSetMatrix = load([paramDir,'paramSets_patternSrch_', dateIn, '_0.dat']);

pIdx = [1,2,5,6,7,8,9,10];    % indices in paramsAll cell array
file = ['pIdx_patternSrch_29Aug19_0.dat'];
paramDir = [pwd , '\Params\'];
csvwrite([paramDir, file] , pIdx)
%pIdx = load([paramDir,'pIdx_patternSrch_' , dateIn , '_0.dat']);

[paramsAll] = genParamStruct();
paramsSub = cell(length(pIdx),1);
lb = [];
ub = [];
for s = 1 : length(pIdx)
    paramsSub{s} = paramsAll{pIdx(s)};
    lb = [lb; paramsSub{s}.lb];
    ub = [ub; paramsSub{s}.ub];
end
lb(8) = lb(8).*10; % re-scale perPartnerHpv to be more similar in scale to other params
ub(8) = ub(8).*10;

initParams = paramSetMatrix(:,1);

%% Cluster information
pc = parcluster('local');    % create a local cluster object
pc.JobStorageLocation = strcat('/gscratch/csde/guiliu' , '/' , getenv('SLURM_JOB_ID'))    % explicitly set the JobStorageLocation to the temp directory that was created in the sbatch script
parpool(pc , str2num(getenv('SLURM_CPUS_ON_NODE')))    % start the pool with max number workers

%% Obtain model output for each set of sampled parameters
options = psoptimset('UseParallel' , true , 'Cache' , 'on' ,...
    'CacheTol' , 0.000000001 , 'CompletePoll' , 'on' , ...
    'TolMesh' , 0.001 , 'MaxIter' , 5 , ... 
    'Display','iter' , 'PlotFcn' , @psplotbestf);
[optiParams , negSumLogL , exitFlag , output] = patternsearch(@calibratorPtrnSrch, initParams , [] , [] , [] , [] , lb , ub , [] , options);

%% Save parameter sets and negSumLogL values
file = ['negSumLogL_patternSrch_' , dateIn , '_' , num2str(t_curr) , '_params.dat'];
paramDir = [pwd , '\Params\'];
dlmwrite([paramDir, file] , optiParams)

file = ['negSumLogL_patternSrch_' , dateIn , '_' , num2str(t_curr) , '_info.dat'];
paramDir = [pwd , '\Params\'];
dlmwrite([paramDir, file] , [1; negSumLogL; exitFlag] , 'delimiter' , ',' , 'roffset' , 1 , 'coffset' , 0 , '-append')

