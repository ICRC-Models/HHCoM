function idParamPosteriors(tstep_abc , date_abc)

t_curr = tstep_abc;
date = date_abc;

%% Load all particles
paramDir = [pwd , '/Params/'];
masterSetMatrix = load([paramDir , 'masterSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent parameter sample
pIdx = load([paramDir , 'pIdx_calib_' , date , '_0' , '.dat']); % load parameter indices
orderedLL = load([paramDir , 'orderedLL_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent ordered log-likelihoods

%% Get indices and parameter values of numBestFits*2  sets
numBestFits = 25;
specIndsList = [1,2,3,6,8,9,11,12,13,15,20,21,22,26,27,32,34,35,38,39,40,41,42,45,47]; 
top50Inds = orderedLL(specIndsList,1); %orderedLL(1:numSets50,1);
top50Params = masterSetMatrix(:,top50Inds);

%% If on Phase 2 of calibration, uncomment the following to plot resampled Ph1 parameters + Ph2 parameters
pIdx = load([paramDir,'pIdx_calib_' , date , '_0_wPh1Resample.dat']);
masterResampleSubsetMatrix = load([paramDir , 'masterResampleSubsetMatrix_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent Ph1 resampled parameters
masterCombinedPhaseMatrix = [masterResampleSubsetMatrix ; masterSetMatrix];
top50Params = masterCombinedPhaseMatrix(:,top50Inds);

%% Set up paramsSub for indexing into paramSet matrix
[paramsAll] = genParamStruct();
paramsSub = cell(length(pIdx),1);
startIdx = 1;
for s = 1 : length(pIdx)
    paramsSub{s}.length = paramsAll{pIdx(s)}.length;
    paramsSub{s}.inds = (startIdx : (startIdx + paramsSub{s}.length - 1));
    startIdx = startIdx + paramsSub{s}.length;
end

%% Print posteriors (mean [lb,ub])
titles = {'partnersM:15-19,hr' , 'partnersM:20-24,hr' , ...    % 22Apr20 calibration Phase 1 + Phase 2
          'partnersM:25-29,hr' , 'partnersM:30-44,hr' , 'partnersM:45-79,hr' , ...    % Resampled Phase 1
          'partnersM:10-14,mr' , 'partnersM:10-14,lr' , ...
          'partnersF:15-19,hr' , 'partnersF:20-24,hr' , ...
          'partnersF:25-29,hr' , 'partnersF:30-44,hr' , 'partnersF:45-79,hr' , ...
          'partnersF:10-14,mr' , 'partnersF:10-14,lr' , ...
          'condUse' , 'epsA' , ...
          'femaleActs:15-19,lr' , 'femaleActs:20-24,lr' , 'femaleActs:25-29,lr' , ...
          'femaleActs:30-44,lr' , 'femaleActs:45-79,lr' , ...
          'HIV transmission' , ...
          'perPartnerHpv' , 'CIN2->CIN3, HIV+' , 'CIN1->CIN2, HIV+' , 'lambdaMultImm' , ...    % Phase 2
          'Inf->CIN1, vax' , 'Inf->CIN1, nv' , 'CIN1->CIN2, vax' , 'CIN1->CIN2, nv' , ...
          'CIN2->CIN3, vax' , 'CIN2->CIN3, nv' , 'CIN3->CC, vax' , 'CIN3->CC, nv' , ...
          'Inf->Well, vax' , 'Inf->Well, nv' , 'CIN1->Inf, vax' , 'CIN1->Inf, nv' , ...
          'CIN2->CIN1, vax' , 'CIN2->CIN1, nv' , 'CIN3->CIN2, vax' , 'CIN3->CIN2, nv' , ...
          'maleHpvClearMult' , 'CIN3->CIN2, HIV+' , 'CIN2->CIN1, HIV+'};
paramSet = top50Params(:,1:25);
for i = 1 : size(top50Params,1)
    disp([titles{i} , ': ' , num2str(mean(top50Params(i,:) , 2)) , ' [' , num2str(min(top50Params(i,:) , [] , 2)) , ', ' , num2str(max(top50Params(i,:) , [] , 2)) , ']'])
end
