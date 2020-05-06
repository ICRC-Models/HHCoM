% Identify parameter ranges of top fifty best-fitting sets and compare to bounds
function [] = idParamRanges(tstep_abc , date_abc)
t_curr = tstep_abc;
date = date_abc;

%% Load all particles
paramDir = [pwd , '/Params/'];
masterSetMatrix = load([paramDir , 'masterSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent parameter sample
pIdx = load([paramDir , 'pIdx_calib_' , date , '_0' , '.dat']); % load parameter indices
orderedLL = load([paramDir , 'orderedLL_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent ordered log-likelihoods
paramSetCurr = load([paramDir , 'paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load next parameter sample
paramSetNext = load([paramDir , 'paramSets_calib_' , date , '_' , num2str(t_curr+1) , '.dat']); % load next parameter sample

%% Get indices and parameter values of top fifty sets --> as applied now, actually all accepted particles
numSets60 = size(orderedLL , 1)*0.60;
top50Inds = orderedLL(1:numSets60,1);
top50Params = masterSetMatrix(:,top50Inds); %masterSetMatrix(:,top50Inds);

%% Find min and max of each parameter and compare to bounds
[paramsAll] = genParamStruct();
paramsSub = cell(length(pIdx),1);
lb = [];
ub = [];
for s = 1 : length(pIdx)
    paramsSub{s} = paramsAll{pIdx(s)};
    lb = [lb; paramsSub{s}.lb];
    ub = [ub; paramsSub{s}.ub];
end

minVals = zeros(size(top50Params,1),1);
maxVals = minVals;
lbVals = minVals;
ubVals = minVals;
for i = 1 : size(top50Params,1)
    minVals(i,1) = min(top50Params(i,:));
    maxVals(i,1) = max(top50Params(i,:));
    lbVals(i,1) = lb(i);
    ubVals(i,1) = ub(i);
end

%% Print ranges
% fileRanges = ['rangesVSbounds_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save ranges compared to bounds
% csvwrite([paramDir, fileRanges] , [minVals , maxVals , lbVals, ubVals]);

%% Plot intermediate parameter posterior distributions
% ***Note: may want to change to top 100 parameter sets for this plot***
figure;
% titles = {'partnersM:10-14,hr' , 'partnersM:15-19,hr' , 'partnersM:20-24,hr' , ...    % 10Jan20 calibration
%          'partnersM:25-29,hr' , 'partnersM:30-44,hr' , 'partnersM:45-79,hr' , ...
%          'partnersM:10-14,mr' , 'partnersM:10-14,lr' , ...
%          'partnersF:10-14,hr' , 'partnersF:15-19,hr' , 'partnersF:20-24,hr' , ...
%          'partnersF:25-29,hr' , 'partnersF:30-44,hr' , 'partnersF:45-79,hr' , ...
%          'partnersF:10-14,mr' , 'partnersF:10-14,lr' , ...
%          'condUse' , 'epsA' , 'perPartnerHpv' , 'lambdaMultImm' , ...
%          'Inf->CIN1, vax' , 'Inf->CIN1, nv' , 'CIN1->CIN2, vax' , 'CIN1->CIN2, nv' , ...
%          'CIN2->CIN3, vax' , 'CIN2->CIN3, nv' , 'CIN3->CC, vax' , 'CIN3->CC, nv' , ...
%          'Inf->Well, vax' , 'Inf->Well, nv' , 'CIN1->Inf, vax' , 'CIN1->Inf, nv' , ...
%          'CIN2->CIN1, vax' , 'CIN2->CIN1, nv' , 'CIN3->CIN2, vax' , 'CIN3->CIN2, nv' , ...
%          'HIV transmission' , 'fert decline 1' , 'fert decline 2'};
% titles = {'partnersM:10-14,hr' , 'partnersM:15-19,hr' , 'partnersM:20-24,hr' , ...    % 10Feb20 calibration
%          'partnersM:25-29,hr' , 'partnersM:30-44,hr' , 'partnersM:45-79,hr' , ...
%          'partnersM:10-14,mr' , 'partnersM:10-14,lr' , ...
%          'partnersF:10-14,hr' , 'partnersF:15-19,hr' , 'partnersF:20-24,hr' , ...
%          'partnersF:25-29,hr' , 'partnersF:30-44,hr' , 'partnersF:45-79,hr' , ...
%          'partnersF:10-14,mr' , 'partnersF:10-14,lr' , ...
%          'condUse' , 'epsA' , 'perPartnerHpv' , 'lambdaMultImm' , 'Inf->Well, vax' , 'Inf->Well, nv' , ...
%          'HIV transmission' , 'fert decline 1' , 'fert decline 2' , 'maleHpvClearMult'};
% titles = {'partnersM:10-14,hr' , 'partnersM:15-19,hr' , 'partnersM:20-24,hr' , ...    % 24Feb20 calibration
%          'partnersM:25-29,hr' , 'partnersM:30-44,hr' , 'partnersM:45-79,hr' , ...
%          'partnersM:10-14,mr' , 'partnersM:10-14,lr' , ...
%          'partnersF:10-14,hr' , 'partnersF:15-19,hr' , 'partnersF:20-24,hr' , ...
%          'partnersF:25-29,hr' , 'partnersF:30-44,hr' , 'partnersF:45-79,hr' , ...
%          'partnersF:10-14,mr' , 'partnersF:10-14,lr' , ...
%          'condUse' , 'epsA' , 'perPartnerHpv' , 'lambdaMultImm' , ...
%          'Inf->CIN1, vax' , 'Inf->CIN1, nv' , 'CIN1->CIN2, vax' , 'CIN1->CIN2, nv' , ...
%          'CIN2->CIN3, vax' , 'CIN2->CIN3, nv' , 'CIN3->CC, vax' , 'CIN3->CC, nv' , ...
%          'Inf->Well, vax' , 'Inf->Well, nv' , 'CIN1->Inf, vax' , 'CIN1->Inf, nv' , ...
%          'CIN2->CIN1, vax' , 'CIN2->CIN1, nv' , 'CIN3->CIN2, vax' , 'CIN3->CIN2, nv' , ...
%          'HIV transmission' , 'fert decline 1' , 'fert decline 2' , 'maleHpvClearMult'};
titles = {'partnersM:15-19,hr' , 'partnersM:20-24,hr' , ...    % 22Apr20 calibration
         'partnersM:25-29,hr' , 'partnersM:30-44,hr' , 'partnersM:45-79,hr' , ...
         'partnersM:10-14,mr' , 'partnersM:10-14,lr' , ...
         'partnersF:15-19,hr' , 'partnersF:20-24,hr' , ...
         'partnersF:25-29,hr' , 'partnersF:30-44,hr' , 'partnersF:45-79,hr' , ...
         'partnersF:10-14,mr' , 'partnersF:10-14,lr' , ...
         'condUse' , 'epsA' , ...
         'femaleActs:15-19,lr' , 'femaleActs:20-24,lr' , 'femaleActs:25-29,lr' , ...
         'femaleActs:30-44,lr' , 'femaleActs:45-79,lr' , ...
         'perPartnerHpv' , 'lambdaMultImm' , ...
         'Inf->Well, vax' , 'Inf->Well, nv' , ...
         'HIV transmission' , 'maleHpvClearMult'};
for i = 1 : size(top50Params,1)
    subplot(4,7,i);
    [N2,edges2] = histcounts(paramSetCurr(i,:) , 100 , 'Normalization' , 'probability');
    edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
    plot(edges2 , N2);
    %[f_curr,xi_curr] = ksdensity(paramSetCurr(i,:) , linspace(lb(i) , ub(i) , 100));
    %plot(xi_curr,f_curr);
    hold all;
    [N3,edges3] = histcounts(paramSetNext(i,:) , 100, 'Normalization' , 'probability');
    edges3 = edges3(2:end) - (edges3(2)-edges3(1))/2;
    plot(edges3 , N3);
    %[f_next,xi_next] = ksdensity(paramSetNext(i,:) , linspace(lb(i) , ub(i) , 100));
    %plot(xi_next,f_next);
    %hold off;
    
    prior = unifpdf(linspace(lb(i) , ub(i) , 100) , lb(i) , ub(i)); %[lb(i):((ub(i)-lb(i))/length(top50Inds)):ub(i)]
    plot(linspace(lb(i) , ub(i) , 100) , prior , 'k');
    hold all;
    [N,edges] = histcounts(top50Params(i,:) , 100 , 'Normalization' , 'probability');
    edges = edges(2:end) - (edges(2)-edges(1))/2;
    plot(edges , N , 'g');
    %hold all;
    %[f,xi] = ksdensity(top50Params(i,:) , linspace(lb(i) , ub(i) , 100));
    %plot(xi,f);
    hold off;
    xlim([lb(i) ub(i)]);
    %xlim([min(top50Params(i,:)) max(top50Params(i,:))]);
    %ylim([0.0 max(f)*1.25]);
    ylim([0.0 max(N)*1.25]);
    %ylim([0.0 max(N3)*1.25]);
    title(titles{i});
    if i == size(top50Params,1)
        %legend('paramSetCurr' , 'paramSetNext');
        %legend('prior' , 't0' , 't1' , 't2' , 't3');
        legend('paramSetCurr' , 'paramSetNext' , 'prior' , 't0');
    end
end
