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

%% Get indices and parameter values of top 50 sets (*0.0089), or all accepted particles (*0.40)
numSets60 = size(orderedLL , 1)*0.0089;
top50Inds = orderedLL(1:numSets60,1);
top50Params = masterSetMatrix(:,top50Inds);

%% If on Phase 2 of calibration, uncomment the following to plot resampled Ph1 parameters + Ph2 parameters
% pIdx = load([paramDir,'pIdx_calib_' , date , '_0_wPh1Resample.dat']);
% masterResampleSubsetMatrix = load([paramDir , 'masterResampleSubsetMatrix_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent Ph1 resampled parameters
% masterCombinedPhaseMatrix = [masterSetMatrix(1:27,:); masterResampleSubsetMatrix(1:19,:); masterSetMatrix(28,:); masterResampleSubsetMatrix(20:22,:); masterSetMatrix(29:34,:)];
% top50Params = masterCombinedPhaseMatrix(:,top50Inds);

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

%% Plot intermediate parameter posterior distributions
figure;
titles = {'partnersMmult:15-29,mr' , 'partnersMmult:30-44,mr' , 'partnersMmult:45-79,mr' , ...    % 28Feb21
          'partnersMmult:15-29,hr' , 'partnersMmult:30-49,hr' , 'partnersMmult:50-79,mr' , ...
          'partnersFmult:15-29,mr' , 'partnersFmult:30-44,mr' , 'partnersFmult:45-79,mr' , ...
          'partnersFmult:15-29,hr' , 'partnersFmult:30-44,hr' , 'partnersFmult:45-79,mr' , ...
          'condUse' , 'epsA' , 'epsR' , ...
          'maleActsMult:10-19' , 'maleActsMult:20-24' , 'maleActsMult:25-34' , ...
          'maleActsMult:35-44' , 'maleActsMult:45-79' , ...
          'femaleActsMult:10-19' , 'femaleActsMult:20-24' , 'femaleActsMult:25-34' , ... 
          'femaleActsMult:35-44' , 'femaleActsMult:45-79' , ...
          'HPV transmission:vax' , 'HPV transmission: nonvax' , ...
          'HIV transmission' , ...
          'partnersTimeMult:inc15-24M' , 'partnersTimeMult:inc15-24F' , 'partnersTimeMult:inc25-49MF' , ...
          'partnersTimeMult:dec15-24M' , 'partnersTimeMult:dec15-24M' , 'partnersTimeMult:dec25-49MF'};
%           'CIN2->CIN3, HIV+' , 'CIN1->CIN2, HIV+' , 'lambdaMultImm' , ...
%           'Inf->CIN1, vax' , 'Inf->CIN1, nv' , 'CIN1->CIN2, vax' , 'CIN1->CIN2, nv' , ...
%           'CIN2->CIN3, vax' , 'CIN2->CIN3, nv' , 'CIN3->CC, vax' , 'CIN3->CC, nv' , ...
%           'Inf->Well, vax' , 'Inf->Well, nv' , 'CIN1->Inf, vax' , 'CIN1->Inf, nv' , ...
%           'CIN2->CIN1, vax' , 'CIN2->CIN1, nv' , 'CIN3->CIN2, vax' , 'CIN3->CIN2, nv' , ...
%           'maleHpvClearMult'
%           'CIN3->CIN2, HIV+' , 'CIN2->CIN1, HIV+'
for i = 1 : size(top50Params,1)
    subplot(6,6,i); 
%     hold all; 
%     [N2,edges2] = histcounts(paramSetCurr(i,:) , 100 , 'Normalization' , 'probability');
%     edges2 = edges2(2:end) - (edges2(2)-edges2(1))/2;
%     plot(edges2 , N2);
%     %[f_curr,xi_curr] = ksdensity(paramSetCurr(i,:) , linspace(lb(i) , ub(i) , 100));
%     %plot(xi_curr,f_curr);
%     hold all;
      
%     [N3,edges3] = histcounts(paramSetNext(i,:) , 100, 'Normalization' , 'probability');
%     edges3 = edges3(2:end) - (edges3(2)-edges3(1))/2;
%     plot(edges3 , N3);
%     %[f_next,xi_next] = ksdensity(paramSetNext(i,:) , linspace(lb(i) , ub(i) , 100));
%     %plot(xi_next,f_next);
%     %hold off;
    
%     prior = unifpdf(linspace(lb(i) , ub(i) , 100) , lb(i) , ub(i)); %[lb(i):((ub(i)-lb(i))/length(top50Inds)):ub(i)]
%     plot(linspace(lb(i) , ub(i) , 100) , prior , 'c');
%     hold all;
%     [N,edges] = histcounts(top50Params(i,:) , 8 , 'Normalization' , 'probability');
%     edges = edges(2:end) - (edges(2)-edges(1))/2;
%     p1 = plot(edges , N , 'k' , 'LineWidth' , 1);
%     p1.Color(4) = 1.0;
    hold all;
    [f,xi] = ksdensity(top50Params(i,:) , linspace(lb(i) , ub(i) , 100));
    plot(xi,f);
    hold off;
    xlim([lb(i) ub(i)]);
    %xlim([min(top50Params(i,:)) max(top50Params(i,:))]);
    ylim([0.0 max(f)*1.25]);
    %ylim([0.0 max(N)*1.25]);
    %ylim([0.0 max(max(f)*1.25 , max(N)*1.25)]);
    title(titles{i});
    if  i == size(top50Params,1) %idx == size(tempIdx,2)
        %legend('paramSet-t0' , 'paramSet-t1' , 'paramSet-t2' , 'paramSet-t3' , 'paramSet-t4' , 'paramSet-t5' , 'paramSet-t6');
        %legend('t0' , 't1' , 't2' , 't3' , 't4' , 't5' , 't6' , 't7' , 't8' , ...
        %    't9' , 't10' , 't11' , 't12' , 't13' , 't14' , 't15' , 't16' , ...
        %    't17' , 't18' , 't19' , 't20' , 't21' , 't22');
        %legend('paramSetCurr' , 'paramSetNext' , 't0'); % 'prior' , 'paramSetNext-fixedWeights');
        legend(['t' , num2str(tstep_abc) , ' ksdensity']);
    end
end
