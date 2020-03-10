% Identify parameter ranges of top fifty best-fitting sets and compare to bounds
function [] = idParamInPrior(tstep_abc , date_abc , paramSetInd)
t_curr = tstep_abc;
date = date_abc;

%% Load all particles
paramDir = [pwd , '/Params/'];
paramSetMatrix = load([paramDir , 'paramSets_calib_' , date , '_' , num2str(t_curr) , '_7867mod17' , '.dat']); % load most recent parameter sample
pIdx = load([paramDir , 'pIdx_calib_' , date , '_0' , '.dat']); % load parameter indices

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

lbVals = zeros(size(paramSetMatrix,1),1);
ubVals = lbVals;
for i = 1 : size(paramSetMatrix,1)
    lbVals(i,1) = lb(i);
    ubVals(i,1) = ub(i);
end

%% Plot intermediate parameter posterior distributions
figure;
titles = {'partnersM:10-14,hr' , 'partnersM:15-19,hr' , 'partnersM:20-24,hr' , ...    % 24Feb20 calibration
         'partnersM:25-29,hr' , 'partnersM:30-44,hr' , 'partnersM:45-79,hr' , ...
         'partnersM:10-14,mr' , 'partnersM:10-14,lr' , ...
         'partnersF:10-14,hr' , 'partnersF:15-19,hr' , 'partnersF:20-24,hr' , ...
         'partnersF:25-29,hr' , 'partnersF:30-44,hr' , 'partnersF:45-79,hr' , ...
         'partnersF:10-14,mr' , 'partnersF:10-14,lr' , ...
         'condUse' , 'epsA' , 'perPartnerHpv' , 'lambdaMultImm' , ...
         'Inf->CIN1, vax' , 'Inf->CIN1, nv' , 'CIN1->CIN2, vax' , 'CIN1->CIN2, nv' , ...
         'CIN2->CIN3, vax' , 'CIN2->CIN3, nv' , 'CIN3->CC, vax' , 'CIN3->CC, nv' , ...
         'Inf->Well, vax' , 'Inf->Well, nv' , 'CIN1->Inf, vax' , 'CIN1->Inf, nv' , ...
         'CIN2->CIN1, vax' , 'CIN2->CIN1, nv' , 'CIN3->CIN2, vax' , 'CIN3->CIN2, nv' , ...
         'HIV transmission' , 'fert decline 1' , 'fert decline 2' , 'maleHpvClearMult'};
for i = 1 : size(paramSetMatrix,1)
    subplot(5,8,i);
    %prior = unifpdf([lb(i):((ub(i)-lb(i))/size(paramSetMatrix,2)):ub(i)] , lb(i) , ub(i));
    %plot([lb(i):((ub(i)-lb(i))/length(top50Inds)):ub(i)] , 'k');
    %hold on;
    plot(paramSetMatrix(i,paramSetInd) , 0.0 , 'k*');
    hold on;
    text(paramSetMatrix(i,paramSetInd) , 0.1 , num2str(paramSetMatrix(i,paramSetInd)));
    hold off;
    xlim([lb(i) ub(i)]);
    xticks([lb(i) ub(i)]); %((ub(i)-lb(i))/2)
    xticklabels({num2str(lb(i)) num2str(ub(i))}); %num2str((ub(i)-lb(i))/2)
    ylim([0.0 0.3]);
    set(gca,'YTick', [])
    text;
    title(titles{i});
end
