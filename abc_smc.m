% ABC-SMC
% Note: "particle" : a set of parameters 

function [] = abc_smc(t , alpha , p_acc_min , date)

t_prev = t-1;
t_curr = t;
t_next = t+1;

%% Load all particles 
paramDir = [pwd , '/Params/'];
paramSetMatrix = load([paramDir , 'paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent parameter sample
paramSetMatrix = paramSetMatrix(2:end,:);
pIdx = load([paramDir , 'pIdx_calib_' , date , '_0' , '.dat']); % load parameter indices
negSumLogLmatrix = load([paramDir , 'negSumLogL_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent log-likelihoods

%% Filter out failed parameter sets (timed-out, etc.)
numSubsets = size(negSumLogLmatrix,1)/17; % calculate number of sub-sets that actually ran (vs. timed-out, failed, etc.)
negS_format = reshape(negSumLogLmatrix , [17,numSubsets]); % first row is paramSetIdx, next 16 rows log-likelihoods for that sub-set
maxV = max(negS_format(1,:)); % find maximum paramSetIdx
setVec = [1:16:maxV];
missingV = [];
keepV = [];
for j = 1 : length(setVec) % identify failed parameter sets
     if ~any(setVec(j) == negS_format(1,:)) 
         missingV = [missingV , [setVec(j) : setVec(j)+15]];
     else
         keepV = [keepV , [setVec(j) : setVec(j)+15]];
     end
end
paramSetMatrix = paramSetMatrix(:,keepV); % filter failed parameter sets from matrix
numFltrdSets = length(keepV);

fileF = ['filteredSets_calib_' , date , '_' , num2str(t_curr) , '.dat'];
csvwrite([paramDir, fileF] , paramSetMatrix);

%% Calculate weights
fileW = ['weights_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save file of sorted log-likelihoods
if t_curr == 0
    % Set w to 1 for all particles in iteration 0
    weights = ones(numFltrdSets,1);
    csvwrite([paramDir, fileW] , weights);
elseif t_curr > 0
    % Load particles from iteration t-1
     paramSetMatrix_prev = load([paramDir , 'alphaParamSets_calib_' , date , '_' , num2str(t_prev) , '.dat']); % load most recent parameter sample
     weights_prev = load([paramDir , 'weights_calib_' , date , '_' , num2str(t_prev) , '.dat']);
     [weights] = calcWeights(paramSetMatrix , paramSetMatrix_prev , weights_prev' , pIdx);
end

%% Specify iteration and concatenate particle matrices
if t_curr == 0
    % Save initial set of particles
    masterSetMatrix = paramSetMatrix;
    fileM = ['masterSets_calib_' , date , '_' , num2str(t_curr) , '.dat'];
    csvwrite([paramDir, fileM] , masterSetMatrix);
elseif t_curr > 0
    % Append particles in iteration t-1 to particles in iteration t
    alphaSets_prev = load([paramDir , 'alphaParamSets_calib_' , date , '_' , num2str(t_prev) , '.dat']);
    masterSetMatrix = [alphaSets_prev , paramSetMatrix];
    fileM = ['masterSets_calib_' , date , '_' , num2str(t_curr) , '.dat'];
    csvwrite([paramDir, fileM] , masterSetMatrix);
    
    weights = [weights_prev ; weights];
end

%% Keep the top alpha-proportion of particles 
% Rho is the likelihood value of parameters given data. Want to maximize this value.
% Order log likelihood values in descending order, and keep top alpha-proportion.
[temp,firstRowOrder] = sort(negS_format(1,:)); % sort by paramSetIdx
negS_ordered = negS_format(:,firstRowOrder);

negS_ordered_flatDat = reshape(negS_ordered(2:end,:),[numFltrdSets,1]); % flatten into vector of log-likelihoods 
if t_curr == 0
    fileMLL = ['masterLL_calib_' , date , '_' , num2str(t_curr) , '.dat'];
    csvwrite([paramDir, fileMLL] , negS_ordered_flatDat);
elseif t_curr > 0
    negS_ordered_flatDat_prev = load([paramDir , 'masterLL_calib_' , date , '_' , num2str(t_prev) , '.dat']);
    negS_ordered_flatDat = [negS_ordered_flatDat_prev ; negS_ordered_flatDat];
end

[vals,inds] = sort(negS_ordered_flatDat,'descend'); % sort log-likelihoods in descending order

fileLL = ['orderedLL_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save file of sorted log-likelihoods
csvwrite([paramDir, fileLL] , [inds , vals]);

%% Save accepted particles from iteration t
alphaSets = masterSetMatrix(:,inds(1:(numFltrdSets*alpha)));
fileAlpha = ['alphaParamSets_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save file of top alpha-proportion of particles
csvwrite([paramDir, fileAlpha] , alphaSets);

%% Normalize weights
alphaWeights = weights(inds(1:(numFltrdSets*alpha)));
normWeights = alphaWeights/sum(alphaWeights);
csvwrite([paramDir, fileW] , normWeights);

%% Calculate epsilon (distance criterion)
eps = max(vals(1:(numFltrdSets*alpha)));
fileEps = ['epsilon_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save epsilon
csvwrite([paramDir, fileEps] , eps);

%% Calculate the proportion of accepted particles and compare to the minimum acceptance criterion
if t == 0
    p_acc = 1;
elseif t > 0
    p_acc = sumall(vals(1:(numFltrdSets*alpha)) < eps) / (numFltrdSets*alpha);
end

%% If the proportion of accepted particles is greater than the minimum acceptance criterion,
% Increase value of t and select new batch of particles to simulate
if p_acc > p_acc_min
    % Select particles for next iteration (t+1). Particles are sampled from iteration t with probability equal to their normalized weight.
    n_new_particles = (numFltrdSets - numFltrdSets*alpha);
    select_particles = datasample(alphaSets, round(n_new_particles) , 2 , 'Weights' , normWeights'); % uniform sample taken at random, with replacement
    
    % Perturn each of the selected particles.
    % Each particle is then moved with a Gaussian perturbation kernel with variance equal to twice the empirical weighted variance of the particles in the previous iteration
    new_particles = select_particles;
    new_particles = [zeros(1 , size(select_particles , 2)) ; select_particles]; % row of bools: all set to false
      
    while any(~new_particles(1,:))
        not_in_prior = find(~new_particles(1,:));
    
        for iS = 1 : length(not_in_prior)
            i = not_in_prior(iS);
            stddev = std(alphaSets,0,2); % standard deviation of sample by parameter (normalized by numPart-1)
            new_particles(2:end,i) = normrnd(select_particles(:,i), stddev);
            
            % Ensure that all perturbations are within the specified priors to ensure weights of subsequent particles > 0.
            [paramsAll] = genParamStruct();
            paramsSub = cell(length(pIdx),1);
            lb = [];
            ub = [];
            for s = 1 : length(pIdx)
                paramsSub{s} = paramsAll{pIdx(s)};
                lb = [lb; paramsSub{s}.lb];
                ub = [ub; paramsSub{s}.ub];
            end
            dBool = 1;
            for s = 1 : size(select_particles,1)
                dBool = dBool * unifpdf(new_particles(s+1,i) , lb(s) , ub(s));
            end
            new_particles(1,i) = (dBool > 0.0);  
        end
    end
    fileD = ['paramSets_calib_' , date , '_' , num2str(t_next) , '.dat'];
    csvwrite([paramDir, fileD] , new_particles(2:end,:));
end
