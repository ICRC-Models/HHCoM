% ABC-SMC
% Note: "particle" : a set of parameters
% Maximizing the positive summed Log-Likelihood even thought variable naming is old and suggests minimizing the negative summed-LL 

function [] = abc_smc(tstep_abc , date_abc)  %(alpha , p_acc_min)
t = tstep_abc;
alpha = 0.6;
p_acc_min = 0.05;
date = date_abc;

t_prev = t-1;
t_curr = t;
t_next = t+1;

%% Load all particles 
paramDir = [pwd , '/Params/'];
paramSetMatrix = load([paramDir , 'paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent parameter sample
pIdx = load([paramDir , 'pIdx_calib_' , date , '_0' , '.dat']); % load parameter indices
negSumLogLmatrix = load([paramDir , 'negSumLogL_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent log-likelihoods

%% Filter out failed parameter sets (timed-out, etc.)
numSubsets = size(negSumLogLmatrix,1)/17; % calculate number of sub-sets that actually ran (vs. timed-out, failed, etc.)
negS_format = reshape(negSumLogLmatrix , [17,numSubsets]); % first row is paramSetIdx, next 16 rows log-likelihoods for that sub-set

[uniqueN uInds v] = unique(negS_format(1,:) , 'first'); % remove duplicate sub-set runs
negS_format = negS_format(:,uInds);
numSubsets = size(negS_format,2); % reset as number unique sub-sets

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
paramSetMatrix = paramSetMatrix(:,keepV); % filter and save only successful parameter sets from matrix
numFltrdSets = length(keepV); % number sets that successfully ran

fileF = ['filteredSets_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save file of successfully run parameter sets
csvwrite([paramDir, fileF] , paramSetMatrix);

%% Calculate weights
if t_curr == 0
    % Set w to 1 for all particles in iteration 0
    masterWeights = ones(numFltrdSets,1);
elseif t_curr > 0
    % Load particles from iteration t-1
     alphaSets_prev = load([paramDir , 'alphaParamSets_calib_' , date , '_' , num2str(t_prev) , '.dat']); % load most recent parameter sample
     alphaWeights_prev = load([paramDir , 'alphaWeights_calib_' , date , '_' , num2str(t_prev) , '.dat']);
     [weights] = calcWeights(paramSetMatrix , alphaSets_prev , alphaWeights_prev' , pIdx);
     masterWeights = [alphaWeights_prev ; weights];
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
    master_negS_ordered_flatDat = negS_ordered_flatDat;
elseif t_curr > 0
    negS_ordered_flatDat_prev = load([paramDir , 'alphaLL_calib_' , date , '_' , num2str(t_prev) , '.dat']);
    master_negS_ordered_flatDat = [negS_ordered_flatDat_prev ; negS_ordered_flatDat];
end
masterNumFltrdSets = length(master_negS_ordered_flatDat);

[vals,inds] = sort(master_negS_ordered_flatDat,'descend'); % sort log-likelihoods in descending order

fileLL = ['orderedLL_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save file of sorted log-likelihoods
csvwrite([paramDir, fileLL] , [inds , vals]);

%% Save accepted particles from iteration t and their associated data
alphaSets = masterSetMatrix(:,inds(1:(masterNumFltrdSets*alpha)));
fileAlpha = ['alphaParamSets_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save file of top alpha-proportion of particles
csvwrite([paramDir, fileAlpha] , alphaSets);

alphaNegS_ordered = master_negS_ordered_flatDat(inds(1:(masterNumFltrdSets*alpha)));
fileALL = ['alphaLL_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save sorted negLogLL for alpha accepted particles
csvwrite([paramDir, fileALL] , alphaNegS_ordered);

alphaWeights = masterWeights(inds(1:(masterNumFltrdSets*alpha)));
%% Normalize weights
normWeights = alphaWeights./sum(alphaWeights);
fileW = ['alphaWeights_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save file of weights of top alpha-proportion of particles
csvwrite([paramDir, fileW] , normWeights);

%% Calculate epsilon (distance criterion)
eps = min(alphaNegS_ordered); % set epsilon as the largest distance/ lowest LL of the alpha-proportion of accelpted particles (which have
                              % the smallest distances/ largest LL between simulated and observed data)
fileEps = ['epsilon_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save epsilon
csvwrite([paramDir, fileEps] , eps);

%% Calculate the proportion of accepted particles and compare to the minimum acceptance criterion
if t == 0
    p_acc = 1;
elseif t > 0
    eps_prev = load([paramDir , 'epsilon_calib_' , date , '_' , num2str(t_prev) , '.dat']);
    p_acc = sum(negS_ordered_flatDat > eps_prev) / (numFltrdSets);
end
filePacc = ['p_acc_calib_' , date , '_' , num2str(t_curr) , '.dat']; % save p_acc 
csvwrite([paramDir, filePacc] , p_acc);

%% If the proportion of accepted particles is greater than the minimum acceptance criterion,
% Increase value of t and select new batch of particles to simulate
if p_acc > p_acc_min
    % Select particles for next iteration (t+1). Particles are sampled from iteration t with probability equal to their normalized weight.
    n_new_particles = (masterNumFltrdSets * (1 - alpha));
    select_particles = datasample(alphaSets, round(n_new_particles) , 2 , 'Weights' , normWeights'); % uniform sample taken at random, with replacement
    
    % Perturb each of the selected particles.
    % Each particle is then moved with a Gaussian perturbation kernel with variance equal to twice the empirical weighted variance of the particles in the previous iteration
    new_particles = select_particles;
    new_particles = [zeros(1 , size(select_particles , 2)) ; select_particles]; % row of bools: all set to false
      
    % Load bounds for comparison
    pIdx_length = length(pIdx);
    [paramsAll] = genParamStruct();
    paramsSub = cell(length(pIdx),1);
    lb = [];
    ub = [];
    for s = 1 : pIdx_length
        paramsSub{s} = paramsAll{pIdx(s)};
        lb = [lb; paramsSub{s}.lb];
        ub = [ub; paramsSub{s}.ub];
    end
 
    while any(~new_particles(1,:)) % while any particles have any parameter values outside bounds
        not_in_prior = find(~new_particles(1,:)); % indices for particles not in prior
        
        select_particles_length = size(select_particles , 1); % number parameters in particle
        notInPrior = zeros((select_particles_length + 1) , length(not_in_prior)); % matrix to hold new particles
   
        parfor iS = 1 : length(not_in_prior)
            i = not_in_prior(iS);
            v = 2 .* var(alphaSets,0,2); % 2x variance of previous accepted particles (normalized by numPart-1)
            stddev = v .^ (1/2);
            %stddev = std(alphaSets,0,2); % standard deviation of sample by parameter (normalized by numPart-1)
            
            new_particlesT = zeros(select_particles_length + 1 , 1);
            new_particlesT(2:end) = normrnd(select_particles(:,i) , stddev);
            
            % Ensure that all perturbations are within the specified priors to ensure weights of subsequent particles > 0.
            dBoolTot = 1;
            for sP = 1 : select_particles_length
                sPT = sP + 1;
                dBool = unifpdf(new_particlesT(sPT) , lb(sP) , ub(sP));
                while (dBool <= 0.0)
                    new_particlesT(sPT) = normrnd(select_particles(sP,i) , stddev(sP)); 
                    dBool = unifpdf(new_particlesT(sPT) , lb(sP) , ub(sP));
                end
                dBoolTot = dBoolTot * dBool;
            end
            new_particlesT(1) = (dBoolTot > 0.0);
            notInPrior(:,iS) = new_particlesT;  
        end
        new_particles(:,not_in_prior) = notInPrior;
    end
    fileD = ['paramSets_calib_' , date , '_' , num2str(t_next) , '.dat'];
    csvwrite([paramDir, fileD] , new_particles(2:end,:));
end
