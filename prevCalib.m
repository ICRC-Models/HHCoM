function objFun = prevCalib(toOptimize)
% choose whether to model hysterectomy
hyst = 'off';
hivOn = 1;
hpvOn = 0;
% choose whether to model HIV

disp(' ')
startYear = 1980;
endYear = 2013;
years = endYear - startYear;
% stepsPerYear = 4; % must reload data if changed
save('settings' , 'years' , 'startYear' , 'endYear')
%% results
% Load parameters and constants for main
load('general')
%% Initial population 
load('popData')
mInit = popInit(: , 1);
MsumInit = sum(mInit);

fInit = popInit(: , 2);
FsumInit = sum(fInit);

MpopStruc = riskDistM;
FpopStruc = riskDistF;

mPop = zeros(age , risk);
fPop = mPop;

for i = 1 : age
    mPop(i , :) = MpopStruc(i, :).* mInit(i) / 1.13;
    fPop(i , :) = FpopStruc(i, :).* fInit(i) / 1.13;
end

dim = [disease , viral , hpvTypes , hpvStates , periods , gender , age ,risk];
initPop = zeros(dim);
initPop(1 , 1 , 1 , 1 , 1 , 1 , : , :) = mPop;
initPop(1 , 1 , 1 , 1 , 1 , 2 , : , :) = fPop;

if hivOn
    toInfectM = (sum(mPop(:)) + sum(fPop(:))) * 0.001 * 0.5;
    toInfectF = (sum(mPop(:)) + sum(fPop(:))) * 0.001 * 0.5;
%     initPop(3 : 4 , 2 : 4 , 1 , 1 , 1 , 1 , 5 , 2) = 1; % initial HIV infected males (acute)
%     initPop(3 : 4 , 2 : 4 , 1 , 1 , 1 , 2 , 4 , 2) = 1; % initial HIV infected females (acute)
    initPop(3 , 2 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = toInfectM / 6; % initial HIV infected male (0.1% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) = ...
        initPop(1 , 1 , 1 , 1 , 1 , 1 , 4 : 6 , 2 : 3) - toInfectM / 6; % moved to HIV infected
    initPop(3 , 2 , 1 , 1 , 1 , 2 , 5 : 6 , 2 : 3) = toInfectF / 4; % initial HIV infected female (0.1% prevalence)
    initPop(1 , 1 , 1 , 1 , 1 , 2 , 5 : 6 , 2 : 3) = ...
        initPop(1 , 1 , 1 , 1 , 1 , 2 , 5 : 6 , 2 : 3) - toInfectF / 4; % moved to HIV infected
end

if hpvOn
    % initPop(1 , 1 , 2 : 4 , 1 , 1 , : , 4 : 9 , :) = 2; % initial HPV hr and lr infecteds (test)
    infected = initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) * 0.5;
    initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) = ...
            initPop(1 , 1 , 1 , 1 , 1 , : , 4 : 9 , :) - infected;
    for h = 2 : hpvTypes
        initPop(1 , 1 , h , 1 , 1 , : , 4 : 9 , :) = infected ./ 3;
    end
end

% Intervention start years
circStartYear = 1990;
vaxStartYear = 2010;

%% Simulation
% disp('Trial run started')
load('general');
load('settings');
load('mixInfectIndices')
load('vlAdvancer')
load('fertMat')
load('hivFertMat')
load('deathMat')
load('circMat')
load('vaxer')
load('mixInfectParams');
load('hpvData')
load('popData')
load('HIVParams')
load('hivIndices')
load('ager')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to be optimized
% epsA = eps(: , 1);
% epsR = eps(: , 2);
partnersM(5 : 9 , 2 : 3) = partnersM(5 : 9 , 2 : 3) .* toOptimize(: , : , 1);
partnersF(5 : 9 , 2 : 3) = partnersF(5 : 9 , 2 : 3) .* toOptimize(: , : , 2);
yr = [1998 ; 2003 ; 2010];
stepsPerYear = 6;
step = 1 / stepsPerYear;
epsA_vec = cell(size(yr , 1) - 1, 1); % save data over time interval in a cell array
epsR_vec = cell(size(yr , 1) - 1, 1); 
for i = 1 : size(yr , 1) - 1
    period = [yr(i) , yr(i + 1)];
    epsA_vec{i} = interp1(period , epsA(i : i + 1 , 1) , ...
        yr(i) : step : yr(i + 1));
    epsR_vec{i} = interp1(period , epsR(i : i + 1 , 1) , ...
        yr(i) : step : yr(i + 1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hivPositiveArtAll = sort(toInd(allcomb(10 , 6 , 1 : hpvStates , 1 : hpvTypes , ...
    1 : periods , 1 : gender  , 1 : age , 1 : risk)));
% profile on
% disp(' ')
% Initialize vectors
timeStep = 1 / stepsPerYear;
years = endYear - startYear;
s = 1 : timeStep : years + 1; % stepSize and steps calculated in loadUp.m
artDistMat = zeros(size(prod(dim) , 20)); % initialize artDistMat to track artDist over past 20 time steps
%performance tracking
runtimes = zeros(size(s , 2) - 2 , 1);
import java.util.LinkedList
artDistList = LinkedList();
popVec = spalloc(years / timeStep , prod(dim) , 10 ^ 8);
popIn = reshape(initPop , prod(dim) , 1); % initial population to "seed" model
newHiv = zeros(length(s) - 1 , gender , age , risk);
newHpv = newHiv;
hivDeaths = zeros(length(s) - 1 , age);
deaths = popVec;
artTreatTracker = zeros(length(s) - 1 , disease , viral , gender , age , risk);
popVec(1 , :) = popIn;
tVec = linspace(startYear , endYear , size(popVec , 1));
k = cumprod([disease , viral , hpvTypes , hpvStates , periods , gender , age]);
artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0
% disp(['Simulating period from ' num2str(startYear) ' to ' num2str(endYear) ...
%     ' with ' num2str(stepsPerYear), ' steps per year.'])
% disp(' ')
% disp('Simulation running...')
% disp(' ')
% Equations for each module in simulation
% Evaluate diff eqs over one time interval. Use output from preceding
% differential equations as input to successive differential equations.
% % Incident infection plot that updates at each time step
% figure()
% xlabel('Year'); ylabel('New cases (%)'); title('Relative HIV Incidence')
try
%     progressbar('Simulation Progress')
    for i = 2 : length(s) - 1
        tic
        year = startYear + s(i) - 1;
        currStep = round(s(i) * stepsPerYear);
%         disp(['current step = ' num2str(startYear + s(i) - 1) ' ('...
%             num2str(length(s) - i) ' time steps remaining until year ' ...
%             num2str(endYear) ')'])
        tspan = [s(i) , s(i + 1)]; % evaluate diff eqs over one time interval
        [~ , pop , newHiv(i , : , : , :) , newHpv(i , : , : , :)] = ...
            ode4xtra(@(t , pop) mixInfect(t , pop , 'vl' , ...
            currStep , naive , coInf, gar , hivSus , hpvSus , toHiv , ...
            toHpv , mCurr , fCurr , mCurrArt , fCurrArt ,epsA_vec , ...
            epsR_vec , yr , modelYr1 ,  circProtect , condProtect , condUse ,...
            actsPer , partnersM , partnersF , beta_hrHPV_val , beta_lrHPV_val , ...
            disease , viral , gender , age , risk , hpvStates , hpvTypes , ...
            k , periods , stepsPerYear , year) , tspan , popIn); 
        if any(pop(end , :) < 0)
            disp('After mixInfect')
            break
        end
        if hivOn
            [~ , pop , hivDeaths(i , :) , artTreat] =...
                ode4xtra(@(t , pop) hiv(t , pop , vlAdvancer , artDist , muHIV , hivPositiveArtAll , ...
                kCD4 , disease , viral , gender , age , risk , k , artOut , hivInds , ...
                stepsPerYear , below200Art_2004 , above200Art_2004 , pie4Vec_2004 , ...
                below200Art_2006 , year) , tspan , pop(end , :));
            artTreatTracker(i , : , : , : , :  ,:) = artTreat;
            if any(pop(end , :) < 0)
                disp('After hiv')
                break
            end
    %         [~ , artTreat] = ode4x(@(t , artDist) treatDist(t , popCopy(end , :) , year) , tspan , artDist); 
            if size(artDistList) >= 20
                artDistList.remove(); % remove earlier artDist matrix more than "20 time steps old"
            else 
                artDistList.add(artTreat); 
            end
            artDist = calcDist(artDistList);
        end
        if hpvOn
            [~ , pop] = ode4x(@(t , pop) hpv(t , pop , hyst) , ...
                tspan , pop(end , :));
            [~ , pop] = ode4x(@cinAdv , tspan , pop(end , :));
            [~ , pop] = ode4x(@(t , pop) hpvTreat(t , pop , hyst , year)...
                , tspan , pop(end , :)); 
        end
        [~ , pop , deaths(i , :)] = ode4xtra(@(t , pop) bornAgeDie(t , pop , ager , year , age , currStep , fertility2k , fertMat , hivFertMat ,...
        deathMat , circMat , vaxer , mtctVec , circStartYear , vaxStartYear , modelYr1 , stepsPerYear) , tspan , pop(end , :));
        if any(pop(end , :) < 0)
            disp('After bornAgeDie')
            break
        end
        % add results to population vector
        popVec(i , :) = pop(end , :);
        popIn = popVec(i , :);
        runtimes(i - 1) = toc;
%     progressbar(i/(length(s) - 1))
    end

%     disp(['Reached year ' num2str(endYear)])
    popVec = sparse(popVec); % compress population vectors
%     save('results' , 'tVec' ,  'popVec' , 'newHiv' , 'newHpv' , 'hivDeaths' , ...
%         'deaths'); 
%     disp(' ')
%     disp('Trial run complete.')
catch
    disp('Trial run aborted.')
    disp(['Error after ' num2str(startYear + s(i) - 1) ' ('...
            num2str(length(s) - i) ' time steps remaining until year ' ...
            num2str(endYear) ')'])
    save('errorWorkspace')
    disp('All variables saved to errorWorkspace for debugging.')
    return
end

% compare output to actual results
hivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : gender , 1 : age , 1 : risk));
hivOut = sum(popVec(: , hivInds) , 2);
artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates, ...
    1 : periods , 1 : gender , 1 : age , 1 : risk));
art = sum(popVec(: , artInds) , 2);

modelHivPos = (hivOut + art) ./ sum(popVec , 2);
yrs = [1990 : 2010 , 2012];
[~ , x , ~] = intersect(round(tVec) , yrs);
dataModelHivPos = modelHivPos(x);
load('actual')
% Evaluate objective function
weights = 2017 - fliplr(yrs)';
objFun = sum((100 * (dataModelHivPos - [overallHivPrev_SA(: , 2) ; overallHivPrev_KZN(: , 2)]) ...
    .* (weights ./ max(weights))) .^ 2);
