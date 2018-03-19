function [cd4Cost , artCost ,  ccHivTreatCost , ccArtTreatCost , ccTreatCost_hNeg] ...
    = costify(popVec , ccTreated , ccCost , hivTreatCost , artTreatCost)
%% Helper functions
annlz = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear)); % finds annual total
annAvg = @(x) sum(reshape(x , stepsPerYear , size(x , 1) / stepsPerYear))... % finds annual average
    / stepsPerYear; 
%% HIV-related costs
% Costs by CD4 count
cd4Cost = zeros(4 , size(popVec , 1));
for d = 1 : 4 % four CD4 categories in model
    cd4Cat = toInd(allcomb((d + 2) , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : gender , 1 : age , 1 : risk));
    cd4Cost(d , :) = sum(annAvg(popVec(: , cd4Cat)) .* hivTreatCost(d) , 2);
end

% ART treatment costs
artInd = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , ...
        1 : gender , 1 : age , 1 : risk));
artCost = sum(annAvg(popVec(: , artInd)) .* artTreatCost , 2);

%% Cervical cancer related costs
% Cervical cancer related costs by region (General)
ccTreatCost = zeros(3 , size(popVec , 1));
for s = 1 : 3
    for h = 2 : hpvTypes
    ccTreatCost(s) = ...
        sum(ccTreated(1 : disease , h , 1 : age , s) .* ccCost(s) , 2);
    end
end

% Cervical cancer related costs in HIV-positive not on ART
ccHivTreatCost = zeros(5 , 3 , size(popVec , 1));
for d = 2 : 6
    for s = 1 : 3
        ccHivTreatCost(d , s) = ...
            sum(ccTreated(d , 2 : hpvTypes , 1 : age , s) .* ccCost(s) , 2);
    end
end

% Cervical cancer related costs in HIV-positive on ART
ccArtTreatCost = zeros(3 , size(popVec , 1));
for s = 1 : 3
    ccArtTreatCost(s) = ...
        sum(ccTreated(10 , 2 : hpvTypes , 1 : age , s) .* artTreatCost , 2);
end

% Cervical cancer related costs in HIV-negative
ccTreatCost_hNeg = zeros(3 , size(popVec , 1));
for s = 1 : 3
    ccTreatCost_hNeg(s) = ...
        sum(ccTreated(1 , 2 : hpvTypes , 1 : age , s) .* artTreatCost , 2);
end