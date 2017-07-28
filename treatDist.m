% Calculates ART treatment distribution
% Accepts:
% 1) Population matrix (pop)
% Returns:
% artTreat, a matrix describing the distribution of individuals who went
% on ART according to their disease and viral load status at the time they
% went on treatment.
function[artTreat] = treatDist(t , pop , year)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load constants and parameters
load('general')
load('settings')
load('HIVParams')
load('hivIndices')
%%%
% pie [d x v x g x a x r]
% values obtained from k4_mat in MainMain.m of HIV model.
pie = zeros(disease , viral , gender , age ,risk);
pie(1 , : , : , 4 : end , :) = 0; % rate of going on PrEP (0.92)
pie(2 , : , : , 4 : end , :) = 0; % ART coverage
pie(3 , : , : , 4 : end , :) = 0; % ART coverage
pie(4 , : , : , 4 : end , :) = 0; % ART coverage
pie(5 , : , : , 4 : end , :) = 0; % ART coverage
pie(6 , : , : , 4 : end, :) = 0; % ART coverage
if year >= 2014
    pie(6 , : , : , 4 : end , :) = 0.06; % Home HTC study
    pie(5 , : , : , 4 : end , :) = 0.02;
    pie(4 , : , : , 4 : end , :) = 0.01;
elseif year >= 2006
    pie(6 , : , : , 4 : end, :) = 0.03; % ART coverage for persons with CD4 < 200
end

% see model notes for index values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
artTreat = zeros(disease , viral , gender , age , risk);
% hsp = ones(prod([hpvTypes , hpvStates , periods]) , 1)';
% sp = ones(hpvStates * periods , 1)';
% hp = ones(hpvTypes * periods, 1)';
% hs = ones(hpvTypes * hpvStates  , 1)';
% allTypes = kron(sp , 1 : hpvTypes);
% allStates = kron(hp , 1 : hpvStates);
% allPeriods = kron(1 : periods , hs);

for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            for v = 1 : viral
%                 acuteInf = toInd({2 * hsp , v * hsp , allTypes , ...
%                     allStates , allPeriods , g * hsp , a * hsp , r * hsp});
%                 acuteInf = toInd(hivInds(2 , v , g , a , r)); % allcomb(2 , v , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r));
%                 artTreat(2 , v , g , a , r) = pie(2 , v , g , a , r) ... % keep track of distribution of people going on ART
%                     .* sumall(pop(acuteInf)); % going on ART
                % CD4 > 500 cells/uL (d = 3)
                % CD4 500-350 cells/uL (d = 4)
                % CD4 350-200 cells/uL (d = 5)
                % CD4 <200 cells/uL (d = 6)
                for d = 3 : 6
%                     cd4Curr = toInd({d * hsp , v * hsp , ...
%                         allTypes , allStates , allPeriods , g * hsp , a * hsp , ...
%                         r * hsp});
                    cd4Curr = toInd(hivInds(d , v , g , a , r)); % allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r));
                    artTreat(d , v , g , a , r) = pie(d , v , g , a , r) ... % keep track of distribution of people going on ART
                        .* sumall(pop(cd4Curr)); % going on ART
                end
            end
        end
    end
end

% Convert to column vector for output to ODE solver
artTreat = reshape(artTreat , [numel(artTreat) , 1]);