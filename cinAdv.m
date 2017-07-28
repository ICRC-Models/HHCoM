% Update duration of CIN status
% Restructures population according to how long certain subgroups have had
% CIN2 or CIN3 status. Advances 1/20th (1 / (period size * steps per year)) of individuals in a given CIN2 or
% CIN3 disease duration group into the subsequent CIN duration group.
% "Remaining fraction" of individuals in given disease duration group "move
% within" current duration group.
function[dPop] = cinAdv(t , pop)

% load general model parameters/values
load('general')
load('settings')
load('cinAdvIndices')
load('cinAdvData')
load('hpvData')
%initialize vector
dPop = zeros(size(pop));

% Advance 1st infection duration group
% inf1 , kAdv loaded from 'cinAdvIndices'

% dPop(inf1) = - pop(inf1) .* kAdv(1); % leaving 1st group
% 
% % Advance into 2 - 8
% 
% for p = 2 : periods - 1
%     % get indices
%     infCurr = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : hpvStates , ...
%         p , 1 : gender , 1 : age , 1 : risk));
%     infPrev = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : hpvStates , ...
%         (p - 1) , 1 : gender , 1 : age , 1 : risk));
%     % calculate
%     dPop(infCurr) = (pop(infPrev) - pop(infCurr)) .* kAdv(p); 
% end
% 
% % Advance into last period
% % get indices
% infLast = toInd(allcomb(1 : disease , 1  : viral , 2 : 4 , 1 : hpvStates ,...
%     periods , 1 : gender , 1 : age , 1 : risk));
% infPrev = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : hpvStates ,...
%     (periods - 1) , 1 : gender , 1 : age , 1 : risk));
% % calculate
% dPop(infLast) = pop(infPrev) .* kAdv(p); 

% CC region progression
dPop(local) = dPop(local) - pCC(1) * pop(local); % local -> regional
dPop(regional) = dPop(regional) + pCC(1) * pop(local) - pCC(2) * pop(regional); % local -> regional -> distant
dPop(distant) = dPop(distant) + pCC(2) * pop(regional); % regional -> distant

% CC duration progression
% for p = 1 : 2
%     cc = toInd(ccInds(p));
%     dPop(cc) = dPop(cc) - kCC(p) * pop(cc);
% end

% CC region and duration progression

dPop = sparse(dPop);