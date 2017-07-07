% Viral load progression
% Simulates viral load progression in HIV infected subpopulation.
% Accepts a population matrix as input and returns dPop, a matrix of
% derivatives that describes the change in the population's subgroups due
% to viral load progression.
function[dPop] = vlAdv(t , pop , vlAdvancer)
if size(pop , 1) ~= size(vlAdvancer , 2) 
    pop = pop';
end
dPop = vlAdvancer * pop;