function [dPop] = agePop(t , pop , ager)
% Aging module
% Ages population.
% 1/5th of the previous age group progresses into the next age group each year.
% Aged cohort is redistributed into new risk groups upon entry into new age
% group. Risk group status in previous age group does not affect risk group
% status in new age group. Risk redistribution is solely age dependent.
% Accepts a population vector as input and returns dPop, a vector of
% derivatives that describes the change in the population's subgroups due
% to aging.
if size(pop , 1) ~= size(ager , 2) 
    pop = pop';
end
dPop = ager * pop;