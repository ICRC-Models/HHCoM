   
tspan = 0 : 0.1 : 100;
artDist = zeros(prod([disease , viral , gender , age , risk]) , 1);
artDist = artDist + 100;
[~ , artTreat] = timeit(ode45(@(t , artDist) treatDist(t , pop1(end , :)) , tspan , artDist)); 

