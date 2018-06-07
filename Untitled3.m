tspan = [2018 : 1/6 : 2025];
load([paramDir , 'ageRiskInds'])
riskDist(: , : , 1) = riskDistM;
riskDist(: , : , 2) = riskDistF;
profile on
[~ , pop2 , ~] = ode4xtra(@(t , pop) ...
    ageRisk(t , pop , disease , hpvTypes , viral , hpvStates , gender , ...
    periods , age , risk ,riskDist , ageInd , riskInd) , tspan , pop(end , :));
if any(pop2(end , :) < 0)
    disp('Negative values in compartments')
    return
end
profile viewer
% Population Size
figure()
plot(tspan , sum(pop2 , 2))
title('Population Size'); xlabel('Year'); ylabel('Population size')

%% risk distribution
gen = 'Male';
for g = 1 : gender
    if g == 2
        gen = 'Female';
    end
    figure()
    for a = 1 : age
        aGroup = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
            1 : hpvStates , 1 : periods , 1 : gender , a , 1 : risk));
        subplot(4 , 4 , a)
        for r = 1 : risk
            rGroup = toInd(allcomb(1 : disease , 1 : viral ,  1 : hpvTypes ,...
                1 : hpvStates , 1 : periods , 1 : gender , a , r));
            popRisk = sum(popVec(: , rGroup) , 2) ./ sum(popVec(: , aGroup) , 2);
            plot(tVec , popRisk)
            title([gen , ' Age Group ' , num2str(a)])
            hold on
        end
    end
    legend('Low Risk' , 'Medium Risk' , 'High Risk')
end
