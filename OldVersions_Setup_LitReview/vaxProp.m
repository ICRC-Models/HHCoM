%% HPV Prevalence
fac = 10 ^ 5;
v90_Vaxd =  zeros(age , length(tVec));
v0_Vaxd = v90_Vaxd;
for a = 3 : age
    % all
    % general
    allFTot = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , ...
        1 : periods , 2 , a , 1 : risk));
    % All vaccinated
    allVax = [toInd(allcomb(1 : disease , 1 : viral , 1 : 4 , 1 , ...
        2 : periods , 2 , a , 1 : risk)); toInd(allcomb(1 : disease , 1 : viral , 1 , 9 , ...
        1 : periods, 2 , a , 1 : risk))];
    
    v90_Vaxd(a , :) = ...
        sum(sum(sum(sum(sum(sum(sum(sum(o90.popVec(: , allVax) ,2),3),4),5),6),7),8),9) ./ ...
        sum(sum(sum(sum(sum(sum(sum(sum(o90.popVec(: , allFTot) ,2),3),4),5),6),7),8),9);
    
    v0_Vaxd( a, :) = ...
        sum(sum(sum(sum(sum(sum(sum(sum(oNo.popVec(: , allVax) ,2),3),4),5),6),7),8),9) ./ ...
        sum(sum(sum(sum(sum(sum(sum(sum(oNo.popVec(: , allFTot) ,2),3),4),5),6),7),8),9);    
end
figure()
ageGroup = {'10 - 14' , '15 - 19' , '20 - 24' , '25 - 29' ,...
    '30 - 34' , '35 - 39' , '40 - 44' , '45 - 49' , '50 - 54' , '55 - 59' , ...
    '60 - 64' , '65 - 69' , '70 - 74' , '75 - 79'};

m1 = mesh(1 : age , tVec , v0_Vaxd' * 100);
set(m1 , 'edgecolor' , 'r')
alpha(0.1)
hold on
m2 = mesh(1 : age , tVec , v90_Vaxd' * 100);
set(m2 , 'edgecolor' , 'b')
alpha(0.1)
set(gca , 'yLim' , [tVec(1) tVec(end)]);
set(gca , 'xtick' , 3 : age , 'xtickLabel' , ageGroup);
xlabel('Age Group'); ylabel('Year'); zlabel('Prevalence (%)')
title('Vaccination')
legend('No Vaccination' , '90% coverage')