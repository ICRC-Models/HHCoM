initEff = 80;
kWane = - 80 / 20;
tVec = 1 : 30;
currEff(1 : 10) = 80 * ones(10 , 1);
currEff(11 : 30) = 80 + kWane * (1 : 20);
currEff2(1 : 15) = 80 * ones(15 , 1);
currEff2(16 : 30) = 80 + kWane * (1 : 15);
currEff3(1 : 20) = 80 * ones(20 , 1);
currEff3(21 : 30) = 80 + kWane * (1 : 10);
figure()
plot(tVec , currEff , tVec , currEff2 , tVec , currEff3)
ylim([0 100])
xlabel('Years After Vaccination'); ylabel('Vaccine Efficacy (%)')
title('Vaccine Waning')
legend('10 year waning' , '15 year waning' , '20 year waning')
