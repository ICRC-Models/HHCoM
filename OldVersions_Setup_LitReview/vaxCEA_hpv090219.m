file = [pwd , '\HHCoM_Results\Vaccine082119_WHOP1_SCE1\P1_SCE1Coverage80_AgeStandHpvPrev_ages9-74.xlsx'];
p1S1 = xlsread(file , 'General' , 'B1:B101');

fileB = [pwd , '\HHCoM_Results\Vaccine082119_WHOP1_SCE1\P1_SCE1Coverage0_AgeStandHpvPrev_ages9-74.xlsx'];
p1S0 = xlsread(fileB , 'General' , 'B1:B101');
years = xlsread(fileB , 'General' , 'A1:A101');
hpvReduct = (p1S1 - p1S0) ./ p1S0 * 100;
plot(years , hpvReduct)
grid on

xlim([2020 2120]);
ylim([-100 0]);
xticks([2020 : 20 : 2120]);
xlabel('Year')
ylabel('Percent Change')
