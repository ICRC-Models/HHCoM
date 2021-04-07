
COMPARTMENTAL MODEL STRUCTURE:
_________________________________________________________________________
COMPARTMENTS:  pop(d , v , h , s , x , p , g , a , r)
d - (disease) HIV disease state
v - (viral) Viral load
h - (hpvVaxStates) Vaccine-type HPV precancer state
s - (hpvNonVaxStates) Non-vaccine-type HPV precancer state
x - (endpoints) Cervical cancer or hysterectomy status
p - (intervens) Vaccination and screening history
g - (gender) Gender
a - (age) Age
r - (risk) Risk

For d:
1 - HIV negative, uncircumcised
2 - HIV negative, circumcised
3 - Acute infection
4 - CD4 > 500 cells/uL
5 - CD4 500-350 cells/uL
6 - CD4 350-200 cells/uL
7 - CD4 <= 200 cells/uL
8 - HIV-positive on ART

For v:
1 - Acute infection (if 2<d<8) , 0 (if d = 1/2)
2 - Asympotomatic (3-4.5 log10 VL)
3 - Pre-AIDS symptomatic (4-5.5 log10 VL)
4 - AIDS (5.5-7 log10 VL)
5 - Late-stage 
6 - On ART, virally suppressed

For h:
1 - Susceptible to vaccine type hrHPV 
2 - Infected with vaccine type hrHPV (16/18 & other 9v types)
3 - CIN1 (vaccine type)
4 - CIN2 (vaccine type)
5 - CIN3 (vaccine type)
6 - Cervical Cancer or hysterectomy
7 - Immune (vaccine type)

For s:
1 - Susceptible to non-vaccine type hrHPV
2 - Infected with non-vaccine type hrHPV (other high risk)
3 - CIN1 (non-vaccine type)
4 - CIN2 (non-vaccine type)
5 - CIN3 (non-vaccine type)
6 - Cervical Cancer or hysterectomy
7 - Immune (non-vaccine type)

For x:
1 - Cervical Cancer (Local) (if h or s = 6)
2 - Cervical Cancer (Regional)
3 - Cervical Cancer (Distant)
4 - Hysterectomy

For p:
1 - Non-vaccinated, non-screened
2 - Vaccinated
3 - Screened
4 - Vaccinated, screened

For g:
1 - Male
2 - Female

For a:
1- 0 (0-4)
2- 1 (5-9)
3- 2 (10-14)
4- 3 (15-19)
5- 4 (20-24)
6- 5 (25-29)
7- 6 (30-34)
8- 7 (35-39)
9- 8 (40-44)
10- 9 (45-49)
...
14- (65-69)
15- (70-74)
16- (75-79)

For r:
1 - Low risk
2 - Medium risk
3 - High risk
*risk is defined by the number of partners per year




NOTES: 
____________________________________________________________________

How to retrieve values for each compartment
toInd(allcomb(d , v , h , s , x , p , g , a , r))
Each of the above variables can be a single number, e.g. "3", or a range of numbers, e.g. "2:5"



SCRIPTS AND FUNCTIONS USED TO RUN THE MODEL:
_________________________________________________________________________

Run model:
1) historicalSim.m 
   - calculates historical trends up to current year
   - saves outputs to the \HHCoM_Results directory:
     \HHCoM_Results\toNow... .mat
   - uses the following functions:
     loadUp2.m (loads parameters)
     ode4xtra.m (numerical solver)
     hpvCCNH.m (HPV natural history)
     hpvScreen.m (Cancer screening and treatment)
     mixInfect.m (HPV and HIV transmission)
     hivNH.m (HIV natural history)
     calcDist.m (calculate distribution of persons initiating ART)
     artMinMax.m (calculate ART initiation/removal to maintain min/max coverage by age) 
     artPopCov.m (calculate ART initiation/removal to maintain population level coverage)
     bornAgeDieRisk.m (demography)
     vmmc.m (voluntary male medical circumcision)
     hpvVaxSchool.m (HPV vaccination - school-based regimen)
     likeFun.m (calculate summed log-likelihood)
   - visualize results using:
     showResults.m

2) futureSim.m 
   - calculates future predictions
   - loads historical results from the \HHCoM_Results directory for initial conditions:
     \HHCoM_Results\toNow... .mat
   - saves outputs to the \HHCoM_Results\Vaccine... directory:
     \HHCoM_Results\Vaccine...\vaxSimResult or \vaxWaneSimResult
   - uses the following functions:
     loadUp2.m (loads parameters)
     ode4xtra.m (numerical solver)
     hpvCCNH.m (HPV natural history)
     hpvScreen.m (Cancer screening and treatment)
     mixInfect.m (HPV and HIV transmission)
     hivNH.m (HIV natural history)
     calcDist.m (calculate distribution of persons initiating ART)
     artMinMax.m (calculate ART initiation/removal to maintain min/max coverage by age) 
     artPopCov.m (calculate ART initiation/removal to maintain population level coverage)
     bornAgeDieRisk.m (demography)
     vmmc.m (voluntary male medical circumcision)
     hpvVaxLmtd.m (HPV vaccination - vaccine limited years)
     hpvVaxSchool.m (HPV vaccination - school-based regimen)
     hpvVaxCU.m (HPV vaccination - catch-up regimen)
     parsave.m (save future results)
   - visualize results using:
     vaxCEA.m

