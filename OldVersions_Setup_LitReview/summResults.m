function summResults(year)
% read in data from model runs
%% Vaccination coverage
%CC Mortality
coverMort_Gen = readtable('General_CCMortality_VaxCover_stand.csv');
coverMort_Neg = readtable('HIVNeg_CCMortality_VaxCover_stand.csv');
coverMort_Art = readtable('VaxCover_ARTMort_stand.csv');
coverMort_500 = readtable('VaxCover_CD4_500Mort_stand.csv');
coverMort_500_350 = readtable('VaxCover_CD4_500_350Mort_stand.csv');
coverMort_350_200 = readtable('VaxCover_CD4_350_200Mort_stand.csv');
coverMort_200 = readtable('VaxCover_CD4_200Mort_stand.csv');

% CC Incidence
coverInc_Gen = readtable('General_Incidence_stand.csv');
coverInc_Neg = readtable('HIVNeg_Reduction_stand.csv');
coverInc_Art = readtable('ART_Incidence_stand.csv');
coverInc_500 = readtable('CD4_500_Incidence_stand.csv');
coverInc_500_350 = readtable('CD4_500_350_Incidence_stand.csv');
coverInc_350_200 = readtable('CD4_350_200_Incidence_stand.csv');
coverInc_200 = readtable('CD4_200_Incidence_stand.csv');

%% Vaccine waning
%CC Mortality
waneMort_Gen = readtable('waning_General_CCMortality_stand.csv');
waneMort_Neg = readtable('waning_hivNegMortality_stand.csv');
waneMort_Art = readtable('ARTMort_stand.csv');
waneMort_500 = readtable('CD4_500Mort_stand.csv');
waneMort_500_350 = readtable('CD4_500_350Mort_stand.csv');
waneMort_350_200 = readtable('CD4_350_200Mort_stand.csv');
waneMort_200 = readtable('CD4_200Mort_stand.csv');

% CC Incidence
waneInc_Gen = readtable('GenInc_Waning_stand.csv');
waneInc_Neg = readtable('HIVNeg_Waning_stand.csv');
waneInc_Art = readtable('ART_Waning_stand.csv');
waneInc_500 = readtable('CD4_500_Waning_stand.csv');
waneInc_500_350 = readtable('CD4_500_350_Waning_stand.csv');
waneInc_350_200 = readtable('CD4_350_200_Waning_stand.csv');
waneInc_200 = readtable('CD4_200_Waning_stand.csv');

%% Generate table
% pick year to analyze
stepsPerYear = 6;
yr = (year - 2018) * stepsPerYear;

%% Vaccine coverage scenarios
% Vaccine coverage: incidence reductions
disp('Relative reductions in CC incidence for vaccination coverage scenarios (90%/70%/50%)')
GenFem = flipud(coverInc_Gen{yr , 6 : 8})';
NegFem = flipud(coverInc_Neg{yr , 2 : 4})';
ArtFem = flipud(coverInc_Art{yr , 2 : 4})';
Fem500 = flipud(coverInc_500{yr , 2 : 4})';
Fem500_350 = flipud(coverInc_500_350{yr , 2 : 4})';
Fem_350_200 = flipud(coverInc_350_200{yr , 2 : 4})';
Fem_200 = flipud(coverInc_200{yr , 2 : 4})';

coverIncidence = table({'90' ; '70' ; '50'}, GenFem, NegFem, ArtFem, Fem500 , ...
    Fem500_350 , Fem_350_200 , Fem_200);
coverIncidence.Properties.VariableNames{1} = 'Coverage';
disp(coverIncidence)

% Vaccine coverage: mortality reductions
disp('Relative reductions in CC mortality for vaccination coverage scenarios (90%/70%/50%)')
GenFem = flipud(coverMort_Gen{yr , 5 : 7})';
NegFem = flipud(coverMort_Neg{yr , 5 : 7})';
ArtFem = flipud(coverMort_Art{yr , 5 : 7})';
Fem500 = flipud(coverMort_500{yr , 5 : 7})';
Fem500_350 = flipud(coverMort_500_350{yr , 5 : 7})';
Fem_350_200 = flipud(coverMort_350_200{yr , 5 : 7})';
Fem_200 = flipud(coverMort_200{yr , 5 : 7})';

coverMortality = table({'90' ; '70' ; '50'}, GenFem, NegFem, ArtFem, Fem500 , ...
    Fem500_350 , Fem_350_200 , Fem_200);
coverMortality.Properties.VariableNames{1} = 'Coverage';
disp(coverMortality)

%% Vaccine waning scenarios
% Vaccine coverage: incidence reductions
disp('Relative reductions in CC incidence for waning scenarios (20/15/10 years)')
GenFem = flipud(waneInc_Gen{yr , 6 : 8})';
NegFem = flipud(waneInc_Neg{yr , 6 : 8})';
ArtFem = flipud(waneInc_Art{yr , 6 : 8})';
Fem500 = flipud(waneInc_500{yr , 6 : 8})';
Fem500_350 = flipud(waneInc_500_350{yr , 6 : 8})';
Fem_350_200 = flipud(waneInc_350_200{yr , 6 : 8})';
Fem_200 = flipud(waneInc_200{yr , 6 : 8})';

waneIncidence = table({'20' ; '15' ; '10'}, GenFem, NegFem, ArtFem, Fem500 , ...
    Fem500_350 , Fem_350_200 , Fem_200);
waneIncidence.Properties.VariableNames{1} = 'WaningPeriod';
disp(waneIncidence)

% Vaccine coverage: mortality reductions
disp('Relative reductions in CC mortality for waning scenarios (20/15/10 years)')
GenFem = flipud(waneMort_Gen{yr , 6 : 8})';
NegFem = flipud(waneMort_Neg{yr , 6 : 8})';
ArtFem = flipud(waneMort_Art{yr , 6 : 8})';
Fem500 = flipud(waneMort_500{yr , 6 : 8})';
Fem500_350 = flipud(waneMort_500_350{yr , 6 : 8})';
Fem_350_200 = flipud(waneMort_350_200{yr , 6 : 8})';
Fem_200 = flipud(waneMort_200{yr , 6 : 8})';

waneMortality = table({'20' ; '15' ; '10'}, GenFem, NegFem, ArtFem, Fem500 , ...
    Fem500_350 , Fem_350_200 , Fem_200);
waneMortality.Properties.VariableNames{1} = 'WaningPeriod';
disp(waneMortality)