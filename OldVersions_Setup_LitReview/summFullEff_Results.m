function summFullEff_Results(year)
% read in data from model runs
%% Vaccination coverage
%CC Mortality
coverMort_Gen = readtable('General_CCMortality_VaxCover_standFull.csv');
coverMort_Neg = readtable('HIVNeg_CCMortality_VaxCover_standFull.csv');
coverMort_Art = readtable('VaxCover_ARTMort_standFull.csv');
coverMort_500 = readtable('VaxCover_CD4_500Mort_standFull.csv');
coverMort_500_350 = readtable('VaxCover_CD4_500_350Mort_standFull.csv');
coverMort_350_200 = readtable('VaxCover_CD4_350_200Mort_standFull.csv');
coverMort_200 = readtable('VaxCover_CD4_200Mort_standFull.csv');

% CC Incidence
coverInc_Gen = readtable('General_Incidence_standFull.csv');
coverInc_Neg = readtable('HIVNeg_Reduction_standFull.csv');
coverInc_Art = readtable('ART_Incidence_standFull.csv');
coverInc_500 = readtable('CD4_500_Incidence_standFull.csv');
coverInc_500_350 = readtable('CD4_500_350_Incidence_standFull.csv');
coverInc_350_200 = readtable('CD4_350_200_Incidence_standFull.csv');
coverInc_200 = readtable('CD4_200_Incidence_standFull.csv');

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