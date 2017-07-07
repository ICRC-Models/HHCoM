load('general')
syms d v h s p g a r
%% general


%% HIV
% on ART
artInds = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , 1: gender , ...
    1 : age , 1 : risk)); 
%artInds = matlabFunction(artInds));

% CD4 < 500
cd4_500Inds = toInd(allcomb(4 : 6 , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%cd4_500Inds = matlabFunction(cd4_500Inds));

% Test and Treat

% All HIV+
hivInds = toInd(allcomb(2 : 6 , 1: viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , ...
    1: gender , 1 : age , 1 : risk));
%hivInds = matlabFunction(hivInds));

%% HPV

% All HPV+
hpvInds = toInd(allcomb(1 : disease , 1 : viral , 2 : 4 , 1 : 7 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%hpvInds = matlabFunction(hpvInds));

% HIV and all HPV+
hivHpvInds = toInd(allcomb(2 : 6 , 1: viral , 2 : 4 , 1 : 7 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%hivHpvInds = matlabFunction(hivHpvInds));

% hrHPV
hrHpvInds = toInd(allcomb(1 : disease , 1 : viral , 2 , 1 : 7 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%hrHpvInds = matlabFunction(hrHpvInds));

% hrHPV and HIV
hrHpvHivInds = toInd(allcomb(2 : 6 , 1 : viral , 2 , 1 : 7 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%hrHpvHivInds = matlabFunction(hrHpvHivInds));

% lrHPV
lrHpvInds = toInd(allcomb(1 : disease , 1 : viral , 3 , 1 : 7 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%lrHpvInds = matlabFunction(lrHpvInds));

% lrHPV and HIV
lrHpvHivInds = toInd(allcomb(2 : 6 , 1 : viral , 3 , 1 : 7 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%lrHpvHivInds = matlabFunction(lrHpvHivInds));

% coHPV
coHpvInds = toInd(allcomb(1 : disease , 1 : viral , 4 , 1 : 7 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%coHpvInds = matlabFunction(coHpvInds));

% coHPV and HIV
coHpvHivInds = toInd(allcomb(2 : 6 , 1 : viral , 4 , 1 : 7 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%coHpvHivInds = matlabFunction(coHpvHivInds));


% CIN1
cin1Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : 4 , 2 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%cin1Inds = matlabFunction(cin1Inds));

% CIN1 and HIV
cin1HivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : 4 , 2 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%cin1HivInds = matlabFunction(cin1HivInds));

% CIN2
cin2Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : 4 , 3 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%cin2Inds = matlabFunction(cin2Inds));

% CIN2 and HIV
cin2HivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : 4 , 3 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%cin2HivInds = matlabFunction(cin2HivInds));

% CIN3
cin3Inds = toInd(allcomb(1 : disease , 1 : viral , 1 : 4 , 4 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%cin3Inds = matlabFunction(cin3Inds));

% CIN3 and HIV
cin3HivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : 4 , 4 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%cin3HivInds = matlabFunction(cin3HivInds));

% Cervical cancer
ccInds = toInd(allcomb(1 : disease , 1 : viral , 1 : 4 , 5 : 7 , 1 : periods , ...
    2 , 1 : age , 1 : risk));
%ccInds = matlabFunction(ccInds));

% Cervical cancer and HIV
ccHivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : 4 , 5 : 7 , 1 : periods , ...
    2 , 1 : age , 1 : risk));
%ccHivInds = matlabFunction(ccHivInds));

% Vaccinated
vaxInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 9 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%vaxInds = matlabFunction(vaxInds));

% Vaccinated and HIV
vaxHivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 9 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%vaxHivInds = matlabFunction(vaxHivInds));

% Hysterectomy
hystInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 8 , 1 : periods , 2 , ...
    1 : age , 1 : risk));
%hystInds = matlabFunction(hystInds));

% Hysterectomy and HIV
hystHivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 8 , 1 : periods , 2 , ...
    1 : age , 1 : risk));
%hystHivInds = matlabFunction(hystHivInds));

% Immune
immInds = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 10 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%immInds = matlabFunction(immInds));

% Immune and HIV
immHivInds = toInd(allcomb(2 : 6 , 1 : viral , 1 : hpvTypes , 10 , 1 : periods , ...
    1 : gender , 1 : age , 1 : risk));
%immHivInds = matlabFunction(immHivInds));

save('groupedInds', 'artInds' , 'cd4_500Inds' , 'hivInds' , 'hpvInds' , 'hivHpvInds' ,  'hrHpvInds' , ...
    'hrHpvHivInds' , 'lrHpvInds' , 'lrHpvHivInds' , 'coHpvInds' , 'coHpvHivInds' , ...
    'cin1Inds' , 'cin1HivInds' , 'cin2Inds' , 'cin2HivInds' , 'cin3Inds' , 'cin3HivInds' , ...
    'ccInds' , 'ccHivInds' , 'vaxInds' , 'vaxHivInds' , 'hystInds' , 'hystHivInds' , ...
    'immInds' , 'immHivInds');

