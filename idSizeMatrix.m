
date = '18July19';
t = 3;

t_prev = t-1;
t_curr = t;
t_next = t+1;

%% Load all particles
paramDir = [pwd , '/Params/'];
paramSetMatrix = load([paramDir , 'paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent parameter sample

dim = size(paramSetMatrix);

fileD = ['matrixSize_calib_' , date , '_' , num2str(t_curr) , '.dat'];
csvwrite([paramDir, fileD] , dim);



