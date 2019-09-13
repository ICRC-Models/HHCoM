function idSizeMatrix(tstep_abc , date_abc)

t_curr = tstep_abc;
date = date_abc;

%% Load all particles
paramDir = [pwd , '/Params/'];
paramSetMatrix = load([paramDir , 'paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat']); % load most recent parameter sample

dim = size(paramSetMatrix,2);

fileD = ['matrixSize_calib_' , date , '_' , num2str(t_curr) , '.dat'];
csvwrite([paramDir, fileD] , dim);



