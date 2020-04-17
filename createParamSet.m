function [] = createParamSet(tstep_abc , date_abc)
t_curr = tstep_abc;
date = date_abc;

%% Load parameters
paramDir = [pwd ,'/Params/'];
paramSetMatrix = load([paramDir,'paramSets_calib_' , date , '_' , num2str(t_curr) , '.dat']);
pIdx = load([paramDir,'pIdx_calib_' , date , '_0.dat']);

% paramSetNew = [paramSetMatrix(1:18 , 672); 0.01; ...       % t=0, 672-mod
%     paramSetMatrix(20 , 672); ...
%     paramSetMatrix(21:24 , 672); ...
%     0.8; 2.0; 1.0; 5.0; ...
%     paramSetMatrix(29:34 , 672); ...
%     1.0; 0.5; 0.0008; ...
%     paramSetMatrix(38 , 672); ...
%     0.85; 2.0];  
% paramSetNew = [paramSetMatrix(1:18 , 2414); 0.01; ...     % t=0, 2414-mod
%     paramSetMatrix(20:39 , 2414); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:22 , 7867); ...           % t=0, 7867-mod
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     paramSetMatrix(29:32 , 7867); ...
%     0.5; 0.5; 1.0; 1.0; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:15 , 7867); 0.4; ...     % t=0, 7867-mod2
%     paramSetMatrix(17:18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     paramSetMatrix(29:32 , 7867); ...
%     0.5; 0.2; 0.5; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1 , 7867).*0.5; ...        % t=0, 7867-mod3
%     paramSetMatrix(2:4 , 7867); ...
%     paramSetMatrix(5:6 , 7867).*1.10; ...
%     paramSetMatrix(7:8 , 7867); ...
%     paramSetMatrix(9:10 , 7867).*1.10; ...
%     paramSetMatrix(11:15 , 7867); 0.4; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     1.2; 0.8; ...
%     paramSetMatrix(31:32 , 7867); ...
%     0.5; 0.2; 0.5; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.5; ...      % t=0, 7867-mod4
%     paramSetMatrix(3:4 , 7867); ...
%     paramSetMatrix(5:6 , 7867).*1.10; ...
%     paramSetMatrix(7:8 , 7867); ...
%     paramSetMatrix(9:11 , 7867).*1.25; ...
%     paramSetMatrix(12:15 , 7867); 0.4; ...
%     0.25; paramSetMatrix(18 , 7867); 0.05; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     1.2; 0.5; ...
%     paramSetMatrix(31:32 , 7867); ...
%     0.5; 0.2; 0.8; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.2; ...        % t=0, 7867-mod5
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.90; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.4; ...
%     0.25; paramSetMatrix(18 , 7867); 0.05; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.2; 0.5; ...
%     paramSetMatrix(31:32 , 7867); ...
%     0.8; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod6
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.75; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.25; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.8; 1.0; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod7
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.4; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.25; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.4; 0.5; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod8
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.15; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.4; 0.5; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod9
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.15; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.2; 0.25; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod10
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.15; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     0.9; 0.2; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod11
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.15; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     0.7; 0.2; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod12
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     0.7; 0.2; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod13
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.8; 1.4; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod14
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.3; 1.1; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod15
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:11 , 7867).*1.50; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.3; 1.1; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod16
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     paramSetMatrix(9:10 , 7867).*1.75; ...
%     paramSetMatrix(11 , 7867).*1.60; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.5; 1.1; ...    
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod17
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     36.0; ...
%     48.0; ...
%     30.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.9; 1.3; ...    
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod18
%     paramSetMatrix(3 , 7867).*1.05; ...
%     paramSetMatrix(4 , 7867).*0.4; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     10.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     36.0; ...
%     48.0; ...
%     33.0; ...
%     21.0; ...
%     paramSetMatrix(13:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.008; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.9; 1.3; ...    
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2; 
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod19
%     paramSetMatrix(3 , 7867).*1.05; ...
%     paramSetMatrix(4 , 7867).*0.4; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     10.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     2.0; ...
%     3.0; ...
%     13.0; ...
%     13.0; ...
%     paramSetMatrix(13:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.008; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.9; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod20
%     paramSetMatrix(3 , 7867).*1.05; ...
%     paramSetMatrix(4 , 7867).*0.4; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     10.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     20.0; ...
%     30.0; ...
%     22.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.008; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.9; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod21
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     36.0; ...
%     48.0; ...
%     30.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.9; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     1.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod22
%     paramSetMatrix(3 , 7867).*1.10; ...
%     paramSetMatrix(4 , 7867).*0.6; ...
%     paramSetMatrix(5 , 7867).*1.10; ...
%     paramSetMatrix(6:8 , 7867).*1.25; ...
%     36.0; ...
%     48.0; ...
%     30.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.9; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod23
%     paramSetMatrix(3 , 7867).*1.10; ...
%     16.0; ...
%     24.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     36.0; ...
%     48.0; ...
%     30.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.9; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod24
%     paramSetMatrix(3 , 7867).*1.10; ...
%     16.0; ...
%     24.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     26.0; ...
%     42.0; ...
%     27.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.9; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod25
%     paramSetMatrix(3 , 7867).*1.10; ...
%     16.0; ...
%     24.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     15.0; ...
%     32.0; ...
%     27.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20:22 , 7867); ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     1.9; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod26
%     paramSetMatrix(3 , 7867).*1.10; ...
%     16.0; ...
%     24.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     15.0; ...
%     32.0; ...
%     27.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.5; 1.0; ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     2.0; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod27
%     paramSetMatrix(3 , 7867).*1.25; ...
%     16.0; ...
%     24.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     7.0; ...
%     24.0; ...
%     27.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.3; 1.3; ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     2.0; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.1; ...        % t=0, 7867-mod28
%     paramSetMatrix(3 , 7867).*1.3; ...
%     22.0; ...
%     24.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.3; 1.3; ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     2.0; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.05; ...        % t=0, 7867-mod29
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     24.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     paramSetMatrix(12:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.3; 1.3; ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     2.0; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod30
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     paramSetMatrix(14:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.3; 1.3; ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     2.0; 1.3; ...
%     paramSetMatrix(31:32 , 7867); ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod31
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     paramSetMatrix(14:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     2.0; 1.3; ...
%     1.3; 0.5; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod32
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     paramSetMatrix(14:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     2.2; 1.3; ...
%     1.3; 0.6; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod33
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     paramSetMatrix(14:15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     2.4; 1.3; ...
%     1.3; 0.7; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod34
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     2.4; 1.3; ...
%     1.3; 0.7; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod35
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 0.8; 5.0; ...
%     2.3; 1.25; ...
%     1.3; 0.7; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod36
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.55; 1.38; ...
%     1.3; 0.7; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod37
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.55; 1.38; ...
%     1.3; 0.7; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37 , 7867); ...
%     0.33; 0.70; ...
%     5.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod38
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.55; 1.38; ...
%     1.3; 0.7; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     4.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod39
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.55; 1.38; ...
%     1.3; 0.7; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     3.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod40
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.55; 1.38; ...
%     1.3; 0.7; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     2.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod41
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.55; 1.38; ...
%     2.0; 1.08; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     2.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod42
%     paramSetMatrix(3 , 7867).*1.3; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     12.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.55; 1.38; ...
%     2.0; 1.08; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     3.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.025; ...        % t=0, 7867-mod43
%     paramSetMatrix(3 , 7867).*1.05; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.5; ...
%     10.0; ...
%     27.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.55; 1.38; ...
%     2.0; 1.08; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     3.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03; ...        % t=0, 7867-mod44
%     paramSetMatrix(3 , 7867).*1.05; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0; ...
%     9.0; ...
%     26.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     2.0; 1.08; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     3.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03; ...        % t=0, 7867-mod45
%     paramSetMatrix(3 , 7867).*1.05; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0; ...
%     9.0; ...
%     26.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     1.8; 0.97; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     3.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03; ...        % t=0, 7867-mod46
%     paramSetMatrix(3 , 7867).*1.05; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0; ...
%     9.0; ...
%     26.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     0.66667; ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     1.8; 0.97; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     3.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03; ...        % t=0, 7867-mod47
%     paramSetMatrix(3 , 7867).*1.05; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0; ...
%     9.0; ...
%     26.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     0.66667; ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     1.8; 0.97; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37 , 7867); ...
%     0.42; 0.70; ...
%     3.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03; ...        % t=0, 7867-mod48
%     paramSetMatrix(3 , 7867).*1.05; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0; ...
%     9.0; ...
%     26.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7876); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     1.8; 0.97; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37 , 7867); ...
%     0.55; 0.48; ...
%     3.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03; ...        % t=0, 7867-mod49
%     paramSetMatrix(3 , 7867).*1.05; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0; ...
%     9.0; ...
%     26.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7876); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     1.8; 0.97; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:38 , 7867); ...
%     0.62; ...
%     3.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03; ...        % t=0, 7867-mod50
%     paramSetMatrix(3 , 7867).*1.05; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0; ...
%     9.0; ...
%     26.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7876); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     1.8; 0.97; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:38 , 7867); ...
%     0.68; ...
%     3.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03; ...        % t=0, 7867-mod51
%     paramSetMatrix(3 , 7867).*1.05; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0; ...
%     9.0; ...
%     26.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     1.8; 0.97; ...
%     1.4; 0.2; 1.0; 0.2; 
%     0.0006;
%     paramSetMatrix(38:39 , 7867); ...
%     3.5];
%paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03.*0.70; ...        % t=0, 7867-mod52
%     paramSetMatrix(3 , 7867).*1.05.*0.70; ...
%     25.0.*0.70; ...
%     25.0.*0.70; ...
%     12.0.*0.70; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0.*0.70; ...
%     9.0.*0.70; ...
%     26.0.*0.70; ...
%     19.0.*0.70; ...
%     14.0.*0.70; ...
%     10.0.*0.70; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     0.65; ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     1.8; 0.97; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     2.0];
%paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03; ...        % t=0, 7867-mod53
%    paramSetMatrix(3 , 7867).*1.05; ...
%    25.0; ...
%    25.0; ...
%    12.0; ...
%    paramSetMatrix(7:8 , 7867).*1.25; ...
%    3.0; ...
%    9.0; ... 
%    26.0; ...
%    19.0; ...
%    14.0; ...
%    10.0; ...
%    paramSetMatrix(15 , 7867); 0.3; ...
%    0.25; paramSetMatrix(18 , 7867); 0.01; ...
%    paramSetMatrix(20 , 7867); ...
%    0.8; 2.0; ...
%    1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%    2.60; 1.41; ...
%    1.8; 0.97; ...
%    1.4; 0.2; 1.0; 0.2;
%    0.000575;
%    paramSetMatrix(38:39 , 7867); ...
%    3.5];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03.*0.75; ...        % t=0, 7867-mod54
%      paramSetMatrix(3 , 7867).*1.05.*0.75; ...
%      25.0.*0.75; ...
%      25.0.*0.75; ...
%      12.0.*0.75; ...
%      paramSetMatrix(7:8 , 7867).*1.25; ...
%      3.0.*0.75; ...
%      9.0.*0.75; ...
%      26.0.*0.75; ...
%      19.0.*0.75; ...
%      14.0.*0.75; ...
%      10.0.*0.75; ...
%      paramSetMatrix(15 , 7867); 0.3; ...
%      0.25; paramSetMatrix(18 , 7867); 0.01; ...
%      0.65; ...
%      0.8; 2.0; ...
%      1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%      2.60; 1.41; ...
%      1.8; 0.97; ...
%      1.4; 0.2; 1.0; 0.2;
%      paramSetMatrix(37:39 , 7867); ...
%      2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03.*0.80; ...        % t=0, 7867-mod55
%      paramSetMatrix(3 , 7867).*1.05.*0.80; ...
%      25.0.*0.80; ...
%      25.0.*0.80; ...
%      12.0.*0.80; ...
%      paramSetMatrix(7:8 , 7867).*1.25; ...
%      3.0.*0.80; ...
%      9.0.*0.80; ...
%      26.0.*0.80; ...
%      19.0.*0.80; ...
%      14.0.*0.80; ...
%      10.0.*0.80; ...
%      paramSetMatrix(15 , 7867); 0.3; ...
%      0.25; paramSetMatrix(18 , 7867); 0.01; ...
%      0.65; ...
%      0.8; 2.0; ...
%      1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%      2.60; 1.41; ...
%      1.8; 0.97; ...
%      1.4; 0.2; 1.0; 0.2;
%      paramSetMatrix(37:39 , 7867); ...
%      2.0];
%paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03.*0.80; ...        % t=0, 7867-mod56
%     paramSetMatrix(3 , 7867).*1.05.*0.80; ...
%     25.0.*0.80; ...
%     25.0.*0.80; ...
%     12.0.*0.80; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0.*0.80; ...
%     9.0.*0.80; ...
%     26.0.*0.80; ...
%     19.0.*0.80; ...
%     14.0.*0.80; ...
%     10.0.*0.80; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; paramSetMatrix(18 , 7867); 0.01; ...
%     0.65; ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     1.8; 0.97; ...
%     1.4; 0.2; 1.0; 0.2;
%     0.0008
%     paramSetMatrix(38:39 , 7867); ...
%     2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03.*0.80; ...        % t=0, 7867-mod57
%      paramSetMatrix(3 , 7867).*1.05.*0.80; ...
%      25.0.*0.80; ...
%      25.0.*0.80; ...
%      12.0.*0.80; ...
%      paramSetMatrix(7:8 , 7867).*1.25; ...
%      3.0.*0.80; ...
%      9.0.*0.80; ...
%      26.0.*0.80; ...
%      19.0.*0.80; ...
%      14.0.*0.80; ...
%      10.0.*0.80; ...
%      paramSetMatrix(15 , 7867); 0.3; ...
%      0.25; paramSetMatrix(18 , 7867); 0.01; ...
%      0.65; ...
%      0.8; 2.0; ...
%      1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%      2.60; 1.41; ...
%      1.8; 0.97; ...
%      1.4; 0.2; 1.0; 0.2;
%      0.0007375
%      paramSetMatrix(38:39 , 7867); ...
%      2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03.*0.80; ...        % t=0, 7867-mod58
%      paramSetMatrix(3 , 7867).*1.05.*0.80; ...
%      25.0.*0.80; ...
%      25.0.*0.80; ...
%      12.0.*0.80; ...
%      paramSetMatrix(7:8 , 7867).*1.25; ...
%      3.0.*0.80; ...
%      9.0.*0.80; ...
%      26.0.*0.80; ...
%      19.0.*0.80; ...
%      14.0.*0.80; ...
%      10.0.*0.80; ...
%      paramSetMatrix(15 , 7867); 0.3; ...
%      0.25; paramSetMatrix(18 , 7867); 0.01; ...
%      0.65; ...
%      0.8; 2.0; ...
%      1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%      2.60; 1.41; ...
%      1.8; 0.97; ...
%      1.4; 0.2; 1.0; 0.2;
%      0.0007
%      paramSetMatrix(38:39 , 7867); ...
%      2.0];
% paramSetNew = [paramSetMatrix(1:2 , 7867).*0.03; ...        % t=0, 7867-mod59
%     paramSetMatrix(3 , 7867).*1.05; ...
%     25.0; ...
%     25.0; ...
%     12.0; ...
%     paramSetMatrix(7:8 , 7867).*1.25; ...
%     3.0; ...
%     9.0; ...
%     26.0; ...
%     19.0; ...
%     14.0; ...
%     10.0; ...
%     paramSetMatrix(15 , 7867); 0.3; ...
%     0.25; 0.3; 0.01; ...
%     paramSetMatrix(20 , 7867); ...
%     0.8; 2.0; ...
%     1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
%     2.60; 1.41; ...
%     1.8; 0.97; ...
%     1.4; 0.2; 1.0; 0.2;
%     paramSetMatrix(37:39 , 7867); ...
%     3.5];
paramSetNew = [paramSetMatrix(2 , 7867).*0.03.*0.5; ...        % t=0, 7867-mod60
    paramSetMatrix(2 , 7867).*0.03; ...       
    paramSetMatrix(3 , 7867).*1.05; ...
    25.0; ...
    25.0; ...
    12.0; ...
    paramSetMatrix(7:8 , 7867).*1.25; ... 
    4.5; ...
    9.0; ...
    26.0; ...
    19.0; ...
    14.0; ...
    10.0; ...
    paramSetMatrix(15 , 7867); 0.3; ...
    0.25; paramSetMatrix(18 , 7867); 0.01; ...
    paramSetMatrix(20 , 7867); ...
    0.8; 2.0; ...
    1.4; 2.0; 1.3; 2.0; 1.0; 5.0; ...
    2.60; 1.41; ...
    1.8; 0.97; ...
    1.4; 0.2; 1.0; 0.2;
    paramSetMatrix(37:39 , 7867); ...
    3.5];




%% Save pIdx and parameter set
file = ['paramSets_calib_' , date , '_' , num2str(t_curr) , '_7867mod60' , '.dat'];
paramDir = [pwd , '/Params/'];
csvwrite([paramDir, file] , paramSetNew)
