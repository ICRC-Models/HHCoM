names = {'male partners mult' , 'female partners mult' , 'condom use' , 'mixing by age' , ...
    'mixing by risk' , 'male acts mult' , 'female acts mult' , 'HPV transmission rate' , 'natural immunity mult' , ...
    'ART transmission mult' , 'CIN3 to CC mult'};

numParams = 11;
t = [1:numParams];
baseline = [1.0 , 1.0 , 0.25 , 0.3 , 0.3 , 1.0 , 1.0, 0.0045 , 1.0 , 1.0 , 1.0];
lb = [0.2 , 0.2 , 0.11 , 0.1 , 0.1 , 0.2 , 0.2 0.001 , 0.5 , 1.0 , 0.2];
ub = [5.0 , 5.0 , 0.85 , 1.0 , 1.0 , 2.0 , 2.0 , 1.0 , 1.01 , 2.32 , 5.0];
min = [0.27 , 0.23 , 0.30 , 0.21 , 0.11 , 0.23 , 0.26 , 0.0045 , 0.56 , 1.19 , 1.11];
max = [4.97 , 4.57 , 0.49 , 0.59 , 0.95 , 1.96 , 1.61 , 0.0271 , 0.92 , 2.14 , 4.76];

%plot(baseline,'k*');
errorbar(t,baseline,baseline-lb,ub-baseline,'k*')
hold all;
plot(min,'ro');
hold all;
plot(max,'ro');
ylim([0 6]);
xlim([0 12])
%set(gca , 'xtick' , 1 : length(numParams) , 'xtickLabel' , names);
