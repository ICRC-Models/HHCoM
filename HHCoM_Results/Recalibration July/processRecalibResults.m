% sympParams = []; 
% 
% for l = 0.001 : 0.002 : 0.01 % changed this
%     for r = 0.014 : 0.002 : 0.02
%         for d = 0.6 : 0.1 : 0.8 
%             for l_r = 0.06 : 0.02 : 0.14
%                 for r_d = 0.025 : 0.02 : 0.1
% 
%                     sympParams = [sympParams; l r d l_r r_d]; 
% 
%                 end
%             end
%         end 
%     end 
% end 

files = dir('*.mat'); 

% Initialize results matrices 
results = zeros(length(files), 3+4+5); 

% For loop for reading in results
for i = 1 : length(files)
    fullFileName = files(i).name;
    curr = load(fullFileName); 

    results(i, 1:3) = curr.stageDist(1:end); 
    results(i, 4:7) = curr.numDxCC(1:end); 
    results(i, 8:end) = [curr.kSymp curr.kRL curr.kDR]; 
end 

% ultimately, you want both the stageDist and numDxCC to be the target numbers 

% turn into array, add labels to each column, export as CSV because it's easier to process in R
results2 = array2table(results, 'VariableNames', {'locDist', 'regDist', 'disDist', ...
        'locCount', 'regCount', 'disCount', 'totCount', ...
        'loc_kSymp', 'reg_kSymp', 'dis_kSymp', 'locToReg', 'regToDis'}); 

writetable(results2, ['recalibration_kSymp_transProb_14Jul23.csv']);

