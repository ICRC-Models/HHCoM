files = dir('*.mat'); 

% Initialize results matrices 
results = []; 

% For loop for reading in results
for i = 1 : length(files)
    fullFileName = files(i).name;
    curr = load(fullFileName); 

    results_temp = zeros(length(curr.numDxCC), 10); 

    results_temp(:, 1:4) = curr.numDxCC(1:end, 1:4); 
    results_temp(:, 5) = curr.numDxCC(1:end, 5); 
    results_temp(:, 6:end) = repmat([curr.kSymp curr.negSumLogL i], size(results_temp, 1), 1); 

    results = [results; results_temp]; 
end 

% ultimately, you want both the stageDist and numDxCC to be the target numbers 

% turn into array, add labels to each column, export as CSV because it's easier to process in R
results2 = array2table(results, 'VariableNames', {'year', 'locNum', 'regNum', 'disNum', ...
        'ccMortNum', ...
        'loc_kSymp', 'reg_kSymp', 'dis_kSymp', 'logLike', 'index'}); 

writetable(results2, ['recalibration_kSymp_transProb_18Jul23.csv']);

