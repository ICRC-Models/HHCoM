%% summary results
load('general')
load('results')

sumall = @(x) sum(x(:));
disp('Saving HIV results')

% hiv (1 = positive, 0 = negative), age (1-12 categories), sex (0 = female, 1 = male),...
% year (decimal, ie, 1990.25), CD4 (0-9, 0 = uninfected --> 5 = < 200 , 6 =
% HIV-, circ, no Prep , 7 = HIV- circ, on PrEP, 8 = HIV-, uncirc, on PrEP, 9 = ART), ...
% VL (0-5, 0 = uninfected, 1 = acute...5 = > 50,000 , 6 = HIV Positive on ART),...
% ART (0 = not on ART, 1 = on ART)
progressbar('Overall' , 'Disease ' , 'Viral');
row = 1;
dVec = [1 : 6 , 10];
% Collapse to HIV and HIV-susceptible population
hivDims = [length(tVec) , length(dVec) , viral , gender , age];
hivResults = zeros(prod(hivDims) , 8);
for i = 1 : length(tVec)
    for d = 1 : length(dVec)
        hivStat = dVec(d) > 1 && dVec(d) < 7 || dVec(d) == 10;
        art = dVec(d) == 10;
        for v = 1 : 6
            for g = 1 : gender
                for a = 1 : age
                    if dVec(d) == 10 && v == 6 % ART
                        group = toInd(allcomb(10 , 6 , 1 : hpvTypes , ...
                            1 : hpvStates , 1 : periods , g , a , 1 : risk));
                        popGroup = sumall(popVec(i , group));
                        hivResults(row , :) = ...
                            [hivStat , a , g == 1, tVec(i) , -1 , ...
                            -1 , art , popGroup];
                    else
                        group = toInd(allcomb(dVec(d) , v , 1 : hpvTypes , ...
                            1 : hpvStates , 1 : periods , g , a , 1 : risk));
                        popGroup = sumall(popVec(i , group));
                        %  resultRow = sub2ind(hivDims , i , d , v , g , a , r);
                        if dVec(d) == 1
                            group = toInd(allcomb(7 : 9 , v , 1 : hpvTypes , ...
                                1 : hpvStates , 1 : periods , g , a , 1 : risk));
                            popGroup = popGroup + sumall(popVec(i , group));
                        end
                        hivResults(row , :) = ...
                            [hivStat , a , g == 1, tVec(i) , dVec(d) - 1 , ...
                            hivStat * v , art , popGroup];
                    end
                    row = row + 1;
                end
            end
            progressbar([] , [] , v / viral)
        end
        progressbar([] , d / length(dVec) , [])
    end
    progressbar(i / length(tVec) , [] , [])
end
file = [date ' NT_HIV_Results.dat'];
csvwrite(file , hivResults);
disp(['Results saved to ' file]);
%%
% Number of deaths by age, sex, year, HIV status
deathOut = zeros(prod([length(tVec) , disease , viral , gender , age]) , 5);
row = 1;
for i = 1 : length(tVec)
    for d = 1 : disease
        hivStat = d > 1 && d < 7 || d == 10;
        for v = 1 : viral
            for g = 1 : gender
                for a = 1 : age
                        curr = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , 1 : risk));
                        deathOut(row , 1) = sumall(deaths(curr));
                        deathOut(row, 2 : end) = [hivStat , a , g == 1 , tVec(i)];
                        row = row + 1;
                end
            end
        end
    end
end
file = [date 'Deaths.dat'];
csvwrite(file , deathOut);
disp(['Results for deaths saved to ' file]);


% Number of new HIV infections
newHivOut = zeros(prod([length(tVec) , gender , age]) , 4);
row = 1;

for i = 1 : length(tVec)
    for g = 1 : gender
        for a = 1 : age
            newHivOut(row , 1) = sum(newHiv(i , g , a , 1 : risk));
            newHivOut(row, 2 : end) = [a , g == 1 , tVec(i)];
            row = row + 1;
        end
    end
end
file = [date 'HIV_Infections.dat'];
csvwrite(file , newHivOut);
disp(['Results for hiv infections saved to ' file]);

% newHpvOut = zeros(prod([length(tVec) , gender , age]) , 4);
% row = 0;
%
% for i = 1 : length(tVec)
%     for g = 1 : gender
%         for a = 1 : age
%             newHpvOut(row + 1 , 1) = sum(newHpv(i , g , a , 1 : risk));
%             newHpvOut(row + 1, 2 : end) = [a , g == 1 , tVec(i)];
%             row = row + 1;
%         end
%     end
% end
% file = 'HPV_Infections.dat';
% csvwrite(file , hivResults);
% disp(['Results for hiv infections saved to ' file]);
%%
% ART uptake by CD4, viral load, gender, age
row = 1;
artTreatTrackerOut = zeros(prod([length(tVec) , 5 , 5 , gender , age]) , 6);
for i = 1 : length(tVec)
  for d = 2 : 6
    for v = 1 : 5
      for g = 1 : gender
        for a = 1 : age
          artTreatTrackerOut(row , 1) = sum(artTreatTracker(i , d , v , g , a  , 1 : risk));
          artTreatTrackerOut(row , 2 : end) = [d - 1 , v , g == 1 , a , tVec(i)];
          row = row + 1;
        end
      end
    end
  end
end
file = [date '_artTreatTracker.dat'];
csvwrite(file , artTreatTrackerOut)
disp(['Results for artTreatTracker saved to ' file])
disp('Done')
disp(' ')
