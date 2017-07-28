function[dPop , extraOuts] = hivART(pop , vlAdvancer , artDist , muHIV , hivPositiveArtAll , ...
    kCD4 , disease , viral , gender , age , risk , k , artOut , hivInds , ...
    stepsPerYear , below200Art_2004 , above200Art_2004 , pie4Vec_2004 , ...
    below200Art_2006 , year)

artOut = 0; % for testing
toInd = @(x) (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + (x(: , 6) - 1) * k(5) ...
    + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + (x(: , 3) - 1) * k(2) ...
    + (x(: , 2) - 1) * k(1) + x(: , 1);
sumall = @(x) sum(x(:));
artDist = reshape(artDist, [disease , viral , gender , age , risk]);
% disease related mortality
%
% pie [d x v x g x a x r]
% values obtained from k4_mat in MainMain.m of HIV model.
treat = zeros(disease , viral , gender , age ,risk);
treat(1 , : , : , 4 : end , :) = 0; % rate of going on PrEP
treat(2 , : , : , 4 : end , :) = 0; % ART coverage
treat(3 , : , : , 4 : end , :) = 0; % ART coverage
treat(4 , : , : , 4 : end , :) = 0; % ART coverage
treat(5 , : , : , 4 : end , :) = 0; % ART coverage
treat(6 , : , : , 4 : end, :) = 0; % ART coverage
kOn = ones(6 , 1);
kOn(6) = 0;

% if year >= 2013
%     treat(6 , : , : , 4 : end, :) = below200Art_2006(end); % ART coverage for persons with CD4 < 200
%     treat(5 , : , : , 4 : end , :) = above200Art_2004(end);
%     treat(4 , : , : , 4 : end , :) = pie4Vec_2004(end);
% elseif year >= 2006
%     treat(6 , : , : , 4 : end, :) = below200Art_2006((year - 2006)  * stepsPerYear + 1); % ART coverage for persons with CD4 < 200
% %     pie(5 , : , : , 4 : end , :) = pie5Vec_2004(end);
% %     pie(4 , : , : , 4 : end , :) = pie4Vec_2004(end);
% elseif year >= 2004
%     treat(6 , : , : , 4 : end, :) = below200Art_2004((year - 2004)  * stepsPerYear + 1); % ART coverage for persons with CD4 < 200
% %     pie(5 , : , : , 4 : end , :) = pie5Vec_2004((year - 2004)  * stepsPerYear + 1);
% %     pie(4 , : , : , 4 : end , :) = pie4Vec_2004((year - 2004)  * stepsPerYear + 1);
% end

if year >= 2013
    maxCover = 0.36;
    if year >= 2020
        maxCover = 0.42;
    end
    for g = 1 : gender
        for a = 1 : age
            for r = 1 : risk
                hivPositiveArt = toInd(hivInds(10 , 6 , g , a , r));
                onArt = sum(pop(hivPositiveArt));
                totHivPos = 0;
                for d = 2 : 6
                    for v = 1 : 5
                        hivPositive = toInd(hivInds(d , v , g , a , r));
                        totHivPos = totHivPos + sum(pop(hivPositive));
                    end
                end
                fracART = onArt / (onArt + totHivPos);
                if fracART < maxCover
                    for d = 2 : 6
                        for v = 1 : 5
                            hivPositive = toInd(hivInds(d , v , g , a , r));
                            hivPos = sum(pop(hivPositive));
                            cover = (maxCover - fracART) ./ (1 - fracART);
                            treat(d , v , g , a , r) = max(cover , 0);                         
                        end
                    end
                end
            end
        end
    end
elseif year >= 2006
    maxCover = 0.35;
    for g = 1 : gender
        for a = 1 : age
            for r = 1 : risk
                hivPositiveArt = toInd(hivInds(10 , 6 , g , a , r));
                onArt = sum(pop(hivPositiveArt));
                totHivPos = 0;
                for d = 2 : 6
                    for v = 1 : 5
                        hivPositive = toInd(hivInds(d , v , g , a , r));
                        totHivPos = totHivPos + sum(pop(hivPositive));
                    end
                end
                fracART = onArt / (onArt + totHivPos);
                if fracART < maxCover
                    for d = 2 : 6
                        for v = 1 : 5
                            hivPositive = toInd(hivInds(d , v , g , a , r));
                            hivPos = sum(pop(hivPositive));
                            cover = (maxCover - fracART) ./ (1 - fracART);
                            treat(d , v , g , a , r) = max(cover , 0);                         
                        end
                    end
                end
            end
        end
    end
end

dPop = zeros(size(pop));
artTreat = zeros(disease , viral , gender , age , risk);
hivDeaths = zeros(age , 1);

% Dropout from PrEP
prepOut = 0; % for now



for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
%             if sum(isnan(pop(:))) > 0
%                 disp(['Gender: ' , num2str(g), ' Age : ' , num2str(a) , ' Risk : ' , num2str(r)])
%             end
            % HIV Negative, uncircumcised (d = 1)
            hivPositiveArt = toInd(hivInds(10 , 6 , g , a , r)); % allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , a , r))

            for v = 1 : 5
                %                 artTreat(2 , v , g , a , r) = pie(2 , v , g , a , r) ... % keep track of distribution of people going on ART
                %                         .* sumall(pop(acuteInf)); % going on ART
                %
                % CD4 > 500 cells/uL (d = 3)
                % CD4 500-350 cells/uL (d = 4)
                % CD4 350-200 cells/uL (d = 5)
                % CD4 <200 cells/uL (d = 6)

                for d = 3 : 6
                    % get indices
                    cd4Curr = toInd(hivInds(d , v , g , a , r));
                    % calculate CD4 changes
                    dPop(cd4Curr) = ...
                        artDist(d , v , g , a , r) ... % Distributed dropouts from ART
                        .* pop(hivPositiveArt)...
                        - (treat(d , v , g , a , r))... % going on ART
                        .* pop(cd4Curr);
                    
                    artTreat(d , v , g , a , r) = treat(d , v , g , a , r) ... % keep track of distribution of people going on ART
                        .* sumall(pop(cd4Curr)); % going on ART

                    % HIV-positive going on ART (d = 10)
                    dPop(hivPositiveArt) = ...
                        dPop(hivPositiveArt)...
                        + treat(d , v , g , a , r)... % rate of going on ART
                        .* pop(cd4Curr); % going on ART
                end
            end

            % ART dropout weighted by proportion of individuals with CD4 above/below 200 receiving ART
            % Dropout from ART (d = 10)
             dPop(hivPositiveArt) = ...
                dPop(hivPositiveArt)...
                - artOut .* pop(hivPositiveArt); % artOut to d = 2:6 as determined by distribution matrix
        end
    end
end
extraOuts{1} = artTreat; %reshape(artTreat , [numel(artTreat) , 1]);
dPop = sparse(dPop);
