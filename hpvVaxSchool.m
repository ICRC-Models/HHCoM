% HPV school-based vaccination (assumes girls are not also screened in vaccination age group)
function[dPop , hpvVaxd , hpvVaxdNonPrev , hpvVaxdPrev] = hpvVaxSchool(pop , disease , viral , risk , ...
    hpvVaxStates , hpvNonVaxStates , endpoints , intervens , vaxG , vaxAge , ...
    vaxRate_vec , toInd , vaxYrs , year , stepsPerYear , gradScaleUp)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
hpvVaxd = 0;
hpvVaxdNonPrev = 0; 
hpvVaxdPrev = 0; 

if gradScaleUp == 1 % note that gradScaleUp has not been set up for future sim, only historical !!!!!!!
        % Vaccination level
        if year >= vaxYrs(1) && year < vaxYrs(2)
            periodInd = 1;
        elseif year >= vaxYrs(2) && year < vaxYrs(3)
            periodInd = 2;
        else
            periodInd = 2; 
        end
        
            dataYr1 = vaxYrs(1);
            dataYrLast = vaxYrs(size(vaxYrs , 1));
            baseYrInd = max(find(year >= vaxYrs , 1, 'last') , 1); % get index of first year <= current year
            baseYr = vaxYrs(baseYrInd);
            vaxRate = vaxRate_vec{periodInd}(1); % vax coverage up to 1st year
            if year < dataYrLast && year > dataYr1 % vax coverage between 1st and last year
                vaxRate = vaxRate_vec{periodInd}(round((year - baseYr) * stepsPerYear) + 1);
            elseif year >= dataYrLast % vax coverage last year and after
                lastInd = size(vaxRate_vec , 1);
                vaxRate = vaxRate_vec{periodInd}(size(vaxRate_vec{periodInd} , 2));
            end 
else 
    vaxRate = vaxRate_vec; 
end 

%% Apply school-based vaccination regimen
for d = 1 : disease
    for v = 1 : viral
        for g = min(vaxG) : max(vaxG) 
            for r = 1 : risk
                for aV = 1 : length(vaxAge)
                    a = vaxAge(aV);
                    fromNonVSus = toInd(allcomb(d , v , 1 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ... 
                        g , a , r));
                    fromNonVImm = toInd(allcomb(d , v , 7 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ... 
                        g , a , r));
                    fromNonVInf = toInd(allcomb(d , v , 2 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ...
                        g , a , r)); %include infected and CIN1
                    fromNonVCin1 = toInd(allcomb(d , v, 3 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ...
                        g , a, r)); 
                    fromNonVCin2 = toInd(allcomb(d , v, 4 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ...
                        g , a, r)); 
                    fromNonVCin3 = toInd(allcomb(d , v, 5 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ...
                        g , a, r)); 
                    fromNonVCc = toInd(allcomb(d , v, 6 , 1 : hpvNonVaxStates , 1 : 3 , 1 , ...
                        g , a, r)); 
                    toVSus = toInd(allcomb(d , v , 1 , 1 : hpvNonVaxStates , 1 : 3 , 2 , ...
                        g , a , r));
                    toVImm = toInd(allcomb(d , v , 7 , 1 : hpvNonVaxStates , 1 : 3 , 2 , ...
                        g , a , r));
                    toVInf = toInd(allcomb(d , v , 2 , 1 : hpvNonVaxStates , 1 : 3 , 5 , ...
                        g , a , r)); %into a vaccinated while HPV+ compartment, assigned to vaxxed in intervens
                    toVCin1 = toInd(allcomb(d , v, 3 , 1 : hpvNonVaxStates , 1 : 3 , 5 , ...
                        g , a, r)); 
                    toVCin2 = toInd(allcomb(d , v, 4 , 1 : hpvNonVaxStates , 1 : 3 , 5 , ...
                        g , a, r)); 
                    toVCin3 = toInd(allcomb(d , v, 5 , 1 : hpvNonVaxStates , 1 : 3 , 5 , ...
                        g , a, r)); 
                    toVCc = toInd(allcomb(d , v, 6 , 1 : hpvNonVaxStates , 1 : 3 , 5 , ...
                        g , a, r)); 
                    otherV = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                        1 : endpoints , 2 , g , a , r));
                    allVNonV = toInd(allcomb(d , v , 1 : hpvVaxStates , 1 : hpvNonVaxStates , ...
                        1 : endpoints , 1 : intervens , g , a , r)); 
                    
                    fracVaxd = sumall(pop(otherV)) / ... % find proportion of population that is currently vaccinated
                        sumall(pop(allVNonV));
                    if vaxRate - fracVaxd > 10 ^ -6 % when proportion vaccinated is below target vaccination level
                        vaxCover = max(0 , (vaxRate - fracVaxd) ./ (1 - fracVaxd)); % vaccinate enough people in age group to reach target
                        vaxdGroupSus = vaxCover .* pop(fromNonVSus);
                        vaxdGroupImm = vaxCover .* pop(fromNonVImm);
                        vaxdGroupInf = vaxCover .* pop(fromNonVInf); % add vax of infected individuals
                        vaxdGroupCin1 = vaxCover .* pop(fromNonVCin1);
                        vaxdGroupCin2 = vaxCover .* pop(fromNonVCin2);
                        vaxdGroupCin3 = vaxCover .* pop(fromNonVCin3);
                        vaxdGroupCc = vaxCover .* pop(fromNonVCc);
                        dPop(fromNonVSus) = dPop(fromNonVSus) - vaxdGroupSus;
                        dPop(fromNonVImm) = dPop(fromNonVImm) - vaxdGroupImm;
                        dPop(fromNonVInf) = dPop(fromNonVInf) - vaxdGroupInf; % remove vax of infected individuals 
                        dPop(fromNonVCin1) = dPop(fromNonVCin1) - vaxdGroupCin1; 
                        dPop(fromNonVCin2) = dPop(fromNonVCin2) - vaxdGroupCin2; 
                        dPop(fromNonVCin3) = dPop(fromNonVCin3) - vaxdGroupCin3; 
                        dPop(fromNonVCc) = dPop(fromNonVCc) - vaxdGroupCc; 
                        dPop(toVSus) = dPop(toVSus) + vaxdGroupSus;
                        dPop(toVImm) = dPop(toVImm) + vaxdGroupImm;
                        dPop(toVInf) = dPop(toVInf) + vaxdGroupInf; % add vax of infected individuals 
                        dPop(toVCin1) = dPop(toVCin1) + vaxdGroupCin1; 
                        dPop(toVCin2) = dPop(toVCin2) + vaxdGroupCin2; 
                        dPop(toVCin3) = dPop(toVCin3) + vaxdGroupCin3; 
                        dPop(toVCc) = dPop(toVCc) + vaxdGroupCc; 
                        hpvVaxd = hpvVaxd + sumall(vaxdGroupSus) + sumall(vaxdGroupImm) + sumall(vaxdGroupInf) + sumall(vaxdGroupCin1) + sumall(vaxdGroupCin2) + sumall(vaxdGroupCin3) + sumall(vaxdGroupCc); % count number of people vaccinated at current time step
                        hpvVaxdNonPrev = hpvVaxdNonPrev + sumall(vaxdGroupSus) + sumall(vaxdGroupImm); 
                        hpvVaxdPrev = hpvVaxdPrev + sumall(vaxdGroupInf) + sumall(vaxdGroupCin1) + sumall(vaxdGroupCin2) + sumall(vaxdGroupCin3) + sumall(vaxdGroupCc); 
                    end
                end
            end
        end
    end
end

%% Convert dPop to a column vector for output to ODE solver
dPop = sparse(dPop);
