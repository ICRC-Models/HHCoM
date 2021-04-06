% HPV school-based vaccination (assumes girls are not also screened in vaccination age group)
function[dPop , hpvVaxd] = hpvVaxSchool(pop , disease , viral , risk , ...
    hpvVaxStates , hpvNonVaxStates , endpoints , intervens , vaxG , vaxAge , ...
    vaxRate , toInd)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
hpvVaxd = 0;

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
                    toVSus = toInd(allcomb(d , v , 1 , 1 : hpvNonVaxStates , 1 : 3 , 2 , ...
                        g , a , r));
                    toVImm = toInd(allcomb(d , v , 7 , 1 : hpvNonVaxStates , 1 : 3 , 2 , ...
                        g , a , r));
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
                        dPop(fromNonVSus) = dPop(fromNonVSus) - vaxdGroupSus;
                        dPop(fromNonVImm) = dPop(fromNonVImm) - vaxdGroupImm;
                        dPop(toVSus) = dPop(toVSus) + vaxdGroupSus;
                        dPop(toVImm) = dPop(toVImm) + vaxdGroupImm;
                        hpvVaxd = hpvVaxd + sumall(vaxdGroupSus) + sumall(vaxdGroupImm); % count number of people vaccinated at current time step
                    end
                end
            end
        end
    end
end

%% Convert dPop to a column vector for output to ODE solver
dPop = sparse(dPop);
