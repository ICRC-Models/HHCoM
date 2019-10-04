% HPV school-based vaccination (assumes girls are not also screened in vaccination age group)
function[dPop , hpvVaxd] = hpvVaxSchool(pop , k , disease , viral , risk , ...
    hpvTypes , hpvStates , periods , vaxG , vaxAge , vaxRate)

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
                    fromNonVSus(:,aV) = toInd(allcomb(d , v , 1 , 1 , 1 , ... 
                        g , a , r));
                    fromNonVImm(:,aV) = toInd(allcomb(d , v , 2 , 10 , 1 , ... 
                        g , a , r));
                    fromNonVImmNonV(:,aV) = toInd(allcomb(d , v , 3 , 10 , 1 , ... 
                        g , a , r));
                    toV(:,aV) = toInd(allcomb(d , v , 1 , 9 , 1 , ...
                        g , a , r));
                    otherV(:,aV) = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , ...
                        2 , g , a , r));
                    allVNonV(:,aV) = toInd(allcomb(d , v , 1 : hpvTypes , 1 :hpvStates , ...
                        1 : periods , g , a , r)); 
                    
                    fracVaxd = (sumall(pop(toV(:,aV))) + sumall(pop(otherV(:,aV)))) / ... % find proportion of population that is currently vaccinated
                        (sumall(pop(allVNonV(:,aV))));
                    if vaxRate - fracVaxd > 10 ^ -6 % when proportion vaccinated is below target vaccination level
                        vaxCover = max(0 , (vaxRate - fracVaxd) ./ (1 - fracVaxd)); % vaccinate enough people in age group to reach target
                        vaxdGroupSus = vaxCover .* pop(fromNonVSus(:,aV));
                        vaxdGroupImm = vaxCover .* pop(fromNonVImm(:,aV));
                        vaxdGroupImmNonV = vaxCover .* pop(fromNonVImmNonV(:,aV));
                        dPop(fromNonVSus(:,aV)) = dPop(fromNonVSus(:,aV)) - vaxdGroupSus;
                        dPop(fromNonVImm(:,aV)) = dPop(fromNonVImm(:,aV)) - vaxdGroupImm;
                        dPop(fromNonVImmNonV(:,aV)) = dPop(fromNonVImmNonV(:,aV)) - vaxdGroupImmNonV;
                        dPop(toV(:,aV)) = dPop(toV(:,aV)) + vaxdGroupSus + vaxdGroupImm + vaxdGroupImmNonV;
                        hpvVaxd = hpvVaxd + sumall(vaxdGroupSus) + sumall(vaxdGroupImm) + sumall(vaxdGroupImmNonV); % count number of people vaccinated at current time step
                    end
                end
            end
        end
    end
end

%% Convert dPop to a column vector for output to ODE solver
dPop = sparse(dPop);
