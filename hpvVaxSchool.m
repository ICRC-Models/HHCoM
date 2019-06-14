% HPV school-based vaccination (assumes girls are not also screened in vaccination age group)
function[dPop , hpvVaxd] = hpvVaxSchool(pop , k , disease , viral , risk , ...
    hpvTypes , hpvStates , periods , vaxG , vaxAge , vaxRate)

%% Set constants and initialize vectors
toInd = @(x) (x(: , 8) - 1) * k(7) + (x(: , 7) - 1) * k(6) + (x(: , 6) - 1) * k(5) ...
    + (x(: , 5) - 1) * k(4) + (x(: , 4) - 1) * k(3) + (x(: , 3) - 1) * k(2) ...
    + (x(: , 2) - 1) * k(1) + x(: , 1);
sumall = @(x) sum(x(:));
hpvVaxd = 0;

dPop = zeros(size(pop));

% Apply school-based vaccination regimen
for d = 1 : disease
    for v = 1 : viral
        for g = min(vaxG) : max(vaxG) 
            for r = 1 : risk
                fromNonVSus = toInd(allcomb(d , v , 1 , 1 , 1 , ... 
                    g , vaxAge , r));
                fromNonVImm = toInd(allcomb(d , v , 2 , 10 , 1 , ... 
                    g , vaxAge , r));
                toV = toInd(allcomb(d , v , 1 , 9 , 1 , ...
                    g , vaxAge , r));
                otherV = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , ...
                    2 , g , vaxAge , r));
                allVNonV = toInd(allcomb(d , v , 1 : hpvTypes , 1 :hpvStates , ...
                    1 : periods , g , vaxAge , r)); 

                fracVaxd = (sumall(pop(toV)) + sumall(pop(otherV))) / ... % find proportion of population that is currently vaccinated
                    (sumall(pop(allVNonV)));
                if vaxRate - fracVaxd > 10 ^ -6 % when proportion vaccinated is below target vaccination level
                    vaxCover = max(0 , (vaxRate - fracVaxd) ./ (1 - fracVaxd)); % vaccinate enough people in age group to reach target
                    vaxdGroupSus = vaxCover .* pop(fromNonVSus);
                    vaxdGroupImm = vaxCover .* pop(fromNonVImm);
                    dPop(fromNonVSus) = dPop(fromNonVSus) - vaxdGroupSus;
                    dPop(fromNonVImm) = dPop(fromNonVImm) - vaxdGroupImm;
                    dPop(toV) = dPop(toV) + vaxdGroupSus + vaxdGroupImm;
                    hpvVaxd = hpvVaxd + sumall(vaxdGroupSus) + sumall(vaxdGroupImm); % count number of people vaccinated at current time step
                end
            end
        end
    end
end

dPop = sparse(dPop);
