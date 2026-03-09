% HPV school-based vaccination (assumes girls are not also screened in vaccination age group)
function[dPop , hpvVaxd] = hpvVaxSchool(pop , disease , viral , risk , ...
    hpvVaxStates , hpvNonVaxStates , endpoints , intervens , vaxG , vaxAge , ...
    vaxRate_vec , toInd , vaxYrs , year , stepsPerYear , gradScaleUp)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
hpvVaxd = 0;
% **Uncomment when want school vax with different vax eff for HIV -/+; This
% sets up tracking number of vaccines given to each group**

%hpvVaxd1 = 0; %HIV neg vaccination 
%hpvVaxd2 = 0; % HIV positive

if gradScaleUp == 1 % note that gradScaleUp has not been set up for historical sim, only future !!!!!!!
        % Vaccination level
        if year >= vaxYrs(1) && year < vaxYrs(2)
            periodInd = 1;
        else
            periodInd = 2; 
        end
        
            dataYr1 = vaxYrs(1);
            dataYrLast = vaxYrs(size(vaxYrs , 1));
            baseYrInd = max(find(year >= vaxYrs , 1, 'last') , 1); % get index of first year <= current year
            baseYr = vaxYrs(baseYrInd);
            if year < dataYrLast && year > dataYr1 % vax coverage between 1st and last year
%                 vaxRate = vaxRate_vec{periodInd}(1); % vax coverage up to 1st year
                vaxRate = vaxRate_vec{periodInd}(round((year - baseYr) * stepsPerYear) + 1);
            elseif year >= dataYrLast % vax coverage last year and after
                lastInd = size(vaxRate_vec , 1);
                vaxRate = vaxRate_vec{2}(size(vaxRate_vec{2} , 2));
            end 
else 
    vaxRate = vaxRate_vec; 
end 

%% Apply school-based vaccination regimen
for d = 1 : disease

%Classify disease state ** UNCOMMENT WHEN RUNNING DIFFERENT VAX EFF FOR HIV -/+ 
    %TRACKED FOR CU
%isHIVneg = (d == 1 || d == 2); 
%isHIVpos = (d >= 3 && d <= 8);

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

                        % **Comment out when tracking hiv-/+ seperate vacciantion
                        hpvVaxd = hpvVaxd + sumall(vaxdGroupSus) + sumall(vaxdGroupImm); % count number of people vaccinated at current time step

%% ** Uncomment when tracking hiv -/+ vaccination, groups vaccines given to
% those + or - **                      
%numVaxd = sumall(vaxdGroupSus) + ...
          %sumall(vaxdGroupImm);

% Update overall vaccination count
%hpvVaxd = hpvVaxd + numVaxd;

 %Track totals by HIV status
%if isHIVneg
  %  hpvVaxd1 = hpvVaxd1 + numVaxd;   % disease states 1–2
%elseif isHIVpos
 % hpvVaxd2 = hpvVaxd2 + numVaxd;   % disease states 3–8
%end
                    end
                end
            end
        end
    end
end

%% Convert dPop to a column vector for output to ODE solver
dPop = sparse(dPop);
