% HPV progression
% Simulates HPV screening, treatment, and hysterectomy on women with HPV,
% precancerous lesions , and cervical cancer.
% Accepts a population matrix as input and returns dPop, a matrix of
% derivatives that describes the change in the population's subgroups due
% to HPV treatment and screening.
function[dPop] = hpvTreat(t , pop , disease , viral , hpvTypes , age , ...
    periods , detCC , hivCC , muCC , ccRInds , ccSusInds , ...
    hystPopInds , screenFreq , screenCover , hpvSens , ccTreat , ...
    cytoSens , cin1Inds , cin2Inds , cin3Inds , normalInds , getHystPopInds ,...
    OMEGA , leep , hystOption , year)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hyst = 0;
if strcmp(hystOption , 'on')
    hyst = 1;
end
% constants
% rCin regression
% kCin progression to CIN from HPV infection
% see model notes for index values
% leep (effective treatment rate by leep)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dPop = zeros(size(pop));
for v = 1 : viral
    for d = 1 : disease
        fScreen = screenFreq(: , 2);
        if d >= 2 && d <= 6 || d == 10 % Infected or low CD4 or HIV+ on ART
            fScreen = screenFreq(: , 1);
        end
        hpvEff = screenCover * fScreen(1) * hpvSens; % effective screen and detect rate by hpv test
        cytoEff = screenCover * fScreen(2) * cytoSens; % effective screen and detect rate by cytology
        for h = 2 : hpvTypes % infected onwards
            for a = 5 : age
%                 cin1 = cin1Inds(d , v , h , a , :);
                cin2 = cin2Inds(d , v , h , a , :);
                cin3 = cin3Inds(d , v , h , a , :);
                normal = normalInds(d , v , 2 , a , :);
                
                % CIN2 group detection and treatment
                dPop(cin2) = dPop(cin2)...
                    - (hpvEff(1) * cytoEff(1) * leep) .* pop(cin2); % detection with hpv test and cytology, and treatment with leep
                
                dPop(normal) = dPop(normal) + (hpvEff(1) * cytoEff(1) * leep) .* pop(cin2);
                
                % CIN3 group
                if year <= 2015 % screening before 2015
                    dPop(cin3) = dPop(cin3) ...
                        - (cytoEff(2) * leep) .* pop(cin3);
                    
                    dPop(normal) = dPop(normal) + (cytoEff(2) * leep) .* pop(cin3);
                    
                else % screening after 2015
                    % 35-55 years
                    if a > 7 && a < 12
                        cin3_35Plus = cin3Inds(d , v , h , a , :);
                        dPop(cin3_35Plus) = dPop(cin3_35Plus) ...
                            - (hpvEff(2) * leep) .* pop(cin3_35Plus); % hpv test screen before leep
                        
                        dPop(normal) = dPop(normal) + (hpvEff(2) * leep) .* pop(cin3_35Plus);
                    end
                    % 25-35 years old
                    if a > 5 && a < 8
                        if d == 1 % HIV negative and 25-35 years old
                            screenTreat = cytoEff(2) * hpvEff(2) * leep; % detect by cytology and hpv test then leep
                        else % HIV positive and 25-35 years old
                            screenTreat = hpvEff(2) * leep; % detect by hpv test then leep
                        end
                        cin3_25_35 = cin3Inds(d , v , h , a , :);
                        dPop(cin3_25_35) = dPop(cin3_25_35) ...
                            - (screenTreat) .* pop(cin3_25_35); % test and leep
                        
                        dPop(normal) = dPop(normal) + (screenTreat) .* pop(cin3_25_35);
                    end
                end
                
                % cervical cancer screen and treat, death
                for s = 5 : 7
                    for p = 1 : periods
                        ccR = ccRInds(d , v , h , s , p , a , :); % s refers to region of cervical cancer (local, regional, distant)
                        dPop(ccR) = dPop(ccR) ...
                            - (detCC(s - 4) .* ccTreat + muCC(p , s - 4) * hivCC(s - 4)) .* pop(ccR); % detection probability and death given region. Add in probability of successful treatment later!
                        
                        % CC -> susceptible
                        cc2Sus = ccSusInds(d , v , h , a , :);
                        dPop(cc2Sus) = dPop(cc2Sus)...
                            + pop(ccR) .* detCC(s - 4) .* ccTreat;
                    end
                end
            end
        end
    end
end

% hysterectomy (age dependent rate)
% if hyst
%     for a = 3 : age % can adjust age group range
%         hystPop = hystPopInds(a , :);
%         getHystPop = getHystPopInds(a , :);
%         hyst = OMEGA(a) .* pop(getHystPop); % women (susceptible --> cervical cancer(distant))
%         dPop(getHystPop) = - hyst;
%         dPop(hystPop) = hyst;
%     end
% end

% % Convert to column vector for output to ODE solver
% dPop = reshape(dPop , [prod(dim) , 1]);