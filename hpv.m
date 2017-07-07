% HPV progression
% Simulates progression through HPV states as well as HPV treatment and
% hysterectomy.
% Accepts a population matrix as input and returns dPop, a vector of
% derivatives that describes the change in the population's subgroups due
% to HPV progression.
function[dPop , extraOut] = hpv(t , pop , immuneInds , infInds , cin1Inds , ...
    cin2Inds , cin3Inds , normalInds , ccInds ,  kInf_Cin1 , kInf_Cin2 , ...
    kCin1_Cin2 , kCin1_Cin3 , kCin2_Cin3 , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , ...
    kCin1_Inf , kCin2_Inf , kCin3_Cin1 , kNormal_Cin1 , kNormal_Cin2 , ...
    rNormal_Inf , hpv_hivClear , c3c2Mults , c2c1Mults , fImm ,...
    disease , viral , age , hpvTypes , hpvStates , hystOption)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumall = @(x) sum(x(:));
hyst = 0;
if strcmp(hystOption , 'on')
    hyst = 1;
end
ccInc = zeros(disease , viral , hpvTypes , age);
% constants
% see model notes for index values
% leep (effective treatment rate by leep)
rImmune = 0.024; % for HPV16, Johnson
%rNormal_Inf = rNormal_Inf * 0.5; % adjust rate of inf -> immunity downward
% fImm(1 : 3) = 1; %0.1;
% fImm(4 : age) = 0.58; % (0.48; 0.27 , 0.69) fraction fully protected by immunity based on RR of natural immunity (Beachler, 2017)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dPop = zeros(size(pop));
for d = 1 : disease
    c3c2Mult = c3c2Mults(1); % multiplier used for CIN2 -> CIN3 in HIV infecteds
    c2c1Mult = c2c1Mults(1); % multiplier used for CIN1 -> CIN2 in HIV infecteds
    c1c2Mult = 1; % CIN2 -> CIN1 regression multiplier
    c2c3Mult = 1; % CIN3 -> CIN2 regression multiplier
    rHivHpvMult = 1; % for HIV negative
    if d > 2 && d < 7 % CD4 > 500 -> CD4 < 200
        c3c2Mult = c3c2Mults(d - 2); % CIN2 -> CIN3 multiplier
        c2c1Mult = c2c1Mults(d - 2); % CIN1 -> CIN2 multiplier
        c1c2Mult = hpv_hivClear(d - 2); % CIN2 -> CIN1 regression multiplier
        c2c3Mult = hpv_hivClear(d - 2); % CIN3 -> CIN2 regression multiplier
        rHivHpvMult = hpv_hivClear(d - 2); % Infection clearance multiplier
    end
    for v = 1 : viral       
        for h = 2 : hpvTypes % infected onwards
            for a = 1 : age
                % Infected group
                immuneM = immuneInds(d , v , h , 1 , a , :);
                immuneF = immuneInds(d , v , h , 2 , a , :);
                infM = infInds(d , v , h , 1 , a , :);
                infF = infInds(d , v  , h , 2 , a , :);
                cin1 = cin1Inds(d , v , h , a , :);
                cin2 = cin2Inds(d , v , h , a , :);
                cin3 = cin3Inds(d , v , h , a , :);
                normalM = normalInds(d , v , 1 , a , :);
                normalF = normalInds(d , v , 2 , a , :);
                % to immune from HPV infected, CIN1, CIN2, CIN3 (remove immunity for now)
                
                dPop(normalF) = dPop(normalF) + rImmune * pop(immuneF) ...
                    + kNormal_Cin2(a , h - 1) * (1 - fImm(a)) .* pop(cin2)... % CIN2 -> Normal
                    + kNormal_Cin1(a , h - 1) * (1 - fImm(a)) .*pop(cin1)... % CIN1 -> Normal
                    + rNormal_Inf(a , h - 1) * (1 - fImm(a)) * rHivHpvMult .* pop(infF); % Inf -> Normal
                
                dPop(normalM) = dPop(normalM) + rImmune * pop(immuneM) ...% Immune -> Normal
                    + rNormal_Inf(a , h - 1) * (1 - fImm(a)) * rHivHpvMult .* pop(infM);  % Infected -> Normal
                
                dPop(immuneF) = dPop(immuneF) + kNormal_Cin2(a , h - 1) * fImm(a) .* pop(cin2)... %CIN2 -> immune
                    + kNormal_Cin1(a , h - 1) * fImm(a) .* pop(cin1)... % CIN1 -> immune
                    + rNormal_Inf(a , h - 1) * fImm(a) * rHivHpvMult .* pop(infF)... % Inf -> immune
                    - rImmune * pop(immuneF); % immune -> normal
                
                dPop(immuneM) = dPop(immuneM) + rNormal_Inf(a , h - 1) ...
                    * fImm(a) * rHivHpvMult * pop(infM)...; % infected -> immune
                    - rImmune * pop(immuneM); % immune -> normal
                
                % infected -> CIN1
                dPop(infF) = dPop(infF) ...
                    + kInf_Cin2(a , h - 1) * pop(cin2)... % CIN2 -> infF
                    + kInf_Cin1(a , h -1) * pop(cin1)... % CIN1 -> infF
                    - (kCin1_Inf(a , h - 1) ... % progression to CIN1 from infected females
                    + kCin2_Inf(a , h - 1) ... % progression to CIN2 from infected females (fast progressors)
                    + rNormal_Inf(a , h - 1) * rHivHpvMult) .* pop(infF); % regression to immune from infected females
                
                dPop(infM) = dPop(infM) - rNormal_Inf(a , h - 1) * rHivHpvMult * pop(infM); % regression to immune from infected males              
                
                % Infection and CIN progression in females only
                % kCin_Inf(stage , hpvType , age group)
                % CIN1 group
                dPop(cin1) = dPop(cin1) + kCin1_Inf(a , h - 1) .* pop(infF)... % infected -> CIN1
                    + kCin1_Cin2(a , h - 1) * c1c2Mult .* pop(cin2)... % CIN2 -> CIN1
                    + kCin1_Cin3(a , h - 1) * c2c3Mult .* pop(cin3) ... % CIN3 -> CIN1 (Fast regressors)
                    - (kCin2_Cin1(a , h - 1) * c2c1Mult + kCin3_Cin1(a , h - 1)... % CIN1 -> CIN2, CIN1 -> CIN3 (Fast Progressors)
                    + kInf_Cin1(a , h - 1) + kNormal_Cin1(a , h - 1)) .* pop(cin1); % CIN1 -> Infected , CIN1 -> Normal/immuneF
                
                % CIN2 group
                dPop(cin2) = dPop(cin2) + kCin2_Inf(a , h - 1) .* pop(infF)... % Infected -> CIN2 (Fast progressors)
                    + kCin2_Cin1(a , h - 1) * c2c1Mult * pop(cin1) ... % CIN1 -> CIN2
                    + kCin2_Cin3(a , h - 1) * c2c3Mult * pop(cin3) ... % CIN3 -> CIN2
                    - (kCin3_Cin2(a , h - 1) * c3c2Mult... % CIN2 -> CIN3
                    + kInf_Cin2(a , h - 1) + kNormal_Cin2(a , h - 1) + kCin1_Cin2(a , h - 1)) .* pop(cin2) ; % CIN2 -> Inf, CIN2 -> Normal/immuneF , CIN2 -> CIN1
                
                % CIN3 group
                dPop(cin3) = dPop(cin3) + kCin3_Cin2(a , h - 1) * c3c2Mult .* pop(cin2)... %CIN2 -> CIN3
                    + kCin3_Cin1(a , h - 1) * pop(cin1) +... % CIN1 -> CIN3
                    - (kCC_Cin3(a , h - 1)... % CIN3 -> CC
                    + kCin1_Cin3(a , h - 1) * c2c3Mult + kCin2_Cin3(a , h - 1))...
                    * c2c3Mult .* pop(cin3); % CIN3 -> CIN1 (Fast regressors), CIN3 -> CIN2
                
                % CC group
                cc = ccInds(d , v , h , a , :);
                dPop(cc) = dPop(cc) + kCC_Cin3(a , h - 1) .* pop(cin3); % CIN3 -> CC
                
                % CC incidence tracker
                ccInc(d , v , h , a) = sumall(dPop(cc));
            end
        end
    end
end
extraOut{1} = ccInc;