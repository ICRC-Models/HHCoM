% HPV progression
% Simulates progression through HPV states as well as HPV treatment and
% hysterectomy.
% Accepts a population matrix as input and returns dPop, a vector of
% derivatives that describes the change in the population's subgroups due
% to HPV progression.
function[dPop , extraOut] = hpvCCdet(t , pop , immuneInds , infInds , cin1Inds , ...
    cin2Inds , cin3Inds , normalInds , ccInds , ccRegInds , ccDistInds , ...
    ccTreatedInds , ccLocDetInds , ccRegDetInds , ccDistDetInds , ...
    kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
    kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf  ,...
    rNormal_Inf , hpv_hivClear , c3c2Mults , ...
    c2c1Mults , fImm , kRL , kDR , muCC , muCC_det , kCCDet , ...
    disease , viral , age , hpvTypes , ...
    rImmuneHiv , vaccinated , hystOption)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% kCCDet = kCCDet .* 0 ; %TESTING!
sumall = @(x) sum(x(:));
hyst = 0;
if strcmp(hystOption , 'on')
    hyst = 1;
end
ccInc = zeros(disease , hpvTypes , age);
ccDeath = ccInc;
ccTreated = zeros(disease , hpvTypes , age , 3); % 3 for cancer stages - local, regional, distant
% constants
% see model notes for index values
% leep (effective treatment rate by leep)
rImmune = 0.024; % for HPV16, Johnson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dPop = zeros(size(pop));
for d = 1 : disease
    c3c2Mult = c3c2Mults(1); % multiplier used for CIN2 -> CIN3 in HIV infecteds
    c2c1Mult = c2c1Mults(1); % multiplier used for CIN1 -> CIN2 in HIV infecteds
    c1c2Mult = 1; % CIN2 -> CIN1 regression multiplier
    c2c3Mult = 1; % CIN3 -> CIN2 regression multiplier
    rHivHpvMult = 1; % for HIV negative
    deathCC = muCC(6 , :); % HIV negative undetected CC mortality
    deathCC_det = muCC_det(6 , :); % HIV negative detected CC mortality
    rHiv = 1; % Multiplier for immunity clearance for HIV+
    rHivHpv_Clear = 1;
    if d > 2 && d < 7 % CD4 > 500 -> CD4 < 200
        c3c2Mult = c3c2Mults(d - 2); % CIN2 -> CIN3 multiplier
        c2c1Mult = c2c1Mults(d - 2); % CIN1 -> CIN2 multiplier
        c1c2Mult = hpv_hivClear(d - 2); % CIN2 -> CIN1 regression multiplier
        c2c3Mult = hpv_hivClear(d - 2); % CIN3 -> CIN2 regression multiplier
        rHivHpvMult = hpv_hivClear(d - 2);%hpvClearMult(d - 2); % Regression multiplier
        rHivHpv_Clear = hpv_hivClear(d - 2); % Infection clearance multiplier
        deathCC = muCC(d - 2 , :); % HIV+ CC death rate
        deathCC_det = muCC_det(d - 2 , :);
        rHiv = rImmuneHiv(d - 2); % Multiplier for immunity clearance for HIV+
    elseif d == 10 % CD4 > 500 multipliers for HIV+ on ART
        c3c2Mult = c3c2Mults(1); % CIN2 -> CIN3 multiplier
        c2c1Mult = c2c1Mults(1); % CIN1 -> CIN2 multiplier
        c1c2Mult = hpv_hivClear(1); % CIN2 -> CIN1 regression multiplier
        c2c3Mult = hpv_hivClear(1); % CIN3 -> CIN2 regression multiplier
        rHivHpvMult = hpv_hivClear(1);%hpvClearMult(1); % Regression multiplier
        rHivHpv_Clear = hpv_hivClear(1); % Infection clearance multiplier
        deathCC_det = muCC_det(5 , :);
        deathCC = muCC(5 , :); % HIV+ CC death rate
        rHiv = rImmuneHiv(1); % Multiplier for immunity clearance for HIV+
    end
    for h = 2 % infected onwards
        for a = 1 : age
            for p = 1 : 2
                % Infected group
                immuneM = immuneInds(d , h , 1 , a , p , :);
                immuneF = immuneInds(d , h , 2 , a , p , :);
                infM = infInds(d , h , 1 , a , p , :);
                infF = infInds(d , h , 2 , a , p , :);
                cin1 = cin1Inds(d , h , a , p , :);
                cin2 = cin2Inds(d , h , a , p , :);
                cin3 = cin3Inds(d , h , a , p , :);
                normalM = normalInds(d , 1 , a , p , :);
                normalF = normalInds(d , 2 , a , p , :);
                % to immune from HPV infected, CIN1, CIN2, CIN3 (remove immunity for now)
                
                dPop(normalF) = dPop(normalF) + rImmune * rHiv * pop(immuneF); % immuneF -> normalF
                
                dPop(normalM) = dPop(normalM) + rNormal_Inf(a , h - 1) * rHivHpv_Clear .* pop(infM);  % infM -> normalM
                
                dPop(immuneF) = dPop(immuneF)...
                    + rNormal_Inf(a , h - 1) * rHivHpv_Clear .* pop(infF)... % infF -> immuneF
                    - rImmune * rHiv * pop(immuneF); % immuneF -> normalF
                
                %                 dPop(immuneM) = dPop(immuneM) + rNormal_Inf(a) ...
                %                     * fImm(a) * rHivHpvMult * pop(infM)...; % infected -> immune
                %                     - rImmune * pop(immuneM); % immune -> normal
                
                % infected -> CIN1
                dPop(infF) = dPop(infF) ...
                    + kInf_Cin1(a , h - 1) * pop(cin1)... % CIN1 -> infF
                    - (kCin1_Inf(a , h - 1) + ... % infF -> CIN1
                    rNormal_Inf(a , h - 1) ...
                    * rHivHpv_Clear) .* pop(infF); % infF -> immuneF
                
                dPop(infM) = dPop(infM) - rNormal_Inf(a , h - 1) * rHivHpv_Clear * pop(infM); % regression to normal from infected males
                
                % Infection and CIN progression in females only
                % kCin_Inf(stage , hpvType , age group)
                % CIN1 group
                dPop(cin1) = dPop(cin1) + kCin1_Inf(a , h - 1) .* pop(infF)... % progression to CIN1 from infected females
                    + rHivHpvMult * kCin1_Cin2(a , h - 1) * c1c2Mult .* pop(cin2)... % CIN2 -> CIN1
                    - (kCin2_Cin1(a , h - 1) * c2c1Mult ... % CIN1 -> CIN2
                    + kInf_Cin1(a , h - 1)) .* pop(cin1); % CIN1 -> Infected
                
                % CIN2 group
                dPop(cin2) = dPop(cin2) ...
                    + kCin2_Cin1(a , h - 1) * c2c1Mult * pop(cin1) ... % CIN1 -> CIN2
                    + kCin2_Cin3(a , h - 1) * c2c3Mult * pop(cin3) ... % CIN3 -> CIN2
                    - (kCin3_Cin2(a , h - 1) * c3c2Mult... % CIN2 -> CIN3
                    + kCin1_Cin2(a , h - 1) * rHivHpvMult * c1c2Mult) .* pop(cin2) ; % CIN2 -> CIN1
                
                % CIN3 group
                dPop(cin3) = dPop(cin3) + kCin3_Cin2(a , h - 1) * c3c2Mult .* pop(cin2)... %CIN2 -> CIN3
                    - (kCC_Cin3(a , h - 1)... % CIN3 -> CC
                    + kCin2_Cin3(a , h - 1) * c2c3Mult)... % CIN3 -> CIN2
                    .* pop(cin3);
                
                % CC group
%                 kCCDet = kCCDet .* 0; % TESTING!!!
%                 deathCC = deathCC .* 0; % TESTING!!!
                
                ccLoc = ccInds(d , h , a , p , :);
                locDet = kCCDet(1) * pop(ccLoc);
                dPop(ccLoc) = dPop(ccLoc) + kCC_Cin3(a , h - 1) .* pop(cin3)... % CIN3 -> CC
                    - kRL * pop(ccLoc)... % local -> regional
                    - deathCC(1) * pop(ccLoc)... % local CC mortality
                    - locDet; % detect local CC -> treat
                
                ccReg = ccRegInds(d , h , a , p , :);
                regDet = kCCDet(2) * pop(ccReg);
                dPop(ccReg) = dPop(ccReg) + kRL * pop(ccLoc)...  % local -> regional
                    - kDR * pop(ccReg)... % regional -> distant
                    - deathCC(2) * pop(ccReg)... % regional CC mortality
                    - regDet; % detect regional CC -> treat
                
                ccDist = ccDistInds(d , h , a , p , :);
                distDet = kCCDet(3) * pop(ccDist);
                dPop(ccDist) = dPop(ccDist) + kDR * pop(ccReg) ... % regional -> distant
                    - deathCC(3) * pop(ccDist)... % distant CC mortality
                    - distDet; % detect distant CC -> treat
                
                
                % CC treated tracker
                ccTreated(d , h , a , 1) = sum(locDet);
                ccTreated(d , h , a , 2) = sum(regDet);
                ccTreated(d , h , a , 3) = sum(distDet);
                
                ccLocDet = ccLocDetInds(d , h , a , :);
                ccRegDet = ccRegDetInds(d , h , a , :);
                ccDistDet = ccDistDetInds(d , h , a , :);
                dPop(ccLocDet) = dPop(ccLocDet) + locDet;
                dPop(ccRegDet) = dPop(ccRegDet) + regDet;
                dPop(ccDistDet) = dPop(ccDistDet) + distDet;
                
                % CC incidence tracker
                ccInc(d , h , a) = ccInc(d , h , a) + sum(kCC_Cin3(a , h - 1) .* pop(cin3));
                
                % CC death tracker
                ccDeath(d , h , a) = ccDeath(d , h , a) + ...
                    sum(deathCC(1) * pop(ccLoc) ...
                    + deathCC(2) * pop(ccReg) + deathCC(3) * pop(ccDist));
            end
            
            % Detected CC deaths and transitions
            dPop(ccLocDet) = dPop(ccLocDet) - kRL * pop(ccLocDet) ...
                - deathCC_det(1) * pop(ccLocDet);
            dPop(ccRegDet) = dPop(ccRegDet) + kRL * pop(ccLocDet) ...
                - kDR * pop(ccRegDet) - deathCC_det(2) * pop(ccRegDet);
            dPop(ccDistDet) = dPop(ccDistDet) + kDR * pop(ccRegDet) ...
                - deathCC_det(3) * pop(ccDistDet);
            
            % Add detected CC deaths to CC death tracker
            ccDeath(d , h , a) = ccDeath(d , h , a) + ...
                sum(deathCC_det(1) * pop(ccLocDet) ...
                + deathCC_det(2) * pop(ccRegDet)...
                + deathCC_det(3) * pop(ccDistDet));
        end
    end
end

extraOut{1} = ccInc;
extraOut{2} = ccDeath;
extraOut{3} = ccTreated;