% HPV progression
% Simulates progression through HPV states as well as HPV treatment and
% hysterectomy.
% Accepts a population matrix as input and returns dPop, a vector of
% derivatives that describes the change in the population's subgroups due
% to HPV progression.
function[dPop , extraOut] = hpv(t , pop , immuneInds , infInds , cin1Inds , ...
    cin2Inds , cin3Inds , normalInds , ccInds , ccRegInds , ccDistInds , ...
    kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
    kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf  ,...
    rNormal_Inf , hpv_hivClear , c3c2Mults , hpvClearMult , ...
    c2c1Mults , fImm , kRL , kDR , muCC , disease , viral , age , hpvTypes , ...
    rImmuneHiv , vaccinated , hystOption)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sumall = @(x) sum(x(:));
hyst = 0;
if strcmp(hystOption , 'on')
    hyst = 1;
end
ccInc = zeros(disease , hpvTypes , age);
ccDeath = ccInc;
% constants
% see model notes for index values
% leep (effective treatment rate by leep)
rImmune = 0.024; % for HPV16, Johnson
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dPop = zeros(size(pop));
for d = 1 : disease
    c3c2Mult = c3c2Mults(1); % multiplier used for CIN2 -> CIN3 in HIV infecteds
    c2c1Mult = c2c1Mults(1); % multiplier
    %multiplier used for CIN1 -> CIN2 in HIV infecteds
    c1c2Mult = 1; % CIN2 -> CIN1 regression multiplier
    c2c3Mult = 1; % CIN3 -> CIN2 regression multiplier
    rHivHpvMult = 1; % for HIV negative
    hivPos = 0;
    deathCC = muCC(6 , :); % HIV negative
    rHiv = 1; % Multiplier for immunity clearance for HIV+
    rHivHpv_Clear = 1;
    if d > 2 && d < 7 % CD4 > 500 -> CD4 < 200
        c3c2Mult = c3c2Mults(d - 2); % CIN2 -> CIN3 multiplier
        c2c1Mult = c2c1Mults(d - 2); % CIN1 -> CIN2 multiplier
        c1c2Mult = hpv_hivClear(d - 2); % CIN2 -> CIN1 regression multiplier
        c2c3Mult = hpv_hivClear(d - 2); % CIN3 -> CIN2 regression multiplier
        rHivHpvMult = hpvClearMult(d - 2); % Regression multiplier
        rHivHpv_Clear = 0.1 .* hpv_hivClear(d - 2); % Infection clearance multiplier
        deathCC = muCC(d - 2 , :); % HIV+ CC death rate
        rHiv = rImmuneHiv(d - 2); % Multiplier for immunity clearance for HIV+
    elseif d == 10 % CD4 > 500 multipliers for HIV+ on ART
        c3c2Mult = c3c2Mults(1); % CIN2 -> CIN3 multiplier
        c2c1Mult = c2c1Mults(1); % CIN1 -> CIN2 multiplier
        c1c2Mult = hpv_hivClear(1); % CIN2 -> CIN1 regression multiplier
        c2c3Mult = hpv_hivClear(1); % CIN3 -> CIN2 regression multiplier
        rHivHpvMult = hpvClearMult(1); % Regression multiplier
        rHivHpv_Clear = 1; % Infection clearance multiplier
        deathCC = muCC(5 , :); % HIV+ CC death rate
        rHiv = rImmuneHiv(1); % Multiplier for immunity clearance for HIV+
    end
    for h = 2 : hpvTypes % infected onwards
        for a = 1 : age
            % Infected group
            immuneM = immuneInds(d , h , 1 , a , :);
            immuneF = immuneInds(d , h , 2 , a , :);
            infM = infInds(d , h , 1 , a , :);
            infF = infInds(d , h , 2 , a , :);
            cin1 = cin1Inds(d , h , a , :);
            cin2 = cin2Inds(d , h , a , :);
            cin3 = cin3Inds(d , h , a , :);
            normalM = normalInds(d , 1 , a , :);
            normalF = normalInds(d , 2 , a , :);
            % to immune from HPV infected, CIN1, CIN2, CIN3 (remove immunity for now)
            
            dPop(normalF) = dPop(normalF) + rImmune * rHiv * pop(immuneF) ... % immuneF -> normalF
                + rNormal_Inf(a , h - 1) * (1 - fImm(a)) * rHivHpvMult .* pop(infF); % infF -> normalF
            
            dPop(normalM) = dPop(normalM) + rImmune * pop(immuneM) ...% immuneM -> normalM
                + rNormal_Inf(a , h - 1) * rHivHpvMult .* pop(infM);  % infM -> normalM
            
            dPop(immuneF) = dPop(immuneF)...
                + rNormal_Inf(a , h - 1) * fImm(a) * rHivHpv_Clear .* pop(infF)... % infF -> immuneF
                - rImmune * rHiv * pop(immuneF); % immuneF -> normalF
            
            %                 dPop(immuneM) = dPop(immuneM) + rNormal_Inf(a) ...
            %                     * fImm(a) * rHivHpvMult * pop(infM)...; % infected -> immune
            %                     - rImmune * pop(immuneM); % immune -> normal
            
            % infected -> CIN1
            dPop(infF) = dPop(infF) ...
                + kInf_Cin1(a , h - 1) * rHivHpvMult * pop(cin1)... % CIN1 -> infF
                - (kCin1_Inf(a , h - 1) + ... % infF -> CIN1
                rNormal_Inf(a , h - 1)...
                * rHivHpv_Clear) .* pop(infF); % infF -> immuneF
            
            
            dPop(infM) = dPop(infM) - rNormal_Inf(a , h - 1) * rHivHpvMult * pop(infM); % regression to normal from infected males
            
            % Infection and CIN progression in females only
            % kCin_Inf(stage , hpvType , age group)
            % CIN1 group
            dPop(cin1) = dPop(cin1) + kCin1_Inf(a , h - 1) .* pop(infF)... % progression to CIN1 from infected females
                + rHivHpvMult * kCin1_Cin2(a , h - 1) * c1c2Mult .* pop(cin2)... % CIN2 -> CIN1
                - (kCin2_Cin1(a , h - 1) * c2c1Mult ... % CIN1 -> CIN2
                + rHivHpvMult * kInf_Cin1(a , h - 1)) .* pop(cin1); % CIN1 -> Infected
            
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
            ccLoc = ccInds(d , h , a , :);
            dPop(ccLoc) = dPop(ccLoc) + kCC_Cin3(a , h - 1) .* pop(cin3)... % CIN3 -> CC
                - kRL * pop(ccLoc)... % local -> regional
                - deathCC(1) * pop(ccLoc); % local CC mortality
            
            ccReg = ccRegInds(d , h , a , :);
            dPop(ccReg) = dPop(ccReg) + kRL * pop(ccLoc)...  % local -> regional
                - kDR * pop(ccReg)... % regional -> distant
                - deathCC(2) * pop(ccReg); % regional CC mortality
            
            ccDist = ccDistInds(d , h , a , :);
            dPop(ccDist) = dPop(ccDist) + kDR * pop(ccReg) ... % regional -> distant
                - deathCC(3) * pop(ccDist); % distant CC mortality
            
            % CC incidence tracker
            ccInc(d , h , a) = sum(kCC_Cin3(a , h - 1) .* pop(cin3));
            
            % CC death tracker
            ccDeath(d , h , a) = sum(deathCC(1) * pop(ccLoc) ...
                + deathCC(2) * pop(ccReg) + deathCC(3) * pop(ccDist));
        end
    end
end


% vaccine waning
% dPop(vaccinated) = -k_wane * pop(vaccinated);
% dPop(waned) = dPop(waned) + k_wane * pop(vaccinated);

extraOut{1} = ccInc;
extraOut{2} = ccDeath;
%extraOut{3} = ccTreatCost;