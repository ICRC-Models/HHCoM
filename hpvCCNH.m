% HPV Natural History
% Simulates progression through HPV states.
% Accepts a population matrix as input and returns dPop, a vector of
% derivatives that describes the change in the population's subgroups due
% to HPV progression.
function[dPop , extraOut] = hpvCCNH(t , pop , ...
    hpv_hivClear , rImmuneHiv , c3c2Mults , c2c1Mults , muCC , ...
    normalhpvVaxInds , immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , ...
    immunehpvNonVaxInds , infhpvNonVaxInds , cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ...
    ccLochpvVaxIndsFrom , ccReghpvVaxInds , ccDisthpvVaxInds , ...
    cin3hpvNonVaxIndsFrom , ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ...
    ccReghpvNonVaxInds , ccDisthpvNonVaxInds , cin1hpvVaxInds , ...
    cin2hpvVaxInds , cin3hpvVaxInds , cin1hpvNonVaxInds , ...
    cin2hpvNonVaxInds , cin3hpvNonVaxInds , kInf_Cin1 , kCin1_Cin2 , kCin2_Cin3 , ...
    kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , kCin1_Inf , rNormal_Inf , ...
    rImmune , fImm , kRL , kDR , maleHpvClearMult , disease , age , hpvVaxStates , ...
    hpvNonVaxStates , hpvTypeGroups)

%% Initialize dPop and output vectors
dPop = zeros(size(pop));
ccInc = zeros(disease , age , hpvTypeGroups);
% cin1Inc = ccInc;
% cin2Inc = ccInc;
% cin3Inc = ccInc;
ccDeath = ccInc;

%% Progress HPV disease states
% Set transition multipliers based on HIV disease state
for d = 1 : disease
    % No multipliers for HIV-negative; set CC-associated mortality
    rHivHpv_Clear = 1;
    rHiv = 1;
    c3c2Mult = c3c2Mults(1); % =1
    c2c1Mult = c2c1Mults(1); % =1
    c1c2Mult = 1;
    rHivHpvMult = 1;
    c2c3Mult = 1;
    deathCC = muCC(6 , :); % HIV-negative CC-associated mortality
    % Multipliers for HIV-positives with CD4>500 --> CD4<200; set increased CC-associated mortality
    if d > 3 && d < 8
        rHivHpv_Clear = hpv_hivClear(d - 3); % Infection clearance multiplier
        rHiv = rImmuneHiv(d - 3); % Multiplier for immunity clearance for HIV+
        c3c2Mult = c3c2Mults(d - 3); % CIN2 -> CIN3 progression multiplier
        c2c1Mult = c2c1Mults(d - 3); % CIN1 -> CIN2 progression multiplier
        c1c2Mult = hpv_hivClear(d - 3); % CIN2 -> CIN1 regression multiplier
        rHivHpvMult = hpv_hivClear(d - 3);%hpvClearMult(d - 2); % Regression multiplier, compounds c1c2Mult
        c2c3Mult = hpv_hivClear(d - 3); % CIN3 -> CIN2 regression multiplier
        deathCC = muCC(d - 3 , :); % HIV-positive CC-associated mortality     
    % Multipliers for HIV-positives on ART equivalent to those for CD4>500; ...
    % set CC-associated mortality equivalent to CD4>500
    elseif d == 8
        rHivHpv_Clear = hpv_hivClear(1); % Infection clearance multiplier
        rHiv = rImmuneHiv(1); % Multiplier for immunity clearance for HIV+
        c3c2Mult = c3c2Mults(1); % CIN2 -> CIN3 progression multiplier
        c2c1Mult = c2c1Mults(1); % CIN1 -> CIN2 progression multiplier
        c1c2Mult = hpv_hivClear(1); % CIN2 -> CIN1 regression multiplier
        rHivHpvMult = hpv_hivClear(1);%hpvClearMult(1); % Regression multiplier, compounds c1c2Mult
        c2c3Mult = hpv_hivClear(1); % CIN3 -> CIN2 regression multiplier
        deathCC = muCC(5 , :); % HIV-positive on ART CC-associated mortality
    end
    
    % Background hysterectomy (age dependent rate) % NOT UPDATED!!!!!!!!!!!!!!!!!
    %if hyst
    %    for a = 11 : age % can adjust age group range
    %        for s = [1:7 , 9:10]
    %            hystSusPop = hystSusInds(d , s , a , :);
    %            hystPop = hystInds(d , a , :);
    %            toHyst = pop(hystSusPop) .* OMEGA(a);        
    %            dPop(hystSusPop) = dPop(hystSusPop) - toHyst;
    %            dPop(hystPop) = dPop(hystPop) + toHyst;
    %        end
    %    end
    %end
    
    % HPV state transitions
    for a = 1 : age
        % prepare indices
        normalMVax = normalhpvVaxInds(d , 1 , a , :);
        normalFVax = normalhpvVaxInds(d , 2 , a , :);
        normalMNonVax = normalhpvNonVaxInds(d , 1 , a , :);
        normalFNonVax = normalhpvNonVaxInds(d , 2 , a , :);
        
        immuneFVax = immunehpvVaxInds(d , 2 , a , :);
        immuneFNonVax = immunehpvNonVaxInds(d , 2 , a , :);
        
        infMVax = infhpvVaxInds(d , 1 , a , :);
        infFVax = infhpvVaxInds(d , 2 , a , :);
        infMNonVax = infhpvNonVaxInds(d , 1 , a , :);
        infFNonVax = infhpvNonVaxInds(d , 2 , a , :);
        
        cin1Vax = cin1hpvVaxInds(d , a , :);
        cin2Vax = cin2hpvVaxInds(d , a , :);
        cin3Vax = cin3hpvVaxInds(d , a , :);
        cin1NonVax = cin1hpvNonVaxInds(d , a , :);
        cin2NonVax = cin2hpvNonVaxInds(d , a , :);
        cin3NonVax = cin3hpvNonVaxInds(d , a , :);
        
        % Adjust compartments
        % normal, immune, and infected transitions
        dPop(normalFVax) = dPop(normalFVax) + fImm(a) * rImmune * rHiv * pop(immuneFVax)... % if fImm(a)=1, immuneF -> normalF 
            + (1 - fImm(a)) * rNormal_Inf(a , 1) * rHivHpv_Clear .* pop(infFVax); % if fImm(a)=0, infF -> normalF
        dPop(normalMVax) = dPop(normalMVax) + maleHpvClearMult * rNormal_Inf(a , 1) * rHivHpv_Clear .* pop(infMVax);  % infM -> normalM
        
        dPop(normalFNonVax) = dPop(normalFNonVax) + fImm(a) * rImmune * rHiv * pop(immuneFNonVax)... % if fImm(a)=1, immuneF -> normalF 
            + (1 - fImm(a)) * rNormal_Inf(a , 2) * rHivHpv_Clear .* pop(infFNonVax); % if fImm(a)=0, infF -> normalF
        dPop(normalMNonVax) = dPop(normalMNonVax) + maleHpvClearMult * rNormal_Inf(a , 2) * rHivHpv_Clear .* pop(infMNonVax);  % infM -> normalM
            

        dPop(immuneFVax) = dPop(immuneFVax)...
            + fImm(a) * rNormal_Inf(a , 1) * rHivHpv_Clear .* pop(infFVax)... % if fImm(a)=1, infF -> immuneF
            - fImm(a) * rImmune * rHiv * pop(immuneFVax); % if fImm(a)=1, immuneF -> normalF
        
        dPop(immuneFNonVax) = dPop(immuneFNonVax)...
            + fImm(a) * rNormal_Inf(a , 2) * rHivHpv_Clear .* pop(infFNonVax)... % if fImm(a)=1, infF -> immuneF
            - fImm(a) * rImmune * rHiv * pop(immuneFNonVax); % if fImm(a)=1, immuneF -> normalF
        
        
        % infected -> CIN1
        dPop(infFVax) = dPop(infFVax) ...
            + kInf_Cin1(a , 1) * pop(cin1Vax)... % CIN1 -> infF
            - (kCin1_Inf(a , 1) + ... % infF -> CIN1
            rNormal_Inf(a , 1) * rHivHpv_Clear) .* pop(infFVax); % infF -> immuneF or normalF
        dPop(infMVax) = dPop(infMVax) - maleHpvClearMult * rNormal_Inf(a , 1) * rHivHpv_Clear * pop(infMVax); % regression to normal from infected males
        
        dPop(infFNonVax) = dPop(infFNonVax) ...
            + kInf_Cin1(a , 2) * pop(cin1NonVax)... % CIN1 -> infF
            - (kCin1_Inf(a , 2) + ... % infF -> CIN1
            rNormal_Inf(a , 2) * rHivHpv_Clear) .* pop(infFNonVax); % infF -> immuneF or normalF
        dPop(infMNonVax) = dPop(infMNonVax) - maleHpvClearMult * rNormal_Inf(a , 2) * rHivHpv_Clear * pop(infMNonVax); % regression to normal from infected males
        

        % Infection and CIN progression in females only
        % CIN1 group
        dPop(cin1Vax) = dPop(cin1Vax) + kCin1_Inf(a , 1) .* pop(infFVax)... % progression to CIN1 from infected females
            + rHivHpvMult * kCin1_Cin2(a , 1) * c1c2Mult .* pop(cin2Vax)... % CIN2 -> CIN1
            - (kCin2_Cin1(a , 1) * c2c1Mult ... % CIN1 -> CIN2
            + kInf_Cin1(a , 1)) .* pop(cin1Vax); % CIN1 -> Infected
        
        dPop(cin1NonVax) = dPop(cin1NonVax) + kCin1_Inf(a , 2) .* pop(infFNonVax)... % progression to CIN1 from infected females
            + rHivHpvMult * kCin1_Cin2(a , 2) * c1c2Mult .* pop(cin2NonVax)... % CIN2 -> CIN1
            - (kCin2_Cin1(a , 2) * c2c1Mult ... % CIN1 -> CIN2
            + kInf_Cin1(a , 2)) .* pop(cin1NonVax); % CIN1 -> Infected
        

        % CIN2 group
        dPop(cin2Vax) = dPop(cin2Vax) ...
            + kCin2_Cin1(a , 1) * c2c1Mult * pop(cin1Vax) ... % CIN1 -> CIN2
            + kCin2_Cin3(a , 1) * c2c3Mult * pop(cin3Vax) ... % CIN3 -> CIN2
            - (kCin3_Cin2(a , 1) * c3c2Mult... % CIN2 -> CIN3
            + kCin1_Cin2(a , 1) * rHivHpvMult * c1c2Mult) .* pop(cin2Vax); % CIN2 -> CIN1
        
        dPop(cin2NonVax) = dPop(cin2NonVax) ...
            + kCin2_Cin1(a , 2) * c2c1Mult * pop(cin1NonVax) ... % CIN1 -> CIN2
            + kCin2_Cin3(a , 2) * c2c3Mult * pop(cin3NonVax) ... % CIN3 -> CIN2
            - (kCin3_Cin2(a , 2) * c3c2Mult... % CIN2 -> CIN3
            + kCin1_Cin2(a , 2) * rHivHpvMult * c1c2Mult) .* pop(cin2NonVax) ; % CIN2 -> CIN1

        
        % CIN3 group
        dPop(cin3Vax) = dPop(cin3Vax) + kCin3_Cin2(a , 1) * c3c2Mult .* pop(cin2Vax)... %CIN2 -> CIN3
            - (kCC_Cin3(a , 1)... % CIN3 -> CC
            + kCin2_Cin3(a , 1) * c2c3Mult)... % CIN3 -> CIN2
            .* pop(cin3Vax);
        
        dPop(cin3NonVax) = dPop(cin3NonVax) + kCin3_Cin2(a , 2) * c3c2Mult .* pop(cin2NonVax)... %CIN2 -> CIN3
            - (kCC_Cin3(a , 2)... % CIN3 -> CC
            + kCin2_Cin3(a , 2) * c2c3Mult)... % CIN3 -> CIN2
            .* pop(cin3NonVax);

        
        % Cervical cancer        
        for s = 1 : hpvNonVaxStates
            % prepare indices
            cin3VaxFrom = cin3hpvVaxIndsFrom(d , s , a , :);
            ccLocVaxTo = ccLochpvVaxIndsTo(d , s , a , :);
            ccLocVaxFrom = ccLochpvVaxIndsFrom(d , s , a , :);
            ccRegVax = ccReghpvVaxInds(d , s , a , :);
            ccDistVax = ccDisthpvVaxInds(d , s , a , :);

            % local
            dPop(ccLocVaxTo) = dPop(ccLocVaxTo) + kCC_Cin3(a , 1) .* pop(cin3VaxFrom); % CIN3 -> CC
            if (s ~= 6) % track CC incidence; do not double count if person already has cancer from a different HPV type
                ccInc(d , a , 1) = ccInc(d , a , 1) + sum(kCC_Cin3(a , 1) .* pop(cin3VaxFrom));
                %cin1Inc(d , a , 1) = cin1Inc(d , a , 1) + sum(kCin1_Inf(a , 1) .* pop(infF)); % !!!!!!infF indices out-dated
                %cin2Inc(d , a , 1) = cin2Inc(d , a , 1) + sum(kCin2_Cin1(a , 1) * c2c1Mult * pop(cin1)); % !!!!!!cin1 indices out-dated
                %cin3Inc(d , a , 1) = cin3Inc(d , a , 1) + sum(kCin3_Cin2(a , 1) * c3c2Mult .* pop(cin2)); % !!!!!!cin2 indices out-dated
            end
            dPop(ccLocVaxFrom) = dPop(ccLocVaxFrom) - kRL * pop(ccLocVaxFrom)... % local -> regional
                - deathCC(1) * pop(ccLocVaxFrom); % local CC mortality
                    
            % regional
            dPop(ccRegVax) = dPop(ccRegVax) + kRL * pop(ccLocVaxFrom)...  % local -> regional
                - kDR * pop(ccRegVax)... % regional -> distant
                - deathCC(2) * pop(ccRegVax); % regional CC mortality
        
            % distant
            dPop(ccDistVax) = dPop(ccDistVax) + kDR * pop(ccRegVax) ... % regional -> distant
                - deathCC(3) * pop(ccDistVax); % distant CC mortality
            
            % CC death tracker
            ccDeath(d , a , 1) = ccDeath(d , a , 1) + ...
                sum(deathCC(1) * pop(ccLocVaxFrom) ...
                + deathCC(2) * pop(ccRegVax) + deathCC(3) * pop(ccDistVax));
        end
       
        for h = 1 : hpvVaxStates
            cin3NonVaxFrom = cin3hpvNonVaxIndsFrom(d , h , a , :);
            ccLocNonVaxTo = ccLochpvNonVaxIndsTo(d , h , a , :);
            ccLocNonVaxFrom = ccLochpvNonVaxIndsFrom(d , h , a , :);
            ccRegNonVax = ccReghpvNonVaxInds(d , h , a , :);
            ccDistNonVax = ccDisthpvNonVaxInds(d , h , a , :);
            
            % local
            dPop(ccLocNonVaxTo) = dPop(ccLocNonVaxTo) + kCC_Cin3(a , 2) .* pop(cin3NonVaxFrom); % CIN3 -> CC
            if (h ~= 6) % track CC incidence; do not double count if person already has cancer from a different HPV type
                ccInc(d , a , 2) = ccInc(d , a , 2) + sum(kCC_Cin3(a , 2) .* pop(cin3NonVaxFrom));
                %cin1Inc(d , a , 2) = cin1Inc(d , a , 2) + sum(kCin1_Inf(a , 2) .* pop(infF)); % !!!!!!infF indices out-dated
                %cin2Inc(d , a , 2) = cin2Inc(d , a , 2) + sum(kCin2_Cin1(a , 2) * c2c1Mult * pop(cin1)); % !!!!!!cin1 indices out-dated
                %cin3Inc(d , a , 2) = cin3Inc(d , a , 2) + sum(kCin3_Cin2(a , 2) * c3c2Mult .* pop(cin2)); % !!!!!!cin2 indices out-dated
            
                dPop(ccLocNonVaxFrom) = dPop(ccLocNonVaxFrom) - kRL * pop(ccLocNonVaxFrom)... % local -> regional
                    - deathCC(1) * pop(ccLocNonVaxFrom); % local CC mortality

                % regional
                dPop(ccRegNonVax) = dPop(ccRegNonVax) + kRL * pop(ccLocNonVaxFrom)...  % local -> regional
                    - kDR * pop(ccRegNonVax)... % regional -> distant
                    - deathCC(2) * pop(ccRegNonVax); % regional CC mortality

                % distant
                dPop(ccDistNonVax) = dPop(ccDistNonVax) + kDR * pop(ccRegNonVax) ... % regional -> distant
                    - deathCC(3) * pop(ccDistNonVax); % distant CC mortality
                
                % CC death tracker
                ccDeath(d , a , 2) = ccDeath(d , a , 2) + ...
                    sum(deathCC(2) * pop(ccLocNonVaxFrom) ...
                    + deathCC(2) * pop(ccRegNonVax) + deathCC(3) * pop(ccDistNonVax));
            end
        end
    end   
end

%% Save outputs and convert dPop to a column vector for output to ODE solver
extraOut{1} = ccInc;
extraOut{2} = ccDeath;

dPop = sparse(dPop);
