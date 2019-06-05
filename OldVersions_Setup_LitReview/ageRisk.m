function [dPop , extraOut] = ageRisk(t , pop , disease , hpvTypes , ...
    viral , hpvStates , gender , ...
    periods , age , risk ,riskDist , ageInd , riskInd)
sumall = @(x) sum(x(:));

dPop = zeros(size(pop));

for g = 1 : gender
    for a = 2 : age
        aPrev = ageInd(g , a - 1 , :);
        aCurr = ageInd(g , a , :);
        
        r1 = riskInd(g , a - 1 , 1 , :);
        r2 = riskInd(g , a - 1 , 2 , :);
        r3 = riskInd(g , a - 1 , 3 , :);
        r1To = riskInd(g , a , 1 , :);
        r2To = riskInd(g , a , 2 , :);
        r3To = riskInd(g , a , 3 , :);
        
        popR1Tot = sumall(pop(r1));
        popR2Tot = sumall(pop(r2));
        popR3Tot = sumall(pop(r3));
        
        % get prospective risk distribution if staying in same risk group when
        % aging
        agedOut = 1/5 .* sumall(pop(aPrev)); % age 1/5th of previous age group
        agedProsp = agedOut + 4/5 .* sumall(pop(aCurr)); % age 1/5th of previous age group into current age group
        riskTarget = agedProsp .* riskDist(a , : , g);
        riskNeed = riskTarget - 4/5 .* [sumall(pop(r1To)) , sumall(pop(r2To)) , sumall(pop(r3To))]; % numbers needed to fill risk groups
        riskAvail = 1/5 .* [popR1Tot , popR2Tot , popR3Tot];
        riskDiff = riskNeed - riskAvail; % difference between numbers needed and available for each risk group

        riskFrac1 = 0;
        riskFrac2 = 0;
        riskFrac3 = 0;

        

        if riskDiff(3) > 0 % if risk 3 deficient
            % start with moving from risk 2 to risk 3
            if riskAvail(2) > 0
                riskFrac2 = min(min(riskDiff(3) , riskAvail(2)) / popR2Tot , 1);
                dPop(r2To) = dPop(r2To) - pop(r2) .* riskFrac2;
                dPop(r3To) = dPop(r3To) + pop(r2) .* riskFrac2;
            end
            % if needed, move from risk 1 to risk 3
            if riskDiff(3) / riskAvail(2) > 1
                riskFrac1 = ...
                    min(min(riskAvail(1) , (riskDiff(3) - riskAvail(2))) / popR1Tot , 1);
                dPop(r1To) = dPop(r1To) - pop(r1) .* riskFrac1;
                dPop(r3To) = dPop(r3To) + pop(r1) .* riskFrac1;
            end
            riskAvail(1) = riskAvail(1) - sum(pop(r1) .* riskFrac1);
            riskAvail(2) = riskAvail(2) - sum(pop(r2) .* riskFrac2);
        end

        if riskDiff(2) > 0 % if risk 2 deficient
            % start with moving from risk 3 to risk 2
            if riskAvail(3) > 0
                riskFrac3 = min(min(riskDiff(2) , riskAvail(3)) / popR3Tot , 1);
                dPop(r3To) = dPop(r3To) - pop(r3) .* riskFrac3;
                dPop(r2To) = dPop(r2To) + pop(r3) .* riskFrac3;
            end
            % if needed, move from risk 1 to risk 2
            if riskDiff(2) / riskAvail(3) > 1
                riskFrac1 =...
                    min(min((riskDiff(2) - riskAvail(3)) , riskAvail(1)) / popR1Tot , 1);
                dPop(r1To) = dPop(r1To) - pop(r1) .* riskFrac1;
                dPop(r2To) = dPop(r2To) + pop(r1) .* riskFrac1;
            end
            riskAvail(1) = riskAvail(1) - sum(pop(r1) .* riskFrac1);
            riskAvail(3) = riskAvail(3) - sum(pop(r3) .* riskFrac3);
        end

        if riskDiff(1) > 0 % if risk 1 deficient
            % start with moving from risk 2 to risk 1
            if riskAvail(2) > 0
                riskFrac2 = min(min(riskDiff(1), riskAvail(2)) / popR2Tot , 1);
                dPop(r2To) = dPop(r2To) - pop(r2) .* riskFrac2;
                dPop(r1To) = dPop(r1To) + pop(r2) .* riskFrac2;
            end
            % if needed, move from risk 3 to risk 1
            if riskDiff(1) / riskAvail(2) > 1
                riskFrac3 = ...
                    min(min((riskDiff(1) - riskAvail(2)) , riskAvail(3)) / popR3Tot , 1);
                dPop(r3To) = dPop(r3To) - pop(r3) .* riskFrac3;
                dPop(r1To) = dPop(r1To) + pop(r3) .* riskFrac3;
            end
            riskAvail(2) = riskAvail(2) - sum(pop(r2) .* riskFrac2);
            riskAvail(3) = riskAvail(3) - sum(pop(r3) .* riskFrac3);
        end
        
        dPop(r1To) = dPop(r1To) + 1/5 .* pop(r1);
        dPop(r2To) = dPop(r2To) + 1/5 .* pop(r2);
        dPop(r3To) = dPop(r3To) + 1/5 .* pop(r3); 
        
        dPop(r1) = dPop(r1) - 1/5 .* pop(r1);
        dPop(r2) = dPop(r2) - 1/5 .* pop(r2);
        dPop(r3) = dPop(r3) - 1/5 .* pop(r3);
    end
    % age last age group
%     dPop(r1To) = dPop(r1To) - 1/5 .* pop(r1To);
%     dPop(r2To) = dPop(r2To) - 1/5 .* pop(r2To);
%     dPop(r3To) = dPop(r3To) - 1/5 .* pop(r3To);
end


