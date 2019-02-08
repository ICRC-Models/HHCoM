function makeMat
close all; clear all; clc
% loadUp(6);
%% Initialize pop vector
disp('Building matrices')
paramDir = [pwd , '\Params\'];
load([paramDir, 'general'])
load([paramDir , 'popData'])
savedir = [pwd , '\Params']; 
pop = spalloc(prod(dim) , 1 , prod(dim));
at = @(x , y) sort(prod(dim)*(y-1) + x); 
%% aging and risk assortment
disp('Building aging matrix')
% aging matrix
ageIn = spalloc(numel(pop) , numel(pop) , numel(pop));
for g = 1 : gender
    for a = 1 : age - 1
        for r = 1 : risk
            fromAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
                1 : hpvStates , 1 : periods , g , a, r));
            toAge = toInd (allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
                1 : hpvStates , 1 : periods , g , a + 1 , r));
            ageIn(at(toAge , fromAge)) = 1/5;
        end
    end
end

% risk sorting matrix
disp('Building risk sorting matrix')
riskSorter = speye(numel(pop) , numel(pop));
riskDist(1 , : , :) = riskDistM;
riskDist(2 , : , :) = riskDistF;
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            for rr = 1 : risk
                fromRisk = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
                    1 : hpvStates , 1 : periods , g , a , r));
                toRisk = toInd (allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
                    1 : hpvStates , 1 : periods , g , a , rr));
                riskSorter(at(toRisk , fromRisk)) = riskDist(g , a , rr);
            end
        end
    end
end
save(fullfile(savedir , 'riskSorter') , 'riskSorter');

disp('Combining age and risk sorting matrices')
ageRiskSorter = riskSorter * ageIn;
save(fullfile(savedir , 'ageRiskSorter') , 'ageRiskSorter')
ageOut = spalloc(numel(pop) , numel(pop) , numel(pop));
for g = 1 : gender
    for a = 1 : age
        for r = 1 : risk
            fromAge = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , ...
                1 : hpvStates , 1 : periods , g , a, r));
            ageOut(at(fromAge , fromAge)) = - 1/5;
        end
    end
end
save(fullfile(savedir ,'ageOut') , 'ageOut')
ager = riskSorter * ageIn + ageOut;
save(fullfile(savedir ,'ager') , 'ager')
disp('Finished building age and risk matrices')
%% hiv
% produces hivTrans, hivDeathMat, artMat, and prepMat

% load([paramDir , 'HIVParams'])
% disp('Building HIV CD4 progression matrix')
% % cd4 progression
% % Time dependent components: ART uptake/dropout (artIn, artOut) , PrEP
% % uptake/dropout (toPrep , prepOut)
% 
% % PrEP uptake/dropout intially 0
% toPrep = 0;
% prepOut = 0;
% artIn = zeros(2 , 1);
% artOut = zeros(2 , 1);
% 
% % Tracks group characteristics at time of ART commencement. Set to 0 here just to fill intial matrix. 
% import java.util.LinkedList
% artDistList = LinkedList();
% artDist = zeros(disease , viral , gender , age , risk); % initial distribution of inidividuals on ART = 0
% 
% 
% hivTrans = spalloc(numel(pop) , numel(pop), numel(pop));
% hivDeathMat = spalloc(numel(pop) , numel(pop), numel(pop));
% artMat = spalloc(numel(pop) , numel(pop), numel(pop));
% prepMat = spalloc(numel(pop) , numel(pop) , numel(pop));
% 
% hivPositiveArt = toInd(allcomb(10 , 6 , 1: hpvTypes , 1 : hpvStates, ...
%     1 : periods , 1 : gender , 1 : age , 1 : risk));
% 
% % Dropout from ART
% above200 = artDist(1 : 5 , 1 : viral , 1 : gender , 1 : age , 1 : risk);
% below200 = artDist(6 , 1 : viral , 1 : gender , 1 : age , 1 : risk);
% 
% % ART dropout weighted by proportion of individuals with CD4 above/below 200 receiving ART
% hivTrans(at(hivPositiveArt , hivPositiveArt)) = ...
%     - (artOut(1) * sum(above200(:)) + artOut(2) * sum(below200(:)));
% 
% hivNegative = toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates, ...
%     1 : periods , 1 : gender , 1 : age , 1 : risk));
% hivNegativePrep = toInd(allcomb(9 , 1 , 1: hpvTypes , 1 : hpvStates, ...
%     1 : periods , 1 : gender , 1 : age , 1 : risk));
% 
% prepMat(at(hivNegativePrep , hivNegative)) = toPrep; % update at specific time steps
% prepMat(at(hivNegative , hivNegativePrep)) = prepOut; % update at specific time steps
% 
% for g = 1 : gender
%     for a = 1 : age
%         for r = 1 : risk
%             for v = 1 : 5
%                 acuteInf = toInd(allcomb(2 , v , 1 : hpvTypes , 1 : hpvStates , ...
%                     1 : periods , 1 : gender , a , 1 : risk));
%                 hivDeathMat(at(acuteInf , acuteInf)) = - muHIV(a , 2); % HIV acute infection mortality
%                 
%                 % build HIV death tally later!
%                 
%                 for d = 3 : 6
%                     cd4Curr = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates, ...
%                         1 : periods , g , 1 : age , 1 : risk));
%                     cd4Prev = toInd(allcomb(d - 1 , v , 1 : hpvTypes , 1 : hpvStates, ...
%                         1 : periods , g , 1 : age , 1 : risk));
%                     artDropOut = artOut(1); % dropout rate when Cd4 > 200
%                     toArt = artIn(1); % art initiation rate when CD4 > 200
%                     if d == 6
%                         artDropOut = artOut(2); % dropout rate when CD4 <= 200
%                         toArt = artIn(2); % art initiation rate when CD4 <= 200
%                     end
%                     
%                     % CD4 progression
%                     hivTrans(at(cd4Curr , cd4Prev)) = kCD4(g , v , d - 2); % previous -> current
%                     
%                     % HIV+ on ART -> dropout. Note: transition out of ART accounted for in weighted dropout.
%                     artMat(at(cd4Curr , hivPositiveArt)) = artDropOut * artDist(d , v , g , a , r); 
%                     hivTrans(at(cd4Curr , cd4Curr)) = -(kCD4(g , v , d - 1) * (d ~= 6)); % progression to next CD4 state
%                     
%                     % deaths due to HIV by CD4
%                     hivDeathMat(at(cd4Curr , cd4Curr)) = - muHIV(a , d);
%                     
%                     % HIV-positive going on ART
%                     artMat(at(cd4Curr , cd4Curr)) = - toArt;
%                     artMat(at(hivPositiveArt , cd4Curr)) = toArt;   
%                 end
%             end
%             
%             % HIV-negative circumcised
%             if g == 1
%                 hivNegCirc = toInd(allcomb(7 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , 1 , a , r));
%                 hivNegCircPrep = toInd(allcomb(8 , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , 1 , a , r));
%             end
%             
%             prepMat(at(hivNegCirc , hivNegCircPrep)) = prepOut; % circumcised and on PrEP -> circumcised
%             prepMat(at(hivNegCirc , hivNegCirc)) = - toPrep; % circumcised -> PrEP
%             
%             prepMat(at(hivNegCircPrep , hivNegCirc)) = toPrep; % circumcised -> circumcised and on PrEP
%             prepMat(at(hivNegCircPrep , hivNegCircPrep)) = - prepOut; % circumcised and on PrEP -> circumcised
%             
%             prepMat(at(hivNegativePrep , hivNegative)) = toPrep; % HIV negative -> HIV negative on PrEP
%             prepMat(at(hivNegativePrep , hivNegativePrep)) = - prepOut; % HIV negative on PrEP -> HIV negative
%         end
%     end
% end
% 
% save(fullfile(savedir , 'hivTrans') , 'hivTrans')
% save(fullfile(savedir , 'hivDeathMat') , 'hivDeathMat')
% save(fullfile(savedir , 'artMat') , 'artMat')
% save(fullfile(savedir , 'prepMat') , 'prepMat')
% disp('Finished building HIV matrices.')
%% Viral load progression (by CD4 count)
disp('Building viral load progression matrix')
load([paramDir , 'HIVParams'])
vlAdvancer = spalloc(numel(pop) , numel(pop) , numel(pop));

for g = 1 : gender
    for d = 2 : 6
        % get vlAcute index
        vlAcute = toInd(allcomb(d , 1 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , 1 : age , 1 : risk));
        
        % for v = 1 %VL acute
        vlAdvancer(at(vlAcute , vlAcute)) = - kVl(g , d - 1 , 1);
        
        % for vl = 2 : 4 % VL < 1,000 ->  1,000 - 10,000 -> 10,000 - 50,000
        for v = 2 : 4
            % get vl index
            vlCurr = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , 1 : age , 1 : risk));
            
            vlPrev = toInd(allcomb(d , (v - 1) , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , 1 : age , 1 : risk));
            
            % calculate vl transitions
            vlAdvancer(at(vlCurr , vlPrev)) = kVl(g , d - 1 , v - 1); % prev -> curr
            vlAdvancer(at(vlCurr , vlCurr)) = - kVl(g , d - 1 , v);   % curr -> (next)
        end
        
        % for v = 5, VL > 50,000
        % get indices
        vl_50k = toInd(allcomb(d , 5 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , 1 : age , 1 : risk));
        
        vl_10kto50k = toInd(allcomb(d , 4 , 1 : hpvTypes , 1 : hpvStates , 1 : periods , g , 1 : age , 1 : risk));
        
        % calculate transition to last vl
        vlAdvancer(at(vl_50k , vl_10kto50k)) = kVl(g , d - 1 , 4);
    end
end

save(fullfile(savedir ,'vlAdvancer') , 'vlAdvancer')
disp('Finished building viral load progression matrix')

%% hpv
% rImmune = 0.024; % for HPV16, Johnson
% %rNormal_Inf = rNormal_Inf * 0.5; % adjust rate of inf -> immunity downward
% % fImm(1 : 3) = 1; %0.1;
% % fImm(4 : age) = 0.58; % (0.48; 0.27 , 0.69) fraction fully protected by immunity based on RR of natural immunity (Beachler, 2017)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dPop = zeros(size(pop));
% for d = 1 : disease
%     c3c2Mult = c3c2Mults(1); % multiplier used for CIN2 -> CIN3 in HIV infecteds
%     c2c1Mult = c2c1Mults(1); % multiplier used for CIN1 -> CIN2 in HIV infecteds
%     c1c2Mult = 1; % CIN2 -> CIN1 regression multiplier
%     c2c3Mult = 1; % CIN3 -> CIN2 regression multiplier
%     rHivHpvMult = 1; % for HIV negative
%     if d > 2 && d < 7 % CD4 > 500 -> CD4 < 200
%         c3c2Mult = c3c2Mults(d - 2); % CIN2 -> CIN3 multiplier
%         c2c1Mult = c2c1Mults(d - 2); % CIN1 -> CIN2 multiplier
%         c1c2Mult = hpv_hivClear(d - 2); % CIN2 -> CIN1 regression multiplier
%         c2c3Mult = hpv_hivClear(d - 2); % CIN3 -> CIN2 regression multiplier
%         rHivHpvMult = hpv_hivClear(d - 2); % Infection clearance multiplier
%     end
%     for v = 1 : viral       
%         for h = 2 : hpvTypes % infected onwards
%             for a = 1 : age
%                 % Infected group
%                 immuneM = immuneInds(d , v , h , 1 , a , :);
%                 immuneF = immuneInds(d , v , h , 2 , a , :);
%                 infM = infInds(d , v , h , 1 , a , :);
%                 infF = infInds(d , v  , h , 2 , a , :);
%                 cin1 = cin1Inds(d , v , h , a , :);
%                 cin2 = cin2Inds(d , v , h , a , :);
%                 cin3 = cin3Inds(d , v , h , a , :);
%                 normalM = normalInds(d , v , 1 , a , :);
%                 normalF = normalInds(d , v , 2 , a , :);
%                 % to immune from HPV infected, CIN1, CIN2, CIN3 (remove immunity for now)
%                 hpvMat(at(normalF , immuneF)) = rImmune;
%                 hpvMat(at(normalF , cin2)) = kNormal_Cin2(a , h - 1) * (1 - fImm(a));
%                 hpvMat(at(normalF , cin1)) = kNormal_Cin1(a , h - 1) * (1 - fImm(a));
%                 hpvMat(at(normalF , infF)) = rNormal_Inf(a , h - 1) * (1 - fImm(a)) * rHivHpvMult;
%                 
%                 
%                 dPop(normalM) = dPop(normalM) + rImmune * pop(immuneM) ...% Immune -> Normal
%                     + rNormal_Inf(a , h - 1) * (1 - fImm(a)) * rHivHpvMult .* pop(infM);  % Infected -> Normal
%                 
%                 dPop(immuneF) = dPop(immuneF) + kNormal_Cin2(a , h - 1) * fImm(a) .* pop(cin2)... %CIN2 -> immune
%                     + kNormal_Cin1(a , h - 1) * fImm(a) .* pop(cin1)... % CIN1 -> immune
%                     + rNormal_Inf(a , h - 1) * fImm(a) * rHivHpvMult .* pop(infF)... % Inf -> immune
%                     - rImmune * pop(immuneF); % immune -> normal
%                 
%                 dPop(immuneM) = dPop(immuneM) + rNormal_Inf(a , h - 1) ...
%                     * fImm(a) * rHivHpvMult * pop(infM)...; % infected -> immune
%                     - rImmune * pop(immuneM); % immune -> normal
%                 
%                 % infected -> CIN1
%                 dPop(infF) = dPop(infF) ...
%                     + kInf_Cin2(a , h - 1) * pop(cin2)... % CIN2 -> infF
%                     + kInf_Cin1(a , h -1) * pop(cin1)... % CIN1 -> infF
%                     - (kCin1_Inf(a , h - 1) ... % progression to CIN1 from infected females
%                     + kCin2_Inf(a , h - 1) ... % progression to CIN2 from infected females (fast progressors)
%                     + rNormal_Inf(a , h - 1) * rHivHpvMult) .* pop(infF); % regression to immune from infected females
%                 
%                 dPop(infM) = dPop(infM) - rNormal_Inf(a , h - 1) * rHivHpvMult * pop(infM); % regression to immune from infected males              
%                 
%                 % Infection and CIN progression in females only
%                 % kCin_Inf(stage , hpvType , age group)
%                 % CIN1 group
%                 dPop(cin1) = dPop(cin1) + kCin1_Inf(a , h - 1) .* pop(infF)... % infected -> CIN1
%                     + kCin1_Cin2(a , h - 1) * c1c2Mult .* pop(cin2)... % CIN2 -> CIN1
%                     + kCin1_Cin3(a , h - 1) * c2c3Mult .* pop(cin3) ... % CIN3 -> CIN1 (Fast regressors)
%                     - (kCin2_Cin1(a , h - 1) * c2c1Mult + kCin3_Cin1(a , h - 1)... % CIN1 -> CIN2, CIN1 -> CIN3 (Fast Progressors)
%                     + kInf_Cin1(a , h - 1) + kNormal_Cin1(a , h - 1)) .* pop(cin1); % CIN1 -> Infected , CIN1 -> Normal/immuneF
%                 
%                 % CIN2 group
%                 dPop(cin2) = dPop(cin2) + kCin2_Inf(a , h - 1) .* pop(infF)... % Infected -> CIN2 (Fast progressors)
%                     + kCin2_Cin1(a , h - 1) * c2c1Mult * pop(cin1) ... % CIN1 -> CIN2
%                     + kCin2_Cin3(a , h - 1) * c2c3Mult * pop(cin3) ... % CIN3 -> CIN2
%                     - (kCin3_Cin2(a , h - 1) * c3c2Mult... % CIN2 -> CIN3
%                     + kInf_Cin2(a , h - 1) + kNormal_Cin2(a , h - 1) + kCin1_Cin2(a , h - 1)) .* pop(cin2) ; % CIN2 -> Inf, CIN2 -> Normal/immuneF , CIN2 -> CIN1
%                 
%                 % CIN3 group
%                 dPop(cin3) = dPop(cin3) + kCin3_Cin2(a , h - 1) * c3c2Mult .* pop(cin2)... %CIN2 -> CIN3
%                     + kCin3_Cin1(a , h - 1) * pop(cin1) +... % CIN1 -> CIN3
%                     - (kCC_Cin3(a , h - 1)... % CIN3 -> CC
%                     + kCin1_Cin3(a , h - 1) * c2c3Mult + kCin2_Cin3(a , h - 1))...
%                     * c2c3Mult .* pop(cin3); % CIN3 -> CIN1 (Fast regressors), CIN3 -> CIN2
%                 
%                 % CC group
%                 cc = ccInds(d , v , h , a , :);
%                 dPop(cc) = dPop(cc) + kCC_Cin3(a , h - 1) .* pop(cin3); % CIN3 -> CC
%                 
%                 % CC incidence tracker
%                 ccInc(d , v , h , a) = sumall(dPop(cc));
%             end
%         end
%     end
% end
%% hpv screening and treatment
% load([paramDir , 'hpvData'])
% load([paramDir , 'hpvIndices'])
% load([paramDir , 'hpvTreatIndices'])
% disp('Building screening and treatment matrix for pre-2015')
% screenTreater = spalloc(numel(pop) , numel(pop) , numel(pop));
% screenTreater2015 = spalloc(numel(pop) , numel(pop) , numel(pop));
% for v = 1 : viral
%     for d = 1 : disease
%         fScreen = screenFreq(: , 2);
%         if d >= 2 && d <= 6 || d == 10 % Infected or low CD4 or HIV+ on ART
%             fScreen = screenFreq(: , 1);
%         end
%         hpvEff = screenCover * fScreen(1) * hpvSens; % effective screen and detect rate by hpv test
%         cytoEff = screenCover * fScreen(2) * cytoSens; % effective screen and detect rate by cytology
%         for h = 2 : hpvTypes % infected onwards
%             for a = 5 : age
%                 %                 cin1 = cin1Inds(d , v , h , a , :);
%                 cin2 = cin2Inds(d , v , h , a , :);
%                 cin3 = cin3Inds(d , v , h , a , :);
%                 normal = normalInds(d , v , 2 , a , :);
%                 
%                 % CIN2 group detection and treatment
%                 screenTreater(at(cin2 , cin2)) = ...
%                     - (hpvEff(1) * cytoEff(1) * leep); % detection with hpv test and cytology, and treatment with leep
%                 screenTreater(at(normal , cin2)) = hpvEff(1) * cytoEff(1) * leep;
%                 
%                 % CIN3 group
%                 %                 if year <= 2015 % screening before 2015
%                 screenTreater(at(cin3 , cin3)) = - (cytoEff(2) * leep);
%                 
%                 screenTreater(at(normal , cin3)) = (cytoEff(2) * leep);
%                 
%                 % cervical cancer screen and treat, death
%                 for s = 5 : 7
%                     for p = 1 : periods
%                         ccR = squeeze(ccRInds(d , v , h , s , p , a , :)); % s refers to region of cervical cancer (local, regional, distant)
%                         screenTreater(at(ccR , ccR)) = ...
%                             - (detCC(s - 4) .* ccTreat + muCC(p , s - 4) * hivCC(s - 4)); % detection probability and death given region. Add in probability of successful treatment later!
%                         
%                         % CC -> susceptible
%                         cc2Sus = squeeze(ccSusInds(d , v , h , a , :));
%                         screenTreater(at(cc2Sus , ccR)) = ...
%                             detCC(s - 4) .* ccTreat;
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% screenTreater2015 = screenTreater;
% disp('Building screening and treatment matrix for 2015 onwards')
% for v = 1 : viral
%     for d = 1 : disease
%         fScreen = screenFreq(: , 2);
%         if d >= 2 && d <= 6 || d == 10 % Infected or low CD4 or HIV+ on ART
%             fScreen = screenFreq(: , 1);
%         end
%         hpvEff = screenCover * fScreen(1) * hpvSens; % effective screen and detect rate by hpv test
%         cytoEff = screenCover * fScreen(2) * cytoSens; % effective screen and detect rate by cytology
%         for h = 2 : hpvTypes % infected onwards
%             for a = 5 : age
%                 %                 cin1 = cin1Inds(d , v , h , a , :);
%                 cin2 = cin2Inds(d , v , h , a , :);
%                 cin3 = cin3Inds(d , v , h , a , :);
%                 normal = normalInds(d , v , 2 , a , :);
%                 % 35-55 years
%                 if a > 7 && a < 12
%                     cin3_35Plus = cin3Inds(d , v , h , a , :);
%                     screenTreater2015(at(cin3_35Plus , cin3_35Plus)) = ...
%                         - (hpvEff(2) * leep); % hpv test screen before leep
%                     
%                     screenTreater2015(at(normal , cin3_35Plus)) = hpvEff(2) * leep;
%                 end
%                 % 25-35 years old
%                 if a > 5 && a < 8
%                     if d == 1 % HIV negative and 25-35 years old
%                         screenTreat = cytoEff(2) * hpvEff(2) * leep; % detect by cytology and hpv test then leep
%                     else % HIV positive and 25-35 years old
%                         screenTreat = hpvEff(2) * leep; % detect by hpv test then leep
%                     end
%                     cin3_25_35 = cin3Inds(d , v , h , a , :);
%                     screenTreater2015(at(cin3_25_35 , cin3_25_35)) = ...
%                         - screenTreat; % test and leep
%                     
%                     screenTreater2015(normal , cin3_25_35) = screenTreat;
%                 end
%             end
%         end
%     end
% end
% save(fullfile(savedir ,'screenTreatMats') , 'screenTreater' , 'screenTreater2015')
% disp('Finished building HPV screening and treatment matrix')
%% bornDie
%% Fertility prior to 1995
load([paramDir ,'HIVParams'])
load([paramDir ,'popData'])
fertMat = spalloc(numel(pop) , numel(pop) , numel(pop));
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix accounting for uninfected mothers
disp('Building fertility matrix for uninfected mothers')
for a = 1 : age
    hivUninf = toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates ,...
        1 : periods , 2 , a , 1 : risk));
    hivPosArt = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates ,...
        1 : periods , 2 , a , 1 : risk));
    fertMat(at(negMaleBirth , hivUninf)) = 0.5 * fertility(a , 1);
    fertMat(at(negFemaleBirth , hivUninf)) = 0.5 * fertility(a , 1);
    fertMat(at(negMaleBirth , hivPosArt)) = 0.5 * fertility(a , 1);
    fertMat(at(negFemaleBirth , hivPosArt)) = 0.5 * fertility(a , 1);
end

% fertility matrix accounting for infected mothers
% time dependent terms: kHiv, update with term in MTCTRate corresponding to
% year
disp('Building fertility matrix for HIV-infected mothers')
hivFertPosBirth = spalloc(numel(pop) , numel(pop) , numel(pop));
hivFertNegBirth = hivFertPosBirth;
posMaleBirth = toInd(allcomb(2 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(2 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% kHiv = MTCTRate(1); % proportion of births from HIV-positive females that result in vertical transmission (before 2004)
%2005: 0.202 , >2008: 0.071. Will change depending on time step. Overwrite
%existing value in hivFertMat at appropriate time step.

for d = 2 : 6 % hiv infected
    for v = 1 : 5 % hiv infected
        for a = 1 : age
%             kHiv = MTCTRate(1);
            hivInfected = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , 1 : periods , 2 , a , 1 : risk));
            hivFertPosBirth(at(posMaleBirth , hivInfected)) = 0.5 * fertility(a , d); % * kHiv; 
            hivFertPosBirth(at(posFemaleBirth , hivInfected)) = 0.5 * fertility(a , d); % * kHiv;
            hivFertNegBirth(at(negMaleBirth , hivInfected)) = 0.5 * fertility(a , d); % * (1 - kHiv);
            hivFertNegBirth(at(negFemaleBirth , hivInfected)) = 0.5 * fertility(a , d);%  * (1 - kHiv);             
        end
    end
end

%fertMat = fertMat + hivFertMat; 
save(fullfile(savedir ,'fertMat') , 'fertMat')
save(fullfile(savedir ,'hivFertMats') , 'hivFertPosBirth' , 'hivFertNegBirth')
%% Fertility from 2005 onwards
load([paramDir ,'HIVParams'])
load([paramDir ,'popData'])
fertMat2 = spalloc(numel(pop) , numel(pop) , numel(pop));
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negFemaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% fertility matrix accounting for uninfected mothers
disp('Building fertility matrix for uninfected mothers for 2005 onwards')
for a = 1 : age
    hivUninf = toInd(allcomb(1 , 1 , 1 : hpvTypes , 1 : hpvStates ,...
        1 : periods , 2 , a , 1 : risk));
    hivPosArt = toInd(allcomb(10 , 6 , 1 : hpvTypes , 1 : hpvStates ,...
        1 : periods , 2 , a , 1 : risk));
    fertMat2(at(negMaleBirth , hivUninf)) = 0.5 * fertility2(a , 1);
    fertMat2(at(negFemaleBirth , hivUninf)) = 0.5 * fertility2(a , 1);
    fertMat2(at(negMaleBirth , hivPosArt)) = 0.5 * fertility2(a , 1);
    fertMat2(at(negFemaleBirth , hivPosArt)) = 0.5 * fertility2(a , 1);
end

% fertility matrix accounting for infected mothers
% time dependent terms: kHiv, update with term in MTCTRate corresponding to
% year
disp('Building fertility matrix for HIV-infected mothers for 2005 onwards')
hivFertPosBirth2 = spalloc(numel(pop) , numel(pop) , numel(pop));
hivFertNegBirth2 = hivFertPosBirth2;
posMaleBirth = toInd(allcomb(2 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
posFemaleBirth = toInd(allcomb(2 , 1 , 1 , 1 , 1 , 2 , 1 , 1));

% kHiv = MTCTRate(1); % proportion of births from HIV-positive females that result in vertical transmission (before 2004)
%2005: 0.202 , >2008: 0.071. Will change depending on time step. Overwrite
%existing value in hivFertMat at appropriate time step.

for d = 2 : 6 % hiv infected
    for v = 1 : 5 % hiv infected
        for a = 1 : age
%             kHiv = MTCTRate(1);
            hivInfected = toInd(allcomb(d , v , 1 : hpvTypes , 1 : hpvStates , 1 : periods , 2 , a , 1 : risk));
            hivFertPosBirth2(at(posMaleBirth , hivInfected)) = 0.5 * fertility2(a , d); % * kHiv; 
            hivFertPosBirth2(at(posFemaleBirth , hivInfected)) = 0.5 * fertility2(a , d); % * kHiv;
            hivFertNegBirth2(at(negMaleBirth , hivInfected)) = 0.5 * fertility2(a , d); % * (1 - kHiv);
            hivFertNegBirth2(at(negFemaleBirth , hivInfected)) = 0.5 * fertility2(a , d);%  * (1 - kHiv);             
        end
    end
end

%fertMat = fertMat + hivFertMat; 
save(fullfile(savedir ,'fertMat2') , 'fertMat2')
save(fullfile(savedir ,'hivFertMats2') , 'hivFertPosBirth2' , 'hivFertNegBirth2')

%% 
% Background deaths
disp('Building death matrix')
deathMat = spalloc(numel(pop) , numel(pop) , numel(pop));
for a = 1 : age
    males = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , 1 , a , 1 : risk));
    females = toInd(allcomb(1 : disease , 1 : viral , 1 : hpvTypes , 1 : hpvStates , 1 : periods , 2 , a , 1 : risk));
    deathMat(at(males , males)) = - mue(a , 1);
    deathMat(at(females , females)) = - mue(a , 2);
end
save(fullfile(savedir ,'deathMat') , 'deathMat')
disp('Death matrix complete')
%%
% Vaccination (model time dependent). Activate when vax year begins.
% ALPHA version: assumes vaccination rate is dependent on age. Female
% vaccination only.
disp('Building vaccination matrix')
vaxer = speye(numel(pop) , numel(pop));
V = zeros(gender , age);
% V(2 , 3) = 1; % turn on vaccination for females in 10 - 14 age group
% for a = 1 : age
%     susFemale = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 : periods , 2 , a , 1 : risk));
%     vaxdFemale = toInd(allcomb(1 : disease , 1 : viral , 1 , 10 , 1 : periods , 2 , a , 1 : risk));
%     vaxer(at(vaxdFemale , susFemale)) = V(2 , a);
%     vaxer(at(susFemale , susFemale)) = -V(2 , a);
%     
%     % for males (future version?)
% %     susMale = toInd(allcomb(1 : disease , 1 : viral , 1 , 1 , 1 : periods , 1 , a , 1 : risk));
% %     vaxdMale = toInd(allcomb(1 : disease , 1 : viral , 5 , 6 , 1 : periods , 1 , a , 1 : risk));
% %     vaxer(vaxdMale , susMale) = V(2 , a);
% %     vaxer(susMale , susMale) = -V(2 , a);
% end
save(fullfile(savedir ,'vaxer') , 'vaxer')
disp('Vaccination matrix complete')
%% Make circumcision matrix before current year (use when circumcision begins in model)
disp('Building circumcision matrix')
negCircMaleBirth = toInd(allcomb(7 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
circMat = spalloc(numel(pop) , numel(pop) , 2);
circMat(at(negCircMaleBirth , negMaleBirth)) = circ(1);
circMat(at(negMaleBirth , negMaleBirth)) = - circ(1);
save(fullfile(savedir ,'circMat') , 'circMat')
disp('Circumcision matrix complete')

%% Make circumcision matrix after 2030 
disp('Building circumcision matrix after 2030')
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negCircMaleBirth = toInd(allcomb(7 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
circMat2 = spalloc(numel(pop) , numel(pop) , 2);
circMat2(at(negCircMaleBirth , negMaleBirth)) = 0.70;
circMat2(at(negMaleBirth , negMaleBirth)) = - 0.70;
save(fullfile(savedir ,'circMat2') , 'circMat2')
disp('Circumcision matrix after 2030 complete')

%% Make circumcision matrix for HIV-only scenarios (use when circumcision begins in model)
disp('Building circumcision matrix 2')
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negCircMaleBirth = toInd(allcomb(7 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
circMale16_29 = toInd(allcomb(7 , 1 : viral , 1 , 1 , 1 , 1 , 4 : 6 , 1 : risk)); % circumcised
circMat2 = spalloc(numel(pop) , numel(pop) , 2);
for d = 1 : disease
    male16_29 = toInd(allcomb(d , 1 : viral , 1 , 1 , 1 , 1 , 4 : 6 , 1 : risk)); % all
    circMat2(at(circMale16_29 , male16_29)) = 0.9;
    circMat2(at(male16_29 , male16_29)) = -0.9;
end
save(fullfile(savedir ,'circMat2') , 'circMat2')
disp('Circumcision matrix 2 complete')

disp('Building circumcision matrix 2B')
negMaleBirth = toInd(allcomb(1 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
negCircMaleBirth = toInd(allcomb(7 , 1 , 1 , 1 , 1 , 1 , 1 , 1));
circMale16_29 = toInd(allcomb(7 , 1 : viral , 1 , 1 , 1 , 1 , 4 : 6 , 1 : risk)); % circumcised
circMat2B = spalloc(numel(pop) , numel(pop) , 2);

for d = 1 : disease
    male16_29 = toInd(allcomb(d , 1 : viral , 1 , 1 , 1 , 1 , 4 : 6 , 1 : risk)); % all
    circMat2B(at(circMale16_29 , male16_29)) = 0.4;
    circMat2B(at(male16_29 , male16_29)) = -0.4;   
end
save(fullfile(savedir ,'circMat2B') , 'circMat2B')
disp('Circumcision matrix 2B complete')
disp(' ')
disp('Matrix construction complete')
disp('All matrices saved to current directory')


