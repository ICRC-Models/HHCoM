%function vaxCEA_multSims_CISNETsces()

%% LOAD PARAMETERS
paramDir = [pwd , '\Params\'];

[stepsPerYear , timeStep , startYear , currYear , endYear , ...
    years , disease , viral , hpvVaxStates , hpvNonVaxStates , endpoints , ...
    intervens , gender , age , risk , hpvTypeGroups , dim , k , toInd , ...
    annlz , ...
    ageSexDebut , mInit , fInit , partnersM , partnersF , maleActs , ...
    femaleActs , riskDist , fertility , fertility2 , fertility3 , fertility4 , ...
    mue , mue2 , mue3 , mue4 , epsA_vec , epsR_vec , ...
    yr , ...
    hivOn , betaHIV_mod , muHIV , kCD4 , ...
    hpvOn , beta_hpvVax_mod , beta_hpvNonVax_mod , fImm , rImmune , ...
    kCin1_Inf , kCin2_Cin1 , kCin3_Cin2 , kCC_Cin3 , rNormal_Inf , kInf_Cin1 , ...
    kCin1_Cin2 , kCin2_Cin3 , lambdaMultImm , hpv_hivClear , rImmuneHiv , ...
    c3c2Mults , c2c1Mults , c2c3Mults , c1c2Mults , muCC , kRL , kDR , artHpvMult , ...
    hpv_hivMult , maleHpvClearMult , ...
    condUse , screenYrs , hpvScreenStartYear , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear , circNatStartYear , vaxStartYear , ...
    baseline , who , spCyto , spHpvDna , spGentyp , spAve , spHpvAve , ...
    circProtect , condProtect , MTCTRate , hyst , ...
    OMEGA , ...
    ccInc2012_dObs , ccInc2018_dObs , cc_dist_dObs , cin3_dist_dObs , ...
    cin1_dist_dObs , hpv_dist_dObs , cinPos2002_dObs , cinNeg2002_dObs , ...
    cinPos2015_dObs , cinNeg2015_dObs , hpv_hiv_dObs , hpv_hivNeg_dObs , ...
    hpv_hivM2008_dObs , hpv_hivMNeg2008_dObs , hivPrevM_dObs , hivPrevF_dObs , ...
    popAgeDist_dObs , totPopSize_dObs , ...
    hivCurr , ...
    gar , hivSus , hpvVaxSus , hpvVaxImm , hpvNonVaxSus , hpvNonVaxImm , ...
    toHiv , vaxInds , nonVInds , hpvVaxInf , hpvNonVaxInf , ...
    hivInds , ...
    cin3hpvVaxIndsFrom , ccLochpvVaxIndsTo , ccLochpvVaxIndsFrom , ...
    ccReghpvVaxInds , ccDisthpvVaxInds , cin3hpvNonVaxIndsFrom , ...
    ccLochpvNonVaxIndsTo , ccLochpvNonVaxIndsFrom , ccReghpvNonVaxInds , ...
    ccDisthpvNonVaxInds , cin1hpvVaxInds , cin2hpvVaxInds , cin3hpvVaxInds , ...
    cin1hpvNonVaxInds , cin2hpvNonVaxInds , cin3hpvNonVaxInds , normalhpvVaxInds , ...
    immunehpvVaxInds , infhpvVaxInds , normalhpvNonVaxInds , immunehpvNonVaxInds , ...
    infhpvNonVaxInds , fromVaxNoScrnInds , fromVaxScrnInds , toNonVaxNoScrnInds , ...
    toNonVaxScrnInds , ageInd , riskInd , ...
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3 , dFertMat3 , deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2(1 , 0 , [] , [] , []);

%% LOAD SAVED RESULTS

% ***SET ME***: name of directory where files are located
dirName = '22Apr20Ph2V11_2v57BaseVax_spCytoScreen_noVMMChpv_currYr2021_CISNET-S0_6_1';

% ***SET ME***: names of teams with data to plot; should match team names in input file names
teamVec = {'UW' , 'DC'};
nTeams = length(teamVec);
% ***SET ME***: date of most recent data to plot; should match dates in input file names from each team, respectively
dateVec = {'07192021' , 'v5'};
% ***SET ME***: scenario numbers to plot; should match scenario names in input file names
% sceVec = {'S0' , 'S1' , 'S2'};
sceVec = {'S0', 'S0b', 'S1', 'S2', 'S3'}; 
nSces = length(sceVec);

% ***SET ME***: legend titles
% plotTits = {'' , 'UW: No HIV' , '      : HIV, no ART' , '      : HIV, observed ART' , ...
%     'DC: No HIV' , '      : HIV, no ART' , '      : HIV, observed ART'};
plotTits = {'' , 'UW: No HIV, no HIV-associated interventions' , '      : No HIV', ...
    '      : HIV, no ART' , '      : HIV, observed ART' , '      : HIV, ART & CC intervention scale-up', ...
    'DC: No HIV, no HIV-associated interventions' , '      : No HIV', ...
    '      : HIV, no ART' , '      : HIV, observed ART', '      : HIV, ART & CC intervention scale-up'};
plotTits2 = {'UW: Total' , '      : HIV-negative' , '      : HIV-positive, all' , ...
    'DC: Total' , '      : HIV-negative' , '      : HIV-positive, all'};
plotTits3 = {'Total' , 'HIV-negative' , 'HIV-positive, all' , 'HIV-positive, untreated' , 'HIV-positive, on ART'};
colorVec = {[0 0.4470 0.7410] , [0.4660 , 0.6740 , 0.1880] , [0.8500 0.3250 0.0980]}; % plot line colors
% styleVec = {'-' , '--' , ':'}; % plot line styles (for plots by scenario)
styleVec = {'-' , '--' , ':', '-.' , '-'}; % plot line styles (for plots by scenario)
widthVec = {4.0 , 0.5 , 2.0}; % pot line widths (for plots by HIV status)

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)


%% Crude HPV prevalence over time
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
t82on = (1982:2121)';
t82onLen = length(t82on);

xline(2021);
for i = 1 : nTeams
    for j = 1 : nSces  
        for dInd = 1 : 1 %5
            % Load results
            fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
            hpvPrevTime = readmatrix(fname , 'Sheet' , 'HPV prevalence' , 'Range' , ['C' , num2str(((dInd-1)*t82onLen+1+3)) , ':C' , num2str(dInd*t82onLen+3)]);
            hold all;
            p = plot(t82on , hpvPrevTime(:,1) , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
        end
    end
end

legend(plotTits{:} , 'Location' , 'northeast');
xlim([1980 2121]); ylim([0 1.0]);
grid on;
xlabel('Year'); ylabel('Crude HPV prevalence among women');

%% Age-standardized HPV prevalence over time
% Note: the age-standardization process shifts the prevalence of the
% last modelled age group to the next age group in the following year.
% However, prevalence is NaN prior to HIV introduction in the
% HIV-positive no ART group, and NaN prior to ART introduction in the
% HIV-positive ART group. Since we have four age groups past the 16 we
% model, a NaN value is present for four years past the introduction of
% HIV/ART, leading to a NaN value for summed HPV infected during these 
% years. We therefore lack data in this four-year interval in the
% saved/plotted results.
worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
    247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
    23477 9261 2155];

fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
t82on = (1982:2121)';
t82onLen = length(t82on);

xline(2021);
for i = 1 : nTeams
    for j = 1 : nSces  
        for dInd = 1 : 1 %5
            % Load results
            fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
            hpv_hivAgeW_dis = readmatrix(fname , 'Sheet' , 'HPV prevalence' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
            
            numHpvTot = zeros(1 , size(hpv_hivAgeW_dis,2));       
            for aInd = 1:age+4
                a = aInd;
                if aInd >= age
                    a = age;
                end
                if aInd <= age    
                    numHpv = hpv_hivAgeW_dis(a , :) .* worldStandard_WP2015(aInd);
                    if (a < 3)
                        numHpv = zeros(1 , size(hpv_hivAgeW_dis,2));
                    end
                elseif aInd > age
                    numHpv = hpv_hivAgeW_dis(a , :);
                    numHpv = cat(2 , (ones(1,aInd-a).*numHpv(1,1)) , numHpv(1 ,1:end-(aInd-a)));
                    numHpv = numHpv .* worldStandard_WP2015(aInd);
                end
                numHpvTot = numHpvTot + numHpv;
            end
            hpvPrevTot = numHpvTot ./ (sum(worldStandard_WP2015(1:age+4)));

            hold all;
            p = plot(t82on , hpvPrevTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
        end
    end
end

legend(plotTits{:} , 'Location' , 'northeast');
xlim([1980 2121]); ylim([0 0.25]);
grid on;
xlabel('Year'); ylabel('AS-HPV prevalence among women');

                            %% Facet wrap the above: No HIV, no HIV associated interventions
                            fig = figure;
                            set(fig,'DefaultAxesFontSize' , 18);
                            t82on = (1982:2121)';
                            t82onLen = length(t82on);
                            
                            xline(2021);
                            for i = 1 : nTeams
                            %     for j = 1 : nSces  
                                for j = 1 % trying to facet wrap
                                    for dInd = 1 : 1 %5
                                        % Load results
                                        fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                                            teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
                                        hpv_hivAgeW_dis = readmatrix(fname , 'Sheet' , 'HPV prevalence' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
                                        
                                        numHpvTot = zeros(1 , size(hpv_hivAgeW_dis,2));       
                                        for aInd = 1:age+4
                                            a = aInd;
                                            if aInd >= age
                                                a = age;
                                            end
                                            if aInd <= age    
                                                numHpv = hpv_hivAgeW_dis(a , :) .* worldStandard_WP2015(aInd);
                                                if (a < 3)
                                                    numHpv = zeros(1 , size(hpv_hivAgeW_dis,2));
                                                end
                                            elseif aInd > age
                                                numHpv = hpv_hivAgeW_dis(a , :);
                                                numHpv = cat(2 , (ones(1,aInd-a).*numHpv(1,1)) , numHpv(1 ,1:end-(aInd-a)));
                                                numHpv = numHpv .* worldStandard_WP2015(aInd);
                                            end
                                            numHpvTot = numHpvTot + numHpv;
                                        end
                                        hpvPrevTot = numHpvTot ./ (sum(worldStandard_WP2015(1:age+4)));
                            
                                        hold all;
                                        p = plot(t82on , hpvPrevTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
                                    end
                                end
                            end
                            
                            % legend(plotTits{:} , 'Location' , 'northeast');
                            legend({'', 'UW: No HIV, no HIV-associated interventions', 'DC: No HIV, no HIV-associated interventions'}, 'Location', 'northeast');
                            % xlim([1980 2121]); ylim([0 0.5]);
                            xlim([1980 2121]); ylim([0 0.25]);
                            grid on;
                            xlabel('Year'); ylabel('AS-HPV prevalence among women');
                
                            %% Facet wrap the above: no HIV
                            fig = figure;
                            set(fig,'DefaultAxesFontSize' , 18);
                            t82on = (1982:2121)';
                            t82onLen = length(t82on);
                            
                            xline(2021);
                            for i = 1 : nTeams
                            %     for j = 1 : nSces  
                                for j = 2 % trying to facet wrap
                                    for dInd = 1 : 1 %5
                                        % Load results
                                        fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                                            teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
                                        hpv_hivAgeW_dis = readmatrix(fname , 'Sheet' , 'HPV prevalence' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
                                        
                                        numHpvTot = zeros(1 , size(hpv_hivAgeW_dis,2));       
                                        for aInd = 1:age+4
                                            a = aInd;
                                            if aInd >= age
                                                a = age;
                                            end
                                            if aInd <= age    
                                                numHpv = hpv_hivAgeW_dis(a , :) .* worldStandard_WP2015(aInd);
                                                if (a < 3)
                                                    numHpv = zeros(1 , size(hpv_hivAgeW_dis,2));
                                                end
                                            elseif aInd > age
                                                numHpv = hpv_hivAgeW_dis(a , :);
                                                numHpv = cat(2 , (ones(1,aInd-a).*numHpv(1,1)) , numHpv(1 ,1:end-(aInd-a)));
                                                numHpv = numHpv .* worldStandard_WP2015(aInd);
                                            end
                                            numHpvTot = numHpvTot + numHpv;
                                        end
                                        hpvPrevTot = numHpvTot ./ (sum(worldStandard_WP2015(1:age+4)));
                            
                                        hold all;
                                        p = plot(t82on , hpvPrevTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
                                    end
                                end
                            end
                            
                            legend({'', 'UW: No HIV', 'DC: No HIV'}, 'Location', 'northeast');
                            xlim([1980 2121]); ylim([0 0.25]);
                            grid on;
                            xlabel('Year'); ylabel('AS-HPV prevalence among women');
                
                            %% Facet wrap the above: HIV, no ART
                            fig = figure;
                            set(fig,'DefaultAxesFontSize' , 18);
                            t82on = (1982:2121)';
                            t82onLen = length(t82on);
                            
                            xline(2021);
                            for i = 1 : nTeams
                            %     for j = 1 : nSces  
                                for j = 3 % trying to facet wrap
                                    for dInd = 1 : 1 %5
                                        % Load results
                                        fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                                            teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
                                        hpv_hivAgeW_dis = readmatrix(fname , 'Sheet' , 'HPV prevalence' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
                                        
                                        numHpvTot = zeros(1 , size(hpv_hivAgeW_dis,2));       
                                        for aInd = 1:age+4
                                            a = aInd;
                                            if aInd >= age
                                                a = age;
                                            end
                                            if aInd <= age    
                                                numHpv = hpv_hivAgeW_dis(a , :) .* worldStandard_WP2015(aInd);
                                                if (a < 3)
                                                    numHpv = zeros(1 , size(hpv_hivAgeW_dis,2));
                                                end
                                            elseif aInd > age
                                                numHpv = hpv_hivAgeW_dis(a , :);
                                                numHpv = cat(2 , (ones(1,aInd-a).*numHpv(1,1)) , numHpv(1 ,1:end-(aInd-a)));
                                                numHpv = numHpv .* worldStandard_WP2015(aInd);
                                            end
                                            numHpvTot = numHpvTot + numHpv;
                                        end
                                        hpvPrevTot = numHpvTot ./ (sum(worldStandard_WP2015(1:age+4)));
                            
                                        hold all;
                                        p = plot(t82on , hpvPrevTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
                                    end
                                end
                            end
                            
                            legend({'', 'UW: HIV, no ART', 'DC: HIV, no ART'}, 'Location', 'northeast');
                            xlim([1980 2121]); ylim([0 0.25]);
                            grid on;
                            xlabel('Year'); ylabel('AS-HPV prevalence among women');
                
                            %% Facet wrap the above: HIV, observed ART
                            fig = figure;
                            set(fig,'DefaultAxesFontSize' , 18);
                            t82on = (1982:2121)';
                            t82onLen = length(t82on);
                            
                            xline(2021);
                            for i = 1 : nTeams
                            %     for j = 1 : nSces  
                                for j = 4 % trying to facet wrap
                                    for dInd = 1 : 1 %5
                                        % Load results
                                        fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                                            teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
                                        hpv_hivAgeW_dis = readmatrix(fname , 'Sheet' , 'HPV prevalence' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
                                        
                                        numHpvTot = zeros(1 , size(hpv_hivAgeW_dis,2));       
                                        for aInd = 1:age+4
                                            a = aInd;
                                            if aInd >= age
                                                a = age;
                                            end
                                            if aInd <= age    
                                                numHpv = hpv_hivAgeW_dis(a , :) .* worldStandard_WP2015(aInd);
                                                if (a < 3)
                                                    numHpv = zeros(1 , size(hpv_hivAgeW_dis,2));
                                                end
                                            elseif aInd > age
                                                numHpv = hpv_hivAgeW_dis(a , :);
                                                numHpv = cat(2 , (ones(1,aInd-a).*numHpv(1,1)) , numHpv(1 ,1:end-(aInd-a)));
                                                numHpv = numHpv .* worldStandard_WP2015(aInd);
                                            end
                                            numHpvTot = numHpvTot + numHpv;
                                        end
                                        hpvPrevTot = numHpvTot ./ (sum(worldStandard_WP2015(1:age+4)));
                            
                                        hold all;
                                        p = plot(t82on , hpvPrevTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
                                    end
                                end
                            end
                            
                            legend({'', 'UW: HIV, observed ART', 'DC: HIV, observed ART'}, 'Location', 'northeast');
                            xlim([1980 2121]); ylim([0 0.25]);
                            grid on;
                            xlabel('Year'); ylabel('AS-HPV prevalence among women');

                            %% Facet wrap the above: HIV, observed ART
                            fig = figure;
                            set(fig,'DefaultAxesFontSize' , 18);
                            t82on = (1982:2121)';
                            t82onLen = length(t82on);
                            
                            xline(2021);
                            for i = 1 : nTeams
                            %     for j = 1 : nSces  
                                for j = 5 % trying to facet wrap
                                    for dInd = 1 : 1 %5
                                        % Load results
                                        fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                                            teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
                                        hpv_hivAgeW_dis = readmatrix(fname , 'Sheet' , 'HPV prevalence' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
                                        
                                        numHpvTot = zeros(1 , size(hpv_hivAgeW_dis,2));       
                                        for aInd = 1:age+4
                                            a = aInd;
                                            if aInd >= age
                                                a = age;
                                            end
                                            if aInd <= age    
                                                numHpv = hpv_hivAgeW_dis(a , :) .* worldStandard_WP2015(aInd);
                                                if (a < 3)
                                                    numHpv = zeros(1 , size(hpv_hivAgeW_dis,2));
                                                end
                                            elseif aInd > age
                                                numHpv = hpv_hivAgeW_dis(a , :);
                                                numHpv = cat(2 , (ones(1,aInd-a).*numHpv(1,1)) , numHpv(1 ,1:end-(aInd-a)));
                                                numHpv = numHpv .* worldStandard_WP2015(aInd);
                                            end
                                            numHpvTot = numHpvTot + numHpv;
                                        end
                                        hpvPrevTot = numHpvTot ./ (sum(worldStandard_WP2015(1:age+4)));
                            
                                        hold all;
                                        p = plot(t82on , hpvPrevTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
                                    end
                                end
                            end
                            
                            legend({'', 'UW: HIV, ART & CC intervention scale-up', 'DC: HIV, ART & CC intervention scale-up'}, 'Location', 'northeast');
                            xlim([1980 2121]); ylim([0 0.25]);
                            grid on;
                            xlabel('Year'); ylabel('AS-HPV prevalence among women');

%% Crude HPV prevalence by age
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79' , '80+'};
year2021 = (2002-1982) + 1; 
t82on = (1982:2121)';
t82onLen = length(t82on);

% Calibration data error bars
meanObs = hpv_hiv_dObs(: , 2);
sdevObs = (hpv_hiv_dObs(: , 3).^(1/2)).*2;
meanNeg = hpv_hivNeg_dObs(: , 2);
sdevNeg = (hpv_hivNeg_dObs(: , 3).^(1/2)).*2;
errorbar(4 : length(meanObs)+4-1 , meanNeg , sdevNeg , ...
    'rs' , 'LineWidth' , 0.5);
hold all;
errorbar(4 : length(meanObs)+4-1 , meanObs , sdevObs , ...
    'rs' , 'LineWidth' , 2.0);
for i = 1 : nTeams
%     for j = 3 : nSces 
    for j = 4 : 4; 
        for dInd = 1 : 3 %5
            % Load results
            fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
            hpvPrevTime = readmatrix(fname , 'Sheet' , 'HPV prevalence' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':T' , num2str(dInd*t82onLen+3)]);
            hold all;
            p = plot([1:length(ageGroup)] , hpvPrevTime(year2021,:) , 'Color' , colorVec{i} , 'LineWidth' , widthVec{dInd});
            set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
        end
    end
end
legend('(McDonald, 2014) HIV-negative, Observed Cape Town: mean, 2SD' , '(McDonald, 2014) HIV-positive, Observed Cape Town: mean, 2SD' , plotTits2{:} , 'Location' , 'northeast');
xlim([0 length(ageGroup)]); ylim([0 1.0]);
grid on;
xlabel('Age Group'); ylabel('Crude HPV prevalence among women'); title('S2: HIV, observed ART; Year: 2002')

%% Crude cervical cancer incidence over time
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
t82on = (1982:2121)';
t82onLen = length(t82on);

xline(2021);
for i = 1 : nTeams
    for j = 1 : nSces  
        for dInd = 1 : 1 %5
            % Load results
            fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
            ccIncTime = readmatrix(fname , 'Sheet' , 'CC incidence rate' , 'Range' , ['C' , num2str(((dInd-1)*t82onLen+1+3)) , ':C' , num2str(dInd*t82onLen+3)]);
            hold all;
            p = plot(t82on , ccIncTime(:,1) , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
        end
    end
end

legend(plotTits{:} , 'Location' , 'northeast');
xlim([1980 2121]); %ylim([0 1.0]);
grid on;
xlabel('Year'); ylabel('Crude cervical cancer incidence');

%% Age-standardized cervical cancer incidence over time
% Note: the age-standardization process shifts the incidence rate of the
% last modelled age group to the next age group in the following year.
% However, CC incidence is NaN prior to HIV introduction in the
% HIV-positive no ART group, and NaN prior to ART introduction in the
% HIV-positive ART group. Since we have four age groups past the 16 we
% model, a NaN value is present for four years past the introduction of
% HIV/ART, leading to a NaN value for summed incidence during these 
% years. We therefore lack data in this four-year interval in the
% saved/plotted results.
fac = 10 ^ 5;
worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
    247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
    23477 9261 2155];

fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
t82on = (1982:2121)';
t82onLen = length(t82on);

xline(2021);
for i = 1 : nTeams
    for j = 1 : nSces  
        for dInd = 1 : 1 %5
            % Load results
            fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
            ccIncAge_dis = readmatrix(fname , 'Sheet' , 'CC incidence rate' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
            
            numCCTot = zeros(1 , size(ccIncAge_dis,2));       
            for aInd = 1:age+4
                a = aInd;
                if aInd >= age
                    a = age;
                end
                if aInd <= age    
                    numCC = ccIncAge_dis(a , :) .* worldStandard_WP2015(aInd);
                    if (a < 3)
                        numCC = zeros(1 , size(ccIncAge_dis,2));
                    end
                elseif aInd > age
                    numCC = ccIncAge_dis(a , :);
                    numCC = cat(2 , (ones(1,aInd-a).*numCC(1,1)) , numCC(1 ,1:end-(aInd-a)));
                    numCC = numCC .* worldStandard_WP2015(aInd);
                end
                numCCTot = numCCTot + numCC;
            end
            ccIncTot = numCCTot ./ (sum(worldStandard_WP2015(1:age+4)));

            hold all;
            p = plot(t82on , ccIncTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
        end
    end
end
hold all;
plot(t82on , ones(1,size(t82on,1)).*4.0 , 'r:')
hold all;
plot(t82on , ones(1,size(t82on,1)).*10.0 , 'r--')

legend(plotTits{:} , 'Elimination: <4/100K' , 'Benchmark: <10/100K' , 'Location' , 'northeast');
xlim([1980 2121]); ylim([0 90]);
grid on;
xlabel('Year'); ylabel('AS-cervical cancer incidence (per 100K)');

                        %% Stratify by scenario: No HIV, no HIV-associated interventions
                        fig = figure;
                        set(fig,'DefaultAxesFontSize' , 18);
                        t82on = (1982:2121)';
                        t82onLen = length(t82on);
                        
                        xline(2021);
                        for i = 1 : nTeams
                            for j = 1 : 1 % stratify by scenario  
                                for dInd = 1 : 1 %5
                                    % Load results
                                    fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                                        teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
                                    ccIncAge_dis = readmatrix(fname , 'Sheet' , 'CC incidence rate' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
                                    
                                    numCCTot = zeros(1 , size(ccIncAge_dis,2));       
                                    for aInd = 1:age+4
                                        a = aInd;
                                        if aInd >= age
                                            a = age;
                                        end
                                        if aInd <= age    
                                            numCC = ccIncAge_dis(a , :) .* worldStandard_WP2015(aInd);
                                            if (a < 3)
                                                numCC = zeros(1 , size(ccIncAge_dis,2));
                                            end
                                        elseif aInd > age
                                            numCC = ccIncAge_dis(a , :);
                                            numCC = cat(2 , (ones(1,aInd-a).*numCC(1,1)) , numCC(1 ,1:end-(aInd-a)));
                                            numCC = numCC .* worldStandard_WP2015(aInd);
                                        end
                                        numCCTot = numCCTot + numCC;
                                    end
                                    ccIncTot = numCCTot ./ (sum(worldStandard_WP2015(1:age+4)));
                        
                                    hold all;
                                    p = plot(t82on , ccIncTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
                                end
                            end
                        end
                        hold all;
                        plot(t82on , ones(1,size(t82on,1)).*4.0 , 'r:')
                        hold all;
                        plot(t82on , ones(1,size(t82on,1)).*10.0 , 'r--')
                        
                        legend({'', 'UW: No HIV, no HIV-associated interventions', 'DC: No HIV, no HIV-associated interventions' , 'Elimination: <4/100K' , 'Benchmark: <10/100K'} , 'Location' , 'northeast');
                        xlim([1980 2121]); ylim([0 90]);
                        grid on;
                        xlabel('Year'); ylabel('AS-cervical cancer incidence (per 100K)');

                        %% Stratify by scenario: No HIV
                        fig = figure;
                        set(fig,'DefaultAxesFontSize' , 18);
                        t82on = (1982:2121)';
                        t82onLen = length(t82on);
                        
                        xline(2021);
                        for i = 1 : nTeams
                            for j = 2:2 % stratify by scenario  
                                for dInd = 1 : 1 %5
                                    % Load results
                                    fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                                        teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
                                    ccIncAge_dis = readmatrix(fname , 'Sheet' , 'CC incidence rate' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
                                    
                                    numCCTot = zeros(1 , size(ccIncAge_dis,2));       
                                    for aInd = 1:age+4
                                        a = aInd;
                                        if aInd >= age
                                            a = age;
                                        end
                                        if aInd <= age    
                                            numCC = ccIncAge_dis(a , :) .* worldStandard_WP2015(aInd);
                                            if (a < 3)
                                                numCC = zeros(1 , size(ccIncAge_dis,2));
                                            end
                                        elseif aInd > age
                                            numCC = ccIncAge_dis(a , :);
                                            numCC = cat(2 , (ones(1,aInd-a).*numCC(1,1)) , numCC(1 ,1:end-(aInd-a)));
                                            numCC = numCC .* worldStandard_WP2015(aInd);
                                        end
                                        numCCTot = numCCTot + numCC;
                                    end
                                    ccIncTot = numCCTot ./ (sum(worldStandard_WP2015(1:age+4)));
                        
                                    hold all;
                                    p = plot(t82on , ccIncTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
                                end
                            end
                        end
                        hold all;
                        plot(t82on , ones(1,size(t82on,1)).*4.0 , 'r:')
                        hold all;
                        plot(t82on , ones(1,size(t82on,1)).*10.0 , 'r--')
                        
                        legend({'', 'UW: No HIV', 'DC: No HIV' , 'Elimination: <4/100K' , 'Benchmark: <10/100K'} , 'Location' , 'northeast');
                        xlim([1980 2121]); ylim([0 90]);
                        grid on;
                        xlabel('Year'); ylabel('AS-cervical cancer incidence (per 100K)');

                        %% Stratify by scenario: HIV, no ART
                        fig = figure;
                        set(fig,'DefaultAxesFontSize' , 18);
                        t82on = (1982:2121)';
                        t82onLen = length(t82on);
                        
                        xline(2021);
                        for i = 1 : nTeams
                            for j = 3:3 % stratify by scenario  
                                for dInd = 1 : 1 %5
                                    % Load results
                                    fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                                        teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
                                    ccIncAge_dis = readmatrix(fname , 'Sheet' , 'CC incidence rate' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
                                    
                                    numCCTot = zeros(1 , size(ccIncAge_dis,2));       
                                    for aInd = 1:age+4
                                        a = aInd;
                                        if aInd >= age
                                            a = age;
                                        end
                                        if aInd <= age    
                                            numCC = ccIncAge_dis(a , :) .* worldStandard_WP2015(aInd);
                                            if (a < 3)
                                                numCC = zeros(1 , size(ccIncAge_dis,2));
                                            end
                                        elseif aInd > age
                                            numCC = ccIncAge_dis(a , :);
                                            numCC = cat(2 , (ones(1,aInd-a).*numCC(1,1)) , numCC(1 ,1:end-(aInd-a)));
                                            numCC = numCC .* worldStandard_WP2015(aInd);
                                        end
                                        numCCTot = numCCTot + numCC;
                                    end
                                    ccIncTot = numCCTot ./ (sum(worldStandard_WP2015(1:age+4)));
                        
                                    hold all;
                                    p = plot(t82on , ccIncTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
                                end
                            end
                        end
                        hold all;
                        plot(t82on , ones(1,size(t82on,1)).*4.0 , 'r:')
                        hold all;
                        plot(t82on , ones(1,size(t82on,1)).*10.0 , 'r--')
                        
                        legend({'', 'UW: HIV, no ART', 'DC: HIV, no ART' , 'Elimination: <4/100K' , 'Benchmark: <10/100K'} , 'Location' , 'northeast');
                        xlim([1980 2121]); ylim([0 90]);
                        grid on;
                        xlabel('Year'); ylabel('AS-cervical cancer incidence (per 100K)');

                        %% Stratify by scenario: HIV, observed ART
                        fig = figure;
                        set(fig,'DefaultAxesFontSize' , 18);
                        t82on = (1982:2121)';
                        t82onLen = length(t82on);
                        
                        xline(2021);
                        for i = 1 : nTeams
                            for j = 4 
                                for dInd = 1 : 1 %5
                                    % Load results
                                    fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                                        teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
                                    ccIncAge_dis = readmatrix(fname , 'Sheet' , 'CC incidence rate' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)])';
                                    
                                    numCCTot = zeros(1 , size(ccIncAge_dis,2));       
                                    for aInd = 1:age+4
                                        a = aInd;
                                        if aInd >= age
                                            a = age;
                                        end
                                        if aInd <= age    
                                            numCC = ccIncAge_dis(a , :) .* worldStandard_WP2015(aInd);
                                            if (a < 3)
                                                numCC = zeros(1 , size(ccIncAge_dis,2));
                                            end
                                        elseif aInd > age
                                            numCC = ccIncAge_dis(a , :);
                                            numCC = cat(2 , (ones(1,aInd-a).*numCC(1,1)) , numCC(1 ,1:end-(aInd-a)));
                                            numCC = numCC .* worldStandard_WP2015(aInd);
                                        end
                                        numCCTot = numCCTot + numCC;
                                    end
                                    ccIncTot = numCCTot ./ (sum(worldStandard_WP2015(1:age+4)));
                        
                                    hold all;
                                    p = plot(t82on , ccIncTot , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
                                end
                            end
                        end
                        hold all;
                        plot(t82on , ones(1,size(t82on,1)).*4.0 , 'r:')
                        hold all;
                        plot(t82on , ones(1,size(t82on,1)).*10.0 , 'r--')
                        
                        legend({'', 'UW: HIV, observed ART', 'DC: HIV, observed ART' , 'Elimination: <4/100K' , 'Benchmark: <10/100K' }, 'Location' , 'northeast');
                        xlim([1980 2121]); ylim([0 90]);
                        grid on;
                        xlabel('Year'); ylabel('AS-cervical cancer incidence (per 100K)');

%% Crude cervical cancer incidence by age
% Load adjusted Globocan 2018 rates for KZN
file = [pwd , '/Config/Reweighted_GlobocanCC_rates.xlsx'];
ccInc2018adjKZN(:,1) = xlsread(file , 'CC rates' , 'AB4:AB15');
% Globocan 2018 calibration error bars
meanObs = ccInc2018_dObs(: , 2);
sdevObs = (ccInc2018_dObs(: , 3).^(1/2)).*2;

fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
ageGroup = {'0-4' , '5-9' , '10-14' , '15-19' , '20-24' , '25-29' ,...
    '30-34' , '35-39' , '40-44' , '45-49' , '50-54' , '55-59' , ...
    '60-64' , '65-69' , '70-74' , '75-79' , '80+'};
year2021 = (2018-1982) + 1; 
t82on = (1982:2121)';
t82onLen = length(t82on);

errorbar(4 : age-1 , meanObs , sdevObs , 'rs' , 'LineWidth' , 1.5);
hold all;
plot(4 : age-1 , ccInc2018adjKZN , 'r*');
hold all;
for i = 1 : nTeams
%     for j = 3 : nSces 
    for j = 4 : 4; 
        for dInd = 1 : 3 %5
            % Load results
            fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
            ccIncTime = readmatrix(fname , 'Sheet' , 'CC incidence rate' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':T' , num2str(dInd*t82onLen+3)]);
            hold all;
            p = plot([1:length(ageGroup)] , ccIncTime(year2021,:) , 'Color' , colorVec{i} , 'LineWidth' , widthVec{dInd});
            set(gca , 'xtick' , 1 : length(ageGroup) , 'xtickLabel' , ageGroup);
        end
    end
end
legend('(Globocan, 2018) Observed SA: mean, 2SD' , 'Estimated KZN, adjusted Globocan 2018' , plotTits2{:} , 'Location' , 'northwest');
xlim([0 length(ageGroup)]); ylim([0 500]);
grid on;
xlabel('Age Group'); ylabel('Crude cervical cancer incidence'); title('S2: HIV, observed ART; Year: 2018')

%% Crude cumulative cervical cancer cases over time
% Note: need to check UW trends (more cumulative cases of CC in S0 compared
% to S1/2). This may be due to drastic difference in population size between
% scenarios due to large impact of HIV-associated mortality and no scale up 
% of ART until the 2000s.
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
t82on = (1982:2121)';
t82onLen = length(t82on);

xline(2021);
for i = 1 : nTeams
    for j = 1 : nSces  
        for dInd = 1 : 1 %5
            % Load results
            fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
                teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
            ccCasesTime = readmatrix(fname , 'Sheet' , 'CC case counts' , 'Range' , ['D' , num2str(((dInd-1)*t82onLen+1+3)) , ':S' , num2str(dInd*t82onLen+3)]);
            ccCasesTotTime = sum(ccCasesTime,2);
            ccCasesCumTime = cumsum(ccCasesTotTime);
            hold all;
            p = plot(t82on , ccCasesCumTime , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
        end
    end
end

legend(plotTits{:} , 'Location' , 'northwest');
xlim([1980 2121]); %ylim([0 10000]);
grid on;
xlabel('Year'); ylabel('Crude cumulative cervical cancer cases among women');

%% Crude HIV prevalence over time
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
t82on = (1982:2121)';
t82onLen = length(t82on);

xline(2021);
%hold all;
%plot(2016 , 0.233 , 'ro');
for i = 1 : nTeams
    for j = 1 : nSces  
        % Load results
        fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
            teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
        hivPrevTime = readmatrix(fname , 'Sheet' , 'HIV prevalence' , 'Range' , ['C' , num2str((1+3)) , ':C' , num2str(t82onLen+3)]);
        hold all;
        p = plot(t82on , hivPrevTime(:,1) , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
    end
end

legend('' , plotTits{2:end} , 'Location' , 'northeast');    %'(SA DHS) Observed SA, ages 15+' , 
xlim([1980 2121]); ylim([0 0.5]);
grid on;
xlabel('Year'); ylabel('Crude HIV prevalence among women');

%% Crude proportion on ART + VS over time
fig = figure;
set(fig,'DefaultAxesFontSize' , 18);
t82on = (1982:2121)';
t82onLen = length(t82on);

xline(2021);
hold all;
plot([(artYr(1:end-1) + 1) ; 2020] , [maxRateF(1:end-1) ; 0.627829] , 'ro');
for i = 1 : nTeams
    for j = 1 : nSces  
        % Load results
        fname = [pwd , '\HHCoM_Results\' , dirName , '\' , ...
            teamVec{i} , '_' , sceVec{j} , '_outcome_template_' , dateVec{i} , '.xlsx'];
        artCovTime = readmatrix(fname , 'Sheet' , 'ART coverage' , 'Range' , ['C' , num2str((1+3)) , ':C' , num2str(t82onLen+3)]);
        hold all;
        p = plot(t82on , artCovTime(:,1) , 'Color' , colorVec{i} , 'LineStyle' , styleVec{j});
    end
end

legend('' , 'Observed KZN, on ART + VS' , plotTits{2:end} , 'Location' , 'northeast');
xlim([1980 2121]); ylim([0 1.0]);
grid on;
xlabel('Year'); ylabel('Proportion on ART + VS among women');


