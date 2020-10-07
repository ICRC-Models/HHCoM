%function vaxCEA_plotSaveIncMort_032020()

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
    condUse , screenYrs , hpvScreenStartYear , waning , ...
    artYr , maxRateM , maxRateF , ...
    artYr_vec , artM_vec , artF_vec , minLim , maxLim , ...
    circ_aVec , vmmcYr_vec , vmmc_vec , vmmcYr , vmmcRate , ...
    hivStartYear , circStartYear , circNatStartYear , vaxStartYear , ...
    baseline , cisnet , who , whob , circProtect , condProtect , MTCTRate , ...
    hyst , OMEGA , ...
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
    infhpvNonVaxInds , ageInd , riskInd , ...
    hivNegNonVMMCinds , hivNegVMMCinds , ...
    vlAdvancer , ...
    fertMat , hivFertPosBirth , hivFertNegBirth , fertMat2 , ...
    hivFertPosBirth2 , hivFertNegBirth2 , fertMat3 , hivFertPosBirth3 , hivFertNegBirth3 , ...
    fertMat4 , hivFertPosBirth4 , hivFertNegBirth4 , ...
    dFertPos1 , dFertNeg1 , dFertMat1 , dFertPos2 , dFertNeg2 , dFertMat2 , ...
    dFertPos3 , dFertNeg3 , dFertMat3 , deathMat , deathMat2 , deathMat3 , deathMat4 , ...
    dDeathMat , dDeathMat2 , dDeathMat3 , dMue] = loadUp2(1 , 0 , [] , [] , []);

%% LOAD SAVED RESULTS
% ***SET ME***: save names of potential scenarios to analyze as variables
baseDirName = 'Vaccine22Apr20Ph2V11_noBaseVax_baseScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_';
dirName_reductBaseline = [baseDirName , 'WHO-SCES012_6_1'];
dirName_reductBaseline2 = ['Vaccine22Apr20Ph2V11_whoS0bBaseVax_baseScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_' , 'WHO-SCES012_6_1'];
dirName_reductBaseline3 = ['Vaccine22Apr20Ph2V11_whoS0cBaseVax_baseScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_' , 'WHO-SCES012_6_1'];
dirName_reductBaseline4 = ['Vaccine22Apr20Ph2V11_noBaseVax_noBaseScreen_hpvHIVcalib_adjFert2_adjCCAgeMults3_KZNCC4_noVMMChpv_' , 'WHO-SCES012_6_1'];
dirName_P1_SCE34 = [baseDirName , 'WHO-SCES34_6_1'];
dirName_P1_SCE56 = [baseDirName , 'WHO-SCES56_6_1'];
dirName_P2_SCE7a7 = [baseDirName , 'WHO-SCES7a7_6_1'];
dirName_P2_SCE7b = [baseDirName , 'WHO-SCES7b_6_1'];
dirName_P2_SCE8 = [baseDirName , 'WHO-SCES8_6_1'];
dirName_P2_SCE9 = [baseDirName , 'WHO-SCES9_6_1'];
dirName_P2_SCE10 = [baseDirName , 'WHO-SCES10_6_1'];
dirName_P2_SCE11 = [baseDirName , 'WHO-SCES11_6_1'];

% ***SET ME***: choose which scenarios you want to save data in Excel for
simVec = {dirName_reductBaseline , dirName_reductBaseline2 , ...
    dirName_reductBaseline3 , dirName_reductBaseline4 , ...
    dirName_reductBaseline , dirName_reductBaseline , ...
    dirName_P2_SCE7a7 , dirName_P2_SCE7b , dirName_P2_SCE7a7 , ...
    dirName_P2_SCE8};
fileVec = {'sim0' , 'sim0' , 'sim0' , 'sim0' , 'sim1' , 'sim2' , 'sim1' , 'sim1' , 'sim2' , 'sim2'}; 
% ***SET ME***: make sure the names here correspond to scenarios in simVec above
fileTits = {'S0 (no vax, baseline screen)' , 'S0b (60% 2v, baseline screen)' , ...
    'S0c (80% 2v, baseline screen)' , 'S0d (no vax, no screen)' , ...
    'P1-S1 (80% 9v, baseline screen)' , 'P1-S2 (90% 9v, baseline screen)' , ...
    'P2-S7a (80% 9v, 50% CU, baseline screen)' , 'P2-S7b (80% 9v, 80% CU, baseline screen)' , ...
    'P2-S7 (90% 9v, 50% CU, baseline screen)' , 'P2-S8 (90% 9v, 90% CU, baseline screen)'};

% % ***SET ME***: choose which scenarios you want to save data in Excel for
% simVec = {dirName_reductBaseline , dirName_reductBaseline4 , ...
%     dirName_reductBaseline , dirName_reductBaseline , ...
%     dirName_P1_SCE34 , dirName_P1_SCE34 , dirName_P1_SCE56 , dirName_P1_SCE56 , ...
%     dirName_P2_SCE9 , dirName_P2_SCE10 , dirName_P2_SCE11};
% fileVec = {'sim0' , 'sim0' , 'sim1' , 'sim2' , 'sim1' , 'sim2' , 'sim1' , 'sim2' , 'sim2' , 'sim2' , 'sim2'};
% % ***SET ME***: make sure the names here correspond to scenarios in simVec above
% fileTits = {'S0 (no vax, baseline screen)' , 'S0d (no vax, no screen)' , ...
%     'P1-S1 (80% 9v, baseline screen)' , 'P1-S2 (90% 9v, baseline screen)' , ...  
%     'P1-S3 (80% 9v, screen@35)' , 'P1-S4 (90% 9v, screen@35)' , ...
%     'P1-S5 (80% 9v, screen@35,45)' , 'P1-S6 (90% 9v, screen@35,45)' , ...
%     'P2-S9 (90% 9v, screen by HIV status)' , 'P2-S10(90% 9v, 50% CU, screen by HIV status)' , ...
%     'P2-S11(90% 9v, 90% CU, screen by HIV status)'};

% Plot settings
reset(0)
set(0 , 'defaultlinelinewidth' , 1.5)

nResults = length(simVec);

for j = 1 : nResults 
    %% CC INCIDENCE - age standardized
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
    
    diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};    
    
    for dInd = 1 : 1 %length(diseaseLabels);
        % Load results
        fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
            'ICC_' , diseaseLabels{dInd} , '_' , fileVec{j} , '.csv'];
        ccIncHivAgeTime = xlsread(fname);
        
        ccIncRefTot = zeros(1 , size(ccIncHivAgeTime,2)-1);       
        for aInd = 1:age+4
            a = aInd;
            if aInd >= age
                a = age;
            end
            
            if aInd <= age    
                ccIncRef = ccIncHivAgeTime(a+1 , 2:end) .* worldStandard_WP2015(aInd);
                if (dInd == 5) && (a < 3)
                    ccIncRef = zeros(1 , size(ccIncHivAgeTime,2)-1);
                end
            elseif aInd > age
                ccIncRef = ccIncHivAgeTime(a+1 , 2:end);
                ccIncRef = [(ones(1,aInd-a).*ccIncRef(1,1)) , ccIncRef(1,1:end-(aInd-a))];
                ccIncRef = ccIncRef .* worldStandard_WP2015(aInd);
            end
            ccIncRefTot = ccIncRefTot + ccIncRef;
        end
        ccInc = ccIncRefTot ./ (sum(worldStandard_WP2015(1:age+4)));
          
        % Plot baseline incidence
        if (j == 1) && (dInd == 1)
            figure;
        end
        hold all;
        plot(ccIncHivAgeTime(1 , 2:end) , ccInc , '-')
        axis([1980 2120 0 80])
        grid on;
        xlabel('Year'); ylabel('AS ICC per 100K'); 
        %legend('General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt');
    end
    
    %% HPV INCIDENCE - age standardized
%     % Note: the age-standardization process shifts the incidence rate of the
%     % last modelled age group to the next age group in the following year.
%     % However, HPV incidence is NaN prior to HIV introduction in the
%     % HIV-positive no ART group, and NaN prior to ART introduction in the
%     % HIV-positive ART group. Since we have four age groups past the 16 we
%     % model, a NaN value is present for four years past the introduction of
%     % HIV/ART, leading to a NaN value for summed incidence during these 
%     % years. We therefore lack data in this four-year interval in the
%     % saved/plotted results.
%
%     % Note: a limitation of these plots is that vax and nonVax HPV
%     % incidence is currently double-counted in coinfected individuals
%     
%     fac = 100;
%     worldStandard_WP2015 = [325428 311262 295693 287187 291738 299655 272348 ...
%         247167 240167 226750 201603 171975 150562 113118 82266 64484 42237 ...
%         23477 9261 2155];
%     
%     riskLabels = {'lr' , 'mr' , 'hr'};
%     riskPlotLines = {'-' , '--' , ':'};
%     diseaseLabels = {'General' , 'HIV_neg' , 'HIV_posAll' , 'HIV_posNoArt' , 'HIV_posArt'};    
%     for r = 1 : risk
%         set(gca,'ColorOrderIndex',1)
%         for dInd = 1 : length(diseaseLabels);
%             % Load results
%             fname = [pwd , '\HHCoM_Results\' , simVec{j} , '\' , ...
%                 'HPVinc_' , diseaseLabels{dInd} , '_' , riskLabels{r} , '_' , fileVec{j} , '.csv'];
%             hpvIncHivAgeTime = xlsread(fname);
% 
%             hpvIncRefTot = zeros(1 , size(hpvIncHivAgeTime,2)-1);       
%             for aInd = 1:age+4
%                 a = aInd;
%                 if aInd >= age
%                     a = age;
%                 end
% 
%                 if aInd <= age    
%                     hpvIncRef = hpvIncHivAgeTime(a+1 , 2:end) .* worldStandard_WP2015(aInd);
%                     if (a < 3)
%                         hpvIncRef = zeros(1 , size(hpvIncHivAgeTime,2)-1);
%                     end
%                 elseif aInd > age
%                     hpvIncRef = hpvIncHivAgeTime(a+1 , 2:end);
%                     hpvIncRef = [(ones(1,aInd-a).*hpvIncRef(1,1)) , hpvIncRef(1,1:end-(aInd-a))];
%                     hpvIncRef = hpvIncRef .* worldStandard_WP2015(aInd);
%                 end
%                 hpvIncRefTot = hpvIncRefTot + hpvIncRef;
%             end
%             hpvInc = hpvIncRefTot ./ (sum(worldStandard_WP2015(1:age+4)));
% 
%             % Plot baseline incidence
%             if (j == 1) && (dInd == 1) && (r == 1)
%                 figure;
%             end
%             hold all;
%             plot(hpvIncHivAgeTime(1 , 2:end) , hpvInc , riskPlotLines{r});
%             hold all;
%             axis([1980 2120 0 30])
%             grid on;
%             xlabel('Year'); ylabel('AS HPV incidence per 100'); 
%             legend('General, lr' , 'HIV-neg' , 'HIV-posAll' , 'HIV-posNoArt' , 'HIV-posArt' , ...
%                 'General, mr' , 'HIV-neg' , 'HIV-posAll' , 'HIV-posNoArt' , 'HIV-posArt' , ...
%                 'General, hr' , 'HIV-neg' , 'HIV-posAll' , 'HIV-posNoArt' , 'HIV-posArt');
%         end
%     end

end

%%
if j == nResults
    hold all;
    plot (ccIncHivAgeTime(1 , 2:end) , ones(size(ccIncHivAgeTime,2)-1,1).*4.0 , 'r:')
    hold all;
    plot (ccIncHivAgeTime(1 , 2:end) , ones(size(ccIncHivAgeTime,2)-1,1).*10.0 , 'r--')
    legend(fileTits{:} , ... %'General' , 'HIV-negative' , 'HIV-positive, all' , 'HIV-positive, untreated' , 'HIV-positive, on ART' , ...
       'Elimination: <4/100K' , 'Benchmark: <10/100K'); %
        %'P1-SCE0 (10sims) - adjusted KZN' , 'P1-SCE1 (10sims) - adjusted KZN' , 'P1-SCE2 (10sims) - adjusted KZN' , ...
end
