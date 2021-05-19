%function vaxCEA_multSims_IPVCsces()

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
% ***SET ME***: save names of potential scenarios to analyze as variables
baseDirName = 'Vaccine22Apr20Ph2V11_baseVax057_baseScreen_shortName_noVMMChpv_discontFxd_screenCovFxd_8ts-2021_';
dirName_reductBaseline = [baseDirName , 'WHO-SCES0b_6_1'];
dirName_reductBaseline2 = [baseDirName , 'WHO-SCES0b_gbV_6_1'];
dirName_SCE2 = [baseDirName , 'WHO-SCES2_6_1'];
dirName_SCE9 = [baseDirName , 'WHO-SCES6_6_1'];
dirName_SCE9gbV = [baseDirName , 'WHO-SCES6_gbV_6_1'];

% ***SET ME***: choose which scenarios you want to save data in Excel for
simVec = {dirName_reductBaseline , dirName_reductBaseline2 , ...
    dirName_SCE2 , dirName_SCE9 , dirName_SCE9gbV};
fileVec = {'sim0' , 'sim0' , 'sim2' , 'sim2' , 'sim2'};
% ***SET ME***: make sure the names here correspond to scenarios in simVec above
fileTits = {'Baseline w/ 57% girls vax' , ...
    '57% boys and girls vax' , ...
    '90% girls multi-cohort vax' , ...
    'WHO guidelines' , ...
    'WHO guidelines + 90% boys multi-cohort vax'};

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
    diseaseTits = {'General' , 'HIV-negative' , 'All WLHIV' , 'WLHIV-untreated' , 'WLHIV-on ART'};    
    
    for dInd = 1 : 3 %length(diseaseLabels)
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
        
        ccIncRefTot_lb = zeros(1 , size(ccIncHivAgeTime,2)-1);       
        for aInd = 1:age+4
            a = aInd;
            if aInd >= age
                a = age;
            end
            if aInd <= age    
                ccIncRef = ccIncHivAgeTime(a+1+age , 2:end) .* worldStandard_WP2015(aInd);
                if (dInd == 5) && (a < 3)
                    ccIncRef = zeros(1 , size(ccIncHivAgeTime,2)-1);
                end
            elseif aInd > age
                ccIncRef = ccIncHivAgeTime(a+1+age , 2:end);
                ccIncRef = [(ones(1,aInd-a).*ccIncRef(1,1)) , ccIncRef(1,1:end-(aInd-a))];
                ccIncRef = ccIncRef .* worldStandard_WP2015(aInd);
            end
            ccIncRefTot_lb = ccIncRefTot_lb + ccIncRef;
        end
        ccInc_lb = ccIncRefTot_lb ./ (sum(worldStandard_WP2015(1:age+4)));
        
        ccIncRefTot_ub = zeros(1 , size(ccIncHivAgeTime,2)-1);       
        for aInd = 1:age+4
            a = aInd;
            if aInd >= age
                a = age;
            end
            if aInd <= age    
                ccIncRef = ccIncHivAgeTime(a+1+age*2 , 2:end) .* worldStandard_WP2015(aInd);
                if (dInd == 5) && (a < 3)
                    ccIncRef = zeros(1 , size(ccIncHivAgeTime,2)-1);
                end
            elseif aInd > age
                ccIncRef = ccIncHivAgeTime(a+1+age*2 , 2:end);
                ccIncRef = [(ones(1,aInd-a).*ccIncRef(1,1)) , ccIncRef(1,1:end-(aInd-a))];
                ccIncRef = ccIncRef .* worldStandard_WP2015(aInd);
            end
            ccIncRefTot_ub = ccIncRefTot_ub + ccIncRef;
        end
        ccInc_ub = ccIncRefTot_ub ./ (sum(worldStandard_WP2015(1:age+4)));
          
        % Plot baseline incidence
        if (j == 1) && (dInd == 1)
            fig = figure;
            set(fig,'DefaultAxesFontSize' , 18);
        end
        subplot(1,3,dInd);
        hold all;
        firstYrInd = (currYear-1-startYear)+2;
        firstYrInd2 = (currYear-1-startYear)+1;
        p = plot(ccIncHivAgeTime(1 , firstYrInd:end) , ccInc(1 , firstYrInd2:end) , '-');
%         hold all;
%         x2 = [ccIncHivAgeTime(1 , firstYrInd:end) , fliplr(ccIncHivAgeTime(1 , firstYrInd:end))];
%         inBetween = [ccInc_ub(1 , firstYrInd2:end) , fliplr(ccInc_lb(1 , firstYrInd2:end))];
%         colorP = get(p,'Color');
%         h = fill(x2 , inBetween , colorP);
%         h.FaceAlpha = 0.3;
%         h.LineStyle = '--';
%         set(h,'EdgeColor', colorP);
        axis([2020 2121 0 200])
        grid on;
        xlabel('Year'); ylabel('AS ICC per 100K');
        title(diseaseTits{dInd});
        
        if j == nResults
            hold all;
            plot (ccIncHivAgeTime(1 , 2:end) , ones(size(ccIncHivAgeTime,2)-1,1).*4.0 , 'k:')
            hold all;
            plot (ccIncHivAgeTime(1 , 2:end) , ones(size(ccIncHivAgeTime,2)-1,1).*10.0 , 'k--')
        end  
    end
    legend(fileTits{:} , 'Elimination: <4/100K' , 'Benchmark: <10/100K');
end
