---------------------------------------------------------
MODEL BRANCHES:
---------------------------------------------------------
master
  -Setting: KwaZulu-Natal, South Africa
  -Base model
  
nonVaxInfections_DoArt
  -Setting: KwaZulu-Natal, South Africa
  -Used for DO ART CEA analysis
  -DoArtOutputs directory includes model outputs for DO ART analysis

nonVaxInfections_DoArt_diagHiv
  -Setting: KwaZulu-Natal, South Africa
  -Used for DO ART CEA analysis
  -Includes wrapper equations to calculate percent of PLWHIV diagnosed/undiagnosed
  -vaxCEA_plotSaveDoArt_071720.m includes code to plot/save DO ART analysis outputs
  
nonVaxInfections_Kenya
  -Setting: Kenya
  -CISNET Kenya model
  -Calibrating with ABC-SMC algorithm (in progress)
  -See Issue #107 for key differences with nonVaxInfections_Kenya_national

nonVaxInfections_Kenya_national
  -Setting: Kenya
  -Gui Liu dissertation model
  -Hand-calibrated
  -Assumes HPV affects HIV natural history

PAF_HIV
  -Setting: KwaZulu-Natal, South Africa
  -Used for PAF analysis

OptScrnAge
  -Setting: KwaZulu-Natal, South Africa
  -Used for optimal screening age analysis
  
  
---------------------------------------------------------
MODEL TAGS:
---------------------------------------------------------
KwaZulu-Natal model, multiple-HPV-types calibration:
Phase 2
  091520_backup_22Apr20Ph2V11calib
  083120_backup_22Apr20Ph2V10calib
  082520_backup_22Apr20Ph2V9calib
  082420_backup_22Apr20Ph2V8calib
  082120_backup_22Apr20Ph2V7calib
  080620_backup_22Apr20Ph2V6calib
  080320_backup_22Apr20Ph2V5calib
  072820_backup_22Apr20Ph2V4calib
  072120_backup_22Apr20Ph2V3calib
  071820_backup_22Apr20Ph2V2calib
  062520_backup_22Apr20Ph2calib
Phase 1
  010821_backup_22Apr20Ph1calib
Phase 0
  041620_backup_handCalibFrom24Feb20calib
Trial calibration attempt
  030420_backup_24Feb20calib
  022420_backup_21Feb20calib
  022120_backup_10Feb20calib

KwaZulu-Natal model outputs (presentation, abstract, result set, etc.):
  100720_output_whoCompare
  062320_output_IPVC2020
  101819_output_5yearAgeGroups-whoReport
  091019_output_uncertaintyAbstract
  090519_output_whoMeetingWupdates
  062719_output_whoMeeting
  053119_output_whoReportP1VaxScreen
  052019_output_cisnetMeetingHivSubgroup
  040919_output_whoReportArtCCexplore
  021919_output_whoReportHpv1910
  020619_output_whoReport

Intermediate KwaZulu-Natal model versions:
Main model predecessor with HPV modeled as single infection
(prior to HPV infection split into 9v and non-9v types, several other model revisions, and re-calibration)
  120319_backup_5yearAgeGroups
Version with single HPV infection, single age groups
  071619_backup_singleAgeVax
Version testing all differential equations solved per timestep
  100419_backup_allFuncsPerTimestep
Additional model backups, mostly before major merges or model additions/fixes
  040621_backup_nonVaxInfections-b4mergeWmaster
  040621_backup_master-b4caughtUpToNonVaxInfections
  012720_backup_b4ARTlims
  121719_backup_nonVaxInfections-retryMergeWnewHIVcd4vlStates
  121319_backup_nonVaxInfections-b4mergeWnewHIVcd4vlStates
  081919_backup_beforeMergeWautoCalibrate
  081619_backup_singleAgeGroups-b4MasterAutoCalMerge

KwaZulu-Natal model version left by programmer Nicholas Tan:
  010219_backup_afterNicholasTan
Version with some changes used for discarded CISNET comparative modeling analysis with Erasmus:
  051420_backup_after_nicholas_Tan_Erasmus_05102020
Stale branches left by programmer Nicholas Tan:
  120118_backup_4v
  120118_backup_omniHPV
  120118_backup_newLambda
  120118_backup_fastMix
  120118_backup_HPV_Smooth
  120118_backup_ART_fix

Intermediate Kenya model versions:
  090120_backup_KenyaWNyanza
  073120_backup_KenyaHandCalib
  031120_backup_Kenya


 
