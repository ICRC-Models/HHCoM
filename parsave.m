function parsave(fname, fivYrAgeGrpsOn , tVec , popVec , newHiv ,...
    newHpvVax , newImmHpvVax , newHpvNonVax , newImmHpvNonVax , ...
    hivDeaths , deaths , ccDeath , ...
    newCC , menCirc , vaxdLmtd , vaxdSchool , vaxdCU , newScreen , artDist , artDistList , ... %artTreatTracker , newTreatImm , newTreatHpv , newTreatHyst ,
    currYear , lastYear , vaxRate , vaxEff , popLast , pathModifier)

savDir = [pwd , '/HHCoM_Results/Vaccine' , pathModifier, '/'];

save(fullfile(savDir , fname) , 'fivYrAgeGrpsOn' , 'tVec' ,  'popVec' , 'newHiv' ,...
    'newHpvVax' , 'newImmHpvVax' , 'newHpvNonVax' , 'newImmHpvNonVax' , ...
    'hivDeaths' , 'deaths' , 'ccDeath' , ...
    'newCC' , 'menCirc' , 'vaxdLmtd' , 'vaxdSchool' , 'vaxdCU' , 'newScreen' , 'artDist' , 'artDistList' , ... %'artTreatTracker' ,'newTreatImm' , 'newTreatHpv' , 'newTreatHyst' , ... 
    'currYear' , 'lastYear' , 'vaxRate' , 'vaxEff' , 'popLast' , '-v7.3') 
end
