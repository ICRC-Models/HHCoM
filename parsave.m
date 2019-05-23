function parsave(fname,tVec ,  popVec , newHiv ,...
    newImmHpv , newVaxHpv , newHpv , deaths , hivDeaths , ccDeath , ...
    newCC , artTreatTracker , vaxdLmtd , vaxdSchool , vaxdCU , newScreen , newTreatImm , newTreatHpv , newTreatHyst , ...
    ccTreated , currYear , lastYear , vaxRate , vaxEff , popLast, pathModifier)

savDir = [pwd , '\HHCoM_Results\Vaccine' , pathModifier,'\'];
save(fullfile(savDir , fname) , 'tVec' ,  'popVec' , 'newHiv' ,...
    'newImmHpv' , 'newVaxHpv' , 'newHpv' , 'deaths' , 'hivDeaths' , ...
    'ccDeath' , 'newCC' , 'artTreatTracker' , 'currYear' , 'lastYear' , ...
    'popLast' , 'ccTreated' , 'vaxdLmtd' , 'vaxdSchool' , 'vaxdCU' , ...
    'newScreen' , 'newTreatImm' , 'newTreatHpv' , 'newTreatHyst' , 'vaxRate' , 'vaxEff' , '-v7.3')

end
